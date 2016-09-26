#include "bbound.hh"
#include "time.hh"

struct time *times;

void evaluate_children(CacheTree* tree, CacheNode* parent, VECTOR parent_not_captured, std::set<size_t> ordered_parent) {
    VECTOR captured, captured_zeros, not_captured, not_captured_zeros;
    int num_captured, c0, c1, captured_correct;
    int num_not_captured, d0, d1, default_correct;
    bool prediction, default_prediction;
    double lower_bound, objective, parent_lower_bound;
    int nsamples = tree->nsamples();
    double c = tree->c();
    rule_vinit(nsamples, &captured);
    rule_vinit(nsamples, &captured_zeros);
    rule_vinit(nsamples, &not_captured);
    rule_vinit(nsamples, &not_captured_zeros);
    size_t i, len_prefix;
    len_prefix = parent->depth() + 1;
    parent_lower_bound = parent->lower_bound();
    double t0 = timestamp();
    for (i = 1; i < tree->nrules(); i++) {
        double t1 = timestamp();
        if (ordered_parent.find(i) != ordered_parent.end())
            continue;
        rule_vand(captured, parent_not_captured, tree->rule(i).truthtable, nsamples, &num_captured);
        rule_vand(captured_zeros, captured, tree->label(0).truthtable, nsamples, &c0);
        c1 = num_captured - c0;
        if (c0 > c1) {
            prediction = 0;
            captured_correct = c0;
        } else {
            prediction = 1;
            captured_correct = c1;
        }
        if (captured_correct < (c * nsamples))
            continue;
        lower_bound = parent_lower_bound + (float)(num_captured - captured_correct) / nsamples + c;
        times->lower_bound_time = timestamp() - t1;
        ++times->lower_bound_num;
        double t2 = timestamp();
        rule_vandnot(not_captured, parent_not_captured, captured, nsamples, &num_not_captured);
        rule_vand(not_captured_zeros, not_captured, tree->label(0).truthtable, nsamples, &d0);
        d1 = num_not_captured - d0;
        if (d0 > d1) {
            default_prediction = 0;
            default_correct = d0;
        } else {
            default_prediction = 1;
            default_correct = d1;
        }
        objective = lower_bound + (float)(num_not_captured - default_correct) / nsamples;
        times->objective_time += timestamp() - t2;
        ++times->objective_num;
        if (objective < tree->min_objective()) {
            printf("min(objective): %1.5f -> %1.5f, length: %zu, cache size: %zu\n",
                   tree->min_objective(), objective, len_prefix, tree->num_nodes());
            tree->update_min_objective(objective);
        }
        if ((lower_bound + c) < tree->min_objective()) {
            double t3 = timestamp();
            tree->insert(i, prediction, default_prediction, lower_bound, objective, parent);
            times->tree_insertion_time += timestamp() - t3;
            ++times->tree_insertion_num;
        }
    }
    times->rule_evaluation_time += timestamp() - t0;
    ++times->rule_evaluation_num;
    if (parent->num_children() == 0) {
        tree->prune_up(parent);
    } else {
        parent->set_done();
        tree->increment_num_evaluated();
    }
    rule_vfree(&captured);
    rule_vfree(&captured_zeros);
    rule_vfree(&not_captured);
    rule_vfree(&not_captured_zeros);
}

std::pair<CacheNode*, std::set<size_t> > stochastic_select(CacheTree* tree, VECTOR not_captured) {
    std::map<size_t, CacheNode*>::iterator iter;
    CacheNode* node = tree->root();
    rule_copy(not_captured, tree->rule(node->id()).truthtable, tree->nsamples());
    int cnt;
    std::set<size_t> ordered_prefix;
    while (node->done()) {
        if ((node->lower_bound() + tree->c()) >= tree->min_objective()) {
            if (node->depth() > 0) {
                CacheNode* parent = node->parent();
                parent->delete_child(node->id());
                tree->delete_subtree(node);
            }
            return std::make_pair((CacheNode*) 0, ordered_prefix);
        }
        if (node->num_children() == 0) {
            tree->prune_up(node);
            return std::make_pair((CacheNode*) 0, ordered_prefix);
        }
        iter = node->random_child();
        node = iter->second;
        ordered_prefix.insert(iter->first);
        rule_vandnot(not_captured, not_captured, tree->rule(iter->first).truthtable, tree->nsamples(), &cnt);
    }
    return std::make_pair(node, ordered_prefix);
}

struct time* bbound_stochastic(CacheTree* tree, size_t max_num_nodes) {
    std::pair<CacheNode*, std::set<size_t> > node_ordered;
    VECTOR not_captured;
    times = (struct time*) calloc(1, sizeof(*times));
    double tot = timestamp();
    size_t num_iter = 0;
    rule_vinit(tree->nsamples(), &not_captured);
    tree->insert_root();
    while ((tree->num_nodes() < max_num_nodes) and (tree->num_nodes() > 0)) {
        double t0 = timestamp();
        node_ordered = stochastic_select(tree, not_captured);
        times->stochastic_select_time += timestamp() - t0;
        ++times->stochastic_select_num;
        if (node_ordered.first) {
            double t1 = timestamp();
            evaluate_children(tree, node_ordered.first, not_captured, node_ordered.second);
            times->evaluate_children_time += timestamp() - t1;
            ++times->evaluate_children_num;
        }
        ++num_iter;
        if ((num_iter % 10000) == 0)
            printf("num_iter: %zu, num_nodes: %zu\n", num_iter, tree->num_nodes());
    }
    times->total_time = timestamp() - tot;
    rule_vfree(&not_captured);
    return times;
}

void bbound_greedy(size_t nsamples, size_t nrules, rule_t *rules, rule_t *labels, 
        size_t max_prefix_length) {
    // Initialize variables
    rule_t greedy_list[max_prefix_length];
    std::vector<rule_t> available_rules (rules, rules + nrules);
    VECTOR captured, captured_zeros, unseen;// not_captured, not_captured_zeros;
    int num_captured, c0, c1, prediction, captured_correct;
    rule_vinit(nsamples, &captured);
    rule_vinit(nsamples, &captured_zeros);
    rule_vinit(nsamples, &unseen);
    make_default(&unseen, nsamples);

    for(size_t i = 0; i < max_prefix_length; ++i) {
        float best_percent_captured = 0.0;
        int best_num_captured = 0;
        int best_index = 0;
        // Iterate over rule array to find best rule
        for(size_t j = 0; j < available_rules.size(); ++j) {
            rule_vand(captured, available_rules[j].truthtable, unseen, nsamples, &num_captured);
            rule_vand(captured_zeros, captured, labels[0].truthtable, nsamples, &c0);
            c1 = num_captured - c0;
            if (c0 > c1) {
                prediction = 0;
                captured_correct = c0;
            } else {
                prediction = 1;
                captured_correct = c1;
            }
            float percent_captured = (float)captured_correct / num_captured;
            if (percent_captured > best_percent_captured || (percent_captured == best_percent_captured && num_captured > best_num_captured)) {
                best_percent_captured = percent_captured;
                best_num_captured = num_captured;
                best_index = j;
            }
        }
        rule_t best_rule = available_rules[best_index];
        // Update unseen with best rule so far
        rule_vor(captured, unseen, best_rule.truthtable, nsamples, &c0);
        available_rules.erase(available_rules.begin() + best_index);
        greedy_list[i] = best_rule;
    }
    rule_print_all(greedy_list, max_prefix_length, nsamples);
}
