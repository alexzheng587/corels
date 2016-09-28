#include "bbound.hh"


BaseNode* base_construct_policy(size_t new_rule, size_t nrules, bool prediction,
                                bool default_prediction, double lower_bound,
                                double objective, BaseNode* parent,
                                int num_not_captured, int nsamples,
                                int len_prefix, double c) {
    return (new BaseNode(new_rule, nrules, prediction, default_prediction, lower_bound, objective, 0, parent));
}

CuriousNode* curious_construct_policy(size_t new_rule, size_t nrules, bool prediction,
                                      bool default_prediction, double lower_bound,
                                      double objective, CuriousNode* parent,
                                      int num_not_captured, int nsamples,
                                      int len_prefix, double c) {
    double curiosity = (lower_bound - c * len_prefix) * nsamples / (float)(nsamples - num_not_captured) + c * len_prefix;
    return (new CuriousNode(new_rule, nrules, prediction, default_prediction, lower_bound, objective, curiosity, parent));
}

template<class T>
void evaluate_children(CacheTree<T>* tree, Node<T>* parent, VECTOR parent_not_captured, construct_signature<T> construct_policy) {
    VECTOR captured, captured_zeros, not_captured, not_captured_zeros;
    int num_captured, c0, c1, captured_correct;
    int num_not_captured, d0, d1, default_correct;
    bool prediction, default_prediction;
    double lower_bound, objective, parent_lower_bound;
    int nsamples = tree->nsamples();
    int nrules = tree->nrules();
    double c = tree->c();
    rule_vinit(nsamples, &captured);
    rule_vinit(nsamples, &captured_zeros);
    rule_vinit(nsamples, &not_captured);
    rule_vinit(nsamples, &not_captured_zeros);
    int i, len_prefix;
    len_prefix = parent->depth() + 1;
    parent_lower_bound = parent->lower_bound();
    for (i = 1; i < nrules; i++) {
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
        lower_bound = parent_lower_bound + (float)(num_captured - captured_correct) / nsamples + c;
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
        if (objective < tree->min_objective()) {
            printf("min(objective): %1.5f -> %1.5f, length: %d, cache size: %zu\n",
                   tree->min_objective(), objective, len_prefix, tree->num_nodes());
            tree->update_min_objective(objective);
        }
        if ((lower_bound + c) < tree->min_objective()) {
            tree->insert(construct_policy(i, nrules, prediction, default_prediction, lower_bound,
                                          objective, parent, num_not_captured, nsamples, len_prefix, c));
            //tree->insert(i, prediction, default_prediction, lower_bound, objective, parent, num_not_captured);
        }
    }
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

template<class T>
Node<T>* stochastic_select(CacheTree<T>* tree, VECTOR not_captured) {
    typename std::map<size_t, Node<T>*>::iterator iter;
    Node<T>* node = tree->root();
    rule_copy(not_captured, tree->rule(node->id()).truthtable, tree->nsamples());
    int cnt;
    while (node->done()) {
        if ((node->lower_bound() + tree->c()) >= tree->min_objective()) {
            if (node->depth() > 0) {
                Node<T>* parent = node->parent();
                parent->delete_child(node->id());
                tree->delete_subtree(node);
            }
            return 0;
        }
        if (node->num_children() == 0) {
            tree->prune_up(node);
            return 0;
        }
        iter = node->random_child();
        node = iter->second;
        rule_vandnot(not_captured, not_captured, tree->rule(iter->first).truthtable, tree->nsamples(), &cnt);
    }
    return node;
}

template<class T>
void bbound_stochastic(CacheTree<T>* tree, size_t max_num_nodes, construct_signature<T> construct_policy) {
    Node<T>* node;
    VECTOR not_captured;
    size_t num_iter = 0;
    rule_vinit(tree->nsamples(), &not_captured);
    tree->insert_root();
    while ((tree->num_nodes() < max_num_nodes) and (tree->num_nodes() > 0)) {
        node = stochastic_select<T>(tree, not_captured);
        if (node)
            evaluate_children<T>(tree, node, not_captured, construct_policy);
        ++num_iter;
        if ((num_iter % 10000) == 0)
            printf("num_iter: %zu, num_nodes: %zu\n", num_iter, tree->num_nodes());
    }
    rule_vfree(&not_captured);
}

template void evaluate_children<bool>(CacheTree<bool>* tree, BaseNode* parent, VECTOR parent_not_captured, construct_signature<bool> construct_policy);

template BaseNode* stochastic_select<bool>(CacheTree<bool>* tree, VECTOR not_captured);

template void bbound_stochastic<bool>(CacheTree<bool>* tree, size_t max_num_nodes, construct_signature<bool> construct_policy);

template void evaluate_children<double>(CacheTree<double>* tree, CuriousNode* parent, VECTOR parent_not_captured, construct_signature<double> construct_policy);

template CuriousNode* stochastic_select<double>(CacheTree<double>* tree, VECTOR not_captured);

template void bbound_stochastic<double>(CacheTree<double>* tree, size_t max_num_nodes, construct_signature<double> construct_policy);

