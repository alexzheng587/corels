#include "bbound.hh"
#include "time.hh"

BaseNode* base_construct_policy(size_t new_rule, size_t nrules, bool prediction,
                                bool default_prediction, double lower_bound,
                                double objective, BaseNode* parent,
                                int num_not_captured, int nsamples,
                                int len_prefix, double c) {
    (void) num_not_captured, nsamples, len_prefix, c;
    return (new BaseNode(new_rule, nrules, prediction, default_prediction,
                         lower_bound, objective, 0, parent));
}

CuriousNode* curious_construct_policy(size_t new_rule, size_t nrules, bool prediction,
                                      bool default_prediction, double lower_bound,
                                      double objective, CuriousNode* parent,
                                      int num_not_captured, int nsamples,
                                      int len_prefix, double c) {
    double curiosity = (lower_bound - c * len_prefix + c) * nsamples / (double)(nsamples - num_not_captured);
    return (new CuriousNode(new_rule, nrules, prediction, default_prediction,
                            lower_bound, objective, curiosity, parent));
}

template<class N>
N* prefix_permutation_insert(construct_signature<N> construct_policy, size_t new_rule,
                             size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                             double objective, N* parent, int num_not_captured, int nsamples, int len_prefix,
                             double c, CacheTree<N>* tree, VECTOR captured, std::vector<size_t> parent_prefix,
                             PrefixPermutationMap* p) {
    typename PrefixPermutationMap::iterator iter;
    parent_prefix.push_back(new_rule);
    std::set<size_t> key(parent_prefix.begin(), parent_prefix.end());
    N* child = NULL;
    iter = p->find(key);
    if (iter != p->end()) {
        std::vector<size_t> permuted_prefix = iter->second.first;
        double permuted_lower_bound = iter->second.second;
        if (lower_bound < permuted_lower_bound) {
            N* permuted_node;
            if ((permuted_node = tree->check_prefix(permuted_prefix)) != NULL) {
                N* permuted_parent = permuted_node->parent();
                permuted_parent->delete_child(permuted_node->id());
                delete_subtree<N>(tree, permuted_node, false);
            }
            child = construct_policy(new_rule, nrules, prediction, default_prediction,
                                       lower_bound, objective, parent,
                                        num_not_captured, nsamples, len_prefix, c);
            iter->second = std::make_pair(parent_prefix, lower_bound);
            //permutation_map_.insert(std::make_pair(key, child));
        }
    } else {
        child = construct_policy(new_rule, nrules, prediction, default_prediction,
                                    lower_bound, objective, parent,
                                    num_not_captured, nsamples, len_prefix, c);
        p->insert(std::make_pair(key, std::make_pair(parent_prefix, lower_bound)));
    }
    return child;
};

static
std::vector<bool> VECTOR_to_bitvector(VECTOR vec, size_t len) {
    std::vector<bool> bitvector;
    bitvector.resize(len);
    for(size_t i = 0; i < len / sizeof(unsigned long); ++i) {
        for(size_t j = 0; j < sizeof(unsigned long); ++j) {
            size_t index = (i * sizeof(unsigned long)) + j;
            size_t bmask = (1 << j) & vec[i];
            if (bmask == 0) {
                bitvector[index] = false;
            } else {
                bitvector[index] = true;
            }
        }
    }
    return bitvector;
}

template<class N>
N* captured_permutation_insert(construct_signature<N> construct_policy, size_t new_rule,
                               size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                               double objective, N* parent, int num_not_captured, int nsamples, int len_prefix,
                               double c, CacheTree<N>* tree, VECTOR captured, std::vector<size_t> parent_prefix,
                               CapturedPermutationMap* p) {
    typename CapturedPermutationMap::iterator iter;
    N* child = NULL;
    std::vector<bool> key = VECTOR_to_bitvector(captured, nsamples);
    iter = p->find(key);
    if (iter != p->end()) {
        std::vector<size_t> permuted_prefix = iter->second.first;
        double permuted_lower_bound = iter->second.second;
        if (lower_bound < permuted_lower_bound) {
            N* permuted_node;
            if ((permuted_node = tree->check_prefix(permuted_prefix)) != NULL) {
                N* permuted_parent = permuted_node->parent();
                permuted_parent->delete_child(permuted_node->id());
                delete_subtree<N>(tree, permuted_node, false);
            }
            child = construct_policy(new_rule, nrules, prediction, default_prediction,
                                       lower_bound, objective, parent,
                                        num_not_captured, nsamples, len_prefix, c);
            iter->second = std::make_pair(parent_prefix, lower_bound);
            //permutation_map_.insert(std::make_pair(key, child));
        }
    } else {
        child = construct_policy(new_rule, nrules, prediction, default_prediction,
                                    lower_bound, objective, parent,
                                    num_not_captured, nsamples, len_prefix, c);
        p->insert(std::make_pair(key, std::make_pair(parent_prefix, lower_bound)));
    }
    return child;
};


BaseNode* base_queue_front(BaseQueue* q) {
    return q->front();
}

CuriousNode* curious_queue_front(CuriousQueue* q) {
    return q->top();
}

template<class N, class Q, class P>
void evaluate_children(CacheTree<N>* tree, N* parent, VECTOR parent_not_captured,
                       std::set<size_t> ordered_parent,
                       construct_signature<N> construct_policy, Q* q, struct time* times,
                       permutation_insert_signature<N, P> permutation_insert, P* p) {
    std::vector<size_t> parent_prefix = parent->get_prefix();
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
    double t0 = timestamp();
    for (i = 1; i < nrules; i++) {
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
        times->lower_bound_time = time_diff(t1);
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
        times->objective_time += time_diff(t2);
        ++times->objective_num;
        if (objective < tree->min_objective()) {
            printf("min(objective): %1.5f -> %1.5f, length: %d, cache size: %zu\n",
                   tree->min_objective(), objective, len_prefix, tree->num_nodes());
            tree->update_min_objective(objective);
        }
        if ((lower_bound + c) < tree->min_objective()) {
            N* n;
            if (p) {
                double t3 = timestamp();
                n = permutation_insert(construct_policy, i, nrules, prediction, default_prediction,
                                       lower_bound, objective, parent, num_not_captured, nsamples,
                                       len_prefix, c, tree, captured, parent_prefix, p);
                times->permutation_map_insertion_time += time_diff(t3);
                ++times->permutation_map_insertion_num;
            }
            else
                n = construct_policy(i, nrules, prediction, default_prediction,
                                    lower_bound, objective, parent,
                                    num_not_captured, nsamples, len_prefix, c);
            if (n) {
                double t4 = timestamp();
                tree->insert(n);
                times->tree_insertion_time += time_diff(t4);
                ++times->tree_insertion_num;

                if (q) q->push(n);
            }
        }
    }
    times->rule_evaluation_time += time_diff(t0);
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

template<class N>
std::pair<N*, std::set<size_t> > stochastic_select(CacheTree<N>* tree, VECTOR not_captured) {
    typename std::map<size_t, N*>::iterator iter;
    N* node = tree->root();
    rule_copy(not_captured, tree->rule(node->id()).truthtable, tree->nsamples());
    int cnt;
    std::set<size_t> ordered_prefix;
    while (node->done()) {
        if ((node->lower_bound() + tree->c()) >= tree->min_objective()) {
            if (node->depth() > 0) {
                N* parent = node->parent();
                parent->delete_child(node->id());
                delete_subtree<N>(tree, node, true);
                if (parent->num_children() == 0)
                    tree->prune_up(parent);
            }
            return std::make_pair((N*) 0, std::set<size_t>{});
        }
        if (node->num_children() == 0) {
            tree->prune_up(node);
            return std::make_pair((N*) 0, std::set<size_t>{});
        }
        iter = node->random_child();
        node = iter->second;
        ordered_prefix.insert(iter->first);
        rule_vandnot(not_captured, not_captured, tree->rule(iter->first).truthtable, tree->nsamples(), &cnt);
    }
    return std::make_pair(node, ordered_prefix);
}

template<class N>
void bbound_stochastic(CacheTree<N>* tree, size_t max_num_nodes,
                       construct_signature<N> construct_policy, struct time* times) {
    std::pair<N*, std::set<size_t> > node_ordered;
    VECTOR not_captured;
    NullQueue<N>* q = NULL;

    double tot = timestamp();
    size_t num_iter = 0;
    rule_vinit(tree->nsamples(), &not_captured);
    tree->insert_root();
    while ((tree->num_nodes() < max_num_nodes) and (tree->num_nodes() > 0)) {
        double t0 = timestamp();
        node_ordered = stochastic_select<N>(tree, not_captured);
        times->node_select_time += timestamp() - t0;
        ++times->node_select_num;
        if (node_ordered.first) {
            double t1 = timestamp();
            evaluate_children<N, NullQueue<N>, PrefixPermutationMap>(tree, node_ordered.first, not_captured,
                                                                     node_ordered.second, construct_policy,
                                                                     q, times, NULL, NULL);
            times->evaluate_children_time += timestamp() - t1;
            ++times->evaluate_children_num;
        }
        ++num_iter;
        if ((num_iter % 10000) == 0)
            printf("num_iter: %zu, num_nodes: %zu\n", num_iter, tree->num_nodes());
    }
    times->total_time = timestamp() - tot;
    rule_vfree(&not_captured);
}

template<class N, class Q>
std::pair<N*, std::set<size_t> >
queue_select(CacheTree<N>* tree, Q* q, N*(*front)(Q*), VECTOR captured) {
    int cnt;

    N* selected_node = front(q); //q->front();
    q->pop();

    N* node = selected_node;
    if (node->deleted()) {  // lazily delete leaf nodes
        tree->decrement_num_nodes();
        delete node;
        return std::make_pair((N*) 0, std::set<size_t>{});
    }

    std::set<size_t> ordered_prefix;
    rule_vclear(tree->nsamples(), &captured);

    while (node != tree->root()) { /* or node->id() != root->id() */
        if (node->deleted()) {
            delete node;
            return std::make_pair((N*) 0, std::set<size_t>{});
        }
        ordered_prefix.insert(node->id());
        rule_vor(captured,
                 captured, tree->rule(node->id()).truthtable,
                 tree->nsamples(), &cnt);
        node = node->parent();
    }
    return std::make_pair(selected_node, ordered_prefix);
}

template<class N, class Q, class P>
void bbound_queue(CacheTree<N>* tree,
                size_t max_num_nodes,
                construct_signature<N> construct_policy,
                Q* q, N*(*front)(Q*),
                struct time* times, 
                permutation_insert_signature<N, P> permutation_insert,
                P* p) {
    int cnt;
    double min_objective = 1.0;
    std::pair<N*, std::set<size_t> > node_ordered;

    VECTOR captured, not_captured;
    rule_vinit(tree->nsamples(), &captured);
    rule_vinit(tree->nsamples(), &not_captured);

    double tot = timestamp();
    size_t num_iter = 0;

    tree->insert_root();
    ++times->tree_insertion_num;
    q->push(tree->root());
    while ((tree->num_nodes() < max_num_nodes) &&
           !q->empty()) {
        double t0 = timestamp();
        node_ordered = queue_select<N, Q>(tree, q, front, captured);
        times->node_select_time += time_diff(t0);
        ++times->node_select_num;
        if (node_ordered.first) {
            double t1 = timestamp();
            min_objective = tree->min_objective();
            /* not_captured = default rule truthtable & ~ captured */
            rule_vandnot(not_captured,
                         tree->rule(tree->root()->id()).truthtable, captured,
                         tree->nsamples(), &cnt);
            evaluate_children<N, Q, P>(tree, node_ordered.first, not_captured,
                                 node_ordered.second, construct_policy, q, times, 
                                 permutation_insert, p);
            times->evaluate_children_time += time_diff(t1);
            ++times->evaluate_children_num;
            if (tree->min_objective() < min_objective) {
                min_objective = tree->min_objective();
                printf("num_nodes before garbage_collect: %zu\n", tree->num_nodes());
                tree->garbage_collect();
                printf("num_nodes after garbage_collect: %zu\n", tree->num_nodes());
            }
        }
        ++num_iter;
        if ((num_iter % 10000) == 0)
            printf("iter: %zu, tree: %zu, queue: %zu, tree inserts: %d, permutation_insert(): %d\n",
                   num_iter, tree->num_nodes(), q->size(), times->tree_insertion_num, times->permutation_map_insertion_num);
    }

    if (q->empty())
        printf("Exited because queue empty\n");
    else
        printf("Exited because max number of nodes in the tree was reached\n"); 
    
    times->total_time = time_diff(tot);

    printf("Deleting queue elements and corresponding nodes in the cache, since they may not be reachable by the tree's destructor\n");
    N* node;
    while (!q->empty()) {
        node = front(q);
        q->pop();
        if (node->deleted()) {
            tree->decrement_num_nodes();
            delete node;
        }
    }
    rule_vfree(&captured);
    rule_vfree(&not_captured);
}

void bbound_greedy(size_t nsamples, size_t nrules, rule_t *rules, rule_t *labels,
                   size_t max_prefix_length) {
    // Initialize variables
    rule_t greedy_list[max_prefix_length];
    std::vector<rule_t> available_rules (rules, rules + nrules);
    VECTOR captured, captured_zeros, unseen; // not_captured, not_captured_zeros;
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

    rule_vfree(&captured);
    rule_vfree(&captured_zeros);
    rule_vfree(&unseen);
}

template<class N>
void delete_subtree(CacheTree<N>* tree, N* node, bool destructive) {
    N* child;
    typename std::map<size_t, N*>::iterator iter;
    if (node->done()) {
        iter = node->children_begin();
        while (iter != node->children_end()) {
            child = iter->second;
            delete_subtree<N>(tree, child, destructive);
            ++iter;
        }
        tree->decrement_num_nodes(); // always delete interior (non-leaf) nodes
        delete node;
    } else {
        if (destructive) {  // only delete leaf nodes in destructive mode
            tree->decrement_num_nodes();
            delete node;
        } else
            node->set_deleted();
    }
}

template BaseNode*
prefix_permutation_insert(construct_signature<BaseNode> construct_policy, size_t new_rule,
                          size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                          double objective, BaseNode* parent, int num_not_captured, int nsamples, int len_prefix,
                          double c, CacheTree<BaseNode>* tree, VECTOR captured, std::vector<size_t> parent_prefix,
                          PrefixPermutationMap* p);

template CuriousNode*
prefix_permutation_insert(construct_signature<CuriousNode> construct_policy, size_t new_rule,
                          size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                          double objective, CuriousNode* parent, int num_not_captured, int nsamples, int len_prefix,
                          double c, CacheTree<CuriousNode>* tree, VECTOR captured, std::vector<size_t> parent_prefix,
                          PrefixPermutationMap* p);

template BaseNode*
captured_permutation_insert(construct_signature<BaseNode> construct_policy, size_t new_rule,
                            size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                            double objective, BaseNode* parent, int num_not_captured, int nsamples, int len_prefix,
                            double c, CacheTree<BaseNode>* tree, VECTOR captured, std::vector<size_t> parent_prefix,
                            CapturedPermutationMap* p);

template CuriousNode*
captured_permutation_insert(construct_signature<CuriousNode> construct_policy, size_t new_rule,
                            size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                            double objective, CuriousNode* parent, int num_not_captured, int nsamples, int len_prefix,
                            double c, CacheTree<CuriousNode>* tree, VECTOR captured, std::vector<size_t> parent_prefix,
                            CapturedPermutationMap* p);

template void
evaluate_children<BaseNode, NullQueue<BaseNode>, PrefixPermutationMap>(CacheTree<BaseNode>* tree,
                                                  BaseNode* parent,
                                                  VECTOR parent_not_captured,
                                                  std::set<size_t> ordered_parent,
                                                  construct_signature<BaseNode> construct_policy,
                                                  NullQueue<BaseNode>* q, struct time*, 
                                                  permutation_insert_signature<BaseNode, PrefixPermutationMap> permutation_insert,
                                                  PrefixPermutationMap* p);

template void
evaluate_children<BaseNode, BaseQueue, PrefixPermutationMap>(CacheTree<BaseNode>* tree,
                                       BaseNode* parent,
                                       VECTOR parent_not_captured,
                                       std::set<size_t> ordered_parent,
                                       construct_signature<BaseNode> construct_policy,
                                       BaseQueue* q, struct time*,
                                       permutation_insert_signature<BaseNode, PrefixPermutationMap> permutation_insert,
                                       PrefixPermutationMap* p);

template void
evaluate_children<BaseNode, BaseQueue, CapturedPermutationMap>(CacheTree<BaseNode>* tree,
                                       BaseNode* parent,
                                       VECTOR parent_not_captured,
                                       std::set<size_t> ordered_parent,
                                       construct_signature<BaseNode> construct_policy,
                                       BaseQueue* q, struct time*,
                                       permutation_insert_signature<BaseNode, CapturedPermutationMap> permutation_insert,
                                       CapturedPermutationMap* p);

template void
evaluate_children<CuriousNode, CuriousQueue, PrefixPermutationMap>(CacheTree<CuriousNode>* tree,
                                             CuriousNode* parent,
                                             VECTOR parent_not_captured,
                                             std::set<size_t> ordered_parent,
                                             construct_signature<CuriousNode> construct_policy,
                                             CuriousQueue* q, struct time*,
                                             permutation_insert_signature<CuriousNode, PrefixPermutationMap> permutation_insert,
                                             PrefixPermutationMap* p);

template std::pair<BaseNode*, std::set<size_t> >
stochastic_select<BaseNode>(CacheTree<BaseNode>* tree, VECTOR not_captured); 

template void
bbound_stochastic<BaseNode>(CacheTree<BaseNode>* tree,
                            size_t max_num_nodes,
                            construct_signature<BaseNode> construct_policy,
                            struct time*);

template std::pair<BaseNode*, std::set<size_t> >
queue_select<BaseNode, BaseQueue>(CacheTree<BaseNode>* tree,
                                  BaseQueue* q,
                                  BaseNode*(*front)(BaseQueue*),
                                  VECTOR captured);

template std::pair<CuriousNode*, std::set<size_t> >
queue_select<CuriousNode, CuriousQueue>(CacheTree<CuriousNode>* tree,
                                        CuriousQueue* q,
                                        CuriousNode*(*front)(CuriousQueue*),
                                        VECTOR captured);

template void
bbound_queue<BaseNode, BaseQueue, PrefixPermutationMap>(CacheTree<BaseNode>* tree,
                                  size_t max_num_nodes,
                                  construct_signature<BaseNode> construct_policy,
                                  BaseQueue* q,
                                  BaseNode*(*front)(BaseQueue*),
                                  struct time*, 
                                  permutation_insert_signature<BaseNode, PrefixPermutationMap> permutation_insert,
                                  PrefixPermutationMap* p);

template void
bbound_queue<BaseNode, BaseQueue, CapturedPermutationMap>(CacheTree<BaseNode>* tree,
                                  size_t max_num_nodes,
                                  construct_signature<BaseNode> construct_policy,
                                  BaseQueue* q,
                                  BaseNode*(*front)(BaseQueue*),
                                  struct time*, 
                                  permutation_insert_signature<BaseNode, CapturedPermutationMap> permutation_insert,
                                  CapturedPermutationMap* p);

template void
bbound_queue<CuriousNode, CuriousQueue, PrefixPermutationMap>(CacheTree<CuriousNode>* tree,
                                        size_t max_num_nodes,
                                        construct_signature<CuriousNode> construct_policy,
                                        CuriousQueue* q,
                                        CuriousNode*(*front)(CuriousQueue*),
                                        struct time*,
                                        permutation_insert_signature<CuriousNode, PrefixPermutationMap> permutation_insert,
                                        PrefixPermutationMap* p);

template void
bbound_queue<CuriousNode, CuriousQueue, CapturedPermutationMap>(CacheTree<CuriousNode>* tree,
                                        size_t max_num_nodes,
                                        construct_signature<CuriousNode> construct_policy,
                                        CuriousQueue* q,
                                        CuriousNode*(*front)(CuriousQueue*),
                                        struct time*,
                                        permutation_insert_signature<CuriousNode, CapturedPermutationMap> permutation_insert,
                                        CapturedPermutationMap* p);

template void
delete_subtree<BaseNode>(CacheTree<BaseNode>* tree, BaseNode* n, bool destructive);

template void
delete_subtree<CuriousNode>(CacheTree<CuriousNode>* tree, CuriousNode* n, bool destructive);
