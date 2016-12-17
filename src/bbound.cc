#include "bbound.hh"
#include "utils.hh"

BaseNode* base_construct_policy(size_t new_rule, size_t nrules, bool prediction,
                                bool default_prediction, double lower_bound,
                                double objective, BaseNode* parent,
                                int num_not_captured, int nsamples,
                                int len_prefix, double c) {
    (void) len_prefix, c;
    size_t num_captured = nsamples - num_not_captured; // number captured by the new rule
    return (new BaseNode(new_rule, nrules, prediction, default_prediction,
                         lower_bound, objective, 0, parent, num_captured));
}

CuriousNode* curious_construct_policy(size_t new_rule, size_t nrules, bool prediction,
                                      bool default_prediction, double lower_bound,
                                      double objective, CuriousNode* parent,
                                      int num_not_captured, int nsamples,
                                      int len_prefix, double c) {
    size_t num_captured = nsamples - num_not_captured; // number captured by the new rule
    double curiosity = (lower_bound - c * len_prefix + c) * nsamples / (double)(num_captured);
    return (new CuriousNode(new_rule, nrules, prediction, default_prediction,
                            lower_bound, objective, curiosity, parent, num_captured));
}

void prefix_map_garbage_collect(PrefixPermutationMap* p, size_t queue_min_length) {
    typename PrefixPermutationMap::iterator iter;
    size_t num_deleted = 0;
    printf("pmap gc for length %zu: %zu -> ", queue_min_length, p->size());
    for (iter = p->begin(); iter != p->end(); ) {
        if (iter->first.size() <= queue_min_length) {
            iter = p->erase(iter);
            ++num_deleted;
        } else {
            ++iter;
        }
    }
    printf("%zu\n", p->size());
    logger.decreasePmapSize(num_deleted);
}

void bfs_prefix_map_garbage_collect(PrefixPermutationMap* p, size_t queue_min_length) {
    size_t pmap_size = p->size();
    printf("bfs pmap gc for length %zu: %zu -> ", queue_min_length, pmap_size);
    p->clear();
    printf("%zu\n", p->size());
    logger.decreasePmapSize(pmap_size);
}

void captured_map_garbage_collect(CapturedPermutationMap* p, size_t min_length) {
}

template<class N>
N* prefix_permutation_insert(construct_signature<N> construct_policy, size_t new_rule,
                             size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                             double objective, N* parent, int num_not_captured, int nsamples, int len_prefix,
                             double c, CacheTree<N>* tree, VECTOR not_captured, std::vector<size_t> parent_prefix,
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
                delete_subtree<N>(tree, permuted_node, false, true);
                logger.incPmapDiscardNum();
            } else {
                logger.incPmapNullNum();
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
        logger.incPmapSize();
    }
    return child;
};

static
std::vector<bool> VECTOR_to_bitvector(VECTOR vec, size_t len) {
    std::vector<bool> bitvector;
    bitvector.resize(len);
    for (size_t index = 0; index < len; index++) {
        size_t i = index / BITS_PER_ENTRY;
        size_t j = (index % BITS_PER_ENTRY);
        size_t bmask = (1 << j) & vec[i];
        if (bmask != 0)
            bitvector[index] = true;
        else
            bitvector[index] = false;
    }
    return bitvector;
}

template<class N>
N* captured_permutation_insert(construct_signature<N> construct_policy, size_t new_rule,
                               size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                               double objective, N* parent, int num_not_captured, int nsamples, int len_prefix,
                               double c, CacheTree<N>* tree, VECTOR not_captured, std::vector<size_t> parent_prefix,
                               CapturedPermutationMap* p) {
    typename CapturedPermutationMap::iterator iter;
    parent_prefix.push_back(new_rule);
    N* child = NULL;
    std::vector<bool> key = VECTOR_to_bitvector(not_captured, nsamples);
    iter = p->find(key);
    if (iter != p->end()) {
        std::vector<size_t> permuted_prefix = iter->second.first;
        double permuted_lower_bound = iter->second.second;
        if (lower_bound < permuted_lower_bound) {
            N* permuted_node;
            if ((permuted_node = tree->check_prefix(permuted_prefix)) != NULL) {
                N* permuted_parent = permuted_node->parent();
                permuted_parent->delete_child(permuted_node->id());
                delete_subtree<N>(tree, permuted_node, false, true);
                logger.incPmapDiscardNum();
            } else {
                logger.incPmapNullNum();
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
        logger.incPmapSize();
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
                       construct_signature<N> construct_policy, Q* q,
                       permutation_insert_signature<N, P> permutation_insert, P* p) {
    auto pp_pair = parent->get_prefix_and_predictions();
    std::vector<size_t> parent_prefix = std::move(pp_pair.first);
    std::vector<bool> parent_predictions = std::move(pp_pair.second);

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
        // captured represents data captured by the new rule
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
        logger.setLowerBoundTime(time_diff(t1));
        logger.incLowerBoundNum();
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
        logger.addToObjTime(time_diff(t2));
        logger.incObjNum();
        if (objective < tree->min_objective()) {
            printf("min(objective): %1.5f -> %1.5f, length: %d, cache size: %zu\n",
                   tree->min_objective(), objective, len_prefix, tree->num_nodes());

            tree->update_min_objective(objective);
            tree->update_opt_rulelist(parent_prefix, i);
            tree->update_opt_predictions(parent_predictions, prediction, default_prediction);

            logger.dumpState(); // dump state when min objective is updated (keep this)
        }
        if ((lower_bound + c) < tree->min_objective()) {
            N* n;
            if (p) {
                double t3 = timestamp();
                n = permutation_insert(construct_policy, i, nrules, prediction, default_prediction,
                                       lower_bound, objective, parent, num_not_captured, nsamples,
                                       len_prefix, c, tree, not_captured, parent_prefix, p);
                logger.addToPermMapInsertionTime(time_diff(t3));
                logger.incPermMapInsertionNum();
            }
            else
                n = construct_policy(i, nrules, prediction, default_prediction,
                                    lower_bound, objective, parent,
                                    num_not_captured, nsamples, len_prefix, c);
            if (n) {
                double t4 = timestamp();
                tree->insert(n);
                logger.addToTreeInsertionTime(time_diff(t4));
                logger.incTreeInsertionNum();
                logger.incPrefixLen(len_prefix);
                if (q) {
                    q->push(n);
                    logger.setQueueSize(q->size());
                    logger.addQueueElement(len_prefix, lower_bound);
                }
            }
        }
    }

    rule_vfree(&captured);
    rule_vfree(&captured_zeros);
    rule_vfree(&not_captured);
    rule_vfree(&not_captured_zeros);

    logger.addToRuleEvalTime(time_diff(t0));
    logger.incRuleEvalNum();
    logger.decPrefixLen(parent->depth());
    logger.removeQueueElement(len_prefix - 1, parent_lower_bound);
    if (parent->num_children() == 0) {
        tree->prune_up(parent);
    } else {
        parent->set_done();
        tree->increment_num_evaluated();
    }
    //logger.dumpState(); // dump state at least once for each call to evaluate_children (possibly too frequent)
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
                delete_subtree<N>(tree, node, true, true);
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
    logger.setCurrentLowerBound(node->lower_bound() + tree->c());
    return std::make_pair(node, ordered_prefix);
}

template<class N, class P>
void bbound_stochastic(CacheTree<N>* tree, size_t max_num_nodes,
                       construct_signature<N> construct_policy,
                       permutation_insert_signature<N, P> permutation_insert,
                       pmap_garbage_collect_signature<P> pmap_garbage_collect,
                       P* p) {
    double start = timestamp();
    logger.setInitialTime(start);
    logger.dumpState();         // initial log record

    double min_objective = 1.0;
    std::pair<N*, std::set<size_t> > node_ordered;
    VECTOR not_captured;
    NullQueue<N>* q = NULL;
    logger.incPrefixLen(0);

    rule_vinit(tree->nsamples(), &not_captured);

    size_t num_iter = 0;
    tree->insert_root();
    logger.incTreeInsertionNum();
    q->push(tree->root());
    logger.dumpState();         // log record for empty rule list
    while ((tree->num_nodes() < max_num_nodes) and (tree->num_nodes() > 0)) {
        double t0 = timestamp();
        node_ordered = stochastic_select<N>(tree, not_captured);
        logger.addToNodeSelectTime(time_diff(t0));
        logger.incNodeSelectNum();
        if (node_ordered.first) {
            double t1 = timestamp();
            min_objective = tree->min_objective();
            evaluate_children<N, NullQueue<N>, P>(tree, node_ordered.first, not_captured,
                                                  node_ordered.second, construct_policy,
                                                  q, permutation_insert, p);
            logger.addToEvalChildrenTime(time_diff(t1));
            logger.incEvalChildrenNum();

            if (tree->min_objective() < min_objective) {
                min_objective = tree->min_objective();
                printf("num_nodes before garbage_collect: %zu\n", tree->num_nodes());
                logger.dumpState();
                tree->garbage_collect();
                logger.dumpState();
                printf("num_nodes after garbage_collect: %zu\n", tree->num_nodes());
            }
        }
        ++num_iter;
        if ((num_iter % 10000) == 0) {
            if (p)
                printf("iter: %zu, tree: %zu, queue: %zu, pmap: %zu, time elapsed: %f\n",
                       num_iter, tree->num_nodes(), q->size(), p->size(), time_diff(start));
            else
                printf("iter: %zu, tree: %zu, queue: %zu, time elapsed: %f\n",
                       num_iter, tree->num_nodes(), q->size(), time_diff(start));
        }
        if ((num_iter % logger.getFrequency()) == 0)
            logger.dumpState();     // want ~1000 records for detailed figures
    }
    logger.dumpState();
    rule_vfree(&not_captured);
}

template<class N, class Q>
std::pair<N*, std::set<size_t> >
queue_select(CacheTree<N>* tree, Q* q, N*(*front)(Q*), VECTOR captured) {
    int cnt;

    N* selected_node = front(q); //q->front();
    q->pop();
    logger.setCurrentLowerBound(selected_node->lower_bound() + tree->c());

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
int bbound_queue(CacheTree<N>* tree,
                size_t max_num_nodes,
                construct_signature<N> construct_policy,
                Q* q, N*(*front)(Q*),
                permutation_insert_signature<N, P> permutation_insert,
                pmap_garbage_collect_signature<P> pmap_garbage_collect,
                P* p, size_t num_iter) {
    double start;
    int cnt;
    double min_objective;
    std::pair<N*, std::set<size_t> > node_ordered;
    VECTOR captured, not_captured;
    rule_vinit(tree->nsamples(), &captured);
    rule_vinit(tree->nsamples(), &not_captured);

    size_t queue_min_length = logger.getQueueMinLen();

    if (tree->num_nodes() == 0) {
        start = timestamp();
        logger.setInitialTime(start);
        logger.initializeState();
        logger.dumpState();         // initial log record

        min_objective = 1.0;
        tree->insert_root();
        logger.incTreeInsertionNum();
        q->push(tree->root());
        logger.setQueueSize(q->size());
        logger.incPrefixLen(0);
        logger.dumpState();         // log record for empty rule list
        logger.initRemainingSpaceSize();
    } else {
        start = logger.getInitialTime();
        min_objective = tree->min_objective();
    }
    while ((tree->num_nodes() < max_num_nodes) &&
           !q->empty()) {
        double t0 = timestamp();
        node_ordered = queue_select<N, Q>(tree, q, front, captured);
        logger.addToNodeSelectTime(time_diff(t0));
        logger.incNodeSelectNum();
        if (node_ordered.first) {
            double t1 = timestamp();
            min_objective = tree->min_objective();
            /* not_captured = default rule truthtable & ~ captured */
            rule_vandnot(not_captured,
                         tree->rule(tree->root()->id()).truthtable, captured,
                         tree->nsamples(), &cnt);
            evaluate_children<N, Q, P>(tree, node_ordered.first, not_captured,
                                       node_ordered.second, construct_policy, q,
                                       permutation_insert, p);
            logger.addToEvalChildrenTime(time_diff(t1));
            logger.incEvalChildrenNum();

            if (tree->min_objective() < min_objective) {
                min_objective = tree->min_objective();
                printf("before garbage_collect. num_nodes: %zu, log10(remaining): %zu\n", tree->num_nodes(), logger.getLogRemainingSpaceSize());
                logger.dumpState();
                tree->garbage_collect();
                logger.dumpState();
                printf("after garbage_collect. num_nodes: %zu, log10(remaining): %zu\n", tree->num_nodes(), logger.getLogRemainingSpaceSize());
            }
        }
        if (queue_min_length < logger.getQueueMinLen()) {
            // garbage collect the permutation map: can be simplified for the case of BFS
            queue_min_length = logger.getQueueMinLen();
            pmap_garbage_collect(p, queue_min_length);
        }
        ++num_iter;
        if ((num_iter % 10000) == 0) {
            if (p)
                printf("iter: %zu, tree: %zu, queue: %zu, pmap: %zu, log10(remaining): %zu, time elapsed: %f\n",
                       num_iter, tree->num_nodes(), q->size(), p->size(), logger.getLogRemainingSpaceSize(), time_diff(start));
            else
                printf("iter: %zu, tree: %zu, queue: %zu, log10(remaining): %zu, time elapsed: %f\n",
                       num_iter, tree->num_nodes(), q->size(), logger.getLogRemainingSpaceSize(), time_diff(start));
        }
        if ((num_iter % logger.getFrequency()) == 0)
            logger.dumpState();     // want ~1000 records for detailed figures
        if (num_iter == 20000) {
            return num_iter;
        }
    }
    logger.dumpState(); // second last log record (before queue elements deleted)
    if (p)
        printf("iter: %zu, tree: %zu, queue: %zu, pmap: %zu, log10(remaining): %zu, time elapsed: %f\n",
               num_iter, tree->num_nodes(), q->size(), p->size(), logger.getLogRemainingSpaceSize(), time_diff(start));
    else
        printf("iter: %zu, tree: %zu, queue: %zu, log10(remaining): %zu, time elapsed: %f\n",
               num_iter, tree->num_nodes(), q->size(), logger.getLogRemainingSpaceSize(), time_diff(start));

    if (q->empty())
        printf("Exited because queue empty\n");
    else
        printf("Exited because max number of nodes in the tree was reached\n");

    char fname[] = "queue.txt";
    ofstream f;
    printf("Writing queue elements to: %s\n", fname);
    f.open(fname, ios::out | ios::trunc);
    f << "lower_bound length num_captured rule_list\n";

    printf("Deleting queue elements and corresponding nodes in the cache, since they may not be reachable by the tree's destructor\n");
    N* node;
    while (!q->empty()) {
        node = front(q);
        q->pop();
        if (node->deleted()) {
            tree->decrement_num_nodes();
            delete node;
        } else {
            auto pp_pair = node->get_prefix_and_predictions();
            std::vector<size_t> prefix = std::move(pp_pair.first);
            std::vector<bool> predictions = std::move(pp_pair.second);
            f << node->lower_bound() << " " << node->depth() << " "
              << node->num_captured() << " ";
            for(size_t i = 0; i < prefix.size(); ++i) {
                f << tree->rule_features(prefix[i]) << "~"
                  << predictions[i] << ";";
            }
            f << "default~" << predictions.back() << "\n";
        }
    }
    f.close();
    logger.dumpState(); // last log record (before cache deleted)

    rule_vfree(&captured);
    rule_vfree(&not_captured);
    return num_iter;
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
void delete_subtree(CacheTree<N>* tree, N* node, bool destructive, bool update_remaining_state_space) {
    N* child;
    typename std::map<size_t, N*>::iterator iter;
    if (node->done()) {
        iter = node->children_begin();
        while (iter != node->children_end()) {
            child = iter->second;
            delete_subtree<N>(tree, child, destructive, update_remaining_state_space);
            ++iter;
        }
        tree->decrement_num_nodes(); // always delete interior (non-leaf) nodes
        delete node;
    } else {
        if (destructive) {  // only delete leaf nodes in destructive mode
            tree->decrement_num_nodes();
            delete node;
        } else {
            logger.decPrefixLen(node->depth());
            if (update_remaining_state_space)
                logger.removeQueueElement(node->depth(), node->lower_bound());
            node->set_deleted();
        }
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
                                                  NullQueue<BaseNode>* q, 
                                                  permutation_insert_signature<BaseNode, PrefixPermutationMap> permutation_insert,
                                                  PrefixPermutationMap* p);

template void
evaluate_children<BaseNode, BaseQueue, PrefixPermutationMap>(CacheTree<BaseNode>* tree,
                                       BaseNode* parent,
                                       VECTOR parent_not_captured,
                                       std::set<size_t> ordered_parent,
                                       construct_signature<BaseNode> construct_policy,
                                       BaseQueue* q,
                                       permutation_insert_signature<BaseNode, PrefixPermutationMap> permutation_insert,
                                       PrefixPermutationMap* p);

template void
evaluate_children<BaseNode, BaseQueue, CapturedPermutationMap>(CacheTree<BaseNode>* tree,
                                       BaseNode* parent,
                                       VECTOR parent_not_captured,
                                       std::set<size_t> ordered_parent,
                                       construct_signature<BaseNode> construct_policy,
                                       BaseQueue* q,
                                       permutation_insert_signature<BaseNode, CapturedPermutationMap> permutation_insert,
                                       CapturedPermutationMap* p);

template void
evaluate_children<CuriousNode, CuriousQueue, PrefixPermutationMap>(CacheTree<CuriousNode>* tree,
                                             CuriousNode* parent,
                                             VECTOR parent_not_captured,
                                             std::set<size_t> ordered_parent,
                                             construct_signature<CuriousNode> construct_policy,
                                             CuriousQueue* q,
                                             permutation_insert_signature<CuriousNode, PrefixPermutationMap> permutation_insert,
                                             PrefixPermutationMap* p);

template std::pair<BaseNode*, std::set<size_t> >
stochastic_select<BaseNode>(CacheTree<BaseNode>* tree, VECTOR not_captured); 

template void
bbound_stochastic<BaseNode, PrefixPermutationMap>(CacheTree<BaseNode>* tree,
                                                  size_t max_num_nodes,
                                                  construct_signature<BaseNode> construct_policy,
                                                  permutation_insert_signature<BaseNode, PrefixPermutationMap> permutation_insert,
                                                  pmap_garbage_collect_signature<PrefixPermutationMap> pmap_garbage_collect,
                                                  PrefixPermutationMap* p);

template void
bbound_stochastic<BaseNode, CapturedPermutationMap>(CacheTree<BaseNode>* tree,
                                                    size_t max_num_nodes,
                                                    construct_signature<BaseNode> construct_policy,
                                                    permutation_insert_signature<BaseNode, CapturedPermutationMap> permutation_insert,
                                                    pmap_garbage_collect_signature<CapturedPermutationMap> pmap_garbage_collect,
                                                    CapturedPermutationMap* p);

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

template int
bbound_queue<BaseNode, BaseQueue, PrefixPermutationMap>(CacheTree<BaseNode>* tree,
                                  size_t max_num_nodes,
                                  construct_signature<BaseNode> construct_policy,
                                  BaseQueue* q,
                                  BaseNode*(*front)(BaseQueue*),
                                  permutation_insert_signature<BaseNode, PrefixPermutationMap> permutation_insert,
                                  pmap_garbage_collect_signature<PrefixPermutationMap> pmap_garbage_collect,
                                  PrefixPermutationMap* p, size_t num_iter);

template int
bbound_queue<BaseNode, BaseQueue, CapturedPermutationMap>(CacheTree<BaseNode>* tree,
                                  size_t max_num_nodes,
                                  construct_signature<BaseNode> construct_policy,
                                  BaseQueue* q,
                                  BaseNode*(*front)(BaseQueue*),
                                  permutation_insert_signature<BaseNode, CapturedPermutationMap> permutation_insert,
                                  pmap_garbage_collect_signature<CapturedPermutationMap> pmap_garbage_collect,
                                  CapturedPermutationMap* p, size_t num_iter);

template int
bbound_queue<CuriousNode, CuriousQueue, PrefixPermutationMap>(CacheTree<CuriousNode>* tree,
                                        size_t max_num_nodes,
                                        construct_signature<CuriousNode> construct_policy,
                                        CuriousQueue* q,
                                        CuriousNode*(*front)(CuriousQueue*),
                                        permutation_insert_signature<CuriousNode, PrefixPermutationMap> permutation_insert,
                                        pmap_garbage_collect_signature<PrefixPermutationMap> pmap_garbage_collect,
                                        PrefixPermutationMap* p, size_t num_iter);

template int
bbound_queue<CuriousNode, CuriousQueue, CapturedPermutationMap>(CacheTree<CuriousNode>* tree,
                                        size_t max_num_nodes,
                                        construct_signature<CuriousNode> construct_policy,
                                        CuriousQueue* q,
                                        CuriousNode*(*front)(CuriousQueue*),
                                        permutation_insert_signature<CuriousNode, CapturedPermutationMap> permutation_insert,
                                        pmap_garbage_collect_signature<CapturedPermutationMap> pmap_garbage_collect,
                                        CapturedPermutationMap* p, size_t num_iter);

template void
delete_subtree<BaseNode>(CacheTree<BaseNode>* tree, BaseNode* n, bool destructive, bool update_remaining_state_space);

template void
delete_subtree<CuriousNode>(CacheTree<CuriousNode>* tree, CuriousNode* n, bool destructive, bool update_remaining_state_space);
