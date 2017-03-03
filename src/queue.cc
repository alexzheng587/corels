#include "queue.hh"
#include "pmap.hh"
#include <algorithm>

extern int ablation;

/*
 * Constructs a node of type BaseNode.
 */
/*
BaseNode* base_construct_policy(unsigned short new_rule, size_t nrules, bool prediction,
                                bool default_prediction, double lower_bound,
                                double objective, BaseNode* parent,
                                int num_not_captured, int nsamples,
                                int len_prefix, double c, double minority) {
    (void) len_prefix, c;
    size_t num_captured = nsamples - num_not_captured;
    return (new BaseNode(new_rule, nrules, prediction, default_prediction,
                         lower_bound, objective, 0, parent, num_captured, minority));
}*/

/*
 * Constructs a node of type CuriousNode.
 */
/*CuriousNode* curious_construct_policy(unsigned short new_rule, size_t nrules, bool prediction,
                                      bool default_prediction, double lower_bound,
                                      double objective, CuriousNode* parent,
                                      int num_not_captured, int nsamples,
                                      int len_prefix, double c, double minority) {
    size_t num_captured = nsamples - num_not_captured;
    double curiosity = (lower_bound - c * len_prefix + c) * nsamples / (double)(num_captured);
    return (new CuriousNode(new_rule, nrules, prediction, default_prediction,
                            lower_bound, objective, curiosity, parent, num_captured, minority));
}*/

BaseQueue::BaseQueue(std::function<bool(Node*, Node*)> cmp)
    : q_(new std::priority_queue<Node*, std::vector<Node*>, 
            std::function<bool(Node*, Node*)> > (cmp)) {}

CuriousQueue::CuriousQueue()
    : q_(new std::priority_queue<Node*, std::vector<Node*>, 
            std::function<bool(Node*, Node*)> > (curious)) {}

LowerBoundQueue::LowerBoundQueue()
    : q_(new std::priority_queue<Node*, std::vector<Node*>, 
            std::function<bool(Node*, Node*)> > (lb)) {}

ObjectiveQueue::ObjectiveQueue()
    : q_(new std::priority_queue<Node*, std::vector<Node*>, 
            std::function<bool(Node*, Node*)> > (objective)) {}

DFSQueue::DFSQueue()
    : q_(new std::priority_queue<Node*, std::vector<Node*>, 
            std::function<bool(Node*, Node*)> > (dfs)) {}

/*
 * Performs incremental computation on a node, evaluating the bounds and inserting into the cache,
 * queue, and permutation map if appropriate.
 * This is the function that contains the majority of the logic of the algorithm.
 *
 * parent -- the node that is going to have all of its children evaluated.
 * parent_not_captured -- the vector representing data points NOT captured by the parent.
 * construct_policy -- function pointer to how to construct nodes.
 * permutation_insert -- function pointer to how to insert into the permutation map.
 * p -- the permutation map
 */
void evaluate_children(CacheTree* tree, Node* parent, VECTOR parent_not_captured,
                       BaseQueue* q, PermutationMap* p) {
    auto pp_pair = parent->get_prefix_and_predictions();
    std::vector<unsigned short> parent_prefix = std::move(pp_pair.first);
    std::vector<bool> parent_predictions = std::move(pp_pair.second);

    VECTOR captured, captured_zeros, not_captured, not_captured_zeros;
    int num_captured, c0, c1, captured_correct;
    int num_not_captured, d0, d1, default_correct;
    bool prediction, default_prediction;
    double lower_bound, objective, parent_lower_bound, lb;
    double parent_minority;
    double minority = 0.;
    int nsamples = tree->nsamples();
    int nrules = tree->nrules();
    double c = tree->c();
    double threshold = c * nsamples;
    rule_vinit(nsamples, &captured);
    rule_vinit(nsamples, &captured_zeros);
    rule_vinit(nsamples, &not_captured);
    rule_vinit(nsamples, &not_captured_zeros);
    VECTOR not_captured_minority;
    int num_not_captured_minority;
    rule_vinit(nsamples, &not_captured_minority);
    int i, len_prefix;
    len_prefix = parent->depth() + 1;
    parent_lower_bound = parent->lower_bound();
    parent_minority = parent->minority();
    double t0 = timestamp();
    for (i = 1; i < nrules; i++) {
        double t1 = timestamp();
        if (std::find(parent_prefix.begin(), parent_prefix.end(), i) != parent_prefix.end())
            continue;
        // captured represents data captured by the new rule
        rule_vand(captured, parent_not_captured, tree->rule(i).truthtable, nsamples, &num_captured);
        // lower bound on antecedent support
        if ((ablation != 1) && (num_captured < threshold)) 
            continue;
        rule_vand(captured_zeros, captured, tree->label(0).truthtable, nsamples, &c0);
        c1 = num_captured - c0;
        if (c0 > c1) {
            prediction = 0;
            captured_correct = c0;
        } else {
            prediction = 1;
            captured_correct = c1;
        }
        // lower bound on accurate antecedent support
        if ((ablation != 1) && (captured_correct < threshold))
            continue;
        lower_bound = parent_lower_bound - parent_minority + (float)(num_captured - captured_correct) / nsamples + c;
        logger.addToLowerBoundTime(time_diff(t1));
        logger.incLowerBoundNum();
        // hierarchical objective lower bound
        if (lower_bound >= tree->min_objective())
            continue;
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

            logger.setTreeMinObj(objective);
            tree->update_min_objective(objective);
            tree->update_opt_rulelist(parent_prefix, i);
            tree->update_opt_predictions(parent_predictions, prediction, default_prediction);

            // dump state when min objective is updated
            logger.dumpState();
        }
        if (tree->meta_size() == 1) {
            rule_vand(not_captured_minority, not_captured, tree->meta(0).truthtable, nsamples, &num_not_captured_minority);
            minority = (float)(num_not_captured_minority) / nsamples;
            lower_bound += minority;
        }
        if (ablation != 2)
            lb = lower_bound + c;
        else
            lb = lower_bound;
        if (lb < tree->min_objective()) {
            Node* n;
            // check permutation bound
            double t3 = timestamp();
            if (p) {
                n = p->insert(i, nrules, prediction, default_prediction,
                                       lower_bound, objective, parent, num_not_captured, nsamples,
                                       len_prefix, c, minority, tree, not_captured, parent_prefix);
                logger.incPermMapInsertionNum();
            }
            else
                n = tree->construct_node(i, nrules, prediction, default_prediction,
                                    lower_bound, objective, parent,
                                    num_not_captured, nsamples, len_prefix, c, minority);
            logger.addToPermMapInsertionTime(time_diff(t3));
            if (n) {
                double t4 = timestamp();
                tree->insert(n);
                logger.incTreeInsertionNum();
                logger.incPrefixLen(len_prefix);
                logger.addToTreeInsertionTime(time_diff(t4));
                if (q) {
                    double t5 = timestamp();
                    q->push(n);
                    logger.setQueueSize(q->size());
                    logger.addQueueElement(len_prefix, lower_bound, true);
                    logger.addToQueueInsertionTime(time_diff(t5));
                }
            }
        } // else:  objective lower bound with one-step lookahead
    }

    rule_vfree(&captured);
    rule_vfree(&captured_zeros);
    rule_vfree(&not_captured);
    rule_vfree(&not_captured_zeros);
    rule_vfree(&not_captured_minority);

    logger.addToRuleEvalTime(time_diff(t0));
    logger.incRuleEvalNum();
    logger.decPrefixLen(parent->depth());
    logger.removeQueueElement(len_prefix - 1, parent_lower_bound, true);
    if (parent->num_children() == 0) {
        tree->prune_up(parent);
    } else {
        parent->set_done();
        tree->increment_num_evaluated();
    }
}

/*
 * Randomly selects a leave node to evaluate by stochastically walking down the tree.
 */
Node* stochastic_select(CacheTree* tree, VECTOR not_captured) {
    Node* node = tree->root();
    rule_copy(not_captured, tree->rule(node->id()).truthtable, tree->nsamples());
    int cnt;
    while (node->done()) {
        if ((node->lower_bound() + tree->c()) >= tree->min_objective()) {
            if (node->depth() > 0) {
                Node* parent = node->parent();
                parent->delete_child(node->id());
                delete_subtree(tree, node, true, true);
                if (parent->num_children() == 0)
                    tree->prune_up(parent);
            }
            return NULL;
        }
        if (node->num_children() == 0) {
            tree->prune_up(node);
            return NULL;
        }
        node = node->random_child();
        rule_vandnot(not_captured, not_captured, tree->rule(node->id()).truthtable, 
                tree->nsamples(), &cnt);
    }
    logger.setCurrentLowerBound(node->lower_bound() + tree->c());
    return node;
}

/*
 * Uses a stochastic search process to explore the search space.
 */
void bbound_stochastic(CacheTree* tree, size_t max_num_nodes,
                       PermutationMap* p) {
    double start = timestamp();
    logger.setInitialTime(start);
    // initial log record
    logger.dumpState();

    double min_objective = 1.0;
    Node* node_ordered;
    VECTOR not_captured;
    NullQueue* q = NULL;
    logger.incPrefixLen(0);

    rule_vinit(tree->nsamples(), &not_captured);

    size_t num_iter = 0;
    tree->insert_root();
    logger.incTreeInsertionNum();
    // log record for empty rule list
    logger.dumpState();
    while ((tree->num_nodes() < max_num_nodes) and (tree->num_nodes() > 0)) {
        double t0 = timestamp();
        node_ordered = stochastic_select(tree, not_captured);
        logger.addToNodeSelectTime(time_diff(t0));
        logger.incNodeSelectNum();
        if (node_ordered) {
            double t1 = timestamp();
            min_objective = tree->min_objective();
            evaluate_children(tree, node_ordered, not_captured, q, p);
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
                printf("iter: %zu, tree: %zu, pmap: %zu, time elapsed: %f\n",
                       num_iter, tree->num_nodes(), p->size(), time_diff(start));
            else
                printf("iter: %zu, tree: %zu, time elapsed: %f\n",
                       num_iter, tree->num_nodes(), time_diff(start));
        }
        if ((num_iter % logger.getFrequency()) == 0) {
            // want ~1000 records for detailed figures
            logger.dumpState();     
        }
    }
    logger.dumpState();
    rule_vfree(&not_captured);
}

/*
 * Selects the node off the front of the queue.
 *
 * front -- determines whether the queue is FIFO or a priority queue.
 */
Node* queue_select(CacheTree* tree, BaseQueue* q, VECTOR captured) {
    int cnt;

    Node* selected_node = q->front();
    q->pop();
    logger.setCurrentLowerBound(selected_node->lower_bound() + tree->c());

    Node* node = selected_node;
    // delete leaf nodes that were lazily marked
    if (node->deleted()) {
        tree->decrement_num_nodes();
        delete node;
        return NULL;
    }

    rule_vclear(tree->nsamples(), captured);

    while (node != tree->root()) {
        if (node->deleted()) {
            delete node;
            return NULL;
        }
        rule_vor(captured,
                 captured, tree->rule(node->id()).truthtable,
                 tree->nsamples(), &cnt);
        node = node->parent();
    }
    return selected_node;
}

/*
 * Explores the search space by using a queue to order the search process.
 * The queue can be ordered by DFS, BFS, or an alternative priority metric (e.g. lower bound).
 */
int bbound_queue(CacheTree* tree, size_t max_num_nodes, BaseQueue* q,
                 PermutationMap* p, size_t num_iter, size_t switch_iter) {
    bool print_queue = 0;
    double start;
    int cnt;
    double min_objective;
    Node* node_ordered;
    VECTOR captured, not_captured;
    rule_vinit(tree->nsamples(), &captured);
    rule_vinit(tree->nsamples(), &not_captured);

    size_t queue_min_length = logger.getQueueMinLen();

    // Set up is different when we initialize the tree as opposed to when we restart execution.
    if (tree->num_nodes() == 0) {
        start = timestamp();
        logger.setInitialTime(start);
        logger.initializeState();
        // initial log record
        logger.dumpState();         

        min_objective = 1.0;
        tree->insert_root();
        logger.incTreeInsertionNum();
        q->push(tree->root());
        logger.setQueueSize(q->size());
        logger.incPrefixLen(0);
        // log record for empty rule list
        logger.dumpState(); 
        logger.initRemainingSpaceSize();
    } else {
        start = logger.getInitialTime();
        min_objective = tree->min_objective();
    }
    while ((tree->num_nodes() < max_num_nodes) &&
           !q->empty()) {
        double t0 = timestamp();
        node_ordered = queue_select(tree, q, captured);
        logger.addToNodeSelectTime(time_diff(t0));
        logger.incNodeSelectNum();
        if (node_ordered) {
            double t1 = timestamp();
            min_objective = tree->min_objective();
            // not_captured = default rule truthtable & ~ captured
            rule_vandnot(not_captured,
                         tree->rule(tree->root()->id()).truthtable, captured,
                         tree->nsamples(), &cnt);
            evaluate_children(tree, node_ordered, not_captured, q, p);
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
        if (q)
            logger.setQueueSize(q->size());
        if (queue_min_length < logger.getQueueMinLen()) {
            // garbage collect the permutation map: can be simplified for the case of BFS
            queue_min_length = logger.getQueueMinLen();
            //pmap_garbage_collect(p, queue_min_length);
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
        if ((num_iter % logger.getFrequency()) == 0) {
            // want ~1000 records for detailed figures
            logger.dumpState();
        }
        if (num_iter == switch_iter)
            return num_iter;
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

    printf("TREE mem usage: %zu\n", logger.getTreeMemory());
    printf("PMAP mem usage: %zu\n", logger.getPmapMemory());

    // Print out queue
    ofstream f;
    if (print_queue) {
        char fname[] = "queue.txt";
        printf("Writing queue elements to: %s\n", fname);
        f.open(fname, ios::out | ios::trunc);
        f << "lower_bound objective length frac_captured rule_list\n";
    }

    // Clean up data structures
    printf("Deleting queue elements and corresponding nodes in the cache,"
            "since they may not be reachable by the tree's destructor\n");
    printf("\nminimum objective: %1.10f\n", tree->min_objective());
    Node* node;
    double min_lower_bound = 1.0;
    double lb;
    size_t num = 0;
    while (!q->empty()) {
        node = q->front();
        q->pop();
        if (node->deleted()) {
            tree->decrement_num_nodes();
            delete node;
        } else {
            lb = node->lower_bound() + tree->c();
            if (lb < min_lower_bound)
                min_lower_bound = lb;
            if (print_queue) {
                auto pp_pair = node->get_prefix_and_predictions();
                std::vector<unsigned short> prefix = std::move(pp_pair.first);
                std::vector<bool> predictions = std::move(pp_pair.second);
                f << node->lower_bound() << " " << node->objective() << " " << node->depth() << " "
                  << (double) node->num_captured() / (double) tree->nsamples() << " ";
                for(size_t i = 0; i < prefix.size(); ++i) {
                    f << tree->rule_features(prefix[i]) << "~"
                      << predictions[i] << ";";
                }
                f << "default~" << predictions.back() << "\n";
                num++;
            }
        }
    }
    printf("minimum lower bound in queue: %1.10f\n\n", min_lower_bound);
    if (print_queue)
        f.close();
    // last log record (before cache deleted)
    logger.dumpState();

    rule_vfree(&captured);
    rule_vfree(&not_captured);
    return num_iter;
}
