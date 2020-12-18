#ifdef VAL
#include "common.hh"
#endif
#include "loss.hh"
#include "queue.hh"
#include <algorithm>
#include <iostream>
#include <stdio.h>

Queue::Queue(std::function<bool(Node*, Node*)> cmp, char const* type)
        : q_(new q(cmp)), type_(type) {}

/*
 * Performs incremental computation on a node, evaluating the bounds and inserting into the cache,
 * queue, and permutation map if appropriate.
 * This is the function that contains the majority of the logic of the algorithm.
 *
 * parent -- the node that is going to have all of its children evaluated.
 * parent_not_captured -- the vector representing data points NOT captured by the parent.
 */
void evaluate_children(CacheTree* tree, Node* parent, tracking_vector<unsigned short, DataStruct::Tree> parent_prefix,
                       VECTOR parent_not_captured, Queue* q, PermutationMap* p) {
    int nsamples = tree->nsamples();
    int nrules = tree->nrules();
    size_t len_prefix = parent->depth() + 1;
    //std::set<std::string> verbosity = logger->getVerbosity();

    Loss* loss = tree->loss();
    ProblemState* state = tree->threadState(0); // TODO update when multithreaded
    loss->initialize_subproblem(state, tree, parent, parent_not_captured);

    double t0 = timestamp();
    // nrules is actually the number of rules + 1 (since it includes the default rule), so the maximum
    // value of i is nrules - 1 instead of nrules
    // the loop bound represents the trivial upper bound on prefix length
    // no other upper bound on prefix length is implemented
    for (int i = 1; i < nrules; i++) {
        double t1 = timestamp();
        // check if this rule is already in the prefix
        if (std::find(parent_prefix.begin(), parent_prefix.end(), i) != parent_prefix.end())
            continue;

        // Check support bounds, compute minimum loss and optimal prediction
        rule_vand(state->captured, parent_not_captured, tree->rule(i).truthtable, tree->nsamples(), &state->n_captured);
        if (tree->support_bound_enabled() && loss->minimum_support_bounded(state, parent, tree))
            continue;

        // Evaluate the two possible predictions and set state variables based on whichever minimizes loss
        rule_vand(state->captured_zeros, state->captured, tree->label(0).truthtable, tree->nsamples(),
                &state->n_zeros_captured);
        // variable initialization for default error calculation
        rule_vandnot(state->not_captured, parent_not_captured, state->captured, nsamples, &state->n_not_captured);
        loss->make_prediction(state, tree->npos(), tree->nneg(), parent);

        if (tree->support_bound_enabled() && loss->accurate_support_bounded(state))
            continue;

        loss->compute_lower_bound(state, tree->npos(), tree->nneg());
        logger->addToLowerBoundTime(time_diff(t1));
        logger->incLowerBoundNum();

        if (loss->min_objective_bounded(state, tree)) // hierarchical objective lower bound
            continue;

        // if (state->rule_lower_bound >= tree->min_objective())
        //     continue;

        double t2 = timestamp();

        loss->make_default_prediction(state, tree->npos(), tree->nneg());
        loss->compute_objective(state);

        logger->addToObjTime(time_diff(t2));
        logger->incObjNum();

        // equivalent points bound
        if (tree->has_minority() && loss->equivalent_points_bounded(state, tree, parent))
            continue;

        // dump state when min objective is updated
        if (state->objective < tree->min_objective()) {
            if (logger->getVerbosity().count("progress")) {
                printf("min(objective): %1.5f -> %1.5f, length: %zu, cache size: %zu\n",
                       tree->min_objective(), state->objective, len_prefix, tree->num_nodes());
            }

            logger->setTreeMinObj(state->objective);
            tree->update_opt_rulelist(state->objective, parent_prefix, i, parent, state->rule_prediction,
                                      state->default_prediction);

            logger->dumpState();
        }

        // only add node to our data structures if its children will be viable
        if (tree->lookahead_bound_enabled() && Loss::lookahead_bounded(state, tree))
            continue;

        double t3 = timestamp();
        // check permutation bound
        if (loss->type() == AUC) {
            AUCLoss *auc_loss = (AUCLoss*) loss;
            auc_loss->compute_new_classes(state, tree, parent, tree->rule(i).truthtable);
        }
        Node* n = p->insert(i, state, parent, tree, parent_prefix);
        logger->addToPermMapInsertionTime(time_diff(t3));
        // n is NULL if this rule fails the permutation bound
        if (!n)
            continue;

        double t4 = timestamp();
        tree->insert(n);
        logger->incTreeInsertionNum();
        logger->incPrefixLen(len_prefix);
        logger->addToTreeInsertionTime(time_diff(t4));

        double t5 = timestamp();
        n->set_in_queue(true);
        q->push(n);
        logger->setQueueSize(q->size());
        if (tree->calculate_size())
            logger->addQueueElement(len_prefix, state->rule_lower_bound, false);

        logger->addToQueueInsertionTime(time_diff(t5));
    }

    logger->addToRuleEvalTime(time_diff(t0));
    logger->incRuleEvalNum();
    logger->decPrefixLen(parent->depth());
    if (tree->calculate_size())
        logger->removeQueueElement(len_prefix - 1, state->parent_lower_bound, false);

    if (parent->num_children() != 0) {
        parent->set_done();
        tree->increment_num_evaluated();
    }
}

/*
 * Explores the search space by using a queue to order the search process.
 * The queue can be ordered by DFS, BFS, or an alternative priority metric (e.g. lower bound).
 */
int bbound(CacheTree* tree, size_t max_num_nodes, Queue* q, PermutationMap* p, double min_objective) {
    size_t num_iter = 0;
    int cnt;
    VECTOR captured, not_captured;
    rule_vinit(tree->nsamples(), &captured);
    rule_vinit(tree->nsamples(), &not_captured);

    size_t queue_min_length = logger->getQueueMinLen();

    double start = timestamp();
    logger->setInitialTime(start);
    logger->initializeState(tree->calculate_size());
    std::set<std::string> verbosity = logger->getVerbosity();
    // initial log record
    logger->dumpState();

    tree->insert_root();
    if (tree->min_objective() > min_objective) {
        tracking_vector<unsigned short, DataStruct::Tree> opt_rl;
        tree->update_opt_rulelist(min_objective, opt_rl, 0, tree->root(), false, false);
    }
    logger->incTreeInsertionNum();
    q->push(tree->root());
    logger->setQueueSize(q->size());
    logger->incPrefixLen(0);
    // log record for empty rule list
    logger->dumpState();
    while ((tree->num_nodes() < max_num_nodes) && !q->empty()) {
        double t0 = timestamp();
        std::pair<Node*, tracking_vector<unsigned short, DataStruct::Tree> > node_ordered = q->select(tree, captured);
        logger->addToNodeSelectTime(time_diff(t0));
        logger->incNodeSelectNum();
        if (node_ordered.first) {
            double t1 = timestamp();
            // not_captured = default rule truthtable & ~ captured
            rule_vandnot(not_captured, tree->rule(0).truthtable, captured, tree->nsamples(), &cnt);
            evaluate_children(tree, node_ordered.first, node_ordered.second, not_captured, q, p);
            logger->addToEvalChildrenTime(time_diff(t1));
            logger->incEvalChildrenNum();

            if (tree->min_objective() < min_objective) {
                min_objective = tree->min_objective();
                if (verbosity.count("progress"))
                    printf("before garbage_collect. num_nodes: %zu, log10(remaining): %zu\n",
                           tree->num_nodes(), logger->getLogRemainingSpaceSize());
                logger->dumpState();
                tree->garbage_collect();
                logger->dumpState();
                if (verbosity.count("progress"))
                    printf("after garbage_collect. num_nodes: %zu, log10(remaining): %zu\n", tree->num_nodes(),
                           logger->getLogRemainingSpaceSize());
            }
        }
        logger->setQueueSize(q->size());
        if (queue_min_length < logger->getQueueMinLen()) {
            // garbage collect the permutation map: can be simplified for the case of BFS
            queue_min_length = logger->getQueueMinLen();
            // pmap_garbage_collect(p, queue_min_length);
        }
        ++num_iter;
        if ((num_iter % 10000) == 0) {
            if (verbosity.count("progress"))
                printf("iter: %zu, tree: %zu, queue: %zu, pmap: %zu, log10(remaining): %zu, time elapsed: %f\n",
                       num_iter, tree->num_nodes(), q->size(), p->size(), logger->getLogRemainingSpaceSize(),
                       time_diff(start));
        }
        if ((num_iter % logger->getFrequency()) == 0) {
            // want ~1000 records for detailed figures
            logger->dumpState();
        }
    }
    logger->dumpState(); // second last log record (before queue elements deleted)
    if (verbosity.count("progress")) {
        printf("iter: %zu, tree: %zu, queue: %zu, pmap: %zu, log10(remaining): %zu, time elapsed: %f\n",
               num_iter, tree->num_nodes(), q->size(), p->size(), logger->getLogRemainingSpaceSize(), time_diff(start));
        if (q->empty())
            printf("Exited because queue empty\n");
        else
            printf("Exited because max number of nodes in the tree was reached\n");
    }

    size_t tree_mem = logger->getTreeMemory();
    size_t pmap_mem = logger->getPmapMemory();
    size_t queue_mem = logger->getQueueMemory();
    if (verbosity.count("progress")) {
        printf("TREE mem usage: %zu\n", tree_mem);
        printf("PMAP mem usage: %zu\n", pmap_mem);
        printf("QUEUE mem usage: %zu\n", queue_mem);
    }

    // Print out queue
    ofstream f;
    if (verbosity.count("log")) {
        char fname[] = "queue.txt";
        if (verbosity.count("progress"))
            printf("Writing queue elements to: %s\n", fname);
        f.open(fname, ios::out | ios::trunc);
        f << "lower_bound objective length frac_captured rule_list\n";
    }

    // Clean up data structures
    if (verbosity.count("progress")) {
        printf("Deleting queue elements and corresponding nodes in the cache,"
               "since they may not be reachable by the tree's destructor\n");
        printf("\nminimum objective: %1.10f\n", tree->min_objective());
    }
    Node* node;
    double min_lower_bound = 1.0;
    double lb;
    size_t num = 0;
    while (!q->empty()) {
        node = q->front();
        q->pop();
        if (node->deleted()) {
            tree->decrement_num_nodes();
            logger->removeFromMemory(sizeof(*node), DataStruct::Tree);
            delete node;
        } else {
            lb = node->lower_bound() + tree->c();
            if (lb < min_lower_bound)
                min_lower_bound = lb;
            if (verbosity.count("log")) {
                std::pair<tracking_vector<unsigned short, DataStruct::Tree>, tracking_vector<bool, DataStruct::Tree> > pp_pair = node->get_prefix_and_predictions();
                tracking_vector<unsigned short, DataStruct::Tree> prefix = std::move(pp_pair.first);
                tracking_vector<bool, DataStruct::Tree> predictions = std::move(pp_pair.second);
                f << node->lower_bound() << " " << node->objective() << " " << node->depth() << " "
                  << (double) node->num_captured() / (double) tree->nsamples() << " ";
                for (size_t i = 0; i < prefix.size(); ++i) {
                    f << tree->rule_features(prefix[i]) << "~"
                      << predictions[i] << ";";
                }
                f << "default~" << predictions.back() << "\n";
                num++;
            }
        }
    }
    if (verbosity.count("progress"))
        printf("minimum lower bound in queue: %1.10f\n\n", min_lower_bound);
    if (verbosity.count("log"))
        f.close();
    // last log record (before cache deleted)
    logger->dumpState();

    rule_vfree(&captured);
    rule_vfree(&not_captured);
    return num_iter;
}
