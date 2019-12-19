#include "queue.hh"
#include <algorithm>
#include <iostream>
#include <mutex>
#include <numeric>
#include <stdio.h>
#include <sys/resource.h>

extern std::mutex min_obj_lk;
extern std::mutex inactive_thread_lk;
std::condition_variable inactive_thread_cv;

Queue::Queue(std::function<bool(EntryType, EntryType)> cmp, char const *type)
    : q_(new q(cmp)), type_(type), cmp_(cmp) {}

/*
 * Performs incremental computation on a node, evaluating the bounds and inserting into the cache,
 * queue, and permutation map if appropriate.
 * This is the function that contains the majority of the logic of the algorithm.
 *
 * parent -- the node that is going to have all of its children evaluated.
 * parent_not_captured -- the vector representing data points NOT captured by the parent.
 */
void evaluate_children(CacheTree *tree, Node *parent,
                       tracking_vector<unsigned short, DataStruct::Tree> parent_prefix,
                       VECTOR parent_not_captured, tracking_vector<unsigned short, DataStruct::Tree> rules,
                       Queue *q, PermutationMap *p, unsigned short thread_id) {

    VECTOR captured, captured_zeros, not_captured, not_captured_zeros, not_captured_equivalent;
    int num_captured, c0, c1, captured_correct;
    int num_not_captured, d0, d1, default_correct, num_not_captured_equivalent;
    bool prediction, default_prediction;
    double lower_bound, objective, parent_lower_bound, lookahead_bound;
    double parent_equivalent_minority;
    double equivalent_minority = 0.;
    int nsamples = tree->nsamples();
    int nrules = tree->nrules();
    double c = tree->c();
    double threshold = c * nsamples;
    rule_vinit(nsamples, &captured);
    rule_vinit(nsamples, &captured_zeros);
    rule_vinit(nsamples, &not_captured);
    rule_vinit(nsamples, &not_captured_zeros);
    rule_vinit(nsamples, &not_captured_equivalent);
    int i, len_prefix;
    len_prefix = parent->depth() + 1;
    parent_lower_bound = parent->lower_bound();
    parent_equivalent_minority = parent->equivalent_minority();
    double t0 = logger->timestamp();

    for (tracking_vector<unsigned short, DataStruct::Tree>::iterator it = rules.begin();
         it != rules.end(); ++it) {

        i = *it;
        double t1 = logger->timestamp();

        // check if this rule is already in the prefix
        if (std::find(parent_prefix.begin(), parent_prefix.end(), i) !=
            parent_prefix.end())
            continue;
        // captured represents data captured by the new rule
        rule_vand(captured, parent_not_captured, tree->rule(i).truthtable, nsamples, &num_captured);
        // lower bound on antecedent support
        if ((tree->ablation() != 1) && (num_captured < threshold))
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
        if ((tree->ablation() != 1) && (captured_correct < threshold))
            continue;
        // subtract off parent equivalent points bound because we want to use pure lower bound from parent
        lower_bound = parent_lower_bound - parent_equivalent_minority + (double)(num_captured - captured_correct) / nsamples + c;
        logger->addToLowerBoundTime(logger->time_diff(t1));
        logger->incLowerBoundNum();
        //min_obj_lk.lock();
        if (lower_bound >= tree->min_objective()) // hierarchical objective lower bound
                                                  //min_obj_lk.unlock();
            continue;
        double t2 = logger->timestamp();
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
        objective = lower_bound + (double)(num_not_captured - default_correct) / nsamples;
        logger->addToObjTime(logger->time_diff(t2));
        logger->incObjNum();
        if (objective < tree->min_objective()) {
            double old_min_obj = tree->min_objective();
            if (tree->update_obj_and_list(objective, parent_prefix, i, parent, prediction, default_prediction)) {
                printf("THREAD %zu: min(objective): %.17g -> %.17g, length: %d, cache size: %zu\n",
                       thread_id, old_min_obj, objective, len_prefix, tree->num_nodes());
                logger->setTreeMinObj(objective);
                // dump state when min objective is updated
                logger->dumpState();
            }
        }
        // calculate equivalent points bound to capture the fact that the minority points can never be captured correctly
        if (tree->has_minority()) {
            rule_vand(not_captured_equivalent, not_captured, tree->minority(0).truthtable, nsamples, &num_not_captured_equivalent);
            equivalent_minority = (double)(num_not_captured_equivalent) / nsamples;
            lower_bound += equivalent_minority;
        }
        if (tree->ablation() != 2)
            lookahead_bound = lower_bound + c;
        else
            lookahead_bound = lower_bound;
        //min_obj_lk.lock();
        // only add node to our datastructures if its children will be viable
        if (lookahead_bound < tree->min_objective()) {
            //min_obj_lk.unlock()
            double t3 = logger->timestamp();
            // check permutation bound
            Node *n = p->insert(i, nrules, prediction, default_prediction,
                                lower_bound, objective, parent, num_not_captured, nsamples,
                                len_prefix, c, equivalent_minority, tree, not_captured, parent_prefix, thread_id);
            logger->addToPermMapInsertionTime(logger->time_diff(t3));
            // n is NULL if this rule fails the permutaiton bound
            if (n) {
                double t4 = logger->timestamp();
                tree->insert(n, thread_id);
                logger->incTreeInsertionNum();
                logger->incPrefixLen(len_prefix);
                logger->addToTreeInsertionTime(logger->time_diff(t4));
                double t5 = logger->timestamp();
                n->set_in_queue(true);
                InternalRoot* iroot = new InternalRoot(n, NULL);
                q->push(iroot);
                logger->setQueueSize(q->size());
                if (tree->calculate_size())
                    logger->addQueueElement(len_prefix, lower_bound, false);
                logger->addToQueueInsertionTime(logger->time_diff(t5));
            }
        } // else:  objective lower bound with one-step lookahead
        // else { min_obj_lk.unlock() }
    }

    rule_vfree(&captured);
    rule_vfree(&captured_zeros);
    rule_vfree(&not_captured);
    rule_vfree(&not_captured_zeros);
    rule_vfree(&not_captured_equivalent);

    logger->addToRuleEvalTime(logger->time_diff(t0));
    logger->incRuleEvalNum();
    logger->decPrefixLen(parent->depth());
    if (tree->calculate_size())
        logger->removeQueueElement(len_prefix - 1, parent_lower_bound, false);
    if (parent->num_children() == 0) {
        //tree->prune_up(parent);
        //parent->set_done();
    } else {
        //parent->set_done();
        tree->increment_num_evaluated();
    }
}

void bbound_init(CacheTree *tree) {

    // Initialize tree and queue
    if (tree->root() == NULL)
        tree->insert_root();
    logger->initializeState(tree->calculate_size());
    logger->incPrefixLen(0);
    logger->incTreeInsertionNum();
    logger->setQueueSize(0);

    // Initialize log
    logger->dumpState();
    logger->initRemainingSpaceSize();
    logger->setInitialTime(logger->timestamp());
}

// Assumes we have thread_inactive_lk when enters, releases lock within this function
void split_work(CacheTree *tree, Queue* q, SharedQueue *shared_q) {
    std::vector<Queue*> other_qs;
    //shared_q->lock();
    /*if (shared_q->empty()) {
       shared_q->unlock(); 
       inactive_thread_lk.unlock();
       return;
    } else {
        */
    for (size_t i = 0; i < shared_q->size(); ++i) {
        other_qs.push_back(shared_q->pop());
    }
   // }
    shared_q->unlock();
    
    size_t q_sz = (q->size() / (other_qs.size() + 1)) + 1;
    for (std::vector<Queue*>::iterator it = other_qs.begin(); it != other_qs.end(); ++it) {
        Queue* other_q = *it;
        /*for(size_t j = 0; j < q_sz / n_splits; ++j) {
            other_q->push(q->front());
        }*/
        q->move(other_q, q_sz);
        //tree->decrement_num_inactive_threads();
    }
    //inactive_thread_lk.unlock();
    inactive_thread_cv.notify_all();
    //tree->wake_n_inactive(tree->num_inactive_threads());
    //for(unsigned short i = 0; i < tree->num_inactive_threads(); ++i)
    //    inactive_thread_cv.notify_one();
    return;
}

/**
 * Returns true if exiting because queue is empty, false if num_nodes reached
 */
bool bbound_loop(CacheTree *tree, size_t max_num_nodes, Queue *q, PermutationMap *p,
                 VECTOR captured, VECTOR not_captured, unsigned short thread_id, SharedQueue *shared_q, double start) {
    size_t num_iter = 0;
    // Not used anywhere
    int cnt;
    double cur_min_objective = tree->min_objective();

    while ((tree->num_nodes() < max_num_nodes) && !q->empty()) {
        //if (tree->num_inactive_threads() > 0) {
        if (!shared_q->empty()) {
            //inactive_thread_lk.lock();
            shared_q->lock();
            //if (tree->num_inactive_threads() > 0) {
            if (!shared_q->empty() && q->size() > tree->num_threads()) {
                split_work(tree, q, shared_q);
            } else {
                shared_q->unlock();
//                inactive_thread_lk.unlock();
            }
        }
        double t0 = logger->timestamp();
        std::pair<EntryType, tracking_vector<unsigned short, DataStruct::Tree> > node_ordered = q->select(tree, captured, thread_id);
        logger->addToNodeSelectTime(logger->time_diff(t0));
        logger->incNodeSelectNum();
        InternalRoot* iroot = node_ordered.first;
        if(iroot) {
            tracking_vector<unsigned short, DataStruct::Tree> parent_prefix = node_ordered.second;
            Node *current_node = iroot->node();
            tracking_vector<unsigned short, DataStruct::Tree>* initialization_rules = iroot->rules();
            // not_captured = default rule truthtable & ~ captured
            rule_vandnot(not_captured,
                         tree->rule(0).truthtable, captured,
                         tree->nsamples(), &cnt);
            if (initialization_rules == NULL || initialization_rules->empty()) {
                initialization_rules = tree->rule_perm();
            }

            double t1 = logger->timestamp();
            evaluate_children(tree, current_node, parent_prefix, not_captured, *initialization_rules, q, p, thread_id);
            logger->addToEvalChildrenTime(logger->time_diff(t1));
            logger->incEvalChildrenNum();

            min_obj_lk.lock();
            // SET CUR MIN OBJECTIVE
            if (tree->min_objective() < cur_min_objective) {
                cur_min_objective = tree->min_objective();
                min_obj_lk.unlock();
                if (featureDecisions->do_garbage_collection()) {
                    printf("THREAD %zu: before garbage_collect. num_nodes: %zu, log10(remaining): %zu\n",
                           thread_id, tree->num_nodes(), logger->getLogRemainingSpaceSize());
                    logger->dumpState();
                    tree->garbage_collect(thread_id);
                    logger->dumpState();
                    printf("THREAD %zu: after garbage_collect. num_nodes: %zu, log10(remaining): %zu\n", thread_id, tree->num_nodes(), logger->getLogRemainingSpaceSize());
                }
            } else {
                min_obj_lk.unlock();
            }
        }
        logger->setQueueSize(q->size());
        ++num_iter;
        if ((num_iter % 10000) == 0) {
            if (logger->getVerbosity() >= 10)
                printf("THREAD %zu: iter: %zu, tree: %zu, queue: %zu, pmap: %zu, log10(remaining): %zu, time elapsed: %f\n",
                       thread_id, num_iter, tree->num_nodes(), q->size(), p->size(), logger->getLogRemainingSpaceSize(), logger->time_diff(start));
        }
        // want ~1000 records for detailed figures
        if ((num_iter % logger->getFrequency()) == 0)
            logger->dumpState();
    }
    logger->dumpState(); // second last log record (before queue elements deleted)
    if (logger->getVerbosity() >= 10)
        printf("THREAD %zu: iter: %zu, tree: %zu, queue: %zu, pmap: %zu, log10(remaining): %zu\n",
               thread_id, num_iter, tree->num_nodes(), q->size(), p->size(), logger->getLogRemainingSpaceSize());
    if (q->empty()) {
        // printf("Exited because queue empty\n");
        return true;
    } else {
        printf("Exited because max number of nodes (%zu) in the tree was reached\n", tree->num_nodes());
        return false;
    }
}

/*
 * Explores the search space by using a queue to order the search process.
 * The queue can be ordered by DFS, BFS, or an alternative priority metric (e.g. lower bound).
 */
int bbound(CacheTree *tree, size_t max_num_nodes, Queue *q, PermutationMap *p,
           unsigned short thread_id, SharedQueue *shared_q) {
    int cnt;
    // These are initialized in q->select
    VECTOR captured, not_captured;
    rule_vinit(tree->nsamples(), &captured);
    rule_vinit(tree->nsamples(), &not_captured);
    rule_vandnot(not_captured,
                 tree->rule(0).truthtable, captured,
                 tree->nsamples(), &cnt);

    double start = logger->timestamp();
    while (!tree->done()) {
        bool queue_empty = bbound_loop(tree, max_num_nodes, q, p, captured, not_captured, thread_id, shared_q, start);
        // printf("THREAD %zu: time elapsed %f\n", thread_id, logger->time_diff(start));
        if (queue_empty) {
            shared_q->lock();
            shared_q->push(q);
            shared_q->unlock();
            // std::mutex inactive_thread_lk = tree->get_inactive_thread_lk();
            std::unique_lock<std::mutex> unique_inactive_thread_lk(inactive_thread_lk);
            //if (tree->num_inactive_threads() == tree->num_threads() - 1) {
            if (shared_q->size() == tree->num_threads()) {
                //tree->increment_num_inactive_threads();
                printf("All threads done\n");
                tree->set_done(true);
                // tree->unlock_inactive_thread_lk();
                // tree->wake_all_inactive();
                inactive_thread_cv.notify_all();
                printf("Thread %zu notifying and exiting\n", thread_id);
                return 0;
            }
            while (q->empty()) {
                //tree->increment_num_inactive_threads();
                // printf("Thread %zu sleeping\n", thread_id);
                // printf("Thread %zu awake\n", thread_id);
                if (tree->done()) {
                    printf("Thread %zu exiting\n", thread_id);
                    //tree->decrement_num_inactive_threads();
                    return 0;
                } else {
                    inactive_thread_cv.wait(unique_inactive_thread_lk);
                }
                //tree->decrement_num_inactive_threads();
            }
            //tree->unlock_inactive_thread_lk();
        } else {
            // max num nodes has been reached
            tree->set_done(true);
            inactive_thread_cv.notify_all();
        }
    }
    // Print out queue
    /*ofstream f;
    if (print_queue) {
        char fname[] = "queue.txt";
        printf("Writing queue elements to: %s\n", fname);
        f.open(fname, ios::out | ios::trunc);
        f << "lower_bound objective length frac_captured rule_list\n";
        f.close();
    }*/
    // last log record (before cache deleted)
    logger->dumpState();

    rule_vfree(&captured);
    rule_vfree(&not_captured);
    return 0;
}
