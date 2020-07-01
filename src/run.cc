#include <stdio.h>
#include <iostream>
#include <set>

#include "queue.hh"
#include "run.hh"
#include "features.hh"
#include "logger.hh"

#define BUFSZ 512

std::mutex log_lk;
std::mutex min_obj_lk;
std::mutex inactive_thread_lk;
std::mutex shared_q_lk;
extern std::atomic<int> emptyQueues;

NullLogger* logger;
FeatureToggle* featureDecisions;

extern "C" {

double run_corels (double c, char* vstring, int curiosity_policy,
                   int map_type, int ablation, int calculate_size, int nrules, int nlabels,
                   int nsamples, rule_t* rules, rule_t* labels, rule_t* meta, int freq, char* log_fname,
                   PermutationMap*& pmap, CacheTree*& tree, Queue*& queue, double& init,
                   int verbosity, int num_threads, int max_num_nodes, int nmeta, int random_seed,
                   std::vector<int>* rulelist, std::vector<int>* classes) {
    printf("log_fname=%s\n", log_fname);
    if (verbosity >= 10)
        print_machine_info();

    if (verbosity >= 1000) {
        printf("\n%d rules %d samples\n\n", nrules, nsamples);
        rule_print_all(rules, nrules, nsamples, 1);

        printf("\nLabels (%d) for %d samples\n\n", nlabels, nsamples);
        rule_print_all(labels, nlabels, nsamples, 1);
    }


    if (verbosity > 1)
        logger = new Logger(c, nrules, verbosity, log_fname, freq);
    else
        logger = new NullLogger();

    featureDecisions = new FeatureToggle(false);

    init = get_timestamp();
    char run_type[BUFSZ];
    strcpy(run_type, "LEARNING RULE LIST via ");
    char const *type = "node";
    std::function<bool(EntryType, EntryType)> cmp;
    if (curiosity_policy == 1) {
        strcat(run_type, "CURIOUS");
        cmp = curious_cmp;
        type = "curious";
    } else if (curiosity_policy == 2) {
        strcat(run_type, "LOWER BOUND");
        cmp = lb_cmp;
    } else if (curiosity_policy == 3) {
        strcat(run_type, "OBJECTIVE");
        cmp = objective_cmp;
    } else if (curiosity_policy == 4) {
        strcat(run_type, "DFS");
        cmp = dfs_cmp;
    } else {
        strcat(run_type, "BFS");
        cmp = base_cmp;
    }

    PermutationMap* p;
    if (map_type == 1) {
        strcat(run_type, " Prefix Map\n");
        PrefixPermutationMap* prefix_pmap = new PrefixPermutationMap;
        p = (PermutationMap*) prefix_pmap;
    } else if (map_type == 2) {
        strcat(run_type, " Captured Symmetry Map\n");
        CapturedPermutationMap* cap_pmap = new CapturedPermutationMap;
        p = (PermutationMap*) cap_pmap;
    } else {
        strcat(run_type, " No Permutation Map\n");
        NullPermutationMap* null_pmap = new NullPermutationMap;
        p = (PermutationMap*) null_pmap;
    }

    tree = new CacheTree(nsamples, nrules, c, num_threads,
            rules, labels, meta, ablation, calculate_size, type, random_seed);
    printf("%s", run_type);

    // Initialize logger
    bbound_init(tree);

    // Set up per-thread queues
    std::thread threads[num_threads];

    Queue* qs[num_threads];
    for(size_t i = 0; i < num_threads; ++i) {
        qs[i] = new Queue(cmp, run_type);
        tracking_vector<unsigned short, DataStruct::Tree> init_rules = tree->get_subrange(i);
        InternalRoot* iroot = new InternalRoot(tree->root(), init_rules);
        qs[i]->push(iroot);
    }

    SharedQueue* shared_q = new SharedQueue();

    // Let the threads loose
    emptyQueues = 0;
    printf("num_threads=%d", num_threads);
    for(size_t i = 0; i < num_threads; ++i) {
        threads[i] = std::thread(bbound, tree, max_num_nodes, qs[i], p, i, shared_q);
    }

    for(size_t i = 0; i < num_threads; ++i) {
        threads[i].join();
    }

    size_t tree_mem = logger->getTreeMemory();
    size_t pmap_mem = logger->getPmapMemory();
    size_t queue_mem = logger->getQueueMemory();
    printf("TREE mem usage: %zu\n", tree_mem);
    printf("PMAP mem usage: %zu\n", pmap_mem);
    printf("QUEUE mem usage: %zu\n", queue_mem);

    printf("final num_nodes: %zu\n", tree->num_nodes());
    printf("final num_evaluated: %zu\n", logger->getTreeNumEvaluated());
    printf("final min_objective: %1.5f\n", tree->min_objective());
    double accuracy = 1 - tree->min_objective() + c*tree->opt_rulelist().size();
    printf("final accuracy: %1.5f\n",
           accuracy);
    print_final_rulelist(tree->opt_rulelist(), tree->opt_predictions(),
                         0, rules, labels, log_fname);

    for(size_t i = 0; i < tree->opt_rulelist().size(); i++) {
        rulelist->push_back(tree->opt_rulelist()[i]);
        classes->push_back(tree->opt_predictions()[i]);
    }
    classes->push_back(tree->opt_predictions().back());

    printf("final total time: %f\n", get_time_diff(init));
    printf("Number of tree acquisitions: %zu\n", tree->n_acc());

    logger->dumpState();
    logger->closeFile();

    printf("delete queue(s)\n");
//    for(size_t i = 0; i < num_threads; ++i) {
//        delete qs[i];
//    }
    printf("delete shared queue\n");
    assert(shared_q->size_approx() == 0);
    delete shared_q;
    printf("delete permutation map\n");
    delete p;
    printf("tree destructors\n");
    // TODO: maybe deadlocking here?
    delete tree;
    delete featureDecisions;
    delete logger;
    if (meta) {
        printf("\ndelete identical points indicator");
        rules_free(meta, nmeta, 0);
    }
    printf("\ndelete rules\n");
    rules_free(rules, nrules, 1);
    printf("delete labels\n");
    rules_free(labels, nlabels, 0);

    tree = nullptr;
    queue = nullptr;
    pmap = nullptr;

    return accuracy;
}

} // end extern C