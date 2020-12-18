#include <stdio.h>
#include <iostream>
#include <set>
#include <string.h>

#include "loss.hh"
#include "queue.hh"
#include "run.hh"

#define BUFSZ 512

NullLogger* logger;

extern "C" {

double run_corels (double c, char* vstring, char* loss_type_str, int curiosity_policy,
                   int map_type, int ablation, int calculate_size, int latex_out, int nrules, int nlabels,
                   int nsamples, rule_t* rules, rule_t* labels, rule_t* meta, minority_class_t* minority_class, int freq, char* log_fname,
                   PermutationMap*& pmap, CacheTree*& tree, Queue*& queue, double& init,
                   int verbosity, int num_threads, int max_num_nodes, int nmeta, int random_seed,
                   std::vector<int>* rulelist, std::vector<int>* classes, double weight) {

    std::set<std::string> verbosity_set;
    char opt_fname[BUFSZ];
    const char *voptions = "rule|rulelist|label|samples|progress|log|silent";

    char *vopt = NULL;
    char *vcopy = strdup(vstring);
    while ((vopt = strsep(&vcopy, ",")) != NULL) {
        if (!strstr(voptions, vopt)) {
            fprintf(stderr, "verbosity options must be one or more of (%s), separated with commas (i.e. -v progress,log)\n", voptions);
            return -1.0;
        }
        verbosity_set.insert(vopt);
    }
    free(vcopy);

    if (verbosity_set.count("samples") && !(verbosity_set.count("rule") || verbosity_set.count("label"))) {
        fprintf(stderr, "verbosity 'samples' option must be combined with at least one of (rule|label)\n");
        return -1.0;
    }
    if (verbosity_set.size() > 2 && verbosity_set.count("silent")) {
        fprintf(stderr, "verbosity 'silent' option must be passed without any additional verbosity parameters\n");
        return -1.0;
    }

    if (verbosity_set.size() == 0) {
        verbosity_set.insert("progress");
    }

    if (verbosity_set.count("silent")) {
        verbosity_set.clear();
    }

    if (verbosity_set.count("log"))
        print_machine_info();

    if (verbosity_set.count("rule")) {
        printf("%d rules %d samples\n\n", nrules, nsamples);
        rule_print_all(rules, nrules, nsamples, (verbosity_set.count("samples")));
        printf("\n\n");
    }

    if (verbosity_set.count("label")) {
        printf("Labels (%d) for %d samples\n\n", nlabels, nsamples);
        rule_print_all(labels, nlabels, nsamples, (verbosity_set.count("samples")));
        printf("\n\n");
    }

    if (verbosity_set.count("log")) {
        logger = new Logger(c, nrules, verbosity_set, log_fname, freq);
    } else {
        logger = new NullLogger();
        logger->setVerbosity(verbosity_set);
    }

    init = timestamp();
    char run_type[BUFSZ];
    Queue* q;
    strcpy(run_type, "LEARNING RULE LIST via ");
    char const *type = "node";
    if (curiosity_policy == 1) {
        strcat(run_type, "CURIOUS");
        q = new Queue(curious_cmp, run_type);
        type = "curious";
    } else if (curiosity_policy == 2) {
        strcat(run_type, "LOWER BOUND");
        q = new Queue(lb_cmp, run_type);
    } else if (curiosity_policy == 3) {
        strcat(run_type, "OBJECTIVE");
        q = new Queue(objective_cmp, run_type);
    } else if (curiosity_policy == 4) {
        strcat(run_type, "DFS");
        q = new Queue(dfs_cmp, run_type);
    } else {
        strcat(run_type, "BFS");
        q = new Queue(base_cmp, run_type);
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

    // this is a provisional solution using fixed loss type
    strcat(run_type, " and loss function type ");
    strcat(run_type, loss_type_str);
    strcat(run_type, "\n");

    Loss* loss;
    if (strcmp(loss_type_str, "bacc") == 0) {
        loss = new BalancedAccuracy();
    } else if (strcmp(loss_type_str, "wacc") == 0) {
        loss = new WeightedAccuracy(weight);
    } else if (strcmp(loss_type_str, "auc") == 0) {
        loss = new AUCLoss();
        delete p;
        p = static_cast<PermutationMap*>(new AUCPermutationMap());
    } else if (strcmp(loss_type_str, "fscore") == 0) {
        loss = new FScore(weight);
        delete p;
        p = static_cast<PermutationMap*>(new FScorePermutationMap());
    } else {
        loss = new Accuracy();
    }
    //printf("tree address before cachetree is: %p\n", tree);
    tree = new CacheTree(loss, nsamples, nrules, c, 1, rules, labels, meta, minority_class, 0, ablation, calculate_size, type);
    if (verbosity_set.count("progress"))
        printf("%s", run_type);
    // runs our algorithm
    bbound(tree, max_num_nodes, q, p, 1.0);

    const tracking_vector<unsigned short, DataStruct::Tree>& r_list = tree->opt_rulelist();

    double accuracy = 1.0 - tree->min_objective() + c*r_list.size();

    if (verbosity_set.count("progress")) {
        printf("final num_nodes: %zu\n", tree->num_nodes());
        printf("final num_evaluated: %zu\n", tree->num_evaluated());
        printf("final min_objective: %1.5f\n", tree->min_objective());
        printf("final accuracy: %1.5f\n", accuracy);
   }

    for(size_t i = 0; i < tree->opt_rulelist().size(); i++) {
        rulelist->push_back(tree->opt_rulelist()[i]);
        classes->push_back(tree->opt_predictions()[i]);
    }
    classes->push_back(tree->opt_predictions().back());
    print_final_rulelist(r_list, tree->opt_predictions(),
                     latex_out, rules, labels, opt_fname, verbosity_set.count("progress"));

    if (verbosity_set.count("progress"))
        printf("final total time: %f\n", time_diff(init));

    logger->dumpState();
    logger->closeFile();
    if (verbosity_set.count("progress")) {
        printf("\ndelete tree\n");
    }
    delete tree;
    if (verbosity_set.count("progress")) {
        printf("\ndelete symmetry-aware map\n");
    }
    delete p;
    if (verbosity_set.count("progress")) {
        printf("\ndelete priority queue\n");
    }
    delete q;

    tree = nullptr;
    queue = nullptr;
    pmap = nullptr;
    /*if (verbosity.count("progress")) {
        printf("\ndelete logger\n");
    }
    delete logger;*/

    return accuracy;
}

} // end extern C