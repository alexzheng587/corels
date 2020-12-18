#include "queue.hh"
#include "loss.hh"
#include <iostream>
#include <stdio.h>
#include <getopt.h>
#include <unistd.h>
#include <string>

#define BUFSZ 512

/*
 * Logs statistics about the execution of the algorithm and dumps it to a file.
 * To turn off, do not pass "log" as verbosity parameter
 */
NullLogger* logger;

int main(int argc, char *argv[]) {
    const char usage[] = "USAGE: %s [-b] "
        "[-l (acc|bacc|wacc|fscore)] [-w weight] "
        "[-n max_num_nodes] [-r regularization] [-v (rule|label|samples|progress|log|silent)] "
        "-c (1|2|3|4) -p (0|1|2) [-k rulelist_save_number] [-f logging_frequency] "
        "-a (0|1|2) [-s] [-L latex_out] "
        "data.out data.label [data.minor]\n\n"
        "%s\n";

    extern char *optarg;
    std::string loss_type_str = "acc";
    double weight = 1.0;
    bool run_bfs = false;
    bool run_curiosity = false;
    int curiosity_policy = 0;
    bool latex_out = false;
    bool use_prefix_perm_map = false;
    bool use_captured_sym_map = false;
    char *vopt, *verb_trim;
    std::set<std::string> verbosity;
    bool verr = false;
    const char *vstr = "rule|label|samples|progress|log|silent";
    int map_type = 0;
    int max_num_nodes = 100000;
    double c = 0.01;
    int k = 1;
    char ch;
    bool error = false;
    char error_txt[BUFSZ];
    int freq = 1000;
    int ablation = 0;
    double min_obj = 1.0;
    bool calculate_size = false;
    char verbstr[BUFSZ];
    int iterno = 0;
    /* only parsing happens here */
    while ((ch = getopt(argc, argv, "bsLc:k:p:v:n:r:f:a:u:i:l:w:")) != -1) {
        switch (ch) {
        case 'b':
            run_bfs = true;
            break;
        case 's':
            calculate_size = true;
            break;
        case 'c':
            run_curiosity = true;
            curiosity_policy = atoi(optarg);
            break;
        case 'L':
            latex_out = true;
            break;
        case 'p':
            map_type = atoi(optarg);
            use_prefix_perm_map = map_type == 1;
            use_captured_sym_map = map_type == 2;
            break;
        case 'v':
            verb_trim = strtok(optarg, " ");
            strcpy(verbstr, verb_trim);
            vopt = strtok(verb_trim, ",");
            while (vopt != NULL) {
                if (!strstr(vstr, vopt)) {
                    verr = true;
                }
                verbosity.insert(vopt);
                vopt = strtok(NULL, ",");
            }
            break;
        case 'n':
            max_num_nodes = atoi(optarg);
            break;
        case 'r':
            c = atof(optarg);
            break;
        case 'f':
            freq = atoi(optarg);
            break;
        case 'a':
            ablation = atoi(optarg);
            break;
        case 'k':
            k = atoi(optarg);
        case 'u':
            min_obj = atof(optarg);
            break;
        case 'i':
            iterno = atoi(optarg);
            break;
        case 'l':
            loss_type_str = optarg;
            break;
        case 'w':
            weight = atof(optarg);
            break;
        default:
            error = true;
            snprintf(error_txt, BUFSZ, "unknown option: %c", ch);
        }
    }
    argc -= optind;
    argv += optind;
    if (loss_type_str != "acc" && loss_type_str != "bacc" && loss_type_str != "wacc" &&
    loss_type_str != "auc" && loss_type_str != "fscore") {
        error = true;
        snprintf(error_txt, BUFSZ, "loss type must be one of (acc|bacc|wacc|auc|fscore)");
    }
    if (weight < 0) {
        error = true;
        snprintf(error_txt, BUFSZ, "weight must be positive");
    }
    if (max_num_nodes < 0) {
        error = true;
        snprintf(error_txt, BUFSZ, "number of nodes must be positive");
    }
    if (c < 0) {
        error = true;
        snprintf(error_txt, BUFSZ, "regularization constant must be postitive");
    }
    if (k < 1) {
        error = true;
        snprintf(error_txt, BUFSZ, "number of rulelists to save must be greater than 0");
    }
    if (map_type > 2 || map_type < 0) {
        error = true;
        snprintf(error_txt, BUFSZ, "symmetry-aware map must be (0|1|2)");
    }
    if ((run_bfs + run_curiosity) != 1) {
        error = true;
        snprintf(error_txt, BUFSZ,
                "you must use exactly one of (-b | -c)");
    }
    if (argc < 2) {
        error = true;
        snprintf(error_txt, BUFSZ,
                "you must specify data files for rules and labels");
    }
    if (access(argv[0], R_OK) < 0 || access(argv[1], R_OK) < 0) {
        error = true;
        snprintf(error_txt, BUFSZ,
                "specified data files for rules and/or labels do not exist");
    }
    if (run_curiosity && !((curiosity_policy >= 1) && (curiosity_policy <= 4))) {
        error = true;
        snprintf(error_txt, BUFSZ,
                "you must specify a curiosity type (1|2|3|4)");
    }
    if (verr) {
        error = true;
        snprintf(error_txt, BUFSZ,
                 "verbosity options must be one or more of (rule|label|samples|progress|log|silent), separated with commas (i.e. -v progress,log)");
    }
    if (verbosity.count("samples") && !(verbosity.count("rule") || verbosity.count("label"))) {
        error = true;
        snprintf(error_txt, BUFSZ,
                 "verbosity 'samples' option must be combined with at least one of (rule|label)");
    }
    if (verbosity.size() > 2 && verbosity.count("silent")) {
        snprintf(error_txt, BUFSZ,
                 "verbosity 'silent' option must be passed without any additional verbosity parameters");
    }

    char* progname = argv[0];

    if (error) {
        fprintf(stderr, usage, progname, error_txt);
        exit(1);
    }

    // default: show progress
    if (verbosity.size() == 0) {
        verbosity.insert("progress");
    }

    if (verbosity.count("silent")) {
        verbosity.clear();
    }

    std::map<int, std::string> curiosity_map;
    curiosity_map[1] = "curiosity";
    curiosity_map[2] = "curious_lb";
    curiosity_map[3] = "curious_obj";
    curiosity_map[4] = "dfs";

    int nrules, nsamples, nlabels, nsamples_chk;
    rule_t *rules, *labels;
    rules_init(argv[0], &nrules, &nsamples, &rules, 1);
    rules_init(argv[1], &nlabels, &nsamples_chk, &labels, 0);

    int nmeta, nsamples_check;
    // Equivalent points information is precomputed, read in from file, and stored in meta
    rule_t *meta;
    // If the objective is auc, then we expect the equivalent points file to contain all
    // classes and their number positives and negatives
    int nclasses, nclass_samples;
    minority_class_t *minor_class;
    if (loss_type_str != "auc") {
        if (argc == 3)
            rules_init(argv[2], &nmeta, &nsamples_check, &meta, 0);
        else
            meta = NULL;
    } else {
        if (argc == 3)
            auc_rules_init(argv[2], &nclasses, &nclass_samples, &minor_class, 0);
        else
            minor_class = NULL;
    }

    if (verbosity.count("log"))
        print_machine_info();
    char froot[BUFSZ];
    char log_fname[BUFSZ];
    char opt_fname[BUFSZ];
    const char* pch = strrchr(argv[0], '/');
    snprintf(froot, BUFSZ, "../logs/for-%s-%s%s-%s-%s-%s-removed=%s-max_num_nodes=%d-c=%.7f-v=%s-f=%d-i=%d",
            pch ? pch + 1 : "",
            run_bfs ? "bfs" : "",
            run_curiosity ? curiosity_map[curiosity_policy].c_str() : "",
            use_prefix_perm_map ? "with_prefix_perm_map" :
                (use_captured_sym_map ? "with_captured_symmetry_map" : "no_pmap"),
             loss_type_str.c_str(),
            meta ? "minor" : "no_minor",
            ablation ? ((ablation == 1) ? "support" : "lookahead") : "none",
             max_num_nodes, c, verbstr, freq, iterno);
    snprintf(log_fname, BUFSZ, "%s.txt", froot);

    if (verbosity.count("rule")) {
        printf("\n%d rules %d samples\n\n", nrules, nsamples);
        rule_print_all(rules, nrules, nsamples, (verbosity.count("samples")));
    }

    if (verbosity.count("label")) {
        printf("\nLabels (%d) for %d samples\n\n", nlabels, nsamples);
        rule_print_all(labels, nlabels, nsamples, (verbosity.count("samples")));
    }

    if (verbosity.count("log")) {
        logger = new Logger(c, nrules, verbosity, log_fname, freq);
    } else {
        logger = new NullLogger();
        logger->setVerbosity(verbosity);
    }
    double init = timestamp();
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
    if (use_prefix_perm_map) {
        strcat(run_type, " Prefix Map");
        PrefixPermutationMap* prefix_pmap = new PrefixPermutationMap;
        p = (PermutationMap*) prefix_pmap;
    } else if (use_captured_sym_map) {
        strcat(run_type, " Captured Symmetry Map");
        CapturedPermutationMap* cap_pmap = new CapturedPermutationMap;
        p = (PermutationMap*) cap_pmap;
    } else {
        strcat(run_type, " No Permutation Map");
        NullPermutationMap* null_pmap = new NullPermutationMap;
        p = (PermutationMap*) null_pmap;
    }

    strcat(run_type, " and loss function type ");
    strcat(run_type, loss_type_str.c_str());
    strcat(run_type, "\n");

    Loss* loss;
    if (loss_type_str == "bacc") {
        loss = new BalancedAccuracy();
    } else if (loss_type_str == "wacc") {
        loss = new WeightedAccuracy(weight);
    } else if (loss_type_str == "auc") {
        loss = new AUCLoss();
        delete p;
        p = static_cast<PermutationMap*>(new AUCPermutationMap());
    } else if (loss_type_str == "fscore") {
        loss = new FScore(weight);
        delete p;
        p = static_cast<PermutationMap*>(new FScorePermutationMap());
    } else {
        loss = new Accuracy();
    }

    CacheTree* tree = new CacheTree(loss, nsamples, nrules, c, k, rules, labels, meta, minor_class, nclasses, ablation, calculate_size, type);
    if (verbosity.count("progress"))
        printf("%s", run_type);
    // runs our algorithm
    bbound(tree, max_num_nodes, q, p, min_obj);

    for(int i = (tree->nrulelists() - 1); i > -1; i--) {
        const tracking_vector<unsigned short, DataStruct::Tree>& r_list = tree->opt_rulelist(i);

        // don't save default rule
        if(r_list.size() == 0)
            continue;

        snprintf(opt_fname, BUFSZ, "%s-%d-opt.txt", froot, i);

        print_final_rulelist(r_list, tree->opt_predictions(i),
                     latex_out, rules, labels, opt_fname, verbosity.count("progress"));

        printf("Rulelist #%d min_objective: %1.5f\n", i, tree->min_objective(i));
        printf("Rulelist #%d accuracy: %1.5f\n",
            i, 1 - tree->min_objective(i) + c*r_list.size());
    }

    if (verbosity.count("progress")) {
        printf("\nfinal num_nodes: %zu\n", tree->num_nodes());
        printf("final num_evaluated: %zu\n", tree->num_evaluated());
        printf("final min_objective: %1.5f\n", tree->min_objective());
        printf("final accuracy: %1.5f\n",
            1 - tree->min_objective() + c*tree->opt_rulelist().size());
    }

    if (verbosity.count("progress"))
        printf("final total time: %f\n", time_diff(init));

    logger->dumpState();
    logger->closeFile();

    if (meta) {
        if (verbosity.count("progress"))
            printf("\ndelete identical points indicator");
        rules_free(meta, nmeta, 0);
    }
    if (verbosity.count("progress"))
        printf("\ndelete rules\n");
    rules_free(rules, nrules, 1);
    if (verbosity.count("progress"))
        printf("delete labels\n");
    rules_free(labels, nlabels, 0);
    if (verbosity.count("progress"))
        printf("tree destructors\n");
    return 0;
}
