#include "features.hh"
#include "queue.hh"
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <getopt.h>
#include <pthread.h>
#include <stdio.h>

#define BUFSZ 512

std::mutex log_lk;
std::mutex min_obj_lk;
int lock_ac = 0;
//pthread_rwlock_t pmap_lk = PTHREAD_RWLOCK_INITIALIZER;

/*
 * Logs statistics about the execution of the algorithm and dumps it to a file.
 * To turn off, pass verbosity <= 1
 */
NullLogger* logger;

FeatureToggle* featureDecisions;

int main(int argc, char *argv[]) {
    const char usage[] = "USAGE: %s [-b] [-t num_threads]"
        "[-n max_num_nodes] [-r regularization] [-v verbosity] "
        "-c (1|2|3|4) -p (0|1|2) [-f logging_frequency]"
        "-a (0|1|2) [-L latex_out] [-k random_seed]"
        "data.out data.label\n\n"
        "%s\n";

    extern char *optarg;
    bool run_bfs = false;
    bool run_curiosity = false;
    int curiosity_policy = 0;
    bool latex_out = false;
    bool use_prefix_perm_map = false;
    bool use_captured_sym_map = false;
    int verbosity = 0;
    int map_type = 0;
    int max_num_nodes = 100000;
    double c = 0.01;
    char ch;
    bool error = false;
    char error_txt[BUFSZ];
    int freq = 1000;
    size_t num_threads = 1;
    int ablation = 0;
    bool calculate_size = false;
    int iterno = 0;
    size_t random_seed = 89;
    /* only parsing happens here */
    while ((ch = getopt(argc, argv, "bsLc:p:v:n:r:f:a:t:i:k:")) != -1) {
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
            verbosity = atoi(optarg);
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
        case 't':
            num_threads = atoi(optarg);
            break;
        case 'i':
            iterno = atoi(optarg);
            break;
        case 'k':
            random_seed = atoi(optarg);
            break;
        default:
            error = true;
            snprintf(error_txt, BUFSZ, "unknown option: %c", ch);
        }
    }
    if (max_num_nodes < 0) {
        error = true;
        snprintf(error_txt, BUFSZ, "number of nodes must be positive");
    }
    if (c < 0) {
        error = true;
        snprintf(error_txt, BUFSZ, "regularization constant must be postitive");
    }
    if (map_type > 2 || map_type < 0) {
        error = true;
        snprintf(error_txt, BUFSZ, "symmetry-aware map must be (0|1|2)");
    }
    if ((run_bfs + run_curiosity) != 1) {
        error = true;
        snprintf(error_txt, BUFSZ,
                "you must use at least and at most one of (-b | -c)");
    }
    if (argc < 2 + optind) {
        error = true;
        snprintf(error_txt, BUFSZ,
                "you must specify data files for rules and labels");
    }
    if (run_curiosity && !((curiosity_policy >= 1) && (curiosity_policy <= 4))) {
        error = true;
        snprintf(error_txt, BUFSZ,
                "you must specify a curiosity type (1|2|3|4)");
    }

    if (error) {
        fprintf(stderr, usage, argv[0], error_txt);
        exit(1);
    }

    std::map<int, std::string> curiosity_map;
    curiosity_map[1] = "curiosity";
    curiosity_map[2] = "curious_lb";
    curiosity_map[3] = "curious_obj";
    curiosity_map[4] = "dfs";

    argc -= optind;
    argv += optind;

    int nrules, nsamples, nlabels, nsamples_chk, errno;
    rule_t *rules, *labels;
    errno = rules_init(argv[0], &nrules, &nsamples, &rules, 1);
    printf("Rules ERRNO: %d\n", errno);
    errno = rules_init(argv[1], &nlabels, &nsamples_chk, &labels, 0);
    printf("Labels ERRNO: %d\n", errno);

    int nmeta, nsamples_check;
    // Equivalent points information is precomputed, read in from file, and stored in meta
    rule_t *meta;
    if (argc == 3)
        rules_init(argv[2], &nmeta, &nsamples_check, &meta, 0);
    else
        meta = NULL;

    if (verbosity >= 10)
        print_machine_info();
    char froot[BUFSZ];
    char log_fname[BUFSZ];
    char opt_fname[BUFSZ];
    const char* pch = strrchr(argv[0], '/');
    snprintf(froot, BUFSZ, "../logs/for-%s-%s%s-%s-%s-removed=%s-t=%lu-max_num_nodes=%d-c=%.7f-v=%d-f=%d-i=%d",
            pch ? pch + 1 : "",
            run_bfs ? "bfs" : "",
            run_curiosity ? curiosity_map[curiosity_policy].c_str() : "",
            use_prefix_perm_map ? "with_prefix_perm_map" :
                (use_captured_sym_map ? "with_captured_symmetry_map" : "no_pmap"),
            meta ? "minor" : "no_minor",
            ablation ? ((ablation == 1) ? "support" : "lookahead") : "none",
            num_threads, max_num_nodes, c, verbosity, freq, iterno);
    snprintf(log_fname, BUFSZ, "%s.txt", froot);
    snprintf(opt_fname, BUFSZ, "%s-opt.txt", froot);

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

    // TOOD: use cmd line params to initialize features
    featureDecisions = new FeatureToggle(false);

    double init = timestamp();
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

    // Create permutation map and synchronization primitives for it.
    PermutationMap* p;
    if (use_prefix_perm_map) {
        strcat(run_type, " Prefix Map\n");
        PrefixPermutationMap* prefix_pmap = new PrefixPermutationMap;
        p = (PermutationMap*) prefix_pmap;
    } else if (use_captured_sym_map) {
        strcat(run_type, " Captured Symmetry Map\n");
        CapturedPermutationMap* cap_pmap = new CapturedPermutationMap;
        p = (PermutationMap*) cap_pmap;
    } else {
        strcat(run_type, " No Permutation Map\n");
        NullPermutationMap* null_pmap = new NullPermutationMap;
        p = (PermutationMap*) null_pmap;
    }

    CacheTree* tree = new CacheTree(nsamples, nrules, c, num_threads,
        rules, labels, meta, ablation, calculate_size, type, random_seed);
    printf("%s", run_type);

    // Initialize logger
    bbound_init(tree);

    // Set up per-thread queues
    std::thread* threads = new std::thread[num_threads];

    Queue qs[num_threads];
    for(size_t i = 0; i < num_threads; ++i) {
        qs[i] = Queue(cmp, run_type);
        tracking_vector<unsigned short, DataStruct::Tree>* init_rules = tree->get_subrange(i);
        InternalRoot* iroot = new InternalRoot(tree->root(), init_rules);
        qs[i].push(iroot);
    }

    SharedQueue* shared_q = new SharedQueue();

	// Let the threads loose
	for(size_t i = 0; i < num_threads; ++i) {
	    threads[i] = std::thread(bbound, tree, max_num_nodes, &qs[i], p, i, shared_q);
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
    printf("final num_evaluated: %zu\n", tree->num_evaluated());
    printf("final min_objective: %1.5f\n", tree->min_objective());
    const tracking_vector<unsigned short, DataStruct::Tree>& r_list = tree->opt_rulelist();
    printf("final accuracy: %1.5f\n",
           1 - tree->min_objective() + c*r_list.size());
    print_final_rulelist(r_list, tree->opt_predictions(),
                         latex_out, rules, labels, opt_fname);

    printf("final total time: %f\n", time_diff(init));
    printf("Number of lock acquistions: %d\n", lock_ac);
    printf("Number of shared_queue acquisitions: %d\n", shared_q->n_acc());

    logger->dumpState();
    logger->closeFile();
    if (meta) {
        printf("\ndelete identical points indicator");
        rules_free(meta, nmeta, 0);
    }
    printf("\ndelete rules\n");
    rules_free(rules, nrules, 1);
    printf("delete labels\n");
    rules_free(labels, nlabels, 0);
    printf("tree destructors\n");
    return 0;
}
