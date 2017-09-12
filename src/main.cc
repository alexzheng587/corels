#include "queue.hh"
#include <iostream>
#include <numeric>
#include <mutex>
#include <string>
#include <thread>
#include <getopt.h>
#include <pthread.h>
#include <stdio.h>

#define BUFSZ 512

std::mutex log_lk;
int lock_ac = 0;
std::mutex lk_;
//pthread_rwlock_t pmap_lk = PTHREAD_RWLOCK_INITIALIZER;

/*
 * Logs statistics about the execution of the algorithm and dumps it to a file.
 * To turn off, pass verbosity <= 1
 */
NullLogger* logger;

int main(int argc, char *argv[]) {
    const char usage[] = "USAGE: %s [-b] [-t num_threads]"
        "[-n max_num_nodes] [-r regularization] [-v verbosity] "
        "-c (1|2|3|4) -p (0|1|2) [-f logging_frequency]"
        "-a (0|1|2) [-s] [-L latex_out]"
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
    /* only parsing happens here */
    while ((ch = getopt(argc, argv, "bsLc:p:v:n:r:f:a:t:")) != -1) {
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
    snprintf(froot, BUFSZ, "../logs/for-%s-%s%s-%s-%s-removed=%s-t=%d-max_num_nodes=%d-c=%.7f-v=%d-f=%d",
            pch ? pch + 1 : "",
            run_bfs ? "bfs" : "",
            run_curiosity ? curiosity_map[curiosity_policy].c_str() : "",
            use_prefix_perm_map ? "with_prefix_perm_map" : 
                (use_captured_sym_map ? "with_captured_symmetry_map" : "no_pmap"),
            meta ? "minor" : "no_minor",
            ablation ? ((ablation == 1) ? "support" : "lookahead") : "none",
            num_threads, max_num_nodes, c, verbosity, freq);
    snprintf(log_fname, BUFSZ, "%s.txt", froot);
    snprintf(opt_fname, BUFSZ, "%s-opt.txt", froot);

    if (verbosity >= 1000) {
        printf("\n%d rules %d samples\n\n", nrules, nsamples);
        rule_print_all(rules, nrules, nsamples);

        printf("\nLabels (%d) for %d samples\n\n", nlabels, nsamples);
        rule_print_all(labels, nlabels, nsamples);
    }

    if (verbosity > 1)
        logger = new Logger(c, nrules, verbosity, log_fname, freq);
    else
        logger = new NullLogger();


    size_t num_nodes = 0;
    size_t num_evaluated = 0;
    //std::vector<unsigned short> r_list;
    //std::vector<unsigned short, cache_alloc<unsigned short> > r_list;
    //std::vector<bool, cache_alloc<bool> > preds;

    double* min_objective = (double*) malloc(sizeof(double*));
    *min_objective = 1.0;
    std::thread* threads = new std::thread[num_threads];
    std::vector<unsigned short> indices(nrules - 1);
    std::iota(indices.begin(), indices.end(), 1);
    std::random_shuffle(indices.begin(), indices.end());
    std::vector<unsigned short>* ranges = new std::vector<unsigned short>[num_threads];
    for(size_t i = 0; i < num_threads; ++i) {
        printf("RANGE: %u-%u\n", indices[(size_t)(i * ((nrules-1)/(float)num_threads))], indices[(size_t)((i + 1) * ((nrules-1)/(float)num_threads))]);
        printf("RANGE INDICES: %zu-%zu\n", (size_t)(i * ((nrules-1)/(float)num_threads)), (size_t)((i + 1) * ((nrules-1)/(float)num_threads)));
        std::vector<unsigned short> range(&indices[(size_t)(i * ((nrules-1)/(float)num_threads))], 
                                                      &indices[(size_t)((i + 1) * ((nrules-1)/(float)num_threads))]);
        ranges[i] = range;
    }
    //pthread_rwlockattr_t attr;
    //pthread_rwlockattr_init(&attr);
    //pthread_rwlock_init(pmap_lk, &attr);

    double init = timestamp();
    char run_type[BUFSZ];

    /*BaseQueue* qs;
    if (curiosity_policy == 1) {
        printf("Curious ");
        CuriousQueue* curious_q = new CuriousQueue[num_threads];
        //for (size_t i = 0; i < num_threads; ++i)
        //    curious_q[i] = new CuriousQueue;
        qs = (BaseQueue*) curious_q;
    } else if (curiosity_policy == 2) {
        printf("Lower Bound ");
        LowerBoundQueue* lb_q = new LowerBoundQueue[num_threads];
        //for (size_t i = 0; i < num_threads; ++i)
        //    lb_q[i] = new LowerBoundQueue;
        qs = (BaseQueue*) lb_q;
    } else if (curiosity_policy == 3) {
        printf("Objective ");
        ObjectiveQueue* obj_q = new ObjectiveQueue[num_threads];
        //for (size_t i = 0; i < num_threads; ++i)
        //    obj_q[i] = new ObjectiveQueue;
        qs = (BaseQueue*) obj_q;
    } else if (curiosity_policy == 4) {
        printf("DFS ");
        DFSQueue* dfs_q = new DFSQueue[num_threads];
        //for (size_t i = 0; i < num_threads; ++i)
        //    dfs_q[i] = new DFSQueue;
        qs = (BaseQueue*) dfs_q;
    } else {
        printf("BFS ");
        qs = new BaseQueue[num_threads];
        //for (size_t i = 0; i < num_threads; ++i)
        //    qs[i] = new BaseQueue;
    }*/

    strcpy(run_type, "LEARNING RULE LIST via ");
    char const *type = "node";
    std::function<bool(Node*, Node*)> cmp;
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
    Queue qs[num_threads];
    for(int i = 0; i < num_threads; ++i) {
        qs[i] = Queue(cmp, run_type);
    }

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

    CacheTree** trees = new CacheTree*[num_threads];
    //CacheTree* tree = new CacheTree(nsamples, nrules, c, rules, labels, meta, ablation, calculate_size, type);
    printf("%s", run_type);
    // runs our algorithm

    //bbound(tree, max_num_nodes, q, p);
    for(size_t i = 0; i < num_threads; ++i) {
        trees[i] = new CacheTree(nsamples, nrules, c, rules, labels, meta, ablation, calculate_size, type);
        trees[i]->open_print_file(i + 1, num_threads);
        threads[i] = std::thread(bbound_init, trees[i], &qs[i],
                              p, ranges[i], min_objective);
        threads[i].join();
        threads[i] = std::thread(bbound, trees[i], max_num_nodes/num_threads,
                              &qs[i], p, min_objective);
    }
    for(size_t i = 0; i < num_threads; ++i) {
        threads[i].join();
        CacheTree* tree = trees[i];
        num_nodes += tree->num_nodes();
        num_evaluated += tree->num_evaluated();
        printf("final num_nodes: %zu\n", tree->num_nodes());
        printf("final num_evaluated: %zu\n", tree->num_evaluated());
        printf("final min_objective: %1.5f\n", tree->min_objective());
        const tracking_vector<unsigned short, DataStruct::Tree>& r_list = tree->opt_rulelist();
        printf("final accuracy: %1.5f\n",
               1 - tree->min_objective() + c*r_list.size());
        print_final_rulelist(r_list, tree->opt_predictions(),
                             latex_out, rules, labels, opt_fname);
    }

    printf("final total time: %f\n", time_diff(init));
    free(min_objective);
    printf("Number of lock acquistions: %d\n", lock_ac);

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
