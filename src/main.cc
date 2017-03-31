#include "queue.hh"
//#include "utils.hh"
//#include "pmap.hh"

#include <gperftools/heap-profiler.h>
#include <iostream>
#include <stdio.h>
#include <getopt.h>
#include <numeric>
#include <thread>
#include <mutex>
#include <pthread.h>

Logger logger;
int ablation = 0;
std::mutex log_lk;
int lock_ac = 0;
pthread_rwlock_t pmap_lk = PTHREAD_RWLOCK_INITIALIZER;

int main(int argc, char *argv[]) {
    const char usage[] = "USAGE: %s [-s] [-b] [-t num_threads]"
        "[-n max_num_nodes] [-r regularization] [-v verbosity] "
        "-c (1|2|3|4) -p (0|1|2) [-f logging_frequency]"
        "-a (1|2|3)"
        "data.out data.label\n\n"
        "%s\n";

    extern char *optarg;
    bool run_stochastic = false;
    bool run_bfs = false;
    bool run_curiosity = false;
    int curiosity_policy = 0;
    bool latex_out = false;
    bool run_pmap = false;
    bool use_prefix_perm_map = false;
    bool use_captured_sym_map = false;
    int verbosity = 1;
    int max_num_nodes = 100000;
    double c = 0.01;
    char ch;
    bool error = false;
    char error_txt[512];
    int freq = 1000;
    size_t switch_iter = 0;
    size_t num_threads = 1;
    /* only parsing happens here */
    while ((ch = getopt(argc, argv, "sbLc:p:v:n:r:f:w:a:t:")) != -1) {
        switch (ch) {
        case 's':
            run_stochastic = true;
            break;
        case 'b':
            run_bfs = true;
            break;
        case 'c':
            run_curiosity = true;
            curiosity_policy = atoi(optarg);
            break;
        case 'L':
            latex_out = true;
            break;
        case 'p':
	        run_pmap = true;
            use_prefix_perm_map = atoi(optarg) == 1;
            use_captured_sym_map = atoi(optarg) == 2;
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
        case 'w':
            switch_iter = atoi(optarg);
            break;
        case 'a':
            ablation = atoi(optarg);
            break;
        case 't':
            num_threads = atoi(optarg);
            break;
        default:
            error = true;
            sprintf(error_txt, "unknown option: %c", ch);
        }
    }
    if ((run_stochastic + run_bfs + run_curiosity) != 1) {
        error = true;
        sprintf(error_txt,
                "you must use at least and at most one of (-s | -b | -c)");
    }
    if (argc < 2 + optind) {
        error = true;
        sprintf(error_txt,
                "you must specify data files for rules and labels");
    }
    if (run_curiosity && !((curiosity_policy >= 1) && (curiosity_policy <= 4))) {
        error = true;
        sprintf(error_txt,
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
    rule_t *meta;
    if (argc == 3)
        rules_init(argv[2], &nmeta, &nsamples_check, &meta, 0);
    else
        meta = NULL;

    print_machine_info();
    char froot[512];
    char log_fname[512];
    char opt_fname[512];
    const char* pch = strrchr(argv[0], '/');
    sprintf(froot, "../logs/for-%s-%s%s%s-%s-%s-removed=%s-t=%d-max_num_nodes=%d-c=%.7f-v=%d-f=%d",
            pch ? pch + 1 : "",
            run_stochastic ? "stochastic" : "",
            run_bfs ? "bfs" : "",
            run_curiosity ? curiosity_map[curiosity_policy].c_str() : "",
            run_pmap ? (use_prefix_perm_map ? "with_prefix_perm_map" : "with_captured_symmetry_map") : "no_pmap",
            meta ? "minor" : "no_minor",
            ablation ? ((ablation == 1) ? "support" : "lookahead") : "none",
            num_threads, max_num_nodes, c, verbosity, freq);
    sprintf(log_fname, "%s.txt", froot);
    sprintf(opt_fname, "%s-opt.txt", froot);

    if (verbosity >= 1000) {
        printf("\n%d rules %d samples\n\n", nrules, nsamples);
        rule_print_all(rules, nrules, nsamples);

        printf("\nLabels (%d) for %d samples\n\n", nlabels, nsamples);
        rule_print_all(labels, nlabels, nsamples);
    }

    logger.setC(c);
    logger.setNRules(nrules);
    logger.initPrefixVec();
    logger.setVerbosity(verbosity);
    logger.setLogFileName(log_fname);
    logger.setFrequency(freq);

	// Variables for tree iteration
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
    BaseQueue* qs;
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
    }

    PermutationMap* p;
    if (use_prefix_perm_map) {
        printf("Prefix Permutation Map\n");
        PrefixPermutationMap* prefix_pmap = new PrefixPermutationMap;
        //prefix_pmap->bucket_count();
        p = (PermutationMap*) prefix_pmap;
    } else if (use_captured_sym_map) {
        printf("Captured Symmetry Map\n");
        CapturedPermutationMap* cap_pmap = new CapturedPermutationMap;
        p = (PermutationMap*) cap_pmap;
    } else {
        printf("No Permutation Map \n");
        NullPermutationMap* null_pmap = new NullPermutationMap;
        p = (PermutationMap*) null_pmap;
    }

    CacheTree** trees = new CacheTree*[num_threads];
    if (run_stochastic) {
        if (use_prefix_perm_map) {
            printf("BBOUND_STOCHASTIC Permutation Map\n");
        } else if (use_captured_sym_map) {
            printf("BBOUND_STOCHASTIC Captured Symmetry Map\n");
        } else {
            printf("BBOUND_STOCHASTIC No Permutation Map \n");
        }
        CacheTree tree(nsamples, nrules, c, rules, labels, meta);
        bbound_stochastic(&tree, max_num_nodes, p);
        printf("final num_nodes: %zu\n", tree.num_nodes());
        printf("final num_evaluated: %zu\n", tree.num_evaluated());
        printf("final min_objective: %1.5f\n", tree.min_objective());
        const std::vector<unsigned short, cache_alloc<unsigned short> >& r_list = tree.opt_rulelist();
        //auto r_list = tree.opt_rulelist();
        printf("final accuracy: %1.5f\n",
               1 - tree.min_objective() + c*r_list.size());
        print_final_rulelist(r_list, tree.opt_predictions(),
                             latex_out, rules, labels, opt_fname);

        logger.dumpState();
    } else {
        //BFS
        for(size_t i = 0; i < num_threads; ++i) {
            trees[i] = new CacheTree(nsamples, nrules, c, rules, labels, meta);
            trees[i]->open_print_file(i + 1, num_threads);
            threads[i] = std::thread(bbound_queue_init, trees[i], &qs[i],
                                  p, ranges[i], min_objective);
            threads[i].join();
            threads[i] = std::thread(bbound_queue, trees[i], max_num_nodes/num_threads,
                                  &qs[i], p, 0, 0, min_objective);
        }
        for(size_t i = 0; i < num_threads; ++i) {
            threads[i].join();
            CacheTree* tree = trees[i];
            printf("final num_nodes: %zu\n", tree->num_nodes());
            printf("final num_evaluated: %zu\n", tree->num_evaluated());
            printf("final min_objective: %1.5f\n", tree->min_objective());
            const std::vector<unsigned short, cache_alloc<unsigned short> >& r_list = tree->opt_rulelist();
            //auto r_list = tree.opt_rulelist();
            printf("final accuracy: %1.5f\n",
                   1 - tree->min_objective() + c*r_list.size());
            print_final_rulelist(r_list, tree->opt_predictions(),
                                 latex_out, rules, labels, opt_fname);
        }
    }
    /*
    } else if (run_curiosity) {
        bbound_queue(tree, max_num_nodes, q, p, 0, 0, NULL);
        CacheTree* tree = new CacheTree(nsamples, nrules, c, rules, labels, meta);
        //CuriousCacheTree tree(nsamples, nrules, c, rules, labels, meta);
        if (curiosity_policy == 1) {
            if (use_prefix_perm_map) {
                printf("Curiosity Prefix Permutation Map\n");
            } else if (use_captured_sym_map) {
                printf("Curiosity Captured Symmetry Map\n");
            } else {
                printf("Curiosity No Permutation Map \n");
            }
            CuriousQueue* curious_q = new CuriousQueue;
            size_t num_iter;
            num_iter = bbound_queue(tree, max_num_nodes,
                                    curious_q, p, 0, switch_iter, NULL);
            if (curious_q->size() > 0) {
                printf("Switching to curious lower bound policy... \n");
                LowerBoundQueue* curious_lb_q = new LowerBoundQueue;
                while(!curious_q->empty()) {
                    CuriousNode* selected_node = (CuriousNode*) curious_q->front();
                    curious_q->pop();
                    curious_lb_q->push(selected_node);
                }
                num_iter = bbound_queue(tree, max_num_nodes,
                                        curious_lb_q, p, num_iter, 0, NULL);
            }
        } else if (curiosity_policy == 2) {
            if (use_prefix_perm_map) {
                printf("CURIOUS LOWER BOUND Prefix Permutation Map\n");
            } else if (use_captured_sym_map) {
                printf("CURIOUS LOWER BOUND Captured Symmetry Map\n");
            } else {
                printf("CURIOUS LOWER BOUND No Permutation Map\n");
            }
            LowerBoundQueue* lb_q = new LowerBoundQueue;
            bbound_queue(tree, max_num_nodes, lb_q, p, 0, 0, NULL);
        } else if (curiosity_policy == 3) {
            if (use_prefix_perm_map) {
                printf("CURIOUS OBJECTIVE Prefix Permutation Map\n");
            } else if (use_captured_sym_map) {
                printf("CURIOUS OBJECTIVE Captured Symmetry Map\n");
            } else {
                printf("CURIOUS OBJECTIVE No Permutation Map\n");
            }
            ObjectiveQueue* objective_q = new ObjectiveQueue;
            bbound_queue(tree, max_num_nodes, objective_q, p, 0, 0, NULL);
        } else if (curiosity_policy == 4) {
            if (use_prefix_perm_map) {
                printf("DFS Prefix Permutation Map\n");
            } else if (use_captured_sym_map) {
                printf("DFS Captured Symmetry Map\n");
            } else {
                printf("DFS No Permutation Map\n");
            }
            DFSQueue* dfs_q = new DFSQueue;
            bbound_queue(tree, max_num_nodes, dfs_q, p, 0, 0, NULL);
        }
        printf("final num_nodes: %zu\n", tree->num_nodes());
        printf("final num_evaluated: %zu\n", tree->num_evaluated());
        printf("final min_objective: %1.5f\n", tree->min_objective());
        const std::vector<unsigned short, cache_alloc<unsigned short> >& r_list = tree->opt_rulelist();
        //auto r_list = tree->opt_rulelist();
        printf("final accuracy: %1.5f\n",
               1 - tree->min_objective() + c*r_list.size());
        print_final_rulelist(r_list, tree->opt_predictions(),
                             latex_out, rules, labels, opt_fname);
    }
        */

/*    printf("final num_nodes: %zu\n", num_nodes);
    printf("final num_evaluated: %zu\n", num_evaluated);
    printf("final min_objective: %1.5f\n", *min_objective);
    printf("final accuracy: %1.5f\n", 1 - *min_objective + c*r_list.size());
    print_final_rulelist(r_list, preds, latex_out, rules, labels, opt_fname);
    */
    logger.dumpState();
    free(min_objective);

    printf("final total time: %f\n", time_diff(init));
    printf("Number of lock acquistions: %d\n", lock_ac);
    logger.dumpState();
    logger.closeFile();
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
