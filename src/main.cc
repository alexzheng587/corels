#include "bbound.hh"
#include "utils.hh"

#include <iostream>
#include <stdio.h>
#include <queue>
#include <getopt.h>

Logger logger;

int main(int argc, char *argv[]) {
    const char usage[] = "USAGE: %s [-s] [-b] "
        "[-n max_num_nodes] [-r regularization] [-v verbosity] "
        "-c (1|2) -p (0|1|2) [-f logging_frequency]"
        "data.out data.label\n\n"
        "%s\n"; // for error

    extern char *optarg;
    bool run_stochastic = false;
    bool run_bfs = false;
    bool run_curiosity = false;
    bool use_curious_cmp = false;
    bool use_lower_bound_cmp = false;
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
    /* only parsing happens here */
    while ((ch = getopt(argc, argv, "sbLc:p:v:n:r:f:")) != -1) {
        switch (ch) {
        case 's':
            run_stochastic = true;
            break;
        case 'b':
            run_bfs = true;
            break;
        case 'c':
            run_curiosity = true;
            use_curious_cmp = atoi(optarg) == 1;
            use_lower_bound_cmp = atoi(optarg) == 2;
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
    if (run_curiosity && ((use_curious_cmp + use_lower_bound_cmp) != 1)) {
        error = true;
        sprintf(error_txt,
                "you must specify a curiosity type (1|2)");
    }

    if (error) {
        fprintf(stderr, usage, argv[0], error_txt);
        exit(1);
    }

    argc -= optind;
    argv += optind;

    int nrules, nsamples, nlabels, nsamples_chk;
    rule_t *rules, *labels;
    rules_init(argv[0], &nrules, &nsamples, &rules, 1);
    rules_init(argv[1], &nlabels, &nsamples_chk, &labels, 0);

    print_machine_info();
    char log_fname[512];
    char opt_fname[512];
    const char* pch = strrchr(argv[0], '/');
    sprintf(log_fname, "../logs/for-%s-%s%s%s-%s-max_num_nodes=%d-c=%.7f-v=%d-f=%d.txt",
            pch ? pch + 1 : "",
            run_stochastic ? "stochastic" : "",
            run_bfs ? "bfs" : "",
            run_curiosity ? (use_curious_cmp ? "curiosity" : "curious_lb") : "",
            run_pmap ? (use_prefix_perm_map ? "with_prefix_perm_map" : "with_captured_symmetry_map") : "no_pmap",
            max_num_nodes, c, verbosity, freq);

    sprintf(opt_fname, "../logs/for-%s-%s%s%s-%s-max_num_nodes=%d-c=%.7f-v=%d-f=%d-opt.txt",
            pch ? pch + 1 : "",
            run_stochastic ? "stochastic" : "",
            run_bfs ? "bfs" : "",
            run_curiosity ? (use_curious_cmp ? "curiosity" : "curious_lb") : "",
            run_pmap ? (use_prefix_perm_map ? "with_prefix_perm_map" : "with_captured_symmetry_map") : "no_pmap",
            max_num_nodes, c, verbosity, freq);

    if (verbosity >= 1000) {
        printf("\n%d rules %d samples\n\n", nrules, nsamples);
        rule_print_all(rules, nrules, nsamples);

        printf("\nLabels (%d) for %d samples\n\n", nlabels, nsamples);
        rule_print_all(labels, nlabels, nsamples);
    }
    if (verbosity >= 100) {
        printf("\nGreedy algorithm:\n\n");
        bbound_greedy(nsamples, nrules, rules, labels, 8);
    }

    logger.setNRules(nrules);
    logger.initPrefixVec();
    logger.setVerbosity(verbosity);
    logger.setLogFileName(log_fname);
    logger.setFrequency(freq);
    double init = timestamp();
    if (run_stochastic) {
        if (use_prefix_perm_map) {
            printf("BBOUND_STOCHASTIC Permutation Map\n");
            CacheTree<BaseNode> tree(nsamples, nrules, c, rules, labels);
            PrefixPermutationMap p;
            bbound_stochastic<BaseNode,
                              PrefixPermutationMap>(&tree, max_num_nodes,
                                                    &base_construct_policy,
                                                    &prefix_permutation_insert,
                                                    &prefix_map_garbage_collect,
                                                    &p);
            printf("final num_nodes: %zu\n", tree.num_nodes());
            printf("final num_evaluated: %zu\n", tree.num_evaluated());
            printf("final min_objective: %1.5f\n", tree.min_objective());
            const std::vector<size_t>& r_list = tree.opt_rulelist();
            printf("final accuracy: %1.5f\n",
                   1 - tree.min_objective() + c*r_list.size());
            print_final_rulelist(r_list, tree.opt_predictions(),
                                 latex_out, rules, labels, opt_fname);

            logger.dumpState();
        } else if (use_captured_sym_map) {
            printf("BBOUND_STOCHASTIC Captured Symmetry Map\n");
            CacheTree<BaseNode> tree(nsamples, nrules, c, rules, labels);
            CapturedPermutationMap p;
            bbound_stochastic<BaseNode,
                              CapturedPermutationMap>(&tree, max_num_nodes,
                                                      &base_construct_policy,
                                                      &captured_permutation_insert,
                                                      &captured_map_garbage_collect,
                                                      &p);
            printf("final num_nodes: %zu\n", tree.num_nodes());
            printf("final num_evaluated: %zu\n", tree.num_evaluated());
            printf("final min_objective: %1.5f\n", tree.min_objective());
            const std::vector<size_t>& r_list = tree.opt_rulelist();
            printf("final accuracy: %1.5f\n",
                   1 - tree.min_objective() + c*r_list.size());
            print_final_rulelist(r_list, tree.opt_predictions(),
                                 latex_out, rules, labels, opt_fname);
            logger.dumpState();
        } else {
            printf("BBOUND_STOCHASTIC No Permutation Map \n");
            CacheTree<BaseNode> tree(nsamples, nrules, c, rules, labels);
            bbound_stochastic<BaseNode,
                              CapturedPermutationMap>(&tree, max_num_nodes,
                                                      &base_construct_policy,
                                                      &captured_permutation_insert,
                                                      &captured_map_garbage_collect,
                                                      NULL);
            printf("final num_nodes: %zu\n", tree.num_nodes());
            printf("final num_evaluated: %zu\n", tree.num_evaluated());
            printf("final min_objective: %1.5f\n", tree.min_objective());
            const std::vector<size_t>& r_list = tree.opt_rulelist();
            printf("final accuracy: %1.5f\n",
                   1 - tree.min_objective() + c*r_list.size());
            print_final_rulelist(r_list, tree.opt_predictions(),
                                 latex_out, rules, labels, opt_fname);

            logger.dumpState();
        }
    }

    if (run_bfs) {
        if (use_prefix_perm_map) {
            printf("BFS Permutation Map\n");        
            CacheTree<BaseNode> tree(nsamples, nrules, c, rules, labels);
            BaseQueue bfs_q;
            PrefixPermutationMap p;
            bbound_queue<BaseNode,
                         BaseQueue,
                         PrefixPermutationMap>(&tree,
                                               max_num_nodes,
                                               &base_construct_policy,
                                               &bfs_q,
                                               &base_queue_front,
                                               &prefix_permutation_insert,
                                               &bfs_prefix_map_garbage_collect,
                                               &p);

            printf("final num_nodes: %zu\n", tree.num_nodes());
            printf("final num_evaluated: %zu\n", tree.num_evaluated());
            printf("final min_objective: %1.5f\n", tree.min_objective());
            const std::vector<size_t>& r_list = tree.opt_rulelist();
            printf("final accuracy: %1.5f\n",
                   1 - tree.min_objective() + c*r_list.size());
            print_final_rulelist(r_list, tree.opt_predictions(),
                                 latex_out, rules, labels, opt_fname);
        } else if (use_captured_sym_map) {
            printf("BFS Captured Symmetry Map\n");        
            CacheTree<BaseNode> tree(nsamples, nrules, c, rules, labels);
            BaseQueue bfs_q;
            CapturedPermutationMap p;
            bbound_queue<BaseNode,
                         BaseQueue,
                         CapturedPermutationMap>(&tree,
                                                 max_num_nodes,
                                                 &base_construct_policy,
                                                 &bfs_q,
                                                 &base_queue_front,
                                                 &captured_permutation_insert,
                                                 &captured_map_garbage_collect,
                                                 &p);

            printf("final num_nodes: %zu\n", tree.num_nodes());
            printf("final num_evaluated: %zu\n", tree.num_evaluated());
            printf("final min_objective: %1.5f\n", tree.min_objective());
            const std::vector<size_t>& r_list = tree.opt_rulelist();
            printf("final accuracy: %1.5f\n",
                   1 - tree.min_objective() + c*r_list.size());
            print_final_rulelist(r_list, tree.opt_predictions(),
                                 latex_out, rules, labels, opt_fname);
        }
        else {
            printf("BFS No Permutation Map\n");        
            CacheTree<BaseNode> tree(nsamples, nrules, c, rules, labels);
            BaseQueue bfs_q;
            bbound_queue<BaseNode,
                         BaseQueue,
                         CapturedPermutationMap>(&tree,
                                                 max_num_nodes,
                                                 &base_construct_policy,
                                                 &bfs_q,
                                                 &base_queue_front,
                                                 &captured_permutation_insert,
                                                 &captured_map_garbage_collect,
                                                 NULL);

            printf("final num_nodes: %zu\n", tree.num_nodes());
            printf("final num_evaluated: %zu\n", tree.num_evaluated());
            printf("final min_objective: %1.5f\n", tree.min_objective());
            const std::vector<size_t>& r_list = tree.opt_rulelist();
            printf("final accuracy: %1.5f\n",
                   1 - tree.min_objective() + c*r_list.size());
            print_final_rulelist(r_list, tree.opt_predictions(),
                                 latex_out, rules, labels, opt_fname);
        }
    }

    if (run_curiosity) {
        if (use_curious_cmp) {
            if (use_prefix_perm_map) {
                printf("CURIOSITY Prefix Permutation Map\n");
                CacheTree<CuriousNode> tree(nsamples, nrules, c, rules, labels);
                CuriousQueue curious_q(curious_cmp);
                PrefixPermutationMap p;
                bbound_queue<CuriousNode,
                             CuriousQueue,
                             PrefixPermutationMap>(&tree,
                                                   max_num_nodes,
                                                   &curious_construct_policy,
                                                   &curious_q,
                                                   &curious_queue_front,
                                                   &prefix_permutation_insert,
                                                   &prefix_map_garbage_collect,
                                                   &p);
                printf("final num_nodes: %zu\n", tree.num_nodes());
                printf("final num_evaluated: %zu\n", tree.num_evaluated());
                printf("final min_objective: %1.5f\n", tree.min_objective());
                const std::vector<size_t>& r_list = tree.opt_rulelist();
                printf("final accuracy: %1.5f\n",
                       1 - tree.min_objective() + c*r_list.size());
                print_final_rulelist(r_list, tree.opt_predictions(),
                                     latex_out, rules, labels, opt_fname);
            } else if (use_captured_sym_map) {
                printf("CURIOSITY Captured Symmetry Map\n");
                CacheTree<CuriousNode> tree(nsamples, nrules, c, rules, labels);
                CuriousQueue curious_q(curious_cmp);
                CapturedPermutationMap p;
                bbound_queue<CuriousNode,
                             CuriousQueue,
                             CapturedPermutationMap>(&tree,
                                                     max_num_nodes,
                                                     &curious_construct_policy,
                                                     &curious_q,
                                                     &curious_queue_front,
                                                     &captured_permutation_insert,
                                                     &captured_map_garbage_collect,
                                                     &p);
                printf("final num_nodes: %zu\n", tree.num_nodes());
                printf("final num_evaluated: %zu\n", tree.num_evaluated());
                printf("final min_objective: %1.5f\n", tree.min_objective());
                const std::vector<size_t>& r_list = tree.opt_rulelist();
                printf("final accuracy: %1.5f\n",
                       1 - tree.min_objective() + c*r_list.size());
                print_final_rulelist(r_list, tree.opt_predictions(),
                                     latex_out, rules, labels, opt_fname);
            }
            else {
                printf("CURIOSITY No Permutation Map\n");
                CacheTree<CuriousNode> tree(nsamples, nrules, c, rules, labels);
                CuriousQueue curious_q(curious_cmp);
                bbound_queue<CuriousNode,
                             CuriousQueue,
                             CapturedPermutationMap>(&tree,
                                                     max_num_nodes,
                                                     &curious_construct_policy,
                                                     &curious_q,
                                                     &curious_queue_front,
                                                     &captured_permutation_insert,
                                                     &captured_map_garbage_collect,
                                                     NULL);
                printf("final num_nodes: %zu\n", tree.num_nodes());
                printf("final num_evaluated: %zu\n", tree.num_evaluated());
                printf("final min_objective: %1.5f\n", tree.min_objective());
                const std::vector<size_t>& r_list = tree.opt_rulelist();
                printf("final accuracy: %1.5f\n",
                       1 - tree.min_objective() + c*r_list.size());
                print_final_rulelist(r_list, tree.opt_predictions(),
                                     latex_out, rules, labels, opt_fname);
            }
        } else if (use_lower_bound_cmp) {
            if (use_prefix_perm_map) {
                printf("CURIOUS LOWER BOUND Prefix Permutation Map\n");
                CacheTree<CuriousNode> tree(nsamples, nrules, c, rules, labels);
                CuriousQueue curious_q(lower_bound_cmp);
                PrefixPermutationMap p;
                bbound_queue<CuriousNode,
                             CuriousQueue,
                             PrefixPermutationMap>(&tree,
                                                   max_num_nodes,
                                                   &curious_construct_policy,
                                                   &curious_q,
                                                   &curious_queue_front,
                                                   &prefix_permutation_insert,
                                                   &prefix_map_garbage_collect,
                                                   &p);
                printf("final num_nodes: %zu\n", tree.num_nodes());
                printf("final num_evaluated: %zu\n", tree.num_evaluated());
                printf("final min_objective: %1.5f\n", tree.min_objective());
                const std::vector<size_t>& r_list = tree.opt_rulelist();
                printf("final accuracy: %1.5f\n",
                       1 - tree.min_objective() + c*r_list.size());
                print_final_rulelist(r_list, tree.opt_predictions(),
                                     latex_out, rules, labels, opt_fname);
            } else if (use_captured_sym_map) {
                printf("CURIOUS LOWER BOUND Captured Symmetry Map\n");
                CacheTree<CuriousNode> tree(nsamples, nrules, c, rules, labels);
                CuriousQueue curious_q(lower_bound_cmp);
                CapturedPermutationMap p;
                bbound_queue<CuriousNode,
                             CuriousQueue,
                             CapturedPermutationMap>(&tree,
                                                     max_num_nodes,
                                                     &curious_construct_policy,
                                                     &curious_q,
                                                     &curious_queue_front,
                                                     &captured_permutation_insert,
                                                     &captured_map_garbage_collect,
                                                     &p);
                printf("final num_nodes: %zu\n", tree.num_nodes());
                printf("final num_evaluated: %zu\n", tree.num_evaluated());
                printf("final min_objective: %1.5f\n", tree.min_objective());
                const std::vector<size_t>& r_list = tree.opt_rulelist();
                printf("final accuracy: %1.5f\n",
                       1 - tree.min_objective() + c*r_list.size());
                print_final_rulelist(r_list, tree.opt_predictions(),
                                     latex_out, rules, labels, opt_fname);
            }
            else {
                printf("CURIOUS LOWER BOUND No Permutation Map\n");
                CacheTree<CuriousNode> tree(nsamples, nrules, c, rules, labels);
                CuriousQueue curious_q(lower_bound_cmp);
                bbound_queue<CuriousNode,
                             CuriousQueue,
                             CapturedPermutationMap>(&tree,
                                                     max_num_nodes,
                                                     &curious_construct_policy,
                                                     &curious_q,
                                                     &curious_queue_front,
                                                     &captured_permutation_insert,
                                                     &captured_map_garbage_collect,
                                                     NULL);
                printf("final num_nodes: %zu\n", tree.num_nodes());
                printf("final num_evaluated: %zu\n", tree.num_evaluated());
                printf("final min_objective: %1.5f\n", tree.min_objective());
                const std::vector<size_t>& r_list = tree.opt_rulelist();
                printf("final accuracy: %1.5f\n",
                       1 - tree.min_objective() + c*r_list.size());
                print_final_rulelist(r_list, tree.opt_predictions(),
                                     latex_out, rules, labels, opt_fname);
            }
        }
    }

    printf("final total time: %f\n", time_diff(init));
    logger.dumpState();
    logger.closeFile();
    printf("\ndelete rules\n");
	rules_free(rules, nrules, 1);
	printf("delete labels\n");
	rules_free(labels, nlabels, 0);
	printf("tree destructors\n");
}
