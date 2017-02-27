#include "queue.hh"
#include "utils.hh"

#include <iostream>
#include <stdio.h>
#include <queue>
#include <getopt.h>

Logger logger;
int ablation = 0;

int main(int argc, char *argv[]) {
    const char usage[] = "USAGE: %s [-s] [-b] "
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
    /* only parsing happens here */
    while ((ch = getopt(argc, argv, "sbLc:p:v:n:r:f:w:a:")) != -1) {
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

    int nrules, nsamples, nlabels, nsamples_chk;
    rule_t *rules, *labels;
    rules_init(argv[0], &nrules, &nsamples, &rules, 1);
    rules_init(argv[1], &nlabels, &nsamples_chk, &labels, 0);

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
    sprintf(froot, "../logs/for-%s-%s%s%s-%s-%s-removed=%s-max_num_nodes=%d-c=%.7f-v=%d-f=%d",
            pch ? pch + 1 : "",
            run_stochastic ? "stochastic" : "",
            run_bfs ? "bfs" : "",
            run_curiosity ? curiosity_map[curiosity_policy].c_str() : "",
            run_pmap ? (use_prefix_perm_map ? "with_prefix_perm_map" : "with_captured_symmetry_map") : "no_pmap",
            meta ? "minor" : "no_minor",
            ablation ? ((ablation == 1) ? "support" : "lookahead") : "none",
            max_num_nodes, c, verbosity, freq);
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
    double init = timestamp();
    if (run_stochastic) {
        PermutationMap<BaseNode>* p;
        if (use_prefix_perm_map) {
            printf("BBOUND_STOCHASTIC Permutation Map\n");
            PrefixPermutationMap<BaseNode>* prefix_pmap = new PrefixPermutationMap<BaseNode>;
            p = (PermutationMap<BaseNode>*) prefix_pmap;
        } else if (use_captured_sym_map) {
            printf("BBOUND_STOCHASTIC Captured Symmetry Map\n");
            CapturedPermutationMap<BaseNode>* cap_pmap = new CapturedPermutationMap<BaseNode>;
            p = (PermutationMap<BaseNode>*) cap_pmap;
        } else {
            printf("BBOUND_STOCHASTIC No Permutation Map \n");
            NullPermutationMap<BaseNode>* null_pmap = new NullPermutationMap<BaseNode>;
            p = (PermutationMap<BaseNode>*) null_pmap;
        }
        CacheTree<BaseNode> tree(nsamples, nrules, c, rules, labels, meta);
        bbound_stochastic<BaseNode>(&tree, max_num_nodes,
                                    &base_construct_policy,
                                    p);
        printf("final num_nodes: %zu\n", tree.num_nodes());
        printf("final num_evaluated: %zu\n", tree.num_evaluated());
        printf("final min_objective: %1.5f\n", tree.min_objective());
        const std::vector<unsigned short>& r_list = tree.opt_rulelist();
        printf("final accuracy: %1.5f\n",
               1 - tree.min_objective() + c*r_list.size());
        print_final_rulelist(r_list, tree.opt_predictions(),
                             latex_out, rules, labels, opt_fname);

        logger.dumpState();
        delete p;
    } else if (run_bfs) {
        PermutationMap<BaseNode>* p;
        if (use_prefix_perm_map) {
            printf("BFS Permutation Map\n");
            PrefixPermutationMap<BaseNode>* prefix_pmap = new PrefixPermutationMap<BaseNode>;
            p = (PermutationMap<BaseNode>*) prefix_pmap;
        } else if (use_captured_sym_map) {
            printf("BFS Captured Symmetry Map\n");
            CapturedPermutationMap<BaseNode>* cap_pmap = new CapturedPermutationMap<BaseNode>;
            p = (PermutationMap<BaseNode>*) cap_pmap;
        } else {
            printf("BFS No Permutation Map \n");
            NullPermutationMap<BaseNode>* null_pmap = new NullPermutationMap<BaseNode>;
            p = (PermutationMap<BaseNode>*) null_pmap;
        }
        CacheTree<BaseNode> tree(nsamples, nrules, c, rules, labels, meta);
        Queue<BaseNode> bfs_q;
        bbound_queue<BaseNode>(&tree,
                               max_num_nodes,
                               &base_construct_policy,
                               &bfs_q,
                               p, 0, 0);

        printf("final num_nodes: %zu\n", tree.num_nodes());
        printf("final num_evaluated: %zu\n", tree.num_evaluated());
        printf("final min_objective: %1.5f\n", tree.min_objective());
        const std::vector<unsigned short>& r_list = tree.opt_rulelist();
        printf("final accuracy: %1.5f\n",
               1 - tree.min_objective() + c*r_list.size());
        print_final_rulelist(r_list, tree.opt_predictions(),
                             latex_out, rules, labels, opt_fname);
        delete p;
    } else if (run_curiosity) {
        PermutationMap<CuriousNode>* p;
        if (use_prefix_perm_map) {
            PrefixPermutationMap<CuriousNode>* prefix_pmap = new PrefixPermutationMap<CuriousNode>;
            p = (PermutationMap<CuriousNode>*) prefix_pmap;
        } else if (use_captured_sym_map) {
            CapturedPermutationMap<CuriousNode>* cap_pmap = new CapturedPermutationMap<CuriousNode>;
            p = (PermutationMap<CuriousNode>*) cap_pmap;
        } else {
            NullPermutationMap<CuriousNode>* null_pmap = new NullPermutationMap<CuriousNode>;
            p = (PermutationMap<CuriousNode>*) null_pmap;
        }
        if (curiosity_policy == 1) {
            if (use_prefix_perm_map) {
                printf("Curiosity Prefix Permutation Map\n");
            } else if (use_captured_sym_map) {
                printf("Curiosity Captured Symmetry Map\n");
            } else {
                printf("Curiosity No Permutation Map \n");
            }
            CacheTree<CuriousNode> *tree = new CacheTree<CuriousNode>(nsamples, nrules, c, rules, labels, meta);
            size_t num_iter;
            CuriousQueue<CuriousNode> curious_q;
            num_iter = bbound_queue<CuriousNode>(tree, max_num_nodes,
                                                   &curious_construct_policy,
                                                   &curious_q,
                                                   p, 0, switch_iter);
            if (curious_q.size() > 0) {
                printf("Switching to curious lower bound policy... \n");
                LowerBoundQueue<CuriousNode> curious_lb_q;
                while(!curious_q.empty()) {
                    CuriousNode* selected_node = curious_q.front();
                    curious_q.pop();
                    curious_lb_q.push(selected_node);
                }
                num_iter = bbound_queue<CuriousNode>(tree, max_num_nodes,
                                                       &curious_construct_policy,
                                                       &curious_lb_q,
                                                       p, num_iter, 0);
            }
            printf("final num_nodes: %zu\n", tree->num_nodes());
            printf("final num_evaluated: %zu\n", tree->num_evaluated());
            printf("final min_objective: %1.5f\n", tree->min_objective());
            const std::vector<unsigned short>& r_list = tree->opt_rulelist();
            printf("final accuracy: %1.5f\n",
                   1 - tree->min_objective() + c*r_list.size());
            print_final_rulelist(r_list, tree->opt_predictions(),
                                 latex_out, rules, labels, opt_fname);
        } else if (curiosity_policy == 2) {
            if (use_prefix_perm_map) {
                printf("CURIOUS LOWER BOUND Prefix Permutation Map\n");
            } else if (use_captured_sym_map) {
                printf("CURIOUS LOWER BOUND Captured Symmetry Map\n");
            } else {
                printf("CURIOUS LOWER BOUND No Permutation Map\n");
            }
            CacheTree<CuriousNode> tree(nsamples, nrules, c, rules, labels, meta);
            LowerBoundQueue<CuriousNode> curious_q;
            bbound_queue<CuriousNode>(&tree,
                                         max_num_nodes,
                                         &curious_construct_policy,
                                         &curious_q,
                                         p, 0, 0);
            printf("final num_nodes: %zu\n", tree.num_nodes());
            printf("final num_evaluated: %zu\n", tree.num_evaluated());
            printf("final min_objective: %1.5f\n", tree.min_objective());
            const std::vector<unsigned short>& r_list = tree.opt_rulelist();
            printf("final accuracy: %1.5f\n",
                   1 - tree.min_objective() + c*r_list.size());
            print_final_rulelist(r_list, tree.opt_predictions(),
                                 latex_out, rules, labels, opt_fname);
        } else if (curiosity_policy == 3) {
            if (use_prefix_perm_map) {
                printf("CURIOUS OBJECTIVE Prefix Permutation Map\n");
            } else if (use_captured_sym_map) {
                printf("CURIOUS OBJECTIVE Captured Symmetry Map\n");
            } else {
                printf("CURIOUS OBJECTIVE No Permutation Map\n");
            }
                CacheTree<CuriousNode> tree(nsamples, nrules, c, rules, labels, meta);
                ObjectiveQueue<CuriousNode> curious_q;
                bbound_queue<CuriousNode>(&tree,
                                           max_num_nodes,
                                           &curious_construct_policy,
                                           &curious_q,
                                           p, 0, 0);
                printf("final num_nodes: %zu\n", tree.num_nodes());
                printf("final num_evaluated: %zu\n", tree.num_evaluated());
                printf("final min_objective: %1.5f\n", tree.min_objective());
                const std::vector<unsigned short>& r_list = tree.opt_rulelist();
                printf("final accuracy: %1.5f\n",
                       1 - tree.min_objective() + c*r_list.size());
                print_final_rulelist(r_list, tree.opt_predictions(),
                                     latex_out, rules, labels, opt_fname);
        } else if (curiosity_policy == 4) {
            if (use_prefix_perm_map) {
                printf("DFS Prefix Permutation Map\n");
            } else if (use_captured_sym_map) {
                printf("DFS Captured Symmetry Map\n");
            } else {
                printf("DFS No Permutation Map\n");
            }
                CacheTree<CuriousNode> tree(nsamples, nrules, c, rules, labels, meta);
                DFSQueue<CuriousNode> curious_q;
                bbound_queue<CuriousNode>(&tree,
                                           max_num_nodes,
                                           &curious_construct_policy,
                                           &curious_q,
                                           p, 0, 0);
                printf("final num_nodes: %zu\n", tree.num_nodes());
                printf("final num_evaluated: %zu\n", tree.num_evaluated());
                printf("final min_objective: %1.5f\n", tree.min_objective());
                const std::vector<unsigned short>& r_list = tree.opt_rulelist();
                printf("final accuracy: %1.5f\n",
                       1 - tree.min_objective() + c*r_list.size());
                print_final_rulelist(r_list, tree.opt_predictions(),
                                     latex_out, rules, labels, opt_fname);
        }
        delete p;
    }

    printf("final total time: %f\n", time_diff(init));
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
