#include "bbound.hh"
#include "utils.hh"

#include <iostream>
#include <stdio.h>
#include <queue>
#include <getopt.h>

Logger logger;

int main(int argc, char *argv[]) {
    const char usage[] = "USAGE: %s [-s] [-b] [-c] [-p] [-v verbosity]\n";
    int nrules, nsamples, nlabels, nsamples_chk;
    rule_t *rules, *labels;

    const char infile[] = "../data/tdata_R.out";
    rules_init(infile, &nrules, &nsamples, &rules, 1);

    const char lfile[] = "../data/tdata_R.label";
    rules_init(lfile, &nlabels, &nsamples_chk, &labels, 0);

    extern char *optarg;
    bool run_stochastic = false;
    bool run_bfs = false;
    bool run_curiosity = false;
    bool use_perm_map = false;
    int verbosity = 0;
    char ch;
    bool error = false;
    /* only parsing happens here */
    while ((ch = getopt(argc, argv, "sbcpv:")) != -1) {
        switch (ch) {
        case 's':
            run_stochastic = true;
            break;
        case 'b':
            run_bfs = true;
            break;
        case 'c':
            run_curiosity = true;
            break;
        case 'p':
            use_perm_map = true;
            break;
        case 'v':
            verbosity = atoi(optarg);
            break;
        default:
            error = true;
        }
    }
    if (!(run_stochastic || run_bfs || run_curiosity || use_perm_map)) {
        run_stochastic = run_bfs = run_curiosity = use_perm_map = true;
    }
    /** handle errors in parsing */
    if (error) {
        fprintf(stderr, usage, argv[0]);
        exit(1);
    }

    if (verbosity >= 100) {
        printf("\n%d rules %d samples\n\n", nrules, nsamples);
        rule_print_all(rules, nrules, nsamples);

        printf("\nLabels (%d) for %d samples\n\n", nlabels, nsamples);
        rule_print_all(labels, nlabels, nsamples);
    }
    if (verbosity >= 10) {
        printf("\nGreedy algorithm:\n\n");
        bbound_greedy(nsamples, nrules, rules, labels, 8);
    }

    double c = 0.001;
    //logger() = new Logger();
    logger.setVerbosity(verbosity); /** need to define verbosity semantics **/
    if (run_stochastic) {
        printf("BBOUND_STOCHASTIC\n");
        CacheTree<BaseNode> tree(nsamples, nrules, c, rules, labels);
        bbound_stochastic<BaseNode>(&tree, 100000,
                                    &base_construct_policy);
        printf("\nnum_nodes: %zu\n", tree.num_nodes());
        printf("num_evaluated: %zu\n", tree.num_evaluated());
        printf("\nmin_objective: %1.5f\n", tree.min_objective());
        logger.dumpState();
        logger.clearState();
    }

    if (run_bfs) {
        if (use_perm_map) {
            printf("\n\n\nBFS Permutation Map\n");        
            CacheTree<BaseNode> tree(nsamples, nrules, c, rules, labels);
            BaseQueue bfs_q;
            PrefixPermutationMap<BaseNode> p;
            bbound_queue<BaseNode,
                         BaseQueue,
                         PrefixPermutationMap<BaseNode> >(&tree,
                                                          100000,
                                                          &base_construct_policy,
                                                          &bfs_q,
                                                          &base_queue_front,
                                                          &p);

            printf("\nnum_nodes: %zu\n", tree.num_nodes());
            printf("num_evaluated: %zu\n", tree.num_evaluated());
            printf("\nmin_objective: %1.5f\n", tree.min_objective());
            logger.dumpState();
            logger.clearState();
        } else {
            printf("\n\n\nBFS No Permutation Map\n");
            CacheTree<BaseNode> tree(nsamples, nrules, c, rules, labels);
            BaseQueue bfs_q_2;
            NullPermutationMap<BaseNode> p2;
            bbound_queue<BaseNode,
                         BaseQueue,
                         NullPermutationMap<BaseNode> >(&tree, 100000,
                                                        &base_construct_policy,
                                                        &bfs_q_2,
                                                        &base_queue_front,
                                                        NULL);
            printf("\nnum_nodes: %zu\n", tree.num_nodes());
            printf("num_evaluated: %zu\n", tree.num_evaluated());
            printf("\nmin_objective: %1.5f\n", tree.min_objective());
            logger.dumpState();
            logger.clearState();
        }
    }

    if (run_curiosity) {
        printf("\n\n\nCURIOUSITY\n");
        CacheTree<CuriousNode> tree(nsamples, nrules, c, rules, labels);
        CuriousQueue curious_q(curious_cmp);
        PrefixPermutationMap<CuriousNode> p3;
        bbound_queue<CuriousNode,
                     CuriousQueue,
                     PrefixPermutationMap<CuriousNode> >(&tree, 100000,
                                                         &curious_construct_policy,
                                                         &curious_q,
                                                         &curious_queue_front,
                                                         &p3);
        printf("\nnum_nodes: %zu\n", tree.num_nodes());
        printf("num_evaluated: %zu\n", tree.num_evaluated());
        printf("\nmin_objective: %1.5f\n", tree.min_objective());
        logger.dumpState();
        logger.clearState();
    }

    printf("\ndelete rules\n");
	rules_free(rules, nrules, 1);
	printf("delete labels\n");
	rules_free(labels, nlabels, 0);
	printf("tree destructors\n");
}
