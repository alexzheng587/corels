#include "bbound.hh"
#include "time.hh"
#include <iostream>
#include <stdio.h>
#include <queue>

int main()
{
    int nrules, nsamples, nlabels, nsamples_chk;
    rule_t *rules, *labels;

    char infile_name[] = "../data/tdata_R.out";
    char* infile = infile_name;
    rules_init(infile, &nrules, &nsamples, &rules, 1);
    printf("\n%d rules %d samples\n\n", nrules, nsamples);
    rule_print_all(rules, nrules, nsamples);

    char lfile_name[] = "../data/tdata_R.label";
    char* lfile = lfile_name;
    rules_init(lfile, &nlabels, &nsamples_chk, &labels, 0);
    printf("\nLabels (%d) for %d samples\n\n", nlabels, nsamples);
    rule_print_all(labels, nlabels, nsamples);

    double c = 0.01;
    CacheTree<BaseNode> tree(nsamples, nrules, c, rules, labels);
    printf("\nGreedy algorithm:\n\n");
    bbound_greedy(nsamples, nrules, rules, labels, 8);

    struct time* times = (struct time*) calloc(1, sizeof(*times));
    bbound_stochastic<BaseNode>(&tree, 1000000,
                                &base_construct_policy,
                                times);
    printf("BBOUND_STOCHASTIC\n");
    printf("\nnum_nodes: %zu\n", tree.num_nodes());
    printf("num_evaluated: %zu\n", tree.num_evaluated());
    printf("\nmin_objective: %1.5f\n", tree.min_objective());
    printf("Total time: %f\n", times->total_time);
    printf("Evaluate children time: %f\n", times->evaluate_children_time);
    printf("Stochastic select time: %f\n", times->stochastic_select_time);
    printf("Rule evaluation time: %f\n", times->rule_evaluation_time);
    printf("Lower bound time: %f\n", times->lower_bound_time);
    printf("Number of lower bound evaluations: %d\n", times->lower_bound_num);
    printf("Objective time: %f\n", times->objective_time);
    printf("Total tree insertion time: %f\n", times->tree_insertion_time);
    printf("Number of tree insertions: %i\n", times->tree_insertion_num);
    clear_time(times);

    CacheTree<BaseNode> tree2(nsamples, nrules, c, rules, labels);
    BaseQueue bfs_q;
    bbound_queue<BaseNode, BaseQueue>(&tree2, 1000000,
                                      &base_construct_policy,
                                      &bfs_q,
                                      times);
    printf("\n\n\nBBOUND_BFS\n");
    printf("\nnum_nodes: %zu\n", tree2.num_nodes());
    printf("num_evaluated: %zu\n", tree2.num_evaluated());
    printf("\nmin_objective: %1.5f\n", tree2.min_objective());
    printf("Total time: %f\n", times->total_time);
    printf("Evaluate children time: %f\n", times->evaluate_children_time);
    printf("Stochastic select time: %f\n", times->stochastic_select_time);
    printf("Rule evaluation time: %f\n", times->rule_evaluation_time);
    printf("Lower bound time: %f\n", times->lower_bound_time);
    printf("Number of lower bound evaluations: %d\n", times->lower_bound_num);
    printf("Objective time: %f\n", times->objective_time);
    printf("Total tree insertion time: %f\n", times->tree_insertion_time);
    printf("Number of tree insertions: %i\n", times->tree_insertion_num);
    free(times);

    printf("\ndelete rules\n");
	rules_free(rules, nrules, 1);
	printf("delete labels\n");
	rules_free(labels, nlabels, 0);
    printf("deconstructor\n");
}
