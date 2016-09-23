#include "bbound.hh"
#include <iostream>
#include <stdio.h>

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
    CacheTree tree(nsamples, nrules, c, rules, labels);
    bbound_stochastic(&tree, 100000000);
    printf("\nnum_nodes: %zu\n", tree.num_nodes());
    printf("num_evaluated: %zu\n", tree.num_evaluated());
    printf("\nmin_objective: %1.5f\n", tree.min_objective());

    printf("\ndelete rules\n");
	rules_free(rules, nrules, 1);
	printf("delete labels\n");
	rules_free(labels, nlabels, 0);
    printf("deconstructor\n");
}
