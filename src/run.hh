#pragma once

#ifdef __cplusplus
#include "rule.h"
#include "queue.hh"
extern "C" {
#endif

double run_corels (double c, char* vstring, int curiosity_policy,
                   int map_type, int ablation, int calculate_size, int nrules, int nlabels,
                   int nsamples, rule_t* rules, rule_t* labels, rule_t* meta, int freq, char* log_fname,
                   PermutationMap*& pmap, CacheTree*& tree, Queue*& queue, double& init,
                   int verbosity, int num_threads, int max_num_nodes, int nmeta, int random_seed,
                   std::vector<int>* rulelist, std::vector<int>* classes);

#ifdef __cplusplus
}
#endif