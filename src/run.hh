#pragma once

#include "rule.h"
#include "queue.hh"
#ifdef __cplusplus
extern "C" {
#endif

double run_corels (double c, char* vstring, char* loss_type_str, int curiosity_policy,
                   int map_type, int ablation, int calculate_size, int latex_out, int nrules, int nlabels,
                   int nsamples, rule_t* rules, rule_t* labels, rule_t* meta, minority_class_t* minority_class, int freq, char* log_fname,
                   PermutationMap*& pmap, CacheTree*& tree, Queue*& queue, double& init,
                   int verbosity, int num_threads, int max_num_nodes, int nmeta, int random_seed,
                   std::vector<int>* rulelist, std::vector<int>* classes, double weight);

#ifdef __cplusplus
}
#endif
