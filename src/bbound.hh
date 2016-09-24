#include "cache.hh"

struct time {
    double total_time;
    double evaluate_children_time;
    int evaluate_children_num;
    double stochastic_select_time;
    int stochastic_select_num;
    double rule_evaluation_time;
    int rule_evaluation_num;
    double lower_bound_time;
    int lower_bound_num;
    double objective_time;
    int objective_num;
    double tree_insertion_time;
    int tree_insertion_num;
};

void evaluate_children(CacheTree* tree, CacheNode* parent, VECTOR parent_not_captured);

CacheNode* stochastic_select(CacheTree* tree, VECTOR not_captured);

struct time* bbound_stochastic(CacheTree* tree, size_t max_num_nodes);

void bbound_greedy(size_t nsamples, size_t nrules, rule_t *rules, rule_t *labels, size_t max_prefix_length);
