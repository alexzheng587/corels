#include "cache.hh"
#include <set>

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

template<class N>
using construct_signature = N* (*)(size_t, size_t, bool, bool, double, double, N* parent, int, int, int, double);

BaseNode* base_construct_policy(size_t new_rule, size_t nrules, bool prediction, bool default_prediction,
                                double lower_bound, double objective, BaseNode* parent, int num_not_captured,
                                int nsamples, int len_prefix, double c);

CuriousNode* curious_construct_policy(size_t new_rule, size_t nrules, bool prediction, bool default_prediction,
                                      double lower_bound, double objective, CuriousNode* parent, int num_not_captured,
                                      int nsamples, int len_prefix, double c);

template<class N>
extern void evaluate_children(CacheTree<N>* tree, N* parent, VECTOR parent_not_captured, std::set<size_t> ordered_parent, construct_signature<N> construct_policy);

template<class N>
extern std::pair<N*, std::set<size_t> > stochastic_select(CacheTree<N>* tree, VECTOR not_captured);

template<class N>
extern struct time* bbound_stochastic(CacheTree<N>* tree, size_t max_num_nodes, construct_signature<N> construct_policy);

void bbound_greedy(size_t nsamples, size_t nrules, rule_t *rules, rule_t *labels, size_t max_prefix_length);
