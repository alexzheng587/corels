#include "cache.hh"

template<class T>
using construct_signature = Node<T>* (*)(size_t, size_t, bool, bool, double, double, Node<T>* parent, int, int, int, double);

BaseNode* base_construct_policy(size_t new_rule, size_t nrules, bool prediction, bool default_prediction,
                                double lower_bound, double objective, BaseNode* parent, int num_not_captured,
                                int nsamples, int len_prefix, double c);

CuriousNode* curious_construct_policy(size_t new_rule, size_t nrules, bool prediction, bool default_prediction,
                                      double lower_bound, double objective, CuriousNode* parent, int num_not_captured,
                                      int nsamples, int len_prefix, double c);

template<class T>
extern void evaluate_children(CacheTree<T>* tree, Node<T>* parent, VECTOR parent_not_captured, construct_signature<T> construct_policy);

template<class T>
extern Node<T>* stochastic_select(CacheTree<T>* tree, VECTOR not_captured);

template<class T>
extern void bbound_stochastic(CacheTree<T>* tree, size_t max_num_nodes, construct_signature<T> construct_policy);
