#include "cache.hh"

template<class N>
extern void evaluate_children(CacheTree<N>* tree, N* parent, VECTOR parent_not_captured);

template<class N>
extern N* stochastic_select(CacheTree<N>* tree, VECTOR not_captured);

template<class N>
extern void bbound_stochastic(CacheTree<N>* tree, size_t max_num_nodes);
