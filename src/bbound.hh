#include "cache.hh"


CacheNode* stochastic_select(CacheTree* tree, VECTOR not_captured);

void bbound_stochastic(CacheTree* tree, size_t max_num_nodes);
