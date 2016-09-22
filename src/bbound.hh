#include "cache.hh"


void evaluate_children(CacheTree* tree, CacheNode* parent, VECTOR parent_not_captured);

CacheNode* stochastic_select(CacheTree* tree, VECTOR not_captured);

void bbound_stochastic(CacheTree* tree, size_t max_num_nodes);
