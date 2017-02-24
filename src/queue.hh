#pragma once
#include "pmap.hh"
#include <functional>
#include <queue>
#include <set>

/*
 * Queue
 */

template<class N>
class NullQueue {
  public:
    void push(N*) {};
    size_t size() {return 0;};
};

typedef std::queue<BaseNode*> BaseQueue;

// lambda function for priority queue metric using curiosity
auto curious_cmp = [](CuriousNode* left, CuriousNode* right) {
    return left->get_storage() > right->get_storage();
};

// lambda function for priority queue metric using lower bound as curiosity
auto lower_bound_cmp = [](CuriousNode* left, CuriousNode* right) {
    return left->lower_bound() > right->lower_bound();
};

// lambda function for priority queue metric using objective as curiosity
auto objective_cmp = [](CuriousNode* left, CuriousNode* right) {
    return left->objective() > right->objective();
};

// lambda function for priority queue metric implementing depth-first
auto depth_first_cmp = [](CuriousNode* left, CuriousNode* right) {
    return left->depth() < right->depth();
};

typedef std::priority_queue<CuriousNode*, std::vector<CuriousNode*>,
                            std::function<bool(CuriousNode*, CuriousNode*)> > CuriousQueue;

BaseNode* base_queue_front(BaseQueue* q);

CuriousNode* curious_queue_front(CuriousQueue* q);

template<class N>
extern N* stochastic_select(CacheTree<N>* tree, VECTOR not_captured);

template<class N, class P>
extern void bbound_stochastic(CacheTree<N>* tree,
                              size_t max_num_nodes,
                              construct_signature<N> construct_policy,
                              permutation_insert_signature<N, P> permutation_insert,
                              pmap_garbage_collect_signature<P> pmap_garbage_collect,
                              P* p);

template<class N, class Q>
extern N*
queue_select(CacheTree<N>* tree, Q* q, N*(*front)(Q*), VECTOR captured);

template<class N, class Q, class P>
extern int bbound_queue(CacheTree<N>* tree,
                         size_t max_num_nodes,
                         construct_signature<N> construct_policy,
                         Q* q, N*(*front)(Q*),
                         permutation_insert_signature<N, P> permutation_insert,
                         pmap_garbage_collect_signature<P> pmap_garbage_collect,
                         P* p, size_t num_iter, size_t switch_iter);

template<class N, class Q, class P>
extern void evaluate_children(CacheTree<N>* tree, N* parent,
                              VECTOR parent_not_captured,
                              construct_signature<N> construct_policy,
                              Q* q,
                              permutation_insert_signature<N, P> permutation_insert,
                              P* p);

void bbound_greedy(size_t nsamples, size_t nrules, rule_t *rules, rule_t *labels, size_t max_prefix_length);
