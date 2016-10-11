#ifndef _BBOUND_H_
#define _BBOUND_H_

#include "cache.hh"
#include <functional>
#include <queue>
#include <map>
#include <set>

/*
 * Queue
 */

template<class N>
class NullQueue {
  public:
    void push(N* node) {
        (void) node;
    };
};

typedef std::queue<BaseNode*> BaseQueue;

// lambda function for priority queue metric using curiosity
auto curious_cmp = [](CuriousNode* left, CuriousNode* right) {
    return left->get_storage() > right->get_storage();
};

typedef std::priority_queue<CuriousNode*, std::vector<CuriousNode*>,
                            std::function<bool(CuriousNode*, CuriousNode*)> > CuriousQueue;

BaseNode* base_queue_front(BaseQueue* q);

CuriousNode* curious_queue_front(CuriousQueue* q);

template<class N>
using construct_signature = N* (*)(size_t, size_t, bool, bool, double, double,
                                   N* parent, int, int, int, double);

BaseNode* base_construct_policy(size_t new_rule, size_t nrules,
                                bool prediction, bool default_prediction,
                                double lower_bound, double objective,
                                BaseNode* parent, int num_not_captured,
                                int nsamples, int len_prefix, double c);

CuriousNode* curious_construct_policy(size_t new_rule, size_t nrules,
                                      bool prediction, bool default_prediction,
                                      double lower_bound, double objective,
                                      CuriousNode* parent, int num_not_captured,
                                      int nsamples, int len_prefix, double c);

template<class N, class Q, class P>
extern void evaluate_children(CacheTree<N>* tree, N* parent,
                              VECTOR parent_not_captured,
                              std::set<size_t> ordered_parent,
                              construct_signature<N> construct_policy,
                              Q* q, struct time*, P* p);

template<class N, class P>
extern std::pair<N*, std::set<size_t> > stochastic_select(CacheTree<N>* tree, VECTOR not_captured, P* p);

template<class N>
extern void bbound_stochastic(CacheTree<N>* tree,
                              size_t max_num_nodes,
                              construct_signature<N> construct_policy,
                              struct time*);

template<class N, class Q>
extern std::pair<N*, std::set<size_t> >
queue_select(CacheTree<N>* tree, Q* q, N*(*front)(Q*), VECTOR captured);

template<class N, class Q, class P>
extern void bbound_queue(CacheTree<N>* tree,
                         size_t max_num_nodes,
                         construct_signature<N> construct_policy,
                         Q* q, N*(*front)(Q*),
                         struct time*, P* p);

/*
 * Permutation Map
 */

template<class N, class K>
class NullPermutationMap {
  public:
    N* permutation_insert(construct_signature<N> construct_policy, size_t new_rule, 
                                                size_t nrules, bool prediction, bool default_prediction, double lower_bound, 
                                                double objective, N* parent, int num_not_captured, int nsamples, int len_prefix, 
                                                double c, CacheTree<N>* tree, std::set<size_t> key)  {
        return NULL;
    };
    std::set<size_t> get_key(std::set<size_t> ordered_prefix, VECTOR captured) {};
    void remove_node(N* node) {};
    std::map<K, N*> permutation_map_;
};



typedef std::set<size_t> PrefixKey;
//typedef std::map<PrefixKey, BaseNode*> PrefixPermutationMap;
template<class N>
class PrefixPermutationMap {
    private:
        typename std::map<PrefixKey, N*>::iterator iter;
    public:
        N* permutation_insert(construct_signature<N> construct_policy, size_t new_rule, 
                                                size_t nrules, bool prediction, bool default_prediction, double lower_bound, 
                                                double objective, N* parent, int num_not_captured, int nsamples, int len_prefix, 
                                                double c, CacheTree<N>* tree, std::set<size_t> key);
        inline PrefixKey get_key(std::set<size_t> ordered_prefix, VECTOR captured);
        inline void remove_node(N* node);
        std::map<PrefixKey, N*> permutation_map_;
};

template<class N>
inline void PrefixPermutationMap<N>::remove_node(N* node) {
    typename std::map<PrefixKey, N*>::iterator iter;
    std::set<size_t> key;
    size_t depth = node->depth();
    for(size_t i = 0; i < depth; ++i) {
        key.insert(node->id());
        node = node->parent();
    }
    iter = permutation_map_.find(key);
    permutation_map_.erase(iter);
}

template <class N>
inline std::set<size_t> PrefixPermutationMap<N>::get_key(std::set<size_t> ordered_prefix, VECTOR captured) {
    (void) captured;
    return ordered_prefix;
}

template<class N, class P>
extern void delete_subtree(CacheTree<N>* tree, N* node, P* p);

void bbound_greedy(size_t nsamples, size_t nrules, rule_t *rules, rule_t *labels, size_t max_prefix_length);

#endif
