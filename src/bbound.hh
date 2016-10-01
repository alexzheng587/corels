#include "cache.hh"
#include <set>
#include <queue>

template<class N>
class NullQueue {
  public:
    void push(N* node) {};
};

typedef std::queue<BaseNode*> BaseQueue;

template<class N>
using construct_signature = N* (*)(size_t, size_t, bool, bool, double, double, N* parent, int, int, int, double);

BaseNode* base_construct_policy(size_t new_rule, size_t nrules, bool prediction, bool default_prediction,
                                double lower_bound, double objective, BaseNode* parent, int num_not_captured,
                                int nsamples, int len_prefix, double c);

CuriousNode* curious_construct_policy(size_t new_rule, size_t nrules, bool prediction, bool default_prediction,
                                      double lower_bound, double objective, CuriousNode* parent, int num_not_captured,
                                      int nsamples, int len_prefix, double c);

template<class N, class Q>
extern void evaluate_children(CacheTree<N>* tree, N* parent,
                              VECTOR parent_not_captured,
                              std::set<size_t> ordered_parent,
                              construct_signature<N> construct_policy,
                              Q* q, struct time*);

template<class N>
extern std::pair<N*, std::set<size_t> > stochastic_select(CacheTree<N>* tree, VECTOR not_captured);

template<class N>
extern void bbound_stochastic(CacheTree<N>* tree,
                              size_t max_num_nodes,
                              construct_signature<N> construct_policy,
                              struct time*);

template<class N, class Q>
extern std::pair<N*, std::set<size_t> >
queue_select(CacheTree<N>* tree, Q* q, VECTOR captured);

template<class N, class Q>
extern void bbound_queue(CacheTree<N>* tree,
                         size_t max_num_nodes,
                         construct_signature<N> construct_policy,
                         Q* q,
                         struct time*);

void bbound_greedy(size_t nsamples, size_t nrules, rule_t *rules, rule_t *labels, size_t max_prefix_length);
