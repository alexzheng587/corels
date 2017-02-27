#pragma once
#include "pmap.hh"
#include <functional>
#include <queue>
#include <set>

/*
 * Queue class -- performs BFS
 */
template<class N>
std::function<bool(N*, N*)> base = [](N* left, N* right) {
    return left->id() > right->id();
};

template<class N>
class Queue {
    public:
        Queue(); 
        N* front() {
            return get_q()->top();
        }
        void pop() {
            get_q()->pop();
        }
        void push(N* node) {
            get_q()->push(node);
        }
        virtual std::priority_queue<N*, std::vector<N*>, std::function<bool(N*, N*)> >* get_q() {
            return q_;
        }
        size_t size() {
            return get_q()->size();
        }
        bool empty() {
            return get_q()->empty();
        }

    private:
        std::priority_queue<N*, std::vector<N*>, std::function<bool(N*, N*)> >* q_;
};

template<class N>
Queue<N>::Queue()
    : q_(new std::priority_queue<N*, std::vector<N*>, std::function<bool(N*, N*)> > (base<N>)) {}

/*
 * CuriousQueue -- orders based on curiosity metric
 */
template<class N>
std::function<bool(N*, N*)> curious = [](N* left, N* right) {
    return left->get_storage() > right->get_storage();
};

template<class N>
class CuriousQueue : public Queue<N> {
    public:
        CuriousQueue(); 
        std::priority_queue<N*, std::vector<N*>, std::function<bool(N*, N*)> >* get_q() override {
            return q_;
        }
    protected:
        std::priority_queue<N*, std::vector<N*>, std::function<bool(N*, N*)> >* q_;
};

template<class N>
CuriousQueue<N>::CuriousQueue()
    : q_(new std::priority_queue<N*, std::vector<N*>, std::function<bool(N*, N*)> > (curious<N>)) {}

/*
 * LowerBoundQueue -- orders based on curiosity metric
 */
template<class N>
std::function<bool(N*, N*)> lb = [](N* left, N* right) {
    return left->lower_bound() > right->lower_bound();
};

template<class N>
class LowerBoundQueue : public Queue<N> {
    public:
        LowerBoundQueue(); 
        std::priority_queue<N*, std::vector<N*>, std::function<bool(N*, N*)> >* get_q() override {
            return q_;
        }
    protected:
        std::priority_queue<N*, std::vector<N*>, std::function<bool(N*, N*)> >* q_;
};

template<class N>
LowerBoundQueue<N>::LowerBoundQueue()
    : q_(new std::priority_queue<N*, std::vector<N*>, std::function<bool(N*, N*)> > (lb<N>)) {}

/*
 * ObjectiveQueue -- orders based on curiosity metric
 */
template<class N>
std::function<bool(N*, N*)> objective = [](N* left, N* right) {
    return left->objective() > right->objective();
};

template<class N>
class ObjectiveQueue : public Queue<N> {
    public:
        ObjectiveQueue(); 
        std::priority_queue<N*, std::vector<N*>, std::function<bool(N*, N*)> >* get_q() override {
            return q_;
        }
    protected:
        std::priority_queue<N*, std::vector<N*>, std::function<bool(N*, N*)> >* q_;
};

template<class N>
ObjectiveQueue<N>::ObjectiveQueue()
    : q_(new std::priority_queue<N*, std::vector<N*>, std::function<bool(N*, N*)> > (objective<N>)) {}

/*
 * DFSQueue -- orders based on curiosity metric
 */
template<class N>
std::function<bool(N*, N*)> dfs = [](N* left, N* right) {
    return left->depth() > right->depth();
};

template<class N>
class DFSQueue : public Queue<N> {
    public:
        DFSQueue(); 
        std::priority_queue<N*, std::vector<N*>, std::function<bool(N*, N*)> >* get_q() override {
            return q_;
        }
    protected:
        std::priority_queue<N*, std::vector<N*>, std::function<bool(N*, N*)> >* q_;
};

template<class N>
DFSQueue<N>::DFSQueue()
    : q_(new std::priority_queue<N*, std::vector<N*>, std::function<bool(N*, N*)> > (dfs<N>)) {}


template<class N>
class NullQueue : public Queue<N> {
  public:
    void push(N*) {};
    size_t size() {return 0;};
};

template<class N>
extern N* stochastic_select(CacheTree<N>* tree, VECTOR not_captured);

template<class N>
extern void bbound_stochastic(CacheTree<N>* tree,
                              size_t max_num_nodes,
                              construct_signature<N> construct_policy,
                              PermutationMap<N>* p);

template<class N>
extern N*
queue_select(CacheTree<N>* tree, Queue<N>* q, VECTOR captured);

template<class N>
extern int bbound_queue(CacheTree<N>* tree,
                         size_t max_num_nodes,
                         construct_signature<N> construct_policy,
                         Queue<N>* q,
                         PermutationMap<N>* p, size_t num_iter, size_t switch_iter);

template<class N>
extern void evaluate_children(CacheTree<N>* tree, N* parent,
                              VECTOR parent_not_captured,
                              construct_signature<N> construct_policy,
                              Queue<N>* q,
                              PermutationMap<N>* p);
