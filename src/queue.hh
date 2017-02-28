#pragma once
#include "pmap.hh"
#include <functional>
#include <queue>
#include <set>

/*
 * Queue class -- performs BFS
 */
std::function<bool(Node*, Node*)> base_cmp = [](Node* left, Node* right) {
    return left->id() > right->id();
};

class BaseQueue {
    public:
        BaseQueue(); 
        Node* front() {
            return get_q()->top();
        }
        inline void pop() {
            get_q()->pop();
        }
        inline void push(Node* node) {
            get_q()->push(node);
        }
        virtual inline std::priority_queue<Node*, std::vector<Node*>, 
                std::function<bool(Node*, Node*)> >* get_q() {
            return q_;
        }
        inline size_t size() {
            return get_q()->size();
        }
        inline bool empty() {
            return get_q()->empty();
        }

    private:
        std::priority_queue<Node*, std::vector<Node*>, std::function<bool(Node*, Node*)> >* q_;
};

BaseQueue::BaseQueue()
    : q_(new std::priority_queue<Node*, std::vector<Node*>, 
            std::function<bool(Node*, Node*)> > (base_cmp)) {}

/*
 * CuriousQueue -- orders based on curiosity metric
 */
std::function<bool(Node*, Node*)> curious = [](Node* left, Node* right) {
    return left->get_storage() > right->get_storage();
};

class CuriousQueue : public BaseQueue {
    public:
        CuriousQueue(); 
        inline std::priority_queue<Node*, std::vector<Node*>, 
               std::function<bool(Node*, Node*)> >* get_q() override {
            return q_;
        }
    protected:
        std::priority_queue<Node*, std::vector<Node*>, std::function<bool(Node*, Node*)> >* q_;
};

CuriousQueue::CuriousQueue()
    : q_(new std::priority_queue<Node*, std::vector<Node*>, 
            std::function<bool(Node*, Node*)> > (curious)) {}

/*
 * LowerBoundQueue -- orders based on curiosity metric
 */
std::function<bool(Node*, Node*)> lb = [](Node* left, Node* right) {
    return left->lower_bound() > right->lower_bound();
};

class LowerBoundQueue : public BaseQueue {
    public:
        LowerBoundQueue(); 
        inline std::priority_queue<Node*, std::vector<Node*>, 
               std::function<bool(Node*, Node*)> >* get_q() override {
            return q_;
        }
    protected:
        std::priority_queue<Node*, std::vector<Node*>, std::function<bool(Node*, Node*)> >* q_;
};

LowerBoundQueue::LowerBoundQueue()
    : q_(new std::priority_queue<Node*, std::vector<Node*>, 
            std::function<bool(Node*, Node*)> > (lb)) {}

/*
 * ObjectiveQueue -- orders based on curiosity metric
 */
std::function<bool(Node*, Node*)> objective = [](Node* left, Node* right) {
    return left->objective() > right->objective();
};

class ObjectiveQueue : public BaseQueue {
    public:
        ObjectiveQueue(); 
        inline std::priority_queue<Node*, std::vector<Node*>, 
               std::function<bool(Node*, Node*)> >* get_q() override {
            return q_;
        }
    protected:
        std::priority_queue<Node*, std::vector<Node*>, std::function<bool(Node*, Node*)> >* q_;
};

ObjectiveQueue::ObjectiveQueue()
    : q_(new std::priority_queue<Node*, std::vector<Node*>, 
            std::function<bool(Node*, Node*)> > (objective)) {}

/*
 * DFSQueue -- orders based on curiosity metric
 */
std::function<bool(Node*, Node*)> dfs = [](Node* left, Node* right) {
    return left->depth() > right->depth();
};

class DFSQueue : public BaseQueue {
    public:
        DFSQueue(); 
        inline std::priority_queue<Node*, std::vector<Node*>, 
               std::function<bool(Node*, Node*)> >* get_q() override {
            return q_;
        }
    protected:
        std::priority_queue<Node*, std::vector<Node*>, std::function<bool(Node*, Node*)> >* q_;
};

DFSQueue::DFSQueue()
    : q_(new std::priority_queue<Node*, std::vector<Node*>, 
            std::function<bool(Node*, Node*)> > (dfs)) {}


class NullQueue : public BaseQueue {
  public:
    inline void push(Node*) {};
    inline size_t size() {return 0;};
};

Node* stochastic_select(CacheTree* tree, VECTOR not_captured);

void bbound_stochastic(CacheTree* tree, size_t max_num_nodes, PermutationMap* p);

Node* queue_select(CacheTree* tree, BaseQueue* q, VECTOR captured);

int bbound_queue(CacheTree* tree, size_t max_num_nodes, BaseQueue* q, 
                 PermutationMap* p, size_t num_iter, size_t switch_iter);

void evaluate_children(CacheTree* tree, Node* parent, VECTOR parent_not_captured,
                       BaseQueue* q, PermutationMap* p);
