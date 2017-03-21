#pragma once

#include "pmap.hh"
#include "alloc.hh"
#include <functional>
#include <queue>
#include <set>

/*
template <class T> 
struct queue_alloc : track_alloc<T> {
    typedef T value_type;
    T* allocate (size_t n) {
        logger.addToQueueMemory(n * sizeof(T));
        return static_cast<T*>(malloc(n*sizeof(T)));
    }   
    void deallocate (T* p, size_t n) {
        logger.removeFromQueueMemory(n * sizeof(*p));
        free(p);
    }   
};
*/

typedef std::priority_queue<Node*, std::vector<Node*, queue_alloc<Node*> >, std::function<bool(Node*, Node*)> > q_alloc;

/*
 * Queue class -- performs BFS
 */
static std::function<bool(Node*, Node*)> base_cmp = [](Node* left, Node* right) {
    return left->id() > right->id();
};

class BaseQueue {
    public:
        BaseQueue(std::function<bool(Node*, Node*)> cmp); 
        BaseQueue() : BaseQueue(base_cmp) {};
        Node* front() {
            return get_q()->top();
        }
        inline void pop() {
            get_q()->pop();
        }
        void push(Node* node) {
            get_q()->push(node);
        }
        virtual q_alloc* get_q() {
            return q_;
        }
        size_t size() {
            return get_q()->size();
        }
        inline bool empty() {
            return get_q()->empty();
        }

    private:
        q_alloc* q_;
};

/*
 * CuriousQueue -- orders based on curiosity metric
 */
static std::function<bool(Node*, Node*)> curious = [](Node* left, Node* right) {
    return left->get_storage() > right->get_storage();
};

class CuriousQueue : public BaseQueue {
    public:
        //CuriousQueue() : BaseQueue(curious) {}; 
        CuriousQueue();
        q_alloc* get_q() override {
            return q_;
        }
    protected:
        q_alloc* q_;
};

/*
 * LowerBoundQueue -- orders based on curiosity metric
 */
static std::function<bool(Node*, Node*)> lb = [](Node* left, Node* right) {
    return left->lower_bound() > right->lower_bound();
};

class LowerBoundQueue : public BaseQueue {
    public:
        LowerBoundQueue();
        //LowerBoundQueue() : BaseQueue(lb) {}; 
        q_alloc* get_q() override {
            return q_;
        }
    protected:
        q_alloc* q_;
};

/*
 * ObjectiveQueue -- orders based on curiosity metric
 */
static std::function<bool(Node*, Node*)> objective = [](Node* left, Node* right) {
    return left->objective() > right->objective();
};

class ObjectiveQueue : public BaseQueue {
    public:
        ObjectiveQueue();
        //ObjectiveQueue() : BaseQueue(objective) {}; 
        q_alloc* get_q() override {
            return q_;
        }
    protected:
        q_alloc* q_;
};

/*
 * DFSQueue -- orders based on curiosity metric
 */
static std::function<bool(Node*, Node*)> dfs = [](Node* left, Node* right) {
    return left->depth() > right->depth();
};

class DFSQueue : public BaseQueue {
    public:
        DFSQueue();
        //DFSQueue() : BaseQueue(dfs) {}; 
        q_alloc* get_q() override {
            return q_;
        }
    protected:
        q_alloc* q_;
};

class NullQueue : public BaseQueue {
  public:
    void push(Node*) {};
    size_t size() {return 0;};
};

extern Node* stochastic_select(CacheTree* tree, VECTOR not_captured);

extern void bbound_stochastic(CacheTree* tree, size_t max_num_nodes, PermutationMap* p);

extern Node* queue_select(CacheTree* tree, BaseQueue* q, VECTOR captured);

extern int bbound_queue(CacheTree* tree, size_t max_num_nodes, BaseQueue* q, 
                 PermutationMap* p, size_t num_iter, size_t switch_iter, double* min_objective);

extern void evaluate_children(CacheTree* tree, Node* parent, VECTOR parent_not_captured,
                       BaseQueue* q, PermutationMap* p, std::vector<unsigned short>& rules, double* min_objective);
extern void bbound_queue_init(CacheTree* tree, BaseQueue* q, PermutationMap* p, std::vector<unsigned short> rules, double* min_objective);
