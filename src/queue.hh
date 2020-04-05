#pragma once

#include "pmap.hh"
#include "alloc.hh"
#include <assert.h>
#include <functional>
#include <queue>
#include <set>

extern std::mutex min_obj_lk;

class InternalRoot {
    public:
        InternalRoot() {
            node_ = NULL;
        }
        InternalRoot(Node* node) {
            node_ = node;
        }
        InternalRoot(Node* node, tracking_vector<unsigned short, DataStruct::Tree> rules) {
            node_ = node;
            rules_ = rules;
        }
        ~InternalRoot() {
            //printf("delete iroot with node %p\n", node_);
            node_ = (Node*) 0xDEADBEEF;
            /*if (node_ != NULL) {
                delete node_;
            }*/
        }
        inline Node* node() { return node_; }
        inline tracking_vector<unsigned short, DataStruct::Tree> rules() { return rules_; }
    protected:
        Node* node_;
        tracking_vector<unsigned short, DataStruct::Tree> rules_;
};

typedef InternalRoot* EntryType;

// pass custom allocator function to track memory allocations in the queue
typedef std::priority_queue<EntryType, tracking_vector<EntryType, DataStruct::Queue>,
        std::function<bool(EntryType, EntryType)> > q;

// orders based on depth (BFS)
static std::function<bool(EntryType, EntryType)> base_cmp = [](EntryType left, EntryType right) {
    return left->node()->depth() >= right->node()->depth();
};

// orders based on curiosity metric.
static std::function<bool(EntryType, EntryType)> curious_cmp = [](EntryType left, EntryType right) {
    return left->node()->get_curiosity() >= right->node()->get_curiosity();
};

// orders based on lower bound.
static std::function<bool(EntryType, EntryType)> lb_cmp = [](EntryType left, EntryType right) {
    return left->node()->lower_bound() >= right->node()->lower_bound();
};

// orders based on objective.
static std::function<bool(EntryType, EntryType)> objective_cmp = [](EntryType left, EntryType right) {
    return left->node()->objective() >= right->node()->objective();
};

// orders based on depth (DFS)
static std::function<bool(EntryType, EntryType)> dfs_cmp = [](EntryType left, EntryType right) {
    return left->node()->depth() <= right->node()->depth();
};

template <class T, class S, class C>
    S& Container(priority_queue<T, S, C>& q) {
        struct HackedQueue : private priority_queue<T, S, C> {
            static S& Container(priority_queue<T, S, C>& q) {
                return q.*&HackedQueue::c;
            }
        };
    return HackedQueue::Container(q);
}

class Queue {
    public:
        Queue(std::function<bool(EntryType, EntryType)> cmp, char const *type);
        // by default, initialize this as a BFS queue
        Queue() : Queue(base_cmp, "BFS") {}
        ~Queue();

        EntryType front() {
            return q_->top();
        }
        inline void pop() {
            q_->pop();
        }
        void push(EntryType node) {
            q_->push(node);
        }
        size_t size() {
            return q_->size();
        }
        bool empty() {
            return q_->empty();
        }
        inline char const * type() {
            return type_;
        }
        /*void move(Queue* other_q, size_t n_elt) {
            tracking_vector<EntryType, DataStruct::Queue> v1 = Container(*q_);
            tracking_vector<EntryType, DataStruct::Queue> v2 = Container(*other_q->q_);
            v2.insert(v2.end(), std::make_move_iterator(v1.begin() + n_elt), 
                                std::make_move_iterator(v1.end()));
            v1.erase(v1.begin() + n_elt, v1.end());
            std::make_heap(v1.begin(), v1.end(), cmp_);
            std::make_heap(v2.begin(), v2.end(), other_q->cmp_);
        }*/
        void move(Queue* other_q, size_t n_elt) {
            bool odd = true;
            std::vector<EntryType> backup_vec;
            for(size_t i = 0; i < n_elt; ++i) {
                EntryType elt = front();
                pop();
                if (odd) {
                    other_q->push(elt);
                } else {
                    backup_vec.push_back(elt);
                }
                odd = !odd;
            }
            for(auto it = backup_vec.begin(); it != backup_vec.end(); ++it) {
                push(*it);
            }
        }

        InternalRoot* select() {
            InternalRoot* iroot = front();
            pop();
            return iroot;
        }

    protected:
        q* q_;
        char const *type_;
        std::function<bool(EntryType, EntryType)> cmp_;
};


#include "shared_queue.hh"
extern int bbound(CacheTree* tree, size_t max_num_nodes, Queue* q, PermutationMap* p,
    unsigned short thread_id, SharedQueue* shared_q);

extern void bbound_init(CacheTree* tree);

extern void evaluate_children(CacheTree* tree, Node* parent,
    tracking_vector<unsigned short, DataStruct::Tree> parent_prefix,
    VECTOR parent_not_captured, std::vector<unsigned short> rules, Queue* q,
    PermutationMap* p, unsigned short thread_id);

extern bool is_valid_node(Node* node, CacheTree* tree);
extern tracking_vector<unsigned short, DataStruct::Tree> get_parent_prefix(CacheTree* tree, Node* node, VECTOR captured);

extern bool bbound_loop(CacheTree* tree, size_t max_num_nodes, Queue* q, PermutationMap* p,
    VECTOR captured, VECTOR not_captured, unsigned short thread_id, SharedQueue *shared_q);


bool bbound_loop_cond(bool max_node_reached, SharedQueue* shared_q, CacheTree* tree);
