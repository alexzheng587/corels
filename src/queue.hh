#pragma once

#include "pmap.hh"
#include "alloc.hh"
#include "shared_queue.hh"
#include <assert.h>
#include <functional>
#include <queue>
#include <set>

// pass custom allocator function to track memory allocations in the queue
typedef std::priority_queue<Node*, tracking_vector<Node*, DataStruct::Queue>, 
        std::function<bool(Node*, Node*)> > q;

// orders based on depth (BFS)
static std::function<bool(Node*, Node*)> base_cmp = [](Node* left, Node* right) {
    return left->depth() >= right->depth();
};

// orders based on curiosity metric.
static std::function<bool(Node*, Node*)> curious_cmp = [](Node* left, Node* right) {
    return left->get_curiosity() >= right->get_curiosity();
};

// orders based on lower bound.
static std::function<bool(Node*, Node*)> lb_cmp = [](Node* left, Node* right) {
    return left->lower_bound() >= right->lower_bound();
};

// orders based on objective.
static std::function<bool(Node*, Node*)> objective_cmp = [](Node* left, Node* right) {
    return left->objective() >= right->objective();
};

// orders based on depth (DFS)
static std::function<bool(Node*, Node*)> dfs_cmp = [](Node* left, Node* right) {
    return left->depth() <= right->depth();
};

class Queue {
    public:
        Queue(std::function<bool(Node*, Node*)> cmp, char const *type);
        // by default, initialize this as a BFS queue
        Queue() : Queue(base_cmp, "BFS") {};
        Node* front() {
            return q_->top();
        }
        inline void pop() {
            q_->pop();
        }
        void push(Node* node) {
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

        std::pair<Node*, tracking_vector<unsigned short, DataStruct::Tree> > select(CacheTree* tree, VECTOR captured, unsigned short thread_id) {
begin:
            int cnt;
            tracking_vector<unsigned short, DataStruct::Tree> prefix;
            Node *selected_node, *node;
            bool valid = true;
            double lb;
            do {
                selected_node = q_->top();
                q_->pop();
                // We can arrive at a situation where a thread is slow to start
                // and therefore the root hasn't been popped off for the last thread
                // and it already has children -- this is not an error and we shouldn't fail the assert
                if (selected_node != tree->root())
                    assert (selected_node->num_children() == 0);
                if (tree->ablation() != 2)
                    lb = selected_node->lower_bound() + tree->c();
                else
                    lb = selected_node->lower_bound();
                logger->setCurrentLowerBound(lb);

                node = selected_node;
                node->set_in_queue(false);
                // delete leaf nodes that were lazily marked
                if (node->deleted() || (lb >= tree->min_objective())) {
                    if(featureDecisions->do_garbage_collection()) {
                        tree->decrement_num_nodes();
                        logger->removeFromMemory(sizeof(*node), DataStruct::Tree);
                        tree->lock(thread_id);
                        if (!node->deleted()) {
                            Node* parent = node->parent();
                            parent->delete_child(node->id());
                        }
                        node->clear_children();
                        delete node;
                        tree->unlock(thread_id);
                    }
                    valid = false;
                } else {
                    valid = true;
                }
            } while (!q_->empty() && !valid); 
            if (!valid) {
                return std::make_pair((Node*)NULL, prefix);
            }

            rule_vclear(tree->nsamples(), captured);
            while (node != tree->root()) {
                // need to delete interior nodes lazily too when parallel
                if(node->deleted()) {
                    if(featureDecisions->do_garbage_collection()) {
                        // if the node is a leaf node, we can physically delete it (destructive mode)
                        // otherwise, call delete_subtree non-destructively
                        Node* parent = node->parent();
                        parent->delete_child(node->id());
                        delete_subtree(tree, node, (node->live_children() == 0), false, thread_id);
                    }
                    goto begin;
                }
                rule_vor(captured,
                         captured, tree->rule(node->id()).truthtable,
                         tree->nsamples(), &cnt);
                prefix.push_back(node->id());
                node = node->parent();
            }
            std::reverse(prefix.begin(), prefix.end());
            return std::make_pair(selected_node, prefix);
        }

    protected:
        q* q_;
        char const *type_;
};


extern int bbound(CacheTree* tree, size_t max_num_nodes, Queue* q, PermutationMap* p,
    unsigned short thread_id, SharedQueue* shared_q);

extern void bbound_init(CacheTree* tree);

extern void evaluate_children(CacheTree* tree, Node* parent,
    tracking_vector<unsigned short, DataStruct::Tree> parent_prefix,
    VECTOR parent_not_captured, std::vector<unsigned short> rules, Queue* q,
    PermutationMap* p, unsigned short thread_id);

extern bool bbound_loop(CacheTree* tree, size_t max_num_nodes, Queue* q, PermutationMap* p,
    VECTOR captured, VECTOR not_captured, unsigned short thread_id, bool special_call);
