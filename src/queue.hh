#pragma once

#include "pmap.hh"
#include "alloc.hh"
#include <assert.h>
#include <functional>
#include <queue>
#include <set>

class InternalRoot {
    public:
        InternalRoot() {
            node_ = NULL;
            rules_ = NULL;
        }
        InternalRoot(Node* node, tracking_vector<unsigned short, DataStruct::Tree>* rules) {
            node_ = node;
            rules_ = rules;
        }
        inline Node* node() { return node_; }
        inline tracking_vector<unsigned short, DataStruct::Tree>* rules() { return rules_; }
    protected:
        Node* node_;
        tracking_vector<unsigned short, DataStruct::Tree>* rules_;
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

class Queue {
    public:
        Queue(std::function<bool(EntryType, EntryType)> cmp, char const *type);
        // by default, initialize this as a BFS queue
        Queue() : Queue(base_cmp, "BFS") {};
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

        std::pair<EntryType, tracking_vector<unsigned short, DataStruct::Tree> > select(CacheTree* tree, VECTOR captured, unsigned short thread_id) {
begin:
            int cnt;
            tracking_vector<unsigned short, DataStruct::Tree> prefix;
            Node *node;
            EntryType iroot;
            tracking_vector<unsigned short, DataStruct::Tree> init_rules;
            bool valid = true;
            double lb;
            do {
                iroot = q_->top();
                node = iroot->node();
                q_->pop();
                // We can arrive at a situation where a thread is slow to start
                // and therefore the root hasn't been popped off for the last thread
                // and it already has children -- this is not an error and we shouldn't fail the assert
                // TODO: selected_node can't have children in your split
                //if (selected_node != tree->root())
                //    assert (selected_node->num_children() == 0);
                if (tree->ablation() != 2)
                    lb = node->lower_bound() + tree->c();
                else
                    lb = node->lower_bound();
                logger->setCurrentLowerBound(lb);

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
                return std::make_pair((EntryType)NULL, prefix);
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
            return std::make_pair(iroot, prefix);
        }

    protected:
        q* q_;
        char const *type_;
};


#include "shared_queue.hh"
extern int bbound(CacheTree* tree, size_t max_num_nodes, Queue* q, PermutationMap* p,
    unsigned short thread_id, SharedQueue* shared_q);

extern void bbound_init(CacheTree* tree);

extern void evaluate_children(CacheTree* tree, Node* parent,
    tracking_vector<unsigned short, DataStruct::Tree> parent_prefix,
    VECTOR parent_not_captured, std::vector<unsigned short> rules, Queue* q,
    PermutationMap* p, unsigned short thread_id);

extern bool bbound_loop(CacheTree* tree, size_t max_num_nodes, Queue* q, PermutationMap* p,
    VECTOR captured, VECTOR not_captured, unsigned short thread_id, SharedQueue *shared_q);
