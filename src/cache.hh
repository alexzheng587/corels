#pragma once

#include "loss.hh"
#include "utils.hh"
#include "alloc.hh"
#include "rule.h"
#include "problemState.hh"
#include <iterator>
#include <map>
#include <vector>
#include <stdlib.h>
#include <memory>

class Loss; // forward declare because of circular dependency

class Node {
public:
    Node(bool default_prediction, double objective, double equivalent_minority);

    Node(unsigned short id, bool prediction, bool default_prediction, double lower_bound, double objective,
         Node* parent, size_t num_captured, double equivalent_minority);

    virtual ~Node() {}

    inline unsigned short id() const;

    inline bool prediction() const;

    inline bool default_prediction() const;

    inline double lower_bound() const;

    inline double objective() const;

    inline bool done() const;

    inline void set_done();

    inline bool deleted() const;

    inline void set_deleted();

    inline bool in_queue() const;

    inline void set_in_queue(bool new_val);

    // Returns pair of prefixes and predictions for the path from this node to the root
    inline std::pair<tracking_vector<unsigned short, DataStruct::Tree>, tracking_vector<bool, DataStruct::Tree> >
    get_prefix_and_predictions();

    inline size_t depth() const;

    inline Node* child(unsigned short idx);

    inline Node* parent() const;

    inline void delete_child(unsigned short idx);

    inline void clear_children();

    inline size_t num_children() const;

    inline size_t live_children();

    inline size_t num_captured() const;

    inline double equivalent_minority() const;

    inline typename std::map<unsigned short, Node*>::iterator children_begin();

    inline typename std::map<unsigned short, Node*>::iterator children_end();

    virtual inline double get_curiosity() {
        return 0.0;
    }

    inline virtual Node* construct_child(unsigned short new_rule, ProblemState* state) {
        return new Node(new_rule, state->rule_prediction, state->default_prediction, state->rule_lower_bound,
                        state->objective, this, state->n_captured, state->equivalent_points_loss);
    }

protected:
    std::map<unsigned short, Node*, std::less<unsigned short>, track_alloc<std::pair<const unsigned short, Node*>, DataStruct::Tree> > children_;
    Node* parent_;
    double lower_bound_;
    double objective_;
    double equivalent_minority_;
    size_t depth_;
    size_t num_captured_;
    unsigned short id_;
    bool prediction_;
    bool default_prediction_;
    bool done_;
    bool deleted_;
    bool in_queue_;

    friend class CacheTree;
};

class CuriousNode : public Node {
public:
    CuriousNode(unsigned short id, bool prediction, bool default_prediction, double lower_bound, double objective,
                double curiosity, CuriousNode* parent, size_t num_captured, double equivalent_minority) :
            Node(id, prediction, default_prediction, lower_bound, objective, (Node*) parent, num_captured,
                 equivalent_minority) {
        curiosity_ = curiosity;
    }

    inline Node* construct_child(unsigned short new_rule, ProblemState* state) override {
        double lower_bound = state->rule_lower_bound;
        double equivalent_points_loss = state->equivalent_points_loss;
        int num_captured = state->n_captured;
        int nsamples = num_captured + state->n_not_captured;
        double curiosity = (lower_bound - equivalent_points_loss) * nsamples / static_cast<double>(num_captured);
        return (Node*) new CuriousNode(new_rule, state->rule_prediction, state->default_prediction,
                                       lower_bound, state->objective, curiosity, this, num_captured,
                                       equivalent_points_loss);
    }

    inline double get_curiosity() override;

protected:
    double curiosity_;
};

class AUCNode : public Node {
public:
    // constructor for root
    AUCNode(bool default_prediction, double objective, double equivalent_minority, int nclasses, minority_class_t* minor_class,
            int npos, int nneg)
            : Node(default_prediction, objective, equivalent_minority) {
                npos_ = npos;
                nneg_ = nneg;
                minority_classes_ = minor_class;
                for (int i = 0; i < nclasses; i++) {
                    class_idx_.push_back(i);
                }
        }

    AUCNode(unsigned short id, bool prediction, bool default_prediction, double lower_bound,
                double objective, size_t npos, size_t nneg,  AUCNode* parent, size_t num_captured,
                double equivalent_minority, minority_class_t* minority, std::vector<int> class_idx) :
            Node(id, prediction, default_prediction, lower_bound, objective, (Node*) parent, num_captured,
                 equivalent_minority) {
                    npos_ = npos;
                    nneg_ = nneg;
                    minority_classes_ = minority;
                    class_idx_ = class_idx;
        }

        inline int get_npos() {
            return npos_;
        };

        inline int get_nneg() {
            return nneg_;
        };

        // Return the minority class
        inline minority_class_t minority_classes(unsigned short idx) const {
            int index = class_idx_[idx];
            return minority_classes_[index];
        }

        inline minority_class_t* get_minority_pointer() {
            return minority_classes_;
        }

        inline int get_nclasses() const {
            return class_idx_.size();
        }

        inline void delete_minority_class(int idx) {
            class_idx_.erase(class_idx_.begin() + idx);
        }

        inline Node* construct_child(unsigned short new_rule, ProblemState* state) override {
        Node* child = (Node*) (new AUCNode(new_rule, state->rule_prediction, state->default_prediction,
                                          state->rule_lower_bound, state->objective,
                                          state->n_ones_captured, state->n_zeros_captured, this, state->n_captured,
                                          state->equivalent_points_loss, state->class_pointer, state->class_idx_list));
        return child;
    }

protected:
    int npos_;
    int nneg_;
    minority_class_t* minority_classes_;    // Pointer to the minority class of the entire tree, only allocated once
    std::vector<int> class_idx_;             // Keeps track of the index of the minority classes this node still has
};

class F1Node : public Node {
public:
    // Constructor for root
    F1Node(bool default_prediction, double objective, double equivalent_minority)
            : Node(default_prediction, objective, equivalent_minority) {
        // Set to 0 in order to accumulate false negatives and false positives as we construct the rule list.
        // If we included the default rule's points, then it would become more difficult to compute bounds
        // before doing all the necessary bitvector operations.
        fns_ = 0;
        fps_ = 0;
    }

    F1Node(unsigned short id, bool prediction, bool default_prediction, double lower_bound, double objective,
           size_t false_positives, size_t false_negatives, F1Node* parent, size_t num_captured,
           double equivalent_minority) :
            Node(id, prediction, default_prediction, lower_bound, objective, (Node*) parent, num_captured,
                 equivalent_minority) {
        fns_ = false_negatives;
        fps_ = false_positives;
    }

    inline size_t get_false_negatives() const { return fns_; }

    inline size_t get_false_positives() const { return fps_; }

    inline Node* construct_child(unsigned short new_rule, ProblemState* state) override {
        Node* child = (Node*) (new F1Node(new_rule, state->rule_prediction, state->default_prediction,
                                          state->rule_lower_bound, state->objective,
                                          state->false_positives, state->false_negatives, this, state->n_captured,
                                          state->equivalent_points_loss));
        return child;
    }

protected:
    size_t fns_; // Cumulative false negatives
    size_t fps_; // Cumulative false positives
};



class CacheTree {
public:
    CacheTree(Loss* loss, size_t nsamples, size_t nrules, double c, size_t k, rule_t* rules, rule_t* labels,
              rule_t* minority, minority_class_t* minority_class, size_t nclasses, int ablation, bool calculate_size, char const* node_type);

    ~CacheTree();

    // Node*
    // construct_node(unsigned short new_rule, size_t nrules, bool prediction, bool default_prediction, double lower_bound,
    //                double objective, Node* parent, int num_not_captured, int nsamples, double equivalent_minority);

    // Node*
    // construct_auc_node(unsigned short new_rule, size_t nrules, bool prediction, bool default_prediction, double lower_bound,
    //                double objective, Node* parent, int num_not_captured, int nsamples, double equivalent_minority,
    //                minority_class_t* minority_class, int npos_captured, int nneg_captured);

    inline double min_objective(int idx) const;

    inline tracking_vector<unsigned short, DataStruct::Tree> opt_rulelist(int idx) const;

    inline tracking_vector<bool, DataStruct::Tree> opt_predictions(int idx) const;

    inline Loss* loss() const;

    inline size_t num_nodes() const;

    inline size_t num_evaluated() const;

    inline rule_t rule(unsigned short idx) const;

    inline char* rule_features(unsigned short idx) const;

    inline rule_t label(unsigned short idx) const;

    inline rule_t minority(unsigned short idx) const;

    inline bool has_minority() const;

    inline size_t nsamples() const;

    inline size_t nneg() const;

    inline size_t npos() const;

    inline size_t nrules() const;

    inline double c() const;

    inline size_t k() const;

    inline size_t nrulelists() const;

    inline Node* root() const;

    inline ProblemState* threadState(int thread_index) const;

    void update_opt_rulelist(double objective, tracking_vector<unsigned short, DataStruct::Tree>& parent_prefix,
                             unsigned short new_rule_id, Node* parent, bool new_pred, bool new_default_pred);

    inline void increment_num_evaluated();

    inline void decrement_num_nodes();

    inline int ablation() const;

    inline bool calculate_size() const;

    void insert_root();

    void insert(Node* node);

    void prune_up(Node* node);

    void garbage_collect();

    void play_with_rules();

    Node* check_prefix(tracking_vector<unsigned short, DataStruct::Tree>& prefix);

    bool support_bound_enabled() const;

    bool lookahead_bound_enabled() const;

protected:
    Loss* loss_;
    Node* root_;
    vector<unique_ptr<ProblemState>> thread_states_;
    size_t nsamples_;
    size_t npos_;
    size_t nneg_;
    size_t nrules_;
    size_t nclasses_;
    double c_;
    size_t k_;

    size_t num_nodes_;
    size_t num_evaluated_;
    int ablation_; // Used to remove support (1) or lookahead (2) bounds
    bool calculate_size_;

    struct optimal_rule {
        double min_objective_;
        tracking_vector<unsigned short, DataStruct::Tree> opt_rulelist_;
        std::vector<bool, track_alloc<bool, DataStruct::Tree> > opt_predictions_;
    };

    tracking_vector<optimal_rule, DataStruct::Tree> optimal_rules; // Ordered from most to least, i.e. the first element is the most optimal

    rule_t* rules_;
    rule_t* labels_;
    rule_t* minority_;
    minority_class_t* minor_class_;

    char const* node_type_;

protected:

    void gc_helper(Node* node);
};

inline unsigned short Node::id() const {
    return id_;
}

inline bool Node::prediction() const {
    return prediction_;
}

inline bool Node::default_prediction() const {
    return default_prediction_;
}

inline double Node::lower_bound() const {
    return lower_bound_;
}

inline double Node::objective() const {
    return objective_;
}

inline bool Node::done() const {
    return done_;
}

inline void Node::set_done() {
    done_ = 1;
}

inline bool Node::deleted() const {
    return deleted_;
}

inline void Node::set_deleted() {
    deleted_ = 1;
}

inline bool Node::in_queue() const {
    return in_queue_;
}

inline void Node::set_in_queue(bool new_val) {
    in_queue_ = new_val;
}

inline std::pair<tracking_vector<unsigned short, DataStruct::Tree>, tracking_vector<bool, DataStruct::Tree> >
Node::get_prefix_and_predictions() {
    tracking_vector<unsigned short, DataStruct::Tree> prefix;
    tracking_vector<bool, DataStruct::Tree> predictions;
    tracking_vector<unsigned short, DataStruct::Tree>::iterator it1 = prefix.begin();
    tracking_vector<bool, DataStruct::Tree>::iterator it2 = predictions.begin();
    Node* node = this;
    for (size_t i = depth_; i > 0; --i) {
        it1 = prefix.insert(it1, node->id());
        it2 = predictions.insert(it2, node->prediction());
        node = node->parent();
    }
    return std::make_pair(prefix, predictions);
}

inline size_t Node::depth() const {
    return depth_;
}

inline Node* Node::child(unsigned short idx) {
    typename std::map<unsigned short, Node*>::iterator iter;
    iter = children_.find(idx);
    if (iter == children_.end())
        return NULL;
    else
        return iter->second;
}

inline void Node::delete_child(unsigned short idx) {
    children_.erase(idx);
}


inline void Node::clear_children() {
    children_.clear();
}

inline size_t Node::num_children() const {
    return children_.size();
}


inline size_t Node::live_children() {
    if (children_.size() == 0) {
        return 0;
    }
    size_t nchildren = 0;
    for (auto child_it = children_begin(); child_it != children_end(); ++child_it) {
        Node* child = child_it->second;
        if (!child->deleted_ && !child->done_) {
            nchildren += 1;
        }
    }
    return nchildren;
}

inline typename std::map<unsigned short, Node*>::iterator Node::children_begin() {
    return children_.begin();
}

inline typename std::map<unsigned short, Node*>::iterator Node::children_end() {
    return children_.end();
}

inline Node* Node::parent() const {
    return parent_;
}

inline size_t Node::num_captured() const {
    return num_captured_;
}

inline double Node::equivalent_minority() const {
    return equivalent_minority_;
}

inline double CuriousNode::get_curiosity() {
    return curiosity_;
}

inline double CacheTree::min_objective(int idx = 0) const {
    return optimal_rules[idx].min_objective_;
}

inline tracking_vector<unsigned short, DataStruct::Tree> CacheTree::opt_rulelist(int idx = 0) const {
    return optimal_rules[idx].opt_rulelist_;
}

inline tracking_vector<bool, DataStruct::Tree> CacheTree::opt_predictions(int idx = 0) const {
    return optimal_rules[idx].opt_predictions_;
}

inline Loss* CacheTree::loss() const {
    return loss_;
}

inline size_t CacheTree::num_nodes() const {
    return num_nodes_;
}

inline size_t CacheTree::num_evaluated() const {
    return num_evaluated_;
}

inline rule_t CacheTree::rule(unsigned short idx) const {
    return rules_[idx];
}

inline char* CacheTree::rule_features(unsigned short idx) const {
    return rules_[idx].features;
}

inline rule_t CacheTree::label(unsigned short idx) const {
    return labels_[idx];
}

inline rule_t CacheTree::minority(unsigned short idx) const {
    return minority_[idx];
}

inline bool CacheTree::has_minority() const {
    return minority_ != NULL;
}

inline size_t CacheTree::nsamples() const {
    return nsamples_;
}

inline size_t CacheTree::npos() const {
    return npos_;
}

inline size_t CacheTree::nneg() const {
    return nneg_;
}

inline size_t CacheTree::nrules() const {
    return nrules_;
}

inline double CacheTree::c() const {
    return c_;
}

inline size_t CacheTree::k() const {
    return k_;
}

inline size_t CacheTree::nrulelists() const {
    return optimal_rules.size();
}

inline Node* CacheTree::root() const {
    return root_;
}

inline int CacheTree::ablation() const {
    return ablation_;
}

inline bool CacheTree::calculate_size() const {
    return calculate_size_;
}

/*
 * Update the optimal rulelist of the tree.
 */
inline void
CacheTree::update_opt_rulelist(double objective, tracking_vector<unsigned short, DataStruct::Tree>& parent_prefix,
                               unsigned short new_rule_id, Node* parent, bool new_pred, bool new_default_pred) {
    optimal_rules.insert(optimal_rules.begin(), optimal_rule());
    if (optimal_rules.size() > k_)
        optimal_rules.erase(optimal_rules.begin() + k_, optimal_rules.end());

    optimal_rules[0].min_objective_ = objective;

    optimal_rules[0].opt_rulelist_.assign(parent_prefix.begin(), parent_prefix.end());
    optimal_rules[0].opt_rulelist_.push_back(new_rule_id);
    logger->setTreePrefixLen(optimal_rules[0].opt_rulelist_.size());

    tracking_vector<bool, DataStruct::Tree> predictions;
    Node* node = parent;
    for (size_t i = parent->depth(); i > 0; --i) {
        predictions.push_back(node->prediction());
        node = node->parent();
    }
    std::reverse(predictions.begin(), predictions.end());
    optimal_rules[0].opt_predictions_.assign(predictions.begin(), predictions.end());
    optimal_rules[0].opt_predictions_.push_back(new_pred);
    optimal_rules[0].opt_predictions_.push_back(new_default_pred);
}

/*
 * Increment number of nodes evaluated after performing incremental computation
 * in evaluate_children.
 */
inline void CacheTree::increment_num_evaluated() {
    ++num_evaluated_;
    logger->setTreeNumEvaluated(num_evaluated_);
}

/*
 * Called whenever a node is deleted from the tree.
 */
inline void CacheTree::decrement_num_nodes() {
    --num_nodes_;
    logger->setTreeNumNodes(num_nodes_);
}

inline bool CacheTree::support_bound_enabled() const {
    return ablation_ != 1;
}

inline bool CacheTree::lookahead_bound_enabled() const {
    return ablation_ != 2;
}

ProblemState* CacheTree::threadState(int thread_index) const {
    return thread_states_[thread_index].get();
}


void delete_subtree(CacheTree* tree, Node* node, bool destructive, bool update_remaining_state_space);