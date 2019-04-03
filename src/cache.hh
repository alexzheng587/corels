#pragma once

#include "utils.hh"
#include "alloc.hh"
#include "rule.h"
#include <iterator>
#include <map>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <memory>
#include <thread>
#include <atomic>
#include <mutex>


class Node {
  public:
    Node(size_t nrules, bool default_prediction, double objective, double equivalent_minority);

    Node(unsigned short id, size_t nrules, bool prediction, bool default_prediction,
         double lower_bound, double objective, Node* parent,
         size_t num_captured, double equivalent_minority);

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

    // Returns pair of prefixes and predictions for the path from this node to the root
    inline std::pair<tracking_vector<unsigned short, DataStruct::Tree>, tracking_vector<bool, DataStruct::Tree> >
      get_prefix_and_predictions();

    inline size_t depth() const;
    inline Node* child(unsigned short idx);
    inline Node* parent() const;
    inline void delete_child(unsigned short idx);
    inline size_t num_children() const;

    inline size_t num_captured() const;
    inline double equivalent_minority() const;

    inline typename std::map<unsigned short, Node*>::iterator children_begin();
    inline typename std::map<unsigned short, Node*>::iterator children_end();
    virtual inline double get_curiosity() {
        return 0.0;
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

    friend class CacheTree;
};

class CuriousNode: public Node {
    public:
        CuriousNode(unsigned short id, size_t nrules, bool prediction, bool default_prediction,
             double lower_bound, double objective, double curiosity, CuriousNode* parent,
             size_t num_captured, double equivalent_minority) : Node(id, nrules, prediction, default_prediction,
                 lower_bound, objective, (Node*) parent, num_captured, equivalent_minority) {
            curiosity_ = curiosity;
        }

        inline double get_curiosity() override;
    protected:
        double curiosity_;
};

size_t nn_helper(Node* node);

class CacheTree {
  public:
    CacheTree() {};
    CacheTree(size_t nsamples, size_t nrules, double c, size_t nthreads,
        rule_t *rules, rule_t *labels, rule_t *minority, int ablation,
	bool calculate_size, char const *type);
    ~CacheTree();

    Node* construct_node(unsigned short new_rule, size_t nrules,
           bool prediction, bool default_prediction,
           double lower_bound, double objective,
           Node* parent, int num_not_captured,
           int nsamples, int len_prefix, double c, double equivalent_minority);

    inline double min_objective() const;
    inline tracking_vector<unsigned short, DataStruct::Tree> opt_rulelist() const;
    inline tracking_vector<bool, DataStruct::Tree> opt_predictions() const;

    inline size_t num_nodes() const;
    inline size_t num_nodes(size_t thread_id);
    inline size_t num_evaluated() const;
    inline size_t num_threads() const;
    inline rule_t rule(unsigned short idx) const;
    inline char* rule_features(unsigned short idx) const;
    inline rule_t label(unsigned short idx) const;
    inline rule_t minority(unsigned short idx) const;
    inline bool has_minority() const;
    inline size_t nsamples() const;
    inline size_t nrules() const;
    inline double c() const;
    inline Node* root() const;

    bool update_min_objective(double objective);
    void update_opt_rulelist(tracking_vector<unsigned short, DataStruct::Tree>& parent_prefix,
                             unsigned short new_rule_id);
    void update_opt_predictions(Node* parent, bool new_pred, bool new_default_pred);

    inline void increment_num_evaluated();
    inline void decrement_num_nodes();
    inline int ablation() const;
    inline bool calculate_size() const;
    inline std::vector<unsigned short> get_subrange(size_t i);

    inline std::vector<unsigned short> rule_perm();
    void insert_root();
    void insert(Node* node, size_t thread_id);
    void prune_up(Node* node);
    void garbage_collect(size_t thread_id);
    void print_tree();
    void close_print_file();
    void open_print_file(size_t thread_num, size_t num_threads);
    void play_with_rules();
    Node* check_prefix(tracking_vector<unsigned short, DataStruct::Tree>& prefix);

    inline void lock(size_t thread_id);
    inline void unlock(size_t thread_id);

  protected:
    std::ofstream t_;
    Node* root_;
    size_t nsamples_;
    size_t nrules_;
    size_t nthreads_;
    double c_;

    size_t num_nodes_;
    size_t num_evaluated_;
    int ablation_; // Used to remove support (1) or lookahead (2) bounds
    bool calculate_size_;

    std::atomic<double> min_objective_;
    tracking_vector<unsigned short, DataStruct::Tree> opt_rulelist_;
    std::vector<bool, track_alloc<bool, DataStruct::Tree> > opt_predictions_;
    std::vector<unsigned short> rule_perm_;
    std::vector<pair<unsigned short,unsigned short>> ranges_;

    rule_t *rules_;
    rule_t *labels_;
    rule_t *minority_;

    char const *type_;
    void gc_helper(Node* node, size_t thread_id);
    std::vector<std::mutex> tree_lks_;
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

inline bool Node::done() const{
    return done_;
}

inline void Node::set_done() {
    done_ = 1;
}

inline bool Node::deleted() const{
    return deleted_;
}

inline void Node::set_deleted() {
    deleted_ = 1;
}

inline std::pair<tracking_vector<unsigned short, DataStruct::Tree>, tracking_vector<bool, DataStruct::Tree> >
    Node::get_prefix_and_predictions() {
    tracking_vector<unsigned short, DataStruct::Tree> prefix;
    tracking_vector<bool, DataStruct::Tree> predictions;
    tracking_vector<unsigned short, DataStruct::Tree>::iterator it1 = prefix.begin();
    tracking_vector<bool, DataStruct::Tree>::iterator it2 = predictions.begin();
    Node* node = this;
    for(size_t i = depth_; i > 0; --i) {
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

inline size_t Node::num_children() const {
    return children_.size();
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

inline double CacheTree::min_objective() const {
    return min_objective_.load(std::memory_order_seq_cst);
}

inline tracking_vector<unsigned short, DataStruct::Tree> CacheTree::opt_rulelist() const {
    return opt_rulelist_;
}

inline tracking_vector<bool, DataStruct::Tree> CacheTree::opt_predictions() const {
    return opt_predictions_;
}

inline size_t CacheTree::num_nodes() const {
    return num_nodes_;
}

inline size_t CacheTree::num_nodes(size_t thread_id) {
    std::vector<unsigned short> range = get_subrange(thread_id);
    size_t ret = 1;
    for (std::vector<unsigned short>::iterator it = range.begin(); it != range.end(); ++it) {
        ret += nn_helper(root_->child(*it));
    }
    return ret;
}

inline size_t CacheTree::num_threads() const {
    return nthreads_;
}

inline size_t CacheTree::num_evaluated() const {
    return num_evaluated_;
}

inline rule_t CacheTree::rule(unsigned short idx) const{
    return rules_[idx];
}

inline char* CacheTree::rule_features(unsigned short idx) const{
    return rules_[idx].features;
}

inline rule_t CacheTree::label(unsigned short idx) const{
    return labels_[idx];
}

inline rule_t CacheTree::minority(unsigned short idx) const{
    return minority_[idx];
}

inline bool CacheTree::has_minority() const {
    return minority_ != NULL;
}

inline size_t CacheTree::nsamples() const {
    return nsamples_;
}

inline size_t CacheTree::nrules() const {
    return nrules_;
}

inline double CacheTree::c() const {
    return c_;
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
 * Update the minimum objective of the tree.
 */
/*inline void CacheTree::update_min_objective(double objective) {
    min_objective_ = objective;
    logger->setTreeMinObj(objective);
}*/

/*
 * Update the minimum objective of the tree.
 * Return true if the operation succeeded (i.e., min_objective doesn't get updated to a lower value than objective by another thread during operation)
 * Return false otherwise
 */
inline bool CacheTree::update_min_objective(double objective) {
    bool succeeded = false;
    double min_objective_load = min_objective_.load(std::memory_order_seq_cst);
    while (min_objective_load > objective) {
        min_objective_load = min_objective_.load(std::memory_order_seq_cst);
        succeeded = min_objective_.compare_exchange_strong(min_objective_load, objective, std::memory_order_seq_cst);
        if (succeeded) {
            logger->setTreeMinObj(objective);
            break;
        }
    }
    return succeeded;
}

/*
 * Update the optimal rulelist of the tree.
 */
inline void
CacheTree::update_opt_rulelist(tracking_vector<unsigned short, DataStruct::Tree>& parent_prefix,
                                  unsigned short new_rule_id) {
    opt_rulelist_.assign(parent_prefix.begin(), parent_prefix.end());
    opt_rulelist_.push_back(new_rule_id);
    logger->setTreePrefixLen(opt_rulelist_.size());
}

/*
 * Update the optimal rulelist predictions of the tree.
 */
inline void
CacheTree::update_opt_predictions(Node* parent, bool new_pred, bool new_default_pred) {
    tracking_vector<bool, DataStruct::Tree> predictions;
    Node* node = parent;
    for(size_t i = parent->depth(); i > 0; --i) {
        predictions.push_back(node->prediction());
        node = node->parent();
    }
    std::reverse(predictions.begin(), predictions.end());
    opt_predictions_.assign(predictions.begin(), predictions.end());
    opt_predictions_.push_back(new_pred);
    opt_predictions_.push_back(new_default_pred);
}

/*
 * Increment number of nodes evaluated after performing incremental computation
 * in evaluate_children.
 */
inline void CacheTree::increment_num_evaluated() {
    ++num_evaluated_;
    logger->setTreeNumEvaluated(num_evaluated_);
}

inline std::vector<unsigned short> CacheTree::get_subrange(size_t i) {
	size_t start, end;
	start = ranges_[i].first;
	end = ranges_[i].second;
	return std::vector<unsigned short>(rule_perm_.begin() + start,
	        rule_perm_.begin() + end);
}

inline std::vector<unsigned short> CacheTree::rule_perm() {
    return rule_perm_;
}

/*
 * Called whenever a node is deleted from the tree.
 */
inline void CacheTree::decrement_num_nodes() {
    --num_nodes_;
    logger->setTreeNumNodes(num_nodes_);
}

inline void CacheTree::lock(size_t thread_id) {
  tree_lks_[thread_id].lock();
}

inline void CacheTree::unlock(size_t thread_id) {
  tree_lks_[thread_id].unlock();
}

void delete_interior(CacheTree* tree, Node* node, bool destructive, bool update_remaining_state_space);
void delete_subtree(CacheTree* tree, Node* node, bool destructive, bool update_remaining_state_space, size_t thread_id);
