#pragma once

#include "utils.hh"
#include "alloc.hh"
#include "rule.h"
#include <atomic>
#include <condition_variable>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <thread>
#include <stdlib.h>
#include <vector>

extern std::mutex min_obj_lk;

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
	bool calculate_size, char const *type, size_t random_seed);
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
    inline size_t num_nodes(unsigned short thread_id);
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

    bool update_obj_and_list(double objective, tracking_vector<unsigned short, DataStruct::Tree>& parent_prefix,
                                unsigned short new_rule_id,
                                Node* parent, bool new_pred, bool new_default_pred);
    void update_opt_rulelist(tracking_vector<unsigned short, DataStruct::Tree>& parent_prefix,
                             unsigned short new_rule_id);
    void update_opt_predictions(Node* parent, bool new_pred, bool new_default_pred);

    inline void increment_num_evaluated();
    inline void decrement_num_nodes();
    inline int ablation() const;
    inline bool calculate_size() const;
    inline tracking_vector<unsigned short, DataStruct::Tree>* get_subrange(size_t i);
    inline tracking_vector<unsigned short, DataStruct::Tree>* rule_perm();
    inline std::vector<tracking_vector<unsigned short, DataStruct::Tree>* > split_rules(size_t n);

    void insert_root();
    void insert(Node* node, unsigned short thread_id);
    void prune_up(Node* node);
    void garbage_collect(unsigned short thread_id);
    void print_tree();
    void close_print_file();
    void open_print_file(size_t thread_num, size_t num_threads);
    void play_with_rules();
    Node* check_prefix(tracking_vector<unsigned short, DataStruct::Tree>& prefix);

    inline void lock(unsigned short thread_id);
    inline void unlock(unsigned short thread_id);


    inline unsigned short num_inactive_threads() const;
    inline void increment_num_inactive_threads();
    inline void decrement_num_inactive_threads();
    inline void lock_inactive_thread_lk();
    inline void unlock_inactive_thread_lk();
    inline size_t n_acc() const;

    inline bool done() const;
    inline void set_done(bool is_done);

    inline void thread_wait();
    inline void wake_all_inactive();
    inline void wake_n_inactive(size_t n);

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

    double min_objective_;
    tracking_vector<unsigned short, DataStruct::Tree> opt_rulelist_;
    tracking_vector<bool, DataStruct::Tree> opt_predictions_;
    tracking_vector<unsigned short, DataStruct::Tree> rule_perm_;
    tracking_vector<pair<unsigned short, unsigned short>, DataStruct::Tree> ranges_;

    rule_t *rules_;
    rule_t *labels_;
    rule_t *minority_;

    unsigned short inactive_threads_;
    bool done_;

    std::mutex inactive_thread_lk_;
    std::condition_variable inactive_thread_cv_;
    std::mutex root_lk_;
    size_t n_acc_; 

    char const *type_;
    void gc_helper(Node* node, unsigned short thread_id);
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
    typename std::map<unsigned short, Node*>::iterator iter = children_.find(idx);
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
    for(auto child_it = children_begin(); child_it != children_end(); ++child_it) {
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

inline double CacheTree::min_objective() const {
    return min_objective_;
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

inline size_t CacheTree::num_nodes(unsigned short thread_id) {
    tracking_vector<unsigned short, DataStruct::Tree>* range = get_subrange(thread_id);
    size_t ret = 1;
    for (tracking_vector<unsigned short, DataStruct::Tree>::iterator it = range->begin(); it != range->end(); ++it) {
        //std::cout << *it << std::endl;
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

inline bool
CacheTree::update_obj_and_list(double objective, tracking_vector<unsigned short, DataStruct::Tree>& parent_prefix,
                                unsigned short new_rule_id,
                                Node* parent, bool new_pred, bool new_default_pred) {
  min_obj_lk.lock();
  if (objective >= min_objective_) {
    min_obj_lk.unlock();
    return false;
  } else {
    min_objective_ = objective;
    logger->setTreeMinObj(objective);
    this->update_opt_rulelist(parent_prefix, new_rule_id);
    this->update_opt_predictions(parent, new_pred, new_default_pred);
  }
  min_obj_lk.unlock();
  return true;
}

/*
 * Update the minimum objective of the tree.
 */
/*inline void CacheTree::update_min_objective(double objective) {
    min_objective_ = objective;
    logger->setTreeMinObj(objective);
}*/

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

inline tracking_vector<unsigned short, DataStruct::Tree>* CacheTree::get_subrange(size_t i) {
	size_t start, end;
	start = ranges_[i].first;
	end = ranges_[i].second;
	tracking_vector<unsigned short, DataStruct::Tree>* subrange = new tracking_vector<unsigned short, DataStruct::Tree>;
    subrange->assign(rule_perm_.begin() + start,rule_perm_.begin() + end);
    return subrange;
}

inline tracking_vector<unsigned short, DataStruct::Tree>* CacheTree::rule_perm() {
    return &rule_perm_;
}


inline std::vector<tracking_vector<unsigned short, DataStruct::Tree>* > CacheTree::split_rules(size_t n) {
    std::vector<tracking_vector<unsigned short, DataStruct::Tree>* > splits;
    /*size_t nrules = rule_perm_.size();
    size_t split_size = (nrules + n - 1) / n;
    for(size_t i = 0; i < n; ++i){
        tracking_vector<unsigned short, DataStruct::Tree> split(rule_perm_.begin() + split_size * i, rule_perm_.begin() + split_size * (i + 1));
        splits.push_back(split);
    }
    return splits;*/

    tracking_vector<unsigned short, DataStruct::Tree> rule_perm = rule_perm_;
    std::iota(rule_perm.begin(), rule_perm.end(), 1);
    std::random_shuffle(rule_perm.begin(), rule_perm.end());

    unsigned short rules_per_thread = (nrules_-1) / n;
    unsigned short inc = (nrules_ - 1) - (rules_per_thread * n);
    unsigned short sIdx = 0;

    for(size_t i = 0; i < n; ++i) {
        unsigned short eIdx = sIdx + rules_per_thread + (i < inc ? 1 : 0);
        // printf("START INDEX: %zu, END INDEX: %zu\n", sIdx, eIdx);
        tracking_vector<unsigned short, DataStruct::Tree>* split = new tracking_vector<unsigned short, DataStruct::Tree>;
        split->assign(rule_perm.begin() + sIdx, rule_perm.begin() + eIdx);
        splits.push_back(split);
        sIdx = eIdx;
    }
    return splits;
}



/*
 * Called whenever a node is deleted from the tree.
 */
inline void CacheTree::decrement_num_nodes() {
    --num_nodes_;
    logger->setTreeNumNodes(num_nodes_);
}

inline void CacheTree::lock(unsigned short thread_id) {
  tree_lks_[thread_id].lock();
}

inline void CacheTree::unlock(unsigned short thread_id) {
  tree_lks_[thread_id].unlock();
}

inline unsigned short CacheTree::num_inactive_threads() const {
    return inactive_threads_;
}

inline void CacheTree::increment_num_inactive_threads() {
    inactive_threads_++;
}

inline void CacheTree::decrement_num_inactive_threads() {
    inactive_threads_--;
}

inline bool CacheTree::done() const {
    return done_;
}

inline void CacheTree::set_done(bool is_done) {
    done_ = is_done;
}

// NOTE: The thread will have the inactive_thread_lk_ when it wakes up
// Also, assume the calling thread has NOT already acquired the inactive_thread_lk_ before claling
inline void CacheTree::thread_wait() {
    std::unique_lock<std::mutex> inactive_thread_lk(inactive_thread_lk_);
    inactive_thread_cv_.wait(inactive_thread_lk);
    inactive_thread_lk.unlock();
}

inline void CacheTree::wake_all_inactive() {
    inactive_thread_cv_.notify_all();
}

inline void CacheTree::wake_n_inactive(size_t n) {
    for(size_t i = 0; i < n; ++i)
        inactive_thread_cv_.notify_one();
}

inline void CacheTree::lock_inactive_thread_lk() {
    n_acc_++;
    inactive_thread_lk_.lock();
}

inline void CacheTree::unlock_inactive_thread_lk() {
    inactive_thread_lk_.unlock();
}

inline size_t CacheTree::n_acc() const {
    return n_acc_;
}


void delete_interior(CacheTree* tree, Node* node, bool destructive, bool update_remaining_state_space);
void delete_subtree(CacheTree* tree, Node* node, bool destructive, bool update_remaining_state_space, unsigned short thread_id);
