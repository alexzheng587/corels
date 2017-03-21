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
#include <scoped_allocator>

//template <class T> class Node;
//template <class N> class CacheTree;

//typedef Node<bool> BaseNode;       // nothing extra
//typedef Node<double> CuriousNode;  // curiosity

/*
template <class T> 
struct cache_alloc : track_alloc<T> { 
	typedef T value_type;
	T* allocate (size_t n) {
		logger.addToTreeMemory(n * sizeof(T));
		return static_cast<T*>(malloc(n*sizeof(T)));
	}
	void deallocate (T* p, size_t n) {
		logger.removeFromTreeMemory(n * sizeof(*p));
		free(p);
	}
};
*/

//template <class T>
class Node {
  public:
    explicit Node(size_t nrules, bool default_prediction, double objective, double minority);

    Node(unsigned short id, size_t nrules, bool prediction, bool default_prediction,
         double lower_bound, double objective, Node* parent,
         size_t num_captured, double minority);

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
    inline std::pair<std::vector<unsigned short, cache_alloc<unsigned short> >, std::vector<bool, cache_alloc<bool> > >
        get_prefix_and_predictions();

    inline size_t depth() const;
    inline Node* child(unsigned short idx);
    inline Node* parent() const;
    inline void delete_child(unsigned short idx);
    inline size_t num_children() const;

    inline size_t num_captured() const;
    inline double minority() const;

    inline typename std::map<unsigned short, Node*>::iterator children_begin();
    inline typename std::map<unsigned short, Node*>::iterator children_end();
    inline Node* random_child(); // FIXME
    virtual inline double get_storage() {
        return 0.0;
    }

  protected:
    std::map<unsigned short, Node*, std::less<unsigned short>, cache_alloc<std::pair<unsigned short, Node*> > > children_;
    Node* parent_;
    double lower_bound_;
    double objective_;
    double minority_;
    size_t depth_;
    size_t num_captured_;
    unsigned short id_;
    bool prediction_;
    bool default_prediction_;
    bool done_;
    bool deleted_;

    // space for something extra, like curiosity or a bit vector
    //T storage_;

    friend class CacheTree;//<Node<T> >;
};

class CuriousNode: public Node {
    public:
        CuriousNode(unsigned short id, size_t nrules, bool prediction, bool default_prediction,
             double lower_bound, double objective, double storage, CuriousNode* parent,
             size_t num_captured, double minority) : Node(id, nrules, prediction, default_prediction,
                 lower_bound, objective, parent, num_captured, minority) {
            storage_ = storage;
            /*return new Node(id, nrules, prediction, default_prediction, lower_bound, objective,
                    parent, num_captured, minority);
                    */
        }

        inline double get_storage();
    protected:
        double storage_;
};

//template<class N>
class CacheTree {
  public:
    CacheTree() {};
    CacheTree(size_t nsamples, size_t nrules, double c, rule_t *rules,
              rule_t *labels, rule_t *meta);
    ~CacheTree();

    CuriousNode* construct_node(unsigned short new_rule, size_t nrules,
           bool prediction, bool default_prediction,
           double lower_bound, double objective,
           Node* parent, int num_not_captured,
           int nsamples, int len_prefix, double c, double minority);

    inline double min_objective() const;
    inline std::vector<unsigned short, cache_alloc<unsigned short> > opt_rulelist() const;
    inline std::vector<bool, cache_alloc<bool> > opt_predictions() const;

    inline size_t num_nodes() const;
    inline size_t num_evaluated() const;
    inline rule_t rule(unsigned short idx) const;
    inline char* rule_features(unsigned short idx) const;
    inline rule_t label(unsigned short idx) const;
    inline rule_t meta(unsigned short idx) const;
    inline size_t meta_size() const;
    inline size_t nsamples() const;
    inline size_t nrules() const;
    inline double c() const;
    inline Node* root() const;

    void update_min_objective(double objective);
    void update_opt_rulelist(std::vector<unsigned short, cache_alloc<unsigned short> >& parent_prefix,
                             unsigned short new_rule_id);
    void update_opt_predictions(std::vector<bool, cache_alloc<bool> >& parent_predictions,
                                bool new_pred,
                                bool new_default_pred);

    inline void increment_num_evaluated();
    inline void decrement_num_nodes();

    void insert_root();
    void insert(Node* node);
    void prune_up(Node* node);
    void garbage_collect();
    void print_tree();
    void open_print_file(size_t thread_num, size_t num_threads);
    void play_with_rules();
    Node* check_prefix(std::vector<unsigned short, cache_alloc<unsigned short> >& prefix);

  protected:
    std::ofstream t_;
    Node* root_;
    size_t nsamples_;
    size_t nrules_;
    double c_;

    size_t num_nodes_;
    size_t num_evaluated_;

    double min_objective_;
    std::vector<unsigned short, cache_alloc<unsigned short> > opt_rulelist_;
    std::vector<bool, cache_alloc<bool> > opt_predictions_;

    std::vector<rule_t> rules_;
    std::vector<rule_t> labels_;
    std::vector<rule_t> meta_;

    void gc_helper(Node* node);
    void print_tree_helper(Node* node, std::vector<short>& rlist);
    void close_print_file();
};

/*
class CuriousCacheTree : public CacheTree {
    public:
//        CuriousCacheTree() {};
        CuriousCacheTree(size_t nsamples, size_t nrules, double c, rule_t *rules,
                                rule_t *labels, rule_t *meta) : CacheTree(nsamples,
                                    nrules, c, rules, labels, meta) {
            curiousity = true;
        }

        CuriousCacheTree(size_t nsamples, size_t nrules, double c, rule_t *rules,
                                rule_t *labels, rule_t *meta) {
			root_ = 0;
            num_nodes_ = 0;
            num_evaluated_ = 0;
			nsamples_ = nsamples;
            nrules_ = nrules;
            c_ = c;
            min_objective_ = 0.5;
            opt_rulelist_ = {};
            opt_predictions_ = {};
			opt_rulelist_.resize(0);
			opt_predictions_.resize(0);
			rules_.resize(nrules);
			labels_.resize(2);
			size_t i;
			for (i = 0; i < nrules; i++) {
			  rules_[i] = rules[i];
			}
			labels_[0] = labels[0];
			labels_[1] = labels[1];
			if (meta) {
			  meta_.resize(1);
			  meta_[0] = meta[0];
			} else {
			  meta_.resize(0);
			}
			//logger.setTreeMinObj(min_objective_);
			//logger.setTreeNumNodes(num_nodes_);
			//logger.setTreeNumEvaluated(num_evaluated_);
        }
        Node* construct_node(unsigned short new_rule, size_t nrules,
           bool prediction, bool default_prediction,
           double lower_bound, double objective,
           Node* parent, int num_not_captured,
           int nsamples, int len_prefix, double c, double minority) override;

        CuriousNode* construct_node(unsigned short new_rule, size_t nrules,
           bool prediction, bool default_prediction,
           double lower_bound, double objective,
           CuriousNode* parent, int num_not_captured,
           int nsamples, int len_prefix, double c, double minority);
    protected:
        bool curiousity;
};
           */

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

inline std::pair<std::vector<unsigned short, cache_alloc<unsigned short> >, std::vector<bool, cache_alloc<bool> > >
    Node::get_prefix_and_predictions() {
    std::vector<unsigned short, cache_alloc<unsigned short> > prefix;
    std::vector<bool, cache_alloc<bool> > predictions;
    auto it1 = prefix.begin();
    auto it2 = predictions.begin();
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

inline Node* Node::random_child() {
    typename std::map<unsigned short, Node*>::iterator iter;
    unsigned short idx;
    iter = children_.begin();
    idx = rand() % (children_.size());
    std::advance(iter, idx);
    return iter->second;
}

inline Node* Node::parent() const {
    return parent_;
}

inline size_t Node::num_captured() const {
    return num_captured_;
}

inline double Node::minority() const {
    return minority_;
}

inline double CuriousNode::get_storage() {
    return storage_;
}


inline double CacheTree::min_objective() const {
    return min_objective_;
}

inline std::vector<unsigned short, cache_alloc<unsigned short> > CacheTree::opt_rulelist() const {
    return opt_rulelist_;
}

inline std::vector<bool, cache_alloc<bool> > CacheTree::opt_predictions() const {
    return opt_predictions_;
}

inline size_t CacheTree::num_nodes() const {
    return num_nodes_;
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

inline rule_t CacheTree::meta(unsigned short idx) const{
    return meta_[idx];
}

inline size_t CacheTree::meta_size() const {
    return meta_.size();
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

/*
 * Update the minimum objective of the tree.
 */
inline void CacheTree::update_min_objective(double objective) {
    min_objective_ = objective;
    //logger.setTreeMinObj(objective);
}

/*
 * Update the optimal rulelist of the tree.
 */
inline void
CacheTree::update_opt_rulelist(std::vector<unsigned short, cache_alloc<unsigned short> >& parent_prefix,
                                  unsigned short new_rule_id) {
    opt_rulelist_.assign(parent_prefix.begin(), parent_prefix.end());
    opt_rulelist_.push_back(new_rule_id);
    //logger.setTreePrefixLen(opt_rulelist_.size());
}

/*
 * Update the optimal rulelist predictions of the tree.
 */
inline void
CacheTree::update_opt_predictions(std::vector<bool, cache_alloc<bool> >& parent_predictions,
                                     bool new_pred,
                                     bool new_default_pred) {
    opt_predictions_.assign(parent_predictions.begin(), parent_predictions.end());
    opt_predictions_.push_back(new_pred);
    opt_predictions_.push_back(new_default_pred);
}

/*
 * Increment number of nodes evaluated after performing incremental computation 
 * in evaluate_children.
 */
inline void CacheTree::increment_num_evaluated() {
    ++num_evaluated_;
    //logger.setTreeNumEvaluated(num_evaluated_);
}

/*
 * Called whenever a node is deleted from the tree.
 */
inline void CacheTree::decrement_num_nodes() {
    --num_nodes_;
    //logger.setTreeNumNodes(num_nodes_);
}

void delete_subtree(CacheTree* tree, Node* node, bool destructive, bool update_remaining_state_space);

/*template<class N>
using construct_signature = N* (*)(unsigned short, size_t, bool, bool, double, double,
                                   N* parent, int, int, int, double, double);

BaseNode* base_construct_policy(unsigned short new_rule, size_t nrules,
                                bool prediction, bool default_prediction,
                                double lower_bound, double objective,
                                BaseNode* parent, int num_not_captured,
                                int nsamples, int len_prefix, double c, double minority);

CuriousNode* curious_construct_policy(unsigned short new_rule, size_t nrules,
                                      bool prediction, bool default_prediction,
                                      double lower_bound, double objective,
                                      CuriousNode* parent, int num_not_captured,
                                      int nsamples, int len_prefix, double c, double minority);
                                      */
