//#include "cache.hh"
#include "bbound.hh"
#include "utils.hh"
#include <memory>
#include <vector>
#include <stdlib.h>

template<class T>
Node<T>::Node(size_t nrules, bool default_prediction, double objective)
    : id_(0), default_prediction_(default_prediction),
      lower_bound_(0.), objective_(objective), done_(0), deleted_(0), depth_(0),
      storage_(0), num_captured_(0) {
}

template<class T>
Node<T>::Node(unsigned short id, size_t nrules, bool prediction,
              bool default_prediction, double lower_bound,
              double objective, T storage, Node<T>* parent, size_t num_captured)
    : id_(id), prediction_(prediction), default_prediction_(default_prediction),
      lower_bound_(lower_bound), objective_(objective),
      done_(0), deleted_(0), depth_(1 + parent->depth_), parent_(parent),
      storage_(storage), num_captured_(num_captured) {
}

template<class N>
CacheTree<N>::CacheTree(size_t nsamples, size_t nrules, double c, rule_t *rules, rule_t *labels)
    : root_(0), nsamples_(nsamples), nrules_(nrules), c_(c),
      num_nodes_(0), num_evaluated_(0), min_objective_(0.5),
      opt_rulelist_({}), opt_predictions_({}) {
    opt_rulelist_.resize(0);
    opt_predictions_.resize(0);
    rules_.resize(nrules);
    labels_.resize(2);
    size_t i;
    for (i = 0; i < nrules; i++)
        rules_[i] = rules[i];
    labels_[0] = labels[0];
    labels_[1] = labels[1];

    logger.setTreeMinObj(min_objective_);
    logger.setTreeNumNodes(num_nodes_);
    logger.setTreeNumEvaluated(num_evaluated_);
}

template<class N>
CacheTree<N>::~CacheTree() {
    if (root_)
        delete_subtree<N>(this, root_, true, false);
}

template<class N>
void CacheTree<N>::insert_root() {
    VECTOR tmp_vec;
    size_t d0, d1;
    bool default_prediction;
    double objective;
    make_default(&tmp_vec, nsamples_);
    d0 = labels_[0].support;
    d1 = nsamples_ - d0;
    if (d0 > d1) {
        default_prediction = 0;
        objective = (float)(d1) / nsamples_;
    } else {
        default_prediction = 1;
        objective = (float)(d0) / nsamples_;
    }
    root_ = new N(nrules_, default_prediction, objective);
    min_objective_ = objective;
    logger.setTreeMinObj(objective);
    ++num_nodes_;
    logger.setTreeNumNodes(num_nodes_);
    opt_predictions_.push_back(default_prediction);
    logger.setTreePrefixLen(0);
}

template<class N>
void CacheTree<N>::insert(N* node) {
    node->parent()->children_.insert(std::make_pair(node->id(), node));
    ++num_nodes_;
    logger.setTreeNumNodes(num_nodes_);
}

template<class N>
void CacheTree<N>::prune_up(N* node) {
    unsigned short id;
    size_t depth = node->depth();
    N* parent;
    while (node->children_.size() == 0) {
        if (depth > 0) {
            id = node->id();
            parent = node->parent();
            parent->children_.erase(id);
            --num_nodes_;
            node->set_deleted();
            node = parent;
            --depth;
        } else {
            --num_nodes_;
            break;
        }
    }
    logger.setTreeNumNodes(num_nodes_);
}

template<class N>
N* CacheTree<N>::check_prefix(std::vector<unsigned short>& prefix) {
    N* node = this->root_;
    for(std::vector<unsigned short>::iterator it = prefix.begin(); it != prefix.end(); ++it) {
        node = node->child(*it);
        if (node == NULL)
            return NULL;
    }
    return node;
}

template<class N>
void CacheTree<N>::gc_helper(N* node) {
    if (!node->done())
        logger.addQueueElement(node->depth(), node->lower_bound());
    N* child;
    std::vector<N*> children;
    for (typename std::map<unsigned short, N*>::iterator cit = node->children_.begin(); cit != node->children_.end(); ++cit)
        children.push_back(cit->second);
    for (typename std::vector<N*>::iterator cit = children.begin(); cit != children.end(); ++cit) {
        child = *cit;
        if ((child->lower_bound() + c_) >= min_objective_) {
            node->delete_child(child->id());
            delete_subtree<N>(this, child, false, false);
        } else
            gc_helper(child);
    }
}

template<class N>
void CacheTree<N>::garbage_collect() {
    logger.clearRemainingSpaceSize();
    gc_helper(root_);
}

template<class N>
inline void CacheTree<N>::update_min_objective(double objective) {
    min_objective_ = objective;
    logger.setTreeMinObj(objective);
}

template<class N>
inline void
CacheTree<N>::update_opt_rulelist(std::vector<unsigned short>& parent_prefix,
                                  unsigned short new_rule_id) {
    opt_rulelist_.assign(parent_prefix.begin(), parent_prefix.end());
    opt_rulelist_.push_back(new_rule_id);
    logger.setTreePrefixLen(opt_rulelist_.size());
}

template<class N>
inline void
CacheTree<N>::update_opt_predictions(std::vector<bool>& parent_predictions,
                                     bool new_pred,
                                     bool new_default_pred) {
    opt_predictions_.assign(parent_predictions.begin(), parent_predictions.end());
    opt_predictions_.push_back(new_pred);
    opt_predictions_.push_back(new_default_pred);
}

template<class N>
inline void CacheTree<N>::increment_num_evaluated() {
    ++num_evaluated_;
    logger.setTreeNumEvaluated(num_evaluated_);
}

template<class N>
inline void CacheTree<N>::decrement_num_nodes() {
    --num_nodes_;
    logger.setTreeNumNodes(num_nodes_);
}


template class Node<bool>; // BaseNode

template class Node<double>; // CuriousNode

template class CacheTree<BaseNode>;

template class CacheTree<CuriousNode>;
