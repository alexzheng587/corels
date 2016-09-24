#include "cache.hh"
#include <memory>
#include <vector>
#include <stdlib.h>


CacheNode::CacheNode(size_t nrules, bool default_prediction, double objective)
    : id_(0), depth_(0), default_prediction_(default_prediction),
      lower_bound_(0.), objective_(objective), done_(0) {
}

CacheNode::CacheNode(size_t id, size_t nrules, bool prediction,
                     bool default_prediction, double lower_bound,
                     double objective, CacheNode* parent)
    : id_(id), depth_(1 + parent->depth_), prediction_(prediction),
      default_prediction_(default_prediction), lower_bound_(lower_bound),
      objective_(objective), done_(0), parent_(parent) {
    parent_->children_.insert(std::make_pair(id, this));
}

CacheNodeCurious::CacheNodeCurious(size_t nrules, bool default_prediction, double objective)
    : CacheNode(nrules, default_prediction, objective), curiosity_(0.) {
}

CacheNodeCurious::CacheNodeCurious(size_t id, size_t nrules, bool prediction,
                     bool default_prediction, double lower_bound,
                     double objective, CacheNodeCurious* parent, double curiosity)
    : CacheNode(id, nrules, prediction, default_prediction, lower_bound, objective, parent),
      curiosity_(curiosity) {
}

template<class N>
CacheTree<N>::CacheTree(size_t nsamples, size_t nrules, double c, rule_t *rules, rule_t *labels)
    : root_(0), nsamples_(nsamples), nrules_(nrules), c_(c), min_objective_(1.),
      num_nodes_(0), num_evaluated_(0) {
    rules_.resize(nrules);
    labels_.resize(2);
    size_t i;
    for (i = 0; i < nrules; i++)
        rules_[i] = rules[i];
    labels_[0] = labels[0];
    labels_[1] = labels[1];
}

template<class N>
CacheTree<N>::~CacheTree() {
    delete_subtree(root_);
}

template<class N>
void CacheTree<N>::insert_root() {
    VECTOR tmp_vec;
    size_t d0, d1;
    bool default_prediction;
    double objective;
    make_default(&tmp_vec, nsamples_);
    d0 = labels_[0].support;
    d1 = labels_[1].support;
    if (d0 > d1) {
        default_prediction = 0;
        objective = (float)(d1) / nsamples_;
    } else {
        default_prediction = 1;
        objective = (float)(d0) / nsamples_;
    }
    root_ = new N(nrules_, default_prediction, objective);
    min_objective_ = objective;
    ++num_nodes_;
}

template<class N>
void CacheTree<N>::insert(size_t new_rule, bool prediction, bool default_prediction,
                          double lower_bound, double objective, N* parent) {
    N* child = new N(new_rule, nrules_, prediction, default_prediction,
                     lower_bound, objective, parent);
    parent->children_.insert(std::make_pair(new_rule, child));
    ++num_nodes_;
}

template<class N>
void CacheTree<N>::prune_up(N* node) {
    size_t id, depth = node->depth();
    N* parent;
    while (node->children_.size() == 0) {
        if (depth > 0) {
            id = node->id();
            parent = node->parent();
            parent->children_.erase(id);
            --num_nodes_;
            delete node;
            node = parent;
            --depth;
        } else {
            --num_nodes_;
            break;
        }
    }
}

template<class N>
void CacheTree<N>::delete_subtree(N* node) {
    N* child;
    typename std::map<size_t, N*>::iterator iter;
    if (node->done()) {
        iter = node->children_.begin();
        while (iter != node->children_.end()) {
            child = iter->second;
            delete_subtree(child);
            ++iter;
        }
    }
    --num_nodes_;    
    //printf("delete node %zu at depth %zu (lb=%1.5f, ob=%1.5f) %zu\n",
    //       node->id(), node->depth(), node->lower_bound(), node->objective(), num_nodes_);
    delete node;
}

template class CacheTree<CacheNode>;
