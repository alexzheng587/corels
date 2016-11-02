//#include "cache.hh"
#include "bbound.hh"
#include <memory>
#include <vector>
#include <stdlib.h>

template<class T>
Node<T>::Node(size_t nrules, bool default_prediction, double objective)
    : id_(0), default_prediction_(default_prediction),
      lower_bound_(0.), objective_(objective), done_(0), deleted_(0), depth_(0), storage_(0) {
}

template<class T>
Node<T>::Node(size_t id, size_t nrules, bool prediction,
           bool default_prediction, double lower_bound,
           double objective, T storage, Node<T>* parent)
    : id_(id), prediction_(prediction), default_prediction_(default_prediction),
      lower_bound_(lower_bound), objective_(objective),
      done_(0), deleted_(0), depth_(1 + parent->depth_), parent_(parent), storage_(storage) {
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
    if (root_)
        delete_subtree<N>(this, root_, true);
    printf("num_nodes: %zu\n", num_nodes_);
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
void CacheTree<N>::insert(N* node) {
    node->parent()->children_.insert(std::make_pair(node->id(), node));
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
N* CacheTree<N>::check_prefix(std::vector<size_t> prefix) {
    typename std::map<size_t, N*>::iterator child_iter;
    N* node = this->root_;
    for(std::vector<size_t>::iterator it = prefix.begin(); it != prefix.end(); ++it) {
        child_iter = node->children_.find(*it);
        if (child_iter != node->children_.end())
            node = node->child(*it);
        else
            return NULL;
    }
    return node;
}

template<class N>
void CacheTree<N>::gc_helper(N* node) {
    N* child;
    std::vector<N*> children;
    for (typename std::map<size_t, N*>::iterator cit = node->children_.begin(); cit != node->children_.end(); ++cit)
        children.push_back(cit->second);
    for (typename std::vector<N*>::iterator cit = children.begin(); cit != children.end(); ++cit) {
        child = *cit;
        if ((child->lower_bound() + c_) >= min_objective_) {
            node->delete_child(child->id());
            delete_subtree<N>(this, child, false);
        } else
            gc_helper(child);
    }
}

template<class N>
void CacheTree<N>::garbage_collect() {
    gc_helper(root_);
}

template class Node<bool>; // BaseNode

template class Node<double>; // CuriousNode

template class CacheTree<BaseNode>;

template class CacheTree<CuriousNode>;
