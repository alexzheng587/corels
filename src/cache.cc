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
      objective_(objective), done_(0), parent_(parent){
    parent_->children_.insert(std::make_pair(id, this));
}

CacheTree::CacheTree(size_t nsamples, size_t nrules, double c, rule_t *rules, rule_t *labels)
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

CacheTree::~CacheTree() {
    delete_subtree(root_);
}

void CacheTree::insert_root() {
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
    root_ = new CacheNode(nrules_, default_prediction, objective);
    min_objective_ = objective;
    ++num_nodes_;
}

void CacheTree::insert(size_t new_rule, bool prediction, bool default_prediction,
                       double lower_bound, double objective, CacheNode* parent) {
    CacheNode* child = new CacheNode(new_rule, nrules_, prediction, default_prediction,
                                     lower_bound, objective, parent);
    parent->children_.insert(std::make_pair(new_rule, child));
    ++num_nodes_;
}

void CacheTree::evaluate_children(CacheNode* parent, VECTOR parent_not_captured) {
    VECTOR captured, captured_zeros, not_captured, not_captured_zeros;
    int num_captured, c0, c1, captured_correct;
    int num_not_captured, d0, d1, default_correct;
    bool prediction, default_prediction;
    double lower_bound, objective, parent_lower_bound;
    rule_vinit(nsamples_, &captured);
    rule_vinit(nsamples_, &captured_zeros);
    rule_vinit(nsamples_, &not_captured);
    rule_vinit(nsamples_, &not_captured_zeros);
    size_t i, len_prefix;
    len_prefix = parent->depth() + 1;
    parent_lower_bound = parent->lower_bound();
    for (i = 1; i < nrules_; i++) {
        rule_vand(captured, parent_not_captured, rules_[i].truthtable, nsamples_, &num_captured);
        rule_vand(captured_zeros, captured, labels_[0].truthtable, nsamples_, &c0);
        c1 = num_captured - c0;
        if (c0 > c1) {
            prediction = 0;
            captured_correct = c0;
        } else {
            prediction = 1;
            captured_correct = c1;
        }
        lower_bound = parent_lower_bound + (float)(num_captured - captured_correct) / nsamples_ + c_;
        rule_vandnot(not_captured, parent_not_captured, captured, nsamples_, &num_not_captured);
        rule_vand(not_captured_zeros, not_captured, labels_[0].truthtable, nsamples_, &d0);
        d1 = num_not_captured - d0;
        if (d0 > d1) {
            default_prediction = 0;
            default_correct = d0;
        } else {
            default_prediction = 1;
            default_correct = d1;
        }
        objective = lower_bound + (float)(num_not_captured - default_correct) / nsamples_;
        if (objective < min_objective_) {
            printf("min(objective): %1.5f -> %1.5f, length: %zu, cache size: %zu\n",
                   min_objective_, objective, len_prefix, num_nodes_);
            min_objective_ = objective;
        }
        if ((lower_bound + c_) < min_objective_)
            insert(i, prediction, default_prediction, lower_bound, objective, parent);
    }
    if (parent->children_.size() == 0) {
        prune_up(parent);
    } else {
        parent->set_done();
        ++num_evaluated_;
    }
}

CacheNode* CacheTree::stochastic_select(VECTOR not_captured) {
    CacheNode* node = root_;
    rule_copy(not_captured, rules_[root_->id()].truthtable, nsamples_);
    int cnt;
    size_t idx = 0;
    std::map<size_t, node_type*>::iterator iter;
    while (node->done()) {
        if ((node->lower_bound() + c_) > min_objective_) {
            if (node->depth_ > 0) {
                CacheNode* parent = node->parent();
                parent->children_.erase(node->id());
                delete_subtree(node);
            }
            return 0;
        }
        if (node->children_.size() == 0) {
            prune_up(node);
            return 0;
        }
        iter = node->children_.begin();
        idx = rand() % (node->children_.size());
        std::advance(iter, idx);
        node = iter->second;
        rule_vandnot(not_captured, not_captured, rules_[iter->first].truthtable, nsamples_, &cnt);
    }
    return node;
}

void CacheTree::toy(size_t max_num_nodes) {
    CacheNode* node;
    VECTOR not_captured;
    size_t num_iter = 0;
    rule_vinit(nsamples_, &not_captured);
    insert_root();
    while ((num_nodes_ < max_num_nodes) and (num_nodes_ > 0)) {
        node = stochastic_select(not_captured);
        if (node)
            evaluate_children(node, not_captured);
        ++num_iter;
        if ((num_iter % 10000) == 0)
            printf("num_iter: %zu, num_nodes: %zu\n", num_iter, num_nodes_);
    }
}

void CacheTree::prune_up(CacheNode* node) {
    size_t id, depth = node->depth();
    CacheNode* parent;
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

void CacheTree::delete_subtree(CacheNode* node) {
    CacheNode* child;
    std::map<size_t, node_type*>::iterator iter;
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
