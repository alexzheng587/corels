# ifdef VAL
#include "common.hh"
# endif 
#include "cache.hh"
#include "utils.hh"
#include <memory>
#include <numeric>
#include <vector>
#include <iostream>
#include <stdlib.h>

Node::Node(size_t nrules, bool default_prediction, double objective, double equivalent_minority)
    : lower_bound_(equivalent_minority), objective_(objective), equivalent_minority_(equivalent_minority), depth_(0),
      num_captured_(0), id_(0), default_prediction_(default_prediction),
      done_(0), deleted_(0) {
          (void) nrules;
}

Node::Node(unsigned short id, size_t nrules, bool prediction,
              bool default_prediction, double lower_bound, double objective,
              Node* parent, size_t num_captured, double equivalent_minority)
    : parent_(parent), lower_bound_(lower_bound), objective_(objective),
      equivalent_minority_(equivalent_minority), depth_(1 + parent->depth_),
      num_captured_(num_captured), id_(id), prediction_(prediction),
      default_prediction_(default_prediction), done_(0), deleted_(0) {
          (void) nrules;
}

CacheTree::CacheTree(size_t nsamples, size_t nrules, double c, size_t nthreads,
    rule_t *rules, rule_t *labels, rule_t *minority, int ablation,
    bool calculate_size, char const *type, size_t random_seed)
    : root_(NULL), nsamples_(nsamples), nrules_(nrules), c_(c),
      num_nodes_(0), ablation_(ablation),
      calculate_size_(calculate_size), min_objective_(0.5),
      opt_rulelist_({}), opt_predictions_({}), rule_perm_(nrules - 1),
      ranges_(nthreads), type_(type), inactive_threads_(0), n_acc_(0), done_(false) {
    opt_rulelist_.resize(0);
    opt_predictions_.resize(0);
    rules_ = rules;
    labels_ = labels;
    nthreads_ = nthreads;

    srand(random_seed);
    std::cout << "RANDOM SEED: " << random_seed << std::endl;
    std::iota(rule_perm_.begin(), rule_perm_.end(), 1);
    std::random_shuffle(rule_perm_.begin(), rule_perm_.end());

    tree_lks_ = std::vector<std::mutex> (nthreads);

    unsigned short rules_per_thread = (nrules-1) / nthreads;
    unsigned short inc = (nrules - 1) - (rules_per_thread * nthreads);
    unsigned short sIdx = 0;

    for(size_t i = 0; i < nthreads; ++i) {
        unsigned short eIdx = sIdx + rules_per_thread + (i < inc ? 1 : 0);
        ranges_[i] = std::make_pair(sIdx, eIdx);
        if (logger->getVerbosity() > 10) {
            printf("RANGE INDICES: %hu-%hu\nRANGE: ", sIdx, eIdx);
            for (tracking_vector<unsigned short, DataStruct::Tree>::iterator it = 
                rule_perm_.begin() + sIdx;
                it != rule_perm_.begin() + eIdx; it++) {
                printf("%hu ", *it);
                printf("\n");
            }
        }
        sIdx = eIdx;
    }

    if (minority) {
        minority_ = minority;
    } else {
        minority_ = NULL;
    }

    logger->setTreeMinObj(min_objective_);
    logger->setTreeNumNodes(num_nodes_);
}

CacheTree::~CacheTree() {
	if (t_)
        close_print_file();
	if (root_)
        delete_subtree(this, root_, true, false, 0);
}


void CacheTree::close_print_file() {
    t_.close();
}

void CacheTree::open_print_file(size_t t_num, size_t num_threads) {
    char tname[64];
    sprintf(tname, "tree_%zu_of_%zu.txt", t_num, num_threads);
    t_.open(tname, ios::out | ios::trunc);
}

Node* CacheTree::construct_node(unsigned short new_rule, size_t nrules, bool prediction,
                         bool default_prediction, double lower_bound, double objective, 
                         Node* parent, int num_not_captured, int nsamples,
                         int len_prefix, double c, double equivalent_minority) {
    size_t num_captured = nsamples - num_not_captured;
    Node* n;
    if (strcmp(type_, "curious") == 0) {
        double curiosity = (lower_bound - equivalent_minority) * nsamples / (double)(num_captured);
        n = (Node*) (new CuriousNode(new_rule, nrules, prediction, default_prediction,
                                lower_bound, objective, curiosity, (CuriousNode*) parent, num_captured, equivalent_minority));
    } else {
        n = (new Node(new_rule, nrules, prediction, default_prediction,
                                lower_bound, objective, parent, num_captured, equivalent_minority));
    }
    logger->addToMemory(sizeof(*n), DataStruct::Tree);
    return n;
}

/*
 * Inserts the root of the tree, setting up the default rules.
 */
void CacheTree::insert_root() {
    VECTOR tmp_vec;
    size_t d0, d1;
    bool default_prediction;
    double objective;
    make_default(&tmp_vec, nsamples_);
    d0 = labels_[0].support;
    d1 = nsamples_ - d0;
    if (d0 > d1) {
        default_prediction = 0;
        objective = (double)(d1) / nsamples_;
    } else {
        default_prediction = 1;
        objective = (double)(d0) / nsamples_;
    }
    double equivalent_minority = 0.;
    if (minority_ != NULL)
        equivalent_minority = (double) count_ones_vector(minority_[0].truthtable, nsamples_) / nsamples_;
    root_ = new Node(nrules_, default_prediction, objective, equivalent_minority);
    min_objective_ = objective;
    logger->setTreeMinObj(objective);
    ++num_nodes_;
    logger->setTreeNumNodes(num_nodes_);
    opt_predictions_.push_back(default_prediction);
    logger->setTreePrefixLen(0);
}

/*
 * Insert a node into the tree.
 */
void CacheTree::insert(Node* node, unsigned short thread_id) {
    if (node == NULL) {
        std::cout << "NULL NODE" << std::endl;
    }
    // TODO: clean up
    /* if (node->parent() == root_) {
    } else {
        tree_lks_[thread_id].lock();
    } */
    // root_lk_.lock();
    node->child_lk_.lock();
    node->parent()->children_.insert(std::make_pair(node->id(), node));
    node->child_lk_.unlock();
    // root_lk_.unlock();
    /* if (node->parent() == root_) {
        root_lk_.unlock();
    } else {
        tree_lks_[thread_id].unlock();
    } */
    ++num_nodes_;
    logger->setTreeNumNodes(num_nodes_);
}

/*
 * Removes nodes with no children, recursively traversing tree towards the root.
 */
void CacheTree::prune_up(Node* node) {
    unsigned short id;
    size_t depth = node->depth();
    Node* parent;
    while (node->children_.size() == 0) {
        if (depth > 0) {
            //id = node->id();
            //parent = node->parent();
            //parent->children_.erase(id);
            // --num_nodes_;
            node->lock();
            node->set_deleted();
            node->unlock();
            //delete node;
            node = parent;
            --depth;
        } else {
            // --num_nodes_;
            break;
        }
    }
    logger->setTreeNumNodes(num_nodes_);
}

/*
 * Checks that the prefix is in the tree and hasn't been deleted.
 * Returns NULL if the prefix isn't in the tree, a pointer to the prefix node otherwise.
 */
Node* CacheTree::check_prefix(tracking_vector<unsigned short, DataStruct::Tree>& prefix) {
    Node* node = this->root_;
    for(tracking_vector<unsigned short, DataStruct::Tree>::iterator it = prefix.begin(); 
            it != prefix.end(); ++it) {
        node = node->child(*it);
        if (node == NULL)
            return NULL;
    }
    return node;
}

/*
 * Recursive helper function to traverse down the tree, deleting nodes with a lower bound greater 
 * than the minimum objective.
 */
void CacheTree::gc_helper(Node* node, unsigned short thread_id) {
    if (calculate_size_ & (!node->done()))
        logger->addQueueElement(node->depth(), node->lower_bound(), false);
    Node* child;
    double lb;
    std::vector<Node*> children;
    for (typename std::map<unsigned short, Node*>::iterator cit = node->children_.begin(); 
            cit != node->children_.end(); ++cit)
        children.push_back(cit->second);
    for (typename std::vector<Node*>::iterator cit = children.begin(); cit != children.end(); ++cit) {
        child = *cit;
        if (ablation_ != 2)
            lb = child->lower_bound() + c_;
        else
            lb = child->lower_bound();
        if (lb >= min_objective_) {
            //node->delete_child(child->id());
            delete_subtree(this, child, true, false, thread_id);
        } else
            gc_helper(child, thread_id);
    }
}

/*
 * Public wrapper function to garbage collect the subset of the
 * tree containing the rules in the given range in the index
 * array.
 */
void CacheTree::garbage_collect(unsigned short thread_id) {
    if (calculate_size_)
        logger->clearRemainingSpaceSize();

    tracking_vector<unsigned short, DataStruct::Tree>* subrange = get_subrange(thread_id);
    for (typename tracking_vector<unsigned short, DataStruct::Tree>::iterator rit = subrange->begin();
        rit != subrange->end(); rit++) {
            Node* cur_child = root_->child(*rit);
            if (cur_child != NULL)
                gc_helper(root_->child(*rit), thread_id);
    }
}

/*
 * Deletes a subtree of tree by recursively calling itself on node's children.
 * node -- the node at the root of the subtree to be deleted.
 * destructive -- booelan flag indicating whether to delete node or just lazily mark it.
 * update_remaining_state_space -- boolean flag indicating whether to update the size of
 * the remaining search space (optional calculation in logger state)
 */
void delete_subtree(CacheTree* tree, Node* node, bool destructive, 
        bool update_remaining_state_space, unsigned short thread_id) {
    Node* child;
    Node* parent = node->parent();
    parent->delete_child(node->id());
    // Interior (non-leaf) node
    if (node->num_children() != 0) {
        // copy the children map so we iterate over something that is not being modified by future delete_subtree calls
        std::map<unsigned short, Node*> children_copy;
        children_copy.insert(node->children_begin(), node->children_end());
        for(std::map<unsigned short, Node*>::iterator iter = children_copy.begin(); 
                iter != children_copy.end(); ++iter) {
            child = iter->second;
            delete_subtree(tree, child, destructive, update_remaining_state_space, thread_id);
        }
        //std::cout << "DELETE INTERIOR" << node->id() << std::endl;
        // delete interior nodes
        if (node->in_queue()) {
            // Must have deleted/cleared child map here
            node->lock();
            node->set_deleted();
            node->unlock();
        } else {
            logger->removeFromMemory(sizeof(*node), DataStruct::Tree);
            tree->lock(thread_id);
            node->clear_children();
            delete node;
            tree->unlock(thread_id);
            tree->decrement_num_nodes(); 
        }
        //node->set_deleted();
    } else {
        // only delete leaf nodes in destructive mode
        if (node->in_queue()) {
            // Must have deleted/cleared child map here
            node->lock();
            node->set_deleted();
            node->unlock();
        } else if (destructive) {  
            //std::cout << "DELETE LEAF " << node->id() << std::endl;
            logger->removeFromMemory(sizeof(*node), DataStruct::Tree);
            tree->decrement_num_nodes();
            tree->lock(thread_id);
            node->clear_children();
            delete node;
            tree->unlock(thread_id);
        } else {
            //std::cout << "DELETE ELSE " << node->id() << std::endl;
            logger->decPrefixLen(node->depth());
            if (update_remaining_state_space)
                logger->removeQueueElement(node->depth(), node->lower_bound(), false);
            node->lock();
            node->set_deleted();
            node->unlock();
        }
    }
}

size_t nn_helper (Node* node) {
    /*if (node == NULL)
        return 0;

    size_t ret = 1;
    for (std::map<unsigned short, Node*>::iterator it = node->children_begin(); it != node->children_end(); it++) {
        ret += nn_helper(it->second);
    }
    return ret;*/
    return 0;
}
