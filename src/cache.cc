# ifdef VAL
#include "common.hh"
# endif 
#include "cache.hh"
#include "loss.hh"
#include "utils.hh"
#include <memory>
#include <vector>
#include <stdlib.h>

// each node has lower bound and objective for incremental computation
Node::Node(bool default_prediction, double objective, double equivalent_minority)
        : parent_(0), lower_bound_(equivalent_minority), objective_(objective), equivalent_minority_(equivalent_minority),
          depth_(0), num_captured_(0), id_(0), default_prediction_(default_prediction), done_(0), deleted_(0) {
}

Node::Node(unsigned short id, bool prediction, bool default_prediction, double lower_bound, double objective,
           Node* parent, size_t num_captured, double equivalent_minority)
        : parent_(parent), lower_bound_(lower_bound), objective_(objective),
          equivalent_minority_(equivalent_minority), depth_(1 + parent->depth_),
          num_captured_(num_captured), id_(id), prediction_(prediction),
          default_prediction_(default_prediction), done_(0), deleted_(0) {
}

CacheTree::CacheTree(Loss* loss, size_t nsamples, size_t nrules, double c, size_t k, rule_t* rules, rule_t* labels,
                     rule_t* minority, minority_class_t* minor_class, size_t nclasses, int ablation, bool calculate_size, char const* node_type)
        : loss_(loss), root_(0), nsamples_(nsamples), nrules_(nrules), c_(c), k_(k),
          num_nodes_(0), num_evaluated_(0), nclasses_(nclasses), ablation_(ablation), calculate_size_(calculate_size), node_type_(node_type) {
    optimal_rules.push_back({0.5, {}, {}});
    rules_ = rules;
    labels_ = labels;
    nneg_ = labels_[0].support;
    npos_ = nsamples_ - nneg_;
    loss->initialize_loss(this);
    if (minority) {
        minority_ = minority;
    } else {
        minority_ = NULL;
    }

    if (minor_class) {
        minor_class_ = minor_class;
    } else {
        minor_class_ = NULL;
    }

    // TODO update when multithreaded
    thread_states_.resize(1);
    for (int i = 0; i < 1; i++) {
        thread_states_[i] = unique_ptr<ProblemState>(new ProblemState(nsamples));
    }

    logger->setTreeMinObj(optimal_rules[0].min_objective_);
    logger->setTreeNumNodes(num_nodes_);
    logger->setTreeNumEvaluated(num_evaluated_);
}

CacheTree::~CacheTree() {
    if (root_)
        delete_subtree(this, root_, true, false);

}

// Node* CacheTree::construct_node(unsigned short new_rule, size_t nrules, bool prediction, bool default_prediction,
//                                 double lower_bound, double objective, Node* parent, int num_not_captured, int nsamples,
//                                 double equivalent_minority) {
//     size_t num_captured = nsamples - num_not_captured;
//     Node* n;
//     // curiosity calculation not relevant for the different loss functions
//     if (strcmp(node_type_, "curious") == 0) {
//         double curiosity = (lower_bound - equivalent_minority) * nsamples / (double) (num_captured);
//         n = (Node * )(new CuriousNode(new_rule, nrules, prediction, default_prediction,
//                                       lower_bound, objective, curiosity, (CuriousNode*) parent, num_captured,
//                                       equivalent_minority));
//     } else {
//         n = new Node(new_rule, nrules, prediction, default_prediction, lower_bound, objective, parent, num_captured,
//                      equivalent_minority);
//     }
//     logger->addToMemory(sizeof(*n), DataStruct::Tree);
//     return n;
// }

// Node* CacheTree::construct_auc_node(unsigned short new_rule, size_t nrules, bool prediction, bool default_prediction,
//                                 double lower_bound, double objective, Node* parent, int num_not_captured, int nsamples,
//                                 double equivalent_minority, minority_class_t* minority_class, int npos_captured, int nneg_captured) {
//         size_t num_captured = nsamples - num_not_captured;
//         Node* n = (Node *)(new AUCNode(new_rule, nrules, prediction, default_prediction,
//                                       lower_bound, objective, (AUCNode*) parent, num_captured,
//                                       equivalent_minority, minority_class, npos_captured, nneg_captured));
//         return n;
// }

/*
 * Inserts the root of the tree, setting up the default rules.
 */
void CacheTree::insert_root() {
    VECTOR tmp_vec;
    VECTOR equivalent_fn, equivalent_fp;
    int num_equivalent_fn = 0, num_equivalent_fp = 0;
    double ld0, ld1;
    bool default_prediction;
    double objective;
    make_default(&tmp_vec, nsamples_);
    // loss calculation for default rule in root node
    ld0 = loss_->compute_loss(0, npos_, npos_, nneg_);
    ld1 = loss_->compute_loss(nneg_, 0, npos_, nneg_);
    // choose prediction that minimizes loss
    if (ld0 < ld1) {
        default_prediction = false;
        objective = ld0;
    } else {
        default_prediction = true;
        objective = ld1;
    }
    double equivalent_minority = 0.;

    // equivalent points bound
    if (minority_ != nullptr && loss_->type() != AUC) {
        rule_vinit(nsamples_, &equivalent_fn);
        rule_vinit(nsamples_, &equivalent_fp);

        rule_vand(equivalent_fn, minority_[0].truthtable, label(1).truthtable, nsamples_, &num_equivalent_fn);
        rule_vand(equivalent_fp, minority_[0].truthtable, label(0).truthtable, nsamples_, &num_equivalent_fp);
        equivalent_minority = loss_->compute_loss(num_equivalent_fp, num_equivalent_fn, npos_, nneg_);

        rule_vfree(&equivalent_fn);
        rule_vfree(&equivalent_fp);
    }

    switch (loss_->type()) {
        case F_SCORE:
            root_ = new F1Node(default_prediction, objective, equivalent_minority);
            break;
        case AUC:
            root_ = new AUCNode(default_prediction, objective, equivalent_minority, nclasses_, minor_class_, npos_, nneg_);
            break;
        case ACCURACY: // Fall through
        case BALANCED_ACCURACY: // Fall through
        case WEIGHTED_ACCURACY: // Fall through
        default:
            root_ = new Node(default_prediction, objective, equivalent_minority);
            break;
    }

    optimal_rules[0].min_objective_ = objective;
    logger->setTreeMinObj(objective);
    ++num_nodes_;
    logger->setTreeNumNodes(num_nodes_);
    optimal_rules[0].opt_predictions_.push_back(default_prediction);
    logger->setTreePrefixLen(0);
}

/*
 * Insert a node into the tree.
 */
void CacheTree::insert(Node* node) {
    node->parent()->children_.insert(std::make_pair(node->id(), node));
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
    logger->setTreeNumNodes(num_nodes_);
}

/*
 * Checks that the prefix is in the tree and hasn't been deleted.
 * Returns NULL if the prefix isn't in the tree, a pointer to the prefix node otherwise.
 */
Node* CacheTree::check_prefix(tracking_vector<unsigned short, DataStruct::Tree>& prefix) {
    Node* node = this->root_;
    for (tracking_vector<unsigned short, DataStruct::Tree>::iterator it = prefix.begin();
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
void CacheTree::gc_helper(Node* node) {
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
        if (lb >= optimal_rules[0].min_objective_) {
            // destructively call delete_subtree
            delete_subtree(this, child, true, false);
        } else {
            gc_helper(child);
        }
    }
}

/*
 * Public wrapper function to garbage collect the entire tree beginning from the root.
 */
void CacheTree::garbage_collect() {
    if (calculate_size_)
        logger->clearRemainingSpaceSize();
    gc_helper(root_);
}

/*
 * Deletes a subtree of tree by recursively calling itself on node's children.
 * node -- the node at the root of the subtree to be deleted.
 * destructive -- booelan flag indicating whether to delete node or just lazily mark it.
 * update_remaining_state_space -- boolean flag indicating whether to update the size of
 * the remaining search space (optional calculation in logger state)
 */
void delete_subtree(CacheTree* tree, Node* node, bool destructive,
                    bool update_remaining_state_space) {
    Node* child;
    Node* parent = node->parent();
    if (parent != NULL)
        parent->delete_child(node->id());
    // Interior (non-leaf) node
    if (node->num_children() != 0) {
        // copy the children map so we iterate over something that is not being modified by future delete_subtree calls
        std::map<unsigned short, Node*> children_copy;
        children_copy.insert(node->children_begin(), node->children_end());
        for (std::map<unsigned short, Node*>::iterator iter = children_copy.begin();
             iter != children_copy.end(); ++iter) {
            child = iter->second;
            delete_subtree(tree, child, destructive, update_remaining_state_space);
        }
        //std::cout << "DELETE INTERIOR" << node->id() << std::endl;
        // delete interior nodes
        if (node->in_queue()) {
            // Must have deleted/cleared child map here
            node->set_deleted();
        } else {
            logger->removeFromMemory(sizeof(*node), DataStruct::Tree);
            node->clear_children();
            delete node;
            tree->decrement_num_nodes();
        }
        //node->set_deleted();
    } else {
        // only delete leaf nodes in destructive mode
        if (node->in_queue()) {
            // Must have deleted/cleared child map here
            node->set_deleted();
        } else if (destructive) {
            //std::cout << "DELETE LEAF " << node->id() << std::endl;
            logger->removeFromMemory(sizeof(*node), DataStruct::Tree);
            tree->decrement_num_nodes();
            node->clear_children();
            delete node;
        } else {
            //std::cout << "DELETE ELSE " << node->id() << std::endl;
            logger->decPrefixLen(node->depth());
            if (update_remaining_state_space)
                logger->removeQueueElement(node->depth(), node->lower_bound(), false);
            node->set_deleted();
        }
    }
}
