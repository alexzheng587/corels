#include "cache.hh"
#include "utils.hh"
#include <memory>
#include <vector>
#include <stdlib.h>

Node::Node(size_t nrules, bool default_prediction, double objective, double minority)
    : lower_bound_(minority), objective_(objective), minority_(minority), depth_(0),
      num_captured_(0), id_(0), default_prediction_(default_prediction),
      done_(0), deleted_(0) {
          (void) nrules;
}

Node::Node(unsigned short id, size_t nrules, bool prediction,
              bool default_prediction, double lower_bound, double objective,
              Node* parent, size_t num_captured, double minority)
    : parent_(parent), lower_bound_(lower_bound), objective_(objective),
      minority_(minority), depth_(1 + parent->depth_),
      num_captured_(num_captured), id_(id), prediction_(prediction),
      default_prediction_(default_prediction), done_(0), deleted_(0) {
          (void) nrules;
}

/*CuriousNode::CuriousNode(unsigned short id, size_t nrules, bool prediction,
              bool default_prediction, double lower_bound, double objective,
              double storage, CuriousNode* parent, size_t num_captured, double minority)
    : parent_(parent), lower_bound_(lower_bound), objective_(objective),
      minority_(minority), depth_(1 + parent->depth()),
      num_captured_(num_captured), id_(id), prediction_(prediction),
      default_prediction_(default_prediction), done_(0), deleted_(0),
      storage_(storage) {
}*/

CacheTree::CacheTree(size_t nsamples, size_t nrules, double c, rule_t *rules,
                        rule_t *labels, rule_t *meta)
    : root_(0), nsamples_(nsamples), nrules_(nrules), c_(c),
      num_nodes_(0), num_evaluated_(0), min_objective_(0.5),
      opt_rulelist_({}), opt_predictions_({}) {
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

    logger.setTreeMinObj(min_objective_);
    logger.setTreeNumNodes(num_nodes_);
    logger.setTreeNumEvaluated(num_evaluated_);
}

CacheTree::~CacheTree() {
    if (root_)
        delete_subtree(this, root_, true, false);
}

/*Node* CacheTree::construct_node(unsigned short new_rule, size_t nrules, bool prediction,
           bool default_prediction, double lower_bound, double objective,
           Node* parent, int num_not_captured, int nsamples,
           int len_prefix, double c, double minority) {
    (void) len_prefix, (void) c;
    size_t num_captured = nsamples - num_not_captured;
    return (new Node(new_rule, nrules, prediction, default_prediction,
                         lower_bound, objective, parent, num_captured, minority));
}*/

CuriousNode* CacheTree::construct_node(unsigned short new_rule, size_t nrules, bool prediction,
                         bool default_prediction, double lower_bound, double objective, 
                         Node* parent, int num_not_captured, int nsamples,
                         int len_prefix, double c, double minority) {
    size_t num_captured = nsamples - num_not_captured;
    double curiosity = (lower_bound - c * len_prefix + c) * nsamples / (double)(num_captured);
    return (new CuriousNode(new_rule, nrules, prediction, default_prediction,
                            lower_bound, objective, curiosity, (CuriousNode*) parent, num_captured, minority));
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
        objective = (float)(d1) / nsamples_;
    } else {
        default_prediction = 1;
        objective = (float)(d0) / nsamples_;
    }
    double minority = 0.;
    if (meta_.size() > 0)
        minority = (float) count_ones_vector(meta_[0].truthtable, nsamples_) / nsamples_;
    root_ = new Node(nrules_, default_prediction, objective, minority);
    min_objective_ = objective;
    logger.setTreeMinObj(objective);
    ++num_nodes_;
    logger.setTreeNumNodes(num_nodes_);
    opt_predictions_.push_back(default_prediction);
    logger.setTreePrefixLen(0);
}

/*
 * Insert a node into the tree.
 */
void CacheTree::insert(Node* node) {
    node->parent()->children_.insert(std::make_pair(node->id(), node));
    ++num_nodes_;
    logger.setTreeNumNodes(num_nodes_);
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

/*
 * Checks that the prefix is in the tree and hasn't been deleted.
 * Returns NULL if the prefix isn't in the tree, a pointer to the prefix node otherwise.
 */
Node* CacheTree::check_prefix(std::vector<unsigned short>& prefix) {
    Node* node = this->root_;
    for(std::vector<unsigned short>::iterator it = prefix.begin(); it != prefix.end(); ++it) {
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
    if (!node->done())
        logger.addQueueElement(node->depth(), node->lower_bound());
    Node* child;
    std::vector<Node*> children;
    for (typename std::map<unsigned short, Node*>::iterator cit = node->children_.begin(); 
            cit != node->children_.end(); ++cit)
        children.push_back(cit->second);
    for (typename std::vector<Node*>::iterator cit = children.begin(); cit != children.end(); ++cit) {
        child = *cit;
        if ((child->lower_bound() + c_) >= min_objective_) {
            node->delete_child(child->id());
            delete_subtree(this, child, false, false);
        } else
            gc_helper(child);
    }
}

/*
 * Public wrapper function to garbage collect the entire tree beginning from the root.
 */
void CacheTree::garbage_collect() {
    logger.clearRemainingSpaceSize();
    gc_helper(root_);
}

/*
 * Deletes a subtree of tree by recursively calling itself on node's children.
 * node -- the node at the root of the subtree to be deleted.
 * destructive -- booelan flag indicating whether to delete node or just lazily mark it.
 * update_remaining_state_space -- affects whether we remove the deleted elements from the queue.
 */
void delete_subtree(CacheTree* tree, Node* node, bool destructive, 
        bool update_remaining_state_space) {
    Node* child;
    typename std::map<unsigned short, Node*>::iterator iter;

    if (node->done()) {
        iter = node->children_begin();
        while (iter != node->children_end()) {
            child = iter->second;
            delete_subtree(tree, child, destructive, update_remaining_state_space);
            ++iter;
        }
        tree->decrement_num_nodes(); // always delete interior (non-leaf) nodes
        delete node;
    } else {
        if (destructive) {  // only delete leaf nodes in destructive mode
            tree->decrement_num_nodes();
            delete node;
        } else {
            logger.decPrefixLen(node->depth());
            if (update_remaining_state_space)
                logger.removeQueueElement(node->depth(), node->lower_bound());
            node->set_deleted();
        }
    }
}

/*
// BaseNode
template class Node<bool>;

// CuriousNode
template class Node<double>; 

template class CacheTree<BaseNode>;

template class CacheTree<CuriousNode>;

template void
delete_subtree<BaseNode>(CacheTree<BaseNode>* tree, BaseNode* n, bool destructive, bool update_remaining_state_space);

template void
delete_subtree<CuriousNode>(CacheTree<CuriousNode>* tree, CuriousNode* n, bool destructive, bool update_remaining_state_space);
*/
