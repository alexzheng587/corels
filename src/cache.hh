#include "string.h"
#include <iterator>
#include <map>
#include <vector>
#include <stdlib.h>

extern "C" {
  #include "rule.h"
}

class CacheNode;
class CacheTree;


class CacheNode {
  public:
    typedef CacheNode node_type;

    explicit CacheNode(size_t nrules, bool default_prediction, double objective);
    inline size_t id() const;
    inline size_t depth() const;
    inline bool prediction() const;
    inline bool default_prediction() const;
    inline double lower_bound() const;
    inline double objective() const;
    inline node_type* child(size_t idx) const;
    inline node_type* parent() const;
    inline bool done() const;
    inline void set_done();

  private:
    size_t id_;
    size_t depth_;
    bool prediction_;
    bool default_prediction_;
    double lower_bound_;
    double objective_;
    bool done_;
    node_type* parent_;
    std::map<size_t, node_type*> children_;

    CacheNode(size_t id, size_t nrules, bool prediction, bool default_prediction,
              double lower_bound, double objective, CacheNode* parent);

    friend class CacheTree;
};

class CacheTree {
  public:
    typedef CacheNode node_type;

    CacheTree(size_t nsamples, size_t nrules, double c, rule_t *rules, rule_t *labels);
    ~CacheTree();

    inline double min_objective() const;
    inline size_t num_nodes() const;
    inline size_t num_interior() const;
    inline size_t num_evaluated() const;

    void insert_root();
    void insert(size_t new_rule, bool prediction, bool default_prediction,
                double lower_bound, double objective, CacheNode* parent);
    void evaluate_children(CacheNode* parent, VECTOR parent_not_captured);
    node_type* stochastic_select(VECTOR not_captured);
    void prune_up(CacheNode* node);
    void toy(size_t max_num_nodes);
    void play_with_rules();

  private:
    node_type* root_;
    size_t nsamples_;
    size_t nrules_;
    double c_;
    double min_objective_;
    size_t num_nodes_;
    size_t num_interior_;
    size_t num_evaluated_;
    std::vector<rule_t> rules_;
    std::vector<rule_t> labels_;

    void delete_subtree(CacheNode* node);
};

inline size_t CacheNode::id() const {
    return id_;
}

inline size_t CacheNode::depth() const {
    return depth_;
}

inline bool CacheNode::prediction() const {
    return prediction_;
}

inline bool CacheNode::default_prediction() const {
    return default_prediction_;
}

inline double CacheNode::lower_bound() const {
    return lower_bound_;
}

inline double CacheNode::objective() const {
    return objective_;
}

inline CacheNode* CacheNode::child(size_t idx) const {
    return children_.find(idx)->second;
}

inline CacheNode* CacheNode::parent() const {
    return parent_;
}

inline bool CacheNode::done() const{
    return done_;
}

inline void CacheNode::set_done() {
    done_ = 1;
}

inline double CacheTree::min_objective() const {
    return min_objective_;
}

inline size_t CacheTree::num_nodes() const {
    return num_nodes_;
}

inline size_t CacheTree::num_interior() const {
    return num_interior_;
}

inline size_t CacheTree::num_evaluated() const {
    return num_evaluated_;
}
