#include "string.h"
#include <iterator>
#include <map>
#include <vector>
#include <stdlib.h>

extern "C" {
  #include "rule.h"
}

class CacheNode;
class CacheNodeCurious;
template<class N> class CacheTree;


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
    inline void delete_child(size_t idx);
    inline size_t num_children() const;
    inline std::map<size_t, CacheNode*>::iterator random_child();

  protected:
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

    friend class CacheTree<CacheNode>;
};

class CacheNodeCurious: public CacheNode {
    public:
        explicit CacheNodeCurious(size_t nrules, bool default_prediction, double objective);

        inline double curiosity() const;

    protected:
        double curiosity_;

        CacheNodeCurious(size_t id, size_t nrules, bool prediction, bool default_prediction,
                  double lower_bound, double objective, CacheNodeCurious* parent, double curiosity);
};

template<class N>
class CacheTree {
  public:
    CacheTree(size_t nsamples, size_t nrules, double c, rule_t *rules, rule_t *labels);
    ~CacheTree();

    inline double min_objective() const;
    inline size_t num_nodes() const;
    inline size_t num_evaluated() const;
    inline rule_t rule(size_t idx) const;
    inline rule_t label(size_t idx) const;
    inline size_t nsamples() const;
    inline size_t nrules() const;
    inline double c() const;
    inline N* root() const;

    inline void update_min_objective(double objective);
    inline void increment_num_evaluated();

    void insert_root();
    void insert(size_t new_rule, bool prediction, bool default_prediction,
                double lower_bound, double objective, N* parent);
    void prune_up(N* node);
    void delete_subtree(N* node);
    void play_with_rules();

  private:
    N* root_;
    size_t nsamples_;
    size_t nrules_;
    double c_;
    double min_objective_;
    size_t num_nodes_;
    size_t num_evaluated_;
    std::vector<rule_t> rules_;
    std::vector<rule_t> labels_;
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

inline void CacheNode::delete_child(size_t idx) {
    children_.erase(idx);
}

inline size_t CacheNode::num_children() const {
    return children_.size();
}

inline std::map<size_t, CacheNode*>::iterator CacheNode::random_child() {
    std::map<size_t, CacheNode*>::iterator iter;
    size_t idx;
    iter = children_.begin();
    idx = rand() % (children_.size());
    std::advance(iter, idx);
    return iter;
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

inline double CacheNodeCurious::curiosity() const {
    curiosity_;
}

template<class N>
inline double CacheTree<N>::min_objective() const {
    return min_objective_;
}

template<class N>
inline size_t CacheTree<N>::num_nodes() const {
    return num_nodes_;
}

template<class N>
inline size_t CacheTree<N>::num_evaluated() const {
    return num_evaluated_;
}

template<class N>
inline rule_t CacheTree<N>::rule(size_t idx) const{
    return rules_[idx];
}

template<class N>
inline rule_t CacheTree<N>::label(size_t idx) const{
    return labels_[idx];
}

template<class N>
inline size_t CacheTree<N>::nsamples() const {
    return nsamples_;
}

template<class N>
inline size_t CacheTree<N>::nrules() const {
    return nrules_;
}

template<class N>
inline double CacheTree<N>::c() const {
    return c_;
}

template<class N>
inline N* CacheTree<N>::root() const {
    return root_;
}

template<class N>
inline void CacheTree<N>::update_min_objective(double objective) {
    min_objective_ = objective;
}

template<class N>
inline void CacheTree<N>::increment_num_evaluated() {
    ++num_evaluated_;
}
