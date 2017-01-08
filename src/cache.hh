#include "string.h"
#include <iterator>
#include <map>
#include <vector>
#include <stdlib.h>

extern "C" {
  #include "rule.h"
}

template <class T> class Node;
template <class N> class CacheTree;

typedef Node<bool> BaseNode;       // nothing extra
typedef Node<double> CuriousNode;  // curiosity

template <class T>
class Node {
  public:
    explicit Node(size_t nrules, bool default_prediction, double objective, double minority);

    Node(unsigned short id, size_t nrules, bool prediction, bool default_prediction,
         double lower_bound, double objective, T storage, Node<T>* parent,
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

    // Returns pair of prefixes and predictions for the path from this
    // node to the root
    inline std::pair<std::vector<unsigned short>, std::vector<bool>>
        get_prefix_and_predictions();

    inline size_t depth() const;
    inline Node<T>* child(unsigned short idx);
    inline Node<T>* parent() const;
    inline void delete_child(unsigned short idx);
    inline size_t num_children() const;

    inline T& get_storage(); // can this be const?
    inline size_t num_captured() const;
    inline double minority() const;

    inline typename std::map<unsigned short, Node<T>*>::iterator children_begin();
    inline typename std::map<unsigned short, Node<T>*>::iterator children_end();
    inline typename std::map<unsigned short, Node<T>*>::iterator random_child(); // FIXME
    // inline typename std::map<unsigned short, Node<T>*>::iterator random_child(PRNG prng);

  private:

    unsigned short id_;
    bool prediction_;
    bool default_prediction_;
    double lower_bound_;
    double objective_;
    bool done_;
    bool deleted_;

    size_t depth_;
    Node<T>* parent_;
    std::map<unsigned short, Node<T>*> children_;

    T storage_;  // space for something extra, like curiosity or a bit vector
    size_t num_captured_;
    double minority_;

    friend class CacheTree<Node<T> >;
};

template<class N>
class CacheTree {
  public:
    CacheTree(size_t nsamples, size_t nrules, double c, rule_t *rules,
              rule_t *labels, rule_t *meta);
    ~CacheTree();

    inline double min_objective() const;
    inline std::vector<unsigned short> opt_rulelist() const;
    inline std::vector<bool> opt_predictions() const;

    inline size_t num_nodes() const;
    inline size_t num_evaluated() const;
    inline rule_t rule(unsigned short idx) const;
    inline char* rule_features(unsigned short idx) const;
    inline rule_t label(unsigned short idx) const;
    inline rule_t meta(unsigned short idx) const;
    inline size_t nsamples() const;
    inline size_t nrules() const;
    inline double c() const;
    inline N* root() const;

    void update_min_objective(double objective);
    void update_opt_rulelist(std::vector<unsigned short>& parent_prefix,
                             unsigned short new_rule_id);
    void update_opt_predictions(std::vector<bool>& parent_predictions,
                                bool new_pred,
                                bool new_default_pred);

    void increment_num_evaluated();
    void decrement_num_nodes();

    void insert_root();
    void insert(N* node);
    void prune_up(N* node);
    void garbage_collect();
    void play_with_rules();
    N* check_prefix(std::vector<unsigned short>& prefix);

  private:
    N* root_;
    size_t nsamples_;
    size_t nrules_;
    double c_;

    size_t num_nodes_;
    size_t num_evaluated_;

    double min_objective_;
    std::vector<unsigned short> opt_rulelist_;
    std::vector<bool> opt_predictions_;

    std::vector<rule_t> rules_;
    std::vector<rule_t> labels_;
    std::vector<rule_t> meta_;

    void gc_helper(N* node);
};

template <class T>
inline unsigned short Node<T>::id() const {
    return id_;
}

template <class T>
inline bool Node<T>::prediction() const {
    return prediction_;
}

template <class T>
inline bool Node<T>::default_prediction() const {
    return default_prediction_;
}

template <class T>
inline double Node<T>::lower_bound() const {
    return lower_bound_;
}

template <class T>
inline double Node<T>::objective() const {
    return objective_;
}

template <class T>
inline bool Node<T>::done() const{
    return done_;
}

template <class T>
inline void Node<T>::set_done() {
    done_ = 1;
}

template <class T>
inline bool Node<T>::deleted() const{
    return deleted_;
}

template <class T>
inline void Node<T>::set_deleted() {
    deleted_ = 1;
}

template <class T>
inline std::pair<std::vector<unsigned short>, std::vector<bool>>
    Node<T>::get_prefix_and_predictions() {
    std::vector<unsigned short> prefix;
    std::vector<bool> predictions;
    auto it1 = prefix.begin();
    auto it2 = predictions.begin();
    Node<T>* node = this;
    for(size_t i = depth_; i > 0; --i) {
        it1 = prefix.insert(it1, node->id());
        it2 = predictions.insert(it2, node->prediction());
        node = node->parent();
    }
    return std::make_pair(prefix, predictions);
}

template <class T>
inline size_t Node<T>::depth() const {
    return depth_;
}

template<class T>
inline Node<T>* Node<T>::child(unsigned short idx) {
    typename std::map<unsigned short, Node<T>*>::iterator iter;
    iter = children_.find(idx);
    if (iter == children_.end())
        return NULL;
    else
        return iter->second;
}

template<class T>
inline void Node<T>::delete_child(unsigned short idx) {
    children_.erase(idx);
}

template<class T>
inline size_t Node<T>::num_children() const {
    return children_.size();
}

template<class T>
inline typename std::map<unsigned short, Node<T>*>::iterator Node<T>::children_begin() {
    return children_.begin();
}

template<class T>
inline typename std::map<unsigned short, Node<T>*>::iterator Node<T>::children_end() {
    return children_.end();
}

template<class T>
inline typename std::map<unsigned short, Node<T>*>::iterator Node<T>::random_child() {
    typename std::map<unsigned short, Node<T>*>::iterator iter;
    unsigned short idx;
    iter = children_.begin();
    idx = rand() % (children_.size());
    std::advance(iter, idx);
    return iter;
}

template<class T>
inline Node<T>* Node<T>::parent() const {
    return parent_;
}

template<class T>
inline T& Node<T>::get_storage() {
    return storage_;
}

template<class T>
inline size_t Node<T>::num_captured() const {
    return num_captured_;
}

template<class T>
inline double Node<T>::minority() const {
    return minority_;
}

template<class N>
inline double CacheTree<N>::min_objective() const {
    return min_objective_;
}

template<class N>
inline std::vector<unsigned short> CacheTree<N>::opt_rulelist() const {
    return opt_rulelist_;
}

template<class N>
inline std::vector<bool> CacheTree<N>::opt_predictions() const {
    return opt_predictions_;
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
inline rule_t CacheTree<N>::rule(unsigned short idx) const{
    return rules_[idx];
}

template<class N>
inline char* CacheTree<N>::rule_features(unsigned short idx) const{
    return rules_[idx].features;
}

template<class N>
inline rule_t CacheTree<N>::label(unsigned short idx) const{
    return labels_[idx];
}

template<class N>
inline rule_t CacheTree<N>::meta(unsigned short idx) const{
    return meta_[idx];
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
