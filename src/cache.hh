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
    explicit Node(size_t nrules, bool default_prediction, double objective);

    Node(size_t id, size_t nrules, bool prediction, bool default_prediction,
         double lower_bound, double objective, T storage, Node<T>* parent);

    inline size_t id() const;
    inline bool prediction() const;
    inline bool default_prediction() const;
    inline double lower_bound() const;
    inline double objective() const;
    inline bool done() const;
    inline void set_done();

    inline size_t depth() const;
    inline Node<T>* child(size_t idx) const;
    inline Node<T>* parent() const;
    inline void delete_child(size_t idx);
    inline size_t num_children() const;

    inline T& get_storage(); // can this be const?

    inline typename std::map<size_t, Node<T>*>::iterator random_child(); // FIXME
    // inline typename std::map<size_t, Node<T>*>::iterator random_child(PRNG prng);

  private:

    size_t id_;
    bool prediction_;
    bool default_prediction_;
    double lower_bound_;
    double objective_;
    bool done_;

    size_t depth_;
    Node<T>* parent_;
    std::map<size_t, Node<T>*> children_;

    T storage_;  // space for something extra, like curiosity or a bit vector

    friend class CacheTree<Node<T> >;
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
    void insert(N* node);
    void prune_up(N* node);
    template <class P>
    void delete_subtree(N* node, P* p = NULL);
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

template <class T>
inline size_t Node<T>::id() const {
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
inline size_t Node<T>::depth() const {
    return depth_;
}

template<class T>
inline Node<T>* Node<T>::child(size_t idx) const {
    return children_.find(idx)->second;
}

template<class T>
inline void Node<T>::delete_child(size_t idx) {
    children_.erase(idx);
}

template<class T>
inline size_t Node<T>::num_children() const {
    return children_.size();
}

template<class T>
inline typename std::map<size_t, Node<T>*>::iterator Node<T>::random_child() {
    typename std::map<size_t, Node<T>*>::iterator iter;
    size_t idx;
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
