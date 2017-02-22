#include "pmap.hh"

void prefix_map_garbage_collect(PrefixPermutationMap* p, size_t queue_min_length) {
    typename PrefixPermutationMap::iterator iter;
    size_t num_deleted = 0;
    printf("pmap gc for length %zu: %zu -> ", queue_min_length, p->size());
    for (iter = p->begin(); iter != p->end(); ) {
        if (iter->first.key[0] <= queue_min_length) {
            iter = p->erase(iter);
            ++num_deleted;
        } else {
            ++iter;
        }
    }
    printf("%zu\n", p->size());
    logger.decreasePmapSize(num_deleted);
}

void bfs_prefix_map_garbage_collect(PrefixPermutationMap* p, size_t queue_min_length) {
    size_t pmap_size = p->size();
    printf("bfs pmap gc for length %zu: %zu -> ", queue_min_length, pmap_size);
    p->clear();
    printf("%zu\n", p->size());
    logger.decreasePmapSize(pmap_size);
}

void captured_map_garbage_collect(CapturedPermutationMap* p, size_t queue_min_length) {
    size_t pmap_size = p->size();
    printf("captured pmap gc for length %zu: %zu -> ", queue_min_length, pmap_size);
    p->clear();
    printf("%zu\n", p->size());
    logger.decreasePmapSize(pmap_size);
}

template<class N>
N* prefix_permutation_insert(construct_signature<N> construct_policy, unsigned short new_rule,
                             size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                             double objective, N* parent, int num_not_captured, int nsamples, int len_prefix,
                             double c, double minority, CacheTree<N>* tree, VECTOR not_captured,
                             std::vector<unsigned short> parent_prefix,
                             PrefixPermutationMap* p) {
    typename PrefixPermutationMap::iterator iter;
    parent_prefix.push_back(new_rule);

    unsigned char* ordered = (unsigned char*) malloc(sizeof(unsigned char) * (len_prefix + 1));
    ordered[0] = (unsigned char)len_prefix;
    auto cmp = [&](int i, int j) { return parent_prefix[i] < parent_prefix[j]; };
    std::sort(&ordered[1], &ordered[len_prefix], cmp);

    std::sort(parent_prefix.begin(), parent_prefix.end());
    
    unsigned short *pre_key = (unsigned short*) malloc(sizeof(unsigned short) * (len_prefix + 1));
    pre_key[0] = (unsigned short)len_prefix;
    memcpy(&pre_key[1], &parent_prefix[0], len_prefix * sizeof(unsigned short));
    
    struct prefix_key key = { pre_key };
    N* child = NULL;
    iter = p->find(key);
    if (iter != p->end()) {
        double permuted_lower_bound = iter->second.first;
        if (lower_bound < permuted_lower_bound) {
            N* permuted_node;
            unsigned char* indices = iter->second.second;
            std::vector<unsigned short> permuted_prefix(parent_prefix.size());
            for (unsigned char i = 0; i < indices[0]; ++i)
                permuted_prefix[i] = parent_prefix[indices[i + 1]];

            if ((permuted_node = tree->check_prefix(permuted_prefix)) != NULL) {
                N* permuted_parent = permuted_node->parent();
                permuted_parent->delete_child(permuted_node->id());
                delete_subtree<N>(tree, permuted_node, false, true);
                logger.incPmapDiscardNum();
            } else {
                logger.incPmapNullNum();
            }
            child = construct_policy(new_rule, nrules, prediction,
                                     default_prediction, lower_bound, objective,
                                     parent, num_not_captured, nsamples,
                                     len_prefix, c, minority);
            iter->second = std::make_pair(lower_bound, ordered);
        }
    } else {
        child = construct_policy(new_rule, nrules, prediction,
                                 default_prediction, lower_bound, objective,
                                 parent, num_not_captured, nsamples, len_prefix,
                                 c, minority);
        unsigned char* ordered_prefix = &ordered[0];
        p->insert(std::make_pair(key, std::make_pair(lower_bound, ordered_prefix)));
        logger.incPmapSize();
    }
    return child;
};

template<class N>
N* captured_permutation_insert(construct_signature<N> construct_policy, unsigned short new_rule,
                               size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                               double objective, N* parent, int num_not_captured, int nsamples, int len_prefix,
                               double c, double minority, CacheTree<N>* tree, VECTOR not_captured,
                               std::vector<unsigned short> parent_prefix,
                               CapturedPermutationMap* p) {
    typename CapturedPermutationMap::iterator iter;
    parent_prefix.push_back(new_rule);
    N* child = NULL;
    struct captured_key key;
    rule_vinit(nsamples, &key.key);
    rule_copy(key.key, not_captured, nsamples);
#ifndef GMP
    key.len = (short) nsamples;
#endif
    iter = p->find(key);
    if (iter != p->end()) {
        std::vector<unsigned short> permuted_prefix = iter->second.first;
        double permuted_lower_bound = iter->second.second;
        if (lower_bound < permuted_lower_bound) {
            N* permuted_node;
            if ((permuted_node = tree->check_prefix(permuted_prefix)) != NULL) {
                N* permuted_parent = permuted_node->parent();
                permuted_parent->delete_child(permuted_node->id());
                delete_subtree<N>(tree, permuted_node, false, true);
                logger.incPmapDiscardNum();
            } else {
                logger.incPmapNullNum();
            }
            child = construct_policy(new_rule, nrules, prediction, default_prediction,
                                       lower_bound, objective, parent,
                                        num_not_captured, nsamples, len_prefix, c, minority);
            iter->second = std::make_pair(parent_prefix, lower_bound);
        }
    } else {
        child = construct_policy(new_rule, nrules, prediction, default_prediction,
                                    lower_bound, objective, parent,
                                    num_not_captured, nsamples, len_prefix, c, minority);
        p->insert(std::make_pair(key, std::make_pair(parent_prefix, lower_bound)));
        logger.incPmapSize();
    }
    return child;
};

template BaseNode*
prefix_permutation_insert(construct_signature<BaseNode> construct_policy, unsigned short new_rule,
                          size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                          double objective, BaseNode* parent, int num_not_captured, int nsamples, int len_prefix,
                          double c, double minority, CacheTree<BaseNode>* tree, VECTOR captured,
                          std::vector<unsigned short> parent_prefix,
                          PrefixPermutationMap* p);

template CuriousNode*
prefix_permutation_insert(construct_signature<CuriousNode> construct_policy, unsigned short new_rule,
                          size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                          double objective, CuriousNode* parent, int num_not_captured, int nsamples, int len_prefix,
                          double c, double minority, CacheTree<CuriousNode>* tree, VECTOR captured,
                          std::vector<unsigned short> parent_prefix,
                          PrefixPermutationMap* p);

template BaseNode*
captured_permutation_insert(construct_signature<BaseNode> construct_policy, unsigned short new_rule,
                            size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                            double objective, BaseNode* parent, int num_not_captured, int nsamples, int len_prefix,
                            double c, double minority, CacheTree<BaseNode>* tree, VECTOR captured,
                            std::vector<unsigned short> parent_prefix,
                            CapturedPermutationMap* p);

template CuriousNode*
captured_permutation_insert(construct_signature<CuriousNode> construct_policy, unsigned short new_rule,
                            size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                            double objective, CuriousNode* parent, int num_not_captured, int nsamples, int len_prefix,
                            double c, double minority, CacheTree<CuriousNode>* tree, VECTOR captured,
                            std::vector<unsigned short> parent_prefix,
                            CapturedPermutationMap* p);
