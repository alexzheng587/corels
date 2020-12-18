#include "pmap.hh"

PrefixPermutationMap::PrefixPermutationMap()
        : pmap(new PrefixMap) {}

CapturedPermutationMap::CapturedPermutationMap()
        : pmap(new CapturedMap) {}

AUCPermutationMap::AUCPermutationMap()
        : pmap(new PrefixMap) {}

Node* PrefixPermutationMap::insert(unsigned short new_rule, ProblemState* state, Node* parent, CacheTree* tree,
                                   tracking_vector<unsigned short, DataStruct::Tree> parent_prefix) {
    int len_prefix = parent->depth() + 1;
    double lower_bound = state->rule_lower_bound;

    logger->incPermMapInsertionNum();
    parent_prefix.push_back(new_rule);

    auto* ordered = (unsigned char*) malloc(sizeof(unsigned char) * (len_prefix + 1));
    ordered[0] = (unsigned char) len_prefix;

    for (int i = 1; i < (len_prefix + 1); i++)
        ordered[i] = i - 1;

    std::function<bool(int, int)> cmp = [&](int i, int j) { return parent_prefix[i] < parent_prefix[j]; };
    std::sort(&ordered[1], &ordered[len_prefix + 1], cmp);

    std::sort(parent_prefix.begin(), parent_prefix.end());
    auto* pre_key = (unsigned short*) malloc(sizeof(unsigned short) * (len_prefix + 1));
    pre_key[0] = (unsigned short) len_prefix;
    memcpy(&pre_key[1], &parent_prefix[0], len_prefix * sizeof(unsigned short));

    logger->addToMemory((len_prefix + 1) * (sizeof(unsigned char) + sizeof(unsigned short)), DataStruct::Pmap);
    prefix_key key = {pre_key};

    Node* child = nullptr;
    auto iter = pmap->find(key);
    if (iter != pmap->end()) {
        double permuted_lower_bound = iter->second.first;
        if (lower_bound < permuted_lower_bound) {
            Node* permuted_node;
            tracking_vector<unsigned short, DataStruct::Tree> permuted_prefix(parent_prefix.size());
            unsigned char* indices = iter->second.second;
            for (unsigned char i = 0; i < indices[0]; ++i)
                permuted_prefix[i] = parent_prefix[indices[i + 1]];

            if ((permuted_node = tree->check_prefix(permuted_prefix)) != nullptr) {
                Node* permuted_parent = permuted_node->parent();
                permuted_parent->delete_child(permuted_node->id());
                permuted_node->set_deleted();
                logger->incPmapDiscardNum();
            } else {
                logger->incPmapNullNum();
            }
            child = parent->construct_child(new_rule, state);
            iter->second = std::make_pair(lower_bound, ordered);
        }
    } else {
        child = parent->construct_child(new_rule, state);
        unsigned char* ordered_prefix = &ordered[0];
        pmap->insert(std::make_pair(key, std::make_pair(lower_bound, ordered_prefix)));
        logger->incPmapSize();
    }
    return child;
}

Node* CapturedPermutationMap::insert(unsigned short new_rule, ProblemState* state, Node* parent, CacheTree* tree,
                                     tracking_vector<unsigned short, DataStruct::Tree> parent_prefix) {
    double lower_bound = state->rule_lower_bound;
    int nsamples = tree->nsamples();

    logger->incPermMapInsertionNum();
    parent_prefix.push_back(new_rule);
    Node* child = nullptr;
    captured_key key{};
    rule_vinit(nsamples, &key.key);
    rule_copy(key.key, state->not_captured, nsamples);
#ifndef GMP
    key.len = (short) nsamples;
#endif
    auto iter = pmap->find(key);
    if (iter != pmap->end()) {
        double permuted_lower_bound = iter->second.first;
        tracking_vector<unsigned short, DataStruct::Tree> permuted_prefix = iter->second.second;
        if (lower_bound < permuted_lower_bound) {
            Node* permuted_node;
            if ((permuted_node = tree->check_prefix(permuted_prefix)) != nullptr) {
                Node* permuted_parent = permuted_node->parent();
                permuted_parent->delete_child(permuted_node->id());
                delete_subtree(tree, permuted_node, false, tree->calculate_size());
                logger->incPmapDiscardNum();
            } else {
                logger->incPmapNullNum();
            }
            child = parent->construct_child(new_rule, state);
            iter->second = std::make_pair(lower_bound, parent_prefix);
        }
    } else {
        child = parent->construct_child(new_rule, state);
        pmap->insert(std::make_pair(key, std::make_pair(lower_bound, parent_prefix)));
        logger->incPmapSize();
    }
    return child;
}

Node* FScorePermutationMap::insert(unsigned short new_rule, ProblemState *state, Node *parent, CacheTree *tree,
                                   tracking_vector<unsigned short, DataStruct::Tree> parent_prefix) {
    Node* child = parent->construct_child(new_rule, state);
    logger->addToMemory(sizeof(*child), DataStruct::Tree);
    return child;
}

Node* AUCPermutationMap::insert(unsigned short new_rule, ProblemState *state, Node *parent, CacheTree *tree,
                                   tracking_vector<unsigned short, DataStruct::Tree> parent_prefix) {
    Node* child = parent->construct_child(new_rule, state);
    logger->addToMemory(sizeof(*child), DataStruct::Tree);
    return child;
}