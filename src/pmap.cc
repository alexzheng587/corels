#ifdef VAL
#include "common.hh"
#endif
#include "pmap.hh"
#include <assert.h>

PrefixPermutationMap::PrefixPermutationMap()
    : pmap(new PrefixMap) {
}

PrefixPermutationMap::~PrefixPermutationMap() {
    for (PrefixMap::iterator iter = pmap->begin(); iter != pmap->end(); ++iter) {
        //prefix_key pkey = iter->first;
        //free(pkey.key);
        unsigned char *old_ordered = std::get<1>(iter->second);
        free(old_ordered);
    }
    delete pmap;
}

CapturedPermutationMap::CapturedPermutationMap()
    : pmap(new CapturedMap) {}

CapturedPermutationMap::~CapturedPermutationMap() {
    delete pmap;
}

Node *PrefixPermutationMap::insert(unsigned short new_rule, size_t nrules, bool prediction,
                                   bool default_prediction, double lower_bound, double objective, Node *parent,
                                   int num_not_captured, int nsamples, int len_prefix, double c, double equivalent_minority,
                                   CacheTree *tree, VECTOR not_captured, tracking_vector<unsigned short, DataStruct::Tree> prefix,
                                   unsigned short thread_id) {
    prefix.push_back(new_rule);
    assert(prefix.size() == len_prefix);
    (void)not_captured;
    logger->incPermMapInsertionNum();

    // Initialization of prefix_key permutation
    std::pair<prefix_key, unsigned char *> prefix_key_and_order = construct_prefix_key_from_prefix_vector(prefix, new_rule, len_prefix);
    prefix_key key = prefix_key_and_order.first;
    unsigned char *ordered = prefix_key_and_order.second;
    Node *child = NULL;

    // Checking membership of prefix_key permutation
    PrefixMap::iterator iter = find_prefix_key_in_pmap(key);
    if (iter != pmap->end()) {
        double permuted_lower_bound = std::get<0>(iter->second);
        if (lower_bound < permuted_lower_bound) {
            remove_existing_node(tree, iter, prefix, thread_id);
            // Create new node and insert
            child = tree->construct_node(new_rule, nrules, prediction,
                                         default_prediction, lower_bound, objective,
                                         parent, num_not_captured, nsamples,
                                         len_prefix, c, equivalent_minority);
            unsigned char *old_ordered = std::get<1>(iter->second);
            free(old_ordered);
            iter->second = std::make_tuple(lower_bound, ordered, thread_id);
        }
    } else {
        // Create new node if permutation doesn't exist
        child = tree->construct_node(new_rule, nrules, prediction,
                                     default_prediction, lower_bound, objective,
                                     parent, num_not_captured, nsamples, len_prefix,
                                     c, equivalent_minority);
        // Need to globally lock map when inserting otherwise iterators in other threads could be
        // invalidated if the unordered_map is resized.
        map_lk.lock();
        pmap->insert(std::make_pair(key, std::make_tuple(lower_bound, ordered, thread_id)));
        map_lk.unlock();
        logger->incPmapSize();
    }

    // Clean up/Wake up any other threads waiting on this entry
    {
        std::lock_guard<std::mutex> key_lk(key_lk_);
        active_keys.erase(key);
    }
    key_cv.notify_all();
    return child;
}

std::pair<prefix_key, unsigned char *> PrefixPermutationMap::construct_prefix_key_from_prefix_vector(tracking_vector<unsigned short, DataStruct::Tree> prefix,
                                                                                                     unsigned short new_rule, int len_prefix) {
    // Initializing permutation constructs for purposes of comparison
    unsigned char *ordered = (unsigned char *)malloc(sizeof(unsigned char) * (len_prefix + 1));
    ordered[0] = (unsigned char)len_prefix;

    for (int i = 1; i < (len_prefix + 1); i++)
        ordered[i] = i - 1;

    std::function<bool(int, int)> cmp = [&](int i, int j) { return prefix[i] < prefix[j]; };
    std::sort(&ordered[1], &ordered[len_prefix + 1], cmp);

    std::sort(prefix.begin(), prefix.end());
    unsigned short *pre_key = (unsigned short *)malloc(sizeof(unsigned short) * (len_prefix + 1));
    pre_key[0] = (unsigned short)len_prefix;
    memcpy(&pre_key[1], &prefix[0], len_prefix * sizeof(unsigned short));

    logger->addToMemory((len_prefix + 1) * (sizeof(unsigned char) + sizeof(unsigned short)), DataStruct::Pmap);
    prefix_key key = {pre_key};
    return std::make_pair(key, ordered);
}

PrefixMap::iterator PrefixPermutationMap::find_prefix_key_in_pmap(prefix_key key) {
    std::unique_lock<std::mutex> key_lk(key_lk_);
    PrefixLocks::iterator key_iter = active_keys.find(key);
    // Wait for other thread to finish with current entry
    while (key_iter != active_keys.end()) {
        // TODO add counter to measure how many collisions there are
        key_cv.wait(key_lk);
        key_iter = active_keys.find(key);
    }
    std::pair<prefix_key, bool> active_key = std::make_pair(key, true);
    active_keys.insert(active_key);
    key_lk.unlock();
    map_lk.lock();
    PrefixMap::iterator iter = pmap->find(key);
    map_lk.unlock();
    return iter;
}

void PrefixPermutationMap::remove_existing_node(CacheTree *tree, PrefixMap::iterator iter, tracking_vector<unsigned short, DataStruct::Tree> prefix,
                                                unsigned short thread_id) {
    // If permutation is better than existing node
    Node *permuted_node;
    // Get existing node and garbage collect
    tree->lock(thread_id);
    tracking_vector<unsigned short, DataStruct::Tree> permuted_prefix(prefix.size());
    unsigned char *indices = std::get<1>(iter->second);
    for (unsigned char i = 0; i < indices[0]; ++i)
        permuted_prefix[i] = prefix[indices[i + 1]];
    if ((permuted_node = tree->check_prefix(permuted_prefix)) != NULL) {
        if (featureDecisions->do_garbage_collection()) {
            Node *permuted_parent = permuted_node->parent();
            permuted_parent->delete_child(permuted_node->id());
        }
        permuted_node->lock();
        permuted_node->set_deleted();
        permuted_node->unlock();
        logger->incPmapDiscardNum();
    } else {
        logger->incPmapNullNum();
    }
    tree->unlock(thread_id);
}

Node *CapturedPermutationMap::insert(unsigned short new_rule, size_t nrules, bool prediction,
                                     bool default_prediction, double lower_bound, double objective, Node *parent, int num_not_captured,
                                     int nsamples, int len_prefix, double c, double equivalent_minority, CacheTree *tree,
                                     VECTOR not_captured, tracking_vector<unsigned short, DataStruct::Tree> prefix, unsigned short thread_id) {
    logger->incPermMapInsertionNum();
    prefix.push_back(new_rule);
    Node *child = NULL;
    captured_key key;
    rule_vinit(nsamples, &key.key);
    rule_copy(key.key, not_captured, nsamples);
#ifndef GMP
    key.len = (short)nsamples;
#endif
    //    lk_.lock();
    CapturedMap::iterator iter = pmap->find(key);
    //    lk_.unlock();
    if (iter != pmap->end()) {
        double permuted_lower_bound = iter->second.first;
        tracking_vector<unsigned short, DataStruct::Tree> permuted_prefix = iter->second.second;
        if (lower_bound < permuted_lower_bound) {
            Node *permuted_node;
            if ((permuted_node = tree->check_prefix(permuted_prefix)) != NULL) {
                if (featureDecisions->do_garbage_collection()) {
                    Node *permuted_parent = permuted_node->parent();
                    permuted_parent->delete_child(permuted_node->id());
                }
                logger->incPmapDiscardNum();
            } else {
                logger->incPmapNullNum();
            }
            child = tree->construct_node(new_rule, nrules, prediction, default_prediction,
                                         lower_bound, objective, parent,
                                         num_not_captured, nsamples, len_prefix, c, equivalent_minority);
            iter->second = std::make_pair(lower_bound, prefix);
        }
    } else {
        child = tree->construct_node(new_rule, nrules, prediction, default_prediction,
                                     lower_bound, objective, parent,
                                     num_not_captured, nsamples, len_prefix, c, equivalent_minority);
        //        lk_.lock();
        pmap->insert(std::make_pair(key, std::make_pair(lower_bound, prefix)));
        //        lk_.unlock();
        logger->incPmapSize();
    }
    return child;
}
