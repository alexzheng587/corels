#include "pmap.hh"

PrefixPermutationMap::PrefixPermutationMap()
    : pmap(new PrefixMap) {
}

CapturedPermutationMap::CapturedPermutationMap()
    : pmap(new CapturedMap) {}

Node* PrefixPermutationMap::insert (unsigned short new_rule, size_t nrules, bool prediction, 
        bool default_prediction, double lower_bound, double objective, Node* parent, 
        int num_not_captured, int nsamples, int len_prefix, double c, double equivalent_minority,
        CacheTree* tree, VECTOR not_captured, tracking_vector<unsigned short, 
        DataStruct::Tree> parent_prefix, unsigned short thread_id) {
    (void) not_captured;
    logger->incPermMapInsertionNum();
    parent_prefix.push_back(new_rule);
    
    unsigned char* ordered = (unsigned char*) malloc(sizeof(unsigned char) * (len_prefix + 1));
    ordered[0] = (unsigned char)len_prefix;

    for (int i = 1; i < (len_prefix + 1); i++)
	    ordered[i] = i - 1;

    std::function<bool(int, int)> cmp = [&](int i, int j) { return parent_prefix[i] < parent_prefix[j]; };
    std::sort(&ordered[1], &ordered[len_prefix + 1], cmp);
    
    std::sort(parent_prefix.begin(), parent_prefix.end());
    unsigned short *pre_key = (unsigned short*) malloc(sizeof(unsigned short) * (len_prefix + 1));
    pre_key[0] = (unsigned short)len_prefix;
    memcpy(&pre_key[1], &parent_prefix[0], len_prefix * sizeof(unsigned short));
    
    logger->addToMemory((len_prefix + 1) * (sizeof(unsigned char) + sizeof(unsigned short)), DataStruct::Pmap);
    prefix_key key = { pre_key };
    
    Node* child = NULL;
    ++lock_ac;
    std::unique_lock<std::mutex> key_lk(key_lk_);
    PrefixLocks::iterator key_iter = active_keys.find(key);
    // Wait for other thread to finish with current entry
    while(key_iter != active_keys.end()) {
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
    if (iter != pmap->end()) {
        double permuted_lower_bound = std::get<0>(iter->second);
        if (lower_bound < permuted_lower_bound) {
            Node* permuted_node;
            tree->lock(thread_id);
            tracking_vector<unsigned short, DataStruct::Tree> permuted_prefix(parent_prefix.size());
            unsigned char* indices = std::get<1>(iter->second);
            for (unsigned char i = 0; i < indices[0]; ++i)
                permuted_prefix[i] = parent_prefix[indices[i + 1]];
            if ((permuted_node = tree->check_prefix(permuted_prefix)) != NULL) {
                Node* permuted_parent = permuted_node->parent();
                //permuted_parent->delete_child(permuted_node->id());
                //delete_subtree(tree, permuted_node, false, tree->calculate_size());
                //permuted_node->set_done();
                permuted_node->set_deleted();
                logger->incPmapDiscardNum();
            } else {
                logger->incPmapNullNum();
            }
            tree->unlock(thread_id);
            child = tree->construct_node(new_rule, nrules, prediction,
                                     default_prediction, lower_bound, objective,
                                     parent, num_not_captured, nsamples,
                                     len_prefix, c, equivalent_minority);
            iter->second = std::make_tuple(lower_bound, ordered, thread_id);
        }
    } else {
        child = tree->construct_node(new_rule, nrules, prediction,
                                 default_prediction, lower_bound, objective,
                                 parent, num_not_captured, nsamples, len_prefix,
                                 c, equivalent_minority);
        unsigned char* ordered_prefix = &ordered[0];
        // Need to globally lock map when inserting otherwise iterators in other threads could be
        // invalidated if the unordere_map is resized.
        map_lk.lock();
        pmap->insert(std::make_pair(key, std::make_tuple(lower_bound, ordered_prefix, thread_id)));
        map_lk.unlock();
        logger->incPmapSize();
    }
    // Wake up any other threads waiting on this entry
    key_lk.lock();
    active_keys.erase(key);
    key_lk.unlock();
    key_cv.notify_all();
    return child;
}

Node* CapturedPermutationMap::insert(unsigned short new_rule, size_t nrules, bool prediction, 
        bool default_prediction, double lower_bound, double objective, Node* parent, int num_not_captured, 
        int nsamples, int len_prefix, double c, double equivalent_minority, CacheTree* tree, 
        VECTOR not_captured, tracking_vector<unsigned short, DataStruct::Tree> parent_prefix, unsigned short thread_id) {
    logger->incPermMapInsertionNum();
    parent_prefix.push_back(new_rule);
    Node* child = NULL;
    captured_key key;
    rule_vinit(nsamples, &key.key);
    rule_copy(key.key, not_captured, nsamples);
#ifndef GMP
    key.len = (short) nsamples;
#endif
//    lk_.lock();
    CapturedMap::iterator iter = pmap->find(key);
//    lk_.unlock();
    if (iter != pmap->end()) {
        double permuted_lower_bound = iter->second.first;
        tracking_vector<unsigned short, DataStruct::Tree> permuted_prefix = iter->second.second;
        if (lower_bound < permuted_lower_bound) {
            Node* permuted_node;
            if ((permuted_node = tree->check_prefix(permuted_prefix)) != NULL) {
                Node* permuted_parent = permuted_node->parent();
                permuted_parent->delete_child(permuted_node->id());
                //delete_subtree(tree, permuted_node, false, tree->calculate_size());
                logger->incPmapDiscardNum();
            } else {
                logger->incPmapNullNum();
            }
            child = tree->construct_node(new_rule, nrules, prediction, default_prediction,
                                       lower_bound, objective, parent,
                                        num_not_captured, nsamples, len_prefix, c, equivalent_minority);
            iter->second = std::make_pair(lower_bound, parent_prefix);
        }
    } else {
        child = tree->construct_node(new_rule, nrules, prediction, default_prediction,
                                    lower_bound, objective, parent,
                                    num_not_captured, nsamples, len_prefix, c, equivalent_minority);
//        lk_.lock();
        pmap->insert(std::make_pair(key, std::make_pair(lower_bound, parent_prefix)));
//        lk_.unlock();
        logger->incPmapSize();
    }
    return child;
}
