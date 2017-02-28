#pragma once
#include "cache.hh"
#include "utils.hh"
#include <unordered_map>

/*
 * Represent prefix canonical order using an array of shorts.
 * The 0th index of the pointer contains the length of the prefix.
 */
struct prefix_key {
    unsigned short *key;
    prefix_key(unsigned short* k) {
        key = k;
    }
    bool operator==(const prefix_key& other) const {
        if (key[0] != other.key[0])
            return false;
        for(size_t i = 1; i <= *key; ++i) {
            if (key[i] != other.key[i])
                return false;
        }
        return true;
    }
};

/*
 * Hash function from: http://www.cse.yorku.ca/~oz/hash.html
 */
struct prefix_hash {
    std::size_t operator()(const prefix_key& k) const {
        unsigned long hash = 0;
        for(size_t i = 1; i <= *k.key; ++i)
            hash = k.key[i] + (hash << 6) + (hash << 16) - hash;
        return hash;
    }
};

/*
 * Represents captured vector using the VECTOR type defined in rule.h
 */
struct captured_key {
    VECTOR key;
#ifndef GMP
    short len;
#endif
    bool operator==(const captured_key& other) const {
#ifdef GMP
        return !mpz_cmp(other.key, key);
#else
        size_t nentries = (len + BITS_PER_ENTRY - 1)/BITS_PER_ENTRY;
        for (size_t i = 0; i < nentries; i++)
            if (other.key[i] != key[i])
                return false;
        return true;
#endif
    }
};

/*
 * Hash function from: http://www.cse.yorku.ca/~oz/hash.html
 */
struct captured_hash {
    std::size_t operator()(const captured_key& k) const{
        unsigned long hash = 0;
		size_t i;
#ifdef GMP
        for(i = 0; i < mpz_size(k.key); ++i)
            hash = mpz_getlimbn(k.key, i) + (hash << 6) + (hash << 16) - hash;
#else
   		size_t nentries = (k.len + BITS_PER_ENTRY - 1)/BITS_PER_ENTRY;
		for(i = 0; i < nentries; ++i)
            hash = k.key[i] + (hash << 6) + (hash << 16) - hash;
#endif
        return hash;
    }
};

/*template<class N>
using permutation_insert_signature = N* (*)(construct_signature<N>, unsigned short, size_t, bool, bool, 
                                            double, double, N* parent, int, int, int, double, double, CacheTree<N>*, VECTOR,
                                            std::vector<unsigned short>, PermutationMap*);
                                            */

class PermutationMap {
    public:
        virtual size_t size() {
            return 0;
        }
        virtual Node* insert (unsigned short new_rule,
                             size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                             double objective, Node* parent, int num_not_captured, int nsamples, int len_prefix,
                             double c, double minority, CacheTree* tree, VECTOR not_captured,
                             std::vector<unsigned short> parent_prefix)
        {
            (void) not_captured, (void) parent_prefix;
            return tree->construct_node(new_rule, nrules, prediction,
                                    default_prediction, lower_bound, objective,
                                    parent, num_not_captured, nsamples,
                                    len_prefix, c, minority);
        }
};

class PrefixPermutationMap : public PermutationMap {
	public:
        PrefixPermutationMap ();
        size_t size() override {
            return pmap->size();
        }
        Node* insert (unsigned short new_rule,
					 size_t nrules, bool prediction, bool default_prediction, double lower_bound,
					 double objective, Node* parent, int num_not_captured, int nsamples, int len_prefix,
					 double c, double minority, CacheTree* tree, VECTOR not_captured,
					 std::vector<unsigned short> parent_prefix) override
		{
            (void) not_captured;
            parent_prefix.push_back(new_rule);

            unsigned char* ordered = (unsigned char*) malloc(sizeof(unsigned char) * (len_prefix + 1));
            ordered[0] = (unsigned char)len_prefix;
            auto cmp = [&](int i, int j) { return parent_prefix[i] < parent_prefix[j]; };
            std::sort(&ordered[1], &ordered[len_prefix], cmp);

            std::sort(parent_prefix.begin(), parent_prefix.end());
            
            unsigned short *pre_key = (unsigned short*) malloc(sizeof(unsigned short) * (len_prefix + 1));
            pre_key[0] = (unsigned short)len_prefix;
            memcpy(&pre_key[1], &parent_prefix[0], len_prefix * sizeof(unsigned short));
            
            prefix_key key = { pre_key };
            Node* child = NULL;
            auto iter = pmap->find(key);
            if (iter != pmap->end()) {
                double permuted_lower_bound = iter->second.first;
                if (lower_bound < permuted_lower_bound) {
                    Node* permuted_node;
                    unsigned char* indices = iter->second.second;
                    std::vector<unsigned short> permuted_prefix(parent_prefix.size());
                    for (unsigned char i = 0; i < indices[0]; ++i)
                        permuted_prefix[i] = parent_prefix[indices[i + 1]];

                    if ((permuted_node = tree->check_prefix(permuted_prefix)) != NULL) {
                        Node* permuted_parent = permuted_node->parent();
                        permuted_parent->delete_child(permuted_node->id());
                        delete_subtree(tree, permuted_node, false, true);
                        logger.incPmapDiscardNum();
                    } else {
                        logger.incPmapNullNum();
                    }
                    child = tree->construct_node(new_rule, nrules, prediction,
                                             default_prediction, lower_bound, objective,
                                             parent, num_not_captured, nsamples,
                                             len_prefix, c, minority);
                    iter->second = std::make_pair(lower_bound, ordered);
                }
            } else {
                child = tree->construct_node(new_rule, nrules, prediction,
                                         default_prediction, lower_bound, objective,
                                         parent, num_not_captured, nsamples, len_prefix,
                                         c, minority);
                unsigned char* ordered_prefix = &ordered[0];
                pmap->insert(std::make_pair(key, std::make_pair(lower_bound, ordered_prefix)));
                logger.incPmapSize();
            }
            return child;
        }
	private:
		std::unordered_map<prefix_key, std::pair<double, unsigned char*>, prefix_hash>* pmap;
};

class CapturedPermutationMap : public PermutationMap {
	public:
        CapturedPermutationMap();
        size_t size() override {
            return cmap->size();
        }
        Node* insert(unsigned short new_rule,
					 size_t nrules, bool prediction, bool default_prediction, double lower_bound,
					 double objective, Node* parent, int num_not_captured, int nsamples, int len_prefix,
					 double c, double minority, CacheTree* tree, VECTOR not_captured,
					 std::vector<unsigned short> parent_prefix) override
		{
            parent_prefix.push_back(new_rule);
            Node* child = NULL;
            captured_key key;
            rule_vinit(nsamples, &key.key);
            rule_copy(key.key, not_captured, nsamples);
#ifndef GMP
            key.len = (short) nsamples;
#endif
            auto iter = cmap->find(key);
            if (iter != cmap->end()) {
                std::vector<unsigned short> permuted_prefix = iter->second.first;
                double permuted_lower_bound = iter->second.second;
                if (lower_bound < permuted_lower_bound) {
                    Node* permuted_node;
                    if ((permuted_node = tree->check_prefix(permuted_prefix)) != NULL) {
                        Node* permuted_parent = permuted_node->parent();
                        permuted_parent->delete_child(permuted_node->id());
                        delete_subtree(tree, permuted_node, false, true);
                        logger.incPmapDiscardNum();
                    } else {
                        logger.incPmapNullNum();
                    }
                    child = tree->construct_node(new_rule, nrules, prediction, default_prediction,
                                               lower_bound, objective, parent,
                                                num_not_captured, nsamples, len_prefix, c, minority);
                    iter->second = std::make_pair(parent_prefix, lower_bound);
                }
            } else {
                child = tree->construct_node(new_rule, nrules, prediction, default_prediction,
                                            lower_bound, objective, parent,
                                            num_not_captured, nsamples, len_prefix, c, minority);
                cmap->insert(std::make_pair(key, std::make_pair(parent_prefix, lower_bound)));
                logger.incPmapSize();
            }
            return child;
		}
 
	private:
		std::unordered_map<captured_key, std::pair<std::vector<unsigned short>, double>, captured_hash>* cmap;
};

class NullPermutationMap : public PermutationMap  {
    public:
        size_t size() {return 0;}
};

//typedef std::unordered_map<struct prefix_key, std::pair<double, unsigned char*>, prefix_hash> PrefixPermutationMap;
//typedef std::unordered_map<struct captured_key, std::pair<std::vector<unsigned short>, double>, captured_hash> CapturedPermutationMap;

//using pmap_garbage_collect_signature = void (*)(PermutationMap*, size_t);

//void bfs_prefix_map_garbage_collect(PermutationMap* p, size_t min_length);
//void prefix_map_garbage_collect(PermutationMap* p, size_t min_length);
//void captured_map_garbage_collect(PermutationMap* p, size_t min_length);
