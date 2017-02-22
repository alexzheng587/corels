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

typedef std::unordered_map<struct prefix_key, std::pair<double, unsigned char*>, prefix_hash> PrefixPermutationMap;
typedef std::unordered_map<struct captured_key, std::pair<std::vector<unsigned short>, double>, captured_hash> CapturedPermutationMap;

template<class P>
using pmap_garbage_collect_signature = void (*)(P*, size_t);

void bfs_prefix_map_garbage_collect(PrefixPermutationMap* p, size_t min_length);
void prefix_map_garbage_collect(PrefixPermutationMap* p, size_t min_length);
void captured_map_garbage_collect(CapturedPermutationMap* p, size_t min_length);

template<class N, class P>
using permutation_insert_signature = N* (*)(construct_signature<N>, unsigned short, size_t, bool, bool, 
                                            double, double, N* parent, int, int, int, double, double, CacheTree<N>*, VECTOR,
                                            std::vector<unsigned short>, P*);

template<class N>
N* prefix_permutation_insert(construct_signature<N> construct_policy, unsigned short new_rule,
                        size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                        double objective, N* parent, int num_not_captured, int nsamples, int len_prefix,
                        double c, double minority, CacheTree<N>* tree, VECTOR not_captured,
                        std::vector<unsigned short>, PrefixPermutationMap* p);

template<class N>
N* captured_permutation_insert(construct_signature<N> construct_policy, unsigned short new_rule,
                        size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                        double objective, N* parent, int num_not_captured, int nsamples, int len_prefix,
                        double c, double minority, CacheTree<N>* tree, VECTOR not_captured,
                        std::vector<unsigned short>, CapturedPermutationMap* p);
