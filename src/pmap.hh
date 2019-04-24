#pragma once
#include "cache.hh"
#include "utils.hh"
#include "alloc.hh"
#include <algorithm>
#include <mutex>
#include <set>
#include <unordered_map>
#include <functional>
#include <condition_variable>
#include <set>

extern int lock_ac;
//extern pthread_rwlock_t pmap_lk_;

/*
 * Represent prefix canonical order using an array of shorts.
 * The 0th index of the pointer contains the length of the prefix.
 */
struct prefix_key {
    unsigned short *key;
    prefix_key(unsigned short* k) {
        key = k;
    }
};

struct prefix_eq {
   bool operator()(const prefix_key& k, const prefix_key& other) const {
        // i = 0 checks for equivalent sizes, i > 0 checks for equivalent prefixes
        for(size_t i = 0; i <= k.key[0]; ++i) {
            if (k.key[i] != other.key[i])
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

// Prefix Map typdefs
typedef std::unordered_map<prefix_key, bool, prefix_hash, prefix_eq, track_alloc<std::pair<const prefix_key, bool>, DataStruct::Pmap> > PrefixLocks;
typedef std::tuple<double, unsigned char*, size_t> prefix_val;
typedef std::unordered_map<prefix_key, prefix_val, prefix_hash, prefix_eq, track_alloc<std::pair<const prefix_key, prefix_val>, DataStruct::Pmap> > PrefixMap;

/*
 * Represents captured vector using the VECTOR type defined in rule.h
 */
struct captured_key {
    VECTOR key;
#ifndef GMP
    short len;
#endif
};

struct cap_eq {
    bool operator()(const captured_key& k, const captured_key& other) const {
#ifdef GMP
        return rule_vector_equal(k.key, other.key, 0, 0);
#else
        return rule_vector_equal(k.key, other.key, k.len, other.len);
#endif
    }
};

/*
 * Hash function from: http://www.cse.yorku.ca/~oz/hash.html
 */
struct captured_hash {
    size_t operator()(const captured_key& k) const{
#ifdef GMP
        return rule_vector_hash(k.key, 0);
#else
        return rule_vector_hash(k.key, k.len);
#endif
    }
};

// Prefix Map typdefs
//typedef std::unordered_map<prefix_key, bool> prefix_locks;
typedef std::pair<double, tracking_vector<unsigned short, DataStruct::Tree> > cap_val;
typedef std::unordered_map<captured_key, cap_val, captured_hash, cap_eq, track_alloc<std::pair<const captured_key, cap_val>, DataStruct::Pmap> > CapturedMap;

class PermutationMap {
    public:
        virtual size_t size() { return 0; }
        virtual Node* insert (unsigned short new_rule,
                             size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                             double objective, Node* parent, int num_not_captured, int nsamples, int len_prefix,
                             double c, double equivalent_minority, CacheTree* tree, VECTOR not_captured,
                             tracking_vector<unsigned short, DataStruct::Tree> parent_prefix, unsigned short thread_id) { return NULL; }
        Node* check_permutation_bound (unsigned short new_rule,
                             size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                             double objective, Node* parent, int num_not_captured, int nsamples, int len_prefix,
                             double c, double equivalent_minority, CacheTree* tree, VECTOR not_captured,
                             tracking_vector<unsigned short, DataStruct::Tree> parent_prefix);
};

class PrefixPermutationMap : public PermutationMap {
	public:
        PrefixPermutationMap ();
        size_t size() override {
            return pmap->size();
        }
        Node* insert (unsigned short new_rule, size_t nrules, bool prediction, 
            bool default_prediction, double lower_bound, double objective, Node* parent, 
            int num_not_captured, int nsamples, int len_prefix, double c, double equivalent_minority,
            CacheTree* tree, VECTOR not_captured, tracking_vector<unsigned short, 
            DataStruct::Tree> parent_prefix, unsigned short thread_id) override;
	private:
		PrefixMap* pmap;
        std::mutex map_lk;
        std::mutex key_lk_;
        std::condition_variable key_cv;
        PrefixLocks active_keys;
};

class CapturedPermutationMap : public PermutationMap {
	public:
        CapturedPermutationMap();
        size_t size() override {
            return pmap->size();
        }
        Node* insert(unsigned short new_rule, size_t nrules, bool prediction, bool default_prediction, 
                double lower_bound, double objective, Node* parent, int num_not_captured, int nsamples, 
                int len_prefix, double c, double equivalent_minority, CacheTree* tree, VECTOR not_captured,
                 tracking_vector<unsigned short, DataStruct::Tree> parent_prefix, unsigned short thread_id) override;
	private:
		CapturedMap* pmap;
        std::mutex key_lk;
//        std::unordered_map<captured_key, bool> active_keys;
};

class NullPermutationMap : public PermutationMap  {
    public:
        size_t size() override {return 0;}
        Node* insert (unsigned short new_rule, size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                        double objective, Node* parent, int num_not_captured, int nsamples, int len_prefix,
                        double c, double equivalent_minority, CacheTree* tree, VECTOR not_captured,
                        tracking_vector<unsigned short, DataStruct::Tree> parent_prefix, unsigned short thread_id) override {
            return tree->construct_node(new_rule, nrules, prediction, default_prediction,
                                            lower_bound, objective, parent,
                                            num_not_captured, nsamples, len_prefix, c, equivalent_minority);
        }
};
