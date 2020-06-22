#pragma once

#include "alloc.hh"
#include "cache.hh"
#include "utils.hh"
#include <algorithm>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <set>
#include <unordered_map>

/*
 * Represent prefix canonical order using an array of shorts.
 * The 0th index of the pointer contains the length of the prefix.
 */
struct prefix_key {
    unsigned short *key;
    // Only for debugging
    // prefix_key(prefix_key const&)=delete;
    //prefix_key(prefix_key&&o):key(o.key){key=nullptr;}

    prefix_key(tracking_vector<unsigned short, DataStruct::Tree> prefix, unsigned int len_prefix) {
        key = (unsigned short *)malloc(sizeof(unsigned short) * (len_prefix + 1));
        key[0] = (unsigned short)len_prefix;
        memcpy(&key[1], &prefix[0], len_prefix * sizeof(unsigned short));
    }
    bool operator==(const prefix_key &other) const {
        // i = 0 checks for equivalent sizes, i > 0 checks for equivalent prefixes
        for (size_t i = 0; i <= key[0]; ++i) {
            if (key[i] != other.key[i])
                return false;
        }
        return true;
    }
    ~prefix_key() {
        free(key);
    }
};

struct prefix_val {
    prefix_val(double lb, unsigned char *ord, size_t ord_sz, size_t tid) {
        lower_bound = lb;
        ordered = (unsigned char *)malloc(sizeof(unsigned char) * ord_sz);
        memcpy(ordered, ord, sizeof(unsigned char) * ord_sz);
        thread_id = tid;
    }
    ~prefix_val() {
        free(ordered);
    }
    double lower_bound;
    size_t thread_id;
    unsigned char *ordered;
};

struct prefix_eq {
    bool operator()(const prefix_key* k, const prefix_key* other) const {
        // i = 0 checks for equivalent sizes, i > 0 checks for equivalent prefixes
        for (size_t i = 0; i <= k->key[0]; ++i) {
            if (k->key[i] != other->key[i])
                return false;
        }
        return true;
    }
};

/*
 * Hash function from: http://www.cse.yorku.ca/~oz/hash.html
 */
struct prefix_hash {
    std::size_t operator()(const prefix_key* k) const {
        unsigned long hash = 0;
        for (size_t i = 1; i <= (k->key)[0]; ++i)
            hash = k->key[i] + (hash << 6) + (hash << 16) - hash;
        return hash;
    }
};

// Prefix Map typdefs
typedef std::vector<prefix_key*> PrefixLocks;
typedef std::unordered_map<prefix_key*, prefix_val*, prefix_hash, prefix_eq, track_alloc<std::pair<prefix_key* const, prefix_val*>, DataStruct::Pmap>> PrefixMap;

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
    bool operator()(const captured_key &k, const captured_key &other) const {
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
    size_t operator()(const captured_key &k) const {
#ifdef GMP
        return rule_vector_hash(k.key, 0);
#else
        return rule_vector_hash(k.key, k.len);
#endif
    }
};

// Prefix Map typdefs
//typedef std::unordered_map<prefix_key, bool> prefix_locks;
typedef std::pair<double, tracking_vector<unsigned short, DataStruct::Tree>> cap_val;
typedef std::unordered_map<captured_key, cap_val, captured_hash, cap_eq, track_alloc<std::pair<const captured_key, cap_val>, DataStruct::Pmap>> CapturedMap;

class PermutationMap {
public:
    virtual ~PermutationMap(){};
    virtual size_t size() { return 0; }
    virtual Node *insert(unsigned short new_rule,
                         size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                         double objective, Node *parent, int num_not_captured, int nsamples, int len_prefix,
                         double c, double equivalent_minority, CacheTree *tree, VECTOR not_captured,
                         tracking_vector<unsigned short, DataStruct::Tree> prefix, unsigned short thread_id) { return NULL; }
    virtual bool prefix_exists_and_is_worse(tracking_vector<unsigned short, DataStruct::Tree> prefix, double lb) { return false; }
};

class PrefixPermutationMap : public PermutationMap {
public:
    PrefixPermutationMap();
    ~PrefixPermutationMap();
    size_t size() override {
        map_lk.lock();
        size_t sz = pmap->size();
        map_lk.unlock();
        return sz;
    }
    Node *insert(unsigned short new_rule, size_t nrules, bool prediction,
                 bool default_prediction, double lower_bound, double objective, Node *parent,
                 int num_not_captured, int nsamples, int len_prefix, double c, double equivalent_minority,
                 CacheTree *tree, VECTOR not_captured, tracking_vector<unsigned short, DataStruct::Tree> prefix, unsigned short thread_id) override;
    bool prefix_exists_and_is_worse(tracking_vector<unsigned short, DataStruct::Tree> prefix, double lb) override;

private:
    prefix_key* construct_prefix_key_from_prefix_vector(tracking_vector<unsigned short, DataStruct::Tree> prefix);
    std::pair<prefix_key *, unsigned char *> construct_prefix_key_and_order_from_prefix_vector(tracking_vector<unsigned short, DataStruct::Tree> prefix,
                                                                                   unsigned short new_rule, int len_prefix);
    PrefixMap::iterator find_prefix_key_in_pmap(prefix_key *key);
    void remove_existing_node(CacheTree *tree, unsigned char *indices, tracking_vector<unsigned short, DataStruct::Tree> prefix,
                              unsigned short thread_id);

    PrefixLocks::iterator find_key_in_active_keys(prefix_key *key);

    PrefixMap *pmap;
    std::mutex map_lk;
    std::mutex key_lk_;
    std::condition_variable key_cv;
    PrefixLocks active_keys;
};

class CapturedPermutationMap : public PermutationMap {
public:
    CapturedPermutationMap();
    ~CapturedPermutationMap();
    size_t size() override {
        return pmap->size();
    }
    Node *insert(unsigned short new_rule, size_t nrules, bool prediction, bool default_prediction,
                 double lower_bound, double objective, Node *parent, int num_not_captured, int nsamples,
                 int len_prefix, double c, double equivalent_minority, CacheTree *tree, VECTOR not_captured,
                 tracking_vector<unsigned short, DataStruct::Tree> prefix, unsigned short thread_id) override;
    // TODO: implement this correctly
    bool prefix_exists_and_is_worse(tracking_vector<unsigned short, DataStruct::Tree> prefix, double lb) override { return false; }

private:
    CapturedMap *pmap;
    std::mutex key_lk;
    //        std::unordered_map<captured_key, bool> active_keys;
};

class NullPermutationMap : public PermutationMap {
public:
    ~NullPermutationMap(){};
    size_t size() override { return 0; }
    Node *insert(unsigned short new_rule, size_t nrules, bool prediction, bool default_prediction, double lower_bound,
                 double objective, Node *parent, int num_not_captured, int nsamples, int len_prefix,
                 double c, double equivalent_minority, CacheTree *tree, VECTOR not_captured,
                 tracking_vector<unsigned short, DataStruct::Tree> prefix, unsigned short thread_id) override {
        return tree->construct_node(new_rule, nrules, prediction, default_prediction,
                                    lower_bound, objective, parent,
                                    num_not_captured, nsamples, len_prefix, c, equivalent_minority);
    }
    bool prefix_exists_and_is_worse(tracking_vector<unsigned short, DataStruct::Tree> prefix, double lb) { return false; }
};
