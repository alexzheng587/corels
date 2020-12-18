#pragma once

#include "cache.hh"
#include "utils.hh"
#include "alloc.hh"
#include <unordered_map>
#include <algorithm>
#include <functional>
#include <set>
#include <vector>

/*
 * Represent prefix canonical order using an array of shorts.
 * The 0th index of the pointer contains the length of the prefix.
 */

struct prefix_key {
    unsigned short* key;

    prefix_key(unsigned short* k) {
        key = k;
    }
};

struct prefix_eq {
    bool operator()(const prefix_key &k, const prefix_key &other) const {
        // i = 0 checks for equivalent sizes, i > 0 checks for equivalent prefixes
        for (size_t i = 0; i <= k.key[0]; ++i) {
            if (k.key[i] != other.key[i])
                return false;
        }
        return true;
    }
};

typedef std::pair<double, unsigned char*> prefix_val;

/*
 * Hash function from: http://www.cse.yorku.ca/~oz/hash.html
 */
struct prefix_hash {
    std::size_t operator()(const prefix_key &k) const {
        unsigned long hash = 0;
        for (size_t i = 1; i <= *k.key; ++i)
            hash = k.key[i] + (hash << 6) + (hash << 16) - hash;
        return hash;
    }
};

typedef std::unordered_map<prefix_key, prefix_val, prefix_hash, prefix_eq, track_alloc<std::pair<const prefix_key, prefix_val>, DataStruct::Pmap> > PrefixMap;

typedef std::vector<F1Node*> fscore_val;

typedef std::unordered_map<prefix_key, fscore_val, prefix_hash, prefix_eq, track_alloc<std::pair<const prefix_key, fscore_val>, DataStruct::Pmap> > FScoreMap;

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

typedef std::pair<double, tracking_vector<unsigned short, DataStruct::Tree> > cap_val;

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

typedef std::unordered_map<captured_key, cap_val, captured_hash, cap_eq, track_alloc<std::pair<const captured_key, cap_val>, DataStruct::Pmap> > CapturedMap;

class PermutationMap {
public:
    virtual size_t size() = 0;

    virtual Node* insert(unsigned short new_rule, ProblemState* state, Node* parent, CacheTree* tree,
           tracking_vector<unsigned short, DataStruct::Tree> parent_prefix) = 0;

    Node* check_permutation_bound(unsigned short new_rule, size_t nrules, bool prediction, bool default_prediction,
                                  double lower_bound, double objective, Node* parent, int num_not_captured,
                                  int nsamples, int len_prefix, double c, double equivalent_minority, CacheTree* tree,
                                  VECTOR not_captured, tracking_vector<unsigned short, DataStruct::Tree> parent_prefix);

    virtual ~PermutationMap() = default;
};

class PrefixPermutationMap : public PermutationMap {
public:
    PrefixPermutationMap();

    size_t size() override {
        return pmap->size();
    }

    Node* insert(unsigned short new_rule, ProblemState* state, Node* parent, CacheTree* tree,
                 tracking_vector<unsigned short, DataStruct::Tree> parent_prefix) override;

    inline PrefixMap* getMap() const {
        return pmap;
    }

private:
    PrefixMap* pmap;
};

class CapturedPermutationMap : public PermutationMap {
public:
    CapturedPermutationMap();

    size_t size() override {
        return pmap->size();
    }

    Node* insert(unsigned short new_rule, ProblemState* state, Node* parent, CacheTree* tree,
                 tracking_vector<unsigned short, DataStruct::Tree> parent_prefix) override;

    inline CapturedMap* getMap() const {
        return pmap;
    }

private:
    CapturedMap* pmap;
};

class NullPermutationMap : public PermutationMap {
public:
    size_t size() override { return 0; }

    Node* insert(unsigned short new_rule, ProblemState* state, Node* parent, CacheTree* tree,
                 tracking_vector<unsigned short, DataStruct::Tree> parent_prefix) override {
        Node* child = parent->construct_child(new_rule, state);
        logger->addToMemory(sizeof(*child), DataStruct::Tree);
        return child;
    }
};

class FScorePermutationMap : public PermutationMap {
public:
    size_t size() override {
        return pmap->size();
    }

    Node* insert(unsigned short new_rule, ProblemState* state, Node* parent, CacheTree* tree,
                 tracking_vector<unsigned short, DataStruct::Tree> parent_prefix) override;

private:
    FScoreMap* pmap;
};

class AUCPermutationMap : public PermutationMap {
public:
    AUCPermutationMap();

    size_t size() override {
        return pmap->size();
    }

    Node* insert(unsigned short new_rule, ProblemState* state, Node* parent, CacheTree* tree,
                 tracking_vector<unsigned short, DataStruct::Tree> parent_prefix) override;

private:
    PrefixMap* pmap;
};