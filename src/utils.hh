#pragma once

#include "rule.h"

#include <cstdlib>
#include <sys/time.h>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <gmpxx.h>

using namespace std;

/*
template <class T>  
    struct track_alloc { 
    typedef T value_type;
    track_alloc() noexcept {}
    template <class U> track_alloc (const track_alloc<U>&) noexcept {}
    T* allocate (size_t n) {
        printf("SZ N: %zu\n", n * sizeof(T));
        return static_cast<T*>(malloc(n*sizeof(T)));
    }   
    void deallocate (T* p, size_t n) {
        printf("DEALLOCATE\n");
        free(p);
    }   
};
*/

class NullLogger {
  public:
    virtual void closeFile() {}
    ~NullLogger() {}

    virtual void setLogFileName(char *fname) {}
    virtual void dumpState() {}
    virtual std::string dumpPrefixLens() {}
    virtual std::string dumpRemainingSpaceSize() {}
    //virtual size_t dumpLogRemainingSpaceSize() {}

    virtual inline void setVerbosity(int verbosity) {}
    virtual inline void setFrequency(int frequency) {}
    virtual inline int getFrequency() {}
    virtual inline void addToLowerBoundTime(double t) {}
    virtual inline void incLowerBoundNum() {}
    virtual inline void addToObjTime(double t) {}
    virtual inline void incObjNum() {}
    virtual inline void addToTreeInsertionTime(double t) {}
    virtual inline void incTreeInsertionNum() {}
    virtual inline void addToRuleEvalTime(double t) {}
    virtual inline void incRuleEvalNum() {}
    virtual inline void addToNodeSelectTime(double t) {}
    virtual inline void incNodeSelectNum() {}
    virtual inline void addToEvalChildrenTime(double t) {}
    virtual inline void incEvalChildrenNum() {}
    virtual inline void setInitialTime(double t) {}
    virtual inline double getInitialTime() {}
    virtual inline void setTotalTime(double t) {}
    virtual inline void addToPermMapInsertionTime(double t) {}
    virtual inline void incPermMapInsertionNum() {}
    virtual inline void setCurrentLowerBound(double lb) {}
    virtual inline void setTreeMinObj(double o) {}
    virtual inline void setTreePrefixLen(size_t n) {}
    virtual inline void setTreeNumNodes(size_t n) {}
    virtual inline void setTreeNumEvaluated(size_t n) {}
    virtual inline void addToTreeMemory(size_t n) {}
    virtual inline void removeFromTreeMemory(size_t n) {}
    virtual inline std::pair<size_t, size_t> getTreeMemory() {}
    virtual inline void addToQueueInsertionTime(double t) {}
    virtual inline void setQueueSize(size_t n) {}
    virtual inline void addToQueueMemory(size_t n) {}
    virtual inline void removeFromQueueMemory(size_t n) {}
    virtual inline std::pair<size_t, size_t> getQueueMemory() {}
    virtual inline void setNRules(size_t nrules) {}
    virtual inline void setC(double c) {}
    virtual inline void initPrefixVec() {}
    virtual inline void incPrefixLen(size_t n) {}
    virtual inline void decPrefixLen(size_t n) {}
    virtual inline size_t sumPrefixLens() {}
    virtual inline void updateQueueMinLen() {}
    virtual inline size_t getQueueMinLen() {}
    virtual inline void incPmapSize() {}
    virtual inline void decreasePmapSize(size_t n) {}
    virtual inline void incPmapNullNum() {}
    virtual inline void incPmapDiscardNum() {}
    virtual inline void addToPmapMemory(size_t n) {}
    virtual inline void removeFromPmapMemory(size_t n) {}
    virtual inline std::pair<size_t, size_t> getPmapMemory() {}
    virtual inline void subtreeSize(mpz_t tot, unsigned int len_prefix, double lower_bound) {}
    virtual inline void approxRemainingSize(mpz_t tot, unsigned int len_prefix) {}
    virtual inline void addQueueElement(unsigned int len_prefix, double lower_bound, bool approx) {}
    virtual inline void removeQueueElement(unsigned int len_prefix, double lower_bound, bool approx) {}
    virtual inline void initRemainingSpaceSize() {}
    virtual inline void clearRemainingSpaceSize() {}
    virtual inline size_t getLogRemainingSpaceSize() {}
    /*
     * Initializes the logger by setting all fields to 0.
     * This allwos us to write a log record immediately.
     */
    inline void initializeState() {
        _state.total_time = 0.;
        _state.evaluate_children_time = 0.;
        _state.evaluate_children_num = 0;
        _state.node_select_time = 0.;
        _state.node_select_num = 0;
        _state.rule_evaluation_time = 0.;
        _state.rule_evaluation_num = 0;
        _state.lower_bound_time = 0.;
        _state.lower_bound_num = 0;
        _state.objective_time = 0.;
        _state.objective_num = 0;
        _state.tree_insertion_time = 0.;
        _state.tree_insertion_num = 0;
        _state.permutation_map_insertion_time = 0.;
        _state.permutation_map_insertion_num = 0;
        _state.current_lower_bound = 0.;
        _state.tree_min_objective = 1.;
        _state.tree_prefix_length = 0;
        _state.tree_num_nodes = 0;
        _state.tree_num_evaluated = 0;
        _state.tree_memory = 0;
        _state.tree_max_memory = 0;
        _state.queue_insertion_time = 0;
        _state.queue_size = 0;
        _state.queue_min_length = 0;
        _state.queue_memory = 0;
        _state.queue_max_memory = 0;
        _state.pmap_size = 0;
        _state.pmap_memory = 0;
        _state.pmap_max_memory = 0;
        _state.pmap_null_num = 0;
        _state.pmap_discard_num = 0;
        mpz_init(_state.remaining_space_size);
        initRemainingSpaceSize();
    }


  protected:
    struct State {
        size_t nrules;
        double c;
        double initial_time;                    // initial time stamp
        double total_time;
        double evaluate_children_time;
        size_t evaluate_children_num;
        double node_select_time;
        size_t node_select_num;
        double rule_evaluation_time;
        size_t rule_evaluation_num;
        double lower_bound_time;
        size_t lower_bound_num;
        double objective_time;
        size_t objective_num;
        double tree_insertion_time;
        size_t tree_insertion_num;
        double permutation_map_insertion_time;
        size_t permutation_map_insertion_num;   // number of calls to `permutation_insert` function
        double current_lower_bound;             // monotonically decreases for curious lower bound
        double tree_min_objective;
        size_t tree_prefix_length;
        size_t tree_num_nodes;
        size_t tree_num_evaluated;
        size_t tree_memory;
        size_t tree_max_memory;
        double queue_insertion_time;
        size_t queue_size;
        size_t queue_min_length;                // monotonically increases
        size_t queue_memory;
        size_t queue_max_memory;
        size_t pmap_size;                       // size of pmap
        size_t pmap_null_num;                   // number of pmap lookups that return null
        size_t pmap_discard_num;                // number of pmap lookups that trigger discard
        size_t pmap_memory;
        size_t pmap_max_memory;
        size_t* prefix_lens;
        mpz_t remaining_space_size;
    };
    State _state;
    int _v;                                     // verbosity
    int _freq;                                  // frequency of logging
    ofstream _f;                                // output file
};

class Logger : public NullLogger {
  public:
    void closeFile() override { if (_f.is_open()) _f.close(); }
    ~Logger() { 
        free(_state.prefix_lens);
        closeFile(); 
    }

    void setLogFileName(char *fname) override;
    void dumpState() override;
    std::string dumpPrefixLens() override;
    std::string dumpRemainingSpaceSize() override;
    //size_t dumpLogRemainingSpaceSize() override;

    inline void setVerbosity(int verbosity) override {
        _v = verbosity;
    }
    inline void setFrequency(int frequency) override {
        _freq = frequency;
    }
    inline int getFrequency() override {
        return _freq;
    }
    inline void addToLowerBoundTime(double t) override {
        _state.lower_bound_time += t;
    }
    inline void incLowerBoundNum() override {
        ++_state.lower_bound_num;
    }
    inline void addToObjTime(double t) override {
        _state.objective_time += t;
    }
    inline void incObjNum() override {
        ++_state.objective_num;
    }
    inline void addToTreeInsertionTime(double t) override {
        _state.tree_insertion_time += t;
    }
    inline void incTreeInsertionNum() override {
        ++_state.tree_insertion_num;
    }
    inline void addToRuleEvalTime(double t) override {
        _state.rule_evaluation_time += t;
    }
    inline void incRuleEvalNum() override {
        ++_state.rule_evaluation_num;
    }
    inline void addToNodeSelectTime(double t) override {
        _state.node_select_time += t;
    }
    inline void incNodeSelectNum() override {
        ++_state.node_select_num;
    }
    inline void addToEvalChildrenTime(double t) override {
        _state.evaluate_children_time += t;
    }
    inline void incEvalChildrenNum() override {
        ++_state.evaluate_children_num;
    }
    inline void setInitialTime(double t) override {
        _state.initial_time = t;
    }
    inline double getInitialTime() override {
        return _state.initial_time;
    }
    inline void setTotalTime(double t) override {
        _state.total_time = t;
    }
    inline void addToPermMapInsertionTime(double t) override {
        _state.permutation_map_insertion_time += t;
    }
    inline void incPermMapInsertionNum() override {
        ++_state.permutation_map_insertion_num;
    }
    inline void setCurrentLowerBound(double lb) override {
        _state.current_lower_bound = lb;
    }
    inline void setTreeMinObj(double o) override {
        _state.tree_min_objective = o;
    }
    inline void setTreePrefixLen(size_t n) override {
        _state.tree_prefix_length = n;
    }
    inline void setTreeNumNodes(size_t n) override {
        _state.tree_num_nodes = n;
    }
    inline void setTreeNumEvaluated(size_t n) override {
        _state.tree_num_evaluated = n;
    }
    inline void addToTreeMemory(size_t n) override{
        _state.tree_memory += n;
        _state.tree_max_memory += n;
    }
    inline void removeFromTreeMemory(size_t n) override{
        _state.tree_memory -= n;
    }
    inline std::pair<size_t, size_t> getTreeMemory() override {
        return std::make_pair(_state.tree_memory, _state.tree_max_memory);
    }
    inline void addToQueueInsertionTime(double t) override {
        _state.queue_insertion_time += t;
    }
    inline void setQueueSize(size_t n) override {
        _state.queue_size = n;
    }
    inline void addToQueueMemory(size_t n) override{
        _state.queue_memory += n;
        _state.queue_max_memory += n;
    }
    inline void removeFromQueueMemory(size_t n) override{
        _state.queue_memory -= n;
    }
    inline std::pair<size_t, size_t> getQueueMemory() override {
        return std::make_pair(_state.queue_memory, _state.queue_max_memory);
    }
    inline void setNRules(size_t nrules) override {
        // the first rule is the default rule
        _state.nrules = nrules - 1;
    }
    inline void setC(double c) override {
        _state.c = c;
    }
    inline void initPrefixVec() override {
        _state.prefix_lens = (size_t*) calloc(_state.nrules, sizeof(size_t));
    }
    inline void incPrefixLen(size_t n) override {
        ++_state.prefix_lens[n];
        if (_state.prefix_lens[n] == 1)
            updateQueueMinLen();
    }
    inline void decPrefixLen(size_t n) override {
        --_state.prefix_lens[n];
        if (_state.prefix_lens[n] == 0)
            updateQueueMinLen();
    }
    /*
     * Returns the size of the logical queue.
     */
    inline size_t sumPrefixLens() override {
        size_t tot = 0;
        for(size_t i = 0; i < _state.nrules; ++i) {
            tot += _state.prefix_lens[i];
        }
        return tot;
    }
    inline void updateQueueMinLen() override {
        // Note: min length is logically undefined when queue size is 0
        size_t min_length = 0; 
        for(size_t i = 0; i < _state.nrules; ++i) {
            if (_state.prefix_lens[i] > 0) {
                min_length = i;
                break;
            }
        }
        _state.queue_min_length = min_length;
    }
    inline size_t getQueueMinLen() override {
        return _state.queue_min_length;
    }
    inline void incPmapSize() override {
        ++_state.pmap_size;
    }
    inline void decreasePmapSize(size_t n) override {
        _state.pmap_size -= n;
    }
    inline void incPmapNullNum() override {
        ++_state.pmap_null_num;
    }
    inline void incPmapDiscardNum() override {
        ++_state.pmap_discard_num;
    }
    inline void addToPmapMemory(size_t n) override {
        _state.pmap_memory += n;
        _state.pmap_max_memory += n;
    }
    inline void removeFromPmapMemory(size_t n) {
        _state.pmap_memory -= n;
    }
    inline std::pair<size_t, size_t> getPmapMemory() override {
        return std::make_pair(_state.pmap_memory, _state.pmap_max_memory);
    }
    inline void subtreeSize(mpz_t tot, unsigned int len_prefix, double lower_bound) override {
        // Theorem 4 (fine-grain upper bound on number of remaining prefix evaluations)
        unsigned int f_naive = _state.nrules - len_prefix;
        unsigned int f = (_state.tree_min_objective - lower_bound) / _state.c;
        if (f_naive < f)
            f = f_naive;
        mpz_set_ui(tot, _state.nrules - len_prefix);
        for (unsigned int k = (_state.nrules - len_prefix - 1); 
                k >= (_state.nrules - len_prefix - f + 1); k--) {
            mpz_addmul_ui(tot, tot, k);
        }
    }
    inline void approxRemainingSize(mpz_t tot, unsigned int len_prefix) override {
        // Proposition 3 (coarse-grain upper bound on number of remaining prefix evaluations)
        size_t K = (size_t) (_state.tree_min_objective / _state.c);
        if (K > _state.nrules)
            K = _state.nrules;

        // sum_{j=0}^M Q_j sum_{k=1}^{K-j} (M - j)! / (M - j - k)!
        mpz_set_ui(tot, _state.nrules - len_prefix);
        for(size_t k = (_state.nrules - len_prefix - 1); k >= (_state.nrules - len_prefix - K + 1); --k)
            mpz_addmul_ui(tot, tot, k);

        // multiply by Qj
        mpz_mul_ui(tot, tot, _state.prefix_lens[len_prefix]);
    }
    inline void addQueueElement(unsigned int len_prefix, double lower_bound, bool approx) override {
        mpz_t tot;
        mpz_init(tot);
        if (approx)
            approxRemainingSize(tot, len_prefix);
        else
            subtreeSize(tot, len_prefix, lower_bound);
        mpz_add(_state.remaining_space_size, _state.remaining_space_size, tot);
        mpz_clear(tot);
    }
    inline void removeQueueElement(unsigned int len_prefix, double lower_bound, bool approx) override {
        mpz_t tot;
        mpz_init(tot);
        if (approx)
            approxRemainingSize(tot, len_prefix);
        else
            subtreeSize(tot, len_prefix, lower_bound);
        mpz_sub(_state.remaining_space_size, _state.remaining_space_size, tot);
        mpz_clear(tot);
    }
    inline void initRemainingSpaceSize() override {
        // Proposition 2 (upper bound on total number of prefix evaluations)
        size_t naive_max_length = 0.5 / _state.c;
        if (naive_max_length < _state.nrules)
            mpz_fac_ui(_state.remaining_space_size, naive_max_length);
        else
            mpz_fac_ui(_state.remaining_space_size, _state.nrules);
    }
    inline void clearRemainingSpaceSize() override {
        mpz_set_ui(_state.remaining_space_size, 0);
    }
    inline size_t getLogRemainingSpaceSize() override {
        // This is approximate.
        return mpz_sizeinbase(_state.remaining_space_size, 10); 
    }
    /*
     * Initializes the logger by setting all fields to 0.
     * This allwos us to write a log record immediately.
    inline void initializeState() override {
        _state.total_time = 0.;
        _state.evaluate_children_time = 0.;
        _state.evaluate_children_num = 0;
        _state.node_select_time = 0.;
        _state.node_select_num = 0;
        _state.rule_evaluation_time = 0.;
        _state.rule_evaluation_num = 0;
        _state.lower_bound_time = 0.;
        _state.lower_bound_num = 0;
        _state.objective_time = 0.;
        _state.objective_num = 0;
        _state.tree_insertion_time = 0.;
        _state.tree_insertion_num = 0;
        _state.permutation_map_insertion_time = 0.;
        _state.permutation_map_insertion_num = 0;
        _state.current_lower_bound = 0.;
        _state.tree_min_objective = 1.;
        _state.tree_prefix_length = 0;
        _state.tree_num_nodes = 0;
        _state.tree_num_evaluated = 0;
        _state.queue_insertion_time = 0;
        _state.queue_size = 0;
        _state.queue_min_length = 0;
        _state.pmap_size = 0;
        _state.pmap_null_num = 0;
        _state.pmap_discard_num = 0;
        mpz_init(_state.remaining_space_size);
        initRemainingSpaceSize();
    }
     */
};

extern Logger logger;

inline double timestamp() {
    struct timeval now;
    gettimeofday(&now, 0);
    return now.tv_sec + now.tv_usec * 0.000001;
}

inline double time_diff(double t0) {
    return timestamp() - t0;
}

#include "alloc.hh"
/* 
 * Prints the final rulelist that CORELS returns.
 * rulelist -- rule ids of optimal rulelist
 * preds -- corresponding predictions of rules (+ default prediction)
 */
void print_final_rulelist(const std::vector<unsigned short, cache_alloc<unsigned short> >& rulelist,
                          const std::vector<bool, cache_alloc<bool> >& preds,
                          const bool latex_out,
                          const rule_t rules[],
                          const rule_t labels[],
                          char fname[]);

void print_machine_info();
