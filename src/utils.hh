#include <sys/time.h>
#include <string.h>
#include <fstream>

#ifndef _UTILS_H_
#define _UTILS_H_

using namespace std;

class Logger {
  public:
    void closeFile() { if (_f.is_open()) _f.close(); }
    ~Logger() { closeFile(); }

    inline void setVerbosity(int verbosity) { _v = verbosity; }

    void setLogFileName(char *fname);
    void dumpState();

    inline void setLowerBoundTime(double t) {
        _state.lower_bound_time = t;
    }
    inline void incLowerBoundNum() {
        ++_state.lower_bound_num;
    }
    inline void addToObjTime(double t) {
        _state.objective_time += t;
    }
    inline void incObjNum() {
        ++_state.objective_num;
    }
    inline void addToTreeInsertionTime(double t) {
        _state.tree_insertion_time += t;
    }
    inline void incTreeInsertionNum() {
        ++_state.tree_insertion_num;
    }
    inline void addToRuleEvalTime(double t) {
        _state.rule_evaluation_time += t;
    }
    inline void incRuleEvalNum() {
        ++_state.rule_evaluation_num;
    }
    inline void addToNodeSelectTime(double t) {
        _state.node_select_time += t;
    }
    inline void incNodeSelectNum() {
        ++_state.node_select_num;
    }
    inline void addToEvalChildrenTime(double t) {
        _state.evaluate_children_time += t;
    }
    inline void incEvalChildrenNum() {
        ++_state.evaluate_children_num;
    }
    inline void setTotalTime(double t) {
        _state.total_time = t;
    }
    inline void addToPermMapInsertionTime(double t) {
        _state.permutation_map_insertion_time += t;
    }
    inline void incPermMapInsertionNum() {
        ++_state.permutation_map_insertion_num;
    }
    inline void setTreeMinObj(double o) {
        _state.tree_min_objective = o;
    }
    inline void setTreeNumNodes(size_t n) {
        _state.tree_num_nodes = n;
    }
    inline void setTreeNumEvaluated(size_t n) {
        _state.tree_num_evaluated = n;
    }
    inline void setQueueSize(size_t n) {
        _state.queue_size = n;
    }

  private:
    struct State {
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
        size_t permutation_map_insertion_num;
        double tree_min_objective;
        size_t tree_num_nodes;
        size_t tree_num_evaluated;
        size_t queue_size;
    };
    State _state;
    int _v; // verbosity
    ofstream _f; // output file
};

inline double timestamp() {
    struct timeval now;
    gettimeofday(&now, 0);
    return now.tv_sec + now.tv_usec * 0.000001;
}

inline double time_diff(double t0) {
    return timestamp() - t0;
}

extern Logger logger;

#endif
