#include <sys/time.h>
#include <string.h>

#ifndef _UTILS_H_
#define _UTILS_H_

class Logger {
  public:
    Logger(int verbosity): _v(verbosity) {

    }

    inline void clearState();
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
  private:
    struct State {
        double total_time;
        double evaluate_children_time;
        int evaluate_children_num;
        double node_select_time;
        int node_select_num;
        double rule_evaluation_time;
        int rule_evaluation_num;
        double lower_bound_time;
        int lower_bound_num;
        double objective_time;
        int objective_num;
        double tree_insertion_time;
        int tree_insertion_num;
        double permutation_map_insertion_time;
        int permutation_map_insertion_num;
    };
    State _state;
    int _v; // verbosity
    // TODO: store file handler here
};
inline void Logger::clearState() {
    memset(&_state, 0, sizeof(_state));
}

inline double timestamp() {
    struct timeval now;
    gettimeofday(&now, 0);
    return now.tv_sec + now.tv_usec * 0.000001;
}

inline double time_diff(double t0) {
    return timestamp() - t0;
}

#endif
