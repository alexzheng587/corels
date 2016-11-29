#include "utils.hh"

#include <stdio.h>
#include <assert.h>


void Logger::setLogFileName(char *fname) {
    if (_v < 1) return;

    printf("writing logs to: %s\n\n", fname);
    _f.open(fname, ios::out | ios::trunc);

    _f << "total_time,evaluate_children_time,node_select_time,"
       << "rule_evaluation_time,lower_bound_time,lower_bound_num,"
       << "objective_time,tree_insertion_time,tree_insertion_num,"
       << "permutation_map_insertion_time,permutation_map_insertion_num,"
       << "tree_min_objective,tree_num_nodes,tree_num_evaluated,"
       << "queue_size, prefix_lengths" << endl;
}

void Logger::dumpState() {
    if (_v < 1) return;

    setTotalTime(time_diff(_state.initial_time));   // update timestamp here

    _f << _state.total_time << ","
       << _state.evaluate_children_time << ","
       << _state.node_select_time << ","
       << _state.rule_evaluation_time << ","
       << _state.lower_bound_time << ","
       << _state.lower_bound_num << ","
       << _state.objective_time << ","
       << _state.tree_insertion_time << ","
       << _state.tree_insertion_num << ","
       << _state.permutation_map_insertion_time << ","
       << _state.permutation_map_insertion_num << ","
       << _state.tree_min_objective << ","
       << _state.tree_num_nodes << ","
       << _state.tree_num_evaluated << ","
       << _state.queue_size << ","
       << dumpPrefixLens().c_str() << endl;
}

std::string Logger::dumpPrefixLens() {
    std::string s = "";
    for(size_t i = 0; i < _state.nrules; ++i) {
        if (_state.prefix_lens[i] > 0) {
            s += std::to_string(i); 
            s += ":";
            s += std::to_string(_state.prefix_lens[i]);
            s += ";";
        }
    }
    return s;
}

void print_final_rulelist(const std::vector<size_t>& rulelist,
                          const std::vector<bool>& preds,
                          const rule_t rules[],
                          const rule_t labels[]) {
    assert(rulelist.size() == preds.size() - 1);

    printf("\nOPTIMAL RULE LIST\n");
    printf("if (%s) then (%s)\n", rules[rulelist[0]].features,
           labels[preds[0]].features);
    for (size_t i = 1; i < rulelist.size(); ++i) {
        printf("else if (%s) then (%s)\n", rules[rulelist[i]].features,
               labels[preds[i]].features);
    }
    printf("else (%s)\n\n", labels[preds.back()].features);
}
