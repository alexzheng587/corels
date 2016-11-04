#include "utils.hh"

#include <stdio.h>

void Logger::dumpState() {
    /** TODO: redirect output to file instead of stdout **/
    printf("Total time: %f\n", _state.total_time);
    printf("Evaluate children time: %f\n", _state.evaluate_children_time);
    printf("Node select time: %f\n", _state.node_select_time);
    printf("Rule evaluation time: %f\n", _state.rule_evaluation_time);
    printf("Lower bound time: %f\n", _state.lower_bound_time);
    printf("Number of lower bound evaluations: %d\n", _state.lower_bound_num);
    printf("Objective time: %f\n", _state.objective_time);
    printf("Total tree insertion time: %f\n", _state.tree_insertion_time);
    printf("Number of tree insertions: %i\n", _state.tree_insertion_num);
    printf("Permutation map insertion time: %f\n", _state.permutation_map_insertion_time);
    printf("Number of calls to permutation_insert(): %i\n", _state.permutation_map_insertion_num);
}
