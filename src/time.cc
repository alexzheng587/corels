#include "time.hh"

void clear_time(struct time* t) {
    t->total_time = 0;
    t->evaluate_children_time = 0;
    t->evaluate_children_num = 0;
    t->node_select_time = 0;
    t->node_select_num = 0;
    t->rule_evaluation_time = 0;
    t->rule_evaluation_num = 0;
    t->lower_bound_time = 0;
    t->lower_bound_num = 0;
    t->objective_time = 0;
    t->objective_num = 0;
    t->tree_insertion_time = 0;
    t->tree_insertion_num = 0;
}
