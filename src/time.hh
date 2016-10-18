#include <sys/time.h>

struct time {
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

void clear_time(struct time* t);

inline double timestamp() {
    struct timeval now;
    gettimeofday(&now, 0);
    return now.tv_sec + now.tv_usec * 0.000001;
}

inline double time_diff(double t0) {
    return timestamp() - t0;
}
