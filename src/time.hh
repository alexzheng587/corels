#include <sys/time.h>

struct time {
    double total_time;
    double evaluate_children_time;
    int evaluate_children_num;
    double stochastic_select_time;
    int stochastic_select_num;
    double rule_evaluation_time;
    int rule_evaluation_num;
    double lower_bound_time;
    int lower_bound_num;
    double objective_time;
    int objective_num;
    double tree_insertion_time;
    int tree_insertion_num;
};

void clear_time(struct time* t);

inline double timestamp() {
    struct timeval now;
    gettimeofday(&now, 0);
    return now.tv_sec + now.tv_usec * 0.000001;
}
