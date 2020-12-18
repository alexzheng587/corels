#pragma once

#include "rule.h"
#include <vector>
#include <algorithm>

struct auc_calculation_t {
    int npos;
    int nneg;
    double positive_fraction;
    // auc_calculation_t(int npos_, int nneg_, double positive_fraction_) 
    // : npos(npos_), nneg(nneg_), positive_fraction(positive_fraction_) {}
};

class ProblemState {
public:
    explicit ProblemState(size_t nsamples) {
        rule_vinit(nsamples, &captured);
        rule_vinit(nsamples, &captured_zeros);
        rule_vinit(nsamples, &not_captured);
        rule_vinit(nsamples, &not_captured_equivalent);
        rule_vinit(nsamples, &not_captured_equivalent_fn);
        rule_vinit(nsamples, &not_captured_equivalent_fp);
    }

    ~ProblemState() {
        rule_vfree(&captured);
        rule_vfree(&captured_zeros);
        rule_vfree(&not_captured);
        rule_vfree(&not_captured_equivalent);
        rule_vfree(&not_captured_equivalent_fn);
        rule_vfree(&not_captured_equivalent_fp);
    }

    VECTOR captured{};
    int n_captured{};
    int n_captured_correct{};

    VECTOR captured_zeros{};
    int n_zeros_captured{};

    int n_ones_captured{};

    int n_equiv_captured{};

    VECTOR not_captured{};
    int n_not_captured{};

    int n_parent_not_captured{};
    int n_parent_zeros_not_captured{};

    double rule_loss{};
    bool rule_prediction{};

    double default_loss{};
    bool default_prediction{};

    double equivalent_points_loss{};

    double rule_lower_bound{};
    double parent_lower_bound{};
    double objective{};
    double regularization_loss{};

    VECTOR not_captured_equivalent{};
    int n_not_captured_equivalent{};

    VECTOR not_captured_equivalent_fn{};
    int n_not_captured_equivalent_fn{};

    VECTOR not_captured_equivalent_fp{};

    int false_positives{};
    int false_negatives{};

    minority_class_t* class_pointer;
    std::vector<int> class_idx_list;

    void insert_prefix(auc_calculation_t new_node) {
        prefix_order_.push_back(new_node);
    }

    static bool compare_positive_fraction(auc_calculation_t const &a, auc_calculation_t const &b) {
        return a.positive_fraction > b.positive_fraction;
    }

    void sort_prefix() {
        std::sort(prefix_order_.begin(), prefix_order_.end(), compare_positive_fraction);
    }

    std::vector<auc_calculation_t> get_prefix() {
        return prefix_order_;
    }

    void clear_prefix() {
        prefix_order_.clear();
    }

protected:
    std::vector<auc_calculation_t> prefix_order_;
};