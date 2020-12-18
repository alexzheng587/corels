#pragma once

#include "rule.h"
#include "cache.hh"
#include "problemState.hh"

class Node;

class CacheTree; // forward declares because of circular dependency

enum LossType {
    ACCURACY,
    BALANCED_ACCURACY,
    WEIGHTED_ACCURACY,
    F_SCORE,
    AUC
};

// Abstract class for computing losses and objectives
class Loss {
public:
    Loss() = default;

    inline LossType type() const { return type_; }

    /**
     * Initializes variables that rely on the input data parameters.
     *
     * Computes support bounds and weights based on the specific loss function
     * and stores them to avoid redundant computation.
     * @param tree
     */
    virtual void initialize_loss(CacheTree* tree) = 0;

    /**
     * Computes the loss incurred by the particular loss function being used.
     * @param nfp: number of false positives.
     * @param nfn: number of false negatives.
     * @param npos: overall number of positive labels.
     * @param nneg: overall number of negative labels.
     * @return: the loss incurred by the particular loss function being used.
     */
    virtual double compute_loss(int nfp, int nfn, int npos, int nneg) const = 0;

    /**
     * Checks whether the number captured by a rule has enough support to be worth considering.
     *
     * 'Worth considering' depends on the loss function, but means that even if we assume that all captured points
     *  would produce the maximum possible loss without this rule, it would still not be worth it.
     * @return true if the rule is provably insufficiently supported to ever be in an optimal rule list.
     */
    virtual bool minimum_support_bounded(ProblemState* state, Node* parent, CacheTree* tree) const = 0;

    /**
     * Checks whether a rule has captured enough points correctly to be worth considering.
     *
     * 'Worth considering' depends on the loss function, but means that even if we assume that all points that we
     * captured and classified correctly would produce the maximum possible loss without this rule,
     * it would still not be worth it.
     *
     * Assumes that state->n_captured is set in make_prediction to be the number of points captured correctly.
     * @return true if the rule has provably not captured enough points correctly to be in an optimal rule list.
     */
    virtual bool accurate_support_bounded(ProblemState* state) const = 0;

    /**
     * Populate the state object with the information that can be precomputed for each rule.
     * @modifies state
     * @param tree
     * @param parent
     * @param parent_not_captured
     */
    virtual void initialize_subproblem(ProblemState* state, CacheTree* tree, Node* parent, VECTOR parent_not_captured);

    /**
     * Evaluates the two choices for the rule under evaluation and populates state with the results.
     * @modifies state
     * @param tree
     * @param parent
     */
    virtual void make_prediction(ProblemState* state, int npos, int nneg, Node* parent);

    /**
     * Computes a lower bound for the current list with the rule under evaluation and populates state with the result.
     * @modifies state
     */
    virtual void compute_lower_bound(ProblemState* state, int npos, int nneg);

    /**
     * After the lower bound is computed, check to see if it is better than the current minimum objective.
     * @param state
     * @param tree
     * @return true if the current list is proven to be non-optimal, false otherwise.
     */
    virtual bool min_objective_bounded(ProblemState* state, CacheTree* tree);

    /**
     * Evaluates the two choices for the default rule after adding the current rule and populates state with the results.
     * @modifies state
     * @param tree
     */
    virtual void make_default_prediction(ProblemState* state, int npos, int nneg);

    /**
     * Compute the final objective using the current rule and populates state with the result.
     * @modifies state
     */
    virtual void compute_objective(ProblemState* state);

    /**
     * Checks that the accuracy gained by adding the current rule is larger than the regularization penalty.
     * @param state
     * @param parent
     * @param regularization
     * @return true if adding this rule is provably non-optimal, false otherwise.
     */
    static bool incremental_accuracy_bounded(ProblemState* state, Node* parent, double regularization);

    /**
     * Evaluates the number of points that are guaranteed to be mispredicted by the default rule,
     * computes how much this will add to our loss,
     * adds the amount to the lower bound of our list,
     * and checks to see if we can no longer beat the current minimum objective.
     * @modifies state
     * @param tree
     * @return true if this rule list is verifiably non-optimal, false otherwise.
     */
    virtual bool equivalent_points_bounded(ProblemState* state, CacheTree* tree, Node* parent);

    /**
     * Looks ahead by one regularization penalty to try to see if we should bother
     * adding this list's children to our data structures.
     * @param state
     * @param tree
     * @return true if this list's children are verifiably non-optimal, false otherwise.
     */
    static bool lookahead_bounded(ProblemState* state, CacheTree* tree);

protected:
    // To avoid dividing number captured by support each time, we multiply lambda by N at the start for accuracy.
    // More generally, for additive loss functions, multiply lambda by the denominator of the loss function.
    double min_support_bound_ = 0;

    LossType type_ = ACCURACY;
};

// (Captured correctly) / (#samples)
class Accuracy : public Loss {
public:
    Accuracy() {
        type_ = ACCURACY;
    }

    double compute_loss(int nfp, int nfn, int npos, int nneg) const override;

    bool accurate_support_bounded(ProblemState* state) const override;

    bool minimum_support_bounded(ProblemState* state, Node* parent, CacheTree* tree) const override;

    void initialize_loss(CacheTree* tree) override;
};

// Average of false positive and false negative rate.
// (1/2)(FN / #positive + FP / #negative)
class BalancedAccuracy : public Loss {
public:
    BalancedAccuracy() {
        type_ = BALANCED_ACCURACY;
    }

    void initialize_loss(CacheTree* tree) override;

    double compute_loss(int nfp, int nfn, int npos, int nneg) const override;

    bool minimum_support_bounded(ProblemState* state, Node* parent, CacheTree* tree) const override;

    bool accurate_support_bounded(ProblemState* state) const override;

protected:
    double largest_support_weight = 1.0; // Use with minimum support bound calculation
    double zeros_support_weight_ = 1.0;
    double ones_support_weight_ = 1.0;
};

// Cost-sensitive average of False Positive and False Negative rate.
// Penalizes false negatives differently according to the weight parameter (< 1 is lighter, > 1 is heavier).
// (FP + weight * FN) / (weight * #positives + #negatives)
class WeightedAccuracy : public Loss {
public:
    explicit WeightedAccuracy(double weight) : weight_{weight} {
        largest_support_weight_ = (weight > 1.0) ? weight : 1.0;
        type_ = WEIGHTED_ACCURACY;
    }

    void initialize_loss(CacheTree* tree) override;

    double compute_loss(int nfp, int nfn, int npos, int nneg) const override;

    bool minimum_support_bounded(ProblemState* state, Node* parent, CacheTree* tree) const override;

    bool accurate_support_bounded(ProblemState* state) const override;

protected:
    const double weight_;
    double largest_support_weight_; // Use with minimum support bounds
};

// Harmonic mean of precision and recall.
// (FP + FN) / (2 * #positive + FP - FN)
class FScore : public Loss {
public:
    explicit FScore(double recall_weight=1.0) : recall_weight_(recall_weight) {
        type_ = F_SCORE;
        beta_sqr_ = recall_weight * recall_weight;
    };

    void initialize_loss(CacheTree* tree) override;

    void initialize_subproblem(ProblemState* state, CacheTree* tree, Node* parent, VECTOR parent_not_captured) override;

    bool minimum_support_bounded(ProblemState* state, Node* parent, CacheTree* tree) const override;

    bool accurate_support_bounded(ProblemState* state) const override;

    void compute_lower_bound(ProblemState* state, int npos, int nneg) override;

    void compute_objective(ProblemState* state) override;

    void make_prediction(ProblemState* state, int npos, int nneg, Node* parent) override;

    void make_default_prediction(ProblemState* state, int npos, int nneg) override;

private:
    double compute_loss(int nfp, int nfn, int npos, int nneg) const override;

protected:
    // with a general recall weight B:
    // FB_score = (1 + B^2) * (precision * recall) / (B^2 * precision + recall)
    const double recall_weight_;
    double beta_sqr_;
};

class AUCLoss : public Loss {
public:
    explicit AUCLoss() {
        type_ = AUC;
    };

    double compute_loss(int nfp, int nfn, int npos, int nneg) const override;

    double compute_auc_loss(int nfp, int nfn, int npos, int nneg, ProblemState* state) const;

    void compute_new_classes(ProblemState* state, CacheTree* tree, Node* parent, VECTOR rule);

    void initialize_loss(CacheTree* tree) override;

    void initialize_subproblem(ProblemState* state, CacheTree* tree, Node* parent, VECTOR parent_not_captured) override;

    bool minimum_support_bounded(ProblemState* state, Node* parent, CacheTree* tree) const override;

    bool accurate_support_bounded(ProblemState* state) const override;

    void compute_objective(ProblemState* state) override;

    void compute_lower_bound(ProblemState* state, int npos, int nneg) override;

    void make_prediction(ProblemState* state, int npos, int nneg, Node* parent) override;

    void make_default_prediction(ProblemState* state, int npos, int nneg) override;

    bool equivalent_points_bounded(ProblemState *state, CacheTree *tree, Node* parent) override;

// protected:
//     std::vector<auc_calculation_t> prefix_order_;
};