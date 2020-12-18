#include "loss.hh"

void Loss::initialize_subproblem(ProblemState* state, CacheTree* tree, Node* parent, VECTOR parent_not_captured) {
    double c = tree->c(); // c is the regularization parameter
    double parent_lower_bound = parent->lower_bound();
    double parent_equivalent_minority = parent->equivalent_minority();
    // subtract off parent equivalent points bound because we want to use pure lower bound from parent
    state->parent_lower_bound = parent_lower_bound - parent_equivalent_minority + c;

    // Store how many zeros we haven't yet captured, to be used within the loop
    // use not_captured as a scratch space vector (will be overwritten)
    state->n_parent_not_captured = tree->nsamples() - parent->num_captured();
    rule_vand(state->not_captured, parent_not_captured, tree->label(0).truthtable, tree->nsamples(),
              &state->n_parent_zeros_not_captured);
    
}

void Loss::make_prediction(ProblemState* state, int npos, int nneg, Node* parent) {
    int n_ones_captured = state->n_captured - state->n_zeros_captured;
    double zero_loss = compute_loss(0, n_ones_captured, npos, nneg);
    double one_loss = compute_loss(state->n_zeros_captured, 0, npos, nneg);

    if (zero_loss < one_loss) {
        state->rule_loss = zero_loss;
        state->rule_prediction = false;
        state->n_captured_correct = state->n_zeros_captured;
    } else {
        state->rule_loss = one_loss;
        state->rule_prediction = true;
        state->n_captured_correct = n_ones_captured;
    }
}

void Loss::compute_lower_bound(ProblemState* state, int npos, int nneg) {
    state->rule_lower_bound = state->parent_lower_bound + state->rule_loss;
}

void Loss::make_default_prediction(ProblemState* state, int npos, int nneg) {
    int default_captured_zeros = state->n_parent_zeros_not_captured - state->n_zeros_captured;
    int default_captured_ones = state->n_not_captured - default_captured_zeros;

    double zero_loss = compute_loss(0, default_captured_ones, npos, nneg);
    double one_loss = compute_loss(default_captured_zeros, 0, npos, nneg);

    if (zero_loss < one_loss) {
        state->default_loss = zero_loss;
        state->default_prediction = false;
    } else {
        state->default_loss = one_loss;
        state->default_prediction = true;
    }
}

void Loss::compute_objective(ProblemState* state) {
    state->objective = state->rule_lower_bound + state->default_loss;
}

bool Loss::min_objective_bounded(ProblemState* state, CacheTree* tree) {
    return state->rule_lower_bound >= tree->min_objective();
}

bool Loss::equivalent_points_bounded(ProblemState* state, CacheTree* tree, Node* parent) {
    rule_vand(state->not_captured_equivalent, state->not_captured, tree->minority(0).truthtable, tree->nsamples(),
              &state->n_not_captured_equivalent);
    rule_vand(state->not_captured_equivalent_fn, state->not_captured_equivalent, tree->label(1).truthtable,
              tree->nsamples(), &state->n_not_captured_equivalent_fn);
    int num_not_captured_equivalent_fp = state->n_not_captured_equivalent - state->n_not_captured_equivalent_fn;

    state->equivalent_points_loss = compute_loss(num_not_captured_equivalent_fp, state->n_not_captured_equivalent_fn,
                                                 tree->npos(), tree->nneg());
    state->rule_lower_bound += state->equivalent_points_loss;

    return min_objective_bounded(state, tree);
}

bool Loss::incremental_accuracy_bounded(ProblemState* state, Node* parent, double regularization) {
    return (parent->objective() - state->objective) <= regularization;
}

bool Loss::lookahead_bounded(ProblemState* state, CacheTree* tree) {
    return state->rule_lower_bound + tree->c() >= tree->min_objective();
}

// ACCURACY

void Accuracy::initialize_loss(CacheTree* tree) {
    min_support_bound_ = tree->c() * tree->nsamples();
}

inline bool Accuracy::minimum_support_bounded(ProblemState* state, Node* parent, CacheTree* tree) const {
    return state->n_captured < min_support_bound_;
}

inline bool Accuracy::accurate_support_bounded(ProblemState* state) const {
    return state->n_captured_correct < min_support_bound_;
}

inline double Accuracy::compute_loss(int nfp, int nfn, int npos, int nneg) const {
    return static_cast<double>(nfp + nfn) / (npos + nneg);
}

// BALANCED ACCURACY

void BalancedAccuracy::initialize_loss(CacheTree* tree) {
    // For support bounds, we want to check later if
    // (0.5)(rule_zeros / total_zeros + rule_ones / total_ones) < lambda
    // => rule_zeros*total_ones + rule_ones*total_zeros < 2 * total_zeros * total_ones * lambda
    min_support_bound_ = 2.0 * tree->nneg() * tree->npos() * tree->c();

    // Duplicate the information in the loss object to avoid changing the API of the support bound checks
    zeros_support_weight_ = tree->npos();
    ones_support_weight_ = tree->nneg();

    // When we don't know which points captured are zeros or ones, multiply by the max of total_ones and total_zeros
    // Cache in object rather than computing each iteration
    largest_support_weight = (tree->nneg() > tree->npos()) ? tree->nneg() : tree->npos();
}

inline double BalancedAccuracy::compute_loss(int nfp, int nfn, int npos, int nneg) const {
    return 0.5 * (static_cast<double>(nfn) / npos + static_cast<double>(nfp) / nneg);
}

inline bool BalancedAccuracy::minimum_support_bounded(ProblemState* state, Node* parent, CacheTree* tree) const {
    return (largest_support_weight * state->n_captured) < min_support_bound_;
}

inline bool
BalancedAccuracy::accurate_support_bounded(ProblemState* state) const {
    return (state->rule_prediction)
           ? (state->n_captured_correct * ones_support_weight_) < min_support_bound_
           : (state->n_captured_correct * zeros_support_weight_) < min_support_bound_;
}

// WEIGHTED ACCURACY

void WeightedAccuracy::initialize_loss(CacheTree* tree) {
    double denominator = weight_ * tree->npos() + tree->nneg();
    min_support_bound_ = tree->c() * denominator;
}

inline double WeightedAccuracy::compute_loss(int nfp, int nfn, int npos, int nneg) const {
    return (nfp + weight_ * nfn) / (weight_ * npos + nneg);
}

inline bool WeightedAccuracy::minimum_support_bounded(ProblemState* state, Node* parent, CacheTree* tree) const {
    return (largest_support_weight_ * state->n_captured) < min_support_bound_;
}

inline bool
WeightedAccuracy::accurate_support_bounded(ProblemState* state) const {
    return (state->rule_prediction)
           ? weight_ * state->n_captured_correct < min_support_bound_
           : state->n_captured_correct < min_support_bound_;
}

// FSCORE

inline double FScore::compute_loss(int nfp, int nfn, int npos, int nneg) const {
    // Precision = true_positives / (true_positives + nfp)
    // Recall = true_positives / (true_positives + nfn)
    // F1_score = 2 * (precision * recall) / (precision + recall)

    // with a general recall weight B, where recall is B times as important as precision:
    // FB_score = (1 + B^2) * (precision * recall) / (B^2 * precision + recall)
    // loss = 1 - FB_score

    // TODO: look into Pauline's formula in gcorels paper, optimize this loss function
    // other formula: (nfp + nfn) / (2*npos + nfp - nfn)
    // tp / (tp + 0.5(fp + fn))

    double tp = npos - nfn; // true positives
    if (tp == 0) { // This will lead to 0 / 0
        if (nfp + nfn == 0) { return 0; } // This is the best possible result, define to be no loss
        else { return 1; } // Else, define to be loss of 1
    }

    double precision = tp / (tp + nfp);
    double recall = tp / (tp + nfn);

    double FB_score = (1 + beta_sqr_)*(precision * recall) / (beta_sqr_*precision + recall);
    return 1 - FB_score;
}

void FScore::initialize_loss(CacheTree* tree) {
    min_support_bound_ = tree->c();
}

void FScore::initialize_subproblem(ProblemState* state, CacheTree* tree, Node* parent, VECTOR parent_not_captured) {
    Loss::initialize_subproblem(state, tree, parent, parent_not_captured);

    state->regularization_loss = tree->c() * static_cast<double>(parent->depth() + 1);
}

bool FScore::minimum_support_bounded(ProblemState* state, Node* parent, CacheTree* tree) const {
    return false;
}

bool FScore::accurate_support_bounded(ProblemState* state) const {
    return false;
}

void FScore::make_prediction(ProblemState* state, int npos, int nneg, Node* parent) {
    int n_ones_captured = state->n_captured - state->n_zeros_captured;
    auto* f1Parent = (F1Node*) (parent);
    int parent_fns = f1Parent->get_false_negatives();
    int parent_fps = f1Parent->get_false_positives();

    int potential_fns = parent_fns + n_ones_captured; // The new number of false negatives if we predict 0
    int potential_fps = parent_fps + state->n_zeros_captured; // If we predict 1

    double zero_loss = compute_loss(parent_fps, potential_fns, npos, nneg);
    double one_loss = compute_loss(potential_fps, parent_fns, npos, nneg);

    if (zero_loss < one_loss) {
        state->rule_loss = zero_loss;
        state->rule_prediction = false;
        state->n_captured_correct = state->n_zeros_captured;
        state->false_negatives = potential_fns;
        state->false_positives = parent_fps;
    } else {
        state->rule_loss = one_loss;
        state->rule_prediction = true;
        state->n_captured_correct = n_ones_captured;
        state->false_negatives = parent_fns;
        state->false_positives = potential_fps;
    }
}

void FScore::compute_lower_bound(ProblemState* state, int npos, int nneg) {
    // Because Fscore is non-additive, assume that the rule loss was computed including the parent's points
    // Additionally, the regularization for the list was not implicitly added to the lower bound,
    // as it was in the case of the additive losses

    // The only lower bound I can think of so far is simply computing the loss without the default rule
    state->rule_lower_bound = state->rule_loss + state->regularization_loss;
}

// TODO: is it valid to choose the default prediction conditional on the rule prediction?
//      Should we be checking every permutation of predictions?
//      I.e., is it possible that we chose to predict negative on a rule
//      because it was optimal without considering the default rule,
//      and now using the resulting false positives and false negatives, we choose to predict positive for default,
//      but we could have a smaller objective the other way around?
void FScore::make_default_prediction(ProblemState* state, int npos, int nneg) {
    // Computes entire objective to choose default rule here because F1 is non-additive
    int default_captured_zeros = state->n_parent_zeros_not_captured - state->n_zeros_captured;
    int default_captured_ones = state->n_not_captured - default_captured_zeros;

    int fns = state->false_negatives;
    int fps = state->false_positives;

    double zero_loss = compute_loss(fps, fns + default_captured_ones, npos, nneg);
    double one_loss = compute_loss(fps + default_captured_zeros, fns, npos, nneg);

    if (zero_loss < one_loss) {
        state->default_loss = zero_loss;
        state->default_prediction = false;
    } else {
        state->default_loss = one_loss;
        state->default_prediction = true;
    }
}

void FScore::compute_objective(ProblemState* state) {
    // Because F1 is non-additive, assume that default loss was computed with all points
    // Additionally, the regularization for the list was not implicitly added to the lower bound
    // as it was in the case of the additive losses
    state->objective = state->default_loss + state->regularization_loss;
}

// AUCLoss

void AUCLoss::initialize_subproblem(ProblemState* state, CacheTree* tree, Node* parent, VECTOR parent_not_captured) {
    double c = tree->c(); // c is the regularization parameter
    double parent_lower_bound = parent->lower_bound();
    double parent_equivalent_minority = parent->equivalent_minority();
    // subtract off parent equivalent points bound because we want to use pure lower bound from parent
    state->parent_lower_bound = parent_lower_bound - parent_equivalent_minority + c;

    state->regularization_loss = tree->c() * static_cast<double>(parent->depth() + 1);

    // Store how many zeros we haven't yet captured, to be used within the loop
    // use not_captured as a scratch space vector (will be overwritten)
    state->n_parent_not_captured = tree->nsamples();
    rule_vand(state->not_captured, parent_not_captured, tree->label(0).truthtable, tree->nsamples(),
              &state->n_parent_zeros_not_captured);

    // initialize parent prefix order
    AUCNode* curr = (AUCNode*) parent;
    // set pointer to minority classes
    state->class_pointer = curr->get_minority_pointer();
    state->clear_prefix();

    while (curr != tree->root())
    {
        auc_calculation_t currentNode = {
            curr->get_npos(),
            curr->get_nneg(),
            (double)curr->get_npos() / (curr->get_npos() + curr->get_nneg()),
        };
        state->insert_prefix(currentNode);
        state->n_parent_not_captured -= curr->get_nneg() + curr->get_npos();
        curr = (AUCNode*) curr->parent();
    }
    state->sort_prefix();
}

void AUCLoss::initialize_loss(CacheTree* tree) {
    min_support_bound_ = std::max(tree->c() * tree->npos(), tree->c() * tree->nneg());
}

bool compare_positive_fraction(auc_calculation_t const &a, auc_calculation_t const &b) {
    return a.positive_fraction > b.positive_fraction;
}

inline double AUCLoss::compute_auc_loss(int n_neg_captured, int n_pos_captured, int npos, int nneg, ProblemState* state) const {
    auc_calculation_t newNode = {
        n_pos_captured,
        n_neg_captured,
        (double) n_pos_captured / (n_neg_captured + n_pos_captured),
    };
    int default_captured_zeros = state->n_parent_zeros_not_captured - n_neg_captured;
    int default_captured_ones = state->n_not_captured - default_captured_zeros;
    auc_calculation_t defaultNode = {
        default_captured_ones,
        default_captured_zeros,
        (double) default_captured_ones / (default_captured_zeros + default_captured_ones),
    };
    std::vector<auc_calculation_t> prefix_order = state->get_prefix();
    prefix_order.push_back(newNode);
    prefix_order.push_back(defaultNode);
    std::sort(prefix_order.begin(), prefix_order.end(), compare_positive_fraction);
    double loss = 0.;
    double positive_leaf_sum = 0.;
    for (int i = 0; i < prefix_order.size(); i++) {
        loss += (prefix_order[i].nneg * (prefix_order[i].npos + positive_leaf_sum));
        positive_leaf_sum += 2.0 * prefix_order[i].npos;
    }
    if (1.0 - (1.0/(2*npos*nneg) * loss) < 0) {
        int qwer = 0;
    }
    return 1.0 - (1.0/(2*npos*nneg) * loss);
}


inline double AUCLoss::compute_loss(int n_neg_captured, int n_pos_captured, int npos, int nneg) const {
    // Compute loss for root node, no prefix_order vector required
    // AUC will be 0.5 since there are only 2 points, one at (0,0) and one at (1,1)
    return 0.5;
}


bool AUCLoss::minimum_support_bounded(ProblemState* state, Node* parent, CacheTree* tree) const {
    return state->n_captured < min_support_bound_;
}

bool AUCLoss::accurate_support_bounded(ProblemState* state) const {
    // TODO
    return false;
}

void AUCLoss::make_prediction(ProblemState* state, int npos, int nneg, Node* parent) {
    // Always predic one since AUC doesn't care about predictions
    state->n_ones_captured = state->n_captured - state->n_zeros_captured;

    double one_loss = compute_auc_loss(state->n_zeros_captured, state->n_ones_captured, npos, nneg, state);

    state->rule_loss = one_loss;
    state->default_loss = one_loss;
    state->rule_prediction = true;
    state->default_prediction = true;

}

void AUCLoss::make_default_prediction(ProblemState* state, int npos, int nneg) {
    return;
}

void AUCLoss::compute_objective(ProblemState* state) {
    // Because AUC is non-additive, assume that default loss was computed with all points
    // Additionally, the regularization for the list was not implicitly added to the lower bound
    // as it was in the case of the additive losses
    state->objective = state->default_loss + state->regularization_loss;
}

void AUCLoss::compute_lower_bound(ProblemState* state, int npos, int nneg) {
    // Because AUC is non-additive, rule loss is computed including the parent's positives and negatives
    // Additionally, the regularization for the list was not implicitly added to the lower bound,
    // as it was in the case of the additive losses

    int default_captured_zeros = state->n_parent_zeros_not_captured - state->n_zeros_captured;
    int default_captured_ones = state->n_not_captured - default_captured_zeros;
    auc_calculation_t leafNode = {
        state->n_ones_captured,
        state->n_zeros_captured,
        (double) state->n_ones_captured / state->n_captured,
    };
    std::vector<auc_calculation_t> prefix_order = state->get_prefix();
    prefix_order.push_back(leafNode);
    std::sort(prefix_order.begin(), prefix_order.end(), compare_positive_fraction);
    double loss = 2 * npos * default_captured_zeros;
    double positive_leaf_sum = 0.;
    for (int i = 0; i < prefix_order.size(); i++) {
        loss += (prefix_order[i].nneg * ((2 * default_captured_ones) + prefix_order[i].npos + positive_leaf_sum));
        positive_leaf_sum += 2.0 * prefix_order[i].npos;
    }
    state->rule_lower_bound = (1.0 - (1.0/(2*npos*nneg) * loss)) + state->regularization_loss;
}

bool AUCLoss::equivalent_points_bounded(ProblemState* state, CacheTree* tree, Node* parent) {
    AUCNode* auc_parent = (AUCNode*) parent;
    int npos = 0;
    int nneg = 0;
    int total_pos = 0;
    int total_neg = 0;
    // get current prefix list
    std::vector<auc_calculation_t> prefix_order = state->get_prefix();
    // add in minority classes
    for (int i = 0; i < auc_parent->get_nclasses(); i++)
    {
        npos = auc_parent->minority_classes(i).npos;
        nneg = auc_parent->minority_classes(i).nneg;
        total_pos += npos;
        total_neg += nneg;
        auc_calculation_t newNode = {
            npos,
            nneg,
            (double) npos / (npos + nneg),
        };
        prefix_order.push_back(newNode);
    }
    // split the remaining samples into two leaves, one containing all positives
    // and one containing all negatives
    auc_calculation_t positiveNode = {
        (state->n_parent_not_captured - state->n_parent_zeros_not_captured) - total_pos,
        0,
        (double) 1,
    };

    auc_calculation_t negativeNode = {
        0,
        state->n_parent_zeros_not_captured - total_neg,
        (double) 0,
    };
    prefix_order.push_back(positiveNode);
    prefix_order.push_back(negativeNode);
    std::sort(prefix_order.begin(), prefix_order.end(), compare_positive_fraction);

    double loss = 0.;
    double positive_leaf_sum = 0.;
    for (int i = 0; i < prefix_order.size(); i++) {
        loss += (prefix_order[i].nneg * ((0.5 * prefix_order[i].npos) + positive_leaf_sum));
        positive_leaf_sum += prefix_order[i].npos;
    }
    state->rule_lower_bound = 1.0 - (1.0/(tree->npos()*tree->nneg()) * loss);

    return min_objective_bounded(state, tree);
}

// Function to determine whether the new rule captures any of the equivalent classes
// If so, remove it from the list of minority classes
void AUCLoss::compute_new_classes(ProblemState* state, CacheTree* tree, Node* parent, VECTOR rule) {
    AUCNode* auc_parent = (AUCNode*) parent;
    VECTOR res;
    rule_vinit(tree->nsamples(), &res);
    std::vector<int> v;
    for (int i = 0; i < auc_parent->get_nclasses(); i++)
    {
        rule_vand(res, rule, auc_parent->minority_classes(i).truthtable, tree->nsamples(), &state->n_equiv_captured);
        if (state->n_equiv_captured == 0) {
            v.push_back(i);
        }
    }
    state->class_idx_list = v;
}