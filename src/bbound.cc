#include "bbound.hh"


CacheNode* stochastic_select(CacheTree* tree, VECTOR not_captured) {
    std::map<size_t, CacheNode*>::iterator iter;
    CacheNode* node = tree->root();
    rule_copy(not_captured, tree->rule(node->id()).truthtable, tree->nsamples());
    int cnt;
    while (node->done()) {
        if ((node->lower_bound() + tree->c()) >= tree->min_objective()) {
            if (node->depth() > 0) {
                CacheNode* parent = node->parent();
                parent->delete_child(node->id());
                tree->delete_subtree(node);
            }
            return 0;
        }
        if (node->num_children() == 0) {
            tree->prune_up(node);
            return 0;
        }
        iter = node->random_child();
        node = iter->second;
        rule_vandnot(not_captured, not_captured, tree->rule(iter->first).truthtable, tree->nsamples(), &cnt);
    }
    return node;
}

void bbound_stochastic(CacheTree* tree, size_t max_num_nodes) {
    CacheNode* node;
    VECTOR not_captured;
    size_t num_iter = 0;
    rule_vinit(tree->nsamples(), &not_captured);
    tree->insert_root();
    while ((tree->num_nodes() < max_num_nodes) and (tree->num_nodes() > 0)) {
        node = stochastic_select(tree, not_captured);
        if (node)
            tree->evaluate_children(node, not_captured);
        ++num_iter;
        if ((num_iter % 10000) == 0)
            printf("num_iter: %zu, num_nodes: %zu\n", num_iter, tree->num_nodes());
    }
    rule_vfree(&not_captured);
}
