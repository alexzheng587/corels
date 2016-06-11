# bbcache
Branch-and-bound algorithm, with caching, for decision lists.

## Dependencies

### Python dependencies

    Python 2.7x
    numpy
    tabular
    matplotlib
    gmpy2

### Dependencies of gmpy2

These can be installed on Mac OS X with `brew install`.

    mpfr
    libmpc

## data/

### tic-tac-toe training dataset

    tdata_R.out :  set of rules
    tdata_R.labels :  labels

Note that the labeled training data is biased in the sense that 63.8% (408/639)
have label 1.  We thus expect short rule lists to predict 1 by default.

By design, the set of rules has 377 elements and does not include all possible
rules a human might think to write down for this problem.  Specifically, it does
not include two of the rules for 'o' winning: {c1=o,c2=o,c3=o} and
{c7=o,c8=o,c9=o}.

It turns out that the greedy algorithm identifies as the first rule a rule for
'o' winning: {c2=o,c5=o,c8=o}.  As a result, the greedy algorithm continues to
add the other rules for 'o' winning.  Since only 6 out of 8 such rules exist in
the set of rules, the greedy algorithm, with maximum rule list length set to 8,
outputs a list that only achieves an accuracy of 0.912 on the training set.
This is fine for us, since we are only using the greedy algorithm to find a
reasonable warm start rule list.

    greedily constructed prefix = (153, 372, 176, 178, 134, 51, 121, 125)

Note that if we initialize the greedy algorithm with a rule for 'x' winning, it
will find a perfect rule list.  Furthermore, the rules in this rule list
perfectly partition the dataset.

### adult dataset

    adult_R.out :  set of rules
    adult_R.labels :  labels

    nrules : 284
    ndata : 30081

    sum(ones) : 22645
    bias : 0.7528

    total pairs : 40186
    commuting pairs : 7235

## code/

### branch_bound.py

This module contains functions and data structures used by variants of the
branch-and-bound algorithm, including an object for a cache element.

### serial.py

This module contains a serial implementation of the branch-and-bound algorithm,
with a cache to support incremental computation.  Prefixes are added to the
queue greedily, which causes the queue to grow exponentially fast.  It's here
because it's a bit easier to understand.

### serial_lazy.py

This module contains a serial implementation of the branch-and-bound algorithm,
with a cache to support incremental computation and a queue that grows lazily.
The queue's length is bounded by the total number of rules.

How many prefixes get stored in the cache?  We empirically measure this, for
different (simulated) warm starts.  Note that to simulate a warm start with a
particular accuracy, we don't actually need to generate a prefix with that
accuracy -- we simply set max_accuracy to the desired level.  For the
tic-tac-toe dataset, we know that a perfect rule list can be generated from the
given rules.

### serial_gc.py

This module contains a serial implementation of the branch-and-bound algorithm,
with a cache to support incremental computation, a queue that grows lazily, and
garbage collection.  The queue's length is bounded by the total number of rules.
This will replace serial_lazy.py, as includes introduces garbage collection in a
modular and optional fashion.

Each round, we track groups of prefixes that are equivalent up to permutation.
Since the prefixes in such a group capture the same data, we only keep one that
has the highest accuracy within the group.

    cache_size[i] + captured_zero[i] + dead_prefix[i] + inferior[i]
    = (nrules - i + 1) * (cache_size[i-1] - dead_prefix_start[i] - stunted_prefix[i])

Update: additional symmetry-based pruning based on (all sets of) rules that commute.

### serial_priority.py

This module contains a serial implementation of the branch-and-bound algorithm,
with:

* A cache that supports incremental computation.  The cache data structure is a
hash map (Python `dict`) with keys encoding prefixes and also has a tree
structure reflecting relationships among prefixes.

* A priority queue to manage different scheduling policies (implemented using
Python's `heapq` module).

* Symmetry-aware garbage collection for sets of prefixes that are equivalent up
to permutation -- we only keep the best.

* Symmetry-aware pruning for equivalence classes of prefixes that contain
(possibly multiple) adjacent pairs of commuting rules -- we only evaluate one.

#### Cache

The cache data structure is a hash map that also encodes a tree.
The keys of the hash map encode prefixes; retrieving a cache entry for a given
prefix is thus an O(1) operation.
The tree structure reflects relationships among prefixes, and these
relationships enable incremental computation in our branch-and-bound algorithm.

Each node of the tree is indexed by a prefix.
The root node is the empty prefix and nodes at depth `d` (below) the root
correspond to prefixes of length `d`.
The children of a node indexed by `prefix` are prefixes starting with `prefix`
and longer by one rule (not in `prefix`).
The parent of a node indexed by `prefix` of length `k > 0` is therefore indexed
by `trim(prefix)`, where this operation returns a prefix whose rules are given
by the first `k - 1` rules of `prefix`.
It is thus efficient to find the parent and other ancestors of a given prefix,
i.e., to traverse the tree up towards the root.
To support efficient traversals in the direction of descendents, a node must
also keep track of its children; our implementation is detailed below, in our
description of cache entry attributes.

The cache's tree structure also encodes the state of computation and contains
two kinds of nodes, `used` and `unused`.
When a cache entry is first inserted, it is `unused` in the sense that it has
not yet been (completely) used in incremental computations related to its children.
If it is ever used for such computations, we call it `used`.
Our `prune_up` routine enforces that that all `used` cache entries are interior
nodes in the tree and all `unused` cache entries are leaves.

An `unused` cache entry has the following attributes:

* **prefix** : n-tuple of integers corresponding to rules (redundant w.r.t. key)

* **prediction** : binary n-tuple of corresponding predictions

* **default_rule** : binary value encoding the default rule

* **accuracy** : real value in [0, 1], accuracy (not necessary)

* **upper_bound** : upper bound on the accuracy for rule lists starting with prefix (not necessary)

* **objective** : objective value of the prefix (not necessary)

* **lower_bound** : lower bound on the objective value for rule lists starting with prefix

* **num_captured** : number of data captured by the prefix

* **num_captured_correct** : number of data captured by the prefix and predicted correctly

* **not_captured** : binary vector encoding which data are not captured by the prefix (represented as an integer using `mpz`)

* **curiosity** : curiosity (not necessary for cache but can be used by priority queue)

* **children** : set of integers encoding child nodes

* **num_children** : number of children

* **reject_list** : tuple of integers encoding rejected rules that should never be appended to `prefix` (if a rule captures insufficient data when appended to a prefix, then it will be insufficient for any rule list that starts with that prefix)

If an `unused` cache entry becomes `used`, then we delete these attributes used
in the computation because they will not be used again:

    prefix, prediction, default_rule, accuracy, objective,
    num_captured, num_captured_correct, not_captured, curiosity

A `used` cache entry thus retains four attributes:

* **lower_bound** : the lower bounds of interior nodes are used to prune the cache when the minimum observed objective decreases (via the `garbage_collect` routine)

* **children** : this attribute encodes the tree structure

* **num_children** : if `num_children` reaches zero, the `used` cache entry is a `dead end` and we delete it (via the `prune_up` routine, which also removes `dead end` ancestors)

* **reject_list** : when an `unused` cache entry is retrieved, this attribute is lazily initialized by copying its parent's `reject_list`

##### Elementary cache operations

* **insert(prefix, cache_entry)** : insert an `unused` cache entry for `prefix`
(and update its parent)

* **delete(prefix)** : delete `prefix` and all of its descendents (and update
its parent, calling `prune_up` if relevant)

##### Cache garbage collection routines

* **insert(prefix, cache_entry)** : symmetry-aware garbage collection happens on cache entry insertion; we discuss this below

* **prune_up(prefix)** : if `prefix` corresponds to a `dead end`, i.e., a `used`
cache entry with no children, then remove cache entries for `prefix` as well as
any `dead end` ancestors

* **garbage_collect(min_objective)** : delete all prefixes with `lower_bound > min_objective)`; called when `min_objective` decreases

### Priority queue

#### Prioritization metrics

Each is a function that maps a prefix to a value (bounded below by zero,
corresponding to the highest priority).

**Prefix length** : Implements breadth-first search

**Objective lower bound** : The lower bound on the objective is the sum of the
lower bound on the number of mistakes and the regularization term.

    priority = (# captured and incorrect) / (# data) + c * (prefix length)

**Curiosity** : Expected objective if all data are captured.

    priority = (prefix misclassification) + c * (expected prefix length)

Where

    prefix misclassification = (# captured and incorrect) / (# captured)

and

    expected prefix length = (prefix length) * (# data) / (# captured)

**Objective** : The objective is the sum of the fraction of mistakes and a
regularization term.
 
    priority = (# incorrect) / (# data) + c * (prefix length)

#### Symmetry-aware pruning

Two rules A and B commute if they capture non-intersecting subsets of data.
If rules A and B commute, then a rule list where A and B are adjacent is
equivalent to another rule list where A and B swap positions.
More generally, a rule list containing possibly multiple, possibly overlapping,
pairs of commuting rules is equivalent to any other rule list that can be
generated by swapping one or more such pairs of rules.
We avoid evaluating multiple such equivalent rule lists by eliminating all but one.

#### Symmetry-aware garbage collection

If two prefixes P and Q are composed of the same rules and equivalent up to a
permutation, then they also capture the same data.
The two corresponding rule lists need not yield the same objective, since the
objective depends on rule order.
Obtain a rule list P' by appending P with some ordered list of unique rules not
contained in P, and Q' by appending Q with the same ordered list.
The performance of P' compared to P will be the same as that of Q' compared to Q, e.g.,

    objective(P') - objective(P) = objective(Q') - objective (Q).

Thus, we can delete whichever of P or Q performs worse than the other.
We call this symmetry-aware garbage collection.
Since a prefix of length `k` belongs to an equivalence class of `k!` prefixes
equivalent up to permutation, this garbage collection dramatically prunes the
search space.

(More generally, when two prefixes capture the same data, we only need to keep
the better of the two; we have not investigated the utility of this idea.)

We perform symmetry-aware garbage collection on cache insertion, which enforces
that the cache never contains multiple prefixes equivalent up to permutation.
We support this via a hash map data structure (Python `dict`) that we call the
**inverse canonical map** (ICM).
A set of all prefixes equivalent up to permutation is represented by any member
of the set; we choose the prefix with rules ordered from least to greatest to
represent this set, and call this prefix **canonical**.
The ICM maps canonical prefixes to cached prefixes (and for convenience,
their objectives) and contains one entry for each cache entry, and lets us
quickly check whether a prefix is equivalent to a cached prefix up to
permutation, and if so, determine which is better.
We couple ICM and cache updates:  cache insertions and deletions trigger
corresponding ICM operations.

## tic-tac-toe results

### tdata, breadth-first, no regularization (c = 0.)

* Aggressive "optimistic (lying) warm start" (initialize min_objective = 0.001)
* Finds a minimum length (8 rules) perfect prefix (< 180 sec)

### tdata, curiosity, no regularization (c = 0.)

* Cold start, quickly (< 2 sec) finds a perfect prefix that is very long (68 rules)

### tdata, curiosity, regularization (c = 0.001)

* Cold start, quickly (< 2 sec) finds a perfect prefix of length 10
* Then certifies that the best (perfect) prefix has length 8 (~250 sec)
* Final cache contains 1776 entries with lower bound < 0.008

## adult results

Subsampling 10% of dataset and `min_captured_correct = max(c, 0.003)` unless
otherwise noted.

### adult, curiosity, no regularization (c = 0.)

* Certifies there are no prefixes with objective < 0.005
* Up to symmetries, ~170,000 prefixes have lower bound < 0.005 (length <= 11)

* Initialize min_objective = 0.01
* Inconclusive, > 3,000,000 prefixes have lower bound < 0.01

### adult, curiosity, regularization (c = 0.003)

* Initialize min_objective = 0.01
* Quickly certifies (< 1 sec) that there are no prefixes with objective < 0.01
* Up to symmetries, 471 prefixes have lower bound < 0.01

* Inconclusive for objective < 0.05
* Up to symmetries, > 3,000,000 prefixes have lower bound < 0.05

### adult, curiosity, regularization (c = 0.01)

* Certifies (< 380 sec) that there are no prefixes with objective < 0.05
* Up to symmetries, ~150,000 prefixes have lower bound < 0.05 (length <= 4)

* Certifies (< 18,000 sec) that there are no prefixes with objective < 0.06
* Up to symmetries, ~1,000,000 prefixes have lower bound < 0.06 (length <= 5)

* Inconclusive for objective < 0.08
* Up to symmetries, > 3,000,000 prefixes have lower bound < 0.08

### adult, curiosity, regularization (c = 0.02)

* Certifies (< 95 sec) there are no prefixes with objective < 0.08
* Certifies (< 240 sec) there are no prefixes with objective < 0.09
<!--* Certifies (< 1400 sec) there are no prefixes with objective < 0.10
* Certifies (< 5500 sec) there are no prefixes with objective < 0.11-->

### adult, curiosity, aggressive regularization (min_captured_correct = c)

* (c = 0.1) Certifies (< 1.5 sec) the best prefix is (51,) on school laptop
* (c = 0.09) Certifies (< 5 sec) the best prefix is (118,)
* (c = 0.08) Certifies (< 10 sec) the best prefix is (118,)
* (c = 0.07) Certifies (< 23 sec) the best prefix is (118,)
* (c = 0.06) Certifies (< 280 sec) the best prefix is (69,)
* (c = 0.05) Certifies (< 1760 sec) the best prefix is (69,)
* (c = 0.04) (> 18000 sec)

### adult, objective, no regularization (c = 0.)

* Full dataset
* Quickly (< 0.1 sec) finds a good prefix (accuracy = 0.83328, length = 10)
* (43, 69, 122, 121, 77, 0, 62, 1, 46, 24)
* Slight improvement after ~200 sec (accuracy = 0.83338, length = 15)
* (43, 69, 122, 121, 77, 0, 62, 1, 46, 24, 22, 237, 39, 145, 56)
* Slight improvement at ~630 sec (accuracy = 0.83342, length = 15)
* (43, 69, 122, 121, 77, 0, 62, 1, 46, 24, 113, 129, 233, 145, 242)
* Slight improvement at ~1,800 sec (accuracy = 0.83345, length = 15)
* (43, 69, 122, 121, 77, 0, 62, 1, 46, 28, 113, 237, 144, 130, 56)
* Cache remains very small (< 2,500 entries)

## implemented but not explained

If `c > 0`, don't add prefix to priority queue or cache if
`c * (len(prefix) + 1) >= min_objective`.

## todo

Insertions and deletions into pdict (and priority_queue?) should happen in
parallel with analogous cache operations.

May want to completely remove dependence of incremental on cache.

Should we skip symmetry-based garbage collection (via `pdict`) when
`len(prefix) == max_prefix_len_check`?

Some of the metrics (`commutes` and `captured_zero`) currently include other
recently added savings -- separate these out properly.

Could track children explicitly as a cache entry attribute.  This would
facilitate proper garbage collection. (Currently, finding ancestors is easy, but
finding children is dumb.)

Not properly updating `num_children` in cache entries or `metrics.cache_size`.

Should we maintain a (lazily / partially materialized?) priority queue ordered
by lower bound (to facilitate garbage collection when the min_objective
decreases)?

More efficient way to encode lists of integers (currently tuples of integers).

Combine objective-based priority with restriction (prevent excessive exploration
of a single subtree).

Restrict search to sub-tree.  (Could be implemented by thresholding `rules`.
May want this as a subroutine to enable smaller cache entries, and perhaps as
part of a parallel scheme.)

Implement back-off. (Something like put the cache entries back on the queue and
increase the lying/optimistic `min_objective`.)

An actual pruning routine.

Rename `max_accuracy` to `best_prefix_accuracy`.

Remove prefix from cache entry (it's already the corresponding key).

Fix the greedy algorithm to stop early based on evaluating the default rule.

Think about useful heuristics to cut down on the size of the search space.

Cynthia's optimization.

We might want to explore notions of approximate equivalence.
