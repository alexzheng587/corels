# bbcache
Branch-and-bound algorithm, with caching, for decision lists.

The main function is `bbound` in [`code/serial_priority.py`](https://github.com/elaine84/bbcache/blob/heapq/code/serial_priority.py)

and the work in the inner loop is done by `incremental` in [`code/branch_bound.py`](https://github.com/elaine84/bbcache/blob/heapq/code/branch_bound.py)

## Paper

    cd paper
    make

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

    nrules : 377
    ndata : 639

    total pairs : 70876
    commuting pairs : 22445

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

#### Overview, sort of

A small amount of work is done once upfront to compute pairs of rules that
**commute globally** due to zero overlap, and to compute pairs of rules where
one **dominates** another, i.e., A dominates B if A captures all the data that
B captures.

Note that the **maximum length prefix** we have to check is

    M = (best observed objective) / c < 1 / c.

Furthermore, if we ever encounter a perfect prefix of length L, then we can set

    M = L - 1.

(If the policy is breadth first, we can just stop.  Otherwise, the algorithm
stops when the queue is empty.)

Suppose we take prefix P of length K off the queue.  That means we've already
evaluated P and placed it in the cache.

First we check whether P is still in the cache (it could have been garbage
collected, e.g., because we found a new objective smaller than its lower bound,
or because we found a permutation with lower objective).  You could say that we
garbage collect the queue **lazily**.  If P is not in the cache, we stop (and
continue to the next prefix in the queue).

Next we construct the list of new rules to consider.  Naively, this would be all
rules not in the prefix.  Let R be the last rule in P.  We eliminate the
following rules from consideration:

* rules that have zero overlap with R and have a larger index

* rules that are dominated by (any rule in) P

* rules that we already know would be rejected by P

A rule is **rejected** by a prefix if it doesn't correctly capture enough data
(it must correctly capture `>= c * ndata`).
Let Q be P's parent.  If Q rejects a rule, then P will also reject that rule.
This 'inheritance' of rejected rules only depends on which data are captured
by Q, and doesn't actually depend on the order of rules in Q.
Let S be the set of rules formed from (K-1) rules of P, in any order.
P inherits rejected rules from any elements of S.
Because of our **symmetry-based garbage collection** of prefixes equivalent up
to a permutation, there are at most K elements of S in the cache;
we can identify these via the **inverse canonical map (ICM)** that maps an
ordered prefix to its permutation in the cache.
We thus **lazily** initialize the list of P's **reject list** of rejected rules.

* Aside: This depends on finding elements of S, which depends on what's in the cache. When is a prefix not in the cache? Either it hasn't yet been evaluated, or it has been partially evaluated and not inserted, or evaluated and (not inserted, or inserted and later deleted). The cache is thus complemented by information that either isn't inserted or gets deleted -- are we throwing away something useful here?

Now we have a list of candidate new rules.  For each candidate rule R, we
compute the following:

(1) We compute which data are captured by R.  If R captures `< (ndata * c)`
data, then R **captured insufficient** data.
We add R to the **reject list** of P and continue to the next candidate rule.

(2) If R captures all data uncaptured by P, then R now behaves like the default
rule, but at the cost of adding a rule.
We add R to the reject list of P and continue.

(3) We compute the majority class of data captured by R,
and count how many data it correctly predicts.
If this number is `< (ndata * c)`, then R **correctly captured insufficient** data.
We add R to the reject list of P and continue.

Let P' be prefix P extended with rule R.
We compute the lower bound of the objective of R

    (number of mistakes) / ndata + c * len(P')

(4) If it's greater than the smallest observed objective,
then P' is a **dead prefix** and we continue to the next candidate.

We compute the default prediction for P', which lets us compute the objective.
If it's smaller than the best objective we've seen, we update that.

(5) The children of P' are one longer that P' (i.e., are length K + 2).
If this length exceeds M, the maximum prefix length we have to check,
then there's no point to pursuing P',
it is a **dead prefix** and we continue to the next candidate.

Also because the children of P' are one longer that P',
we can give them a tighter lower bound,

    (number of mistakes) / ndata + c * len(P') + c

(6) If this is greater than the smallest observed objective,
then there's no point in pursuing P',
it is a **dead prefix** and we continue to the next candidate.

(7) We check whether there's a permutation of P' in the ICM.

* If not, then add an entry for P' to the ICM.
* Otherwise, call the permutation T.
  * If the objective of P' is lower than T, then T is **inferior** and we replace it with P'.
  * Otherwise, P' is **inferior** and we continue to the next candidate.

If P' reaches this far, then we construct a cache entry for it.

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

* **reject_list** : when an `unused` cache entry for prefix P of length L is retrieved, this attribute is lazily initialized by taking the union of `reject_list` elements from any cached prefix in S, where S is the set of prefixes formed by taking subsets of (L-1) rules from P; this is efficient due to our symmetry-aware garbage collection and an associated data structure.

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

* **garbage_collect(min_objective)** : delete all prefixes with `lower_bound > min_objective`; called when `min_objective` decreases

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

#### Symmetry-aware pruning I

If two rules A and B capture non-intersecting subsets of data, they
**commute globally (type I)**.
If rules A and B commute, then a rule list where A and B are adjacent is
equivalent to another rule list where A and B swap positions.
More generally, a rule list containing possibly multiple, possibly overlapping,
pairs of commuting rules is equivalent to any other rule list that can be
generated by swapping one or more such pairs of rules.
We avoid evaluating multiple such equivalent rule lists by eliminating all but one.
We achieve this by restricting how we grow a prefix (by appending a rule),
as informed by a hash map `cdict` that maps each rule `R_i` to a set of rules
`S_i = {R_j}` such that `R_i` commutes with every `R_j` in `S_i` and `j > i`.
When growing a prefix that ends with rule `R_i`, we only append a rule `R_j` if
it is **not** in `S_i`.

For `tdata`, `cdict` maps each rule to a set that on average has 59 rules
(15% of 377 total).

For `adult`, `cdict` maps each rule to a set that on average has 25 rules
(8% of 284 total).

Note that we could also define a local version of (type I) commuting.

#### Symmetry-aware pruning II

If two rules A and B are adjacent in a rule list and the both predict "the same
way" (i.e., both predict "0" or both predict "1" for captured data), then they
**commute locally (II)**.
The way we handle this currently could be made more efficient.

#### Relation-aware pruning

Rule A **dominates** rule B if the data captured by B is a subset of the data
captured by A.
Thus, rule B should never follow rule A in a prefix, because it will never
capture additional data.
Similar to symmetry-aware pruning, we achieve this by restricting how we grow a
prefix, as informed by a hash map `rdict` that maps a rule `R_i` to a set of
rules `T_i = {R_k}` such that `R_i` dominates every `R_k` in `T_i`.
When growing a prefix that ends with rule `R_i`, we only append a rule `R_k` if
it is **not** in `T_i`.

Notice that the intersection of `S_i` and `T_i` is the empty set, thus the
mappings represented by `cdict` and `rdict` can easily be combined into a single
mapping.

For `tdata`, `rdict` maps each rule to a set that on average has 3 rules.
For `adult`, `rdict` maps each rule to a set that on average has 2 rules.

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

## small datasets with rule expansion

The last three columns report the number of rules mined for (max cardinality, min support)

| dataset | # data | # 0 | # 1 | f. 0 | f. 1 | # dim | (2, 0.01) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| bcancer | 683 | 444 | 239 | 0.65 | 0.35 | 28 | 1,336 |
| cars | 1,728 | 1,210 | 518 | 0.70 | 0.30 | 22 | 792 |
| haberman | 306 | 81 | 225 | 0.26 | 0.74 | 16 | 334 |
| monks1 | 432 | 216 | 216 | 0.5 | 0.5 | 18 | 396 |
| monks2 | 432 | 290 | 142 | 0.67 | 0.33 | 18 | 396 |
| monks3 | 432 | 204 | 228 | 0.47 | 0.53 | 18 | 396 |
| votes | 435 | 168 | 267 | 0.39 | 0.61 | 17 | 512 |
| compas | 7214 | 3743 | 3471 | 0.52 | 0.48 | 30 | 1,037 |

## implemented but not explained

If `c > 0`, don't add prefix to priority queue or cache if
`c * (len(prefix) + 1) >= min_objective`.

## todo

Semantics should buy us more than we're already achieving.  Let's optimize the crap out of tic-tac-toe!

From Cynthia:  If the rules make the same prediction on all of their overlapping points, then they commute.

Compute commutes as rule is added!

Curiosity doesn't seem to be helping `adult` with `c = 0.01` find the best
prefix of length 3, compared to `objective` and `breadth_first`.
Run `breadth_first` on all prefixes of at least length 4.

Does the FP growth code have ideas that would be useful to us?

Conditional commutativity?

Insertions and deletions into `pdict` (and `priority_queue`?) should happen in
parallel with analogous cache operations.

May want to completely remove dependence of incremental on cache.

Should we skip symmetry-based garbage collection (via `pdict`) when
`len(prefix) == max_prefix_len_check`?

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
