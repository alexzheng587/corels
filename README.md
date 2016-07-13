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

#### Symmetry-aware pruning

Two rules A and B **commute** if they capture non-intersecting subsets of data.
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

## small datasets with rule expansion (max cardinality = 2, min/max support = 10/90%)

| dataset | # data | # 0, # 1 | 0, 1 | # dim | # rules |
| --- | --- | --- | --- | --- | --- |
| bcancer | 683 | 444, 239 | 0.65, 0.35 | 28 | 1343 |
| cars | 1728 | 1210, 518 | 0.70, 0.30 | 22 | 993 |
| haberman | 306 | 81, 225 | 0.26, 0.74 | 16 | 340 |
| monks1 | 432 | 216, 216 | 0.5, 0.5 | 18 | 612 |
| monks2 | 432 | 290, 142 | 0.67, 0.33 | 18 | 621 |
| monks3 | 432 | 204, 228 | 0.47, 0.53 | 18 | 665 |
| votes | 435 | 168, 267 | 0.39, 0.61 | 17 | 645 |

## small datasets without rule expansion (with varying amounts of regularization)

### bcancer (c = 0.01, min_captured_correct = 0.01)

    nrules: 27
    ndata: 683

    seconds: 29.18559

    if {a6>7} then predict 1
    else if {a2<5} then predict 0
    else predict 1

    accuracy: 0.9516837482
    objective: 0.0683162518
    seconds: 1.14501

### bcancer (c = 0.003, min_captured_correct = 0.)

    seconds: 44.94337

    if {a1>7} then predict 1
    else if {a8>7} then predict 1
    else if {a6>7} then predict 1
    else if {a2<5} then predict 0
    else predict 1

    accuracy: 0.9648609078
    objective: 0.0471390922

### bcancer (c = 0.001, min_captured_correct = 0.)

    seconds: 119.64070

    if {a7>7} then predict 1
    else if {a2.5-7} then predict 1
    else if {a4>7} then predict 1
    else if {a6>7} then predict 1
    else if {a1<5} then predict 0
    else if {a9.5-7} then predict 1
    else if {a8>7} then predict 1
    else if {a6.5-7} then predict 1
    else if {a1.5-7} then predict 0
    else if {a7.5-7} then predict 1
    else if {a9<5} then predict 0
    else predict 1

    accuracy: 0.9765739385
    objective: 0.0344260615

### bcancer (c = 0., min_captured_correct = 0.)

    seconds: 1009.41659

    if {a7>7} then predict 1
    else if {a2.5-7} then predict 1
    else if {a4>7} then predict 1
    else if {a6>7} then predict 1
    else if {a1<5} then predict 0
    else if {a9>7} then predict 1
    else if {a9.5-7} then predict 1
    else if {a8>7} then predict 1
    else if {a5>7} then predict 0
    else if {a6.5-7} then predict 1
    else if {a3>7} then predict 0
    else if {a2>7} then predict 1
    else if {a8.5-7} then predict 1
    else if {a5<5} then predict 0
    else if {a7.5-7} then predict 1
    else if {a1.5-7} then predict 1
    else predict 0

    accuracy: 0.9780380673
    objective: 0.0219619327

### cars (c = 0.01, min_captured_correct = 0.01)

    nrules: 21
    ndata: 1728

    seconds: 0.21441

    if {persons=2} then predict 0
    else if {safety=low} then predict 0
    else if {buying=low} then predict 1
    else if {buying=med} then predict 1
    else if {maint=vhigh} then predict 0
    else if {safety=high} then predict 1
    else if {lug-boot=small} then predict 0
    else predict 1

    accuracy: 0.9386574074
    objective: 0.1313425926

###  cars (c = 0.003, min_captured_correct = 0.)

    seconds: 0.17232

    if {persons=2} then predict 0
    else if {safety=low} then predict 0
    else if {buying=low} then predict 1
    else if {buying=med} then predict 1
    else if {maint=vhigh} then predict 0
    else if {buying=high} then predict 1
    else if {maint=high} then predict 0
    else if {safety=high} then predict 1
    else if {lug-boot=small} then predict 0
    else predict 1

    accuracy: 0.9479166667
    objective: 0.0790833333

###  cars (c = 0.001, min_captured_correct = 0.)

    seconds: 0.08900

    if {persons=2} then predict 0
    else if {safety=low} then predict 0
    else if {buying=low} then predict 1
    else if {buying=med} then predict 1
    else if {maint=vhigh} then predict 0
    else if {buying=high} then predict 1
    else if {maint=high} then predict 0
    else if {lug-boot=big} then predict 1
    else if {safety=high} then predict 1
    else if {lug-boot=small} then predict 0
    else if {doors=2} then predict 0
    else predict 1

    accuracy: 0.9502314815
    objective: 0.0607685185

###  cars (c = 0., min_captured_correct = 0.)

    seconds: 0.04759

    if {persons=2} then predict 0
    else if {safety=low} then predict 0
    else if {buying=low} then predict 1
    else if {buying=med} then predict 1
    else if {maint=vhigh} then predict 0
    else if {safety=high} then predict 1
    else if {lug-boot=small} then predict 0
    else if {buying=high} then predict 1
    else if {maint=high} then predict 0
    else if {lug-boot=big} then predict 1
    else if {doors=2} then predict 0
    else if {persons=more} then predict 1
    else if {doors=3} then predict 0
    else predict 1

    accuracy: 0.9513888889
    objective: 0.0486111111

### haberman (c = 0.01, min_captured_correct = 0.01)

    nrules: 15
    ndata: 306

    seconds: 4.09679

    if {age<40} then predict 1
    else if {nodes10-19} then predict 0
    else predict 1

    accuracy: 0.7614379085
    objective: 0.2585620915

### haberman (c = 0.003, min_captured_correct = 0.)

    seconds: 35.42373

    if {age<40} then predict 1
    else if {nodes0} then predict 1
    else if {year>65} then predict 1
    else if {nodes10-19} then predict 0
    else if {age>69} then predict 0
    else if {nodes1-9} then predict 1
    else if {year62-63} then predict 0
    else if {age60-69} then predict 1
    else predict 0

    accuracy: 0.7941176471
    objective: 0.2298823529

### haberman (c = 0.001, min_captured_correct = 0.)

    seconds: 46.61386

    same solution as above

### haberman (c = 0., min_captured_correct = 0.)

    ~same as above

### monks1 (c = 0.01, min_captured_correct = 0.01)

    nrules: 17
    ndata: 432

    seconds: 0.06728

    if {a5=1} then predict 1
    else if {a1=2} then predict 0
    else if {a1=1} then predict 0
    else if {a2=3} then predict 1
    else predict 0

    accuracy: 0.8333333333
    objective: 0.2066666667

### monks1 (c = 0.003, 0.001, 0., min_captured_correct = 0.)

    ~same as above

###  monks2 (c = 0.01, min_captured_correct = 0.01)

    nrules: 17
    ndata: 432

    seconds: 5.79060

    if {a1=1} then predict 0
    else if {a3=2} then predict 0
    else if {a2=1} then predict 0
    else if {a4=1} then predict 0
    else if {a6=1} then predict 1
    else if {a5=1} then predict 1
    else predict 0

    accuracy: 0.7268518519
    upper_bound: 0.7268518519

### monks2 (c = 0.003, 0.001, 0., min_captured_correct = 0.)

    ~same as above

### monks3 (c = 0.01, min_captured_correct = 0.01)

    nrules: 17
    ndata: 432

    seconds: 0.01394

    if {a5=4} then predict 0
    else if {a2=3} then predict 0
    else predict 1

    accuracy: 0.9722222222
    objective: 0.0477777778

### monks3 (c = 0.003, min_captured_correct = 0.)

    seconds: 0.00667

    if {a5=4} then predict 0
    else if {a2=2} then predict 1
    else if {a2=1} then predict 1
    else if {a4=3} then predict 0
    else if {a4=2} then predict 0
    else if {a5=3} then predict 1
    else predict 0

    accuracy: 1.0000000000
    objective: 0.0180000000

### monks3 (c = 0.001, min_captured_correct = 0.)

    seconds: 0.00750

    if {a5=4} then predict 0
    else if {a2=2} then predict 1
    else if {a2=1} then predict 1
    else if {a4=3} then predict 0
    else if {a4=2} then predict 0
    else if {a5=3} then predict 1
    else predict 0

    accuracy: 1.0000000000
    objective: 0.0060000000

### monks3 (c = 0., min_captured_correct = 0.)

    ~same as above

### votes (c = 0.01, min_captured_correct = 0.01)

    nrules: 16
    ndata: 435

    if {V4} then predict 0
    else predict 1

    accuracy: 0.9563218391
    objective: 0.0536781609
    seconds: 0.00188

### votes (c = 0.003, min_captured_correct = 0.)

    ~same as above

### votes (c = 0.001, min_captured_correct = 0.)

    seconds: 0.38524

    if {V4} then predict 0
    else if {V3} then predict 1
    else if {V15} then predict 1
    else if {V5} then predict 1
    else if {V11} then predict 1
    else if {V1} then predict 0
    else if {V12} then predict 0
    else if {V9} then predict 1
    else predict 0

    accuracy: 0.9655172414
    objective: 0.0424827586

### votes (c = 0., min_captured_correct = 0.)

    ~same as above

## tic-tac-toe results

### tdata, breadth-first, no regularization (c = 0.)

* Aggressive "optimistic (lying) warm start" (initialize min_objective = 0.001)
* Finds a minimum length (8 rules) perfect prefix (< 180 sec)

### tdata, curiosity, no regularization (c = 0.)

* Cold start, quickly (< 2 sec) finds a perfect prefix that is very long (68 rules)

### tdata, curiosity, regularization (c = 0.001)

* Cold start, quickly (< 2 sec) finds a perfect prefix of length 10
* Then quickly finds perfect prefixes of length 9, then 8 (< 3 sec)
* Then certifies that the best (perfect) prefix has length 8 (< 225 sec)
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

### adult, breadth-first, aggressive regularization (min_captured_correct = c)

Reporting cache entries by prefix length (for certified runs, these are entries with lower bound less than the best objective).  Also reporting "growth rate."

* (c = 0.06) Certifies (~180 sec) the best prefix is (69,), max cache size ~130,000, max length checked = 3

    [1, 204, 13513, 114479] -> [204, 66.2, 8.47]

* (c = 0.05) Certifies (~900 sec) the best prefix is (69,), max cache size ~580,000, max length checked = 4

    [1, 220, 18561, 455669, 104964] -> [220, 84.3, 24.5, 0.23]

* (c = 0.04) (>1,000 sec) the best known prefix is (69,), max cache size >5,000,000, max length checked >= 4 (certified there are no prefixes of length 2 or 3 that are better)

    [1, 236, 22893, 957038, >4019841] -> [236, 97.0, 41.8, ...]

* (c = 0.01) the best known prefix (43, 122, 121), objective = 0.20088, max cache size >3,000,000 max length checked >=4 (certified this is the best prefix of length 3)

    [1, 263, 32400, 2478819, >488533 ...] -> [263, 123, 76.5, ...]

### adult, (min_captured_correct = c = 0.01), max cache size = 2,600,000

Reporting number of cache leaves by prefix length (things in priority queue) starting at 1

* breadth-first: 0.20088, completed depth = 3

    [0, 0, 2476944, 88552]

* curiosity: 0.22148, explores too deep too fast, e.g., many prefixes of length 1 unexplored, including (43,)

    [237, 7272, 238092, 1449381, 891552]

* lower bound: 0.20420, explores up to depth 4

    [126, 21855, 512618, 2052639]

* objective: 0.20088, explores up to depth 6

    [262, 228, 1137, 4851, 197921, 2366233]

### adult, curiosity, aggressive regularization (min_captured_correct = c)

* (c = 0.1) Certifies (< 1.5 sec) the best prefix is (51,) -> 1.2 sec
* (c = 0.09) Certifies (< 5 sec) the best prefix is (118,) -> 5.3 sec
* (c = 0.08) Certifies (< 10 sec) the best prefix is (118,) -> 9.6 sec
* (c = 0.07) Certifies (< 25 sec) the best prefix is (118,) -> 23 sec
* (c = 0.06) Certifies (< 210 sec) the best prefix is (69,)
* (c = 0.05) Certifies (< 1,060 sec) the best prefix is (69,)
* (c = 0.04) (< 40,000 sec) the best prefix is (69,)

* (c = 0.03) (cache reaches 3 x 10^6 entries in 170 sec)
* (c = 0.02) (cache reaches 3 x 10^6 entries in 131 sec)

* (c = 0.01) (cache reaches 3 x 10^6 entries in 136 sec)

    (0,)    0.26432 0.74568
    (11,)   0.26199 0.74801
    (23,)   0.25601 0.75399
    (43,)   0.22742 0.78258
    (128, 43, 69)   0.22282 0.80718
    (33, 43, 69)    0.22182 0.80818
    (41, 43, 69)    0.22149 0.80851

* (c = 0.0) (cache reaches 3 x 10^6 entries in 102 sec)
* `lower_bound` priority metric performs similarly

    (0,)    0.25432 0.74568
    (11,)   0.25199 0.74801
    (23,)   0.24601 0.75399
    (43,)   0.21742 0.78258
    (114, 59)   0.21676 0.78324
    (114, 69)   0.20180 0.79820
    (114, 38, 69)   0.20146 0.79854
    (114, 236, 69)  0.19648 0.80352
    (114, 236, 38, 69)  0.19614 0.80386
    (114, 236, 38, 87, 43, 69)  0.19249 0.80751
    (114, 236, 38, 87, 43, 267, 69) 0.19149 0.80851
    (114, 236, 38, 87, 43, 34, 69)  0.19116 0.80884
    (43, 38, 87, 267, 34, 69)   0.19016 0.80984

* (c = 0.01)

  * Certifies (< 12 sec) that there are no prefixes with objective < 0.04
  * Up to symmetries, 0 prefixes have lower bound < 0.04

  * Certifies (< 60 sec) that there are no prefixes with objective < 0.05
  * Up to symmetries, 42 prefixes have lower bound < 0.05 (length <= 2)

  * Certifies (< 650 sec) that there are no prefixes with objective < 0.06
  * Up to symmetries, 569 prefixes have lower bound < 0.06 (length <= 3)

### adult, objective

* (c = 0.) quickly (< 20 sec) finds a good prefix
* (objective = 0.15525, accuracy = 0.84475, length = 39)
* (43, 69, 122, 121, 0, 206, 77, 1, 81, 20, 26, 49, 38, 134, 54, 57, 58, 130, 75, 136, 160, 73, 240, 67, 138, 217, 87, 91, 189, 34, 243, 35, 47, 140, 30, 76, 153, 21, 253)
* No improvement after 500 sec

* Full dataset (c = 0.)
* quickly (~1.0 sec) finds a good prefix (accuracy = 0.83431, length = 36)
* (43, 69, 122, 121, 77, 0, 62, 1, 46, 97, 20, 75, 21, 47, 32, 33, 162, 48, 50, 266, 49, 54, 58, 120, 66, 123, 185, 76, 129, 84, 260, 53, 95, 139, 244, 52)
* No improvement after 500 sec

* Full dataset (c = 0.0001)
* (43, 69, 122, 121, 77, 0, 62, 3, 7, 46, 97)
* (objective = 0.16675, accuracy = 0.83435, seconds < 1)

* Full dataset (c = 0.001)
* (43, 69, 122, 77, 46, 118, 97)
* (objective = 0.17169, accuracy: 0.83531, seconds < 100)

* (c = 0.003) quickly (< 1 sec) finds a good prefix
* (43, 69, 122, 121) (objective = 0.17656, accuracy = 0.83544)
* cache reaches 3 x 10^6 entries in 201 sec

* (c = 0.01) quickly (< 4 sec) finds a good prefix
* (43, 122, 121) (objective = 0.20088, accuracy = 0.82912)
* cache reaches 3 x 10^6 entries in 216 sec

* (c = 0.02) quickly (< 1 sec) finds a good prefix
* (122, 121) (objective = 0.23049, accuracy = 0.80951)
* cache reaches 3 x 10^6 entries in 205 sec

### adult, objective, small subsample (sample = 0.017)

* (c = 0.) (< 14 sec, objective = 0.09198, accuracy = 0.90802, length = 54)
* (122, 121, 43, 69, 59, 0, 4, 2, 12, 77, 8, 16, 49, 73, 32, 33, 54, 57, 130, 81, 75, 76, 162, 26, 91, 87, 128, 174, 175, 67, 254, 65, 144, 150, 6, 18, 1, 96, 102, 147, 97, 94, 168, 219, 244, 243, 156, 143, 85, 120, 196, 276, 105, 182)
* No improvement after 500 sec, cache stays small (< 2000)

### adult, curiosity, warm start

* (c = 0.) (0.15525)
* (c = 0.003) (0.17656)
* (c = 0.01) (0.20088)

## implemented but not explained

If `c > 0`, don't add prefix to priority queue or cache if
`c * (len(prefix) + 1) >= min_objective`.

## todo

Curiosity doesn't seem to be helping `adult` with `c = 0.01` find the best
prefix of length 3, compared to `objective` and `breadth_first`.
Run `breadth_first` on all prefixes of at least length 4.

Conditional commutativity?

Symmetry wrt not_captured?

Insertions and deletions into `pdict` (and `priority_queue`?) should happen in
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
