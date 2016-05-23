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

* A cache to support incremental computation.

* A priority queue to manage different scheduling policies (implemented using
Python's `heapq` module).

* Symmetry-aware garbage collection for sets of prefixes that are equivalent up
to permutation -- only keep the best.

* Symmetry-aware pruning for equivalence classes of prefixes that contain
(possibly multiple) adjacent pairs of commuting rules -- only need to evaluate one.

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

## tic-tac-toe results

### tdata, breadth-first, no regularization (c = 0.)

* Aggressive "optimistic (lying) warm start" (initialize min_objective = 0.001)
* Finds a minimum length (8) perfect prefix (~200 sec)

### tdata, curiosity, no regularization (c = 0.)

* Cold start, quickly (< 1 sec) finds a perfect prefix that is very long (82)

### tdata, curiosity, regularization (c = 0.001)

* Cold start, quickly finds a perfect prefix of length 10
* Then certifies that the best (perfect) prefix has length 8 (~600 sec -- check this number)

## adult results

Subsampling 10% of dataset unless otherwise noted.

### adult, curiosity, no regularization (c = 0.)

* Certifies there are no prefixes with objective < 0.005
* Up to symmetries, ~170,000 prefixes have lower bound < 0.005 (length <= 11)

* Initialize min_objective = 0.01
* Inconclusive, > 3,000,000 prefixes have lower bound < 0.01

### adult, curiosity, regularization (c = 0.003)

* Initialize min_objective = 0.01
* Quickly certifies (< 2 sec) that there are no prefixes with objective < 0.01
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

## todo

Combine objective-based priority with restriction (prevent excessive exploration of a single subtree).

Restrict search to sub-tree.  (Could be implemented by thresholding `rules`.)

Implement back-off.

An actual pruning routine.

If prefix P is optimal, is breadth-first the best way to certify?

Rename `max_accuracy` to `best_prefix_accuracy`.

Remove prefix from cache entry (it's already the corresponding key).

Add notes about commuting pairs.

Fix the greedy algorithm to stop early based on evaluating the default rule.

Think about useful heuristics to cut down on the size of the search space.

Cynthia's optimization.

We might want to explore notions of approximate equivalence.
