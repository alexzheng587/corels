# bbcache
Branch-and-bound algorithm, with caching, for decision lists.

Python dependencies
-------------------

    Python 2.7x
    numpy
    tabular
    matplotlib
    gmpy2

Dependencies of gmpy2
---------------------

These can be installed on Mac OS X with `brew install`.

    mpfr
    libmpc

tic-tac-toe training dataset
----------------------------

    data/tdata_R.out :  set of rules
    data/tdata_R.labels :  labels

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

code/branch_bound.py
--------------------

This module contains functions and data structures used by variants of the
branch-and-bound algorithm, including an object for a cache element.

code/serial.py
--------------

This module contains a serial implementation of the branch-and-bound algorithm,
with a cache to support incremental computation.  Prefixes are added to the
queue greedily, which causes the queue to grow exponentially fast.  It's here
because it's a bit easier to understand.

code/serial_lazy.py
-------------------

This module contains a serial implementation of the branch-and-bound algorithm,
with a cache to support incremental computation and a queue that grows lazily.
The queue's length is bounded by the total number of rules.

How many prefixes get stored in the cache?  We empirically measure this, for
different (simulated) warm starts.  Note that to simulate a warm start with a
particular accuracy, we don't actually need to generate a prefix with that
accuracy -- we simply set max_accuracy to the desired level.  For the
tic-tac-toe dataset, we know that a perfect rule list can be generated from the
given rules.

<!-- These are somewhat wrong
            warm=0.999  warm=0.99   warm=0.91   maximum
    len=0   1           1           1           1
    len=1   14          26          351         377
    len=2   171         592         97309       141752
    len=3   1856        12563       ?           53157000
    len=4   20061       259766      ?           19880718000
    len=5   243503      ?           ?           7415507814000

            warm=0.999  warm=0.99   warm=0.91   maximum
    len=0   1           1           1           1
    len=1   14          26          351         377
    len=2   x13         x23         x278        x376
    len=3   x11         x22         ?           x375
    len=4   x11         x21         ?           x374
    len=5   x12         ?           ?           x373
 -->

code/serial_gc.py
-----------------

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

    ############################################################################
    froot = 'tdata_R'
    warm_start = True
    max_accuracy = 0.999
    max_prefix_length = 8
    garbage_collect = True

    cache size: [1, 14, 92, 416, 1746, 8431, 49459, 361698, 520176]
    gc size: [0, 0, 78, 583, 2740, 12990, 99944, 798329, 33948]
    seconds: [0.0, 0.64, 0.1, 0.67, 2.98, 14.05, 57.75, 288.19, 1931.9]

    if {c4=x,c5=x,c6=x} then predict 1
    else if {c7=x,c8=x,c9=x} then predict 1
    else if {c1=x,c2=x,c3=x} then predict 1
    else if {c2=x,c5=x,c8=x} then predict 1
    else if {c3=x,c6=x,c9=x} then predict 1
    else if {c1=x,c4=x,c7=x} then predict 1
    else if {c3=x,c5=x,c7=x} then predict 1
    else if {c1=x,c5=x,c9=x} then predict 1
    else predict 0
    prefix: (53, 304, 19, 314, 281, 110, 333, 29)
    prediction: (1, 1, 1, 1, 1, 1, 1, 1)
    accuracy: 1.000
    upper_bound: 1.000
    num_captured: 408
    num_captured_correct: 408
    sum(not_captured): 231
    curiosity: 0.000

    ############################################################################
    froot = 'tdata_R'
    warm_start = True
    max_accuracy = 0.999
    max_prefix_length = 4
    garbage_collect = False

    cache size: [1, 14, 170, 1842, 19890]
    seconds: [0.0, 0.02, 0.15, 2.07, 21.1]

    if {c3=o,c5=o,c7=o} then predict 0
    else if {c1=o,c5=o,c9=o} then predict 0
    else if {c3=o,c6=o,c9=o} then predict 0
    else if {c4=o,c5=o,c6=o} then predict 0
    else predict 1

    ############################################################################
    froot = 'adult_R'
    warm_start = True
    froot = 'adult_R'
    max_accuracy = 0.834
    max_prefix_length = 3
    garbage_collect = True
    seed = 0
    sample = 0.1    # 10% of the dataset

    cache size: [1, 257, 32217, 2584429]
    dead prefix start: [0, 0, 0, 0]
    caputed zero: [0, 0, 581, 143107]
    stunted prefix: [0, 0, 0, 0]
    dead prefix: [0, 27, 8516, 1387673]
    inferior: [0, 0, 31417, 4969985]
    seconds: [0.0, 0.0, 1.22, 163.95]

    if {capital.gain=7298LessThancapital-gain,capital.loss=capital-lossEQ0} then predict 0
    else if {marital.status=Married,occupation=Prof-specialty} then predict 0
    else if {marital.status=Married,occupation=Exec-managerial} then predict 0
    else predict 1

    prefix: (168, 106, 199)
    prediction: (0, 0, 0)
    accuracy: 0.829
    upper_bound: 0.955
    num_captured: 519
    num_captured_correct: 385
    sum(not_captured): 2489
    curiosity: 0.258


    ############################################################################
    froot = 'adult_R'
    warm_start = True
    froot = 'adult_R'
    max_accuracy = 0.834
    max_prefix_length = 3
    garbage_collect = True
    # full dataset

    cache size: [1, 258, 32232, 2594146]
    dead prefix start: [0, 0, 0, 0]
    caputed zero: [0, 0, 577, 140941]
    stunted prefix: [0, 0, 0, 0]
    dead prefix: [0, 26, 8736, 1348177]
    inferior: [0, 0, 31469, 5006160]
    seconds: [0.0, 4.91, 4.16, 640.08]

    prefix: (168, 106, 199)
    prediction: (0, 0, 0)
    accuracy: 0.826
    upper_bound: 0.953
    num_captured: 5028
    num_captured_correct: 3622
    sum(not_captured): 25053
    curiosity: 0.280


thoughts
--------

We currently do garbage collection at the end of each round, but we could do it as we're
adding prefixes to the cache.  As prefixes grow in length, so does the number of prefixes
in each group of equivalent prefixes.

We might want to explore notions of approximate equivalence.

todo
----

Fix the greedy algorithm to stop early based on evaluating the default rule.

We probably want a similar stopping condition in the branch-and-bound algorithm?

Think about useful heuristics to cut down on the size of the search space.

Cynthia's optimization.

Cynthia suggested using a regularized upper bound, since this will be tighter.
