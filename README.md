# bbcache
Branch-and-bound algorithm, with caching, for decision lists.

dependencies
------------

    Python 2.7x
    numpy
    tabular
    matplotlib

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

thoughts
========

We've been ignoring symmetries so far, and perhaps it's time to start thinking
about them.  If two prefixes are similar in the sense that they capture
approximately the same data and have approximately the same accuracy, then it
doesn't really make sense to pursue both of them.
