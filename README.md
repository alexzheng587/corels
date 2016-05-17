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


tdata with commuting
--------------------

    froot = 'tdata_R'
    min_objective = 1.
    c = 0.
    max_prefix_length = 8
    garbage_collect = True

    cache size: [1, 14, 92, 416, 1746, 8431, 49474, 361820, 983133]
    dead prefix start: [0, 0, 0, 0, 0, 0, 0, 0, 0]
    caputed zero: [0, 0, 0, 0, 27, 478, 4743, 55155, 901311]
    stunted prefix: [0, 0, 0, 0, 0, 0, 0, 0, 0]
    commutes: [0, 0, 1366, 11025, 51660, 208056, 913480, 4870817, 32341564]
    dead prefix: [0, 363, 3791, 22925, 101220, 427543, 2111344, 12507000, 99515029]
    inferior: [0, 0, 15, 134, 931, 6750, 57291, 560062, 132363]
    seconds: [0.0, 0.0, 0.04, 0.18, 0.74, 3.06, 15.74, 104.71, 718.27]
    growth: [14.0, 6.57, 4.52, 4.2, 4.83, 5.87, 7.31, 2.72]

    if {c7=x,c8=x,c9=x} then predict 1
    else if {c1=x,c4=x,c7=x} then predict 1
    else if {c4=x,c5=x,c6=x} then predict 1
    else if {c1=x,c5=x,c9=x} then predict 1
    else if {c3=x,c6=x,c9=x} then predict 1
    else if {c3=x,c5=x,c7=x} then predict 1
    else if {c1=x,c2=x,c3=x} then predict 1
    else if {c2=x,c5=x,c8=x} then predict 1
    else predict 0

    prefix: (359, 64, 264, 73, 219, 211, 47, 147)
    prediction: (1, 1, 1, 1, 1, 1, 1, 1)
    accuracy: 1.0000000000
    upper_bound: 1.0000000000
    objective: 0.0000000000
    lower_bound: 0.0000000000
    num_captured: 408
    num_captured_correct: 408
    sum(not_captured): 231
    curiosity: 0.000

tdata w/o commuting
-------------------

    cache size: [1, 14, 92, 416, 1746, 8431, 49459, 361698, 520176]
    dead prefix start: [0, 0, 0, 0, 0, 0, 0, 0, 0]
    caputed zero: [0, 0, 13, 158, 1000, 6038, 36048, 261984, 2664696]
    stunted prefix: [0, 0, 0, 0, 0, 0, 0, 0, 0]
    dead prefix: [0, 363, 5081, 33343, 150098, 623799, 2950881, 16927278, 130609440]
    inferior: [0, 0, 78, 583, 2740, 12990, 99944, 798329, 33948]
    seconds: [0.0, 0.01, 0.04, 0.34, 1.3, 5.56, 20.83, 112.93, 702.97]

tdata w/o garbage collection
----------------------------

    garbage_collect = False

    cache size: [1, 14, 170, 1842, 19890]
    seconds: [0.0, 0.02, 0.15, 2.07, 21.1]

    if {c3=o,c5=o,c7=o} then predict 0
    else if {c1=o,c5=o,c9=o} then predict 0
    else if {c3=o,c6=o,c9=o} then predict 0
    else if {c4=o,c5=o,c6=o} then predict 0
    else predict 1

adult w/10% of dataset
----------------------

    froot = 'adult_R'
    warm_start = True
    froot = 'adult_R'
    max_accuracy = 0.835
    max_prefix_length = 3
    garbage_collect = True
    seed = 0
    sample = 0.1    # 10% of the dataset

    cache size: [1, 257, 32129, 2571505]
    dead prefix start: [0, 0, 0, 0]
    caputed zero: [0, 0, 326, 72215]
    stunted prefix: [0, 0, 0, 0]
    commutes: [0, 0, 7232, 963103]
    dead prefix: [0, 27, 8041, 1313416]
    inferior: [0, 0, 25003, 4140139]
    seconds: [0.0, 0.01, 1.31, 148.49]
    growth: [257, 125, 80]

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


adult w/full dataset is similar
-------------------------------

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

code/serial_priority.py
-----------------------


adult w/objective-based prioritization
--------------------------------------

    method = 'objective'
    max_cache_size = 3000000
    froot = 'adult_R'
    max_accuracy = None
    min_objective = None
    c = 0.
    max_prefix_length = 70  # not an actual constraint
    seed = 0
    sample = 0.1

    if {capital.gain=7298LessThancapital-gain,capital.loss=capital-lossEQ0} then predict 0
    else if {education=Bachelors,marital.status=Married} then predict 0
    else if {marital.status=Married,occupation=Prof-specialty} then predict 0
    else if {marital.status=Married,occupation=Exec-managerial} then predict 0
    else if {age=Middle-aged,capital.gain=capital-gainEQ0} then predict 1
    else if {relationship=Husband,hours.per.week=Over-time} then predict 0
    else if {education=Grad-school,marital.status=Married} then predict 0
    else if {age=Middle-aged,capital.loss=capital-lossEQ0} then predict 1
    else if {education=Grad-school,sex=Male} then predict 0
    else if {age=Senior,capital.gain=capital-gainEQ0} then predict 1
    else if {age=Senior,occupation=Exec-managerial} then predict 1
    else if {capital.gain=capital-gainEQ0,hours.per.week=Part-time} then predict 1
    else if {age=Young,relationship=Own-child} then predict 1
    else if {marital.status=Never-married,relationship=Own-child} then predict 0
    else if {capital.loss=capital-lossEQ0,hours.per.week=Part-time} then predict 1
    else if {education=Assoc-degree,capital.gain=capital-gainEQ0} then predict 1
    else if {education=Assoc-degree,capital.loss=capital-lossEQ0} then predict 1
    else if {marital.status=Never-married,hours.per.week=Full-time} then predict 1
    else if {education=Grad-school,capital.gain=capital-gainEQ0} then predict 1
    else if {marital.status=Never-married,sex=Male} then predict 1
    else if {occupation=Exec-managerial,capital.gain=capital-gainEQ0} then predict 0
    else if {education=Bachelors,sex=Male} then predict 0
    else if {sex=Male,hours.per.week=Over-time} then predict 0
    else if {education=Bachelors,hours.per.week=Full-time} then predict 1
    else if {marital.status=Not-married-anymore,capital.gain=capital-gainEQ0} then predict 1
    else if {relationship=Not-in-family,sex=Male} then predict 0
    else if {education=HS-grad,marital.status=Never-married} then predict 1
    else if {education=HS-grad,sex=Female} then predict 1
    else if {race=Black,capital.gain=capital-gainEQ0} then predict 0
    else if {age=Young,hours.per.week=Full-time} then predict 1
    else if {workclass=Gov,capital.gain=capital-gainEQ0} then predict 0
    else if {age=Young,marital.status=Never-married} then predict 1
    else if {capital.gain=capital-gainEQ0,hours.per.week=Full-time} then predict 1
    else if {marital.status=Not-married-anymore,hours.per.week=Full-time} then predict 1
    else if {age=Senior,workclass=Private} then predict 0
    else if {education=Grad-school,capital.loss=capital-lossEQ0} then predict 0
    else if {occupation=Craft-repair,capital.gain=capital-gainEQ0} then predict 0
    else if {age=Senior,capital.loss=capital-lossEQ0} then predict 1
    else if {workclass=Private,capital.loss=capital-lossEQ0} then predict 1
    else predict 0
    prefix: (43, 69, 122, 121, 0, 206, 77, 1, 81, 20, 26, 49, 38, 134, 54, 57, 58, 130, 75, 136, 160, 73, 240, 67, 138, 217, 87, 91, 189, 34, 243, 35, 47, 140, 30, 76, 153, 21, 253)
    prediction: (0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1)
    accuracy: 0.8447473404
    upper_bound: 0.8447473404
    objective: 467.0000000000
    lower_bound: 467.0000000000
    num_captured: 3006
    num_captured_correct: 2539
    sum(not_captured): 2
    curiosity: 0.155

    Evaluated on the full dataset:
    prediction: (0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1)
    accuracy: 0.8327515708
    upper_bound: 0.8328513015
    objective: 5031.0000000000
    lower_bound: 5028.0000000000
    num_captured: 30064
    num_captured_correct: 25036
    sum(not_captured): 17
    curiosity: 0.167

thoughts
--------

We currently do garbage collection at the end of each round, but we could do it as we're
adding prefixes to the cache.  As prefixes grow in length, so does the number of prefixes
in each group of equivalent prefixes.

We might want to explore notions of approximate equivalence.

todo
----

Add to logs -- max accuracy, min objective, best prefix

Remove prefix from cache entry (it's already the corresponding key).

Add notes about commuting pairs.

Fix the greedy algorithm to stop early based on evaluating the default rule.

We probably want a similar stopping condition in the branch-and-bound algorithm?

Think about useful heuristics to cut down on the size of the search space.

Cynthia's optimization.

Cynthia suggested using a regularized bound, since this will be tighter.

from Hongyu
-----------

A rule list from the latest `sbrl` R pkg :

    sbrl_model <- sbrl(data, iters=20000, pos_sign="1", neg_sign="0", rule_minlen=1, rule_maxlen=2, minsupport_pos=0.10, minsupport_neg=0.10, lambda=20.0, eta=1.0, nchain=25)
    print(sbrl_model)

The rules list is :

    If      {capital.gain=7298LessThancapital-gain} (rule[46]) then positive probability = 0.01617922
    else if {education=Grad-school,marital.status=Married} (rule[76]) then positive probability = 0.24616368
    else if {age=Young,marital.status=Never-married} (rule[36]) then positive probability = 0.99606078
    else if {age=Young,capital.loss=capital-lossEQ0} (rule[34]) then positive probability = 0.92338388
    else if {education=Bachelors,marital.status=Married} (rule[68]) then positive probability = 0.33762434
    else if {marital.status=Married,occupation=Exec-managerial} (rule[123]) then positive probability = 0.46758105
    else if {education=HS-grad,marital.status=Married} (rule[85]) then positive probability = 0.70224325
    else if {education=Some-college,marital.status=Married} (rule[98]) then positive probability = 0.58744545
    else if {relationship=Own-child,hours.per.week=Full-time} (rule[224]) then positive probability = 0.98528471
    else if {marital.status=Married,occupation=Prof-specialty} (rule[124]) then positive probability = 0.27289720
    else if {education=Assoc-degree,marital.status=Married} (rule[59]) then positive probability = 0.57581069
    else if {capital.loss=capital-lossEQ0,hours.per.week=Part-time} (rule[55]) then positive probability = 0.98307380
    else if {occupation=Other-service,capital.loss=capital-lossEQ0} (rule[172]) then positive probability = 0.97980381
    else if {occupation=Prof-specialty,sex=Male} (rule[183]) then positive probability = 0.73821990
    else if {occupation=Adm-clerical,sex=Female} (rule[154]) then positive probability = 0.98090186
    else if {occupation=Prof-specialty} (rule[184]) then positive probability = 0.86271186
    else if {age=Middle-aged,education=HS-grad} (rule[4]) then positive probability = 0.96911197
    else if {education=Grad-school} (rule[81]) then positive probability = 0.61656442
    else if {age=Middle-aged,hours.per.week=Full-time} (rule[6]) then positive probability = 0.94713161
    else if {education=Bachelors,native.country=N-America} (rule[69]) then positive probability = 0.73352034
    else if {sex=Male,native.country=N-America} (rule[244]) then positive probability = 0.86269540
    else  (default rule)  then positive probability = 0.94477711

Elaine seems to have a slightly different rule list, without a rule for `{occupation=Adm-clerical,sex=Female}` :

    line numbers: (46, 78, 36, 34, 70, 122, 87, 98, 222, 123, 60, 55, 170, 181, ???, 182, ...)

    zero-indexed: (45, 77, 35, 33, 69, 121, 86, 97, 221, 122, 59, 54, 169, 180)


There is another example rule list in our earlier paper http://arxiv.org/pdf/1602.08610v1.pdf on page 3:

    if capital-gain>$7298.00 then probability to make over 50K = 0.986
    else if Young,Never-married, then probability to make over 50K = 0.003
    else if Grad-school,Married, then probability to make over 50K = 0.748
    else if Young,capital-loss=0, then probability to make over 50K = 0.072
    else if Own-child,Never-married, then probability to make over 50K = 0.015
    else if Bachelors,Married, then probability to make over 50K = 0.655
    else if Bachelors,Over-time, then probability to make over 50K = 0.255
    else if Exec-managerial,Married, then probability to make over 50K = 0.531
    else if Married,HS-grad, then probability to make over 50K = 0.300
    else if Grad-school, then probability to make over 50K = 0.266
    else if Some-college,Married, then probability to make over 50K = 0.410
    else if Prof-specialty,Married, then probability to make over 50K = 0.713
    else if Assoc-degree,Married, then probability to make over 50K = 0.420
    else if Part-time, then probability to make over 50K = 0.013
    else if Husband, then probability to make over 50K = 0.126
    else if Prof-specialty, then probability to make over 50K = 0.148
    else if Exec-managerial,Male, then probability to make over 50K = 0.193
    else if Full-time,Private, then probability to make over 50K = 0.026
    else (default rule) then probability to make over 50K = 0.066

    line numbers: (46, 36, 78, 34, 135, 70, 69, 122, 87, 83, 98, 123, 60, 114, 211, 182, 167, 258)

    zero-indexed: (45, 35, 77, 33, 134, 69, 68, 121, 86, 82, 97, 122, 59, 113, 210, 181, 166, 257)
