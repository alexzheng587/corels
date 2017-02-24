# To do

## Code review ([clean-up branch](https://github.com/elaine84/bbcache/tree/clean-up))

### General

- [ ] Would we benefit from some refactoring?
- [ ] Are we fully leveraging the rule library?
- [ ] How should we include or point to the rule library?
- [ ] Are we displaying useful output during execution and at the end?
- [ ] We should add the ability to generate predictions for a test dataset on a learned model
- [ ] Would we like to support a more "plug-in" style framework (for scheduling policies, node types, and/or algorithms)?
- [ ] Can we simplify anything about how we're using templates (are we not really leveraging them in places)?

### cache.cc and cache.hh

- [ ] Anything that we should strip out of our Node types?  (Or perhaps introduce "lite" versions?)
- [ ] Any extra cruft being stored in the trie?  (See below regarding interactions with logger)

### queue.cc and queue.hh

### pmap.cc and pmap.hh

- [ ] TODO: Need to deal with the fact that CapturedKey doesn't contain length of associatd prefix.

### utils.cc and utils.hh

- [ ] Are we ok with the fact that the logger object is a global variable?
- [ ] The logger and the tree both replicate some state -- how should we eliminate redundancies?
- [ ] Can we turn it off completely or are we just not writing to a file?
- [ ] Default setting should probably turn the logger off
- [ ] We should measure the logger's overhead
- [ ] The remaining state space calculation (tighter bound) probably adds noticeable overhead
- [ ] The logger's actions are in some places quite fine-grained -- evaluate pros and cons

### main.cc

- [ ] Can we tidy this up with a more modular structure?
- [ ] It would be nice to support more custom scheduling policies

### GNUmakefile

- [ ] Eliminate machine-specific conditions
- [ ] How do we improve our makefile?

### README.md (start fresh)

- [ ] Explain all external dependencies
- [ ] Succinctly describe our algorithm and point to a paper
- [ ] We should include a small dataset as a working example

### Other

- [ ] Pick a license
- [ ] Should we make a new repo `corels` (with or without our git history)?
- [ ] Should we make R bindings?
- [ ] Should we include any scripts

### Aesthetic things

- [ ] Use `equivalent` or `identical` in names instead of `minority` or `minor`

## Short paper

- [ ] Check language regarding best-first search, priority queue and metric
- [ ] Clarify that our experiments with the priority queue use the lower bound
- [ ] Delete references (in the intro and implementation) to things we haven't done or don't evaluate?
- [ ] Make room for and add acknowledgements
- [ ] Provide github link for our code
- [ ] Revisit Johann's comments

## Long paper

### Definitions

- [ ] Stop overloading the notation for misclassification error (takes 3 arguments always)
- [ ] Swap the arguments in the notation for captures

### Bounds

- [ ] Write missing proofs for equivalent points bounds

### Data

- [ ] Of those stopped, who is frisked (or searched)?
- [ ] Maybe something about hospital readmissions
- [ ] Kaggle had a maybe interesting housing dataset
- [ ] Look at the datasets from Stat 220
- [ ] Figure out what is going on with the nursery dataset -- maybe ping Hongyu again (Elaine)
- [ ] How far did we get with adult (gap between lower bound and objective)? (Elaine)

### Implementation (writing)

- [ ] Extended version of our short paper's implementation section (Nicholas)
- [ ] Add a more detailed description of our algorithm as pseudocode, that includes all bounds and data structures (Nicholas)

### Implementation (wish list)

- [ ] Add BitVector node and a version of our algorithm (probably using new functions) that uses it

### Experiments

- [ ] Describe stop-and-frisk dataset in greater detail
- [ ] Measurements with and without our symmetry map using captured bit vectors
- [ ] Should we output confusion matrices?
- [ ] Provide more examples of rule lists that we find (Elaine)

### Evaluation (analyzing our algorithms and data structures)

- [ ] Detailed analysis of our log files.
      Start by visualizing all logged quantities as functions of time,
      using logs from our ablation study.  Explore when linear versus log scales
      make sense (e.g., matplotlib functions plot, semilogx, semilogy, loglog).
      Plot data for a single cross-validation folds separately, and together.
      Do line plots always make sense? (Nicholas)
- [ ] Analyze our timing measurements -- where does our algorithm spend its time?
      Does this change during execution?  See `eval/scratch.py` for preliminary analysis.

### General

- [ ] Curiosity is missing -- it gives significant improvement for tic-tac-toe, but what else?
- [ ] Do we want a more thorough description of how we do rule mining, and our bit vector representation and operations?

-----

**Stuff below hasn't been organized recently, but contains some useful thoughts and notes**

## bbound improvements

- [ ] Add stopping condition based on the number of iterations

- [ ] Add policy that switches from curiosity to BFS

- [ ] Measure our logging overhead and determine a useful heuristic threshold
      (e.g., "writing a log entry every 50 iterations incurs about 1% overhead on tdata")

- [ ] Make R/Python binding so it is easier to run experiments. (Daniel)
      Link to SBRL cran repo: https://github.com/cran/sbrl. We can follow this format.

- [ ] What else should we measure?  E.g., think about time spent deleting nodes
      from the cache, garbage collection triggered by a new best objective value,
      etc.

- [ ] Consolidate algorithm state in the logger object

- [ ] The lower bound check before the call to `evaluate_children` does NOT require
      a trie lookup in the special case of a curious queue ordered by lower bound,
      in which case it could be considerably faster -- maybe want to specialize `queue_select`?

- [ ] Garbage collect the permutation map of prefixes longer than the maximum prefix length

- [x] How many times do we look up an item in the hash map that is not in the tree --
      added to logger at `_state.pmap_null_num`

- [x] When we do look up an item, how many times do we end up discarding
      something we would have been unable to discard had we not done this --
      added to logger at `_state.pmap_discard_num`, counts number of lookups
      that trigger a discard operation, not the total number of deleted nodes

- [ ] If garbage collecting the cache and queue would reduce the size of the queue
      by at least some factor (e.g., 10%) then do so -- can't iterate over
      `std::priority_queue` so this would require a different data structure,
      e.g., a custom subclass

- [ ] Make permutation map garbage collection optional

- [x] Estimate the size (in memory) of the permutation map compared to the cache --
      compare experiments for a fixed `max_num_nodes` with and without the
      permutation map.  Should we consider an alphabetical tree?

- [ ] Estimate the size (in memory) of the queue -- should we try to eliminate it?
      The queue is especially annoying if it retains many prefixes marked as deleted
      and is difficult to garbage collect.  However, it seems small compared to the cache,
      e.g., consider the following policies (curious lower bound queue with permutation map
      vs. stochastic with no queue or permutation map vs. bfs queue without pmap) --
      not truly a fair comparison because the distributions of prefix lengths differ

      `./bbcache -c 2 -p 1 -r 0.001 -n 100000000 ../data/tdata_R.out ../data/tdata_R.label` (15 GB)

      `./bbcache -s -r 0.001 -n 22460425 ../data/tdata_R.out ../data/tdata_R.label` (4 GB)

      `./bbcache -b -p 0 -r 0.001 -n 22460425 ../data/tdata_R.out ../data/tdata_R.label` (4.2 GB)

- [ ] Eliminate curious lower bound queue by propagating lower bounds up the tree;
      select next prefix by traversing path to node with smallest lower bound
      (implement as a separate algorithm)

- [ ] Can we eliminate the BFS queue by implementing BFS via tree traversal?

- [ ] Implement back-off and test on adult dataset (Nicholas)

- [ ] The above bound seems to incur about a 10% overhead (in time),
      so add an option to switch this off and instead use the coarse-grain bound (Elaine)

- [ ] Calculate above bound in a separate process

- [ ] Properly calculate log10 of the remaining search space --
      see `getLogRemainingSpaceSize()` in `utils.hh` --
      you'll probably need to either write your own function or switch from `mpz` to `mpfr`

## Algorithms and data structures

A place to note things we haven't implemented, but might

- [ ] Framework to remember rejected antecedents
- [ ] Framework aware of rules that commute
- [ ] Framework aware of dominates relationships
- [ ] Ability to switch between scheduling policies
- [ ] Something like Thompson sampling using curiosity
- [ ] Priority metric that blends between breadth-first and curiosity (non-stochastic)
- [ ] Enforce that the output optimal rule list is the simplest
