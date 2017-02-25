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
      State that is vital for the execution of the algorithm (e.g. num_nodes) should be in tree. O/w in logger.
- [ ] Can we turn it off completely or are we just not writing to a file?
- [ ] Default setting should probably turn the logger off
- [ ] The logger's actions are in some places quite fine-grained -- evaluate pros and cons

- [ ] The remaining state space calculation (fine-grain bound) probably adds noticeable overhead
      (Elaine noted 10% at some point?), so add an option to switch this off and
      instead use the coarse-grain bound -- this was the original thing implemented,
      but I think just isn't used currently

- [ ] Measure logging overhead and determine a useful heuristic threshold
      (e.g., "writing a log entry every 50 iterations incurs about 1% overhead on tdata")
      Can test this by just making all logger functions just return right away.

- [ ] Either use inheritance or macros to allow for turning logger off

### main.cc

- [ ] Can we tidy this up with a more modular structure?
- [ ] It would be nice to support more custom scheduling policies

### templates vs inheritance
- [ ] Measure the space-time tradeoff of using templates vs inheritance
- [ ] Rewrite Queue class using inheritance and benchmark changes
- [ ] Extend inheritance to other data structures and rewrite main
- [ ] Check size of executable vs runtime

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
