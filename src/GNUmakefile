# IF DEBUGGING, CHANGE TO -O0
# IF using valgrind, add -DVAL
CFLAGS := -g -W -Wall -Wno-unused-parameter -I. -I../lib -O0 -DGMP -fno-omit-frame-pointer # @Cdefs @Pdefs
CXXFLAGS := $(CFLAGS) -std=c++11
CFLAGS += -std=c99
LIBS += -lpthread -lgmpxx -lgmp # -lprofiler -ltcmalloc
LDFLAGS := -L../lib

ifneq ($(shell uname),Darwin)
    LIBS += -lbsd
endif

CC = gcc
CXX = g++

DEPSDIR := .deps

CLEAN = corels *~ *.o *.so

all: corels

corels: run.o main.o corels.o pmap.o cache.o rulelib.o utils.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LIBS)

evaluate-test: evaluate-test.o evaluate.o corels.o pmap.o cache.o rulelib.o utils.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LIBS)

main.o: main.cc $(DEPSDIR)/stamp
	$(CXX) $(CXXFLAGS) $(DEPCFLAGS) -c $< -o $@

corels.o: corels.cc $(DEPSDIR)/stamp
	$(CXX) $(CXXFLAGS) $(DEPCFLAGS) -c $< -o $@

pmap.o: pmap.cc $(DEPSDIR)/stamp
	$(CXX) $(CXXFLAGS) $(DEPCFLAGS) -c $< -o $@

utils.o: utils.cc $(DEPSDIR)/stamp
	$(CXX) $(CXXFLAGS) $(DEPCFLAGS) -c $< -o $@

cache.o: cache.cc $(DEPSDIR)/stamp
	$(CXX) $(CXXFLAGS) $(DEPCFLAGS) -c $< -o $@

rulelib.o: rulelib.c $(DEPSDIR)/stamp
	$(CC) $(CFLAGS) $(DEPCFLAGS) -c $< -o $@

evaluate-test.o: evaluate-test.cc $(DEPSDIR)/stamp
	$(CXX) $(CFLAGS) $(DEPCFLAGS) -c $< -o $@

evaluate.o: evaluate.cc $(DEPSDIR)/stamp
	$(CXX) $(CFLAGS) $(DEPCFLAGS) -c $< -o $@

run.o: run.cc $(DEPSDIR)/stamp
	$(CXX) $(CXXFLAGS) $(DEPCFLAGS) -c $< -o $@

acceptance:
	cd ../tests_acceptance && ./test.sh

margo:
	cd .. && ./tests_acceptance/margo/test.sh

$(DEPSDIR)/stamp:
	mkdir -p $(dir $@)
	touch $@

clean:
	$(RM) $(wildcard $(CLEAN))
	$(RM) -r $(DEPSDIR)

DEPFILES := $(wildcard $(DEPSDIR)/*.d)
ifneq ($(DEPFILES),)
include $(DEPFILES)
endif

.PHONY: all clean acceptance margo
.SUFFIXES:
