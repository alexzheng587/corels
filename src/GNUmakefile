# CFLAGS := -g -W -Wall -Wno-unused-parameter -I. -DGMP -O2 -fno-omit-frame-pointer
CFLAGS := -g -W -Wall -Wno-unused-parameter -I. -DGMP -O2 -fno-omit-frame-pointer -fPIC

ifneq (,$(filter $(shell uname -n),elaines-MacBook-Pro.local Nicholass-Air-5.wireless.yale.internal elaines-MBP))
    LIBS := -L/usr/local/lib
    CFLAGS += -I/usr/local/include
else
    LIBS :=
endif

LIBS += -L/usr/local/lib -lgmpxx -lgmp

ifneq ($(shell uname -n),beepboop)
    LDFLAGS :=
    CFLAGS += -fno-omit-frame-pointer
    ifneq ($(shell uname -r),15.6.0)
#        LIBS += -lprofiler -ltcmalloc
    endif
else
    LDFLAGS :=
endif

CXXFLAGS := $(CFLAGS) -std=c++14

ifneq ($(shell uname),Darwin)
    LIBS += -lbsd
endif

CC = gcc
CXX = g++ -std=c++11
#CC = clang
#CXX = clang++ -std=c++11

DEPSDIR := .deps

CLEAN = corels *~ *.o *.so *.a

all: corels libcorels.so

corels: main.o corels.o loss.o pmap.o cache.o rulelib.o utils.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LIBS)

libcorels.so: run.o corels.o loss.o pmap.o cache.o rulelib.o utils.o
	$(CXX) -shared -o $@ $^ $(LDFLAGS) $(LIBS)

libcorels.a: run.o corels.o loss.o pmap.o cache.o rulelib.o utils.o
	ar rcs $@ $^

main.o: main.cc $(DEPSDIR)/stamp
	$(CXX) $(CXXFLAGS) $(DEPCFLAGS) -c $< -o $@

run.o: run.cc $(DEPSDIR)/stamp
	$(CXX) $(CXXFLAGS) $(DEPCFLAGS) -c $< -o $@

corels.o: corels.cc $(DEPSDIR)/stamp
	$(CXX) $(CXXFLAGS) $(DEPCFLAGS) -c $< -o $@

loss.o: loss.cc $(DEPSDIR)/stamp
	$(CXX) $(CXXFLAGS) $(DEPCFLAGS) -c $< -o $@

pmap.o: pmap.cc $(DEPSDIR)/stamp
	$(CXX) $(CXXFLAGS) $(DEPCFLAGS) -c $< -o $@

utils.o: utils.cc $(DEPSDIR)/stamp
	$(CXX) $(CXXFLAGS) $(DEPCFLAGS) -c $< -o $@

cache.o: cache.cc $(DEPSDIR)/stamp
	$(CXX) $(CXXFLAGS) $(DEPCFLAGS) -c $< -o $@

rulelib.o: rulelib.c $(DEPSDIR)/stamp
	$(CC) $(CFLAGS) $(DEPCFLAGS) -c $< -o $@

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

.PHONY: all clean
.SUFFIXES:
