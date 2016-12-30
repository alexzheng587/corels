CFLAGS := -g -W -Wall -O$(O) -I.
LIBS := -lgmpxx -lgmp

ifneq ($(shell uname -n),beepboop)
    LDFLAGS := -Wl,-no_pie
    CFLAGS += -fno-omit-frame-pointer
    ifneq ($(shell uname -n),elaines-MacBook-Air.local)
        LIBS += -lprofiler -ltcmalloc
    endif
else
    LDFLAGS := -W
endif

ifneq ($(shell uname),Darwin)
    LIBS += -lbsd
endif

CC = gcc
CXX = g++ -std=c++11

DEPSDIR := .deps

CLEAN = bbcache *~ *.o *.so

all: bbcache

bbcache: main.o bbound.o cache.o rulelib.o utils.o memtrack.o
	$(CXX) $(CFLAGS) $^ -o $@ $(LDFLAGS) $(LIBS)

main.o: main.cc $(DEPSDIR)/stamp
	$(CXX) $(CFLAGS) $(DEPCFLAGS) -c $< -o $@

bbound.o: bbound.cc $(DEPSDIR)/stamp
	$(CXX) $(CFLAGS) $(DEPCFLAGS) -c $< -o $@

utils.o: utils.cc $(DEPSDIR)/stamp
	$(CXX) $(CFLAGS) $(DEPCFLAGS) -c $< -o $@

cache.o: cache.cc $(DEPSDIR)/stamp
	$(CXX) $(CFLAGS) $(DEPCFLAGS) -c $< -o $@

memtrack.o: memtrack.cc $(DEPSDIR)/stamp
	$(CC) $(CFLAGS) $(DEPCFLAGS) -c $< -o $@

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
