CXX = g++ -std=c++11

CFLAGS := -g -W -Wall -Wno-unused-parameter -I. -O$(O)

ifndef disable-gmp
	CFLAGS += -DGMP
endif

LIBS := -L/usr/local/lib -lgmp -lgmpxx -L../src/ -lcorels

tests: test-queue.cc tests-main.cc
	$(CXX) $(CFLAGS) $^ -o $@ $(LIBS)