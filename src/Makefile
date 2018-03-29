CXX = g++
CXXFLAGS = -Wall -O2 -Wno-unused-function -std=c++11
LDFLAGS = -lrocksdb -lz 
OBJS = kad
HEADERS=kstring.h kseq.h

.PHONY: all

all: kad

#%.o: %.c
kad: kad.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm $(OBJS)
