CXX = g++
CXXFLAGS = -Wall -lrocksdb -lz -O2 -Wno-unused-function -std=c++11
OBJS = kad
HEADERS=kstring.h kseq.h

.PHONY: all

all: $(OBJS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm $(OBJS)
