CXX = g++
CXXFLAGS = -Wall -lleveldb -ljansson -lz
OBJS = kad
HEADERS=kstring.h kseq.h
OBJECTS=kstring.o

.PHONY: all

all: $(OBJS)

kad: kad.cc $(HEADERS) $(OBJECTS)
# kad: kad.cc
# 	g++ -Wall -std=c++11 -o $@ $< -lleveldb

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
