CXX = g++
CXXFLAGS = -Wall -lrocksdb -ljansson -lz
OBJS = kad
HEADERS=kstring.h kseq.h
OBJECTS=kstring.o

.PHONY: all

all: $(OBJS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm $(OBJECTS) $(OBJS)
