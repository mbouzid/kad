CXX = g++
CXXFLAGS = -Wall -lrocksdb -lz
OBJS = kad
HEADERS=kstring.h kseq.h

.PHONY: all

all: $(OBJS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm $(OBJS)
