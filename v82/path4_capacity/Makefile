CC = gcc
CFLAGS = -O3 -march=native -std=c11 -Wall -Wextra -fopenmp
LDFLAGS = -fopenmp -lm

all: bin/capacity_pic

bin/capacity_pic: src/capacity_pic.c
	mkdir -p bin results
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f bin/capacity_pic

.PHONY: all clean phaseA
phaseA: bin/capacity_pic
	./bin/capacity_pic phaseA results
