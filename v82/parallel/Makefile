CC = gcc
CFLAGS = -O3 -std=c11 -Wall -Wextra -lm

all: bin/three_approaches

bin/three_approaches: three_approaches.c
	mkdir -p bin results
	$(CC) $(CFLAGS) -o $@ $< -lm

clean:
	rm -f bin/three_approaches

.PHONY: all clean run
run: all
	./bin/three_approaches all results
