# v81/path1_locks — monist PIC sandbox
CC      ?= gcc
CFLAGS  ?= -O3 -march=native -std=c11 -Wall -Wextra -fopenmp
LDFLAGS ?= -fopenmp -lm

SRC     := src/locks_pic.c
BIN     := bin/locks_pic

.PHONY: all clean run run-l0 run-l1 run-l2 run-l3

all: $(BIN)

$(BIN): $(SRC) | bin
	$(CC) $(CFLAGS) -o $@ $(SRC) $(LDFLAGS)

bin:
	mkdir -p bin results

clean:
	rm -f $(BIN)

run: $(BIN)
	./$(BIN) all results

run-l0: $(BIN)
	./$(BIN) l0 results

run-l1: $(BIN)
	./$(BIN) l1 results

run-l2: $(BIN)
	./$(BIN) l2 results

run-l3: $(BIN)
	./$(BIN) l3 results
