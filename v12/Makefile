CC = gcc
CFLAGS = -O3 -Wall -Wextra -fopenmp -march=native
LDFLAGS = -lm -lfftw3 -lfftw3_threads -fopenmp

SRC = src/bi_geon.c
BIN = bin/bi_geon

all: $(BIN)

$(BIN): $(SRC)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(BIN)

.PHONY: all clean
