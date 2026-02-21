CC = gcc
CFLAGS = -O3 -Wall -Wextra -fopenmp -march=native
LDFLAGS = -lm -fopenmp
BINDIR = ../bin

all: $(BINDIR)/bi_hopfion

$(BINDIR)/bi_hopfion: bi_hopfion.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(BINDIR)/bi_hopfion

.PHONY: all clean
