CC = gcc
CFLAGS = -O2 -Wall -Wextra -std=c99
LDFLAGS = -lm
BINDIR = ../bin

.PHONY: all clean

all: $(BINDIR)/strong_geon

$(BINDIR)/strong_geon: strong_geon.c | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -f $(BINDIR)/strong_geon
