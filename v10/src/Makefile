CC = gcc
CFLAGS = -O2 -Wall -Wextra -std=c99
LDFLAGS = -lm
BINDIR = ../bin

.PHONY: all clean

all: $(BINDIR)/tov_geon

$(BINDIR)/tov_geon: tov_geon.c | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -f $(BINDIR)/tov_geon
