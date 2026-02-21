CC = gcc
CFLAGS = -O2 -Wall -Wextra -std=c11
LDFLAGS = -lm

BINDIR = ../bin
DATADIR = ../data

all: $(BINDIR)/soler

$(BINDIR)/soler: soler.c
	@mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -o $@ soler.c $(LDFLAGS)

clean:
	rm -f $(BINDIR)/soler

.PHONY: all clean
