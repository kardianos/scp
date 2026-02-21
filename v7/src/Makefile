CC = gcc
CFLAGS = -O2 -Wall -Wextra -std=c99
LDFLAGS = -lm
BINDIR = ../bin

.PHONY: all clean breather

all: $(BINDIR)/selftrap $(BINDIR)/breather

breather: $(BINDIR)/breather

$(BINDIR)/selftrap: selftrap.c | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

$(BINDIR)/breather: breather.c | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -f $(BINDIR)/selftrap $(BINDIR)/breather
