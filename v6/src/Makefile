CC = gcc
CFLAGS = -O2 -Wall -Wextra -std=c99
LDFLAGS = -lm
BINDIR = ../bin

.PHONY: all clean depletion

all: depletion

depletion: $(BINDIR)/depletion

$(BINDIR)/depletion: depletion.c | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -f $(BINDIR)/depletion
