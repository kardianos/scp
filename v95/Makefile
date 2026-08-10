# v92 build — C kernel of record + Go tools (harness, experimental kernel, viewer)
CC      ?= gcc
CFLAGS  ?= -O2 -Wall -Wextra

all: freecell go

freecell: kernel/freecell.c
	$(CC) $(CFLAGS) -o $@ $< -lm -lzstd

go:
	go build -o fabrun ./cmd/fabrun
	go build -o battery ./cmd/battery
	go build -o volview ./cmd/volview

battery-run: all
	./battery

clean:
	rm -f freecell fabrun battery volview

.PHONY: all go battery-run clean
