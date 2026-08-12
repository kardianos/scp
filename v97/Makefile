# v97 build — C kernel of record + Go tools (carried from v95)
CC      ?= gcc
CFLAGS  ?= -O2 -Wall -Wextra

all: freecell go

freecell: kernel/freecell.c
	$(CC) $(CFLAGS) -o $@ $< -lm -lzstd

go:
	go build -o fabrun ./cmd/fabrun
	go build -o battery ./cmd/battery
	go build -o volview ./cmd/volview
	go build -o fcsdump ./cmd/fcsdump
	go build -o scpan ./cmd/scpan
	go build -o dipole ./cmd/dipole

battery-run: all
	./battery

clean:
	rm -f freecell fabrun battery volview fcsdump scpan dipole

.PHONY: all go battery-run clean
