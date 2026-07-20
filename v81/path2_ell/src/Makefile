# v81/path2_ell — CPU sandbox
CC = gcc
CFLAGS = -O2 -Wall -Wextra -std=c11 -lm
# -lm must be at end for some linkers
LDFLAGS = -lm

BIN = path2_ell
SRC = path2_ell.c

.PHONY: all clean run

all: $(BIN)

$(BIN): $(SRC)
	$(CC) -O2 -Wall -Wextra -std=c11 -o $(BIN) $(SRC) $(LDFLAGS)

run: $(BIN)
	./$(BIN) all ../out

clean:
	rm -f $(BIN)
