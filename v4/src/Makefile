CC = gcc
CFLAGS = -O2 -Wall -Wextra -std=c11
LDFLAGS = -lm

BINDIR = ../bin
DATADIR = ../data

all: $(BINDIR)/qball $(BINDIR)/ym_breather $(BINDIR)/boson_star $(BINDIR)/soler $(BINDIR)/charged_soler $(BINDIR)/ed_soler $(BINDIR)/edm_soler $(BINDIR)/ebid

$(BINDIR)/qball: qball.c
	@mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -o $@ qball.c $(LDFLAGS)

$(BINDIR)/ym_breather: ym_breather.c
	@mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -o $@ ym_breather.c $(LDFLAGS)

$(BINDIR)/boson_star: boson_star.c
	@mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -o $@ boson_star.c $(LDFLAGS)

$(BINDIR)/soler: soler.c
	@mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -o $@ soler.c $(LDFLAGS)

$(BINDIR)/charged_soler: charged_soler.c
	@mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -o $@ charged_soler.c $(LDFLAGS)

$(BINDIR)/ed_soler: ed_soler.c
	@mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -o $@ ed_soler.c $(LDFLAGS)

$(BINDIR)/edm_soler: edm_soler.c
	@mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -o $@ edm_soler.c $(LDFLAGS)

$(BINDIR)/ebid: ebid.c
	@mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -o $@ ebid.c $(LDFLAGS)

run_scan: $(BINDIR)/qball
	@mkdir -p $(DATADIR)
	../bin/qball -outdir $(DATADIR)

run_profile: $(BINDIR)/qball
	@mkdir -p $(DATADIR)
	../bin/qball -omega 0.5 -outdir $(DATADIR)
	../bin/qball -omega 0.8 -outdir $(DATADIR)
	../bin/qball -omega 0.95 -outdir $(DATADIR)

clean:
	rm -f $(BINDIR)/qball $(BINDIR)/ym_breather $(BINDIR)/boson_star $(BINDIR)/soler $(BINDIR)/charged_soler $(BINDIR)/ed_soler $(BINDIR)/edm_soler $(BINDIR)/ebid

.PHONY: all clean run_scan run_profile
