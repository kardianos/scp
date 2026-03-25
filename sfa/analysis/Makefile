CC = gcc
CFLAGS = -O3 -fopenmp
LDFLAGS = -lzstd -lm
SFA_DEF = -DSFA_IMPLEMENTATION

TOOLS = analyze_sfa spatial_analysis cluster_profile stable_vs_unstable shell_analysis modify_sfa freq_phase accel_analysis breathing

all: $(TOOLS)

%: %.c ../format/sfa.h
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(TOOLS)

.PHONY: all clean
