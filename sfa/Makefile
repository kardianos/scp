CC = gcc
NVCC = nvcc
CFLAGS = -O3 -fopenmp
NVFLAGS = -O3 -arch=sm_70
LDFLAGS = -lzstd -lm
SDL_FLAGS = $(shell sdl2-config --cflags --libs 2>/dev/null)

# ---- Analysis tools ----

ANALYSIS_TOOLS = analyze_sfa spatial_analysis cluster_profile stable_vs_unstable \
                 shell_analysis modify_sfa freq_phase accel_analysis breathing sfa_info

ANALYSIS_BINS = $(addprefix analysis/, $(ANALYSIS_TOOLS))

analysis/%: analysis/%.c format/sfa.h
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

# ---- Seed generators ----

SEED_TOOLS = gen_braid gen_oscillon
SEED_BINS = $(addprefix seed/, $(SEED_TOOLS))

seed/%: seed/%.c format/sfa.h
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

# ---- Simulation kernels ----

sim/scp_sim: sim/scp_sim.c format/sfa.h
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

sim/scp_multi: sim/scp_multi.c format/sfa.h
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

# CUDA targets (only if nvcc available)
sim/scp_sim_cuda: sim/scp_sim.cu format/sfa.h
	$(NVCC) $(NVFLAGS) -o $@ $< $(LDFLAGS)

sim/scp_multi_cuda: sim/scp_multi.cu format/sfa.h
	$(NVCC) $(NVFLAGS) -o $@ $< $(LDFLAGS)

# ---- Viewer ----

viewer/volview: viewer/volview.c format/sfa.h
	$(CC) $(CFLAGS) -o $@ $< $(SDL_FLAGS) $(LDFLAGS)

# ---- Targets ----

.PHONY: all analysis seeds sim viewer cuda clean help

all: analysis seeds sim viewer

analysis: $(ANALYSIS_BINS)

seeds: $(SEED_BINS)

sim: sim/scp_sim sim/scp_multi

cuda: sim/scp_sim_cuda sim/scp_multi_cuda

viewer: viewer/volview

clean:
	rm -f $(ANALYSIS_BINS) $(SEED_BINS)
	rm -f sim/scp_sim sim/scp_multi sim/scp_sim_cuda sim/scp_multi_cuda
	rm -f viewer/volview

help:
	@echo "SFA Build System"
	@echo ""
	@echo "  make            Build all CPU targets (analysis + seeds + sim + viewer)"
	@echo "  make analysis   Build analysis tools (10 tools in analysis/)"
	@echo "  make seeds      Build seed generators (gen_braid, gen_oscillon)"
	@echo "  make sim        Build CPU simulation kernels (scp_sim, scp_multi)"
	@echo "  make cuda       Build CUDA simulation kernels (requires nvcc)"
	@echo "  make viewer     Build volume viewer (requires SDL2)"
	@echo "  make clean      Remove all binaries"
	@echo ""
	@echo "Individual tools:  make analysis/sfa_info  make seed/gen_braid  etc."
