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

SEED_TOOLS = gen_braid gen_oscillon gen_phase_confined gen_proton_analytical gen_composite
SEED_BINS = $(addprefix seed/, $(SEED_TOOLS))

seed/%: seed/%.c format/sfa.h
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

# ---- Tools ----

TOOL_BINS = tools/sfa_render

tools/sfa_render: tools/sfa_render.c format/sfa.h
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS) -lwebp -lwebpmux -lsharpyuv

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

# ---- Runner (Go) ----

RUNNER_SRC = $(wildcard runner/*.go)

runner/scp-runner: $(RUNNER_SRC)
	cd runner && go build -o scp-runner .

# ---- Targets ----

BINDIR = ../bin

.PHONY: all analysis seeds sim tools runner viewer cuda install clean help

all: analysis seeds sim tools runner

analysis: $(ANALYSIS_BINS)

seeds: $(SEED_BINS)

sim: sim/scp_sim sim/scp_multi

tools: $(TOOL_BINS)

runner: runner/scp-runner

cuda: sim/scp_sim_cuda sim/scp_multi_cuda

viewer: viewer/volview

install: all
	@mkdir -p $(BINDIR)
	@for t in $(ANALYSIS_BINS); do cp $$t $(BINDIR)/; done
	@for t in $(SEED_BINS); do cp $$t $(BINDIR)/; done
	@cp sim/scp_sim $(BINDIR)/
	@for t in $(TOOL_BINS); do cp $$t $(BINDIR)/; done
	@cp runner/scp-runner $(BINDIR)/
	@echo "Installed $$(ls $(BINDIR) | wc -l) binaries to $(BINDIR)/"

clean:
	rm -f $(ANALYSIS_BINS) $(SEED_BINS) $(TOOL_BINS)
	rm -f sim/scp_sim sim/scp_multi sim/scp_sim_cuda sim/scp_multi_cuda
	rm -f viewer/volview runner/scp-runner
	rm -rf $(BINDIR)

help:
	@echo "SFA Build System"
	@echo ""
	@echo "  make            Build all CPU targets (analysis + seeds + sim + tools + runner)"
	@echo "  make install    Build all and copy binaries to ../bin/"
	@echo "  make analysis   Build analysis tools (10 tools in analysis/)"
	@echo "  make seeds      Build seed generators (5 generators)"
	@echo "  make sim        Build CPU simulation kernels"
	@echo "  make tools      Build visualization tools (sfa_render)"
	@echo "  make runner     Build the MCP simulation runner (Go)"
	@echo "  make cuda       Build CUDA simulation kernels (requires nvcc)"
	@echo "  make viewer     Build volume viewer (requires SDL2)"
	@echo "  make clean      Remove all binaries"
	@echo ""
	@echo "Individual:  make analysis/sfa_info  make seed/gen_braid  make runner  etc."
