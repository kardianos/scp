# v93/theory — progressive / shear probe suite (Go, no Python)
# Run from v93/:  make -C theory probe-run
# Or:             make -C theory probe && ./theory/bin/probe -j 6

V93     := $(abspath ..)
BIN     := bin
PROBE   := $(BIN)/probe
GO      ?= go
J       ?= 6
SEED    ?= 20260802

.PHONY: all probe probe-run probe-shear probe-prog clean

all: probe

probe: $(PROBE)

$(PROBE): probe/*.go
	@mkdir -p $(BIN)
	cd $(V93) && $(GO) build -o theory/$(PROBE) ./theory/probe

# full matrix (shear + progressive + controls)
probe-run: $(PROBE)
	cd $(V93) && ./theory/$(PROBE) -j $(J) -seed $(SEED)

# shear / field-organization arms only
probe-shear: $(PROBE)
	cd $(V93) && ./theory/$(PROBE) -j $(J) -seed $(SEED) -only shear,ctrl

# progressive dense / Bohmian-method arms only
probe-prog: $(PROBE)
	cd $(V93) && ./theory/$(PROBE) -j $(J) -seed $(SEED) -only prog

# long-horizon sustain isolation (default T=3000, Go fabrun)
T       ?= 3000
KERNEL  ?= fabrun
sustain: $(PROBE)
	cd $(V93) && GOMAXPROCS=4 ./theory/$(PROBE) -suite sustain -kernel $(KERNEL) \
		-T $(T) -j $(J) -seed $(SEED) -snap-every 5000

clean:
	rm -rf $(BIN)
