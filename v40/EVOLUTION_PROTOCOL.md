# V40 Evolution Protocol — Binding Reference for All Agents

## Core Principle

This is an EVOLUTIONARY search, NOT a parameter sweep. Each generation builds
on the OUTPUTS of the previous generation. Seeds come from simulation results,
not from fresh initialization.

## The Evolutionary Loop

```
1. RUN promising configurations at long T (T=200+)
2. SELECT survivors (structures that maintain binding through T=200)
3. MODIFY the last frame of survivors (add sub-structures, perturbations)
4. RESUME simulation from the modified state
5. SCORE: did the modification improve or degrade stability?
6. SELECT the best modified structures → go to step 3
```

## What We Learned (Gen 0-3): Exploration Phase

Gen 0-3 were parameter sweeps testing FRESH initial configurations:
- Higher amplitude → deeper binding (but may not be stable long-term)
- Counter-chirality pairing → constructive V(P) interference
- UDD 3-braid → genuine 3D composite binding
- Larger grid (N=192, L=25) eliminates boundary artifacts
- Wider tube radius (R=4) gives better triple-product overlap

**Top configurations from exploration:**
1. CB15: Counter-braid 1.5× (S=3.87, E_pot=-697)
2. S20: Braid 2.0× (S=3.37, E_pot=-612)
3. UDD_R4: 3-braid UDD A=0.5 R=4 (S=3.16, E_pot=-537)

These scores are from T=30 SHORT PROBES. We do NOT know if they survive T=200.

## Gen 4: Stability Testing (Current Phase)

### Step 1: Long-term survival test

Run the top 3 from Gen 3 at T=200 on Tesla_V100 (N=192, L=25, f32 output):
- CB15: Counter-braid 1.5×
- S20: Braid 2.0×
- UDD_R4: 3-braid UDD A=0.5 R=4

Monitor: E_pot, P_int, theta_rms, cluster count over full T=200.
A structure SURVIVES if:
- P_int stays above 20% of initial at T=200
- E_pot remains negative (binding persists)
- phi_max stays above 2× background (not dispersed)

### Step 2: Take last frame of survivors

For each survivor, the T=200 output SFA's last frame becomes the seed for
the next step. Load it via `init=sfa, init_frame=-1`.

### Step 3: Modify and test

For each surviving last frame, generate 3-4 variants using modify_sfa:
- Add a small perturbation (--perturb 0.05)
- Add an oscillon near the core (--add-oscillon 0 0 0 0.3 2.0)
- Add a rotated braid component (--add-braid cx cy cz angle)
- Scale a sub-region (if tool supports it)

Resume each variant for T=50. Score: does S improve or degrade?

### Step 4: Select and iterate

Take the top 2 variants from step 3. Run them to T=200. Take their last frames.
Repeat steps 3-4 until convergence (S stops improving for 2 consecutive rounds).

## GPU Execution Protocol

### Remote setup
```bash
# Search for Tesla_V100
vastai search offers 'gpu_name=Tesla_V100 num_gpus=1 rentable=true disk_space>=20' -o 'dph' --limit 5

# Create instance
vastai create instance <ID> --image nvidia/cuda:12.2.0-devel-ubuntu22.04 --disk 20 --ssh

# Upload source + compile on remote (DO NOT generate seeds locally and upload)
scp scp_sim.cu sfa.h gen_3braid.c gen_braid.c modify_sfa.c analyze_sfa.c remote:~/
ssh remote "mkdir -p format && cp sfa.h format/ && \
  sed -i 's|#include.*sfa\.h.*|#include \"format/sfa.h\"|' *.c *.cu && \
  apt-get update -qq && apt-get install -y -qq libzstd-dev > /dev/null 2>&1 && \
  nvcc -O3 -arch=sm_70 -o scp_sim_cuda scp_sim.cu -lzstd -lm && \
  gcc -O3 -o gen_3braid gen_3braid.c -lzstd -lm && \
  gcc -O3 -o gen_braid gen_braid.c -lzstd -lm && \
  gcc -O3 -o modify_sfa modify_sfa.c -lzstd -lm && \
  gcc -O3 -fopenmp -o analyze_sfa analyze_sfa.c -lzstd -lm && \
  echo BUILD_OK"
```

### Generate seeds remotely
```bash
ssh remote "./gen_3braid -N 192 -L 25 -A 0.5 -R 4 -chirality UDD -o seed.sfa"
```

### Run simulation
```bash
ssh remote "./scp_sim_cuda config.cfg"
```

### Analyze
```bash
ssh remote "./analyze_sfa output.sfa --json analysis.json"
```

### Download results (f32 SFA for viewing, JSON for scoring)
```bash
scp remote:~/output.sfa ./local_results/
scp remote:~/analysis.json ./local_results/
```

### Resume from output (modify last frame, then continue)
```bash
ssh remote "./modify_sfa output.sfa --frame -1 --perturb 0.05 -o modified_seed.sfa"
ssh remote "./scp_sim_cuda modified_config.cfg"  # init_sfa = modified_seed.sfa
```

## Key Constraints

1. **DO NOT modify scp_sim.c, scp_sim.cu, or sfa.h** without explicit user approval
2. **DO NOT generate seeds locally and upload** — compile seed generators remotely
3. **ALL simulation MUST run on GPU, not CPU** (CLAUDE.md policy)
4. **Seeds come from SIMULATION OUTPUT**, not fresh generation (after Gen 3)
5. **Monitor GPU utilization** — if a sim takes too long, check nvidia-smi
6. **Download SFA files** so the user can view results in volview
7. **Keep the GPU instance alive** until ALL runs are verified
8. **Destroy instance** only after downloading all needed files

## Stability Scoring

```
S(t) = 0.4 × energy_score + 0.3 × coherence_score + 0.2 × theta_score + 0.1 × shape_score

energy_score = max(0, -E_pot / 80)         # normalized to single-braid E_pot
coherence_score = P_int(t) / P_int(0)      # triple product retention
theta_score = min(1, theta_rms / 0.05)     # theta equilibration
shape_score = 1 / (1 + aspect - 1)         # penalize extreme anisotropy
```

A structure is DEAD when:
- P_int < 10% of initial (coherence lost)
- E_pot > 0 (no binding)
- phi_max < 2× background (dispersed)
- Any of these sustained for T > 20

## File Locations

- Protocol: `/home/d/code/scp/v40/EVOLUTION_PROTOCOL.md` (THIS FILE)
- Postulates: `/home/d/code/scp/v40/POSTULATES.md`
- Plan: `/home/d/code/scp/v40/PLAN.md`
- Gen summaries: `/home/d/code/scp/v40/GEN{1,2,3}_SUMMARY.md`
- State: `/home/d/code/scp/v40/state/generation.json`
- Tools: `/home/d/code/scp/v40/tools/` (analyze_sfa, modify_sfa, gen_3braid, spatial_analysis)
- Sim kernel: `/home/d/code/scp/sfa/sim/scp_sim.cu` (DO NOT MODIFY)
- SFA header: `/home/d/code/scp/sfa/format/sfa.h` (DO NOT MODIFY)
- Seed generators: `/home/d/code/scp/sfa/seed/gen_braid.c` (DO NOT MODIFY)
- Project rules: `/home/d/code/scp/CLAUDE.md`
