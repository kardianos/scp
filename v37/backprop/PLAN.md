# Backprop Field Search: Full Operational Plan

## Objective
Find a self-sustaining, localized field configuration by training the raw field
values (not parameters) from a random spherical seed. Like training an LLM from
random weights — the field IS the model.

## Tools to Build

### 1. `field_search.py` — Main training loop (JAX/GPU)
- Random sphere initialization (smoothed noise × spherical envelope)
- PDE forward pass with absorbing BC in every timestep
- Curriculum training: T_opt ramps from 5 → 20 → 50 → 100 steps
- Multiple independent seeds (16-32) trained in parallel via vmap
- Loss: -|E_pot(T)| + λ_compact * R_rms² + λ_smooth * grad_energy
- Adam optimizer on the raw field values
- Checkpointing: save best candidates as .npy every 100 epochs
- Output: training log (TSV), best field snapshots

### 2. `validate_candidate.py` — Validation at higher resolution (JAX/GPU)
- Load N=32 candidate, upscale to N=48 or N=64 via trilinear interpolation
- Run forward for T=200 with absorbing BC
- Measure E_pot retention, R_rms, cluster count at checkpoints
- Output: validation report (pass/fail + metrics)

### 3. `candidate_to_sfa.py` — Convert candidate to SFA for C simulation (CPU)
- Load .npy candidate, upscale to N=128
- Write initial condition as first frame of an SFA file
- OR: write a raw binary that v37_knot.c can load as initial condition

### 4. `run_candidate_sim.sh` — Full C simulation of a candidate
- Takes a candidate .npy, converts to N=128 initial condition
- Runs v37_knot.c (or new v37_fieldseed.c) in full C simulation
- Outputs SFA file for visualization
- Runs sfa_frag + sfa_structure analysis

### 5. `download_best.sh` — SFA download pipeline
- While next experiment runs, download SFA from previous best
- Compress with zstd if needed

## Validation Tests (run before full search)

### V1: Gradient check
- Perturb a single grid point by ε, compute loss change
- Compare with JAX autodiff gradient at that point
- Must agree to 4+ digits

### V2: Energy conservation (no absorbing BC)
- Run forward 100 steps with no damping
- Total energy must be conserved to <1% drift
- Confirms the JAX PDE matches the C simulation

### V3: Absorbing BC works
- Initialize a Gaussian pulse at center
- Run forward with absorbing BC
- Energy must decrease monotonically (waves absorbed)
- Compare: same run WITHOUT absorbing should conserve energy

### V4: Known-good seed
- Load braid3(z) seed, run forward 100 steps
- E_pot should remain negative (known to work short-term)
- Gradient should be nonzero (optimizer can improve it)

## Execution Steps

### Phase 0: Build and upload (local)
1. Write field_search.py with all components
2. Write validate_candidate.py
3. Write validation tests as functions in field_search.py
4. Spin up V100 GPU instance on Vast.ai
5. Upload all scripts

### Phase 1: Validate (GPU, ~5 min)
6. Run V1-V4 validation tests
7. Debug any failures
8. Confirm: gradients correct, energy conserved, absorbing works, braid seed works

### Phase 2: Initial search (GPU, ~1 hour)
9. Run 32 random seeds × 1500 epochs with curriculum
10. Monitor: which seeds develop negative E_pot?
11. Save top 5 candidates

### Phase 3: Validate winners (GPU, ~30 min)
12. Run validate_candidate.py on top 5 at N=48
13. Which survive T=200 with absorbing BC?
14. Download SFA of any survivors

### Phase 4: Iterate (GPU, up to 8 hours)
15. If Phase 2 found stable candidates: increase grid to N=48, train longer
16. If Phase 2 found nothing: adjust loss function, try different random scales
17. Run multiple rounds, each time refining the approach
18. Every time a stable candidate is found:
    a. Generate SFA (run JAX sim, save f32 frames)
    b. Download SFA while next experiment runs
    c. Run sfa_frag locally to verify

### Phase 5: Final validation (local CPU)
19. Best candidate: upscale to N=128, run full C simulation T=500
20. sfa_frag + sfa_structure analysis
21. Compare with braid3 reference

## Agent Roles

### Builder Agent
- Writes code, deploys to GPU, runs experiments
- Optimistic: tries things, iterates fast

### Skeptic Agent
- Reviews builder's code for correctness
- Checks: are gradients right? Is the loss measuring what we think?
- Distrusts results: "is this real or an artifact?"
- Validates independently

## Budget
- V100 at ~$0.30-0.50/hr × 8 hours = ~$2.40-4.00
- Well within the Vast.ai credit balance

## Success Criteria
- A field configuration at N≥64 that:
  - Starts from random (not braid-seeded)
  - Survives T≥200 with absorbing BC
  - Maintains ≥30% of initial E_pot
  - Stays localized (R_rms < 0.5 × L)
  - Remains 1 cluster throughout (sfa_frag verified)
