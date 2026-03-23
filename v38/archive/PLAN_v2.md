# V38 Experiment Plan v2 (Corrected)

## Items Removed (Known Wrong or Premature)
- ~~Braid = electron identification~~ (size hierarchy fails: braid 10× Bohr radius)
- ~~Sub-braid quark speculation~~ (no mechanism, unfalsifiable at current resolution)
- ~~Electron-based recalibration~~ (doesn't fix ℏ problem, creates size problem)
- ~~Fine structure constant calculation~~ (eta is free parameter, no prediction possible)

## What IS Established
1. The 6-field Cosserat system (V34+) produces CHARGE-DEPENDENT FORCES
2. Winding number is topologically conserved (integer, robust)
3. The theta field is a COMPONENT of EMF, not the full photon field
4. The braid self-reconstructs through motion (V26)
5. The braid only exists in periodic BC — it is NOT compact
6. ℏ_code = 38.6 (scale-invariant), ℏ_eff = 22,727 (590× too large)
7. V28 CMA-ES parameters (m²=2.25, μ=-41.345, κ=50) were optimized for
   braid stability, NOT for matching any physical mass

## Open Questions (What We Don't Know)

### Q1: What physical scale do the V33+ parameters correspond to?
The V2 Skyrme parameters were fitted to the proton. The V33+ Cosserat parameters
were found by CMA-ES. What mass does the braid NATURALLY correspond to?
- The braid energy is 5350 code units
- The mass gap m=1.5 sets a natural energy scale
- WITHOUT external calibration, what does 1 code unit of energy map to?
- Could the V28 CMA-ES have converged on parameters that correspond to a
  completely different energy scale than the Skyrme model?

### Q2: Why does the braid self-reconstruct and what is the minimum for survival?
V26 Discovery 1: the traveling wave continuously remakes V(P). But:
- What is the minimum propagation speed for self-reconstruction?
- Is there a minimum amplitude? Minimum tube radius?
- If stopped (vel=0), how quickly does it die?
- The MECHANISM is the key — once understood, we can design compact variants

### Q3: What is the theta field's actual radial structure?
V34 showed charge-dependent force, but:
- What is the EXACT radial profile θ(r) around the braid?
- Is it 1/r (Coulomb-like), r^(-0.5) (RESEARCH found this), or something else?
- How does it depend on braid speed?
- The theta field is a COMPONENT of EMF — which component exactly?

### Q4: What is the 590× ℏ discrepancy?
ℏ_eff = E_braid/ω = 22,727 vs ℏ_physical = 38.6.
- Is ℏ_eff the wrong definition? (Maybe ℏ ≠ E/ω for this system)
- Does the braid have internal modes whose action IS 38.6?
- Is the braid a MULTI-quantum object (containing ~590 quanta)?
- This is the most quantitative test available

## Phase 1: Characterize the Braid Theta Field (GPU, 30 min)

**Goal**: Measure θ(r) precisely around a propagating braid at high resolution.

**Method**:
- Run braid at N=256 on V100 using CUDA (`braid_analyze_cuda.cu`)
- Every 10 steps: compute azimuthally-averaged θ(r) at braid center
- Fit to power law: θ ~ r^α. Is α = -1 (Coulomb)? -0.5? Other?
- Measure θ vs braid speed (run at multiple velocities)
- Compare with Biot-Savart: B = μ₀I/(2πr) → θ should ~ 1/r if EM-like

**Output**: `data/braid_theta_profile.tsv`, SFA file for visualization

## Phase 2: Self-Reconstruction Threshold (GPU, 1 hr)

**Goal**: Find the minimum conditions for braid survival.

**Method**:
- Start braid at various speeds: v = 0, 0.1, 0.3, 0.5, 0.7, 0.9
- For each: does the braid survive T=200? At what speed does it die?
- Start braid with various amplitudes: A = 0.2, 0.4, 0.6, 0.8, 1.0
- For each: survival threshold?
- Start braid and STOP it (set vel=0 at t=50): how fast does it die?

**Output**: Phase diagram of (speed, amplitude) → survival

## Phase 3: Multi-Braid Interaction (GPU, 2 hr)

**Goal**: Observe two braids interacting with theta coupling.

**Method**:
- Two parallel braids at separation D=5, 10, 15, 20
- Both same winding (UU) and opposite winding (UD)
- Run T=300 at N=256
- Measure: force vs D, theta field between braids, energy transfer

**Output**: Force law F(D), comparison with 1/r² Coulomb

## Phase 4: ℏ Investigation (analytical + short GPU runs, 1 hr)

**Goal**: Understand the 590× discrepancy.

**Method**:
- Compute the braid's internal mode spectrum (breathing, wobble, etc.)
- Each mode has energy E_n and frequency ω_n
- Compute E_n/ω_n for each mode — is any mode's action = 38.6?
- Count the number of modes with E_n/ω_n ~ 38.6 — this is the braid's "quantum number"
- If the braid contains ~590 quanta of action 38.6, it IS a multi-quantum object

## Budget
- V100 GPU: ~4 hours total (~$0.50)
- Must use CUDA, NOT CPU (see CLAUDE.md GPU policy)
- Monitor GPU every 10s with monitor.sh
