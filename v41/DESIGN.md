# V41: First-Principles Seed Construction from Stability Signatures

## Goal

Construct initial field configurations that satisfy the three stability signatures
discovered in V40 Gen 4 analysis, then predict their stability BEFORE simulation,
and run only the most promising candidates.

## The Three Stability Signatures (from V40 STRUCTURAL_FINDINGS.md)

### Signature 1: θ Confinement (outer/inner < 0.7)
- θ must be concentrated in the core, not radiating outward
- Unstable clusters have θ increasing outward (outer/inner > 1.0)
- **Construction approach**: Initialize θ with a Gaussian envelope that matches
  the φ core, NOT at zero. Pre-loading θ at its equilibrium pattern prevents
  the initial transient where θ grows from zero and radiates.

### Signature 2: Velocity Structure
- Single braids: |v| increases outward (breathing shell)
- 3-braid composites: |v| DECREASES outward (contraction)
- **Construction approach**: Initialize velocities with a radial profile that
  matches the target structure. For 3-braids, set v·r̂ < 0 (inward) in
  the outer shell.

### Signature 3: |P| Concentration (inner/outer > 10×)
- Binding concentrated in core, nearly zero at edges
- Critical P value ≈ 0.082 (maximum binding force from V40 RESEARCH.md)
- **Construction approach**: Tune per-braid amplitudes so that the OVERLAP region
  (where all three braids contribute to P) has P_peak near 0.082.

## Seed Construction Method

### The field structure

Three perpendicular braids (along x, y, z axes) with chirality assignment:
- UUD: braids 1,2 have chirality +1 (Up), braid 3 has chirality -1 (Down)
- UDD: braid 1 has chirality +1 (Up), braids 2,3 have chirality -1 (Down)

Each braid contributes to all three field components φ₀, φ₁, φ₂.

### Parameters to sweep (36 seeds per chirality)

| Parameter | Values | Count | Description |
|-----------|--------|-------|-------------|
| A (amplitude) | 0.3, 0.5, 0.7 | 3 | Per-braid amplitude |
| R_tube | 2.5, 3.5, 4.5 | 3 | Tube radius |
| theta_init | 0, 0.5, 1.0 | 3 (relative) | Pre-initialized θ fraction |
| v_profile | breathing, contracting, zero, mixed | 4 | Velocity initialization |

3 × 3 × 4 = 36 combinations per chirality. Total: 72 seeds.

### Pre-simulation stability prediction

For each seed, BEFORE running the simulation, compute:
1. **P_peak** in the overlap region → should be near 0.082
2. **θ outer/inner** ratio → should be < 0.7
3. **v outer/inner** ratio → should be > 1.3 (breathing) or < 0.7 (contracting)
4. **ρ inner/outer** ratio → should be > 5
5. **E_pot total** → more negative = stronger initial binding
6. **|P| concentration** = P_int(r<R_half) / P_int(total) → should be > 0.8
7. **Force balance** at core: |Σforces| / |largest_force| → should be < 0.5

Combined stability prediction score:
```
S_pred = w1*f(P_peak) + w2*f(θ_ratio) + w3*f(v_ratio) + w4*f(ρ_ratio) + w5*f(E_pot) + w6*f(P_conc) + w7*f(force_bal)
```

### Ranking and selection

1. Generate 36 UUD seeds + 36 UDD seeds = 72 total
2. Compute S_pred for each
3. Rank within each chirality
4. Select top 3 UUD + top 3 UDD = 6 candidates
5. Run simulations at N=192, L=25, T=200 on Tesla_V100

## Implementation

### Tool: `v41/construct_seed.c`

A single C program that:
1. Takes command-line parameters (A, R_tube, theta_init, v_profile, chirality)
2. Constructs the 3-braid field with pre-loaded θ and shaped velocities
3. Computes the 7 stability metrics on the constructed field
4. Writes the seed SFA and outputs the metrics to stdout

### Tool: `v41/sweep.sh`

Shell script that:
1. Compiles construct_seed.c
2. Runs all 72 parameter combinations
3. Collects metrics into a TSV
4. Ranks and selects top 3 per chirality
5. Copies the 6 winning seeds to a run directory

### Simulation: Run on Tesla_V100

Upload 6 seeds + scp_sim.cu, run T=200, download f16 SFAs + analysis JSONs.
