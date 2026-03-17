# V23-F: Critical Gravity Phase 2 Redo (Equilibrated Profiles)

## Thesis

V23-C Phase 2 was inconclusive because the Gaussian initial conditions
caused violent shedding that overwhelmed the inter-oscillon interaction.
This investigation redoes Phase 2 properly: first equilibrate single
oscillons, save their profiles, then use those profiles as initial
conditions for two-oscillon runs.

## Background

V23-C Phase 1 found:
- All oscillons survive μ ∈ [-20, -5] with gap margins 5-12%
- Correlation length ξ ≈ 5-11 (does not diverge)
- κ_tail ≈ 0.1-0.2 (roughly constant)

V23-C Phase 2 found:
- Shedding transients dominate the dynamics
- F_avg measurements unreliable (non-parabolic separation vs time)
- Two merger events at mu=-20/D=50 and mu=-12/D=20, rest is repulsion/separation
- Need pre-equilibrated profiles to isolate the true interaction

## Method

### Step 1: Equilibrate Single Oscillons

For each μ ∈ {-20, -14, -10, -8, -6}:

1. Initialize a single symmetric oscillon (Gaussian, A=1.0, σ=3.0).
2. Evolve for t_equil = 10000 with Nx=4000, xmax=80.
3. At t=10000: save the FULL profile {φ₁(x), φ₂(x), φ₃(x), v₁(x), v₂(x), v₃(x)}
   (both field and velocity at a specific phase of the breathing cycle).
4. Also save the envelope profile and the breathing frequency ω.
5. Record the peak amplitude, energy, and core fraction.

### Step 2: Two-Oscillon Runs with Equilibrated Profiles

For each μ from Step 1, at separations D = {15, 20, 25, 30, 40, 50}:

1. Initialize two oscillons by superposing the saved profile:
     φ_a(x) = f_eq(x - D/2) + f_eq(x + D/2)
     v_a(x) = v_eq(x - D/2) + v_eq(x + D/2)
   (Both in-phase — same breathing phase.)

2. Evolve for t_run = 5000 with Nx=8000, xmax=200.

3. Track the center-of-energy of each oscillon vs time.

4. Measure:
   a. Initial acceleration a₀ = d²(sep)/dt² at t ~ 0-100
      (before any significant displacement)
   b. Separation vs time curve
   c. Whether oscillons attract, repel, or stay put
   d. Energy conservation (should be much better than Gaussian init)

### Step 3: Extract Interaction Potential

1. From the initial accelerations a₀(D) at each separation:
     F(D) = (M/2)·a₀(D)

2. Fit F(D) to candidate forms:
   a. Yukawa: F = F₀·exp(-D/λ)
   b. Power law: F = F₀/D^n
   c. Yukawa + constant: F = F₀·exp(-D/λ) + F₁

3. Extract the interaction range λ. Compare with ξ from V23-C Phase 1.

4. Integrate F(D) to get V(D). Check for equilibrium (V'(d_eq) = 0).

### Step 4: μ Dependence of Interaction Range

1. Plot λ vs μ. Does the range increase as |μ| decreases (gap closes)?
2. Plot λ vs ξ_pred = 1/√(m²-ω²). Is λ ∝ ξ_pred?
3. This directly tests whether the critical gravity mechanism works:
   longer ξ → longer interaction range → toward 1/r at criticality.

## Reference Code

- v23/critical: `/home/d/code/scp/v23/critical/src/critical1d.c`
  (Phase 2 code — modify to load profiles instead of Gaussian init)
- v21 1D solver: `/home/d/code/scp/v21/src/triad1d.c` (base code)

## Output Structure

- `src/critical2.c` — equilibration + two-oscillon interaction solver
- `data/profile_mu{N}.dat` — equilibrated profiles
- `data/interact_mu{N}_D{D}_ts.tsv` — separation vs time
- `RESULTS.md` — results and analysis

## Parameters

Equilibration: μ ∈ {-20, -14, -10, -8, -6}, κ=20, m=1.0
  Nx=4000, xmax=80, t_equil=10000
Interaction: D ∈ {15, 20, 25, 30, 40, 50}
  Nx=8000, xmax=200, t_run=5000

Compile: `gcc -O3 -Wall -o critical2 src/critical2.c -lm`
