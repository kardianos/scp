# V24-S1b: Asymmetric Oscillon Proca Force — RESULTS

## Parameters

mu=-20, kappa=20, mass=1.0, lambda=0.99, A=0.8, sigma=3.0
Proca mass: m_A = sqrt(1 - 0.99) = 0.1, range 1/m_A = 10.0
Equilibration: Nx=4000, xmax=80, t=10000
Interaction: Nx=8000, xmax=300, t=2000

## Phase 1: Equilibration Summary

| Option | omega  | E_total  | A_peak | S      | A2      | alive |
|--------|--------|----------|--------|--------|---------|-------|
| SYM    | 1.6020 | +11.72   | 0.751  | 2.977  | 0.000   | YES   |
| 180    | 0.001  | -54.70   | 2.164  | 0.533  | 26.639  | YES   |
| UUD    | ---    | NaN      | 0.000  | ---    | ---     | NO    |
| SEED   | 1.6020 | +11.55   | 0.749  | 2.953  | 0.001   | YES   |

### Key observations:

1. **180-degree oscillon** survived equilibration with MASSIVE antisymmetric content
   (A2/S = 50). The anti-phase structure is a stable attractor — phi3 grows to -2.2
   while phi1=phi2 stay at ~1.1. Energy is deeply negative (E=-54.7) due to large
   P = phi1*phi2*phi3 < 0 contribution. Frequency ~0 (not oscillating in the usual
   sense; the fields are large and quasi-static with slow modulation).

2. **UUD option** (m3=0.95, all same sign) blew up to NaN. The pairwise coupling at
   lambda=0.99 makes m2_anti = 0.01 near-tachyonic. With all-positive initial data
   and slightly different m3, the configuration is unstable. DEAD.

3. **SEED option** (phi3 scaled by 0.9): the 10% asymmetry was completely radiated
   away during equilibration. By t=3000, phi1=phi2=phi3 to 4 decimal places.
   Final A2/S = 0.0004. The symmetric oscillon is a strong attractor —
   antisymmetric perturbations are expelled, consistent with the Proca mass gap.

## Phase 2: Force Table

| D  | F(SYM)      | F(180)      | F(SEED)     |
|----|-------------|-------------|-------------|
| 20 | +2.14e-06   | -2.15e-03   | -2.36e-06   |
| 30 | +8.15e-06   | -1.96e-03   | +5.34e-06   |
| 40 | -1.71e-06   | -1.79e-03   | +5.91e-07   |
| 60 | -1.77e-05   | -1.44e-03   | -1.00e-05   |
| 80 | +5.42e-06   | -1.13e-03   | +2.59e-06   |

### Yukawa fits:

| Option | F0       | xi_fit   | 1/m_A | sign     |
|--------|----------|----------|-------|----------|
| SYM    | 2.2e-06  | -56.5    | 10.0  | mixed    |
| 180    | 2.7e-03  | 93.8     | 10.0  | attract  |
| SEED   | 2.0e-06  | -124.2   | 10.0  | mixed    |

## Phase 3: Key Comparison

**Ratio |F(asym)/F(sym)| at each D:**

| D  | 180     | SEED   |
|----|---------|--------|
| 20 | 1007    | 1.10   |
| 30 | 240     | 0.65   |
| 40 | 1042    | 0.35   |
| 60 | 82      | 0.57   |
| 80 | 209     | 0.48   |

## CRITICAL CAVEAT: The 180 result is NOT a clean Proca measurement

The 180-degree oscillon is a fundamentally different beast:
- Energy E=-54.7 (vs +11.7 for SYM) — 5x larger magnitude, opposite sign
- Peak amplitude 2.16 (vs 0.75) — nearly 3x larger
- The "force" measured is 100-1000x larger, but the initial separation measurement
  is wrong: at D=20, the energy centroid starts at sep=76 (not 20), meaning the
  oscillons are so extended that the centroids can't be cleanly separated
- The Yukawa fit gives xi=94, not 10. This is NOT Proca exchange; it's the massive
  nonlinear overlap of two very large structures

**The 180 force is dominated by direct field overlap, not Proca exchange.**

The SEED option cleanly tests the Proca hypothesis: start with a small asymmetry,
see if the Proca channel carries any long-range force. Answer: **NO**. The asymmetry
is radiated away during equilibration (A2/S drops from 0.05 to 0.0004), and the
resulting force is indistinguishable from the symmetric control.

## Verdict

**Proca channel NOT detected as a force mediator between oscillons.**

The result is actually deeper than "not detected" — it's structurally impossible:

1. **Asymmetric perturbations are unstable**: The Proca mass gap (m_A = 0.1) means
   antisymmetric fluctuations around the symmetric oscillon are massive. Any small
   A2 content is expelled on timescale ~1/m_A = 10. The SEED option confirms this.

2. **The 180-degree oscillon** retains asymmetry but it's a completely different
   soliton branch (negative energy, non-oscillating). Its force is nonlinear overlap,
   not Proca exchange. The fitted range (94) bears no relation to 1/m_A = 10.

3. **For Proca exchange to matter**, we would need:
   - A stable oscillon with persistent A2 content
   - The A2 tail to decay as exp(-m_A*r) (range 10), much longer than the
     symmetric tail exp(-m_tail*r) (range ~2)
   - But the dynamics expel A2 content — there is no stable source

The symmetric oscillon is a perfect S-wave state; it simply does not couple to the
antisymmetric channel. This is a selection rule, not a sensitivity issue.
