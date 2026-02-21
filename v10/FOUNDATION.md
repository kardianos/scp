# V10: Two-Potential TOV Geon — Does Decoupling g_tt and g_rr Resolve Q9?

## Motivation: The V9 Redshift Paradox (Q9)

V9 proved that the f_2(1270) strong gravity geon produces hadron-scale
objects (3-quark geon at 887 MeV = 94.5% M_p). However, Q9 found a
**fundamental tension**: either the well is deep enough for trapping
(but E ~ 45 MeV due to redshift), or the energy matches M_p (but the
well is too shallow for a barrier).

V9 used a single metric function f(r):
```
ds² = -f(r) dt² + dr²/f(r) + r² dΩ²
```
This forces g_tt · g_rr = -1 (Schwarzschild gauge). As f → 0, BOTH
time slows down AND space stretches. The mode's frequency redshifts
as √f while the mode stretches as 1/f, causing catastrophic energy loss.

## V10 Hypothesis

Inside an extended matter distribution, g_tt and g_rr are independent:
```
ds² = -A(r) dt² + dr²/C(r) + r² dΩ²
```

In the Fierz-Pauli massive spin-2 theory (for EM sources with T = 0):
- A(r): sourced by energy density ρ (always positive → attractive)
- C(r): sourced by radial pressure p_r (changes sign → can be repulsive!)

For EM standing waves, p_r = ρ_E - ρ_B:
- Near origin (magnetic dominated): p_r < 0 → C > 1 (space compresses)
- At mode peak (electric dominated): p_r > 0 → C < 1 (space stretches)

**Hypothesis**: If C stays close to 1 while A drops (deep temporal well),
the mode energy is preserved because the effective diffusion coefficient
P = A·C dies more slowly than f² = A².

## Equations

### Two Yukawa BVPs
```
M_t'' - μ² M_t = -κ r ρ        (temporal, from energy density)
M_s'' - μ² M_s = -κ r ρ_s      (spatial, from mixed source)
```

where ρ_s = α_s·ρ + (1-α_s)·p_r:
- α_s = 1: V9 recovery (A = C, single potential)
- α_s = 0: Fierz-Pauli (spatial sourced by radial pressure)

### EM stress-energy in two-potential metric
```
ρ   = N_l/(8πr²) [ω²u²/A + C(u')²]     (energy density)
p_r = N_l/(8πr²) [ω²u²/A - C(u')²]     (radial pressure)
```

### Wave equation
The Regge-Wheeler equation for EM TM_l modes on `ds² = -A dt² + B dr² + r² dΩ²`
with B = 1/C gives the symmetrized eigenvalue problem:
```
-P v'' + W v = ω² v
```
where:
- P = A·C (diffusion coefficient; replaces f² in V9)
- v = (AC)^{1/4} · u (symmetrized variable; replaces √f·u in V9)
- W = A·l(l+1)/r² + P(β'/4 + β²/16) (effective potential)
- β = (ln P)' = A'/A + C'/C

**Verification**: At α_s = 1 (A = C = f), this reduces exactly to V9's
equation -f²v'' + [f·l(l+1)/r² + f·f''/2 - f'²/4]v = ω²v.

## Numerical Results

### Q1: V9 Recovery (α_s = 1) — CONFIRMED

At κ=0.06, μ=6.47, Rmax=0.5, l=1:
```
V10 (α_s=1): ω² = 47.76, A_min = C_min = 0.473, well = 0.527
V9 result:   ω² = 47.78, f_min = 0.473, well = 0.527
```
Agreement to 0.04%. Both mass functions identical: M_t = M_s.

### Q2: Fierz-Pauli decoupling (α_s = 0) — MODEST IMPROVEMENT

**Padé metric (no floor), κ=0.06, Rmax=0.5:**

| α_s | ω²    | A_min | C_min | C_max | E (MeV) | Comment           |
|-----|-------|-------|-------|-------|---------|-------------------|
| 1.0 | 57.8  | 0.687 | 0.687 | 1.000 | 1499    | V9 recovery       |
| 0.5 | 59.3  | 0.671 | 0.794 | 1.000 | 1520    | moderate decouple |
| 0.0 | 61.3  | 0.643 | 0.881 | 1.057 | 1544    | Fierz-Pauli       |

**Key findings:**
1. ω² increases by 6% from α_s=1 to α_s=0 (less redshift)
2. C stays closer to 1 at α_s=0 (less spatial curvature)
3. C_max = 1.057 at α_s=0: **spatial repulsion confirmed** near origin
4. Effect is quantitative (~6%), not qualitative

### Q3: Spatial repulsion physics — CONFIRMED

At Rmax=3.0, κ=35, μ=6.47, l=1, α_s=0 (subcritical):
```
M_s_min = -0.010  (REPULSIVE spatial mass function)
C at origin = 1.071  (7% spatial compression at deepest temporal well)
A_min = 0.929  (7% temporal well)
```

The spatial mass function IS negative near the origin, creating C > 1.
Physics: magnetic pressure (p_r < 0) from the EM standing wave pushes
space inward, partially compensating the temporal gravitational well.

### Q4: Phase transition — FIERZ-PAULI MAKES COLLAPSE WORSE

**Linear metric (floor), α_s scan at κ=0.06:**

| α_s | Status      | ω²    | A_min | C_min | E (MeV) |
|-----|------------|-------|-------|-------|---------|
| 1.0 | subcritical | 47.8  | 0.473 | 0.473 | 1363    |
| 0.65| subcritical | 45.3  | 0.321 | 0.514 | 1327    |
| 0.60| COLLAPSED  | 0.145 | 0.010 | 0.010 | 75      |
| 0.0 | COLLAPSED  | 0.170 | 0.010 | 0.010 | 81      |

**Sharp transition at α_s ≈ 0.62.** Below this, both potentials collapse.
The Fierz-Pauli decoupling makes collapse EASIER because the temporal
potential A drops faster without C holding it up.

### Q5: Q9 Test — REDSHIFT PARADOX PERSISTS

**Rmax scan at physical f_2 parameters (κ=35, μ=6.47, Padé):**

| Rmax | α_s=1 E (MeV) | α_s=1 A_min | α_s=0 E (MeV) | α_s=0 A_min | α_s=0 C_max |
|------|---------------|-------------|---------------|-------------|-------------|
| 1.50 | 80 (collapsed) | 0.024      | 4 (collapsed)  | 0.000       | 100 (guard) |
| 2.00 | 379           | 0.669       | 4 (collapsed)  | 0.000       | 100 (guard) |
| 2.25 | —             | —           | 367            | 0.608       | 2.23        |
| 2.50 | 337           | 0.866       | 340            | 0.810       | 1.24        |
| 4.00 | 220           | 0.985       | 221            | 0.984       | 1.02        |

**Critical observations:**
1. α_s=0 collapse extends to LARGER Rmax (2.0 vs 1.5 for V9)
2. Above transition, energies are nearly identical (±1%)
3. The Q9 sharp transition persists at α_s=0
4. No Rmax gives both deep well AND proton-mass energy

### Q6: 3-Quark Geon — FIERZ-PAULI COLLAPSES PROTON

**3×l=1 modes at κ=35, μ=6.47, Padé:**

| Rmax | α_s=1 E (MeV) | α_s=0 E (MeV) | α_s=0 C_max |
|------|---------------|---------------|-------------|
| 2.70 | 886 (94.5% M_p) | 0 (collapsed) | —          |
| 2.85 | —             | 937* (99.9% M_p) | 100 (guard) |
| 2.90 | —             | 864           | 2.11        |
| 3.00 | —             | 843           | 1.50        |

(*) Rmax=2.85 at α_s=0: A_min = 0.0001 (at Padé singularity). Solution
is numerically marginal — right at the collapse boundary. The E = 937 MeV
value is a coincidence of the sharp transition passing through M_p.

**The V9 proton-mass matching (886 MeV at Rmax=2.7) COLLAPSES at α_s=0.**
The critical Rmax shifts from 2.65→2.85. Above transition, energies
are comparable but the transition is equally sharp.

## Why the Two-Potential Approach Fails

### The 1/A runaway mechanism

The EM energy density contains the term ρ_E = ω²u²/A. As the temporal
well deepens (A → 0), ρ_E diverges. This creates a self-reinforcing
positive feedback loop:

```
deeper A → larger ρ_E → larger ∫ρ → larger M_t → deeper A
```

This loop is INDEPENDENT of the spatial curvature C. Whether C stays
close to 1 (Fierz-Pauli) or tracks A (V9), the temporal potential
collapses by its own feedback.

### Why spatial repulsion doesn't help

At moderate coupling, C > 1 near the origin (spatial repulsion from
magnetic pressure). But the magnitude is small:
- A drops by ~35% (well_A = 0.35 at κ=0.06, Padé, α_s=0)
- C rises by only 6% (C_max = 1.057)

This mismatch occurs because ρ > |p_r| always (ρ = ρ_E + ρ_B while
|p_r| = |ρ_E - ρ_B| < ρ). The temporal source is always stronger than
the spatial source.

### The P = A·C argument

In V9: P = f² → as f → 0, P → 0 as f².
In V10: P = A·C → if C compensated (C ~ 1/A), then P ~ 1 (constant).

But C does NOT compensate by 1/A. At best, C rises by ~6-10% while A
drops by ~35%. So P ≈ 0.65 × 1.06 ≈ 0.69, not much better than
P = 0.65² = 0.42 (V9). The improvement exists but is insufficient
for the Q9 scenario where A → 0.

## Conclusion

**V10 is a DEFINITIVE NEGATIVE RESULT for the hypothesis that decoupling
g_tt and g_rr resolves the V9 redshift paradox (Q9).**

The Fierz-Pauli massive spin-2 theory does produce:
1. Spatial repulsion near the origin (C > 1, up to 50% at transition)
2. A modest increase in mode energy (~6%)
3. Less spatial curvature (C stays closer to 1 than A)

But it does NOT:
1. Resolve the Q9 tension (sharp transition persists)
2. Prevent temporal collapse (1/A runaway is self-reinforcing)
3. Improve the phase transition structure (actually makes it worse)

**The fundamental issue is that ρ_EM ∝ 1/g_tt**, which drives the temporal
potential into self-reinforcing collapse regardless of the spatial metric.
The spatial curvature is a passive response, too weak to compensate.

### Implications for Options B and C

Per the V9 NEXT.md analysis:

**Option B (Nonlinear Topological Geon)**: The V10 result strengthens the
case for nonlinear EM. The 1/A runaway occurs because ρ ∝ 1/A in LINEAR
Maxwell theory. In Born-Infeld or other nonlinear theories, the energy
density is capped (ρ_BI ~ b² √(1 + ρ_M/b²)), which could prevent the
runaway. Combined with topological stabilization, this could break the
Q9 tension.

**Option C (Bag Model)**: The box IS the stabilizer in V9/V10. Making the
box dynamical (via a confining scalar field) would make Rmax a self-
consistent variable. But V9 Q9 showed that no Rmax gives both barrier
AND proton mass, so a bag model would face the same tension unless the
bag energy itself contributes to the mass.

## Code

- `src/tov_geon.c`: Two-potential Yukawa geon Hartree solver
  - `-hartree`: Two-potential Hartree solve
  - `-scan`: Coupling constant scan
  - `-alphascan`: α_s scan (0=Fierz-Pauli, 1=V9)
  - `-rmaxscan`: Rmax scan (Q9 test)
  - `-lifetime`: WKB lifetime
  - `-multimode`: Multi-mode Hartree (3-quark geon)
  - `-alpha <A>`: Spatial source mixing parameter
  - `-nlmetric N`: 0=linear, 1=exp, 2=Padé
  - Other parameters same as V9
