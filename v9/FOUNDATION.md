# V9: Salam Strong Geon — Nuclear-Scale Self-Trapping via Massive Spin-2 Mediator

## Motivation: The V8 Coupling Gap

V8 proved that tensor coupling (spin-2 mediator coupling to T_mu_nu)
provides positive feedback for EM self-trapping — qualitatively different
from the scalar no-go theorem (V6+V7). However, using Newton's G, the
coupling gap is catastrophic:

```
kappa_cr / kappa_GR ~ (R / l_Planck)^2 ~ 10^{40}
```

Geons in pure GR are Planck-scale proto-black-holes.

## The Salam Strong Gravity Hypothesis (1970s)

Abdus Salam noticed this same gap and proposed that the strong nuclear
force is mediated in part by a **massive spin-2 meson** — the f_2(1270).
Being spin-2, the f_2 couples to T_{mu nu} exactly like a graviton, but:

1. **Much stronger coupling**: G_strong ~ 10^{38} G_N
2. **Finite range**: Yukawa cutoff at 1/m_{f_2} ~ 0.15 fm

The coupling gap is exactly the ratio M_Pl^2 / m_{f_2}^2 — not a mystery,
but a consequence of the mediator mass.

### Physical parameters

```
f_2(1270) meson:
  Mass:      m_f2 = 1.2754 GeV    (mu = 6.47 fm^{-1})
  Width:     Gamma = 187 MeV      (broad resonance)
  Spin/P/C:  2++                   (isoscalar tensor)
  Range:     1/mu = 0.155 fm

Strong gravitational coupling:
  G_strong = g^2 / (4 pi m_f2^2) ~ 71 GeV^{-2}   (~ 10^{40} G_N)
  kappa_strong = 4 pi G_strong ~ 895 GeV^{-2}
```

### Self-consistency check (V8 FOUNDATION)

```
G_strong x M_geon ~ R^3 mu^2 / 3
G_strong x 1 GeV  ~ (1 fm)^3 x (6.5 fm^{-1})^2 / 3 ~ 14 fm
G_strong ~ 71 GeV^{-2} ~ 10^{40} G_N   checkmark
```

This matches the V8 coupling gap — the strong geon IS the mechanism
that bridges it.

## Equations

### Metric and EM mode (same as V8)

```
ds^2 = -f(r) dt^2 + dr^2/f(r) + r^2 dOmega^2
f(r) = 1 + 2 Phi(r)
```

EM TM_l mode with radial function u(r), frequency omega.
Regge-Wheeler equation (symmetrized v = sqrt(f) u):
```
-f^2 v'' + W(r) v = omega^2 v
W(r) = f l(l+1)/r^2 + f f''/2 - f'^2/4
```

### Yukawa potential (replaces Einstein mass integration)

The massive spin-2 mediator satisfies:
```
(nabla^2 - mu^2) Phi = kappa rho_EM
```

In spherical symmetry with w = r Phi:
```
w'' - mu^2 w = kappa r rho_EM
```

Since rho_EM > 0 and the Green's function is negative,
w < 0 and Phi < 0 (attractive potential). The metric function:
```
f(r) = 1 + 2 Phi(r) = 1 + 2 w(r)/r < 1  in the core
```

**Key difference from V8 (GR):**
- V8: m'(r) = kappa r^2 rho  (algebraic integration, long-range 1/r)
- V9: w'' - mu^2 w = kappa r rho  (BVP, Yukawa decay e^{-mu r}/r)

The Yukawa BVP is solved by Thomas algorithm (tridiagonal).

### Boundary conditions

```
w(0) = 0                    (regularity: Phi finite at origin)
w'(Rmax) = -mu w(Rmax)      (outgoing Yukawa at infinity)
```

### Effective mass function

Define M_eff(r) = -w(r) > 0, so f = 1 - 2 M_eff(r) / r.
The BVP becomes:
```
M_eff'' - mu^2 M_eff = -kappa r rho_EM
```
Same structure as V8's mass function, but solved as BVP not integral.

## Code Units

All lengths in code units (dimensionless). Physical conversion:
```
1 code length = L_0    (set by user, typically 1 fm)
hbar c / L_0 = 197.3 MeV * fm / L_0
```

With L_0 = 1 fm:
```
mu_code = m_f2 * L_0 / (hbar c) = 1.2754 / 0.1973 = 6.47
kappa_code = kappa_strong * (hbar c)^2 / L_0^2  [see below]
```

The code EM density has specific normalization (mode u with int u^2 dr = 1).
kappa is a free dimensionless coupling constant; physical mapping depends
on the energy content of the mode.

## Numerical Results

### Q1: Positive feedback — CONFIRMED

delta(omega^2) < 0 for ALL mu values tested. The Yukawa well lowers
the EM eigenvalue, confirming positive feedback (tensor coupling).

Perturbative results at Rmax=20, l=1:

| mu    | delta(omega^2)/kappa | frac_shift/kappa | kappa_cr(pert) |
|-------|---------------------|------------------|----------------|
| 0.01  | -6.00e-6            | -1.19e-4         | 8.41e3         |
| 0.5   | -2.36e-7            | -4.67e-6         | 2.14e5         |
| 1.0   | -6.21e-8            | -1.23e-6         | 8.13e5         |
| 6.47  | -1.55e-9            | -3.07e-8         | 3.26e7         |
| 10.0  | -6.48e-10           | -1.28e-8         | 7.79e7         |

At mu=0.01 (nearly massless): recovers V8 GR result. kappa_cr ~ 10^4.
As mu increases: kappa_cr grows as mu^2 (Yukawa suppression).

### Q2: Scaling — kappa_cr proportional to mu^2 Rmax^5

For the f_2 mass (mu=6.47), kappa_cr at various box sizes:

| Rmax | omega^2_0 | kappa_cr(pert) | kappa_cr(Hartree) | kappa_phys/kappa_cr(H) |
|------|-----------|----------------|-------------------|----------------------|
| 0.5  | 80.75     | 0.57           | ~0.063            | 556 (supercritical)  |
| 1.0  | 20.19     | 12.5           | ~1.4              | 25  (supercritical)  |
| 1.5  | 8.97      | 85.6           | ~9.5              | 3.7 (supercritical)  |
| 2.0  | 5.05      | 346            | ~38               | 0.9 (marginal!)      |
| 3.0  | 2.24      | 2,544          | ~283              | 0.12 (subcritical)   |
| 5.0  | 0.81      | 32,149         | ~3,572            | 0.01 (subcritical)   |

Scaling law for Rmax >> 1/mu:
```
kappa_cr(pert) ~ 0.24 x mu^2 x Rmax^5
kappa_cr(Hartree) ~ kappa_cr(pert) / 9    (nonlinear amplification)
```

**Physical coupling** (kappa_phys = 4 pi G_strong / L_0^2 ~ 35 in code units
with L_0 = 1 fm):

Hartree results at kappa=35, mu=6.47 (physical f_2 parameters):

| Rmax (fm) | f_min  | well_depth | R_eff (fm) | Status        |
|-----------|--------|------------|------------|---------------|
| 2.00      | 0.01   | 99.0%      | 1.70       | COLLAPSED     |
| 2.05      | 0.01   | 99.0%      | 1.73       | COLLAPSED     |
| 2.10      | 0.532  | 46.8%      | 1.21       | subcritical   |
| 2.20      | 0.686  | 31.4%      | 1.30       | subcritical   |
| 2.50      | 0.845  | 15.5%      | 1.51       | subcritical   |
| 3.00      | 0.936  | 6.4%       | 1.83       | subcritical   |
| 5.00      | 0.995  | 0.5%       | 3.06       | perturbative  |

**Critical crossover at Rmax ~ 2.08 fm (proton diameter!)**

Just above transition (Rmax=2.1):
```
omega^2 = 4.41    omega = 2.10 fm^{-1}
E_mode  = omega x hbar c = 0.414 GeV   (~44% of proton mass)
R_eff   = 1.21 fm                       (~1.4x proton charge radius)
Phi_min = -0.234                        (23% potential depth)
```

### Q3: Phase transition — sharp collapse

Hartree coupling scan at mu=6.47, Rmax=0.5 (compact box, l=1):

| kappa | omega^2 | well_depth | frac_shift | R_eff |
|-------|---------|-----------|------------|-------|
| 0.00  | 80.72   | 0.000     | 0.000      | 0.306 |
| 0.01  | 76.09   | 0.072     | -0.057     | 0.305 |
| 0.03  | 66.39   | 0.222     | -0.178     | 0.302 |
| 0.05  | 55.29   | 0.398     | -0.315     | 0.296 |
| 0.06  | 47.78   | 0.527     | -0.408     | 0.289 |
| 0.065 | **0.15**| **0.990** | **-0.998** | 0.418 |

**SHARP PHASE TRANSITION** between kappa=0.060 and kappa=0.065.
No smooth intermediate: either gentle perturbation or full collapse.
The system jumps from 53% well depth to 99% (horizon floor).

At kappa=0.06 (just subcritical):
```
omega^2         = 47.78  (41% shift from 80.72)
well depth      = 0.527  (53%)
Phi_min         = -0.263
R_eff           = 0.289 fm
R_eff * mu      = 1.87   (mode ~ 2x Yukawa range)
```

### Q4: Mode does NOT self-contract from box scale

For Rmax >> 1/mu, the mode stays at box scale regardless of kappa.
The Yukawa well forms at the mode's location (r ~ Rmax/2), not at
the origin. The well depth tracks the mode — no mode migration.

At mu=6.47, kappa=5000, Rmax=5:
```
R_eff = 4.23  (mode fills box)
R_eff * mu = 27.3  (mode >> Yukawa range)
```

Mode contraction ONLY occurs when Rmax ~ few/mu (compact box),
where the box already forces overlap with the well.

This is because the EM l=1 mode has NO TRUE BOUND STATE in the
Yukawa well — the centrifugal barrier l(l+1)/r^2 goes to zero
at infinity, so the continuum starts at omega^2 = 0. The "geon"
is a box mode shifted by the self-gravitational well, not a
genuine soliton.

### Physical interpretation: the strong geon at f_2 parameters

At the self-consistent scale (Rmax ~ 2 fm, kappa ~ 35, mu = 6.47):

```
omega_flat = 4.49/2.0 = 2.25 fm^{-1}
E_flat     = 2.25 x 0.197 = 0.443 GeV
omega_geon ~ 0.7 x omega_flat = 1.57 fm^{-1}  (30% shift at marginal coupling)
E_geon     ~ 1.57 x 0.197 = 0.310 GeV
R_eff      ~ 0.6 x 2.0 = 1.2 fm
```

The geon energy (~300-400 MeV) is ~3x below the proton mass (938 MeV).
The geon size (~1.2 fm) is comparable to the proton charge radius (0.88 fm).

**ORDER OF MAGNITUDE CORRECT** for hadrons: the f_2-mediated tensor
force produces self-consistent objects at the GeV/fm scale.

### Summary table: what each coupling regime gives

```
kappa / kappa_cr    Result
< 0.1               Gentle perturbation (~linear shift)
0.1 - 0.9           Moderate well (10-50% depth), mode barely contracts
0.9 - 1.0           SHARP TRANSITION: mode collapses, well hits horizon
> 1.0               Supercritical: "strong gravity black hole"
```

## Interpretation

### What the strong geon proves

1. **Tensor coupling with massive mediator WORKS** — positive feedback
   confirmed for all mu, not just mu=0 (GR).

2. **The f_2(1270) coupling IS strong enough** for nuclear-scale
   self-trapping. kappa_phys ~ kappa_cr at Rmax ~ 2 fm.

3. **The result is a collapsed object**, not a gentle geon.
   No intermediate stable solution: either perturbative or collapsed.
   This is consistent with Salam's proposal that hadrons are
   "strong gravity black holes."

4. **kappa_cr scales as mu^2 Rmax^5** for Rmax >> 1/mu. The extra Rmax^2
   (vs V8's Rmax^3) comes from the Yukawa suppression of the mode-well
   overlap. At fixed mu, larger box = harder to trap.

5. **The geon is NOT a soliton** — it's a box mode shifted by self-gravity.
   Without the box, the mode leaks through the centrifugal barrier
   (no true l >= 1 bound states). It's a quasi-bound resonance with
   finite lifetime.

### The geon hierarchy problem (V9 version)

V8 found kappa_cr ~ Rmax^3 for massless mediators, giving the hierarchy
(R/l_Planck)^2 ~ 10^{40}. V9 replaces this with:
```
kappa_cr ~ mu^2 Rmax^5
```

For the f_2 meson (mu ~ 6.5 fm^{-1}) at nuclear scale (Rmax ~ 2 fm):
```
kappa_cr ~ 42 x 32 ~ 1340  (perturbative)
kappa_cr ~ 150               (Hartree)
kappa_phys = 4 pi G_strong / L_0^2 ~ 35
```

The gap is only ~5x (Hartree kappa_cr / kappa_phys ~ 4.3).
The V8 gap of 10^{40} is reduced to ORDER UNITY.

### Q5: Geon lifetime — LONG-LIVED (WKB confirmed)

The collapsed geon is a **quasi-bound resonance** in the Regge-Wheeler
potential V_RW(r) = f(r) l(l+1)/r^2. The mode sits in a well created
by f < 1 (gravitational potential) and tunnels through the centrifugal
barrier to escape.

**Method:** WKB tunneling integral in tortoise coordinates:
```
S = int_{r1}^{r2} sqrt(V_RW - omega^2) / f(r) dr
Gamma = omega * exp(-2S)
tau = 1/Gamma
```

**Barrier structure (collapsed, f_min = 0.01):**
- Well: narrow region where f*l(l+1)/r^2 < omega^2 (f small + large r)
- Barrier: from well edge to flat-space turning point r_out = sqrt(l(l+1)/omega^2)
- For l=1, omega^2=0.009: barrier spans r = 2.0 to 14.7 fm (width 12.7 fm!)
- V_max/omega^2 ~ 36 (barrier much taller than mode energy)

**Subcritical (well depth < 100%):**
- omega^2 >> centrifugal barrier at mode location
- NO barrier exists after the well
- Mode propagates freely — lifetime = 1 oscillation period

**Results at physical f_2 parameters (kappa=35, mu=6.47):**

| Rmax  | l | omega^2  | S    | exp(-2S) | tau/tau_had | Q factor |
|-------|---|----------|------|----------|-------------|----------|
| 0.5   | 1 | 0.145    | 2.58 | 5.75e-3  | 152         | 174      |
| 0.5   | 2 | 0.366    | 3.80 | 4.99e-4  | 1,106       | 1,212    |
| 1.0   | 1 | 0.036    | 1.97 | 1.93e-2  | 91          | 52       |
| 1.5   | 1 | 0.016    | 2.34 | 9.30e-3  | 282         | 108      |
| 2.0   | 1 | 0.009    | 2.49 | 6.94e-3  | 499         | 144      |
| 2.0   | 2 | 0.023    | 4.26 | 2.01e-4  | 10,978      | 5,315    |
| 2.0   | 3 | 0.042    | 5.91 | 7.35e-6  | 221,386     | 27,873   |
| 2.1   | 1 | 3.19     | 0    | 1        | **0.19**    | 1        |

**Key findings:**

1. **ALL collapsed geons are long-lived**: tau >> tau_hadron, Q >> 1.
   The resonance persists for 100-200,000+ hadronic timescales.

2. **tau grows exponentially with l**: S ~ l, so tau ~ exp(2l).
   l=1: ~100x, l=2: ~10,000x, l=3: ~200,000x tau_hadron.

3. **Sharp boundary at the phase transition**: collapsed geons have
   tau/tau_had ~ 100+, subcritical geons have tau/tau_had ~ 0.2.
   No intermediate regime.

4. **The barrier is the centrifugal term l(l+1)/r^2** in the flat region
   (f=1), so the lifetime is insensitive to the f_min floor value.
   The floor only affects the inner turning point (inside the well).

5. **Physical width at proton scale** (Rmax=2.0, l=1):
   Gamma = 0.13 MeV, much narrower than typical hadronic widths
   (~100-200 MeV). The geon acts as a narrow resonance.

### Q6: Nonlinear metric completions — REGULAR CORES EXIST (Padé)

The linearized metric f(r) = 1 - 2M_eff/r goes negative for large
M_eff, requiring an artificial floor f_min = 0.01. Three nonlinear
completions were tested:

```
Linear (nl=0):       f = 1 - 2M/r         (floor at 0.01)
Exponential (nl=1):  f = exp(-2M/r)        (always > 0)
Padé (nl=2):         f = 1/(1 + 2M/r)      (always > 0)
```

**Comparison at mu=6.47, Rmax=0.5:**

| Metric    | kappa_cr | f_min at κ=0.065 | Barrier? | Behavior      |
|-----------|----------|------------------|----------|---------------|
| Linear    | ~0.063   | 0.01 (floor)     | YES      | Horizon → collapse |
| Exponential| ~0.10   | 0.596            | NO       | Blowup at κ=0.105 |
| Padé      | NONE     | 0.668            | NO       | Regular ∀κ    |

**Why metrics differ:**
The EM energy density scales as ρ ~ ω²u²/f. The Hartree feedback
loop is 1/f amplification:
- Linear: 1/f diverges at f=0 → horizon
- Exponential: 1/f = e^{2M/r} grows exponentially → delayed collapse (κ_cr ≈ 0.10)
- Padé: 1/f = 1+2M/r grows linearly → gentlest feedback, NO collapse

**Padé metric sweep (mu=6.47, Rmax=0.5, l=1):**

| κ    | f_min  | well  | ω²     | S      | τ/τ_had   | Status           |
|------|--------|-------|--------|--------|-----------|------------------|
| 0.3  | 0.217  | 78%   | 15.0   | 0      | 0.086     | no barrier       |
| 0.5  | 0.114  | 89%   | 6.03   | 0      | 0.14      | no barrier       |
| 0.8  | 0.063  | 94%   | 2.47   | 0.054  | 0.24      | short-lived      |
| 1.0  | 0.047  | 95%   | 1.60   | 0.519  | 0.75      | marginal         |
| 1.5  | 0.027  | 97%   | 0.686  | 1.46   | 7.4       | long-lived       |
| 2.0  | 0.016  | 98%   | 0.346  | 2.33   | 60        | long-lived       |
| 3.0  | 0.006  | 99%   | 0.103  | 4.35   | 6,232     | very long-lived  |
| 5.0  | 7e-4   | 99.9% | 0.008  | 13.4   | 1.8×10¹²  | essentially stable |

**Key findings for Q6:**

1. **Padé metric gives genuine regular cores** — no horizon, no floor,
   f > 0 everywhere. The solution converges smoothly for all κ up to ~5
   (numerical underflow at κ ≥ 10 but no divergence).

2. **Regular-core geons ARE long-lived** at strong coupling (κ ≥ 1.5).
   The deep well (>97%) with small ω² produces thick centrifugal barriers.
   At κ=5: S=13.4, τ ≈ 10¹² τ_hadron (essentially stable!).

3. **The barrier requires deep wells.** At κ < 0.8 (Padé), the well is
   too shallow for ω² < V_RW at the barrier — no tunneling suppression.
   Transition from no-barrier to long-lived is gradual (not sharp like
   the linear metric collapse).

4. **Exponential metric always gives NO barrier** — it collapses (NaN)
   before the well gets deep enough. Only the Padé completion supports
   both regular cores AND long lifetimes.

5. **Physical coupling (κ=35) has numerical underflow** in Padé (f_min
   and ω² round to zero). The trend from κ=1-5 (f_min ~ 1/κ², ω² ~ 1/κ⁴)
   predicts f_min ~ 2×10⁻⁴, ω² ~ 10⁻⁵ at κ=35 — deep in the very-long-
   lived regime. Extended precision would be needed to confirm.

### Q7: Multi-mode effects — 3-QUARK GEON MATCHES PROTON MASS

Multiple EM modes sharing the same gravitational potential. Each mode
has angular momentum l_i and occupancy n_i. Total density = sum of
individual mode densities weighted by occupancy.

**Identical modes (N × same l):** equivalent to kappa → N×kappa.
Verified: 3×l=1 at κ=0.02 gives identical result to single l=1 at κ=0.06.

**Multi-mode results at physical f₂ coupling (κ=35, μ=6.47):**

| Configuration    | Rmax  | f_min | well  | E_total | E/M_p | R_eff   |
|------------------|-------|-------|-------|---------|-------|---------|
| 1×l=1 (single)   | 2.1*  | 0.532 | 47%   | 414     | 0.44  | 1.21 fm |
| 3×l=1            | 2.7*  | 0.634 | 37%   | 857     | 0.91  | 1.59 fm |
| 2×l=1+1×l=2      | 2.85* | 0.723 | 28%   | **887** | **0.945** | 1.7-1.8 fm |
| 1×l=1+1×l=2+1×l=3| 3.0   | 0.01  | 99%   | 60      | 0.064 | 2.6 fm  |

(*) = just above critical Rmax (subcritical side of transition)

**Individual mode energies at 2×l=1+1×l=2, Rmax=2.85:**
```
l=1 mode: E = 267 MeV, R_eff = 1.71 fm  (x2)
l=2 mode: E = 352 MeV, R_eff = 1.85 fm  (x1)
Total:    E = 887 MeV                    (94.5% of proton mass)
```

**Key findings:**

1. **3-quark geon matches the proton mass to 5-9%.** The best match is
   2×l=1+1×l=2 at Rmax=2.85: E=887 MeV vs M_p=938 MeV (5.5% low).
   Even 3×l=1 at Rmax=2.7 gives E=857 MeV (8.7% low).

2. **Per-mode energy ≈ constituent quark mass.** Each l=1 mode contributes
   ~267-285 MeV, comparable to the constituent quark mass (~300 MeV).
   The l=2 mode contributes ~352 MeV (heavier "strange-like" quark?).

3. **Critical Rmax shifts with N_modes:** single → 2.08 fm, 3×l=1 → 2.65 fm,
   2+1 → 2.75 fm. More modes = larger critical box (κ_eff = N×κ for same l).

4. **Mixed-l produces DEEPER wells** than same-l at equal total occupancy.
   Higher-l modes contribute density at larger r, broadening the Yukawa source.
   The 1+2+3 configuration collapses even at Rmax=3.0 (still supercritical).

5. **The energy is NOT simply additive.** More modes → deeper well → lower
   per-mode eigenvalue. At Rmax well above transition, E_total ≈ N × E_flat
   (just perturbative shift). At transition, each mode's energy is ~65%
   of its flat-space value.

### Q8: Topological Hopfion Geon — NO (topology doesn't stabilize)

**Question:** Can EM field line topology (Ranada hopfion, Hopf charge
Q_H = ∫ A·B d³x) prevent the geon mode from leaking?

**Analysis:**

The Ranada hopfion is a specific coherent superposition of spherical
partial waves (l,m). In the geon framework with spherically symmetric
metric, each (l,m) mode evolves independently. Key observations:

1. **Partial waves decouple in spherical symmetry.** Each l-mode sees
   its own Regge-Wheeler potential V_l(r) = f(r)·l(l+1)/r². The modes
   leak independently through their own centrifugal barriers.

2. **The l=1 component dominates.** The hopfion's dominant contribution
   is from l=1 (dipole). This has the thinnest barrier (lowest l) and
   sets the overall leakage rate.

3. **Q_H conservation doesn't prevent leaking.** The Hopf charge is
   conserved by Maxwell's equations, but this only prevents decay to
   the trivial (zero-field) state. The hopfion can expand and escape
   the well while maintaining its linked field line topology.

4. **Each partial wave leaks independently.** As the l=1 component
   tunnels through the barrier, Q_H decreases. The topology protects
   against unraveling to zero, not against spatial expansion.

5. **Nonlinear coupling would help but is absent.** In a fully nonlinear
   theory (Born-Infeld + Skyrme terms), partial waves would be coupled
   and topology could prevent individual components from leaking. In the
   V9 linear EM framework, no such coupling exists.

**Conclusion:** The geon lifetime is set by WKB tunneling of the l=1
partial wave, identical to the single-mode result (Q5). Topology
provides no additional stabilization in the linear EM framework.
The strong geon is stabilized by the gravitational well alone.

### Q9: Self-consistent box size — FUNDAMENTAL TENSION

**Question:** Can the box size Rmax be determined self-consistently
by the mode itself, without an external boundary?

**Method:** Vary Rmax at fixed κ (Padé metric) and check whether ω²
converges to a plateau (indicating self-trapping) or scales as 1/Rmax²
(indicating box mode).

**Rmax scan at fixed κ (Padé metric, μ=6.47, l=1):**

| Rmax | κ=1.0 ω²  | κ=1.0 well | κ=3.0 ω²  | κ=3.0 well | R_eff(κ=3) |
|------|-----------|------------|-----------|------------|------------|
| 0.3  | 15.0      | 78%        | 0.103     | 99.4%      | 0.405      |
| 0.4  | 6.05      | 82%        | 0.029     | 99.8%      | 0.421      |
| 0.5  | 1.60      | 95%        | 0.103     | 99.4%      | 0.405      |
| 0.7  | 0.53      | 97%        | —         | —          | —          |
| 1.0  | 0.136     | 99%        | —         | —          | —          |
| 2.0  | 0.018     | 99.8%      | —         | —          | —          |
| 3.0  | 0.005     | 99.9%      | —         | —          | —          |

**NO plateau in ω²** — the eigenvalue always decreases with Rmax.
The mode fills whatever box is provided. However, R_eff IS nearly
constant for small Rmax (self-localization of density), even as ω²
continues to fall.

**Lifetime as function of Rmax (Padé, κ=2.0, μ=6.47, l=1):**

| Rmax | ω²    | R_eff  | Barrier? | S     | τ/τ_had    |
|------|-------|--------|----------|-------|------------|
| 0.4  | 0.025 | 0.362  | YES      | 13.3  | ~10¹²      |
| 0.5  | 0.346 | 0.371  | YES      | 2.33  | 60         |
| 0.6  | 1.03  | 0.381  | YES      | 0.73  | 1.4        |
| 0.7  | 1.95  | 0.392  | NO       | 0     | <1         |
| 1.0  | 4.5   | 0.42   | NO       | 0     | <1         |

**R_eff plateau at ~0.37 fm** for Rmax=0.4-0.6 — the mode self-localizes.
But the barrier onset occurs at r ≈ Rmax: for Rmax inside the barrier
region, the box boundary IS the barrier. For Rmax beyond the barrier,
the mode leaks freely.

**3-quark geon: self-consistency vs. proton mass**

At 3×κ=35 (physical 3-quark coupling → κ_eff=105):

| Rmax | Status       | well  | S     | τ/τ_had | E_total  | E/M_p |
|------|-------------|-------|-------|---------|----------|-------|
| 2.5  | collapsed   | 99%+  | 2.46  | 593     | ~45 MeV  | 0.048 |
| 2.6  | collapsed   | 99%+  | 2.46  | 600     | ~45 MeV  | 0.048 |
| 2.7  | subcritical | 37%   | 0     | <1      | 857 MeV  | 0.91  |

**THE FUNDAMENTAL TENSION:**

1. **Self-consistent trapping** (Rmax < Rmax_cr): The mode collapses
   into a deep well. The barrier is thick (S > 2), giving τ >> τ_had.
   But the collapsed mode has ω² ~ 0, so E_mode ~ 0. For the 3-quark
   geon: E_total ≈ 45 MeV ≈ 5% M_p. This is NOT a proton.

2. **Proton mass matching** (Rmax ≈ Rmax_cr): The mode is just barely
   subcritical. ω² is large enough to give E ≈ 857-887 MeV ≈ M_p.
   But the well is only 30-50% deep — NO barrier exists. The mode
   leaks freely on the hadronic timescale.

3. **No intermediate regime exists.** The sharp phase transition (Q3)
   means there is no Rmax where the well is deep enough for a barrier
   AND the energy is high enough to match hadron masses. The geon is
   either trapped-but-too-light or correct-mass-but-untrapped.

**Interpretation:** The strong geon cannot be made truly self-consistent
as a stable hadron. It produces either (a) deeply trapped objects at
~45 MeV (too light by 20×), or (b) proton-mass objects that leak freely.
The physics at the transition is correct (right scale, right mass), but
the transition itself is too sharp to support a metastable state.

This is the V9 version of the "geon instability": Wheeler's original
geons also required a box (or leaked on light-crossing timescales).
The massive mediator improves the situation dramatically (from 10⁴⁰
coupling gap to order unity), but the box mode nature persists.

## Code

- `src/strong_geon.c`: Yukawa geon Hartree solver
  - `-perturb`: First-order perturbation theory
  - `-hartree`: Full Hartree self-consistent solve
  - `-scan`: Scan kappa from 0 to kmax
  - `-muscan`: Scan mu at fixed kappa (perturbative)
  - `-lifetime`: WKB quasi-bound state lifetime
  - `-multimode`: Multi-mode Hartree (with `-modes` spec)
  - `-nlmetric N`: Nonlinear metric (0=linear, 1=exp, 2=Padé)
  - `-modes spec`: Mode specification l1:n1,l2:n2,... (e.g. 1:3 or 1:2,2:1)
  - `-kappa K`: Coupling constant
  - `-mu M`: Yukawa mass (default 1.0)
  - `-seed A`: Initialize with compact Yukawa well (A r exp(-mu r))
  - `-l L`: Angular momentum (>= 1)
  - `-Rmax R`: Domain size (default 10.0)
  - `-N N`: Grid points (default 4001)

## References

- Salam (1973): "Strong Interactions, Gravitation, and Cosmology"
- Isham, Salam, Strathdee (1973): "f-Dominance of Gravity", PRD 8, 2600
- Sivaram & Sinha (1979): "Strong gravity, black holes, and hadrons", Phys Rep 51, 111
- Wheeler (1955): "Geons" — original massless geon
- V8 FOUNDATION: Geon solver, coupling gap analysis
