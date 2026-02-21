# V7: Disformal Self-Trapping and Causal Inertia

## Core Idea

Particles are self-trapped massive scalar waves in a Bekenstein disformal metric.
The wave field ψ sources a massless scalar σ, which modifies the effective metric
in which ψ propagates. The self-consistent solution confines ψ via a potential
well V(r) = m²A(r) < m², where A = e^{2ασ} < 1 (conformal factor depressed
by negative σ). The σ field decays as D/r at infinity, mediating a long-range
1/r scalar force between solitons.

## Model

**Fields**: Massive scalar ψ (wave), massless scalar σ (inertia/mediator)

**Disformal metric** (Bekenstein 1993):
```
g̃_μν = e^{2ασ} η_μν + β ∂_μσ ∂_νσ
```
where α = conformal coupling, β = disformal coupling.

**Wave equation** (ψ propagates in g̃):
```
□̃ψ + m²ψ = 0   →   ü = v²u'' + (C_r'/C_t)(u'-u/r) - l(l+1)u/r² - m²Au
```
with u = rψ, v² = A/(A+βσ'²), A = e^{2ασ}.

**Sigma equation** (σ sourced by ψ "Lagrangian density"):
```
□σ = α_s · (∂ψ)²   ≈   -α_s · m² · ψ²   (for bound states)
```
The source ⟨(∂ψ)²⟩ ≈ -m²ψ² < 0 for bound states (kinetic ≈ mass energy),
giving σ < 0 → A < 1 → potential well.

**Key sign**: Source pushes σ NEGATIVE, creating A < 1 (well), NOT positive
(barrier). This works for the massive scalar proxy where ⟨(∂ψ)²⟩ ≈ -m²ψ² < 0.
However, for massless EM: ⟨B²-E²⟩ > 0 in slow-light media (see below).

## Self-Consistent (Hartree) Solution

Iterative scheme:
1. Trial ψ → compute source S(r) = -α_s·m²·ψ²
2. Solve σ BVP: σ'' + 2σ'/r = S(r) with Neumann BC (σ → D/r)
3. Compute V_eff(r) = m²A(σ) + l(l+1)/r²
4. Find eigenstate ψ in V_eff (implicit imaginary time)
5. Under-relax σ, repeat until convergence

**Result**: CONVERGES for ALL tested coupling strengths α_s ∈ [0.5, 20].

## Key Numerical Results

### Coupling scan (α_c=0.3, m_ψ=1, μ_σ=0, l=1):

| α_s | ω² | E_bind | σ_min | A_min | D |
|------|-------|--------|-------|-------|-----|
| 0.5 | 0.998 | 0.002 | -0.04 | 0.98 | -0.5 |
| 1.0 | 0.989 | 0.011 | -0.09 | 0.95 | -1.0 |
| 2.0 | 0.963 | 0.037 | -0.18 | 0.90 | -2.0 |
| 5.0 | 0.826 | 0.174 | -0.54 | 0.72 | -5.0 |
| 10 | 0.598 | 0.402 | -1.10 | 0.52 | -10 |
| 20 | 0.348 | 0.652 | -2.10 | 0.28 | -20 |

### Exact result: D = -α_s

The dilaton charge D = -α_s·m²·∫ψ²d³x/(4π) = -α_s (for normalized ψ
with ∫u²dr = 1, giving ∫ψ²4πr²dr = 4π). This is EXACT by Gauss's law.

### Strong coupling (α_c=0.5, α_s=10, l=1):
- ω² = 0.409, E_bind = 0.591 (59% of rest mass)
- σ_min = -1.50, A_min = 0.223 (wave speed reduced to 0.47c at center)
- V_min = 0.319 (deep potential well)

### l-dependence (α_c=0.3, α_s=5):
- l=0: ω² = 0.653, E_bind = 0.347 (no centrifugal → stronger binding)
- l=1: ω² = 0.826, E_bind = 0.174

## Force Between Solitons

Two solitons at separation R interact via their σ fields:
```
V(R) = D₁·D₂ / (4πR) = α_s² / (4πR)   [attractive for same-sign D]
```
This gives effective gravitational coupling:
```
G_eff = D²/(4πM²) = α_s²/(4πM²)
```

**Problem**: G_eff is a free parameter (α_s). For gravity: need α_s ~ 10⁻²⁰.
The mechanism works but doesn't PREDICT the coupling strength.

## What Works

1. **Self-trapping**: ψ creates σ well, σ confines ψ — fully self-consistent ✓
2. **Long-range force**: σ ~ D/r at infinity (massless mediator) ✓
3. **Universal attraction**: all solitons have D < 0, same-sign scalar = attractive ✓
4. **Disformal refraction**: additional confinement from v_r < c in σ' shell ✓
5. **Bound states for all couplings**: no minimum threshold ✓

## What Fails / Limitations

1. **Scalar, not tensor**: σ mediates spin-0 force (no gravitational waves) ✗
2. **Coupling is free**: G_eff = α_s²/(4πM²) not determined by theory ✗
3. **Same conclusion as V6**: 1/r scalar force via massless mediator
4. **Time evolution unstable**: parametric resonance for μ_σ=0 (oscillating
   source at 2ω drives σ). Static Hartree works but full dynamics crashes.
5. **No tensor gravity**: would need spin-2 massless field (not in this model)

## Comparison with Previous Investigations

| Model | Source | Feedback | 1/r? | Trap? | Lesson |
|-------|--------|----------|------|-------|--------|
| V2 BLV | Skyrme L₂+L₄ | N/A (topological) | No (Yukawa) | Topological | Core-scale only |
| V2 L₆ | Baryon density | N/A (topological) | No (Yukawa) | Topological | Core-scale only |
| V6 | |ω|² → □δρ | Depletion (neg.) | Yes | No | Scalar 1/r, 10⁴⁰× too strong |
| V7 massive ψ | -m²ψ² | Bypassed by m² | Yes | Yes (m² input) | Mass must be input |
| **V7 massless F** | **B²-E²** | **Negative (fatal)** | **No** | **No** | **No-go theorem** |
| V4 EBId | T_μν (GR) | Positive (Shapiro) | Yes | Yes (GR) | Tensor required |

The progression V6 → V7 scalar → V7 vector establishes the **scalar no-go
theorem**: any scalar mediator coupled to the EM Lagrangian F² = B²-E²
suffers negative feedback (depletion destroys the well). Only tensor
coupling to T_μν ∝ B²+E² provides the positive feedback needed for
self-trapping of massless waves.

## Massless Vector Attempt (breather.c) — NEGATIVE RESULT

### Motivation

The scalar proxy uses massive ψ with source □σ = -α_s·m²·ψ². But "matter is light"
requires massless F_μν with source □σ = α_s·F² = α_s·(B²-E²). Mass should emerge
from self-trapping energy, not as a Lagrangian parameter.

### Implementation (breather.c)

Three-phase Hartree solver:
1. **Bootstrap**: massive scalar source (reproduces selftrap.c, ω²=0.826 at α_s=5)
2. **Switch to F² source**: keep m² in eigenvalue equation, change source to BI-saturated F²
3. **Remove mass**: set m²=0, pure massless vector

F² for TM_l standing wave: ⟨F²⟩ = [(u'-u/r)² - ω²u²]/(2r²)

BI saturation: source = f_sign·α_s·F²/√(1+|F²|/(2b²))

### Results: DEFINITIVE NEGATIVE

Phase 3 (m²=0) always converges to **delocalized box mode**:

| α_s | Rmax | ω² | D_charge | ω²·Rmax² |
|-----|------|----|----------|----------|
| 20  | 20   | 0.0505 | -0.195 | 20.2 |
| 50  | 20   | 0.0505 | -0.488 | 20.2 |
| 100 | 20   | 0.0505 | -0.977 | 20.2 |
| 20  | 5    | 0.709  | -2.015 | 17.7 |
| 20  | 8    | 0.315  | -1.211 | 20.2 |

ω²·Rmax² ≈ 20.2 = (4.49)² = first zero of spherical Bessel j₁.
This is the box ground state, NOT a localized bound state.

Mass quench (F² source, α_s=5, l=1, Rmax=40):

| m² | ω² | D_charge |
|----|-----|----------|
| 1.0 | 0.762 | -3.85 |
| 0.1 | 0.096 | -1.23 |
| 0.01 | 0.015 | -0.30 |
| 0.0 | 0.013 | -0.04 |

Smooth transition from bound (m²=1) to delocalized (m²=0). No discontinuity.

Fixed-well test (selftrap.c, m_psi=0): wave leaks with half-life ≈ 9 code units
(quasi-bound resonance, not true bound state).

### Root Cause: F² Sign in Slow-Light Media

**In a disformal well with v < c**:
- Wave speed v = c·√A < c → spatial wavenumber k = ω/v > ω/c
- Therefore k² > ω²/c² always → B² > E² always
- F² = B²-E² > 0 everywhere inside the well

**For massive scalar**: (∂ψ)² = k² - ω²_eff, where ω²_eff = ω² - m²A(r).
Since m²A < m², we get ω²_eff > ω² - m², so (∂ψ)² CAN be negative.
This breaks the B-E cancellation and allows self-trapping.

**For massless vector**: k > ω always in slow medium → F² > 0 always.
The near-cancellation B ≈ E makes the source inherently weak (~10⁻⁴).
No self-consistent localized solution exists.

### Conclusion: Scalar No-Go Theorem

Massless vector disformal self-trapping via F² coupling does NOT work.
The mechanism requires mass as an input parameter — mass does not emerge
from the geometry alone. This is a fundamental obstruction, not a numerical
artifact.

## The Optical Refraction No-Go Theorem (V6 + V7 Synthesis)

The V7 breather failure is not an isolated result — it is the numerical
manifestation of a fundamental negative-feedback loop that connects the
V6 density depletion result with the V7 wave-trapping attempt.

### The Paradox

1. **Trapping requirement**: Self-trapped waves need a converging lens at the
   core — a slow-light medium (v_core < v_background). Waves refract toward
   the slower medium.

2. **Depletion reality (V6)**: The twist/energy of the field depletes the
   local vacuum density ρ. The core is thinner than the background.

3. **Speed consequence**: Lower density → less inertia → faster wave
   propagation. Depleted core = fast-light medium.

4. **Fatal feedback**: A fast-light core is a diverging lens. The wave bends
   outward and scatters instead of self-focusing.

### The Negative Feedback Loop

```
Wave tries to slow down (trap itself)
  → B > E in slow medium (spatial bunching)
  → F² = B²-E² > 0
  → Source depletes σ (or ρ)
  → Depleted density INCREASES local wave speed
  → "Well" becomes "hill"
  → Wave expelled to box boundary
  → ω² = (4.49/Rmax)² — the box mode
```

The V7 code didn't fail — it mathematically proved this feedback loop.
The `20.2` result IS the theorem.

### Why the Massive Scalar Evades This

The massive scalar proxy (selftrap.c) works because m² provides a
negative-definite contribution to the source: ⟨(∂ψ)²⟩ ≈ -m²ψ² < 0.
This overwhelms the B²-E² cancellation by construction. But this puts
mass in by hand — the trapping relies on an input mass, not geometry.

### The No-Go Statement

**No scalar density field can self-trap a massless wave**, because the
coupling is through F² = B²-E², which has the wrong sign in slow-light
media. The act of creating a slow-light well generates F² > 0, which
depletes the well, creating a negative feedback that always wins.

This applies to:
- V6 (density conservation): depletion → 1/r⁶ locally, 1/r only via □
- V7 (disformal scalar): F² > 0 → anti-well → box mode
- Any scalar mediator coupled to the EM Lagrangian density

### The Tensor Resolution

| Coupling | Source | Sign in well | Feedback |
|----------|--------|--------------|----------|
| Scalar (F²) | B²-E² | Positive (wrong) | **Negative** — well depletes itself |
| Tensor (T_μν) | B²+E² | Always positive | **Positive** — energy deepens well |

General Relativity couples to the energy-momentum tensor T_μν ∝ B²+E²,
which is a sum of squares — always strictly positive. A concentration
of wave energy creates a Shapiro time delay (geometric slow-light well)
that DEEPENS with more energy, not depletes. This is positive feedback:
more energy → deeper well → tighter confinement → more energy density.

The no-go theorem forces the theory from scalar density mechanics to
tensor geometric curvature as the trapping mechanism.

## Code

- `src/selftrap.c`: Massive scalar proxy solver (~960 lines)
  - `-fixed`: Fixed σ well, evolve ψ only
  - `-coupled`: Coupled time evolution (unstable for μ=0)
  - `-hartree`: Hartree self-consistent solver (WORKS)
  - `-veff`: Print effective potential
  - `-standing`: Standing wave initialization
- `src/breather.c`: Massless vector breather solver (~810 lines)
  - `-hartree`: Three-phase Hartree (bootstrap → F² → massless) — NEGATIVE
  - `-quench`: Mass quench scan m²: 1→0
  - `-coupled`: Coupled BI time evolution — UNSTABLE
  - `-fsign +/-1`: F² source sign convention
- `src/Makefile`: Build with `make`

## Files
```
v7/
├── FOUNDATION.md       (this file)
├── src/
│   ├── selftrap.c      (massive scalar proxy — WORKS)
│   ├── breather.c      (massless vector — NEGATIVE)
│   └── Makefile
├── bin/
│   ├── selftrap
│   └── breather
└── results/            (output data)
```
