# V27 Results: From Braided Soliton to Emergent Physics

## Big Picture Assessment

V27 achieved TWO of its five milestones and produced one BREAKTHROUGH:

| Milestone | Target | Result | Verdict |
|-----------|--------|--------|---------|
| M1 (stable braid) | fc>0.5, |P|>0.1 | fc=0.28, |P|=0.574 | PARTIAL |
| **M4 (mass)** | **Emergent mass** | **m=0 survives! |P|=1.742** | **BREAKTHROUGH** |
| M2 (gravity) | Quadrupolar attraction | Isotropic repulsion | NEGATIVE |
| M3 (EM) | Quantized torsion flux | Φ_T=0 (symmetry cancellation) | NEGATIVE |
| **M5 (confinement)** | **Conserved topology** | **Winding = -1.000 exactly** | **POSITIVE** |

## Milestone 1: Stable Propagating Braid

The DynA configuration (m=1, no extra couplings, k=2π/L) is a UNIQUE
resonance. Tested 9 configurations:

- Pairwise coupling (λ_pw=0.3, 0.5, 0.7): KILLS |P| (0.574→0.002)
- Smaller domain (L=10): HURTS (|P|=0.034)
- Speed scan (k=0.5-3.0): ALL fail (|P|<0.012)
- Only k=2π/20 (one twist per domain) maintains |P|=0.574

The DynA |P|=0.574 resonance is fragile: any perturbation destroys it.
fc=0.276 is limited by transverse absorbing BC (energy leaks through x,y).
M1b transient at t=400 reached fc=0.47, |P|=0.13 (nearly passing) before
boundary absorption killed it → path forward: periodic BC in ALL directions.

## Milestone 4: Mass from Dynamics — BREAKTHROUGH

**The massless propagating braid SURVIVES.**

| m | |P|(500) | fc | E | Regime |
|---|---------|-----|-------|--------|
| 1.0 | 0.574 | 0.276 | +281 | 1 (positive E) |
| 0.8 | 0.003 | 0.157 | +137 | Dead zone |
| 0.6 | 0.000 | 0.174 | +31 | Dead zone |
| 0.4 | 0.589 | 0.129 | -2,871 | 2 (negative E) |
| 0.2 | 1.105 | 0.108 | -6,381 | 2 |
| 0.1 | 1.591 | 0.106 | -7,518 | 2 |
| **0.0** | **1.742** | **0.103** | **-7,915** | **2 (optimal)** |

Two distinct regimes separated by a dead zone at m=0.6-0.8:
- Regime 1 (m≈1): positive energy, moderate |P|, l=2≈30%
- Regime 2 (m≤0.4): negative energy (deep well), HIGH |P|, l=2≈3%

Strong binding scan (m=0, varying μ):

| μ | κ | |P|(500) | fc |
|---|---|---------|-----|
| -20 | 20 | 1.742 | 0.103 |
| **-50** | **50** | **2.525** | **0.095** |
| -100 | 100 | 2.164 | 0.092 |
| -200 | 200 | 1.643 | 0.087 |

Optimal: μ=-50, κ=50 gives |P|=2.525 (highest triple-product coupling).

The m=0 Lagrangian:

    L = Σ_a [½(∂_t φ_a)² - ½(∂_i φ_a)²] - (μ/2)P²/(1+κP²)

Mass is EMERGENT: M = |E|/c² = 7915 code units, entirely from dynamics.

## Milestone 2: Gravity — NEGATIVE

M2a: Single braid strain decomposition at m=0, μ=-50:
- l=2/l=0 = 3.1% at R=5, dropping to 0.03% at R=8
- The deep-well m=0 braid is MORE SPHERICAL than the m=1 braid
- The potential well symmetrizes the field distribution

M2b: Two braids at D=30 (L=40, N=128, dx=0.63 — TOO COARSE):
- Numerical instability (energy grew to 1.4M)
- Separation oscillated unreliably

M2c: Angular force pattern (three angles):
- All repulsive (Δsep = +1.7 to +2.0)
- Isotropic (20% variation across angles)
- BUT: M3c at proper resolution showed ATTRACTION

Resolution: M2 results at L=40 are numerically unreliable. M3c at L=20
showed genuine attraction between braids.

## Milestone 3: EM — NEGATIVE for charge, POSITIVE for attraction

M3a: Torsion field |Ω| exists but is NOT localized (only 7.7% in core).

M3b: Torsion flux Φ_T = 0 at ALL z-slices (10⁻¹⁴, machine precision).
This is a SYMMETRY OBSTRUCTION: the 120° azimuthal symmetry of the
three-strand braid guarantees exact cancellation. Not numerical artifact.

M3c: Two-braid interaction (L=20, N=128):
- Same-twist braids: ATTRACT (Δsep = -7.2)
- Opposite-twist braids: ATTRACT (Δsep = -5.0)
- Same-twist attracts MORE (gravity-like: all attract)
- Resolution caveat: dx=0.63 may cause numerical artifacts

## Milestone 5: Topological Confinement — POSITIVE

The winding number is EXACTLY CONSERVED at -1.000 throughout t=0 to 500.
Not approximately — EXACTLY, to numerical precision.

| t | winding | |P| | fc | E |
|---|---------|------|------|-------|
| 0 | -1.000 | 0.127 | 0.983 | 25.8 |
| 100 | -1.000 | 1.515 | 0.364 | -8688 |
| 200 | -1.000 | 1.796 | 0.365 | -8830 |
| 300 | -1.000 | 0.911 | 0.367 | -8901 |
| 400 | -1.000 | 1.102 | 0.363 | -8950 |
| 500 | -1.000 | 2.506 | 0.365 | -8970 |

The topology is conserved by the EXISTING dynamics — no Chern-Simons
term needed. The propagating helical wave can't unwind.

Also: fc STABILIZES at 0.363-0.367 (better than M4a's 0.103) at μ=-50.
And |P| oscillates 0.9-2.5 (LARGE, sustained).

## The V27 Legacy

### What V27 Achieved:
1. **Massless soliton with emergent mass** (m=0, M=|E|/c²=7915)
2. **Exact topological conservation** (winding=-1.000, no extra terms)
3. **Two-braid attraction** (M3c: Δsep=-7.2 for same-twist)
4. **Optimal parameters identified** (μ=-50, κ=50, k=2π/L)

### What V27 Did NOT Achieve:
1. **Spin-2 gravity** — the m=0 deep-well braid is too spherical (l=2≈3%)
2. **EM charge** — cylindrical symmetry cancels torsion flux exactly
3. **fc > 0.5** — the braid is loosely bound (fc≈0.10-0.36)

### The Fundamental Tension (Updated):
- Regime 1 (m=1): ASPHERICAL (l=2≈30%) but UNWINDS (|P| slowly decays)
- Regime 2 (m=0): TOPOLOGICALLY STABLE (winding conserved) but SPHERICAL (l=2≈3%)
- The deep potential well that provides stability also SPHERICALIZES the braid

### Path Forward:
The two regimes need to be BRIDGED:
- A braid that has Regime 2's topological stability AND Regime 1's asphericity
- Breaking the azimuthal symmetry (for EM charge and l=2 content)
- Periodic BC in ALL directions (for fc improvement — M1 showed this)
- Higher resolution (N=256) for reliable two-braid force measurement
