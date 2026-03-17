# Combo 1+2+5+7: EOS + Inertia + Condensate + Self-Ref = Full GR Analog

## Thesis

Combine ALL the working pieces into one simulation:
- Test A (EOS): P/ρ≈0 in core (pressureless dust, virial balance)
- Test B (Inertia): a=F/M, quadratic deformation, F_crit≈2e-5
- Test E (Lattice): λ=0.5 stabilizes chain, c_s≈0.9, phonons real
- Test F (Self-ref): α=0.0001 gives Φ=-0.006, self-consistent

Build the stable pairwise-coupled lattice AND evolve it with the
self-consistent gravitational potential Φ from Poisson. Perturb and
measure whether the combined system produces:
1. A modified phonon spectrum (from the Φ correction)
2. A self-consistent lattice metric
3. Any coupling between the acoustic branch and the metric

## Setup

The lattice (from Test E) with pairwise coupling λ=0.5, PLUS the
self-consistent metric (from Test F) with coupling α.

At each timestep:
1. Compute total energy density ρ(x) = Σ grid contributions
2. Solve Poisson: ∂²Φ/∂x² = α·ρ(x) (tridiagonal Thomas algorithm)
3. Modify dynamics: m²_eff = m²(1+2Φ), c²_eff = 1+4Φ
4. Evolve fields with modified dynamics

The key question: does the Φ field COUPLE to the phonon modes?

In a uniform medium: Φ = const (no gradient) → no force on phonons.
In a NON-uniform medium (lattice): Φ has structure (periodic from the
oscillon array) → phonons feel the Φ gradient → gravitational interaction.

## Method

### Phase 1: Self-consistent lattice

1. Build 8-oscillon chain at λ=0.5, d=16 (from Test E)
2. Turn on self-consistent Φ with α=0.0001 (from Test F, safe value)
3. Evolve t=10000
4. Does the lattice survive? Does Φ converge to a periodic profile?
5. Measure Φ(x) — does it have a periodic structure matching the lattice?

### Phase 2: Phonon spectrum with gravity

6. At equilibrium: add small random displacements to oscillon positions
7. Compute phonon spectrum (normal mode decomposition + DFT)
8. Compare with Test E spectrum (no gravity):
   - Is c_s different? (gravity modifies the effective EOS → different sound speed)
   - Are new modes present? (gravitational waves?)
   - Does the antisymmetric branch change?

### Phase 3: Gravitational wave test

9. Apply a localized perturbation to one oscillon
10. Track: does Φ(x,t) show a PROPAGATING wave?
11. If Φ propagates: what speed? (Should be related to c²_eff)
12. Is the Φ wave correlated with the density (phonon) wave?

### Phase 4: Equivalence principle check

13. Perturb oscillons #3 and #6 with DIFFERENT amplitudes
14. Do they respond to the gravitational field with the SAME acceleration?
15. If a₃ = a₆ (same Φ gradient): equivalence principle holds
16. If a₃ ≠ a₆: mass-dependent coupling → no equivalence

## Key Predictions

If the system is truly GR-like:
- Phonon speed should depend on Φ: c_s(Φ) = c_s(0)·√(1+4Φ)
- The Φ wave should propagate at c_eff (or faster)
- The equivalence principle should hold (same a for different M)

If NOT GR-like:
- c_s unchanged by Φ (gravity doesn't affect sound)
- No propagating Φ wave (gravity is instantaneous/Poisson)
- Equivalence violated (coupling depends on oscillon properties)

## Reference Code

- v24/fundamental/testE_lattice/src/lattice.c (lattice chain)
- v24/fundamental/testF_selfref/src/selfref.c (Poisson solver + metric)
- v24/fundamental/combo_25/src/combo25.c (perturbation + phonon analysis)

## Output

- `src/combo1257.c`, `data/`, `RESULTS.md`

## Parameters

λ=0.5, α=0.0001, μ=-20, κ=20, m=1.0
N_osc=8, d=16, periodic BC
Nx = 8*320 = 2560
t_equil=5000, t_gravity=10000

Compile: `gcc -O3 -Wall -o combo1257 src/combo1257.c -lm`
