# Test E: Pairwise-Coupled Oscillon Lattice (Condensate + Time Crystal)

## Thesis

V23-D showed the oscillon lattice is too weakly bound (well depth 4e-4).
V24-ME showed pairwise coupling stabilizes oscillons and creates long-range
mediators. Combine: a lattice of oscillons WITH pairwise coupling.

The pairwise coupling should:
1. Deepen the inter-oscillon potential well (stronger binding)
2. Provide a long-range Proca mediator for lattice cohesion
3. Phase-lock the oscillons into a time crystal (triple product synchronizes)

## Method

1. Equilibrate a single oscillon WITH pairwise coupling λ=0.5 (moderate,
   well within stability). Save the profile.
2. Build an 8-oscillon chain (periodic BC, spacing d=16) from the
   equilibrated profile. Same as V23-D Phase 2 but with pairwise coupling.
3. Evolve for t=10000.
4. Track all 8 positions X_n(t).
5. Measure: does the chain survive (max displacement < d/2)?
6. If stable: compute the phonon spectrum via normal modes.
7. Compare with V23-D (no pairwise): is the well depth larger?
   Is the phonon speed different? Are the modes longer-lived?

## Key Prediction

At λ=0.5: m_A = 0.71, Proca range = 1.4. The inter-oscillon potential
should have contributions from BOTH the direct tail overlap (range ~1)
AND the Proca mediator (range ~1.4). The total well depth should be
deeper than V23-D's 4e-4.

At λ=0.9: m_A = 0.32, range = 3.1. Even deeper well, longer-range
cohesion. But the oscillon profile is larger (more pairwise energy),
so the equilibrium spacing may shift.

## Also Test: Phase Coherence (Time Crystal)

8. At equilibrium: are all oscillons breathing at the SAME frequency?
   (Expected yes, from V24-DG's finding that triple product enforces
   phase locking)
9. Is there a definite phase relationship between oscillons?
10. If yes: this IS a spatiotemporal crystal — periodic in both space
    and time.

## Reference Code

- v23/phonon/src/chain1d.c (periodic chain code)
- v24/maxwell_e (pairwise coupling)

## Output

- `src/lattice.c`, `data/`, `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0
λ scan: {0.0 (control), 0.5, 0.9}
N_osc=8, d=16, periodic BC
Nx = N_osc * 320 per spacing
t_equil=5000 (single), t_chain=10000

Compile: `gcc -O3 -Wall -o lattice src/lattice.c -lm`
