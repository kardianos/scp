# V23-D: Phonon-Graviton from Oscillon Gas

## Thesis

Individual oscillons interact via Yukawa tails (short-range, massive). But a
GAS or LATTICE of many oscillons has collective excitations — phonons — that
can be massless (acoustic branch with ω = c_s k). If these phonons exist and
couple universally to energy density, they are the graviton analog.

This investigation creates a 1D chain of oscillons, measures the collective
mode spectrum, and checks for a gapless acoustic branch.

## Background

In condensed matter: individual atoms interact via short-range potentials
(Lennard-Jones, Yukawa). The atomic lattice has acoustic phonons that are
automatically massless (ω → 0 as k → 0). The speed of sound c_s = √(K/ρ)
where K is the bulk modulus and ρ is the mass density.

For oscillons:
- Individual mass: E ~ 1.3 (μ=-20) to ~6.3 (μ=-5) code units
- Interaction: Yukawa-like with range ξ ~ 5-11 (from V23-C)
- Equilibrium spacing: unknown — need to find it
- Phase dependence: in-phase → attractive, anti-phase → repulsive (from V22)

## Mathematical Setup

### 1D Oscillon Chain

Place N oscillons at positions x_n = n·d (spacing d) on a 1D grid. Each
oscillon consists of three fields φ₁, φ₂, φ₃ initialized from the
equilibrated 1D profile.

The collective coordinate for oscillon n is its center-of-energy position:

    X_n(t) = ∫ x·ρ_n(x,t) dx / ∫ ρ_n(x,t) dx

where ρ_n is the energy density attributed to oscillon n.

### Phonon Spectrum

For small displacements δX_n = X_n - n·d from equilibrium:

    M·d²(δX_n)/dt² = -K(2·δX_n - δX_{n-1} - δX_{n+1})

where M is the oscillon mass and K = -V''(d) is the effective spring constant
from the inter-oscillon potential V(d).

This gives the phonon dispersion:

    ω²(k) = (4K/M)·sin²(kd/2)

which is GAPLESS: ω → c_s·k as k → 0, with c_s = d·√(K/M).

### What Could Go Wrong

1. **No equilibrium spacing**: if oscillons always merge (at all d) or always
   repel (at all d), there's no lattice and no phonons.
2. **Phase decoherence**: if breathing phases randomize, the in-phase/anti-phase
   interaction averages to zero or net repulsion → no lattice.
3. **Instability**: the chain might be dynamically unstable (Jeans instability)
   even if a static equilibrium exists.

## What to Compute

### Phase 1: Inter-Oscillon Potential (prerequisite)

1. Equilibrate a SINGLE oscillon (μ=-20, κ=20, m=1.0) for t=5000.
   Save the profile f_eq(x).

2. Initialize TWO oscillons from f_eq at separations D = 8, 10, 12, 14, 16,
   18, 20, 25, 30 (both in-phase, i.e., same breathing phase).

3. Evolve for t=2000 (short — we want to measure the initial force before
   significant displacement occurs).

4. Measure the acceleration of the separation: a(D) = d²(sep)/dt² at t ~ 100.
   This gives the force F(D) = M·a(D)/2.

5. Integrate F(D) to get the inter-oscillon potential V(D).

6. Find the equilibrium spacing d_eq where V'(d_eq) = 0 (if it exists).
   Compute K = V''(d_eq).

### Phase 2: Oscillon Chain (if Phase 1 finds equilibrium)

7. Place N = 8-16 oscillons at spacing d_eq, initialized from f_eq.

8. Add small random displacements δX_n ~ 0.1 to each oscillon position.

9. Evolve for t=10000.

10. Track X_n(t) for each oscillon.

11. Compute the Fourier transform of X_n(t) in both space (k) and time (ω)
    to get the phonon dispersion ω(k).

12. Check: is there a gapless acoustic branch (ω → 0 as k → 0)?

### Phase 3: Sound Speed and G_eff

13. From the dispersion: extract c_s (slope of ω(k) at k=0).

14. Compute the effective gravitational constant:
    G_eff = 1/(ρ_chain · c_s²)
    where ρ_chain = M/d_eq is the linear mass density.

15. Compare with the microscopic interaction parameters.

## Reference Code

- v21 1D solver: `/home/d/code/scp/v21/src/triad1d.c` (base code)
- v22 two-oscillon: `/home/d/code/scp/v22/src/two_oscillon.c` (centroid tracking)
- v23/critical: `/home/d/code/scp/v23/critical/src/critical1d.c`
  (Phase 2 two-oscillon code — adapt for N oscillons)

## Output Structure

- `src/phonon1d.c` — main code (Phases 1-3)
- `data/` — output data files
- `RESULTS.md` — results and analysis

## Parameters

μ=-20, κ=20, m=1.0 (v21 production, strongest coupling for clearest signal)
Grid: Nx=8000, xmax=200 (need room for N=8-16 oscillons)
Phase 1: t_equil=5000 (single), t_measure=2000 (pair)
Phase 2: t_chain=10000

Start with Phase 1. Only proceed to Phase 2 if an equilibrium spacing exists.

Compile: `gcc -O3 -Wall -o phonon1d src/phonon1d.c -lm`
