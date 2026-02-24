# V12 Conclusion: The Non-Linear Strong Geon

## 1. Theory

V12 tested whether combining three previously isolated mechanisms could
produce a stable, proton-mass electromagnetic geon:

1. **Binding force** (V9): Yukawa massive tensor from f_2(1270), providing
   ~1 fm confinement range with coupling kappa=35.

2. **Energy cap** (V11 + new): Born-Infeld nonlinearity, reformulated in
   the D-field representation to prevent the 1/A runaway that creates
   the Q9 horizon paradox.

3. **Topology** (V11): Hopf magnetic knot, providing topological
   confinement without an artificial Rmax boundary.

The central claim was that BI nonlinearity cuts off the divergence that
causes the Q9 paradox, allowing a self-consistent equilibrium at
proton-mass energies without horizon formation.

## 2. Method

### Phase 1: BI-Modified 1D Hartree (scalar BI, V9 code)

Used V9's existing `strong_geon` solver with the `-bBI` parameter at
physical coupling (kappa=35, mu=6.47, Pade metric). Scanned b_BI from
0.01 to 1e10 at Rmax values spanning the Q9 collapse transition.

No code modifications required — V9 already implemented the scalar BI
formula. See IMPLEMENTATION.md Section "Phase 1" for equation traceability.

### Phase 3: Full 3D BI-Hopfion in Yukawa Metric (new code)

Built `v12/src/bi_geon.c` combining V11's 3D BI evolution with V9's
Yukawa gravity. Three key innovations over V11:

1. **D-field BI formulation**: Analytical constitutive relations
   E = dH/dD, H = dH/dB — no Newton solver, ~10x faster. The D-field
   Hamiltonian rho_BI = b^2(sqrt(R)-1) is finite for ALL field strengths,
   eliminating the NaN crashes that would occur with E-field formulation
   when E approaches b.

2. **FFT Yukawa solver**: 3D Yukawa equation (nabla^2 - mu^2)Phi = -kappa*rho
   solved via Green's function in Fourier space: Phi(k) = -kappa*rho(k)/(k^2+mu^2).
   O(N^3 log N) per metric update, amortized over multiple timesteps.

3. **ADM evolution with correct lensing**: dB/dt = -curl(alpha*E),
   dD/dt = curl(alpha*H), where the lapse alpha goes INSIDE the curl.
   The product rule curl(alpha*E) = alpha*curl(E) + grad(alpha) x E
   generates the gravitational lensing term that bends light toward the
   well center. The Pade lapse alpha = 1/sqrt(1-2Phi) is strictly positive
   for all Phi, preventing coordinate singularities.

Phase 2 (static Hopfion in Yukawa) was skipped per PLAN.md Section 7.

See IMPLEMENTATION.md for full equation-to-code traceability.

## 3. Results

### Phase 1: Q9 Resolution — POSITIVE

Born-Infeld nonlinearity **completely eliminates** the Q9 horizon collapse:

- Without BI: At Rmax=1.0 (kappa=35, mu=6.47), the Hartree iteration
  diverges to M_peak = 2e19, f_min = 0 (full horizon collapse).

- With BI (b_BI=0.5): At the same parameters, the iteration converges
  smoothly to f_min = 0.317 (68% well depth). No collapse.

- The Q9 sharp transition between "collapsed" and "bound" is replaced by
  a continuous crossover. f_min is a smooth, monotonic function of both
  b_BI and Rmax.

- The mechanism is exactly as predicted: BI caps the energy density growth
  from 1/A (divergent) to 1/sqrt(A) (sub-critical). The Hartree feedback
  loop saturates at finite well depth.

### Phase 3: 3D Hopfion Confinement — NEGATIVE at Physical Parameters

The Hopfion disperses at physical parameters (kappa=35, mu=6.47):

- Initial gravitational well: Phi_min = -0.233 (23% well depth)
- R_eff grows from 1.01 to 4.75 in T=5 (comparable to flat-space 4.8)
- Root cause: Yukawa range (1/mu = 0.15 fm) is 7x shorter than Hopfion
  core (a = 1 fm). Only 2.4% of the field volume feels significant pull.

At unphysical parameters (kappa=1000, mu=1.0), a quasi-stable state
IS observed:

- R_eff stabilizes at 1.42 from T=17 to T=20 (0.1% variation)
- B_peak oscillates between 0.1 and 3.1 (EM standing mode in well)
- Phi deepens from -42.8 to -62.4 (dissipative contraction from sponge)
- alpha_min = 0.089 > 0 (deep well, no horizon)

This demonstrates the confinement mechanism works in principle but
requires a mediator with mass ~200 MeV (not f_2(1270) at 1270 MeV)
and coupling ~1000 (not the physical 35).

### Numerical Validation

- Flat-space (gravity OFF) matches V11 exactly: E=0.584, disperses at c
- Energy conservation: dE/E = 3e-4 at N=64 (flat), degraded to ~0.5%/code-time
  with gravity (sponge + metric update artifacts)
- divB stable at 1.09e-2 throughout (4th-order truncation of initial data)
- Pade lapse stays positive for Phi as deep as -62 (no coordinate singularity)
- D-field BI formulation produces zero NaN across all runs

## 4. Verdict

**V12 proves the concept but fails at physical parameters.**

The BI energy cap resolves Q9 (Phase 1 positive), but the Hopfion cannot
be confined by the f_2(1270) Yukawa force (Phase 3 negative). The geometric
mismatch a*mu >> 1 is fatal: the attractive force decays exponentially
before reaching the Hopfion surface.

The quasi-stable state at kappa=1000, mu=1.0 shows that a sufficiently
long-range, strongly-coupled gravitational force CAN confine a BI Hopfion.
This is the first demonstration of a stable 3D electromagnetic geon in
the entire V8-V12 program. However, it requires parameters that do not
correspond to any known physical interaction.

### What V12 Adds to the Overall Picture

| Mechanism | Status | Fatal Flaw |
|-----------|--------|------------|
| V9 Yukawa alone | Collapses to horizon (Q9) | 1/A runaway |
| V11 BI alone | Disperses in flat space | No binding force |
| V12 BI + Yukawa (physical) | Prevents collapse but disperses | Yukawa range too short |
| V12 BI + Yukawa (unphysical) | **QUASI-STABLE GEON** | Parameters nonphysical |

The remaining path would require identifying a physical mechanism that
provides O(1 fm) range attraction with coupling ~10^3, which is not
available in the Standard Model below the electroweak scale.
