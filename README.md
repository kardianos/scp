# SCP Field Theory

A classical field theory based on three real scalar fields with a
triple-product potential, extended by three massless angle fields
(Cosserat coupling). The theory produces stable particle-like
structures (braids), gravitational attraction, charge-dependent
electromagnetic forces, composite baryons with confinement, and
nuclear binding from a single Lagrangian.

All results are obtained by numerical simulation of the equation
of motion on a 3D grid, with no fitting to experimental data.
The parameters (m, mu, kappa, eta) are fixed across all experiments.

## Documents

[CONCEPT.md](CONCEPT.md) describes the current state of the theory,
including the field equation, particle structure, gravity mechanism,
electromagnetism, composite particles, nuclear binding, and all
confirmed and open results.

[DISCOVERIES.md](DISCOVERIES.md) is the chronological record of
simulation results from V23 through V42, including confirmed
findings, negative results, and key numerical measurements.

[FUTURE.md](FUTURE.md) lists open questions, proposed experiments,
and resolved or abandoned directions. Items are prioritized by
their impact on the viability of the theory.

[EM_THEORY.md](EM_THEORY.md) contains the electromagnetic sector
in detail: the mapping of angle fields to the vector potential,
the derivation of Maxwell's equations from the linearized Cosserat
equation, the mechanism by which oscillating magnetic dipoles
produce an effective Coulomb force through radiation scattering,
and open questions regarding the fine structure constant and
composite charge structure.

## Code

The simulation kernel, analysis tools, seed generators, and volume
viewer are in the [sfa/](sfa/) directory. See [sfa/README.md](sfa/README.md)
for build instructions and usage. The Makefile at [sfa/Makefile](sfa/Makefile)
builds all targets.

## Simulation History

Versioned experiment directories (v28 through v42) contain per-run
plans, results, analysis, and generated data. Each version builds
on the findings of the previous.
