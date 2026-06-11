# SCP Field Theory

A classical field theory of three complex scalar fields with a
triple-product potential, three Cosserat angle fields, and a gauged
diagonal U(1) symmetry. The theory supports absolutely stable charged
particles (Q-balls) with quark-like internal substructure, a measured
massless Coulomb force, a short-range phase-coherent nuclear force,
fused composite nuclei with mass defects, a stable flavor multiplet,
and collective nuclear modes (a giant-dipole-resonance analog).

All results are obtained by numerical simulation of the equations of
motion on a 3D grid, with conserved-quantity bookkeeping at machine
precision. The structural parallels to particle physics are measured;
no quantitative mapping to real physical constants has been
established.

Earlier formulations used real fields (v28–v53); those produced only
long-lived breathing structures (oscillons, braids) that are now
understood to be necessarily mortal — the U(1) complexification
(v66) is what first yields true particles. See CONCEPT.md §2 and §9.

## Documents

- [CONCEPT.md](CONCEPT.md) -- The theory as currently understood.
  Written as a cohesive description, not a lab notebook.
- [DISCOVERIES.md](DISCOVERIES.md) -- Chronological lab notebook.
- [FUTURE.md](FUTURE.md) -- Open questions, proposed experiments,
  and resolved or abandoned directions.
- [EM_THEORY.md](EM_THEORY.md) -- Electromagnetic sector (historical
  θ-polariton analysis; the current photon analog is the gauge field).

## Code

The simulation kernel (real / complex / gauged modes in one
config-driven binary), analysis tools, seed generators, and volume
viewer are in [sfa/](sfa/). The runner MCP server manages local and
remote (GPU) execution. See [CLAUDE.md](CLAUDE.md) for build
instructions and conventions.

## Simulation History

Versioned experiment directories (v28/ ... v71/) contain per-run
plans, results, analysis, and generated data. Each version tests
specific physics questions -- equations, parameters, and initial
conditions vary between experiments by design. Key eras: real-field
phenomenology (v28–v53), algebraic side-tracks (v54–v63), stability
crisis (v64–v65), U(1) Q-balls (v66–v68), gauged kernel-v3 (v69),
verification + substructure/nuclei/flavor (v70–v71).
