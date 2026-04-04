# SCP Field Theory

A classical field theory based on three real scalar fields with a
triple-product potential, extended by three angle fields (Cosserat
coupling). The theory produces stable particle-like structures,
gravitational attraction, charge-dependent electromagnetic forces,
composite baryons with confinement, and nuclear binding.

All results are obtained by numerical simulation of the equations
of motion on a 3D grid. The equations themselves evolve across
experiments as the theory is refined -- each version tests specific
hypotheses about the field dynamics.

## Documents

- [CONCEPT.md](CONCEPT.md) -- The theory as currently understood.
  Written as a cohesive description, not a lab notebook.
- [FUTURE.md](FUTURE.md) -- Open questions, proposed experiments,
  and resolved or abandoned directions.
- [EM_THEORY.md](EM_THEORY.md) -- Electromagnetic sector in detail.

## Code

The simulation kernel, analysis tools, seed generators, and volume
viewer are in [sfa/](sfa/). The runner MCP server manages local and
remote (GPU) execution. See [CLAUDE.md](CLAUDE.md) for build
instructions and conventions.

## Simulation History

Versioned experiment directories (v28/, v34/, ..., v52/) contain
per-run plans, results, analysis, and generated data. Each version
tests specific physics questions -- equations, parameters, and
initial conditions vary between experiments by design.
