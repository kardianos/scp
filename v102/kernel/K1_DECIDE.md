# v102 kernel step, K1 — decision (the reference simulator)

**Written:** 2026-10-08, after `cas/out/k1_*.out`. Gates:
`K1_PREREG.md`.

## Ending: **Go**

Formation is one quantum jump, and the kernel can implement it as an
event. Binding needs a sink, and cold light at the boundary is one. Its
rate is now known, and it is very slow.

## K1-1. Formation as a quantum jump: **PASS**

Leading-order n = 4 lock, resonant photons, `cas/out/k1_form_M*.out`:

| fabric | holes made per photon | fidelity, jump against exact | space-lead flux = absorption | weight in the dressed nucleon's ground multiplet |
|---|---:|---:|---|---|
| 3+2 | 12.9% | **1.0000** | 1.421548 = 1.421548 | 1.000 (jump) / 1.000 (exact) |
| W | 5.6% | **1.0000** | 0.617187 = 0.617187 | 1.000 / 1.000 |

- **Why the jump is exact here.** The space band (±2t_s = ±0.12) is
  ten times wider than the lock's internal level spread. So absorption
  does not select among the depleted lock's states: it removes one
  quantum, and the space quantum takes the energy difference.
- The jump carries the photon's (S^z, I^z) by construction.
- The ledger closes. The incoming energy (lock ground + photon kinetic
  ε_res) equals the depleted lock's energy plus the leaving space
  quantum's, which sits inside its band.

## K1-2. Two holes in an isolated lock: **as registered (PASS)**

`cas/out/k1_closed_*.out`; T = 200/J0; norm drift ≤ 7e-12.

| lock | jumps at | bound-pair weight at t = 0 | time mean | energy above the bound pair |
|---|---|---:|---:|---:|
| 3+2 ring, n = 10 | 0, 3 | 0.048 | 0.048 | 0.27 J0 |
| 3+2 ring, n = 10 | 0, 5 | 0.033 | 0.033 | 0.29 J0 |
| W ring, n = 6 | 0, 3 | 0.163 | 0.163 | 0.91 J0 |

**A closed lock never binds two holes.** The binding energy has nowhere
to go: the lock's own quantum costs more than the binding.

## K1-3. Capture by cold light (a report)

3+2 ring, n = 6, two mobile holes. The start is an unbound S = 1
two-hole state 0.356 J0 above the pair. The bound pair is the three
S = 0 states at separation 2.99. `cas/out/k1_capture_d*.out`:

| cold source δ | in J0 | final bound-pair population | photons to ≥ 0.9 | a photon un-binds the pair |
|---:|---:|---:|---:|---:|
| 0.01 | 1.04 | 0.25 (stationary) | never | 1e-8 |
| 0.002 | 0.21 | **1.00** | ≤ 3.4×10¹⁰ | 1e-16 |
| 0.0005 | 0.05 | **1.00** | ≤ 3.4×10¹⁰ | 0 |

- **Cold light does capture.** Light colder than the pair's level
  structure brings two holes into their bound pair, and never back out.
  The pattern is radiative capture: binding needs an emitted quantum to
  carry the energy.
- **It is very slow:** about 10¹⁰ photons per pair. A cold photon is
  far off the boundary resonance, so it scatters only through a weak
  virtual exchange.
- **W was not run.** The six-cell two-hole W lock is too large for the
  full S-matrix (about 2.4 million states); a deviation. W's pair binds
  four times deeper and its lock gap is larger, so the same mechanism
  applies.

## K1-4. The classical field as quanta

Stated in `K1_PREREG.md`, not tested. K2 uses it as written.

## What K2 must therefore carry

1. **Formation:** an event at a lock boundary cell.
   - The trigger: a neighbouring field atom within the resonance
     window.
   - The rate: the K1-1 probability per resonant atom.
   - The action: the jump Σ_μ P[μ, ν] b_iμ on the lock state, one atom
     removed from the field, and Ω plus the detuning added to the
     neighbour's `Es`.
2. **No binding inside the lock without a sink.**
3. **Capture:** a boundary scattering rule for cold light. Its rate per
   photon is about 10⁻¹⁰ or less (K1-3), so at kernel light levels a
   pair forms only over very long times. It is reported, not added
   yet: K2a starts without it, and adding it is a separate decision.
