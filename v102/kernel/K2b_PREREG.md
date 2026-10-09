# v102 kernel step, K2b — pre-registration (the widened field)

**Written:** 2026-10-08, before any K2b code ran. Plan: `PLAN.md`.
Previous step: `K2a_DECIDE.md` (Go; §5 is why K2b is needed).
**Authority:** user, 2026-10-08: "Continue after that if it makes sense
to do so, committing between each step."

## 1. What K2b builds

**The widening.** The kernel field becomes `fmodes` complex amplitudes
per cell, ψ_ν with ν = 0 … F−1.

- **F = 1** is the default and today's field: ψ_0 = fa1 + i·fa2.
- **F = 2** gives the 3+2 modes: spin up, spin down.
- **F = 4** gives the W modes, in the order of `graphm.MODES_W`:
  (+,+), (+,−), (−,+), (−,−) in (2S^z, 2I^z).
- The extra modes live in new arrays, allocated only when F > 1.

**The field rules are mode-blind.**

- **Transport.** Every mode has the same on-site rotation (w1e·dt)
  and the same hop exp(−iτσ_x) on every link (pass F).
- **Energy.** Ee = Σ_ν |ψ_ν|², summed pairwise for F = 4:
  (m0 + m1) + (m2 + m3).
- **Conversions that take field energy** (condensation, `field_take`,
  `pattern_cond`, `pattern_evap`) scale every mode by the same factor,
  so a cell's mode content is kept.
- **Emission into a cell that has field** (`field_inject`) also scales
  every mode by the same factor.
- **Emission into an empty cell** (|ψ|² ≤ 1e-20) puts amplitude
  √(dE/F) at the cell phase th2 into **every** mode. Classical matter
  carries no spin or isospin, so the emitted light favours no mode.
- **Not widened.** Mode-0-specific options are refused when F > 1:
  `amp_door`, `qp_phase`, `door_g`, `bath_xb`, `bh_r`. The FCS snapshot
  and the loop-phase meters still read mode 0 only; this is stated as a
  limit.

**Seeding.** `fseed_mode` acts after all seeds and before the ledger
`E0`:

- `fseed_mode` = m ≥ 0 moves the seeded field (written to mode 0 by
  every seed) into mode m;
- `fseed_mode` = −1 splits it equally: ψ_ν = ψ_0/√F, which is exact
  (ψ_0/2) for F = 4.

**The lock under F > 1** (`lockq`, with `lockq_M` = F required).

- **Candidates.** A candidate is now (i, o, ν), with n_ν = floor(|ψ_ν(o)|²/ε_o).
- **The jump** is J_ν = Σ_μ P[μ, ν] b_iμ, into the sector shifted by
  mode ν's charges. One next sector is cached per mode.
- **The event** removes one atom from mode ν of o: ψ_ν(o) is scaled,
  and Ee_o drops by ε_o.
- **Window, rate and draw** are as in K2a after its correction: the
  exact per-jump energy change, the golden rule, and a Poisson draw.
- **At F = 1** the module behaves exactly as K2a.

## 2. Gates (all must pass before the K2b commit)

| | gate | bar |
|---|---|---|
| **B1** | **Inert at default.** The full battery, both kernels; kernel stdout for `exp=bath`, `blob`, `slit` and `pulse` (T = 20), against the K2a binary | 87 bars ALL GREEN, bar table identical; stdout byte-identical |
| **B2** | **Exact reduction.** F = 4 with `fseed_mode=-1`, against F = 1, for `exp=pulse` (amp 1, σ 2.2) and `exp=blob` (L = 16), T = 20, V3a law | every `t=` row and every `# RESULT` line byte-identical. (F = 2 with the split is not exact, since 1/√2 rounds; its largest relative deviation in the `t=` rows is reported, not gated) |
| **B3** | **Mode permutation is a symmetry.** F = 4, `exp=pulse` (amp 1, σ 2.2), `fseed_mode` = 0, 1, 2, 3 | every `t=` row and `# RESULT` line byte-identical across the four. The energy outside the seeded mode (from empty-cell emission) is reported |
| **B4** | **Per-mode events.** `lockq` on, F = `lockq_M` (3+2 n = 8, F = 2; W n = 6, F = 4), `exp=pulse` amp 4 σ 3, light seeded in each mode ν in turn, three module seeds | every event takes exactly one atom from mode ν, and moves the lock by mode ν's charges (an in-kernel check counts violations: 0); \|drift_rel\| ≤ 1e-12; no event outside the band. Also run with `fseed_mode=-1`: the hole charges made are reported, and whether a zero-charge two-hole sector (the nucleus's sector) is reached |
| **B5** | **Per-mode A/B.** A seeded jump with photon mode ν (W ν = 1, 2, 3; 3+2 ν = 1), 2000 steps, `lockq_rate=0` | hole, S^z and I^z profiles against `cas/k2a_ab.py` (extended with ν), to 1e-6 |

## 3. Not in K2b

- No MPS lock (K2c).
- No capture rule.
- No change to any law constant.
- `v97/fab` is not touched: the Go kernel stays single-mode, and the
  battery's C ↔ Go checks run at F = 1.
