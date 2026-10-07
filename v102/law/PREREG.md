# v102 law step — pre-registration

**Written:** 2026-10-07, after `../phase1/DECIDE.md` and before any
number in this directory.
**Authority:** user, this session: "Commit and continue to 'Next'."
`../phase1/DECIDE.md` §5 named three items:

1. fix the law constants;
2. decide what "forced" means (the user's call; not done here);
3. re-register the two defective bars as a new experiment.

This page covers 1 and 3.

---

## Part L — does the law fix t, g, and Δ?

The derived lock (Phase 1) needs all three:

- **Δ > 0:** pair conversion into space is off-resonant, and uphill.
- **t/Δ < t_c(g/Δ):** about 0.16–0.26 for g/Δ = 0.05–0.2
  (`../phase1/cas/out/p1_tc.out`).
- g small enough for the spin-chain picture.

**Sources allowed:** `v89/PRINCIPLE.md`; the laws_V2g/V3a table and
its prose (`v89/battery/laws_V2g.cfg`); the rulings R1, R2, R5, R6;
and the kernel only as the table's operational definition. Nothing
before v89.

**Each constant is classed** as one of:

- **fixed:** a table constant or a principle sentence gives it, with
  no reading;
- **reading:** it follows only under a named interpretation (each
  reading is listed and computed);
- **undefined:** the table has no counterpart.

**Endings** (exactly one is claimed in `DECIDE.md`):

| ending | condition | consequence |
|---|---|---|
| A | t, g, Δ all fixed, inside the lock phase | ξ is predicted from the Phase 1 P5 table, with its spread |
| B | fixed, or fixed under every listed reading, and outside the lock phase | under the current table the derived lock does not exist |
| C | at least one is undefined, and no listed reading puts the fiber inside the phase | name the constant and the sentence that would fix it; state what the table's constants imply if the fiber inherits them |
| D | the table's conversion cannot be put in Hamiltonian form at all | the derived model needs its own conversion law |

**Registered expectation (not a gate):** g undefined (the table's
conversion is a rate law); Δ reading-dependent; the 3-D transport
constant much larger than any detuning the rungs provide. That
expectation points at C.

## Part R — the two defective bars, re-registered on new cases

These are new experiments on cases not computed before this page.

- **R-P1′.** Microscopic ring n = 6 (Σc ≤ 2, S^z = 0), g/Δ = 0.1,
  t/Δ ∈ {0.05, 0.15}. Neither point was computed in Phase 1.
  Compare the first 12 excitation energies *with multiplicity* to the
  spin-1 Heisenberg ring, with one fitted `J_eff`.
  - Pass: every miss < 5%, ground state total spin 0.
  - The t/Δ = 0.15 point lies inside the lock phase (n = 6 threshold
    in (0.18, 0.20]), so it is a fair test of the mapping and not of
    the phase boundary.
- **R-P6c′.** A new slit geometry: separation ±16, half-width 2, a
  different source. The eraser intensity equals the coherent sum of
  the single-slit fields, relative difference < 1e-10. The both-open
  scalar field is reported, not gated (it includes scattering between
  the slits).

## What does not count

- Choosing a reading because it lands inside the lock phase.
- Tuning `field_J`, `w1`, `w2`, `q_detune`, or `cap`. They are the
  table's and are read, not moved.
- Reading ξ as a cell count.
