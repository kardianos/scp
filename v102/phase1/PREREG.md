# v102 Phase 1 — pre-registration

**Written:** 2026-10-07, after R6 was ruled and before any Phase 1
code or number. Phase 1 opened under R1, R2, R5, and R6
(`../README.md`).

---

## 1. The model

Phase 1 is an instrument. The calculator is not the fabric.

**Microscopic chain.** Cells and links alternate.

- **A cell** holds two bosonic modes `b_{i,1}`, `b_{i,2}`, with total
  occupation capped at its number of link ends (R6). That cap is 2 in
  the bulk or on a ring, and 1 on an end cell of an open chain.
- **A link** carries one scalar space mode `c_l`, truncated at 2.
- **Hamiltonian:**

  ```
  H = Δ Σ_l c_l†c_l
    + g Σ_l (c_l† A_l + h.c.),     A_l = b_{i1}b_{j2} − b_{i2}b_{j1}
    + t Σ_l Σ_a (b_{i,a}† b_{j,a} + h.c.)
  ```

  Transport that would exceed a cell's cap has no matrix element.
- `g/Δ` and `t/Δ` are not fixed by the law. They are measured as
  sensitivities, never chosen.

**Effective chain.** The derived leading-order model: a spin-1
Heisenberg chain with `J = 2g²/Δ` (`../gates/WEDGE.md`). On an open
chain, its end cells are spin ½ (R6, degree 1).

## 2. Bars

**P1. Mapping, with transport** (exact diagonalization, microscopic
ring).

- n = 4 with the full capped space, and n = 6 with Σc ≤ 2.
- Couplings: g/Δ = 0.1, and t/Δ ∈ {0, 0.1, 0.2}.
- Pass requires all of:
  - ground state total spin 0, at full occupation;
  - with one fitted `J_eff`, every one of the first 8 distinct
    excitation energies within 5% of the spin-1 Heisenberg ring's.

**P2. Instrument validation** (DMRG, effective open chain with
spin-½ ends).

- Energy per bulk bond within 1e-4 of −1.401484.
- Gap, extrapolated in L, within 2% of 0.41050 J.
- ξ from bulk correlations within 5% of 6.03.

If P2 misses, Phase 1 is **inconclusive**, not failed.

**P3. ξ from the microscopic chain** (DMRG, cells plus links, degree-1
end cells).

- Points: (g/Δ, t/Δ) = (0.1, 0), (0.2, 0), (0.1, 0.1), (0.2, 0.2).
- Pass: at (0.1, 0), ξ_micro is within 5% of the P2 value.
- The dependence across the other points is reported as a table and
  is not used to pick anything.

**P4. One line, no second line.**

- (a) In the microscopic ring of P1 (n = 6, g/Δ = 0.1, t/Δ = 0.1),
  the lowest level carrying a real space quantum (⟨space⟩ > 0.5)
  lies above twice the magnon gap.
- (b) In the effective chain, the lowest S = 0 excitation and the
  lowest S = 2 state lie at or above 0.98 × 2Δ.
  - On rings of n = 8, 10, 12 (exact diagonalization).
  - On the open chain via DMRG: the lowest `S^z = 2` state against
    twice the lowest `S^z = 1` excitation.

  A level below that bound would be a second line, meaning a bound
  pair.

**P5. Corrections.** `J_eff(g/Δ)` and ξ_micro across the P3 points,
reported. Not a gate beyond P3's limit.

**P6. The slit, from the unlocked linearized field** (lattice
calculator).

- (a) The screen intensity of the two-component field equals the
  scalar field's, relative difference < 1e-10.
- (b) With orthogonal polarizers on the two slits, fringe visibility
  < 0.01.
- (c) Projecting onto the diagonal polarization (the eraser) restores
  the scalar visibility to within 1e-10.

## 3. Advance

Phase 1 **passes** if P2 validates the instrument and P1, P3, P4, and
P6 pass. A pass names Phase 2 (the same chain on cells that are born
and die) without authorizing it. Any miss on P1, P3, P4, or P6 stops
Phase 1, and `DECIDE.md` names it.

## 4. What does not count

- Fitting `g/Δ` or `t/Δ` so that ξ lands anywhere.
- Reading ξ as a cell count, or taking `n = 6` from it.
- Treating the DMRG chain or the slit lattice as the substance.
- Using P2's literature values as anything but a check on the
  instrument.
