# v102 Phase 1 — decision

**Written:** 2026-10-07, after the outputs in `cas/out/`.
Bars: `PREREG.md`. Rulings: R1, R2, R5, R6 (`../README.md`).
**Status:** this page closes Phase 1. It does not name Phase 2 as
passed. It does not edit a kernel.

---

## 1. Instrument (P2): **validated**

The DMRG in `cas/dmrg.py` was written for this phase and checked on
the effective chain: a spin-1 Heisenberg chain with spin-½ end cells,
which is what R6 makes at a degree-1 end.

| quantity | measured | literature | miss | bar |
|---|---:|---:|---:|---:|
| energy per bulk bond | −1.401484 | −1.401484 | 3.4e-8 | 1e-4 |
| gap, extrapolated in L | 0.41317 | 0.41050 | 0.65% | 2% |
| ξ, r ∈ [6, 24], L = 80 | 6.137 | 6.03 | 1.8% | 5% |

`cas/out/p2_validate.out`.

## 2. Results by bar

### P1. Mapping, with transport: **MISS as registered**

The miss has two separate causes.

**(a) A defect in the bar's design (n = 4).** The 4-site Heisenberg
ring is `(S₁+S₃)·(S₂+S₄)`, so it has large accidental degeneracies.
The microscopic model splits them at fourth order, which misaligns
the comparison of the "first 8 distinct levels." Compared with
multiplicity instead (`p1_diag.out`; a post-registration diagnostic,
labelled as such), the miss is 2.1% at t/Δ = 0, 3.1% at 0.05, 7.1% at
0.1, and 21% at 0.15.

**(b) Physics: the lock is a phase, bounded by transport.** At
t/Δ = 0.2 the ground state is no longer the antiferromagnetic singlet.

| ring | t/Δ = 0.2 ground state | ⟨space⟩ |
|---|---|---:|
| n = 4 | S = 3 | 0.86 |
| n = 6 | S = 5 | 0.98 |

The conversion holes become mobile, and a ferromagnetic background
lets them hop freely (Nagaoka-type). The estimate: each hole hops with
amplitude `2t`, so the converted state sits near `Δ − 8t`. Measured
thresholds (`p1_tc.out`, `p1_tc6.out`):

| ring | g/Δ | t_c/Δ |
|---|---:|---:|
| 4 | 0.05 | 0.16 |
| 4 | 0.10 | 0.19 |
| 4 | 0.20 | 0.26 |
| 6 | 0.10 | in (0.18, 0.20] |

Inside the phase, the n = 6 ring passes: 1.9% at t/Δ = 0 and 1.8% at
t/Δ = 0.1, with a singlet ground state (`p1_micro.out`).

### P3. ξ from the microscopic chain: **MISS** (by a replacement instrument)

**The registered route did not converge.** The microscopic DMRG gave a
higher energy at χ = 48 than at χ = 32, so it is not trustworthy. A
first round also used a chemical potential that made space free, which
put it in the wrong sector. Both sets are kept, unused, in
`cas/out/p3_bad_mu05/` and `cas/out/p3_unconverged/`.

**The replacement.** The exact cluster effective Hamiltonian
(des Cloizeaux projection of 2-cell and 3-cell microscopic clusters,
exact through fourth order), followed by the validated DMRG on the
resulting chain, with L = 48 and r ∈ [6, 20] (`p3_cluster.out`):

| g/Δ | t/Δ | J / (2g²/Δ) | K/J (biquadratic) | three-site remainder / J | ξ |
|---:|---:|---:|---:|---:|---:|
| 0.05 | 0 | 0.990 | −0.00001 | 0.024 | 6.06 |
| 0.10 | 0 | 0.963 | −0.00016 | 0.086 | **5.79** |
| 0.20 | 0 | 0.877 | −0.0016 | 0.243 | 5.07 |
| 0.10 | 0.10 | 0.999 | +0.0355 | 0.099 | 5.97 |
| 0.20 | 0.20 | 0.982 | +0.105 | 0.391 | 5.22 |
| 0.30 | 0 | 0.781 | −0.0047 | 0.376 | 4.39 |
| Heisenberg reference (same L, window) | | | | | 6.16 |

At (0.2, 0.2) and (0.3, 0), the three-site remainder is about 0.4 J,
and the 3-cell gap to the space-carrying states falls to about 0.3 Δ.
Fourth-order truncation is not reliable there. Those two rows show the
trend only.

The single-link interaction stays almost pure Heisenberg, and its
SU(2) residual is about 1e-16. ξ moves through the genuinely two-link
(three-site) terms, which grow as (g/Δ)². At the registered point
(0.1, 0), ξ = 5.79 is 5.6% below P2's 6.137 and 6.0% below the
same-window reference (6.16), against a 5% bar. At g/Δ = 0.05 it
would have passed (1.3%).

### P4. One line, no second line: **PASS**

- (a) The first level carrying a real space quantum is index 141,
  right after all 141 spin-manifold states. It lies at 22.9× twice the
  magnon gap.
- (b) On rings of n = 8, 10, 12, the lowest S = 0 excitation is at
  2.26, 2.14, 1.94 × 2Δₙ, and the lowest S = 2 state at 1.47, 1.43,
  1.38 × 2Δₙ. On the open chain, `E(S^z = 2)` sits at 1.075, 1.043,
  1.028 × twice the gap.

There is no bound pair, so no second species, below the two-magnon
threshold.

### P5. Corrections: reported (the table in P3)

ξ depends continuously on g/Δ (roughly quadratically) and, more
weakly, on t/Δ. The single-link `J` approaches `2g²/Δ` as g → 0.

### P6. The slit: **(a), (b) PASS; (c) MISS as registered**

- (a) The two-component field equals the scalar field to 4e-16.
- (b) With orthogonal polarizers on the two slits, the interference
  fraction is 5e-16.
- (c) The eraser's visibility differs from the both-open scalar's by
  6.4e-3, against a registered 1e-10.

**(c) is a defect in the bar's design.** The eraser equals the
coherent sum of the two single-slit fields, to 4e-16. The both-open
field is not that sum: waves scatter between the slits, a 3.2e-3
effect (`p6_slit.out`). The bar compared against the wrong object.

**Consequence:** a lone quantum cannot convert. The conversion
operator needs two quanta, one per link end (W1). Detecting a click
therefore needs matter at the screen to pair with. That is a new open
regression item. The fringe loci carry over exactly.

## 3. Decision

**Phase 1 stops.** Under the advance rule, P1, P3, and P6 each miss.

| bar | how it missed | design defect or physics? |
|---|---|---|
| P1 at n = 4, t ≤ 0.1 | accidental degeneracy misaligns the levels | **design defect**: the bar compared distinct levels on a degenerate ring |
| P6(c) | compared with the both-open field | **design defect**: wrong reference |
| P1 at t/Δ = 0.2 | the ground state is a converted ferromagnet | **physics**: the lock has a transport threshold |
| P3 at (0.1, 0) | ξ = 5.79, 5.6% below 6.137 | **physics**: three-site terms move ξ at order (g/Δ)² |

## 4. What is now true, and what is not

**True:**

1. The derived lock exists only as a phase. Transport must be weak
   against the conversion detuning, with `t/Δ` below about 0.16–0.26
   depending on g/Δ. Beyond that, conversion holes hop, the field
   converts, and the ground state is ferromagnetic.
2. Inside the phase, the fiber is a spin-1 Haldane chain:
   - singlet ground state;
   - one isolated triplet line;
   - no bound pair;
   - the first space-carrying level far above.
3. ξ ≈ 6 (in grains) is the leading-order value. At finite `g/Δ` and
   `t/Δ` it moves continuously: 6.06, 5.79, and 5.07 at g/Δ = 0.05,
   0.1, and 0.2. It is forced only in the limit of weak conversion and
   weak transport.
4. The linear slit is recovered exactly, component by component.
   Detection needs a partner quantum.

**Not true:**

- ξ is forced at finite coupling. It is a function of `g/Δ` and
  `t/Δ`, the law constants R1 and R2 left unspecified.
- `n = 6`, or any ring size, is selected. ξ is a decay length, and
  `E₀(n)` stays concave.
- The microscopic DMRG measured anything. It did not converge.
- Phase 2 is named.

## 5. What would move past this page

*Update, same day:* item 1 was done as a new experiment on new
cases, and both bars passed (`../law/DECIDE.md` Part R). Item 3 was
examined in `../law/`, which ended at C.


1. **Re-register the two defective bars, as a new experiment.**
   - P1: compare with multiplicity, inside the lock phase.
   - P6(c): compare against the coherent single-slit sum.

   Neither would rescue P3.
2. **Decide what "forced" means.** ξ is a leading-order number plus
   corrections that depend continuously on two law ratios. That is
   weaker than "no continuous coefficient," and stronger than v98's
   `ρ`, which moved `n*` at order one. Whether that counts as a size
   law is the user's call (`v99/reopen/OPTIONS.md` C).
3. **Fix the law constants.** R1 and R2 change the field and its
   transport, but the law table does not yet give `g`, `Δ`, or `t`.
   The lock needs `t/Δ` below threshold. Whether the law supplies that
   is the next real question, before any Phase 2.
