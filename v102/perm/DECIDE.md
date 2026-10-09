# v102 perm — decision (re-readings and decisions, alone and in combination)

**Written:** 2026-10-08, after the outputs in `cas/out/`. The factors,
the bars, the tiers and the algebra are in `PREREG.md`. Deviations are
listed in §8.
**Status:** a measurement in models. It rules nothing.

---

## 1. The answer

**One cell reaches the Full tier: W + T + S.**

- **W:** v98's own widening, 3+2 → 3+4. Each quantum carries spin ½
  and isospin ½.
- **T:** bounded transport inside the lock.
- **S:** space may leave the lock's boundary link.

In that cell:

- The fiber is a lock. Its gap extrapolates to 1.32 J0, three times
  3+2's.
- One open slot is a firm (½, ½) object: p↑, p↓, n↑, n↓. Both p and n
  have spin ½, and they are one object at one size: decay length
  2.37 → 2.30 grains between 24 and 32 cells.
- It walks under T, and light makes it under S (≈ 5% per resonant
  photon).

**One caveat stays on that "Full".** The size of the *moving* W hole
was not measured; its static size is used. In 3+2, the moving hole kept
a length-independent cloud at every transport tried. But at t = J0
(T's registered point) that cloud is not a single exponential, and it
fails the registered decay-length bar (§3, T).

**What no cell has:**

- **A ladder.** In every cell, every excited level above the hole is
  the hole plus one free lock quantum, a little above threshold, with a
  repulsive channel.
- **A nucleus.** Two holes either convert into a space quantum (t = 0)
  or collapse into one empty cell that cuts the fiber (t > 0).

**What each piece buys and costs** (§6 has the full list):

- **W** is the only re-reading that gives spin ½ and p/n together. It
  is a law change: four modes per cell instead of two.
- **T** buys walking. It costs firmness of shape as t rises, and it
  lets holes collapse in pairs.
- **S** buys formation. It is an exception to V at the boundary.
- **I** gives p/n but removes spin. Every I cell fails.

## 2. The twelve cells

Marks: ✔ pass; ✘ fail; (a) a verdict fixed by `PREREG.md` §3's
algebra and confirmed numerically; ✘* fails the registered decay-length
bar although the profiles at the two lengths agree to 1e-3 (§3, T);
✔p the static value carried as a proxy (§8.3). T is scored at its
registered t = J0, with t = 0.1 J0 in brackets where it differs.

| cell | K | N1 | N2 | N3 spin ½ | N4 | N5 ladder | N6 walks | N7 p and n | N8 forms | tier |
|---|---|---|---|---|---|---|---|---|---|---|
| R0 | ✔ | ✔ | ✔ | ✔ | ✔ | ✘ | ✘ (a) | ✘ (a) | ✘ | Partial |
| R0 + T | ✔ | ✔ | ✘* (✔) | ✔ | ✔ | ✘ | ✔ | ✘ (a) | ✘ | Fail (Partial) |
| R0 + S | ✔ | ✔ | ✔ | ✔ | ✔ | ✘ | ✘ (a) | ✘ (a) | ✔ | Partial |
| R0 + T + S | ✔ | ✔ | ✘* (✔) | ✔ | ✔ | ✘ | ✔ | ✘ (a) | ✔ | Fail (Partial: only N7) |
| I | ✔ | ✔ | ✔ | ✘ (a) | ✔ | ✘ | ✘ (a) | ✔ | ✘ | Fail |
| I + T | ✔ | ✔ | ✘* (✔) | ✘ (a) | ✔ | ✘ | ✔ | ✔ | ✘ | Fail |
| I + S | ✔ | ✔ | ✔ | ✘ (a) | ✔ | ✘ | ✘ (a) | ✔ | ✔ | Fail |
| I + T + S | ✔ | ✔ | ✘* (✔) | ✘ (a) | ✔ | ✘ | ✔ | ✔ | ✔ | Fail: everything but spin |
| W | ✔ | ✔ | ✔ | ✔ | ✔ | ✘ | ✘ (a) | ✔ | ✘ | Partial |
| W + T | ✔ | ✔ | ✔p | ✔ | ✔ | ✘ | ✔ | ✔ | ✘ | Partial (only N8) |
| W + S | ✔ | ✔ | ✔ | ✔ | ✔ | ✘ | ✘ (a) | ✔ | ✔ | Partial (only N6) |
| **W + T + S** | ✔ | ✔ | ✔p | ✔ | ✔ | ✘ | ✔ | ✔ | ✔ | **Full** (N2 by proxy) |

I + W is W, exactly (X6).

## 3. Each factor alone

### I: the fiber's SU(2) read as isospin (X6)

- Every number equals R0's (P1).
- The hole's two states are exactly degenerate. In the n = 8 effective
  ring the lowest six levels agree to 4e-16, and the profiles mirror to
  1e-14. That is one object with two occupancies (which plane holds
  the open slot), at the same size, with a p/n split of exactly 0.
  **N7 passes.**
- Nothing in the fiber is left to carry spin. **N3 fails.**

### W: the widening 3+4 (X1, X2, X6; `w_lock.out`, `analyze_tp.out`)

- **The pair.**
  - On two full cells, A†A takes the values {0, 6, 10} (3+2:
    {0, 4, 6}).
  - The first conversion still costs Δ = 0.479. Energies depend only
    on how many quanta a cell holds.
- **The lock (K).**
  - Rings n = 4 and 6: the ground state is a unique (0, 0) singlet.
    The first excitation is a 15-fold multiplet (1,1) ⊕ (1,0) ⊕ (0,1).
  - The microscopic n = 4 ring agrees with the effective one.
  - Long chains (TeNPy, χ = 128): the gap is 1.411, 1.365 and 1.347 J0
    at L_bulk = 16, 24, 32, **extrapolating to 1.32 J0**. 3+2's is
    0.42–0.44 J0.
  - At L = 24 the (1,0) member is degenerate with the (1,1) member to
    1e-7.
- **The hole.**
  - Its ground multiplet is (½, ½): four states, exactly degenerate.
  - Its cloud has the hole cell at −0.091 and the neighbours at +0.656
    (the same for S^z and I^z, to 1e-5), decaying fast.
  - **N2: the decay length is 2.365 → 2.299 grains between
    L_bulk = 24 and 32 (2.8%), and 2.600 → 2.600 on a common range.**
    The 90% radius is 5 grains at 16, 24 and 32.
  - With charge read as `I^z + ½`, N3 and N7 both pass.
  - The (3/2, 3/2) level sits above the hole-plus-quantum threshold:
    1.429 against 1.365 J0 (L_bulk = 24), and 1.387 against 1.347
    (32). The margin shrinks as about 1/L, the signature of scattering,
    not binding. **No ladder.**
- **I + W = W.** The spectra of the mirror sectors agree to 6e-16.

### T: transport inside the lock, so the hole can walk (X3)

Rings, band width of the hole's positional states at t = 0.1 J0:

| ring | 3+2 | W |
|---|---:|---:|
| n = 4 | — | 0.436 J0 |
| n = 6 | 0.145 J0 | 0.388 J0 |
| n = 8 | 0.123 J0 | — |
| n = 10 | 0.112 J0 | — |

- At t = 0 the band is exactly zero. At t = J0 it is wider than the
  threshold.
- **N6 passes in both fabrics, with the W hole about three times as
  mobile.**
- The ground state stays spin ½ (3+2) or (½, ½) (W) at every t tried.

Long chains, 3+2, TeNPy (χ = 64), L_bulk = 24 and 40:

| t | E(3/2) − E(½) against the gap | decay length of the hole-centred profile | change | profiles at 24 and 40 |
|---|---|---|---|---|
| 0.1 J0 | 0.512 against 0.479 (24): above | 5.14 → 5.32 | **3.4% ✔** | agree to 3e-3 |
| 1.0 J0 | 0.543 / 0.479 (24); 0.474 / 0.441 (40): above | 6.90 → 5.71 | **17% ✘** | agree to 4e-3 at every d ≤ 20 |
| 3.0 J0 | 0.528 / 0.479; 0.469 / 0.441: above | 3.29 → 5.69 | **73% ✘** | agree to 4e-3 |

- **No ladder under T at any t.**
- **The cloud changes character with t:**
  - at 0.1 J0 it is a staggered cloud like the static hole's;
  - at 1 J0 it is weaker and staggered;
  - at 3 J0 it is a uniform polarization, the precursor of the Nagaoka
    state.
- **The size result splits.** The registered decay-length fit fails at
  1 and 3 J0, but the hole-centred profiles at the two lengths are the
  same to three decimals out to d = 20. The cloud does not grow with
  the chain. It stops being a single exponential, so a fit over each
  chain's own range drifts. A range-free measure, the radius holding
  90% of the spin, gives 7 and 7 grains at 1 J0, 9 and 11 at 0.1 J0,
  and 5 and 7 at 3 J0.
- I keep the registered verdict (✘*) and show the profiles beside it.

### S: space may leave a lock's boundary link (X4)

| fabric | g/Δ = 0.1 | g/Δ = 0.2 | method |
|---|---:|---:|---|
| 3+2 | 12.1% | 24.9% | exact (`../fabric/` F2a) |
| 3+2 | 12.9% | 29.2% | the truncated method: the calibration |
| W | 5.6% (≈ 5.0% after the 3+2 calibration) | 9.7% (≈ 8.3%) | the truncated method |

These are holes made per resonant photon.

- All four W photon states give the same rate.
- The hole made carries the photon's charges (P5). Under I or W, the
  photon's state therefore decides p or n.
- **N8 passes in both fabrics.**

## 4. What the combinations add

- **T + S (3+2):** light makes holes, and the holes then walk.
  - At t = 0.1 J0 only N7 fails, and that is parity.
  - At t = J0, N2 also fails its registered fit.
- **I + T + S:** N6, N7 and N8 all pass, and only spin is missing.
  Everything a nucleon does except carry spin ½.
- **W + T + S:** every criterion passes. N2 for the moving W hole is
  the static value. N5 (no ladder) is reported and fails, as it does in
  every cell.
- **A hazard that only T + S shows:**
  - Under T, two holes that meet collapse into one empty cell, which
    cuts the fiber (§5).
  - Under S, light keeps making holes: broadband light about 0.2–1.4%
    per photon, resonant light 4–25%. Holes will meet.
  - **So a fiber in light, with T on, fills with holes, and pairs of
    them cut it.** Neither factor alone shows this.

## 5. Two holes (X5; reported, not gated)

Binding B = E(two holes) + E(no hole) − 2E(one hole), in J = J_h2.

| fabric | holes | d = 1 | d = 2 | d = 3 | d = 4 |
|---|---|---:|---:|---:|---:|
| 3+2, n = 8 | pinned (t = 0) | −3.44 | −0.05 | −0.27 | −0.02 |
| W, n = 4 | pinned (t = 0) | −5.24 | −0.13 | | |
| 3+2, n = 6 | free (t = J) | −20.4 | | | |
| W, n = 4 | free (t = J) | −16.5 | | | |

- **Adjacent holes bind by converting.** A (1,1) pair sits below its
  space quantum (Δ₁₁ = −0.024), so the binding uses real link
  occupancy, which V forbids.
- **Separated holes interact weakly.** The interaction oscillates with
  sublattice.
- **With transport, two holes merge into one empty cell next to a full
  one.** E_f(2) + E_f(0) − 2E_f(1) = −0.25, about −27 J. That is not a
  nucleus: it cuts the fiber.
- v98's nuclear ladder (deuteron barely bound, then the α-stairs) does
  not appear in either fabric.

## 6. What works and what does not

**Works:**

1. **The widening W.**
   - Its lock is stronger than 3+2's: the gap is 1.32 against 0.44 J0.
   - Its hole is more compact: 2.3 against 5.8 grains.
   - It gives spin ½ and the p/n pair at once. The "one open slot"
     object is a (½, ½) quartet, so p and n are the same object at the
     same size, as v98 asks.
   - This is the only reading tested that removes `../fabric/`'s
     parity wall without giving up spin.
2. **T at small t.** At t = 0.1 J0 the hole walks (band ≈ 0.11 J0 in
   3+2, ≈ 0.4 J0 in W) and keeps its shape (N2 passes).
3. **S.** About one resonant photon in 4–8 makes a hole (3+2), or one
   in 10–20 (W), and the photon's state picks p or n.
4. **W + T + S** passes every gated criterion.

**Does not work:**

1. **I alone.** The relabelling trades spin for p/n. It never has both.
2. **A ladder.** There is no bound excited level in any of the 12 cells,
   at any t and any length tried. Every excitation is the hole plus a
   free lock quantum, scattering repulsively. v98's discrete ladder
   (Δ-like states) needs an attractive channel, and nothing here
   supplies one.
3. **Nuclei.** Two holes do not form a bound pair of nucleons. They
   convert (t = 0, against V) or collapse into a fiber cut (t > 0).
4. **T at larger t.**
   - By t ≈ J0 the walking hole's cloud loses its single-exponential
     form, so the registered size test fails, though the cloud does not
     grow.
   - At 3 J0 it is a ferromagnetic polarization, the Nagaoka
     precursor; Nagaoka proper is near 9–10 J.
   - Walking and firmness pull against each other.
5. **T + S together.** Holes made by light meet and cut the fiber.

## 7. What would move past this page (each is the user's call)

1. **Rule W (the 3+4 widening).** It is v98's named fallback, and it is
   the reading this campaign found working. It is a law change, from two
   field modes per cell to four. A kernel step would follow from it, and
   needs authorization.
2. **Rule T with a bound.** For example: transport inside a lock is
   allowed, with t ≲ 0.1–0.3 J0 (≈ 0.001–0.003 in rate units), small
   enough to keep the hole's shape. Its value would be a new constant.
3. **Rule S as a boundary-only exception to V.**
4. **The two open physics problems** (not decisions):
   - a ladder needs an attractive hole–quantum channel;
   - nuclei need two holes that bind without converting or collapsing.

   Both point at the same missing ingredient: something in the S-law
   that holds two objects at a distance. Candidates to test next, on
   the same instruments:
   - the 4th-order terms at larger g;
   - W's own multiplet structure (an attractive (½, ½) ⊗ (1,1) channel
     in a different charge sector);
   - S-law terms beyond leading order around a moving hole.

## 8. Deviations from `PREREG.md`

1. **W chains.**
   - The first W chains (my DMRG at χ = 32, L_bulk = 12–28) were not
     converged. The spin and isospin profiles of the (½, ½) hole
     differed, though the model is exactly symmetric. On a 4-cell chain
     χ = 32 was 0.009 J0 off the exact energy, with discarded weight
     3e-4.
   - A truncation that never splits a degenerate block made no
     difference.
   - My DMRG at χ = 64 matched exact (1e-5) but was too slow. With the
     user's approval, a dependency was added: **TeNPy 1.1.1**
     (charge-conserving DMRG, in `~/.venvs/scp-tenpy`, outside the
     repo).
   - It was validated against:
     - the exact 4-cell W hole (1e-5);
     - my χ = 64 12-cell W hole (4e-6 lower);
     - my 3+2 chains (2e-10).
   - The W chains then ran **at the registered lengths (16, 24, 32),
     at χ = 128**.
   - The χ = 32 outputs are in `cas/out/chi32_superseded/` and are not
     used. My χ = 64 outputs (`x12_W_*_chi64.out`) agree with TeNPy to
     about 1e-3 and are kept as a cross-check.
2. **Sector penalties.**
   - My DMRG pinned sectors with a 1 J0 penalty. That is weaker than
     W's gap, so a run aimed at S^z = 1 alone stayed at S^z = 0.07;
     those runs were discarded.
   - TeNPy conserves the charges exactly, so no penalties remain. The
     gap is taken from the (1,1) member of the first multiplet, and
     (1,0) agrees to 1e-7.
3. **W with T** ran on exact rings only (n = 4, 6). A mobile W site has
   14 states, and was not attempted in DMRG. N2 under W + T is the
   static W hole's (marked ✔p).
4. **N2 under T.**
   - The registered quantity was the decay length of the hole-centred
     correlation C(d). On an open chain, C(d) carries the hole's
     particle-in-a-box weight, so a conditional (weight-normalized)
     profile was added.
   - Both fail at t = 1 and 3 J0, and both pass at 0.1 J0. The
     registered verdict is kept.
   - The profiles at the two lengths and the range-free 90% radius are
     reported beside it. They show a cloud that does not grow but is not
     a single exponential.
5. **X4 for W** uses the truncated-channel method. Calibrated on 3+2,
   it runs 7–17% high, and the W numbers are reported raw and
   corrected.
6. **Units.** X3's t is in J0 = 2g²/Δ (0.009586), not J_h2 (0.009231).
7. **The 3+2 mobile chains** ran at χ = 64 (TeNPy). My own runs of the
   same chains did not finish in 5.5 hours and were stopped. One of
   them finished (`x3c_M2_mob_40_3_t0.1`) and is kept as a cross-check.
8. **Background-job hygiene.** Two scheduler scripts survived a kill
   (bash had exec'd them under another command line). They relaunched
   χ = 32/48 runs for about an hour before being found and stopped.
   Nothing from those runs is used.

## Reproduce

From `v102/perm/cas/`. `PY` is `~/.venvs/scp-tenpy/bin/python` (TeNPy
1.1.1, numpy, scipy); everything else runs on the system `python3`.

| step | command | output |
|---|---|---|
| W lock, rings | `python3 -I w_lock.py` | `out/w_lock.out` |
| relabelling (I, I+W) | `python3 -I x6_relabel.py` | `out/x6_relabel.out` |
| T, rings | `python3 -I x3_rings.py {2,4} n` | `out/x3_rings_M*_n*.out` |
| S, formation | `python3 -I x4_form.py {2,4} {0.1,0.2}` | `out/x4_M*_g*.out` |
| two holes | `python3 -I x5_pairs.py` | `out/x5_pairs.out` |
| long chains | `./run_tp.sh` (TeNPy), then `$PY -I analyze_tp.py` | `out/tp_*.out`, `out/analyze_tp.out` |
| cross-check, own DMRG | `python3 -I lo_dmrg.py ...`, `python3 -I analyze_chains.py` | `out/x12_W_*_chi64.out`, `out/x3c_*.out`, `out/analyze_chains_ownDMRG.out` |

Instruments:

- `graphm.py`: the S-law with M modes per cell.
- `weff.py`: all-order effective clusters for M modes.
- `lo.py`: the leading-order model, with mobile holes and exact
  sectors.
- `lo_dmrg.py`, `dmrg_sym.py`: my DMRG, kept as the cross-check.
- `tp_chain.py`: the TeNPy chains.
