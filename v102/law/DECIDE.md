# v102 law step — decision

**Written:** 2026-10-07, after `cas/out/law.out` and
`cas/out/rereg.out`. Gates: `PREREG.md`.
**Status:** this page closes the law step. It does not change the law
table. It does not open Phase 2.

---

## Part L — does the law fix t, g, Δ? **Ending C**

### t: transport on the fiber (two readings)

| reading | source | t (rate units) |
|---|---|---:|
| T-a | the fiber inherits 3-D transport, `field_J·ŵ` with ŵ = ½ on a degree-2 ring | 0.90 |
| T-b | v98 §1, "locked field cannot travel as itself," read as no hop inside a lock | 0 |

T-b was added after the first computation, to complete the list. Its
source was already discussed in `../gates/WEDGE.md` W4. It was not
chosen for its result.

### Δ: the pair-conversion detuning (four readings)

Field quanta load their cell's pitch (`w/(1 + q_detune·x)`). At the
R6 ceiling, a full cell holds 2 quanta at pitch 1.336, and converting
one quantum from each of two full cells frees 2.421 (rate units).

| reading | space quantum's pitch | Δ | lock? |
|---|---|---:|---|
| R-a | bare dense rung w2; field also bare | −0.400 | no: the converted state is the ground state |
| R-b | bare dense rung w2; field loaded | +0.479 | needs t/Δ < t_c |
| R-b′ | w2 loaded by the space quantum itself | −0.025 | no |
| **R-c** | **no pitch: space is a continuum** | 0, with decay | **no stationary lock** |

**R-c is the reading that matches the current representation.** v99
§1 retired the dense fluid. Under `pattern_only`, conversion is field
↔ space. The kernel's pitch load never includes `Es`
(`freecell.c:1575`), so space has no pitch. A field pair always finds
a resonant space state, so the would-be lock decays at the
golden-rule rate. R-a, R-b, and R-b′ all borrow the retired dense
rung for the space quantum.

### g: **undefined**

The table converts energy as a rate law: `f_conv` of the excess per
step, `k_dep`, atoms at the source's pitch, and a credit that lapses
at two atoms (`freecell.c:2740`). There is no coherent pair-to-space
amplitude to read.

### Combined

| | R-a | R-b | R-b′ | R-c |
|---|---|---|---|---|
| T-a (t = 0.9) | no | **outside**: t/Δ = 1.88; t would need to be 7–12× smaller | no | no |
| T-b (t = 0) | no | inside, *if* a small g is supplied | no | no |

Only one combination, T-b × R-b, hosts the lock. It uses the strong
reading of v98 §1, and the dense rung that v99 retired. It still
needs a `g` the table does not have. Per `PREREG.md` it is reported,
not chosen.

**Ending C:** g is undefined, and the reading that matches the current
representation (R-c) has no stationary lock. Ending B would need every
reading outside. T-b × R-b keeps that from being claimed.

### The sentences that would fix it

Each would be a law statement, so a law-table change under the
ratchet, and the user's ruling:

- **S1.** A link's space mode is a discrete oscillator with its own
  pitch `Ω`. Without it, R-c holds and any lock decays.
- **S2.** `Ω` exceeds the pair it converts, so Δ = Ω − 2.421 > 0. Under
  the bare dense rung, Ω = 2.9 and Δ = 0.479.
- **S3.** Conversion on a link is a coherent pair amplitude `g`, with
  `g/Δ` small.
- **S4.** Transport inside a locked fiber is slow: t/Δ below about
  0.16–0.26. Either t = 0 (v98 §1 read strongly), or a fiber transport
  constant ≤ 0.08–0.12 in rate units, which is 7–12× below the 3-D
  optics value.

With S1–S4 (Ω = w2, t = 0), the fiber is the Phase 1 lock: Δ = 0.479,
and ξ ≈ 6 grains minus an O((g/Δ)²) drift that `g` would set
(`../phase1/DECIDE.md` P3/P5).

## Part R — the two defective bars, re-registered: **PASS**

| bar | case | result | bar |
|---|---|---|---|
| R-P1′ | n = 6, t/Δ = 0.05 | max miss 1.2% (12 levels with multiplicity), singlet | 5% |
| R-P1′ | n = 6, t/Δ = 0.15 | max miss 4.9%, singlet | 5% |
| R-P6c′ | slits ±16, half-width 2, new source | eraser = coherent single-slit sum, 5e-16 | 1e-10 |

Both new cases were run only after this page's gates were written. The
Phase 1 misses at P1 (n = 4) and P6(c) were defects in how those bars
were written. They do not change Phase 1's stop: P3 (ξ's drift) and
the transport threshold remain physics.

## What is now true, and what is not

**True:**

1. Under the current law, read in its current representation, the
   derived lock does not exist. A field pair is resonant with a
   pitchless space and decays.
2. Borrowing the retired dense rung gives Δ > 0 only when the field is
   loaded (R-b). Even then, 3-D transport is about 10× too fast.
3. The table has no coherent conversion coupling.
4. The lock needs four new law statements (S1–S4), not a reading.

**Not true:**

- The table predicts ξ.
- R-b is the law's reading. It borrows a retired mode.
- Any constant was tuned.
