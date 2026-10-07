# G-two — decision

**Written:** 2026-10-07, after `WEDGE.md` and `cas/out/two.out`.
Gate: `../PREREG.md` §3.
**Status:** this page closes G-two. It does not edit G-wedge.

---

## 1. What the premises say about `N`

- With transport refused at the ceiling, a conversion returns its
  quanta to the cells they came from. For each cell, `N_i` plus the
  space quanta on its two links is conserved. `N` is a label of the
  state.
- Quanta hoard: the energy per quantum falls with `N` at both ring
  sizes.

| N | S | e(N), n = 6 | e(N), n = 8 | gapped? |
|---:|---:|---:|---:|---|
| 1 | ½ | −1.434 | −1.413 | no (LSM) |
| 2 | 1 | −2.436 | −2.417 | yes |
| 3 | 3/2 | −3.433 | −3.411 | no (LSM) |
| 4 | 2 | −4.430 | −4.408 | yes |

  (units g²/Δ; `two.out`)

So cells fill to the ceiling, and **`N = N_c`**. The premises fix that
the ceiling sets `N`. They do not fix the ceiling's value.

## 2. What the record offers for the value

| source | what it says | can it set the value? |
|---|---|---|
| two-atom credit | "a credit that lapses at two atoms" (v89 law 1): a cap on accumulated sub-atom credit | no. It is not an occupation (PREREG §5) |
| `cap = 2.5` (laws_V3a) | an energy | no. `cap/ε(ω)` depends on pitch and is generally not an integer |
| v98 §1 | a ceiling exists | no value |
| R5 (a gapped lock) | odd `N_c` gives half-integer `S`, which is gapless (Lieb–Schultz–Mattis; v99/reopen E2 at S = ½) | parity only: 2, 4, 6, … all remain, with ξ ≈ 6.03, ≈ 50, … |
| "take the smallest even" | — | a minimality preference the record does not hold as a principle (PREREG §3) |

## 3. Decision

**G-two fails.** `N` equals a ceiling whose value the premises do not
fix.

## 4. The sentence that would fix it

> **The occupancy ceiling of a lock is one quantum per link end.**

This is motivated by a reading, not derived. Each conversion takes
exactly one quantum from each end of its link (W1). If each link
carries one actuation (v89 §3 says "a link … actuates"), a link end
needs one quantum and has room for one.

On the fiber ring (degree 2) this gives `N_c = 2` and `S = 1`. On a
cell of degree `d` it gives `S = d/2`, which is the AKLT structure.

If this sentence is ruled, both open items close. W4's ceiling exists,
and `N = 2`. The fiber is then the spin-1 antiferromagnetic Heisenberg
ring with `J = 2g²/Δ`:

- gap ≈ 0.41 J;
- ξ ≈ 6.03 grains at leading order;
- one triplet line;
- corrections of order `(g/Δ)²`.

## 5. What is now true, and what is not

**True:**

1. The conversion bond hoards quanta, so a cell's occupation equals
   the lock's occupancy ceiling.
2. Odd ceilings give no gapped lock.
3. Nothing in the record sets the ceiling's value.

**Not true:**

- `N = 2` is derived.
- The two-atom credit is an occupation.
- The smallest even ceiling is preferred by a principle.
