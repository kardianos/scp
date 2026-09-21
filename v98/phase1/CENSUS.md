# Phase 1A — graph census

**Status:** executed 2026-08-13. Paper only.
Rule: `CLOSURE.md`. Enumerator: `census.py` (NetworkX atlas, all
simple graphs on `n ≤ 7`). Spectra: `cas/spectra.input` (FriCAS)
and the factorisation recorded below. Raw TSV:
`runs/census_atlas.tsv`.

---

## 0. Verdict in four lines

1. Under C1–C3 there is a **unique smallest** lock candidate:
   `K₃ = C₃`.
2. After that the set is a **zoo** (1, 3, 11, 60, 502 types at
   `n = 3..7`). Closure alone does not mint a species.
3. The measured filter F-odd (odd rings π-frustrated, v90) moves
   the unique smallest to **`C₄`**. Still a zoo at larger `n`.
4. The strictest named pair, F-odd + F-3, has a unique smallest
   **`K_{3,3}`** (`n=6`, `e=9`) and only one type at `n=7`
   (`K_{3,4}`). That is the nearest thing this census has to a
   line. It is not a proton.

Question 5 of `CLOSURE.md` is answered: **raw graph closure does
not force a monodisperse species.** 1B is load-bearing.

---

## 1. Counts

Atlas graphs scanned: 1252 (`n ≤ 7`). Cycle-closed: **577**.

| filter | n=3 | 4 | 5 | 6 | 7 | unique smallest |
|---|---:|---:|---:|---:|---:|---|
| C1–C3 | 1 | 3 | 11 | 60 | 502 | `K₃` (`e=3`, `κ=2`, girth 3) |
| F-cut | 1 | 3 | 10 | 56 | 468 | `K₃` |
| F-3 | — | 1 | 3 | 17 | 136 | `K₄` (`e=6`, `κ=3`, girth 3) |
| F-odd | — | 1 | 1 | 5 | 9 | `C₄` (`e=4`, `κ=2`, girth 4) |
| F-odd + F-cut | — | 1 | 1 | 5 | 8 | `C₄` |
| F-odd + F-3 | — | — | — | 1 | 1 | **`K_{3,3}`** (`e=9`, `κ=3`, girth 4) |

F-odd closed types identified:

| n | e | graph | κ | notes |
|---|---:|---|---:|---|
| 4 | 4 | `C₄` | 2 | unique smallest even lock |
| 5 | 6 | `K_{2,3}` | 2 | |
| 6 | 6 | `C₆` | 2 | v93 vacuum ring |
| 6 | 7 | `C₆+chord` | 2 | |
| 6 | 8 | `K_{2,4}` | 2 | |
| 6 | 8 | prism-like cubic-ish | 2 | |
| 6 | 9 | **`K_{3,3}`** | **3** | unique F-odd+F-3 at n=6 |
| 7 | 8–12 | nine types | 1–3 | only `K_{3,4}` has `κ=3` |

The atlas stops at 7. At `n=8` the cube `Q₃` and `K_{4,4}` both
enter F-odd+F-3, so the near-line does not survive the next
integer without a further rule.

---

## 2. Exact spectra (FriCAS)

Characteristic polynomials of the combinatorial Laplacian
`L = D − A`, factored over `ℤ`. Traces equal `2|E|`.

| graph | `det(tI − L)` | eigenvalues |
|---|---|---|
| `K₃` | `t(t−3)²` | `0, 3, 3` |
| `C₄` | `t(t−2)²(t−4)` | `0, 2, 2, 4` |
| `K₄−e` | `t(t−2)(t−4)²` | `0, 2, 4, 4` |
| `K₄` | `t(t−4)³` | `0, 4, 4, 4` |
| `C₆` | `t(t−1)²(t−3)²(t−4)` | `0, 1, 1, 3, 3, 4` |
| `K_{3,3}` | `t(t−3)⁴(t−6)` | `0, 3⁴, 6` |

`C₆` recovers the CHARTEXP / qgraph book `{0,1,1,3,3,4}` exactly.
That is a check, not an adoption of chart write-back.

These shells are **diagnostics**. They are not painted onto hoards.
qgraph already said: magic numbers are occupancy of a *given*
graph. This census asked the other question — *which graphs are
even candidates* — and the answer is “many.”

---

## 3. What the smallest objects are, and are not

- **`K₃`.** Unique smallest simple cycle-closed graph. Odd.
  v90 says odd rings are π-frustrated. So `K₃` is the
  combinatorial minimum and, on this substrate, a bad lock.
- **`C₄`.** Unique smallest even lock. Edge-minimal. `b₁ = 1`.
  Loses cycle-closure if any edge is deleted; loses F-3 if any
  two vertices are deleted (`κ=2`). A thin lock.
- **`K₄`.** Unique smallest 3-connected. Odd triangles. Same
  F-odd objection as `K₃`.
- **`C₆`.** The programme’s vacuum ring. Even, not 3-connected.
  Holds winding in vacuum; is a loop, not a species (v97 H0/H1).
- **`K_{3,3}`.** Unique smallest even *and* 3-connected.
  Non-planar. That is compatible with “the lock’s geometry is
  not our 3-space” (the fiber / 1B). It is still just a graph.

None of these is a nucleon. Unique-*smallest* is not
unique-*size*. The zoo at the next `n` is the hoard, combinatorially.

---

## 4. What this does to Phase 1C

1A alone cannot authorize a substrate. It can constrain one:

- Any lock story that produces a continuum of sizes has ignored
  C1–C3 (already the Face C miss).
- Any lock story that needs an odd triangle as the ground object
  has to defeat the v90 F-odd measurement.
- If 1C wants one smallest robust even lock, the census names
  `K_{3,3}` and then demands a rule that kills `K_{3,4}`, the
  cube, and `K_{4,4}`. That rule is not C1–C3. It has to come
  from 1B (a trough in `n`) or from occupancy (`Δ/T` on one
  chart — qgraph M3, which selects *N* on a given *G*, not *G*).

Reproduce:

```
python3 v98/phase1/census.py
fricas -nosman -noht < v98/phase1/cas/spectra.input
```
