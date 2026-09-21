# Inner-chart probe — note for incorporation (2026-08-13)

**Fork only.** Code and logs: `v97/chartexp/`. Live `fab/` was not
changed. Go kernel only. Not a principle-face. Not a law-table
change. Do not adopt `chart=*` into the carried kernel.

---

## Recommendation

**Incorporate the motif 3+1 projection as a hunt/ANLZ instrument
only.** Do not feed chart pitches back into `w2e`, the beat, hops,
or door `ε`. Do not open the chart-pitch face on this evidence.
Do not fit `{w1,w2}` to the book.

If you add anything to the live tree: compute the Laplacian of the
**tagged / small motif** (n ≤ 24), keep the three most occupied
shells plus residual, project `ψ_m` onto them, print
`{ω_a, n_a, res, W}`. Default off. Physics path untouched. That is
the same contract as `promo_i`/`promo_o` when `promo=0`.

Do **not** copy `chart=2/3/4` write-path hooks from the fork.

---

## Why

P2 says pitch is `ω = √λ` of the current contact graph. The live
kernel carries two table numbers, load-stretched:

`w1e = w1 / (1+q x)`, `w2e = w2 / (1+q x)`.

A cell cannot be a mixture; a painted `seed_mw=2` is a phase, not
an occupation of `χ_2`. The question was whether a cheap 3+1 book
on a small chart (a) tracks that composition and (b) enables more
correct door/hop actions.

**(a) yes, as a meter. (b) no, not with anything tried.**

The C6 object the hunt already uses is exactly four shells
(`m = 0, ±1, ±2, 3`). Motif + unweighted L puts a `seed_mw=2`
ring on the degenerate `ω = √3` pair (`n0 ≈ n1`, residual ~0) and
keeps it for T=200 in vacuum. In a lit bath the same book
scrambles in the same window winding dies. That is worth reading.
It is not worth steering.

---

## What was tried

Fork: `v97/chartexp/fab` (copy of `fab/`) + `chartrun`.
Tests: Jacobi on C6 recovers `{0,1,1,3,3,4}`.
Logs: `chartexp/runs/`. Narrative: `chartexp/RESULTS.md`.

Apparatus (all waves): `exp=ring ring_n=6 seed_mw=2 amp_nat=2
amp_logate=1 hop_order=1 seed=20260802 k_rad=0.05`. Vacuum unless
`noise_amp=0.4`.

| knob | meaning |
|---|---|
| `chart=0` | off |
| `chart=1` | project 3+1; no write-path change |
| `chart=2` | door `ε` from occupancy-weighted `√λ` (batch 2: `w2e` untouched) |
| `chart=3` | overwrite `w2e` (hops + beat + door) |
| `chart=4` | mode door: `ε = A0\|ω0−w1e\|/2π`, residual priced; beat untouched |
| `chart_scope` | `star` (cell+neighbours) or `motif` (connected live-Em, n>24 → star) |
| `chart_acct` | `energy` (unweighted L), `mass` (√Em), `area` (slot area) |

Batch 1: tiny T=40 / medium T=200 vacuum; lit T=80. On vs off,
three accountings, star vs motif, grain/mix writeback.

Batch 2: clean grain (ε only); mode door; composition H0 time
series (vacuum T=200 every 8 t.u.; lit T=80 every 4 t.u.).

---

## What failed

**1. Writing a scalar `weff` into `w2e` (batch 1 grain/mix).**
Lit ring: cond 0.65 → 0, W held at 2. Looked like composition
protecting winding. It was the beat. Pass 6 does
`cbeta += (w1e − w2e) dt`. Table `w2e ≈ 2.09` at x=0.325; star
`weff ≈ 1.36`; `w1e ≈ 1.19`. Gap collapses, condensation never
reaches an atom. Same fake hold as `f_conv=0`.

**2. Clean grain (batch 2, `chart=2`).** `ε` from weff, table
`w2e` left on the beat. Vacuum and lit logs **byte-identical**
to off (CHART lines stripped). Door grain alone does not move
this apparatus. Confirms (1) was the leak.

**3. Mode door (batch 2, `chart=4`).** Usefulness bar: W held
*with* cond > 0, or cond that stays in the √3 pair. Miss. Lit
motif: cond 0.41 (not 0, not 0.65), z_live → 0 (worse than
off’s 0.67), W → 1 and stuck. Star mode: cond 0.25, also
disconnects. Resizing the atom and taxing residual is not a
selection rule. The vertex door still dumps occupation.

**4. Area accounting.** `weff ~ 1e-6` on this ring. Eigenvalues
scale with raw `sA`. Unusable unless weights are normalized.
Not a theory variant.

**5. Fit `{w1,w2}` to the book.** Pre-registered as allowed only
if a mode door did a real job and then collapsed to a cheaper
two-pitch law. Clean grain is a null. Mode door is a worse
melt. Fitting `weff → w2e` would re-enact (1). **Not run.**

**6. Treating star and motif as interchangeable.** Motif recovers
C6 (`weff = √3`). Star is a path-3 (`weff ≈ 1.36`). Different
charts. Do not average; do not “tune” one toward the other.

Conservation was never the issue: |drift| ≤ 6e-14 on every arm.

---

## What to incorporate (if anything)

Into `hunt/` or an ANLZ column, on tagged/small motifs only:

- unweighted Laplacian, 3 most-occupied shells + residual
- per-snap `{ω0,ω1,ω2,n0,n1,n2,res}` and ring W
- H0-shaped checks already passed in the fork:
  - vacuum C6 + `seed_mw=2`: `ω0=ω1=√3`, `n0≈n1`, res ≈ 0, W = 2
    for the run (slow n0↔n1 skew inside the pair is allowed)
  - lit melt: that pair scrambles in the same window W drops

Do not put this in the step function. Do not let it touch `w2e`.
Do not open chart-pitch until a door exists that holds
composition *while converting* (promotion `a→b` that chooses
the receiving shell, not the atom size).

The other agent should treat this page as the verdict. Details
and tables: `chartexp/RESULTS.md`.
