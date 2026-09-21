# chartexp results — inner-chart on/off (2026-08-13)

Fork only (`v97/chartexp/fab`, binary `chartexp/chartrun`). Live
`v97/fab` untouched. Go tests: `go test ./chartexp/fab/` green
(C6 Jacobi recovers `{0,1,1,3,3,4}`).

## What ran

All arms: `exp=ring ring_n=6 seed_mw=2 amp_nat=2 amp_logate=1
hop_order=1 seed=20260802 k_rad=0.05`. Vacuum unless noted.

| wave | T | field | arms |
|---|---|---|---|
| tiny | 40 | dark | off, proj, grain, mix × {energy, mass}, area, motif-proj |
| medium | 200 | dark | off, proj, grain, mix, mix-mass, mix-motif |
| lit | 80 | `noise_amp=0.4` | off, proj, grain, mix, mix-motif, mix-mass |

`chart=1` (proj) is diagnostic. `chart=2` (grain) overwrites `w2e`
at the start of pass 6 (door + beat). `chart=3` (mix) overwrites
`w2e` after pass 1 (hops + door + beat).

## 1. Numerical differences

### Diagnostic is inert

`t_off` ≡ `t_proj` and `l_off` ≡ `l_proj` after stripping `# CHART`
/ `# chartexp` lines. Drift, cond, Em, geometry: identical.
**Storing the 3+1 book does not move the simulation.**

### Vacuum (no door traffic)

Every arm: cond = rad = 0, Em_tot = 4.875, winding held (W = 2
where printed). Differences:

| arm | drift | weff | n0+n1 (mean) | notes |
|---|---|---|---|---|
| off / proj / grain | 1.735e-15 (T=40) | 1.365 (star) | 2.44 | grain: same Em; qcnvF 0.054→0 and gg 0.007→0.0001 |
| mix energy | 1.72e-14 | 1.366 | 2.44 | hops see weff; drift ×10, still FP-floor |
| mix mass | 2.57e-14 | 1.232 | 2.44 | |
| motif proj | 1.735e-15 | **1.732** | 4.87 | √3 = C6 λ=3 shell |
| area proj | 1.735e-15 | **1.4e-6** | 2.44 | eigenvalues scale with slot area; pitch collapses |

Medium T=200: same pattern, mix drift 2.4–3.4e-13 vs off 9.6e-14.
Winding still 2. No conversion on any vacuum arm (atom grain, as
p2_ring).

**Motif + energy** is the only book that matches the theory object:
painted `m=2` occupies the degenerate √3 pair (`n0≈n1≈2.44`,
residual ~0). A **star** (cell + 2 ring neighbours) is a path-3,
weff≈1.36 — a different chart. That is not a bug; it is what
“inner cell chart” *means*.

**Area** is not a usable accounting here unless weights are
normalized. Raw `sA` on this ring is tiny, so `ω → 0`.

### Lit ring (doors can fire)

| arm | cond | Em_tot | z_live | W | weff |
|---|---|---|---|---|---|
| off / proj | **0.6518** | 5.625 | 0.67 | **0** | 1.21 (3 cells left) |
| mix-motif | 0.6518 | 5.625 | 0.67 | 0 | 1.37 (remnant) |
| grain / mix / mix-mass | **0** | 4.875 | 2.00 | **2** | 1.23–1.38 (C6 held) |

Replacing `w2e` with the star pitch **shuts condensation** and
the ring keeps its winding. Off/proj convert noise→matter
(cond=0.65), the cycle breaks (2 live cells), W dies.

Cause (not a mystery): pass 6 does `cbeta += (w1e − w2e) dt`.
Table `w2e ≈ 2.09` at x=0.325; occupied star `weff ≈ 1.36`;
`w1e ≈ 1.19`. Beat gap 0.90 → 0.17. Fewer beat fires → cond
never reaches an atom (qcnvD leftover 0.42 vs 1.18). Grain was
*intended* as ε-only; the hook overwrites `w2e` for the whole
door pass, so the beat changes too. Mix does that on purpose.

Motif mix (`weff = √3 ≈ 1.73` while C6 lives) keeps a larger
gap to `w1e` and condenses like off — then melts the same way.

Conservation: |drift| ≤ 6e-14 on every lit arm. No ledger leak
of the Em+Ee+Es kind.

## 2. Are the differences valid?

**Valid as arithmetic of this hook**

- proj vs off: null. The book is a spectator.
- vacuum mix vs off: real, small, in the unitary hop (pass U
  uses `w2e` for τ). Energy conserved; winding not touched.
- lit grain/mix vs off: real and large (cond on/off, W 2 vs 0).
  Mechanism is the **beat**, not a selection rule on `Δm`.
- motif vs star: real difference in the *chart you chose*.
  Motif recovers the C6 shell the seed painted. Star does not
  claim to.

**Not yet valid as “more correct composition physics”**

1. Occupations `n_a` never enter the door as a promotion `a→b`.
   Only a scalar `weff = Σ n_a ω_a / Σ n_a` is written back.
   That is a one-number fit, not 3+1 accounting in the write
   path.
2. Grain leaked into the beat. A clean ε-only test was not run.
3. Shutting the door by collapsing `w2 − w1` protects winding
   the way turning `f_conv=0` would. That is not evidence the
   *composition* was respected; it is evidence conversion was
   skipped.
4. Area accounting is a unit error, not a theory variant.
5. Single seed, NC=6, no bath foam. Nothing here is H1.

## 3. Batch 2 — clean grain, mode door, composition H0

`chart=2` is now **clean grain** (`ε` from weff in `atomsFire`
only; `w2e` and the beat stay on the table). `chart=4` is the
mode door (`ε = A0|ω0−w1e|/2π`, residual priced as
`(n0+n1)/tot`, refuse if nothing occupied). Script:
`chartexp/run_batch2.sh`.

### (1) Clean grain — reconverged. Beat leak confirmed.

| pair | stripped logs |
|---|---|
| vacuum off vs grain (star or motif) | **IDENTICAL** |
| lit off vs grain (star or motif) | **IDENTICAL** |

Lit both: cond = 0.651795, z_live → 0.67, W → 0. The batch-1
“grain/mix hold W” was entirely the pass-6 overwrite of `w2e`
collapsing the beat. Door `ε` alone does not move this apparatus.

### (2) Mode door — not useful (yet)

| arm | cond | z_live | W | chart on |
|---|---|---|---|---|
| lit off | 0.652 | 0.67 | 0 | (n/a) |
| lit mode motif | **0.413** | **0** | 1 then stuck | 0 (isolated) |
| lit mode star | **0.251** | **0** | 1→0 | 0 |
| vacuum mode motif | 0 | 2.00 | 2 | 6, √3 pair held |

Usefulness bar was: W held *with* cond > 0, or cond that stays
in the √3 pair. **Miss.** Cond is smaller (gap ε ≈ 0.10 vs table
~0.38, plus residual price) but the ring still dies, and dies
*harder* (z=0, no contacts). Vacuum mode is almost off (one
credit line differs: radiance leftover sees the smaller ε).

A promotion `a→b` that only resizes the atom and taxes residual
is not a selection rule. It does not keep occupation on the
painted shell once the vertex door fires.

### (3) Composition H0 — PASS as a meter

**Vacuum** (`b2_h0vac`, motif+energy, T=200, every 8 t.u.):

- `o = (1.732, 1.732)` and W = 2 on every frame
- `n0/n1` starts 2.4375/2.4375, slow skew to 2.556/2.314
- residual 0 → 0.0027 (< 0.06% of Em)

The book names the painted `m=2` shell and keeps it. The slow
n0↔n1 drift is partner mixing inside the same degeneracy, not
leak to λ=0 or λ=4.

**Lit** (`b2_h0lit`, same book, table dynamics):

| t | W | weff | n0, n1, n2, res | z |
|---|---|---|---|---|
| 0–20 | 2 | 1.732 | 2.44 / 2.43 / 0 / 0 | 2.00 |
| 24 | 2 | 1.585 | 2.25 / 0.76 / 0.74 / **0.31** | 1.33 |
| 32 | **1** | 1.463 | scramble, o=(1.41, 1.85, 0) | 1.00 |
| 36+ | **0** | 1.38 | remnant 3-cell, o=(1.73, 1.00) | 0.67 |

The √3 pair empties in the same window W dies. The book is
tracking winding death, not a frozen label. That is what H0
asked.

### (4) Star vs motif — still not a fit knob

Unchanged. Motif recovers C6. Star is path-3. Batch 2’s mode
door is worse on both; that does not license averaging them.

### (5) Fit — not valid, not run

A fit was allowed only if 3+1 + a mode door collapsed to a
cheaper two-pitch law *after* doing a real job. Clean grain is
a null. Mode door is a worse melt. The only working object is
the **diagnostic book**. Writing `weff` back into `w2e` is the
move batch 1 already showed fakes a hold via the beat. Fitting
that map would re-enact the leak.

No `w2e(n_a)` or `ε(ω0)` law is adopted.

## 4. Bottom line (after batch 2)

- **Keep:** motif+energy 3+1 projection as an instrument. It is
  free (proj ≡ off), it holds the √3 pair in vacuum, and it
  tracks W-death in a lit bath.
- **Drop:** feeding a scalar weff into `w2e` (beat leak);
  area accounting; treating star and motif as interchangeable.
- **Not yet a face:** a mode door that only changes `ε` and
  taxes residual. Selection has to act on *which* shell
  receives the atom, not on how big the atom is.
- **Do not fit** the table pitches to this book until a door
  exists that holds composition *while converting*.
