# FREECELL — state of the free-cell substrate programme

**Date:** 2026-08-02. **Status:** exploratory; no kernel change made.
**Governing docs:** `PRINCIPLE.md`, `ROADMAP.md` §8, `DESIGN_2026-08-02.md`
(the running record this distils), `battery/README.md` (the ratchet rule
and the law-dependency map), `LIVEFAB.md` (the earlier dynamic-substrate
design).

This document is a handoff. It states what a "free cell" is, what has
been measured, what has been refuted, what was *claimed* and then
withdrawn, and what stands between here and running the unification
battery on free cells instead of a frozen foam. Every claim carries the
file and log that backs it.

---

## 1. What a free cell is, and why

`cellfab.c` today holds a **frozen foam**: `NC` and `NL` are fixed after
init, every array is `malloc(NC * ...)`, no cell or channel is ever born
or dies, and one global `dt` advances everything synchronously. That is a
permanent index set in space *and* in time — the construction
`PRINCIPLE.md` constraint 2 forbids. GLM's review named it as the
programme's weakest point: the no-background claim has only ever been
tested on an instrument that still has one.

A **free cell** is the proposed replacement: a cell whose extent (and
possibly shape) is a dynamical consequence of its own energy, whose
contacts are re-derived from current geometry rather than fixed at init,
and which advances on its own local clock rather than a global tick.

The programme splits into two independent halves, and they have very
different maturity:

| half | question | status |
|---|---|---|
| **execution** | can the substrate be advanced without a global clock, conservatively and deterministically? | **solved** |
| **substance** | can a free-cell fabric hold a localised object together? | **unsolved** |

---

## 2. WHAT WORKS

### 2.1 Local clocks — READY (`localclock.c`, `charge/runs/localclock.log`)

Four-criterion readiness suite on a degree-7 graph, N=96, energy moving
on channels as paired antisymmetric transfers.

| criterion | result |
|---|---|
| **R1 conservation** | drift +1.4e-14 on 120.68 = **1.2e-16 relative**, identical to synchronous, and *independent of lookahead* |
| **R2 convergence** | against a common dt/64 reference: sync falls 2.03/2.07/2.14/2.33×, async **2.05/1.94/2.03/2.07×** — first order, same solution |
| **R3 determinism** | **0.000e+00** under reversed scan, rotated origin, and both; arrival order differs by 1.2e-1 |
| **R4 width** | ~all events eligible at once for lookahead ≥ 0.05 |

**Four conditions, each a measured requirement, not a preference:**

1. **Order events by a total function of state** — `(t, kind, index)`.
   Never by arrival. The tie-break must be *total*: breaking ties by scan
   position instead of index leaked back in at 3e-6, which is small and
   not byte-identical, and the ratchet compares bytes.
2. **Bound skew in LOCAL TIME**, never tick count. Cells take different
   steps, so equal tick counts are not equal time; a tick bound distorts
   the ordering it exists to protect (1.3e-2 vs 1.1e-4).
3. **Channels own their transfers and their own clock.** Firing a
   transfer once per *endpoint advance* double-counts — each edge moves
   `gE(h_i+h_j)` per round instead of `gE` per unit time, an O(1) rate
   error that does **not** vanish as dt→0. Channel ownership was worth
   **160×** on transport accuracy (1.8e-2 → 1.1e-4).
4. **Never order by the tick counter.** It diverges without bound by
   construction — that divergence *is* the dilation. Instantaneous
   neighbour tick skew reached 1266 while local times stayed inside the
   lookahead. A mod-M byte is fine as a cyclic index; it cannot carry
   causality. If used for anything cyclic, assert `K < M/2` — the wrap
   failure is **silent** (the comparison stops ordering and the bound
   evaporates with no symptom).

### 2.2 Batch-parallel execution — bit-identical to serial

Events sharing no cell touch disjoint memory and commute exactly. An
eligible event joins the batch iff it is the minimum `(t, index)` over
its conflict neighbourhood — a local maximal independent set, a pure
function of state.

* **Batch is bit-identical to serial** (max|dθ| exactly 0). It is not an
  approximation of the schedule; it *is* the schedule, executed wide.
* **Bit-identical at 2, 4, 8, 16 threads.**
* Conflict rule (this was a bug, invisible below 8 threads): a cell event
  **reads its neighbours' published phase**, so two *adjacent cell
  events* race even though they write disjoint memory. The neighbourhood
  is incident channels **plus adjacent cells**.

### 2.3 Geometric pressure bounds size (`geopress.c` §G1)

`r = E^(1/3)`, contact repulsion `(k/2)δ²`, energy conserved and flowing
down `μ = (1/3)E^(-2/3)P`. **No maximum size anywhere in the model.**

| box L | φ | z | r_sd/r_mean | r_max/r_mean | dE/E |
|---|---|---|---|---|---|
| 22 | 0.104 | 3.31 | 0.0448 | 1.0759 | +2.1e-16 |
| 18 | 0.191 | 3.76 | 0.0683 | 1.1577 | −4.2e-16 |
| 14 | 0.405 | 4.47 | 0.0957 | 1.1997 | −6.3e-16 |
| **12** | **0.643** | **7.33** | 0.1318 | **1.2424** | +0.0e+00 |

Size stays bounded at every coordination — no cell exceeds 1.2424× the
mean — with conservation at the FP floor. (The `dE/E` column is quoted
from the current `geopress.log`; earlier revisions of the file printed
slightly different FP-floor values because the contact term now routes
through `kij`/`pow`. The physics and the `r` columns are unchanged.) Only the L=12 row is at jamming
(φ_J ≈ 0.64); the others are dilute, which matters enormously in §3.

**Why this works where consonance failed:** consonance failed by
*counting* (d simultaneous commensurability conditions at degree d,
generically unsatisfiable). A pressure balance is **one scalar equation
per cell** — nothing to frustrate.

### 2.4 Shape quantisation inside a cell (`shapelock.c`)

A cell with semi-axes (a,b,c) has three internal intervals obeying an
**identity**, `(a/b)(b/c)(c/a) = 1` — the comb-closure condition of
`charge_sym.mac`, with only two of three intervals free. One cell, one
closed triangle, **no neighbours to disagree with**.

**S1 — the shape spectrum** (exact integer enumeration; no solver, so no
search artifact):

| Tenney ceiling | admissible shapes | complete? |
|---|---|---|
| 0 | 1 (`1:1:1`, the sphere) | yes |
| 1 | 3 | yes |
| 2 | 8 | yes |
| 3 | **26** | yes |
| 4 | 64 (n≤12) vs 95 (n≤32) | **no, truncated** |
| 6 | ≥538 | no |

Through ceiling 3 the count is identical at both enumeration windows —
**finite and complete**. At v89's `comb_limit=6` the set is ≥538 and
unconverged, large enough that "selection" is doubtful.

**S2 — a genuine staircase** on the integer overtone series:
`1/1` (t 0–0.75) → `2/1` (0.90–1.20) → `3/1` (1.35–1.65) → `4/1` (1.80)
→ `5/1` (1.95–2.10) → 6/1, 7/1, 8/1, 10/1, 11/1, 13/1. Flat plateaux,
converged at max|grad| ~1e-8, plateau width shrinking with ratio exactly
as log-space rung density predicts.

**S3 — control is clean**: comb collapsed to unison gives 1.00000 →
1.00120, linear, no steps.

**S4 — the comb span IS the aspect-ratio bound**: a:b saturates at
precisely the outermost rung at every ceiling (2, 4, 8, 16, 32, 64, 128
for ceilings 1–7), residual shrinking 3.1e-3 → 6.6e-4. A bound on shape
made entirely of intervals, **no length anywhere**.

**The quantised object this programme was looking for is a SHAPE, not a
size.** Every attempt to quantise size failed or reduced to an
idealisation; shape quantises cleanly.

---

## 3. WHAT DOES NOT WORK

Distinguish two categories carefully. Some things are **refuted**; others
were *claimed refuted* and are actually **untested**, because the
configuration could not have shown the effect.

### 3.1 REFUTED — consonance between cells (`harmfrust.c`)

Random graphs, degree 2–12, every cell given its own preferred size, `kE`
scaled with degree (at fixed `kE` every degree collapses to the unison
crystal and nothing is measured):

* locked fraction **0.05–0.13 at every degree including the ring**, no
  trend; no trend with comb span either (0.02–0.10 across ceilings 1–8);
* surviving locks use high-Tenney-height rungs (mean_H ~3–4), not the
  clean series;
* `max_domain` 2–6 with no structure — **no consonant-domain scale**.

Degree 2 with generic targets gives 0.062, which bounds `harmgrow`'s ring
staircase retroactively: that was **one impurity in a uniform medium**,
not a generic property.

### 3.2 REFUTED — harmonic demixing as a boundary mechanism

Proposal: overlap cost depending on the harmonic relation between cells
(`U ∝ D(w_i/w_j)·V_overlap`) gives differential repulsion — the
Flory–Huggins/Ising demixing mechanism, which at zero temperature demixes
*by construction* and yields a firm interface.

**Killed on counting, from data already in hand:**

* A domain needs consonant-bond probability above the bond percolation
  threshold `p_c ≈ 1/(z−1) ≈ 0.167` at z=7. Measured **scalar** lock
  fraction is 0.05–0.13 — already at or below threshold. Requiring three
  axes to match simultaneously gives ~p³ ~ **1e-3**, two orders below; a
  cell then has ~0.007 consonant neighbours. **No connected domain
  exists, so there is nothing for a boundary to enclose.**
* Signal-to-noise: the harmonic term modulates contact stiffness at scale
  β·D·ε_contact, while the foam's geometric disorder is σ_d = 3%
  annealed, 18–19% frozen. The preference sits under the disorder it
  would have to organise against.

### 3.3 The dimensional counting that unifies all of §3

| | constraint structure | outcome |
|---|---|---|
| **within** a cell | ONE closed triangle, 2 free intervals | satisfiable — 26 shapes |
| **between** cells, size | d simultaneous scalar conditions | 5–13%, sub-percolation |
| **between** cells, 3-axis | 3×3 relations across d neighbours, ~p^(3d) | ~1e-3 |

**Harmonics are a per-cell structure, full stop.** They cannot mediate
between cells in 3D — not as a bound, not as a boundary. Going to three
axes makes the inter-cell case *worse*, because it cubes a matching
requirement that already failed.

### 3.4 NOT REFUTED — WITHDRAWN. Localisation has never been tested.

Three rejections of localisation are **withdrawn** because all of them
ran below jamming:

| L | φ | jammed? | used by |
|---|---|---|---|
| 20 | **0.139** | **no** | G2, G4 (N16), G5, G6 (N17) |
| 12 | 0.643 | yes | G1 only |

At φ=0.139 there is ~4.6× more free volume than at jamming; a hot core
disperses whatever the contact law is, because there is room to.
Therefore:

* **G2** "pressure cannot localise — convexity" — *not sound as measured*
* **N16** "`cap` saturation opens no coexistence" — *not sound as measured*
* **N17** "α<2 blows the packing apart" — unsupported; also used the wrong
  function family (a power law with α<2 has force → constant k as δ→0,
  i.e. a hard shell, not a yielding contact)

The **arguments** survive independently — convexity does imply a single
phase, and the derivation below is analytic — but no experiment has
tested localisation in a regime where it was possible.

### 3.5 The one analytic result worth keeping from that line

For a cell in a fixed cage with a power-law contact `U = Σ(k/α)δ^α`:

    mu = (1/3) E^(-2/3) P,  P = sum_j k delta^(alpha-1)
    mu falls with E   <=>   delta/r > (alpha - 1)/2

The non-convexity localisation needs is **already present in geometric
pressure** — it is out of reach. Harmonic contacts (α=2) need 50%
overlap against a jamming overlap of 1–10%.

This also refutes N16 *on paper*: the criterion contains α and nothing
else, so `cap` softening — which rescales `k` — could not move the
crossover at any value.

### 3.6 No GPU speedup (`localclock.cu`, `charge/runs/localclock_cuda.log`)

V100, correctness confirmed (max|gpu−cpu| = 8.882e-16 = 4×DBL_EPSILON, a
`sin()` library difference, not a scheduling one). Speedup 1.6× at
N=1536, 6.7× at N=6144 — then it stops.

**Batch width saturates at ~100**: 58.5 → 86.7 → 98.6 → 102.4 for
N=1536…98304. My earlier CPU extrapolation (~N^0.8 to ~1e4 at N=1e6) was
**wrong** — that was the rising part of a curve that plateaus. At
N=98304 the run is 199,198 *sequential* rounds each with three kernel
launches and two blocking readbacks; at ~30 µs/round that is ~20 s of the
measured 22.5 s. **Latency-bound, not throughput-bound.**

The limiter is the **local-time eligibility window**: with lookahead 0.05
and steps ~0.02, only events in a narrow time band are eligible at once,
and that band holds a bounded count however large the graph. The fix is a
modelling question (lookahead is physical — the rate limit), not a kernel
one.

---

## 4. PORTING THE BATTERY — what stands in the way

The battery is 20 experiments + the ħ-linearity cross-check, currently
**19/20 gated + 1 recorded** on `laws_V2g` (`battery/README.md`;
`e3b_blob_tilt` moved to recorded-not-gated on 2026-08-02 under the seed
panel).

### 4.1 The blocker, stated plainly

**On a frozen foam, a dense blob is held together by the lattice.** Cells
cannot move, so a dense region stays where it is. On free cells nothing
holds it, and **no boundary mechanism exists**: harmonics are refuted
(§3.1–3.3), uniform pressure is convex (§3.5, argument sound), and
localisation at jammed density is untested (§3.4).

Every battery experiment that depends on a self-maintaining dense object
is therefore blocked. That is the single gating item.

### 4.2 Tiered port

**Tier A — portable now** (no self-maintaining dense object required):

| exp | what it needs |
|---|---|
| `e1_conserve` | the ledger only. Local clocks already give FP-floor conservation. |
| `e2_pulse` | field-sector propagation on a dynamic contact graph |
| `e5_bell` | CHSH on a pair — small, local |
| `d1_slit`, `t1_tonomura`, `q2_eraser`, `t4_hom` | optics; need a stable *field* path, not a dense object |
| `qt_lo`, `qt_hi` | conversion threshold at a point |
| `p1_beam` | momentum of light |
| `LIN` | ħ-linearity: an invariant over whatever fires |

Risk for the optics group: they depend on the **connectivity floor**
(F2: mean degree ≥ ~7 for optics). A free-cell substrate must maintain
that, and livefab run 1 densified catastrophically (degree 8.6 → 47) when
links were re-derived under a spring — so link existence must be a
**contact rule with energy-dependent radii relaxed to jamming**, not a
spring over a cutoff (`LIVEFAB.md` §6).

**Tier B — blocked on localisation:**

| exp | why |
|---|---|
| `e3a_blob` | a heavy blob must *seal* — i.e. hold itself |
| `e3b_blob_tilt` | tilted blob must translate without dispersing |
| `g1_footprint` | a mass must *maintain* an extended depression |
| `g3_shadow` | matter must be opaque, so matter must persist |
| `g4_throughput` | a blob's space flux over time |
| `e4_curve` | curvature response to a mass |

**Tier C — uncertain, test early:**
`e6_pairs`, `e7_tune`, `e8_comma`, `e9_fifth` need dense *pairs* that
persist long enough to tune. They are small (individual tuned cells
rather than blobs), so they may survive without a full boundary
mechanism — but `gamma_res_m=0` (the dense rim seal) already destroys the
pair sector wholesale on the frozen foam (tuned pairs 31 → 7, gg 0.684 →
0.126, e7 fails on **all five** panel seeds, comma unpaid at exactly
0.000). If the seal cannot be maintained on free cells, Tier C falls with
Tier B.

### 4.3 Kernel changes required (all need explicit user authorization)

Per `LIVEFAB.md` §5, extended by this thread's findings:

1. **Apparatus key `live_links`, default 0 = byte-identical legacy.**
2. **Birth/death ledger** replacing the per-step `if (A > 0)` liveness
   test. A channel whose overlap drops to zero is marked dying, finishes
   its in-flight cycles under **rule (α) flush-on-death** (a channel may
   die only at `e_mid = 0`), then its slot recycles. Birth is free
   (`e_mid = 0, φ = 0`).
3. **Per-step spatial-hash neighbour query** replacing the init-time dart
   throw.
4. **Local-clock execution** with the four conditions of §2.1 —
   critically, **channels become first-class event holders with their own
   clock**, which is a structural change to how transfers are fired, not
   a scheduling detail.
5. **Cell extent from energy** (`r ∝ E^(1/3)`), and — if shape is adopted
   — a 3-axis shape per cell with volume from E and axis ratios free.

### 4.4 Acceptance criteria

* `laws_V2g` **19/20 gated + 1 recorded**, and **byte-identical output
  with `live_links=0`** (the ratchet's rule 1 comparison).
* Conservation at the FP floor unconditionally — this is the One Law and
  is not negotiable.
* Bit-identical reruns at any thread count (the property the battery's
  baseline comparison depends on).
* Mean degree ≥ ~7 maintained (F2 optics floor), σ_d not worse than the
  frozen baseline it replaces.

---

## 5. WHAT WOULD HAVE TO BE OVERCOME

In dependency order.

**(a) Test localisation at jammed density.** This has never been done.
`geopress` with φ ≥ 0.64, z ≈ 6–7, retention measured at **two
relaxation times** with a convergence witness and a **packing-fraction
column printed alongside every result** so a dilute run announces itself.
Until this is run, "pressure cannot localise" is an argument, not a
measurement.

**(b) Give cells shape.** A scalar radius is a volume, not a morphology.
Under anisotropic load such a cell can only change volume, i.e. shed or
absorb energy — **dispersal is the only channel it has**. Three axes
provide a constant-volume soft mode and, decisively, **shear
resistance**; without shear resistance no assembly can hold against
directional stress. `shapelock.c` shows the shape comb selects discretely
and bounds aspect ratio; it does **not** show a network of such cells
holds together.

**(c) Find a boundary mechanism, or prove none is needed.** Ruled out:
harmonics (percolation, §3.2), uniform pressure (convexity, §3.5).
Remaining candidates are **structural rather than harmonic** — contact
topology, or **v89's own gated-transport machinery**, which is where the
dense sector's physics actually lives and which this entire thread never
touched. Note the law-dependency map: `gamma_res_m` (the dense rim seal)
is what lets two dense objects hold a relationship at all, and `w2` (the
dense/field pitch separation) is the most load-bearing knob of the eight
(14/20). Neither has been examined as a cohesion mechanism.

**(d) Recover parallel width, or accept serial.** The batch saturates at
~100 because of the local-time eligibility window. Raising the lookahead
is a modelling decision (it is the physical rate limit), not tuning.
Removing the two per-round device syncs is worth ~10× on latency but
cannot create absent width.

---

## 6. STANDING METHOD RULES (earned the hard way)

Six times in this thread a sweep or a proposal preceded a cheap analytic
check that would have settled it. Recorded so a successor does not repeat
them:

1. **Derive the criterion, confirm the configuration satisfies it, then
   sweep.** The percolation check (§3.2) is one line against numbers
   already in `harmfrust`'s log. The packing fraction (§3.4) is one line.
   The α-independence that refutes N16 (§3.5) is two lines of algebra.
2. **Print the configuration alongside every result** — packing fraction,
   mean degree, contact count. A dilute run must announce itself.
3. **A plateau over 2× is not convergence.** N16 read 0.192 → 0.194
   across a 2× increase in relaxation time and then collapsed to −0.262
   at 8×. Report the **convergence witness** (max|f|, max|grad|)
   alongside, not the plateau: at 960 cycles the converged case sat at
   6.2e-14 while the false one was at 2.3e-3, ten orders away.
4. **Rugged min-over-rungs landscapes need rung-seeded search**, not
   random restarts. Random restarts report the basin nearest the start
   as if it were the answer, and produced two spurious staircases here.
5. **A drive must carry a scale.** A linear/scale-free drive has
   unbounded benefit and always runs to the outermost rung, so the rung
   structure selects nothing. This error was diagnosed, fixed in
   `harmgrow.c`, and then reintroduced in `shapelock.c`.
6. **Verify the negative control actually runs.** `build_rungs(climb<0)`
   returned *zero* rungs, making the energy 3e300 and the relaxation a
   no-op; the control read a perfectly flat response and looked like a
   pass. In `harmgrow.c` the same bug was invisible because `RUNGS` is
   static and zero-initialised, so `RUNGS[0]` *is* the unison and only
   the gradient was consulted. **A control that passes by accident is
   worse than no control.**
7. **Booleans hide artifacts.** The first GPU run printed `EXACT/DIFFERS`
   and was uninformative; printing the magnitude showed 8.882e-16 = 4 ulp,
   a library difference rather than a bug.

---

## 7. FILE AND LOG INDEX

| artifact | what it establishes |
|---|---|
| `localclock.c` / `charge/runs/localclock.log` | local-clock readiness suite R1–R7; batch-parallel bit-identity; thread invariance |
| `localclock.cu` / `charge/runs/localclock_cuda.log` | GPU batch; correctness to 4 ulp; width saturation; latency bound |
| `geopress.c` / `charge/runs/geopress.log` | G1 size bound (**valid**); G2/G4/G5/G6 localisation (**withdrawn — sub-jamming**) |
| `shapelock.c` / `charge/runs/shapelock.log` | S1 shape spectrum (exact); S2 shape staircase; S3 control; S4 aspect bound |
| `harmfrust.c` / `charge/runs/harmfrust.log` | consonance between cells **refuted**, degrees 2–12 |
| `harmgrow.c` / `charge/runs/harmgrow.log` | ring size staircase — real but **one impurity in a uniform bath** |
| `DESIGN_2026-08-02.md` | the full running record, including every withdrawal |
| `ROADMAP.md` §8 | standing decisions: local clocks are the execution model; consonance rejected |
| `LIVEFAB.md` | dynamic-substrate design; contact rule; death rule (α); cellfab integration path |
| `battery/README.md` | ratchet rule, seed panel, the complete 8-knob law-dependency map |

---

## 8. ONE-PARAGRAPH SUMMARY FOR A MODEL STEPPING IN

The execution half is solved: local clocks pass a four-criterion
readiness suite (conservation 1.2e-16 relative, first-order convergence,
bit-identical determinism, full event width) and batch-parallel execution
is provably identical to serial. The substance half is not: a free-cell
fabric has **no mechanism to hold a localised object together**.
Geometric pressure bounds cell size robustly (that result is solid) but
cannot localise by a convexity argument; consonance is refuted between
cells by percolation counting but works *within* a cell, giving a finite
shape spectrum and a shape staircase. Crucially, **localisation itself
has never been tested at jammed density** — three prior rejections all
ran at φ=0.139 against φ_J≈0.64 and are withdrawn. The first task is to
run that test properly, with cells that have three axes (hence shear
resistance) rather than a scalar radius, and with the packing fraction
printed next to every number.
