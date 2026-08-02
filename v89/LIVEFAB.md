# LIVEFAB — the dynamic-topology substrate (S3 design + prototype)

**Design document and prototype, opened 2026-08-01 by GLM during review.**
Subordinate to `PRINCIPLE.md` and `SUBSTRATE.md`. The keystone next step:
this is what converts the program's central claim ("no background") from
*tested on a stage-sanctioned instrument* to *tested background-free*.

Kernel policy respected: `cellfab.c` is **not** modified. The prototype is
`livefab_proto.c`, a standalone file. Any `cellfab` integration is flagged
for explicit user authorization.

---

## 1. The problem, restated sharply

`CONSTRUCTION.md` §0 forbids any of:

1. labels that persist while values change,
2. a permanent reference geometry,
3. **fixed connectivity not itself produced by the current energy structure**,
4. a quantity whose conservation requires fixed ambient space.

`cellfab` violates **(3)**: the candidate channel set and channel *lengths*
`ld[l]` are generated once at init and persist. Only the contact *areas*
`lA[l]` are live (functions of current radii). The README and `CELLFAB.md`
§8 flag this as a "sanctioned residue"; `SUBSTRATE.md` makes the measured
argument that it cannot stand:

> **F1 — the disorder is topological.** Frozen-graph spring relaxation
> stalls at σ_d ≈ 18% regardless of sweeps: a dense random graph cannot be
> equalized positionally. Uniformity requires either window control (S1's
> contact shell) or **dynamic topology (S3).** This is livefab's measured
> justification.

S1 (annealed contact-shell, σ_d 19%→3.0%) is a frozen-window control: it
removes the disorder but keeps the topology frozen. It passes the battery
at 18/21-effective and needs apparatus re-derivation per substrate (F3:
addresses are substrate-dependent). **S1 does not satisfy constraint (3).
livefab does.**

## 2. The livefab thesis

Drop the frozen channel set entirely. A channel is **born** when two
parcels' areas of effect overlap *now* (a function of the current energy
state, via the live radii `r = r0·(E_s/e_s0)^{1/3}`), and **dies** when
they cease to overlap, or when resonance fails persistently, or when an
endpoint parcel ceases. The channel set at every instant is a product of
the current energy complex — nothing frozen.

The channel *length* is then the current geometric distance between two
existing parcels, recomputed when the channel exists. There is no `ld[l]`
array that persists; the length is a property of the current pair.

**Falsifiable sub-thesis (what the prototype tests):** a system whose
constraint graph is re-derived each step from the current state can
equalize where a frozen-graph spring stalls. If live-links relaxation
also stalls at ~18%, the livefab thesis is in trouble and S3 needs a
different mechanism. If it reaches the S1 ~3% band *and* stays connected,
the keystone works.

## 3. Conservation across topology change — the load-bearing detail

This is where livefab either lives or dies, and it is independent of the
σ_d question. When a channel **dies** with `e_mid > 0` (energy in transfer
mode T, mid-cycle), that energy cannot vanish (PRINCIPLE §1) and cannot
strand on a dead channel (no ghost). Three options, in order of elegance:

* **(α) Flush-on-death.** A channel may die only at `e_mid = 0`. Any
  channel with in-flight energy is kept alive until its cycle completes or
  its energy returns to an endpoint. Death is *gated* by the cycle gate,
  not imposed. **Cleanest; preferred.** Cost: a few zombie channels linger
  finishing their cycles; bounded by the in-flight population.
* **(β) Return-to-source.** On death, `e_mid` returns to the source
  parcel as space (mode S). A paired ledger move; conserves exactly.
  Cheap but physically odd (the energy "un-deposits").
* **(γ) Split to both ends.** `e_mid` splits half/half to the two
  endpoints. Conserves; symmetric. Odder still.

(α) is the choice because it makes topology change *subordinate* to the
cycle gate — the same gate that already governs everything else. No new
conservation rule. **Channel birth** is free (a new channel starts at
`e_mid = 0, φ = 0`); only death needs the gate.

The integer ledger (`int_ledger/`) already moves `e_mid` on per-channel
slots; under (α) a channel's slot is simply not recycled until `e_mid`
hits zero. The books close unconditionally.

## 4. The prototype — `livefab_proto.c`

Standalone; no `cellfab` dependency. Three arms on identical initial
conditions, measuring σ_d and degree vs sweeps:

| arm | mechanism | what it tests |
|---|---|---|
| **A** `frozen-spring` | fix neighbor graph at init; spring each edge toward `dtar`; iterate | reproduce F1's 18% floor (the failure case) |
| **B** `live-links` | each sweep, re-derive links from current positions (`d < rcut`); spring each **live** link toward `dtar`; over-stretched links die | **the livefab thesis** — does dynamic topology equalize? |
| **C** `s1-repulse` | pure repulsion below `dtar` (no attractive term), then trim to first shell | reproduce S1's ~3% (the control) |

If **B** reaches σ_d well below A's floor and near C's band **while keeping
mean degree ≥ ~7** (the F2 connectivity floor for optics), livefab works
mechanically. The prototype does not carry the full gate physics — it
isolates the one contested claim (can dynamic topology equalize?).

Build: `gcc -O2 -o livefab_proto livefab_proto.c -lm`. Run: `./livefab_proto`.

## 5. After the prototype (flagged for user authorization)

If the prototype is positive, the path into `cellfab` is:

1. A new apparatus key `live_links` (default 0 = byte-identical legacy).
2. Replace the per-step `if (A > 0)` channel-liveness test (already
   present — `lA[l]` is recomputed each step from current radii) with a
   **birth/death ledger**: a channel whose overlap drops to zero is marked
   dying; it finishes its in-flight cycles under (α); its slot recycles.
   A channel between two currently-overlapping parcels that is not yet in
   the table is born at `e_mid = 0`.
3. The candidate-pair search (currently the init-time dart throw) becomes
   a per-step spatial-hash neighbor query over the *current* positions.
4. Re-run the full battery. The success bar is `laws_V2g` 21/21 on the
   live substrate (SUBSTRATE's standing bar).

This is a kernel edit and **requires explicit user sign-off** before it is
written into `cellfab.c`. The prototype above is the evidence-gathering
step that does not need it.

---

## 6. Prototype verdict (run 2026-08-01) — the thesis holds

`livefab_proto.c` loads the S1 annealed substrate (NC=5039, σ_d=8.4% over
the 1.15·(rᵢ+rⱼ) link rule), inserts a matter cluster (shrink radii in a
central sphere — "space converted to pattern"), and relaxes under pure
repulsion comparing a FROZEN (vacuum-graph) substrate against a LIVE
(re-derive links each sweep) substrate. Steady state at sweep 300:

| matter strength | region | FROZEN σ_d | LIVE σ_d |
|---|---|---|---|
| strong (r×0.30) | far field | 14.51% (grew from 8.4%) | **4.71%** (below baseline) |
| strong (r×0.30) | core | 24.79% | 21.24% |
| mild (r×0.55) | far field | 14.56% | **4.68%** |
| mild (r×0.55) | core | 23.58% | 19.16% |

**Reading.** The livefab thesis is mechanically confirmed and the
discriminator is clean: the LIVE far field re-jams to *lower* disorder than
the frozen baseline (8.4%→4.7%) — operationally "the remaining space
accounts for the converted region" (PRINCIPLE §4.3) — while the FROZEN
substrate frustrates and its far field degrades (8.4%→14.5%). Live
re-derivation wins in both regions; the dramatic margin is the far field
(~3×). The core is necessarily disturbed in both arms (matter is a real
disturbance: "where there is matter there is less space"), and is not the
discriminating quantity.

Run 1 (recorded as a finding, not a failure) established the negative
result that *defines* what livefab must be: a "live geometric-cutoff links
+ attractive spring" reading **densifies** (degree 8.6→47, σ_d rises),
reproducing SUBSTRATE finding #1. Livefab's link existence must therefore
be a **contact rule with energy-dependent radii** (d < 1.15·(rᵢ+rⱼ), r ∝
E_s^{1/3}), relaxed by **pure repulsion to jamming** — not a spring over a
geometric cutoff. Under that reading the vacuum is continuous-S1 and the
defining capability is the matter-response measured above.

**Status of the keystone.** The mechanism works in isolation. The
unresolved question is whether carrying it into the full gate physics
(cycle-gated transport, the dense comb, in-flight energy on dying
channels under rule (α)) stays stable and passes the battery at 21/21.
That is a `cellfab` edit and awaits explicit user authorization per §5.

