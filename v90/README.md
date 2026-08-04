# v90 — the carried set

v90 carries forward from v89 **exactly three working methods and nothing
else** — the free-cell substrate, local-clock execution, and the battery
discipline — with **C as the kernel of record** and **Go replacing python**
for the harness, viewer, and orchestration (plus a verified Go-kernel
experiment beside the C one). New physics happens on this carried set.
Everything else in v89 stays where it is, as the record and the reference.
The pre-v89 ban is unchanged and permanent.

| carried method | v89 evidence | v90 home |
|---|---|---|
| **the free-cell substrate** — cohesion from the standing law on cells whose geometry is state | `v89/freecell.c`, `v89/FREECELL.md` §9–§10 (substance solved in two regimes; FCQ) | `kernel/freecell.c` (**C, kernel of record**) + `fab/` (Go, verified experiment) |
| **local-clock execution** — no global tick, conservative, deterministic, batch-parallel bit-identical | `v89/localclock.c`, `v89/FREECELL.md` §2 (four measured conditions) | the P2 scheduler spec (below) |
| **the battery discipline** — one law table, apparatus-only configs, gated bars, the ratchet | `v89/battery/` (laws_V2g 20/20; python harness) | `cmd/battery` (**Go harness — python replaced**) |

## The law

**laws_V2g, VERBATIM, unchanged.** `v89/battery/laws_V2g.cfg` remains the
canonical table; `fab/params.go` embeds the same constants as defaults and
the battery purity-checks that they have not drifted. v90 introduces **no
law change**. Apparatus keys (geometry, seeds, meters) are per-experiment;
law keys are the table.

**The ratchet rule governs v90 exactly as v89:** every kernel or law-table
modification runs the FULL v90 battery before commit; experiments that pass
and defend a claim join the suite and gate all future modifications; bars
encode claims (sharpen by measurement, never soften to pass); a green test
leaves the suite only by explicit user decision. The v90 battery compares
bytes where it claims determinism.

## Languages and roles (user decision, 2026-08-03)

**C is the kernel of record** (`kernel/freecell.c`, carried from
`v89/freecell.c`): production physics runs here, and the GPU path (CUDA)
stays open through it — Go does not run on GPUs. New apparatus lands in the
C kernel first, battery-gated.

**Go replaces python as the harness language** (the priority): the battery
(`cmd/battery`), the viewer (`cmd/volview`), and all orchestration are Go.
No python in v90.

**The Go kernel (`fab/`, `cmd/fabrun`) is a sanctioned experiment**, not the
production kernel: a 1:1 port kept verified against C (`VERIFY.md` — the
trig-kernel path is bit-identical arithmetic by construction; whole bath and
blob runs reproduce byte-identically). What it buys: a second independent
implementation that catches kernel bugs by divergence, run-to-run
byte-determinism in a memory-safe language, and the natural host for
prototyping the localclock goroutine scheduler (P2) before it is written
into C. Measured cost: ~1.3× C wall time on the blob box.

## Layout

| path | what |
|---|---|
| `kernel/freecell.c` | **the C kernel of record** (carried from `v89/freecell.c`; the v89 original stays frozen as evidence). New apparatus lands here first. |
| `fab/` | the Go kernel experiment — a 1:1 port (passes 0,S,1,G0,2,3,4,5,D,F,6,7; slot ledger; rule-α birth/death). Deviations listed in `VERIFY.md`. |
| `cmd/fabrun` | Go kernel CLI, same `key=value` config surface and log-line formats as the C instrument (so runs diff textually). |
| `cmd/battery` | **the battery harness (Go)**: kernel-agnostic — runs `./freecell` (gate) and `./fabrun` (cross-check), parses `# RESULT` lines, evaluates gated bars, writes `runs/*.log`. `./battery` must end `ALL GREEN`. Heavier experiments also stream `.fcs` to `runs/streams/` (regenerable). |
| `cmd/volview` | volumetric viewer for `.fcs` snapshots (basis: `sfa/volview`). Headless renders (`-view ... -out f.png`, `-avg` time-average), **`-info`** (stream inventory, CFG provenance, ANLZ tables as TSV), and the **interactive mode `-i`**: orbit rotation, time scrub/play, channels by key (es em ee x r **phase fa1 thd** + any schema column), link rendering, tag tint, per-frame/global normalization, screenshots; **`-follow`** attaches to a stream a kernel is still writing (live view, pin-to-latest). `volview -i run.fcs`, `volview -i -follow live.fcs`. |
| `FCS.md` | the snapshot format spec (v3 chunked streams; v1/v2 legacy) and why it is not SFA (the background assumption, spelled out). |
| `runs/` | committed evidence logs (`BATTERY.log` = latest full-suite verdict). |
| `VERIFY.md` | the A/B record: C↔Go agreement protocol and measurements; byte-determinism. |
| `DS.md` | the double-slit-on-free-cells campaign record (P1). |

Build: `cd v90 && make all` (gcc + go ≥ 1.27). No external dependencies.

## Program

**PRIMARY FOCUS (user directive, 2026-08-03): the composite programme —
`COMPOSITE.md`.** Many-celled quarks (composites whose collective mode
carries flavor), validated by many-cell experiments, with **Mass, Charge,
and EMF as parallel lanes** reading the same runs. First execution done:
the composite bond works (and carries its energy as flight), the
macro-parity selection rule quantizes composite geometry (nv ≡ 0 mod 12),
the U/D molecule/droplet asymmetry is real, and all three lanes converge
on the flight-load-vs-tongue frontier at three scales. P0–P4 below feed
this focus (P2 scaling = bigger composites; P3 shells = composite
response structure).

**P0 — port + verify (the gate for everything else).**
Faithful Go port of the kernel; A/B against the C reference: t=0 INIT
identity, short-horizon trajectory agreement at the libm-ulp floor,
long-run structural agreement (conservation floor, pair rung, ring parity,
pulse speed, blob retention); Go run-to-run byte identity. Record: `VERIFY.md`.

**P1 — verify the standing FREECELL states, then re-host the reality ladder.**
The v90 battery carries the measured v89 free-cell claims as gated bars
(conservation, determinism, pair rung both-sides, ring6 vs ring5 parity,
pulse v/C, blob retention live≥frozen, UUD closed vs UDD open-chain). Then
the first NEW verification: **the double slit on the free substrate** —
tier 0 (wave fringes at parameter-free loci, `DOUBLESLIT.md` conventions)
on a jammed free-cell medium with a carved vacuum barrier; frozen-scaffold
control beside it. Honest scope: tiers that need single-quantum clicks ride
the atoms machinery (carried); anything needing sharp flavored matter waits
on the S2 amplitude completion (v89 FREECELL §10.4).

**P2 — larger areas.**
The spatial hash is O(N) and the box is periodic — the port already scales
mechanically; what larger areas *need* is the local-clock scheduler so quiet
regions do not pay global-tick cost. Implementation spec = FREECELL §2's four
measured conditions, verbatim: (1) order events by a total function of state
`(t, kind, index)`, never arrival; (2) bound skew in LOCAL TIME, never tick
count; (3) channels own their transfers and their own clock; (4) never order
by the tick counter (assert `K < M/2` on any cyclic byte). Batch-parallel =
minimum over conflict neighbourhood, bit-identical to serial by measurement.
Acceptance: same battery, byte-identical where claimed, plus a scaling table
(cells vs ms/t.u. vs RSS).

**P3 — emergent electron shells (the north star).**
Shells = bound response harmonics around a bound core (the v85 shell thread's
surviving form), not an imported species. On free cells: a rung-bonded core
object embedded in the live substrate; look for stationary field-mode
response structure (radial nodes in time-averaged Ee) that appears without
being seeded and survives perturbation. Prerequisite honesty: sharp charge
(flavored composites) is amplitude-completion territory; unison-core response
structure is testable at rate level first.

**P4 — numeric speedups, only after stability.**
The v89 bench discipline (BENCH.md): standard boxes, measured per-rung
speedups, byte-identical outputs across thread counts or the rung is
rejected. The C kernel is the speedup target (OpenMP, then CUDA — the GPU
path is why C remains the kernel of record); the Go kernel prototypes the
localclock batch schedule where goroutines make it cheap to try. `simd`
package on hot loops is a Go-side candidate. No approximation rungs — the
law is not a knob.

## Stretch goals (user-requested, 2026-08-03)

* **Pauli exclusion — design recorded (2026-08-03), blocked on S2.**
  v89 tier-4 HOM measured fermionic exchange signs (g_b 0.42 < 0.5 <
  g_f 0.58, exchange registry). Exclusion is stronger: identical
  excitations refusing one state, distinct from cap saturation (which
  refuses everyone). **The design:** (1) observable = joint occupancy
  statistics at a Y-junction cavity — P(both grains in one cell) for
  *identical* vs *distinguishable* (different-pitch) grain pairs, all
  else equal; Pauli = identical pairs anti-bunch in OCCUPANCY, not just
  in coincidence. (2) Prerequisite: two-quantum amplitudes with an
  exchange sign — the S2 completion's coherent channel plus an exchange
  registry; a classical two-plane ψ cannot exclude (fields add).
  (3) Rate-level precursor, runnable TODAY: verify cap-saturation
  refusal is distinguishability-BLIND (same refusal for identical and
  different-pitch arrivals) — the negative control that later separates
  exclusion from fullness. **RUN AND GATED (2026-08-03, battery
  `pauli0`, 8 bars):** on a frozen contact channel, a receiver AT cap
  admits EXACTLY zero from an identical-pitch sender and from a
  fifth-consonant sender alike (the headroom factor multiplies the
  transport want to zero before the consonance gate enters — and the
  door is FP-sticky); NEAR cap (headroom 0.04·cap) both kinds are
  admitted, throttled (intake 0.106 vs 4.0e-3 — the 26× is consonance
  gate quality, the transport's business, not the door's). Saturation
  refusal carries no identity information at rate level. Apparatus
  lesson: "distinguishable" arrivals must still be CONSONANT (a fifth,
  w_snd = 1.5·w_rcv via load) — arbitrary detunes do not transport at
  all; and the contact channel needs an open lens (d < rᵢ+rⱼ), since
  loads shrink radii (Es = 1 − Em·s_pull/(1+s_pull)).
* **Shells (incomplete harmonics), visualized.** Time-averaged field
  structure around a bound core, rendered volumetrically: `.fcs` sequences
  → volview time-average mode → radial node structure, if it exists, seen
  directly. Pairs with P3 (the shell probe) and needs no new physics to
  *attempt* — a bound core in a live bath and patience.

  **First look (2026-08-03, `runs/shell/`):** the toolchain works
  end-to-end (`volview -avg`), and the rate-level physics sorts cleanly:
  a **unison ring6** in a bath radiates *exactly zero* field (det=0 ⇒
  roughness silent; and sub-atom demands never fire — the atoms machinery
  keeps consonant matter dark) while leaking dense energy into a visible
  Em atmosphere (`shell_ring6_avg_em.png`); the **seeded-exact UUD
  triangle is also silent** (its fifth is exact by construction —
  consonant); the **blob radiates** (rough 1.93 in 60 t.u. — many random
  detunings) and shows a smooth time-averaged **field halo**
  (`shell_blob_avg_ee.png`) — diffuse, no radial nodes at this apparatus.
  So at rate level the field atmosphere is sourced by *detune*, not by
  binding; node structure, if it exists, likely needs a driven/response
  setup (P3) or the amplitude completion.

## Working rules (carried verbatim)

Derive the criterion before the sweep; print the configuration beside every
result; controls must verifiably run (print their config too); magnitudes,
not booleans; drives carry scales; pre-register predictions before first
run. State where every model's state lives and check it against constraint
2 (no background) before writing it.

## Status

**The standing queue and interpretation live in `NEXT.md`** — one page:
hard results index, what they mean, and the prioritized next steps.

| item | state |
|---|---|
| P0 port + verify | **done** — VERIFY.md: bath/blob/blob2 byte-identical, physics digits identical everywhere; queue-#6 ports restore the FULL A/B surface (slit/rings/blob2 + all meters; battery `abx`) |
| P1 battery | **GREEN — 85 bars**, both kernels, through 14 apparatus additions (runs/BATTERY.log) |
| P1 double slit tier 0 | **CONFIRMED on the live substrate** (DS.md; battery `ds` 8 bars) |
| DS tier 1 | **CLAIMED** (DS.md; battery `ds1` 8 bars) — 362 whole-grain clicks, fringe-phase-locked s̄ 0.683 vs controls 0.380, minima 7.7× dark, R3 atomicity exact |
| FCS + viewer | v3 chunked streams, interleaved ANLZ, lossless codec, live `-follow`, interactive orbit/scrub (FCS.md) |
| MOTION baseline | measured: bounded wander, plastic-lock collisions, blob merging, driven-object null (MOTION.md) |
| **XSEC cross-section** | **gated (15 bars)**: headroom=absorber (+7.27, clean cond), flat-top-cap=EMITTER (−7.2, side-glow), inert load=reflector; monotone impact-parameter profile; absorption is clock-rate, not mass (MOTION.md XSEC) |
| **COMPOSITE (primary focus)** | rounds 1–2 executed: macro-parity rule, T6 U/D asymmetry, bond=flight, compensation 8×, H-stiffness→flux-machine interior, **bath protects**, nucleon residual force (COMPOSITE.md) |
| P2 local-clock scheduler | **Go prototype done** (fab/localclock.go; battery `p2lc` 8 bars): four conditions on the real substrate, batch==serial exact at 1–8 workers (pending-min rule), quiet-region economy 6.87× at L=64 (P2.md); C production engine open |
| P3 shells | first-look done (runs/shell/); driven-response probe queued |
| P4 speedups | not started (after stability); P2 production engine is the first rung (P2.md checklist) |
| **S2 frontier** | **entry measured + design memo** (S2.md): kr=1 fails — kills the fifth everywhere load-bearing; completion = coherent dense channel, quantified acceptance surface; [LAW]-level, needs sign-off |
| **REALITY audit** | **done** (REALITY.md): correlation-to-reality retrospective with degrees (EXACT/QUANTITATIVE/STRUCTURAL/ABSENT/ANTI) and the prioritized adjustment list — the refreshed queue in NEXT.md derives from it |
