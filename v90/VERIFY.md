# VERIFY — the C↔Go kernel A/B record

**Claim:** `fab/` (Go) is a faithful 1:1 port of the kernel of record
`kernel/freecell.c` (C, = `v89/freecell.c`). "Faithful" is defined and
measured below, not asserted.

## Construction rules that make the comparison sharp

* Same pass order (0, S, 1, G0, 2, 3, 4, 5, D, F, 6, 7), same serial
  semantics, same canonical orderings (cells ascending, neighbour lists
  sorted, CSR sorted by other-endpoint, hops applied from the i side).
* Same RNG (xorshift64, identical warm-up and draw order — the draw
  *sequence* is part of the port; e.g. the skipped th1 draw is preserved).
* The gate/phase trig (`lut_cos`/`lut_sin`) is pure double arithmetic
  (fmod + minimax kernels, `fab/lut.go` = `int_ledger/trig_lut.inc`), so it
  is **bit-identical by construction**. The C reference is built `-O2` with
  no `-march` (baseline x86-64 has no FMA instruction), and Go's compiler
  does not contract on amd64 — so neither side fuses.
* Same log-line formats (C `%g` ↔ Go `%.6g`), so runs diff textually.
  Line 1 identifies the implementation; comparators skip it.

## Known, bounded divergence sources

libm calls where glibc and Go's `math` may differ in the last ulp:
`cbrt` (pass-1 radii), `sincos` (pass F), `pow` (the P15 bond gates and
non-integer `p_gate`), `log`/`cos` (Box–Muller: init + tumble),
`exp`/`sin`/`cos` (seed construction), `acos` (gyration diagnostics only).
`sqrt`, `fmod`, `floor` are exactly rounded on both sides and cannot
diverge. There are no other sources.

## Measurements (2026-08-03, this machine, seed 20260802)

| test | result |
|---|---|
| `exp=bath T=10` and `T=40` | **byte-identical** C vs Go after line 1 (incl. 2320/4611 channel births; the battery's `cross C==Go` bar re-checks T=40 every run) |
| `exp=blob T=20 L=12` and `L=16` | **byte-identical** after line 1 |
| `exp=pair bath=0 T=40/60` | every physics digit identical (`RESULT pair d_final` to all 6 decimals); only the conservation-drift column differs, at ±1–4·10⁻¹⁶ (the Kahan sum sees the libm ulps) |
| `exp=ring`, `exp=tri` | same pattern: physics digits identical; drift column at the 1e-16 floor; one −0.0000 vs 0.0000 print in a ~1e-17 gyration eigenvalue (acos, diagnostics only) |
| `exp=pulse L=24 T=10` | all columns identical incl. births/deaths/candidates and `Ee_front`; drift column ±5·10⁻¹⁶; `v_over_C` identical to 4 digits (battery cross bar: measured 0) |
| Go determinism | two identical invocations → **byte-identical logs** (battery bar, sha256) |
| law purity | both kernels print the laws_V2g header byte-exactly (battery bar) |
| cost | blob L=16 T=20: C 20.9 s / 20 MB RSS; Go 27.3 s / 48 MB RSS → Go ≈ 1.30× wall, 2.4× memory |

Verdict: the port is faithful at the floating-point floor. Cross-language
divergence is confined to the last ulp of a handful of libm calls and only
ever surfaces in the conservation column at ~1e-16 — never in a physics
digit at printed precision on any run tried.

## Queue-#6 ports (2026-08-04): slit / rings / blob2 + the meter apparatus

`exp=slit` (2D bath, carve, pinned wall fixture, packet, screen +
SCREENCELL, `slit_obj` + `obj_y`, `slit_clicks`/scond), `exp=rings`
(full composite machinery: molecule/droplet, rint charts, Z3 chains,
xcomp), `exp=blob2`, the `p1_meter` (5 instrumented sites), the XSEC
sector meter (`sect_*`, E + occupancy), and the tag-split conversion
ledger (`# CONVTAG`/`# RESULT convtag`) are now in BOTH kernels.
Battery experiment **`abx`** (6 bars) gates the pairs every run:

| test | result |
|---|---|
| `exp=blob2 L=16 T=80 p1_meter=1` (the battery p1 geometry) | **byte-identical** C vs Go after line 1, INCLUDING the drift column and the full `# P1`/`# RESULT p1` meter output |
| `exp=slit` tier-0 (ds_m0 args) | every physics digit identical — all diag counts, the full 96-bin `# SCREEN` table, every `# SCREENCELL` row, exposure to 6 decimals; only the drift column differs (±2e-15 C vs ±4e-14 Go — the packet seed's libm `exp/cos/sin` ulps; the optics field is LINEAR, so ulps never amplify) |
| `exp=slit` ds1 harvest args (`slit_clicks=1 amp=4`) | 22 clicks in both kernels, every `# CLICK t/y/e` row byte-equal, e_sum equal; drift column only (Go floor 1.17e-13 — the abx bar bounds it at 1e-12) |
| `exp=slit` XSEC hdr arm (law medium, occulter, sector meter) | `# RESULT convtag` net **7.274023 identical**, `# RESULT sect` Etot **961.008277 identical**, every SECT row byte-equal — the grain-quantized law absorbs ulp dust instead of amplifying it, even in the nonlinear fog |
| `exp=rings kind=1 nv=6 bath=0 T=100` | every physics digit identical (IR diag column, `# RESULT rint*`, `ringq*`); drift column at the 1e-16 floor; the known −0.0000 gyration zero-sign |

The `abx` comparator masks exactly the allowed divergence (the drift
column, the `RESULT drift_rel` line, negative-zero prints) and demands
byte-equality of everything else; blob2 is gated UNMASKED.

## Deviations (complete list)

* `.fcs` snapshot writer implemented in Go (`snap_every`/`snap_dir` keys —
  vestigial in the v89 C original, now implemented in the v90 C kernel as
  well). Pure output; consumes no RNG; cannot perturb determinism.
  Format **FCS v3** (`FCS.md`): both kernels write identical chunk
  streams (CFG/SCHM/CELL/LINK/ANLZ; compressed chunks differ in zstd
  bytes across kernels — parity is of decoded content).
* Log line 1 identifies the implementation.

## Standing rule

Any change to `fab/` or `kernel/freecell.c` reruns `./battery` (the cross
bars gate the pair). If the two kernels ever disagree in a physics digit,
that is a bug in one of them — the disagreement is the detector.
