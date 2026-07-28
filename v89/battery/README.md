# The unification battery — ROADMAP §6.5, operational

One **law table** (`laws_*.cfg`), byte-shared by every experiment; one
**apparatus file** per experiment (`apparatus/*.cfg`) that may not touch a
law key (the runner refuses to start otherwise); **physics acceptance**,
not byte comparison. Law variants compete as whole tables — switching a
law per experiment is structurally impossible here.

```
python3 battery.py --laws laws_V2.cfg [--jobs 10] [--only e6_pairs ...] [--skip-run]
```

Outputs `runs/<variant>/<experiment>.log` + `summary.tsv`, prints the
pass/fail table, exits nonzero on any fail.

## The experiments (16 runs + 1 cross-check)

e1 conservation (all mechanisms live) · e2 field packet ≥0.3C · e3a heavy
blob seals · e3b light tilted blob translates · e4 curvature linearity ·
e5 CHSH · e6 rung tongue in the foam · e7 tuning curve (P4) · e8 the comma
is paid, monotone · e9 the fifth (3:2) lives · d1 double-slit fringes ·
t1 Tonomura clicks · q2 eraser ± · t4 HOM ordering · qt_lo sub-threshold
(nothing condenses) · qt_hi above-threshold (condensation) · LIN: every
fired grain across every log sits on the eps(w) = A0·w/2π grid (ħ-linearity
as a battery invariant).

## Variant verdict (2026-07-27, S1)

| variant | law | result |
|---|---|---|
| V1 | source atom, per-cycle floor | 14/17 — floor kills slow channels (e4, e8 roughness dead) |
| **V2** | **source atom, credit (lapse at 2)** | **16/17 — adopted** |
| V3 | destination atom, floor | 14/17 — same floor failures |
| V4 | destination atom, credit | 15/17 — big dense atoms over-tax pair exchange (e6 collapses) |

The one V2 failure is e7: tuned pairs unpin from the rung by occupancy
drift (frequency-correction gap, ROADMAP open area A1 — pre-named, now
load-bearing; the uniform route to it is the S2 dense-amplitude rewrite).

Law iteration history: L1 = the §6.4 proposal (Γ_m 0.02); L2 = Γ_m → 0.10
(under the unified pitch landscape w2=2.9, q=1.2, per-link detunes are ~5×
the old landscape's; 0.02 froze the blob rim and starved tuned pairs).
