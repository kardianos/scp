# The unification battery — ROADMAP §6.5, operational

One **law table** (`laws_*.cfg`), byte-shared by every experiment; one
**apparatus file** per experiment (`apparatus/*.cfg`) that may not touch a
law key (the runner refuses to start otherwise); **physics acceptance**,
not byte comparison. Law variants compete as whole tables — switching a
law per experiment is structurally impossible here.

## THE RATCHET RULE (standing, 2026-07-28)

1. **Every modification of the kernel or of a law table runs the FULL
   battery before it is committed.** No exceptions for "obviously inert"
   changes — inert is what the green run proves. (S2 exists because a
   zero-ablation was skipped once.)
2. **Experiments that pass and defend a claim are added to the battery**
   and become part of the gate for every future modification. The suite
   only grows.
3. **Bars encode defended claims, never snapshots of current output.**
   A bar may be sharpened by a measurement; it may not be softened to
   make a change pass.
4. **A green test leaves the suite only by explicit user decision** —
   never as a side effect of making something else work.

```
python3 battery.py --laws laws_V2z.cfg [--jobs 10] [--only e6_pairs ...] [--skip-run]
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

## Variant verdict (2026-07-28, after §6.10 S2 retirement)

| variant | law | result |
|---|---|---|
| V1 | source atom, per-cycle floor | 15/17 — the floor cannot pay the comma (e8) or fit e4 |
| V2 | source atom, credit, kappa_freq=0.6 | 17/17 — the §6.9 table; the bias buys margin, not passes |
| **V2z** | **source atom, credit, NO bias (kf=0, kr=0)** | **17/17 — the standing law table (one fewer free constant)** |
| V2d | source atom, credit, kappa_reac=1 (derived reactive term) | 16/17 — e3b decoheres: rate transport cannot host the raw current; kept as the S2 exhibit and S2-full acceptance criterion |
| V3 | destination atom, floor | 16/17 — e8 (comma needs memory) |
| V4 | destination atom, credit | 16/17 — e3b (dense-atom tax perturbs the front) |

Annealing history (each value applied to ALL variants — mechanisms are
variant-orthogonal; full record in ROADMAP §6.9/§6.10):
L1 = the §6.4 proposal (Γ_m 0.02); L2 = Γ_m → 0.10 (unified landscape
has ~5× per-link detunes); L3 = kappa_freq dispersive bias (choir's
correction, A1); L4 = flight-loads-pitch (zero-constant kernel law —
the shelter was detuning its own pair); L5 = bias windowed by the
acceptance (blob rims feel nothing); L6 = bias × mutual gate closure
+ mob_floor 0.004 (pair bleed vs front recruitment both pass);
L7 = S2 derivation (CONSONANCE Part VII): the bias is real amplitude
physics, derived in v89/s2/ — and unnecessary at rate level after L4;
kappa_freq retired (the zero ablation had never been re-run post-L4;
the battery caught the vestige).
