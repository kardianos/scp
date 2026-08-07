# v93 L1 — first findings (face 1 & 2 of the unitary dense hop)

Opened 2026-08-07. The unitary dense hop (pass U) is implemented behind
`amp_nat` in both kernels, byte-inert at `amp_nat=0`, C≡Go (battery ALL
GREEN 87). This file records what `amp_nat>0` actually does to the dense
sector, measured against the e3b (the L1-A acceptance).

## L1-B (conservation) — GREEN, the theorem holds

Bath T=40 at the V3a table, `amp_nat` sweep:

| amp_nat | drift_rel |
|--------:|----------:|
| 0       | 0.000e+00 |
| 3       | 0.000e+00 |

Conservation is a theorem of the update (pairwise Givens norm-preservation
+ the unchanged door): Σ(Em+Ee+Es) is invariant at machine precision, with
no ledger to patch. This is the design's central claim (II.4), confirmed.

## e3b baseline reproduction — byte-exact

`exp=blob bath=1 L=16 T=80 amp=0.5 sigma=2.5 kx=1.1 wf_on=1` (the v92
documented e3b config), `amp_nat=0`, 3 seeds:

| seed     | speed    | cos_to_kdir |
|----------|----------|-------------|
| 111      | 0.003138 | -0.9394     |
| 20260802 | 0.001544 | -0.2798     |
| 314159   | 0.001556 | -0.9286     |

Byte-identical to the v92 `L1_FINDINGS` addendum-2 baseline. The unitary
lane is genuinely inert at `amp_nat=0`.

## Face 1 (tau_s = amp_nat·base·gsym·sqrt(head_i·head_j)) — FAILED

The blob loads cells to ~cap=2.5 ⇒ head≈0 in the dense core ⇒ tau_s≈0
there ⇒ the coherent tilt (which lives in the core) is FROZEN. Only the
random-phase surface transports. e3b cos went seed-variant/incoherent
(-0.95/+0.88/+0.28 at amp_nat=0.25). **The head factor is wrong**: it
suppresses exactly the region whose phase gradient should drive the
current. (Anticipated by the design: head was not mandated, only the
closure gate.)

## Face 2 (tau_s = amp_nat·base·gsym, head dropped) — ENGAGES, partial

The closure gate (cos^p) survives as the angle envelope; the door (pass 6)
enforces cap. 3-seed e3b:

| amp_nat | seed 111 cos | seed 20260802 cos | seed 314159 cos |
|--------:|-------------|-------------------|-----------------|
| 2       | +0.54        | +0.81             | +0.93           |
| 3       | -0.87        | -0.87             | -0.45           |
| 5       | +0.79        | -0.83             | (mixed)         |

speed (all ~0.0016–0.0046, scaling roughly linearly with amp_nat).

**The positive result:** within each amp_nat the SIGN of cos is seed-
robust (all + at 2, all − at 3) — a coherent directional bias exists,
driven by the phase-tilt. The current J = 2·tau·Im(ψi*·ψj) (the cross
term the additive ledger rejected) IS doing something: this is the first
time dense transport has shown a seed-robust direction in this programme.

**The negative result:** |cos| is 0.5–0.9 (not →1) and the sign REVERSES
between amp_nat=2 (+) and 3 (−). A centroid trace (amp_nat=2.5, diag
every 0.5 t.u.) shows com_dev grows only 0.017→0.12 over 42 t.u. with
clear wobble (0.06→0.10→0.07), and Em_tag slowly melts (81.7→78.8, 3.6%).
So the blob translates SLOWLY with a WOBBLY direction — not the clean
ballistic drift L1-A wants. bath=0 gives exactly zero motion (the medium
is required).

## Diagnosis

The wobble + amp_nat-dependent sign reversal is the signature of a free
tight-binding packet: the unitary hop is purely Hamiltonian (reversible),
so the coherent current sloshes rather than damping into a steady drift,
and the localized packet diffracts/spreads into the bath. The additive
model keeps the blob bound (cap + inflight + headroom) and drifts it
slowly; the bare unitary hop exchanges amplitude freely and the packet
spreads. Two suspected Trotter/error sources to test next:

1. **The dense clock is applied OUT of pass U** (pass 6 advances th2 by
   w2e·dt AFTER the hop), whereas pass F applies the local precession
   BEFORE the hops, in-pass. The hop-then-precess vs precess-then-hop
   split is a first-order Trotter error that accumulates over 4000 steps
   and could drive the wobble. Face 3 moves the precession into pass U.
2. The packet needs a binding mechanism to stay localized while its
   coherent current translates it (the design's "rotations exchange
   reversibly" does not, by itself, bind).

## Status

- L1-B (conservation): **GREEN** (the theorem).
- L1-A (translation): **OPEN** — coherent direction exists (seed-robust
  within amp_nat, a first) but |cos|≠1 and the blob wobbles/melts rather
  than translating ballistically. Not yet the speed≥2.6e-3 cos→1 bar.
- L1-C (anti-ignition): bath geometry at amp_nat=3 (phi_nom 1.7044,
  z_live 16.89) close to amp_nat=0 — preliminary benign; to be measured
  properly.
- Invariant surface (87 V3a bars): **held** at amp_nat=0.

Face 3 (in-pass precession, mirroring pass F exactly) is next.
