# v93 — FACE A: the symmetrized hop schedule (hop_order)

2026-08-08. RESUME §7-A / RING_DNLS §A.4. The sequential i<j Givens sweep of
pass U is a deterministic C6-breaking Trotter asymmetry — the named source of
BOTH the vacuum retention ceiling (frozen scaffold slips ~t=50-100 with all
noise off) AND item 1's dt-invariance failure (same family). This face adds a
byte-inert knob `hop_order` selecting the pass-U schedule, in both kernels.

## The knob (byte-inert at 0; battery ALL GREEN 87)

- `hop_order=0` (default) — sequential canonical sweep (the prior behaviour,
  byte-faithful: verified against the pre-edit kernel at amp_nat=2).
- `hop_order=1` — **Strang symmetric**: forward tau/2 + reverse tau/2. Time-
  reversal-invariant; cancels the chiral sweep bias.
- `hop_order=2` — **randomized**: per-step Fisher-Yates shuffle of the
  canonical slots, full angle (separate xorshift stream; physics draws
  untouched).

Both kernels mirror exactly (C `apply_uhop` / Go `applyUhop`; collect active
canonical slots into a flat list, then dispatch by schedule).

## A. The ring retention ceiling — sweep order FIXED; condensation is separate

Vacuum C6 ring, seed_mw=2, amp_nat=2 amp_logate=1, T=1000. W = cycle winding;
R2 = seeded-m coherence; Em_min = the |psi|>0-on-cycle condition.

| scaffold | hop_order | W@1000 | R2@1000 | hold-until | note |
|---|---|---|---|---|---|
| **live geom** | 0 seq | +1.0 | 0.63 | degraded | sweep chirality breaks C6 |
| **live geom** | 1 strang | **+2.0** | **0.978** | **>=1000** | HELD |
| **live geom** | 2 rand | **+2.0** | **0.993** | **>=1000** | HELD |
| frozen | 0 seq | 0.0 | 0.45 | ~t=20 | phase wanders; Em holds |
| frozen | 1 strang | (+1 walk) | — | **~t=108** | clean W=+2 R2=1.0 to 108 |
| frozen | 2 rand | 0.0 | 0.72 | ~t=16 | stochastic noise seeds slip |

1. **The sweep-order artifact is real and Strang/random fix it** — in LIVE
   geometry (the physical case) sequential degrades to W=+1 while both
   symmetrized schedules hold W=+2 (R2 0.978/0.993) to t=1000. The frozen-
   scaffold Strang arm holds a clean W=+2 R2=1.000 to t=108 (5x sequential),
   proving the chirality is gone.
2. **The frozen t=1000 bar is NOT met — and the reason is a SEPARATE
   mechanism.** Frozen Strang slips at ~t=108 when Em_min->0 (a voice drains,
   detunes, decouples — the q_detune rich-get-richer condensation of route B,
   not sweep order). Live geometry's bond feedback re-equalizes edges and
   suppresses this, so live holds. The frozen ceiling is now OWNED BY
   CONDENSATION (face C), no longer by sweep order.
3. **Random is worse than Strang for retention** — the per-step shuffle
   injects stochastic phase noise that seeds dephasing faster (frozen slip
   ~t=16 < Strang ~t=108). Random trades the systematic chirality for noise.

## B. e3b dt-invariance — both symmetrized schedules improve it

e3b blob, amp_nat=2 FIXED, vary dt; fd = dense-hop current (the +x transport).

| hop_order | fd @dt=0.02 | fd @dt=0.01 | spread |
|---|---|---|---|
| 0 seq | 299.9 | 287.8 | **4.0%** |
| 1 strang | 285.5 | 289.1 | **1.3%** |
| 2 rand | 308.6 | 309.5 | **0.3%** |

Sequential's 4% dt-spread (item 1's Trotter artifact) drops to 1.3% (Strang)
and 0.3% (random). **Random is the most dt-invariant for transport.**

## Synthesis — complementary schedules

| | retention (ring) | transport dt-invariance (e3b) |
|---|---|---|
| sequential | worst (chirality breaks C6) | worst (4%) |
| **Strang** | **best** (deterministic, no noise; live holds to 1000) | good (1.3%) |
| **random** | poor (noise seeds dephasing) | **best** (0.3%) |

- The sweep-order Trotter artifact is CONFIRMED and ADDRESSED: the live ring
  holds W=+2 to t=1000 under Strang/random where sequential degrades; the e3b
  dt-spread falls 3-13x.
- **Strang (hop_order=1) is the recommended default for the retention face**
  (best retention, good transport, deterministic). Random is the transport-
  dt-invariance champion but a retention liability.
- The frozen-scaffold t=1000 bar reveals that condensation (route B) is the
  remaining frozen-ceiling owner — the natural bridge to face C.

## Files
- Kernel: `kernel/freecell.c` (struct `hop_order`; `apply_uhop`; collect +
  dispatch in pass U; `actslot_` scratch; `hop_rng_s`). Go mirror: `fab/step.go`
  (`applyUhop` + dispatch), `fab/params.go` (HopOrder), `fab/sim.go` (actslot,
  hopRng), `fab/run.go` (seed + header).
- Script: `runs/faceA.sh` -> `runs/faceA/{ringA_*,e3b_*}`.
