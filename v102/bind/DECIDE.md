# v102 bind — decision (pricing the (1,1) link; parking the binding; size)

**Written:** 2026-10-09, after the registered runs (M2, M3) in `cas/out/`.
Gates: `PREREG.md`. Analysis: `cas/analyze_m2.py`
(`cas/out/analyze_m2.out`).
**Status:** a measurement in models. It rules nothing. The user ruled
R12 (singlet core, = option A1s) during this step, on the measurements
reported here. The kernel steps that followed (`../kernel/K2d_*`,
`../kernel/K2e_*`) are summarized in §6.

## Ending: **S**

- Three field-reading pricings reach "Works, with size" with the
  formation space quantum (P1) as the parking place, in both fabrics:
  zero (A0), core (A1) and upper (A3).
  - The cluster size is **N\* = 2**: pairs.
  - Pairs are bound and compact at both lengths. The third hole is free;
    a fourth hole makes a second pair.
- A1s (score, now R12) gives the same numbers at the registered
  length (the pre-registration ran score at one length only). Its second
  length, run as a labelled supplement after R12 was ruled, also gives
  N\* = 2 in both fabrics (§Supplementary).
- A2 (lower, the space branch) clumps and violates V. It is reported
  beside the ending (§4).

## 1. M1: the (1,1) link in the exact S-law (reported, not gated)

Two pinned holes at d = 1 (`cas/out/m1_*.out`). B is in J0; "link"
is the space on the shared link.

| | W ring n = 4, sector (0, 0) | 3+2 ring n = 8, sector 0 |
|---|---|---|
| exact ground | B = −4.91, link 0.48 (the space branch) | B = −3.31, link 0.47 |
| exact field-like level | B = +0.92, link 0.08 | B = +2.18, link 0.11 |
| zero: its ground's largest overlap | 0.51 with the field-like level (0.18 with the ground) | 0.37 field-like (0.21 ground) |
| score | 0.65 field-like (0.05 ground) | 0.59 field-like (0.03 ground) |
| upper | 0.60 field-like (0.10 ground) | 0.53 field-like (0.10 ground) |
| lower | 0.43 with the exact ground | 0.39 with the exact ground |
| d = 2 | exact ground B = −0.19; every lo pricing 0.00 | exact B = −0.05 |

- **The exact (1,1) ground is the space branch,** with about half a
  space quantum on the shared link. V excludes it. Lower represents it.
- **The field readings represent the field-like level.** It sits 0.9 J0
  (W) and 2.2 J0 (3+2) above two free holes. Score represents it best:
  overlap 0.59–0.65.
- As the pre-registration disclosed, M1 was seen while the instrument
  was built. It is the check on what each pricing means, not a gate.

## 2. M2: cluster energetics (TeNPy chains, t = 0.1 J0)

The binding of the N-th hole, b(N) = B(N) − B(N − 1), in J0, at the two
lengths. ⟨d⟩ is the mean hole–hole distance.

**3+2** (n_bulk = 26 / 42; bars: bound ≤ −0.05, free ≥ −0.02):

| pricing | b(2) | b(3) | b(4) | b(5) | b(6) | ⟨d⟩ of the pair |
|---|---|---|---|---|---|---|
| zero | −0.130 / −0.129 | −0.016 / +0.000 | −0.114 / −0.128 | −0.034 / −0.017 | −0.103 / −0.115 | 2.99 / 2.99 |
| core | −0.129 / −0.129 | +0.003 / −0.000 | −0.131 / −0.127 | −0.013 / −0.015 | −0.117 / −0.005 | 2.99 / 2.99 |
| score † | −0.129 / −0.129 | −0.016 / −0.001 | −0.113 / −0.127 | −0.011 / −0.016 | −0.117 / −0.002 | 2.99 / 2.99 |
| upper | −0.130 / −0.129 | −0.016 / +0.002 | −0.116 / −0.129 | −0.032 / −0.016 | −0.103 / −0.125 | 2.99 / 2.99 |
| lower | +4.249 / +4.248 | −0.051 / −0.000 | +1.375 / +1.324 | +0.536 / +0.541 | +0.787 / +3.708 | 1.00 / 1.00 |

**W** (n_bulk = 18 / 26; bars: bound ≤ −0.10, free ≥ −0.03):

| pricing | b(2) | b(3) | b(4) | b(5) | ⟨d⟩ of the pair |
|---|---|---|---|---|---|
| zero | −0.504 / −0.503 | −0.029 / +0.005 | −0.471 / −0.508 | −0.025 / −0.028 | 3.06 / 3.06 |
| core | −0.500 / −0.499 | −0.026 / +0.001 | −0.468 / −0.500 | −0.020 / −0.028 | 3.08 / 3.08 |
| score † | −0.503 / −0.501 | −0.028 / +0.005 | −0.470 / −0.506 | −0.023 / −0.028 | 3.07 / 3.07 |
| upper | −0.502 / −0.501 | −0.029 / −0.006 | −0.470 / −0.495 | −0.024 / −0.028 | 3.07 / 3.07 |
| lower | +3.413 / +3.412 | −0.042 / −0.000 | +1.182 / +1.402 | −0.041 / −0.066 | 1.01 / 1.01 |

† Score's second length (3+2 n_bulk = 42, W n_bulk = 26) is supplementary:
run after R12 was ruled, not pre-registered (§Supplementary).

- **Pairs, and only pairs.** For every field reading, in both fabrics:
  - b(2) is bound at both lengths, and the pair's ⟨d⟩ is the same at
    both (3.0 in 3+2, 3.1 in W);
  - b(3) is free at the longer length, where the third hole's ⟨d⟩
    grows by 50–190%;
  - b(4) ≈ b(2): the fourth hole makes a second pair, and the pattern
    repeats.
- **Size forced (i): N\* = 2** for zero, core and upper, in both fabrics.
  The (N\* + 1)-th hole is free.
- **No adjacency.** Pr(1) = 0.00 per pair in every field-reading ground
  state, so the pricing of the (1,1) link does not enter the ground
  states at all. That is why zero, core, score and upper agree to about
  0.01 J0.
- **The pair is the all-singlet channel.** In 3+2 the two-hole S = 1
  sector sits +0.133 J0 above S = 0. In W, (2,0) and (0,2) sit +0.531
  J0 above (0,0). The same ordering holds at N = 3–6 in 3+2 (+0.10 to
  +0.16 J0) and at N = 3–5 in W (+0.53 to +0.66 J0), under zero and
  upper (`analyze_m2.out`, "multiplets").
- Truncation: ≤ 4e-5 (3+2), ≤ 8e-4 (W). The W bar (−0.10) is twice the
  worst unconverged error seen in `../attract/`.

## 3. M3: the parking places

**M3a: P1, the formation space quantum** (`cas/out/m3a_*`, exact rings:
3+2 n = 10, W n = 6). The first hole is at cell 0. The second, of
opposite charges, is at distance d. The table gives the jump state's
weight on bound levels (B < −0.05 J0), and whether the event is inside
the window:

| d | 3+2 (zero, score, upper, core) | W (zero, score, upper, core) |
|---|---|---|
| 1 | 0.000; out of the window for core and score (mean-window rule) | 0.000–0.0003; out for core and score |
| 2 | 0.001 | 0.012 |
| **3** | **0.662** | **0.517–0.519** |
| 4 | 0.001 | — (n = 6 ring) |
| 5 | 0.631 | — |

- **P1 is usable:** at d = 3 the second hole is born into the pair with
  weight 0.66 (3+2) and 0.52 (W), inside the window. Every field
  reading agrees, because the pair never touches a (1,1) link.
- The formation space quantum carries away the binding, which is far
  inside the ±12.5 J0 window.
- Lower: d = 1 is born bound (0.86 in 3+2, 0.81 in W) into the space
  branch, which V excludes.

**M3b: P3 (a third hole) and M3c: P4 (a lock quantum)**, for the
bound step b(2):

| | \|b(2)\| | one-hole band W1 | gap | P3 open (\|b\| < W1) | P4 open (\|b\| > gap) |
|---|---:|---:|---:|---|---|
| 3+2 | 0.13 | 0.112 | 0.44 | no | no |
| W | 0.50 | 0.388 | 1.32 | no | no |

- Neither sink is open for pairs, in either fabric. P3 opens only for
  zero and upper's b(6) in 3+2 (|b(6)| = 0.103 at n_bulk = 26, against
  0.112): the third pair's step, past N\*.
- **P2 (cold light)** is carried from K1-3 at about 10¹⁰ photons per
  pair, so it is not usable.

## 4. Scoring (PREREG §4)

| cell | C1 bound | C2 parking | C3 compact | C4 V | C5 no collapse | C6 size | tier |
|---|---|---|---|---|---|---|---|
| A0 zero × P1 | ✓ | ✓ (0.66 / 0.52) | ✓ | ✓ | ✓ | ✓ N\* = 2 | **Works, with size** |
| A1 core × P1 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ N\* = 2 | **Works, with size** |
| A3 upper × P1 | ✓ | ✓ | ✓ | ✓ (space fraction × adjacency = 0) | ✓ | ✓ N\* = 2 | **Works, with size** |
| A1s score × P1 | ✓ | ✓ | ✓ (second length †) | ✓ | ✓ | ✓ N\* = 2 (second length †) | **Works, with size** (†) |
| any × P0, P2, P3, P4 | ✓ | ✗ (closed or ~10¹⁰ photons) | | | | | Partial |
| A2 lower × P1 (reported beside) | ✗ (b(2) = +4.2 / +3.4 in the chains) | (✓ at d = 1) | ✗ | ✗ | ✗ (adjacent) | ✗ | — |

- A2 note: in the exact rings, lower binds adjacent pairs by about
  −3.7 J0 (3+2) and −3.2 J0 (W) (M3a's bound levels, ⟨d⟩ = 1.0). The
  TeNPy chains under lower are trapped (B(2) = +4.2 J0 at d = 1, with
  holes caught at the ends). Either way it clumps, and it violates V.
  It is a contrast only.

## 5. What this does and does not show

- **Energetically, the lock forces a size: pairs.** It does so under
  every reading that respects V, with the binding parked in the space
  quantum that leaves at formation.
- **It does not show that light builds pairs.** P1 needs the second
  hole to form at d ≈ 3, in the complementary mode, while nothing else
  is nearby. Only the kernel steps can say whether light does that.
- **The (1,1) pricing matters only off the ground state.** None of the
  field-reading ground states holds an adjacent pair. The pricing
  decides what *formation* can make (§6).

## 6. What came after (summary)

- **K2d** (`../kernel/K2d_DECIDE.md`): every pricing became a kernel
  option. Under K2a's mean-energy window, core and score kept holes
  apart in light. The 3+2 lock jammed at three holes, with
  complementary pairs at d = 3 born bound (0.51 under score).
- **R12** (user): singlet core, A1s.
- **K2e** (`../kernel/K2e_DECIDE.md`): R12 made exact in the formation
  event (projection). Bright light then fills the lock with adjacent
  holes in the allowed channels, far from bound (B = +15 to +17 J0 in
  3+2). K2d's jam was the mean window's doing, not R12's.
- **K3** (`../kernel/K3_PREREG.md`): a lock in light under R12, over
  long times.

## Supplementary (labelled; run after R12 was ruled)

- **What was added.** The pre-registration ran score at one length
  (3+2 n_bulk = 26, W n_bulk = 18). After R12 made score the ruled
  pricing, its second length was run with the same instrument, bond
  dimensions and bars (`cas/run_supp.sh`): 3+2 n_bulk = 42 (χ = 64),
  W n_bulk = 26 (χ = 128). Every entry marked † in §2 and §4 comes
  from these runs. They gate nothing; the ending (S) stands on zero,
  core and upper alone.
- **Result.** Score behaves as zero, core and upper do:
  - 3+2: b(2) = −0.129 at both lengths, the pair's ⟨d⟩ = 2.99 at both;
    b(3) = −0.001 (free; the third hole's ⟨d⟩ grows by 96%);
    b(4) = −0.127 (a second pair).
  - W: b(2) = −0.501, ⟨d⟩ = 3.07 at both lengths; b(3) = +0.005 (free;
    ⟨d⟩ grows by 63%); b(4) = −0.506.
  - **N\* = 2, size forced (i)**, in both fabrics. P3 and P4 are closed
    at every N. Pr(1) = 0.00 throughout.
  - Truncation: ≤ 2e-5 (3+2), ≤ 6e-4 (W), inside the §2 figures.
- **b(6) in 3+2.** Score, like core, gives b(6) bound at n_bulk = 26
  (−0.117) and free at 42 (−0.002); zero and upper give it bound at
  both. Whether the 42-site value is the ground state or a DMRG trap was
  not checked. It is past N\* and changes nothing above; it is recorded,
  not explained.
