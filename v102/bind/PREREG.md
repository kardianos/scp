# v102 bind — pre-registration (the two walls before K3; does a mechanism force a size)

**Written:** 2026-10-09, before the registered runs (M2, M3). The M1
instrument was run while it was being built (§3, M1), so M1 is
reported, not gated.
**Authority:** user, 2026-10-09: "For K3 stops, for both questions,
These are beyond me for the most part. If pairing has energy that
needs to be parked, that is a good sign because that means field is
gathering (by definition) but the mechanism for stability is
intrinsically not something I have an answer to. Try different
mechanisms, see what works; ideally it would be something that in a
given configuration, it would also force a specific size."
**Status:** a measurement in models. It rules nothing.
- No kernel edit and no law constant change in this step.
- Models: the leading-order lock (`../perm/cas/lo.py`), TeNPy chains
  (`../perm/cas/tp_chain.py`), exact rings, and the microscopic S-law
  (`../perm/cas/graphm.py`).

---

## 1. The two walls (`../kernel/K2b_DECIDE.md` §5)

1. **The (1,1) link.** Two N = 1 cells side by side (two adjacent
   holes, or a hole next to an end cell) share a link whose conversion
   is close to resonant: Δ₁₁ = Ω − 2E_f(1) = −0.0243, while the pair's
   coupling is g|P| = 0.068 (3+2), 0.096 (W). The leading-order lock
   gives that link no term. The kernel's lock does the same. Under
   bright light the kernel makes such pairs.
2. **The sink.** Two holes bind in the all-singlet channel by
   B = −0.50 J0 (W), −0.13 J0 (3+2). A closed lock cannot shed B (its
   own quantum costs the gap, 1.32 / 0.44 J0). So pairs never form
   from holes that meet, unless the energy is parked somewhere.

## 2. The options

### Wall 1: pricing the (1,1) link (`lo.P11`)

The link's exact two-level problem is the singlet pair against one
space quantum: E± = Δ₁₁/2 ± √(Δ₁₁²/4 + g²|P|²). The non-singlet
channels do not couple to the link at all.

| | option | singlet pair on the link | other channels | status |
|---|---|---|---|---|
| **A0** | `zero` | no term | no term | today's lock (lo.py, `lockq.c`) |
| **A1** | `core` | excluded (penalty 1.0 ≈ 104 J0) | excluded | a new rule: "two N = 1 cells may not share a link" |
| **A1s** | `score` | excluded | no term | V taken strictly: both branches hold about half a space quantum (below), so V admits no singlet (1,1) state; the other channels hold none |
| **A2** | `lower` | E₋: −8.45 J0 (3+2), −11.35 J0 (W); 59% / 56% space on the link | no term | the space branch. **Violates V** (real link occupancy). It is the lowest exact state. Contrast only |
| **A3** | `upper` | E₊: +5.91 J0 (3+2), +8.81 J0 (W); 41% / 44% space | no term | the field branch: the one continuous with "no space quantum" as g → 0. No new constant |

A2 and A3 apply the branch energy as (E/|P|²)·A†A on the (1,1) block.
That operator is exact on the isolated link.

### Wall 2: where the binding is parked

| | parking place | how it is measured |
|---|---|---|
| **P0** | nowhere (closed lock) | `../kernel/K1_DECIDE.md` K1-2: the bound weight stays at its value at formation |
| **P1** | **the formation space quantum.** The hole that completes a bound object forms directly into it. The space quantum leaving through the boundary carries the binding: the band window, ±2s_k ≈ ±12.5 J0, covers any binding here | M3a: the yield per completing event. It is already in K2b's event rule (the coherent jump) |
| **P2** | light: inelastic scattering of cold light at the boundary | carried from K1-3 (≈ 3.4×10¹⁰ photons per pair). That run started from an S = 1 state because the six-ring has no S = 0 unbound state below the gap. A re-run on a larger ring is too large for that S-matrix method; this is stated, not re-run |
| **P3** | **a third hole**, which takes the binding as kinetic energy | M3b: energetically open only if the binding is smaller than the band a free hole can absorb (the one-hole band width) |
| **P4** | **a lock quantum** (a magnon walking off along the fiber) | M3c: energetically open only if a cluster's marginal binding exceeds the lock's gap |

### Size

"Size" here means the number of holes N* in the largest bound cluster.
Two ways it could be forced:

- **(i) Energetic saturation:** clusters up to N* are bound and compact,
  and the (N* + 1)-th hole is free.
- **(ii) Dynamic saturation:** the largest cluster that the open sinks
  can actually build (P3 or P4 open at every step up to N*, and closed
  for the next).

N* is read, never placed. It is not mapped onto any nucleus.

## 3. Measurements

**M1 (reported, not gated): the (1,1) link in the exact S-law.**
`cas/m1_link.py`, run as `cas/run_m3.sh`'s first two lines.

- Pinned holes on a W n = 4 ring (cmax = 2) and a 3+2 n = 8 ring
  (cmax = 3).
- The exact d = 1 spectrum, for each level:
  - its energy as a binding B;
  - its space on the shared link;
  - its overlap with each lo pricing's ground state placed in the
    space-free configurations.
- **Disclosure:** the instrument was run on both rings while being
  written (eigsh with 60 levels for the 3+2 ring; the registered run is
  dense). Seen then:
  - the exact ground is the space branch;
  - a field-like level with about 10% link space sits above it, at
    +0.92 J0 (W) and +2.18 J0 (3+2);
  - zero, upper and score place their ground states 0.05–0.9 J0 from
    it;
  - lower sits within 15% of the exact ground.

  So M1 cannot be gated. It is reported as the check on what each
  pricing represents.

**M2: cluster energetics** (`cas/run_m2.sh`; TeNPy, `tp_chain.py`
extended to K mobile holes and to `p11`).

| | W | 3+2 |
|---|---|---|
| lock | open chain, R6 end cells, n_bulk = 18 and 26 mobile sites | n_bulk = 26 and 42 |
| transport | t = 0.1 J0 (R9) | t = 0.1 J0 |
| holes | N_h = 0 … 5 at both lengths; 6 at 18 | N_h = 0 … 6 |
| pricing | zero, core, lower, upper at both lengths; score at 18 | the same; score at 26 |
| χ | 128 | 64 |
| sector | minimal: (0, 0) for even N_h, (1, 1) for odd | 0 or 1 |
| multiplet (short length, zero and upper) | (2, 0), (0, 2) even; (3, 1), (1, 3) odd | 2 even; 3 odd |

Readouts, in J0:

- **B(N)** = E(N) − E(0) − N·[E(1) − E(0)]: the total binding against
  N free holes.
- **b(N)** = B(N) − B(N − 1): the binding of the N-th hole to the
  cluster of N − 1.
- **⟨d⟩:** the mean hole–hole distance over all pairs, from the
  pair-distance distribution Pr(d).
- **Adjacency:** Pr(1), the number of adjacent hole pairs.
- **The cluster's (S, I) or S,** from the multiplet sectors.
- The truncation error is reported.

**M3: the sinks.**

- **M3a (P1).** `cas/m3_yield.py`, exact rings: 3+2 n = 10, W n = 6;
  all five pricings; t = 0.1 J0.
  - The first jump (mode 0) is at cell 0. The second jump, of opposite
    charges, is at d = 1 … n/2.
  - Reported: the bound weight of the jump state (levels with
    B < −0.05 J0), the mean B, and the event's δ against the window.
- **M3b (P3), energetic.** Compare each bound b(N) from M2 with the
  one-hole band width W_1. Take W_1 from `../perm/` X3 (rings): 3+2
  n = 10 is 0.112 J0, W n = 6 is 0.388 J0. P3 is open for a step if
  |b(N)| < W_1.
- **M3c (P4), energetic.** Compare each bound b(N) with the lock's gap:
  W 1.32 J0, 3+2 0.44 J0 (`../perm/` and `../fabric/`, extrapolated).
  P4 is open for a step if |b(N)| > gap.
- **Dynamics for P3 and P4** are registered here only as follow-ups. A
  step that is energetically open gets its own closed-lock run in the
  next page.

## 4. Bars

**Binding of the N-th hole** (both lengths):

| | W | 3+2 |
|---|---|---|
| bound | b(N) ≤ −0.10 J0 | b(N) ≤ −0.05 J0 |
| free | b(N) ≥ −0.03 J0 at the longer length | b(N) ≥ −0.02 J0 at the longer length |

Anything between is "marginal". The W bar is twice the worst
unconverged error that `../attract/` saw at χ = 128.

**Compact cluster:** ⟨d⟩ changes by less than 15% between the two
lengths. A free hole makes ⟨d⟩ grow with the length.

**Size forced (i):** for some N* ≥ 2, at both lengths:

- every b(2 … N*) is bound;
- the N*-hole cluster is compact;
- b(N* + 1) is free, or the (N* + 1)-hole ⟨d⟩ grows by 15% or more.

**P1 usable:** a bound weight of at least 0.01 at some d whose event is
inside the window.

**Per-cell scoring.** A cell is (pricing × parking place). Each cell is
scored on six criteria:

| | criterion |
|---|---|
| C1 | a cluster with N ≥ 2 is bound |
| C2 | the parking place is open and usable for it (P1 usable; P3 or P4 energetically open; P2 carried as about 10⁻¹⁰ per photon, so not usable) |
| C3 | the cluster is compact |
| C4 | V is respected: A2 fails; A3 is reported with its space fraction times its adjacency |
| C5 | no collapse: the floor holds, and no adjacent block forms (Pr(1) per pair below 0.5) |
| C6 | a size is forced (i or ii) |

| tier | criteria |
|---|---|
| Works, with size | C1–C6 |
| Works | C1–C5 |
| Partial | anything less |

## 5. Endings (exactly one is claimed in `DECIDE.md`)

| ending | condition |
|---|---|
| **S** | some field-reading cell (A0, A1, A1s, A3) reaches "Works, with size" |
| **N** | some field-reading cell reaches "Works", but none forces a size |
| **X** | no field-reading cell reaches "Works" |

A2's cells are reported beside the endings, never as one of them.

## 6. What comes after

If ending S or N, the next step (its own pre-registration) is the
kernel:

- `lockq_p11` in `lockq.c`, off by default and byte-inert;
- the full battery;
- A/B against `lo.py` under the chosen pricing;
- lock-in-light runs to see what a given configuration builds.

That step needs no new law constant. A1 (`core`) would be a new rule,
and needs the user.
