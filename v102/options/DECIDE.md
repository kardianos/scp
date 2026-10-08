# v102 options — decision (where the lock fails, option by option)

**Written:** 2026-10-08, after the outputs in `cas/out/`. Options,
experiments and the scoring rule: `PREREG.md`. Deviations from it are
listed in §9.
**Status:** this page is a measurement. It rules nothing. No kernel and
no law table was edited. The battery runs with changed constants were
made from scratch builds through the battery's `-extra` flag; they are
not adoptions.

---

## 1. The answer

**One bundle makes the most sense:** **B-label-cool**.

| question | choice | the measured reason |
|---|---|---|
| ceiling at the boundary | **C-lock**: a lock cell's cap counts its lock links | both label-free ceilings fail (§2) |
| transport | **T0**: no hop on a link touching a lock cell | T4 equals T0 alone, but collapses under any space flow (§3) |
| boundary conversion | **X1**: S3 on every link whose two cells hold field | the only cooling channel (§7) |
| space flow on lock links | **none**, by the same label; coherent flow along the fiber is tolerated to `2t_s < 0.2275` | rate-law flow kills the lock and creates energy (§3) |
| door | **closed on lock cells**: the work function already does this where a lock carries no space deficit | D2 costs 10 bars and makes detection weak and line-selective (§4) |
| g | **unset**; the window is `g/Δ ≤ 0.1` | §5 |
| size | **not selected** | §6 |

This bundle needs **exactly one label**: which cells belong to a lock.
No local rule replaced it in any test.

- The label-free ceilings either leak the lock (C-geo) or freeze free
  light (C-held).
- Given the label, every other choice is local or free. T0, the space
  freeze on lock links and the door exemption all reuse the same label.

What still fails under the best bundle:

1. **Formation needs light colder than the lock's gap**, and about 10⁸
   photons per lock at g/Δ = 0.1. Resonant or broadband light heats the
   lock. With space flow at the boundary it also destroys lock quanta.
2. **The lock is never the energy minimum against free light.** It
   costs 1.325 per quantum, against −0.34 (1-D) or −3.94 (cubic) for
   free light at its band bottom. It lives by the label's refusal to
   move it, and the same refusal means transport cannot fill a lock.
3. **Touching locks merge.** With X1, a cube beats two separate 4-rings
   by 6.8 J. Keeping locks apart would need a second label: lock
   identity, not just membership.
4. **g is unset,** and the window that keeps the derived line is narrow.

## 2. Q1 — the boundary

### E1. Transport inside the ring

n = 6, caps 2, Σc ≤ 3, S^z = 0, g/Δ = 0.1 (`e1_transport.out`):

| rule | t/Δ | ground S | ⟨Σc⟩/n | gap/J | L1 miss | K1 |
|---|---:|---:|---:|---:|---:|---|
| T0 | 0 | 0 | 0.044 | 0.744 | 3.2% | PASS |
| T1 | 0.05 | 0 | 0.046 | 0.740 | 2.7% | PASS |
| T1 | 0.10 | 0 | 0.052 | 0.729 | 4.1% | PASS |
| T2 | 1.88 | **4** | 0.333 | — | 63% | FAIL |
| T3 | 1.88 | 0 | **0.492** | — | 220% | FAIL |
| T4 | 1.88 | 0 | 0.044 | 0.744 | 3.2% | PASS, identical to T0 |

- **T2 makes the converted ferromagnet,** as in L5(b).
- **T3 lets a converted cell, no longer full, hop,** so the ground
  state converts.
- **T4 refuses every hop** while R6's local law holds, so it is T0
  with no label (algebra in `PREREG.md` §2, confirmed).

### E2. Leak into empty outside

n = 4 ring at the R7 ground state, an empty lead from cell 0, t = 0.9.
Probability that no quantum has escaped by T = 300 (`e2_*.out`):

| rule | C-lock | C-geo | C-held |
|---|---:|---:|---:|
| T0 | **1.000000** | **1.000000** | **1.000000** |
| T2 | 0.0001 | 0.0001 | 0.0003 |
| T3 | 0.0004 | 0.0000 | 0.136 |
| T4 | **1.000000** | 0.0000 | **1.000000** |

- A leaking pair loses its first quantum at rate 0.006–0.025 per time
  unit, which is 0.4–1.5/J.
- **T3 leaks even under C-lock.** A cell depleted by virtual conversion
  is below its cap, so its remaining quantum walks out.
- **T4 × C-geo leaks.** Cell 0 has degree 3, so N + Σc = 2 < 3 never
  freezes it.
- C-held forbids 0.35% of the R7 ground state's weight (states with two
  space quanta on one link). Dropping that weight raises the ring's
  energy by 0.0019 (0.2 J), measured with R7's Hamiltonian.

### E3. A photon arrives

Run only on the five pairs that keep the lock in E2 (`e3_*.out`):

| pair | cold photon (ε = −1.72) | resonant photon (ε = 0.2275), X1 |
|---|---|---|
| T0 or T4, C-lock | reflected; ring fidelity 1.000000; excitation 5e-9 | 5% held in the boundary resonance at the read time; ring fidelity 0.9964 |
| T0, C-geo / C-held | the same | the same, within 0.5% |
| **T4, C-held** | the photon never lands next to the lock | **no interaction: X1 is inert** |

Under T4 × C-held, the cell beside the lock gets cap 1 the moment a
quantum arrives. A quantum arriving would freeze it, and since the rule
is checked in both directions, the hop is refused. **So that pair has
no cooling channel.**

### E3b. What each pair does to free light (added; see §9)

Two photons crossing on an open chain, with no lock anywhere
(`e3b.out`):

| rule | C-lock | C-geo | C-held |
|---|---|---|---|
| T2 | pass 0.9945 | pass 0.9945 | pass 0.9945 |
| T4 | pass 0.9945 | **cannot pass** (hard-core at occupancy z) | **cannot pass; an adjacent pair is frozen forever** (1 reachable state) |

**C-held × T4 stops light.** Any two adjacent quanta freeze. A beam is
many adjacent quanta.

### E4. The geometric ceiling

Under C-geo a cell at its ceiling holds z quanta, its degree, so its
spin is z/2 (`e4_a.out`, `e4_c.out`, `e4b_rings.out`, `e4_b.out`).

| z | spin | Δ_z | J (two-cell cluster) | spin-z/2 chain |
|---:|---:|---:|---:|---|
| 2 | 1 | 0.479 | 0.00932 | gap 0.4105 J |
| 3 | 3/2 | 0.788 | 0.00573 | **gapless** (Lieb–Schultz–Mattis) |
| 4 | 2 | 1.002 | 0.00452 | gap ≈ 0.089 J, ξ ≈ 49 (literature) |
| 6 | 3 | 1.287 | 0.00352 | gap ≈ 0.010 J, ξ ≈ 640 (literature) |

- **Mixed degrees make a ferrimagnet.** An alternating (z = 2, 3)
  ring, solved exactly for n = 6, has ground-state spin 1.5. That is
  Lieb–Mattis's n/2 × |1 − 3/2|, with gapless spin waves.
- Exact rings up to n = 8 cannot yet separate S = 1, 3/2 and 2: gap × n
  is 4.0–4.7 for all three. Open S = 3/2 chains do separate: the gap
  falls from 0.207 (L = 8) to 0.085 (L = 16). The S = 2 and S = 3 gaps
  are literature values.
- So C-geo fails twice:
  - at N = 2 the lock leaks (E2);
  - at N = z the lock loses its gap, or becomes a ferrimagnet.

## 3. Q2 — space flow on lock links

| option | test | result |
|---|---|---|
| F0 | — | holds (R7) |
| F0 globally (`s_k = 0`) | battery | **3 bars move**: p1 COE ratio, xsec absorption 6.849 (bar 6.9–7.65), xsec inert-lens floor 0.666 (bar 0.70–0.88) |
| F1 along the fiber, with T0 (n = 6) | `t_s` = 0.03, 0.06 | PASS: L1 miss 3.3%, 3.8%; ⟨Σc⟩/n 0.045, 0.049 |
| | `t_s` = 0.12, 0.24 | FAIL: ⟨Σc⟩/n 0.12, 0.41; L1 miss 41%, 84% |
| **F1 along the fiber, with T4** (n = 4) | any `t_s ≥ 0.03` | **collapse**: ⟨Σc⟩/n 0.55; L1 miss about 100%; the whole ground state has T4's predicate open |
| F1 lead (L5a, carried) | 1-D lead or tree z = 12.6 | holds (+0.36, +0.07 above the pair) |
| | lattice-edge estimate | decays |
| F2 rate law, γ = s_k | n = 4, 6, 8 | lifetime 0.83, 0.58, 0.44 /J |
| | γ = s_k·s_disp | 2.8, 1.9, 1.5 /J |
| | energy per jump | **creates +0.61 to +0.77 (1.3–1.6 Δ)** |
| F3 | — | identical to S1 wherever conversion can act, since S3 needs field at both ends |

- **On a ring with space flow along it, the threshold is not the band
  edge.** It is `2t_s ≈ Ω − E_f(2) = 0.2275`, the cost of each further
  adjacent conversion. Mobile space makes conversion strings cheaper than
  their cost. The `s_k` reading sits at 0.12, a margin of 1.9.
- **T4 needs no space flow at all on lock links.** Space flow moves a
  converted quantum off its link and breaks R6's local law. T4's frozen
  predicate then opens, and t = 0.9 converts the ring. T0 holds under
  the same flow. **This is the combination that separates T0 from T4.**
- **F2 is excluded by the one law.** Each jump leaves 1.3–1.6 Δ
  unaccounted, and the lock dies in under 3/J.

## 4. Q3 — the door

| option | measured |
|---|---|
| D0, work function open (the cell sits in a space deficit) | halves the lock cell's field in **3 steps** (0.06 time units, 5e-4/J) |
| D0, work function closed (no deficit): threshold `wf_far = 99` against 2.67 | **does not fire.** The lock is exempt with no label |
| D2 globally (`f_conv = 0`) | **10 bars move**: all 7 ds1 click bars, pauli0 at-cap admission, xsec absorption (−0.37), xsec monotone cross-section |
| D2 detection, screen link F0 | absorption 0 |
| D2 detection, screen link F1 (`t_s = s_k`) | peak 7–30% (m = 1, 2; g/Δ = 0.1, 0.2), FWHM 0.21–0.22 around ε = −0.024 or +0.2275; band average 0.4–1.7% |
| D2 detection, screen link R-c (continuum) | band average 0.35–2.7%, peak at most 3.8% |

`pattern_only`'s work function already answers the two-conversion-laws
problem for a lock that carries no space deficit. Whether a lock carries
one depends on how the S-law's link oscillator maps onto the kernel's
scalar `Es`, and that is a kernel-step question. D2 turns detection into
a weak spectral filter, which the click claims cannot survive.

## 5. Q4 — g

n = 6 ring and the S-law clusters (`e8_static_*.out`, `e8_margin_*.out`):

| g/Δ | J / (2g²/Δ) | K/J | r3/J | ⟨Σc⟩/n | L1 miss | hold, 2nd half | t_c/Δ |
|---:|---:|---:|---:|---:|---:|---:|---|
| 0.025 | 0.998 | −0.001 | 0.004 | 0.003 | 0.2% | 0.995 | (0.10, 0.15] |
| 0.05 | 0.990 | −0.005 | 0.015 | 0.012 | 0.9% | 0.978 | (0.10, 0.15] |
| **0.10** | 0.963 | −0.019 | 0.052 | 0.044 | **3.4%** | 0.922 | (0.15, 0.18] |
| 0.125 | 0.944 | −0.027 | 0.078 | 0.065 | **5.06%** | 0.890 | (0.18, 0.20] |
| 0.15 | 0.923 | −0.036 | 0.108 | 0.086 | 6.9% | 0.859 | (0.20, 0.25] |
| 0.2 | 0.877 | −0.051 | 0.175 | 0.128 | 10.7% | **0.804** | (0.20, 0.25] |
| 0.3 | 0.781 | −0.072 | 0.299 | 0.194 | 18.3% | 0.739 | (0.30, 0.40] |
| 0.5 | 0.618 | −0.084 | 0.459 | 0.267 | 29.4% | 0.690 | (0.40, 0.50] |

- The ground state is a singlet at every g. K/J stays far from −1, so
  the chain is in the Haldane phase throughout.
- The registered K1 window is **g/Δ ≤ 0.1**. At 0.125 the line misses
  by 0.06%; at 0.2 the hold fails.
- The transport margin grows with g.
- **Formation pulls the other way** (§7). Cooling needs light colder
  than the gap, and the gap is ∝ g². The photon count falls roughly as
  g⁻⁴: about 3×10⁸ at 0.1 and 1.7×10⁷ at 0.2. At g/Δ = 0.05 no source
  tried (δ ≥ 0.005) was cold enough.

The readings against the window:

| reading | g/Δ | place |
|---|---:|---|
| `mob_floor` | 0.008 | inside; cooling would need light colder than about 10⁻⁵ and about 10¹³ photons |
| `lock_floor` | 0.010 | inside; the same |
| `s_k` | 0.125 | **at the edge**: the line misses by 0.06% |
| `gamma_res` | 0.52 | outside |
| ħΩ | 1.11 | outside |
| `field_J/2` | 1.88 | outside |

No reading is chosen.

## 6. Q5 — size

| option | measured (`e9_*.out`) | reading |
|---|---|---|
| Z0, free rings | E₀/n lowest at n = 4 (L7) | the smallest even cycle |
| **Z1, coalescence with X1** | cube −0.24514 against two 4-rings −0.18298: **−0.0622 (6.8 J)**. The cube is a singlet, gap 0.0071 | **joining wins: locks hoard** |
| Z1 at leading order | per cell −1.50 (ring), −1.95 (cube), −1.44 (6-ring), −1.91 (hexagonal prism); gap 1.00, 0.78, 0.72, 0.56 | more links per cell, lower energy, smaller gap. A 3-D aggregate orders as a Néel solid with no gap |
| Z2, twist θ = π | cost 0.0056, 0.0028, 0.0016 for n = 4, 6, 8, about 1/n² | moves the boundary only if small rings twist and large ones do not. A uniform twist favours larger rings. No interior minimum |
| Z3, per quantum | rings 1.3248–1.3253; open chains 1.385 (n = 3) down to 1.339 (n = 10); a lone full cell 1.3363; free light at its band bottom −0.338 (1-D), −3.938 (cubic) | rings beat chains, and all are far above free light |

No option selects an interior size.

- With X0 between locks, the ring sits on the smallest even cycle.
- With X1, touching locks merge without limit.

## 7. Q6 — formation (E5)

n = 4 lock, one boundary link, exact multichannel S-matrix, from L6's
hot Néel start (ground population 0.17). Unitarity holds to 1e-15
(`e5_g*_*.out`).

| g/Δ | ring gap | photon source | outcome per photon | photons to ground ≥ 0.9 | destruction by then |
|---:|---:|---|---|---:|---:|
| any | — | **X0** | the identity: nothing happens | never | 0 |
| 0.05 | 0.0024 | cold, δ = 0.005 | the stationary ring is warm (p0 = 0.48) | never | 0 |
| 0.1 | 0.0095 | cold, δ = 0.005 | cools | ≤ 2.7×10⁸ | 0 (F0); 1e-8 (F1) |
| 0.1 | | cold, δ = 0.01 | cools (stationary 0.93) | ≤ 5.4×10⁸ | ~0 |
| 0.1 | | δ = 0.02 | stationary 0.47 | never | 0 |
| 0.1 | | **δ = 0.05 (registered)** | stationary 0.09 | **never** | 0 |
| 0.2 | 0.0352 | cold, δ ≤ 0.02 | cools | ≤ 0.8–1.7×10⁷ | ≤ 2e-9 |
| 0.2 | | δ = 0.05 | stationary 0.69 | never | 0 |
| 0.05–0.2 | | **resonant** (ε_res ± 0.05) | excites 9–40% per photon; with F1 also destroys a quantum 4–25% per photon | never: heats into converted levels (E − E0 = 1.0–1.7) | — |
| 0.05–0.2 | | broadband | heats (p0 → 0.04 or less); F1 destroys 0.2–1.4% per photon | never | — |

**A lock forms only with boundary conversion (X1), and only from light
colder than its own gap.** That is about 10⁸ photons at g/Δ = 0.1, with
negligible destruction. Light near the boundary resonance heats it.
With boundary space flow, resonant light destroys it. Under the
registered δ = 0.05, K2 fails at every g.

## 8. The scorecard

### Options in isolation

| option | K1 lives | K2 forms | K3 table | K4 cost | K5 margin |
|---|---|---|---|---|---|
| T0 | ✔ | — | label confines it | 1 label | — |
| T1, 0.05–0.10 Δ | ✔ | — | **16–18 bars** if global | 1 constant (+ a label to confine it) | t_c/t = 1.7–3.4 |
| T2, T3 | ✘ | | | | |
| T4 | ✔ (with C-lock or C-held, no space flow) | | with C-lock ✔; with C-geo/C-held ✘ (E3b) | 0 | **0 under any space flow** |
| C-lock | ✔ | ✔ with X1 | ✔ | 1 label | |
| C-geo | ✘ (leak; spin z/2) | | ✘ under T4 | 0 | |
| C-held | ✔ with T4 | ✘ (photon cannot land) | ✘ under T4 (light freezes) | 0 | |
| X0 | ✔ | ✘ | ✔ | needs the label to refuse | |
| X1 | ✔ | ✔ (cold light) | ✔ | 0 | merges touching locks |
| F0 on lock links | ✔ | | 3 bars if global; label if local | label | |
| F1 at `s_k` along the fiber (T0) | ✔ | | ✔ | 0 | 1.9 |
| F1, lattice-edge | ✘ | | | | |
| F2 | ✘ (≤ 2.8/J; creates energy) | | | | |
| D0, deficit | ✘ (3 steps) | | | | |
| D0, no deficit | ✔ | | ✔ | 0 | depends on representation |
| D2 | ✔ | | **10 bars; detection ≤ 30% on a line** | | |
| g/Δ ≤ 0.1 | ✔ | ✔ at 0.1, with δ ≤ 0.01 | | 1 constant | 1.0–1.25 to the L1 edge |

### Bundles

| bundle | K1 | K2 | K3 | K4 | K5 | verdict |
|---|---|---|---|---|---|---|
| **B-label-cool** (T0, C-lock, X1, F0 on lock links, door exempt) | ✔ | ✔ with cold light; ✘ at the registered δ | unmeasured (the label confines every rule) | 1 label | space flow ≤ 1.9 × `s_k`; g at the window's edge | **admissible, rank 1** |
| B-local-lock (T4, C-lock, X1, F0 inside) | ✔ | ✔ the same | unmeasured | 1 label | **no tolerance of space flow on lock links** | admissible, rank 2 |
| B-label (T0, C-lock, X0) | ✔ | ✘ (identity map) | unmeasured | 1 label | as rank 1 | admissible, rank 3 |
| B-local (T4, C-held, X1) | ✔ ring | ✘ | **✘ light freezes** | 0 | — | not admissible |
| B-geo (T4, C-geo, X1, F1) | ✘ | | ✘ | 0 | — | not admissible |
| B-table (T2, C-geo, X1, F2, D0) | ✘ ×4 | | | 0 | — | not admissible |

The ranking follows `PREREG.md` §5:

- K2 puts B-label below the two cooling bundles.
- K4 ties those two at one label each.
- K5 separates them. T4 tolerates no space flow on lock links, and T0
  tolerates `s_k` with a margin of 1.9.

## 9. Deviations from the pre-registration

1. **E5's cold source.** `PREREG.md` registered one cold source,
   δ = 0.05. Before E5 ran, I added δ = 0.005, 0.01 and 0.02, because
   δ = 0.05 is 5 J, hotter than the n = 4 lock's gap (0.0095). This was
   added in the code, not written into `PREREG.md`. K2 is reported at
   both: **it fails at the registered δ for every g, and passes with
   the colder sources at g/Δ = 0.1 and 0.2.**
2. **E5's numerics.**
   - The photon count first used an eigendecomposition. That was
     unstable at 10¹³ photons, so it was replaced by repeated squaring
     of the augmented (substochastic) map.
   - Without boundary space flow, columns are renormalized to 1. The
     largest correction is 2e-15 after averaging; it absorbs numerical
     unitarity slips near infinitely narrow in-band resonances.
   - Every printed number is from the final code.
3. **E3 bug, fixed.** The first E3 runs launched the photon away from
   the ring, at the wrong end of the band. The engine's hop is +t, so
   ε = +2t cos k. Those outputs were deleted, and every E3 number above
   comes from the corrected rerun.
4. **E2 bug, fixed.** Under C-held, seeds that violate the ceiling were
   first admitted, which broke Hermiticity (norm 1.16). Seeds are now
   filtered. The forbidden weight (0.35%) is dropped and printed.
5. **E3b was added.** It was not registered. It tests free light under
   each pair, and it is reported under K3.
6. **E4(b).** The open-chain DMRG was too slow for its question. Exact
   rings (n ≤ 8) cannot separate the spins, so the gapless and small-gap
   entries rest on the Lieb–Schultz–Mattis and Lieb–Mattis theorems and
   on literature values. They are labelled as such. A reduced DMRG of
   open S = 3/2 chains is in `e4_b.out`: the gap is 0.207 at L = 8 and
   0.085 at L = 16, which is falling at least as fast as 1/L. The S = 1
   chain with spin-½ ends (`../lock/` L4) holds 0.45–0.50 at these
   lengths.
7. **E6(a)'s T4 combination** ran on n = 4 (`e6_t4f1.py`). The n = 6
   run was too slow and was stopped after its t_s = 0 row.
8. **E8 added g/Δ = 0.125.** That is the `s_k` reading, which E8's
   "report which G-read values fall inside" required.
9. **E9(a)'s hexagonal prism** was run at leading order (spin-1
   Heisenberg), not with the microscopic S-law. The microscopic version
   exceeds 10⁶ states. The cube is microscopic and exact within Σc ≤ 3.
10. **The table side** (E6d, E7b, E10) ran at `pattern_only = 0`, the
    battery's default representation.

## 10. What is now true, and what is not

**True (measured here):**

1. **No label-free ceiling keeps a lock and leaves light alone.**
   - C-geo leaks the lock, or makes spin z/2 cells (ferrimagnetic or
     near-gapless).
   - C-held freezes any two adjacent quanta of free light.
   - The lock needs one label: which cells are in it.
2. **Given that label, t = 0 on lock links (T0) is the robust choice.**
   The local rule T4 is identical when nothing else moves, but it
   collapses under any space flow on lock links.
3. **A lock forms only through boundary conversion, and only from light
   colder than its gap,** in 10⁷–10⁹ photons at g/Δ = 0.1–0.2. Resonant
   light heats it, and with boundary space flow it destroys it.
4. **The lock is metastable.** Every lock sits far above free light at
   its band bottom, per quantum. It exists because transport is refused,
   and that refusal also means transport cannot fill it.
5. **Space may flow on a lock only coherently, and slowly:**
   `2t_s < Ω − E_f(2) = 0.2275` along the fiber. Rate-law flow creates
   energy and kills the lock within 3/J.
6. **The door is harmless to a lock that carries no space deficit.**
   Removing the door costs 10 bars and leaves detection weak and
   line-selective.
7. **The registered window for g is g/Δ ≤ 0.1.** The table's `s_k`
   reading sits on its edge.
8. **With boundary conversion, touching locks merge.** No option
   selects a size.

**Not true:**

- Any of this runs in a kernel, or the battery passes with it.
- A label-free lock exists under the options tested.
- g is set.
- A size is selected.
- A lock forms from the table's own light. Its spectrum and temperature
  were not part of this test.

## 11. What would move past this page (each is the user's call)

1. **Rule the label.** For example: "a lock is the set of cells at
   their ceiling joined by converting links", carried as a state flag.
   Without some such statement, no tested rule set keeps the lock.
2. **Rule what happens when locks touch.** Merging is what X1 does.
   Keeping them apart needs lock identity, a second label.
3. **Formation.** Decide whether a light source colder than the gap is
   part of the world, or test lock birth in place: a pair converting out
   of space, which T0 does not forbid. Birth in place was not tested.
4. **A kernel step** (needs authorization). The door's work function
   and the label's representation decide whether D0 spares the lock.

## 12. After the campaign: "cold" is a layer mismatch (user, 2026-10-08)

User: "as I suspect for the link which must be cold; heat is not an
aspect to it if you modify it, it is underlying energy, not capable of
having energy, it is a part of what we think of energy, layer
mismatch."

**My reading. Nothing here is measured, and it is not a ruling.**

1. **The fiber and its links are the layer that energy is made of.**
   They cannot be hot. Asking them to be cold puts a property of the
   energy layer onto the layer beneath it.
2. **By R5, the lock is the fiber's gapped quantum.** A fiber above its
   ground state is not a hot lock. It is a fiber carrying lock quanta.
   - So K2 as registered (ground population ≥ 0.9) measured emptying
     the fiber of its quanta, not making a lock.
   - Two statements on this page are about emptying the fiber, and
     should be read that way: §1's "formation needs light colder than
     the lock's gap", and §10 item 3. The same goes for
     `../lock/DECIDE.md` L6.
3. **The model's own layer mismatch is S1's form.** The link's space
   mode was given its own energy ladder, Ω per quantum, that can be
   really occupied. Every failure in this campaign that runs through
   real link occupation sits on that side:
   - resonant light pumping the boundary oscillator (E5: it heats, and
     with F1 it destroys);
   - space flow moving a real link quantum (E6: F1 above 0.2275, T4's
     collapse, F2's jumps and the energy they leave unaccounted);
   - L6's 10–20% weight in space-carrying levels after a sudden start.

   The exchange needs the link only virtually.
4. **Unchanged:**
   - the label result (§2);
   - the door (§4);
   - the g window (§5);
   - size (§6);
   - the fiber's metastability against free light (§10 item 4).
5. **Testing it would take two things:**
   - a link with no real occupation. The effective chain h2 + r3
     (`../lock/cas/clusters.py`) is that model to fourth order in g;
   - a formation test that counts lock quanta, not ground population.

   Neither has been run.

## Reproduce

From `v102/options/cas/`. The battery variants are built into scratch
and run with `-extra`.

| step | command | output |
|---|---|---|
| validate | `python3 -I val_graph.py` | prints matches against `../../lock/cas/slaw.py` (≤ 2e-14) |
| E1 | `python3 -I e1_transport.py` | `out/e1_transport.out` |
| E2, E3 | `./run_e2.sh`, `./run_e3.sh` | `out/e2_*.out`, `out/e3_*.out` |
| E3b | `python3 -I e3b_freelight.py` | `out/e3b.out` |
| E4 | `python3 -I e4_ceiling.py a`, `... c`, `python3 -I e4b_rings.py`, `python3 -I e4_ceiling.py b` | `out/e4_*.out`, `out/e4b_rings.out` |
| E5 | `./run_e5.sh` | `out/e5_g*_*.out` |
| E6 | `python3 -I e6_space.py a`, `... c`, `python3 -I e6_t4f1.py` | `out/e6_*.out` |
| E7 | `python3 -I e7_detect.py` | `out/e7.out` |
| E8 | `./run_e89.sh`, plus `e8_g.py {static,margin} 0.125` | `out/e8_*.out` |
| E9 | `./run_e9.sh` | `out/e9_*.out` |
| table | `battery -kernel c -extra "<key>=<value>"` for `s_k=0`, `f_conv=0`, `field_J=0.18`, `field_J=0.09` | `out/bat_*.out` |

Instruments:

- `graph.py`: the S-law on any graph, with the option rules (ceilings,
  symmetric transport rules, twist, space flow).
- `../../lock/cas/clusters.py`, `mpo3.py`: carried.
