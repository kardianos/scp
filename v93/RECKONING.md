# v93 — the reckoning (R1/R2/R3): Face C is a TRANSIENT, not matter

2026-08-08. RECONCILIATION §1 short-term experiments, executed. All three
gates return the negative: the Face C condensation is a **transient self-
trapping that peaks at ~t=300 and relaxes away** — phase-incoherent and
energy-unbound. Face C caught the transient peak, not a steady state.

## R1 — Endurance (the gate): transient, both corners
Broad packet, Strang (hop_order=1), T=3000. Em_max / PR / rms:

| corner | t=0 | t=300 (Face C) | t=1000 | t=3000 |
|---|---|---|---|---|
| deep (amp2, qd3.6) Em_max | 2.25 | **7.24** | 4.87 | **2.67** |
| deep PR / rms | 456 / 4.30 | 250 / 4.51 | 273 / 5.02 | **345 / 5.85** |
| law (amp0.5, qd1.2) Em_max | 0.56 | **2.17** | 1.96 | **1.81** |
| law PR / rms | 456 / 4.30 | 203 / 6.01 | 255 / 6.76 | **336 / 7.38** |

Both corners: **peak at t≈300, then monotonic decay + spread.** Deep corner
decays faster (Em_max 7.24→2.67, −63%; evap=130 bleeds over-cap energy to the
field). Law corner decays slower, purely by diffusion (cond=evap=0; Em_max
2.17→1.81, −17%; PR 203→336; rms 6.0→7.4). **Neither holds.** t=300 was the
transient maximum, exactly where Face C measured.

## R2 — Isolation (redefined; bath=0 gives NC=0 in this substrate): tight droplet melts
A tight seed (σ=1, amp=2) self-traps to PR≈22 at t=20 — but it **melts** over
T=1000: Em_max 2.03→1.69, rms 1.73→4.90 (×2.8), PR 29→37. The early tightness
was transient; no natural soliton scale exists (tightness tracks the seed:
PR 22 for σ=1 vs 250 for σ=2.5). Low-energy tight seeds (amp=0.5) don't
condense at all (DNLS threshold).

## R3 — Phase audit (the hinge): phase-random
Natural condensate (seed_mw=0, uniform phase), union-find clustered: ONE
connected condensate (187 cells — the "189 hoards" of FACE_C were a cell-
occupancy histogram, not a particle spectrum). Within-hoard R2_m0 = 0.001
(deep) / 0.17 (law) — **at or below the bath noise floor** (0.023 / 0.11).
No hidden phase order, no meaningful winding. **q_detune is phase-blind; the
condensate is an energy-only object.** Confirms B.4 in the stronger form.

## Verdict — "record and pivot" (both reviewers' named negative)
The q_detune load-detuning **transiently concentrates energy** (a real effect
the additive law cannot reproduce — RING_DNLS B.5), peaking ~t=300, but it is:
- **not stable** (R1/R2: diffuses + evaporates back to seed levels by t=3000),
- **not phase-coherent** (R3),
- **not particulate** (one connected condensate, not a gas of particles).

So Face C is **creation-adjacent relaxation physics**, not a particle factory.
The programme's standing open problem — **coherent, stable binding of local
structure (the retention wall)** — remains THE problem. This reckoning closes
the "is the condensate matter?" question: no.

## Recommended next face (opencode path first, per user directive)
opencode's path (RECONCILIATION §3) names, given R3 phase-random, the **J-
dependent nonlinearity / dynamical Hodge star (v93 §II.5)** as the most
defensible first kernel change — already in the design language ("space curves
because energy is conserved and space is convertible"): the link metric
responds to the link current J = 2·Im(ψ_i*ψ_j), making the binding phase-
sensitive. **R1/R2 raise the bar for it:** the new mechanism must produce a
*genuinely stable* bound state (non-diffusing), not merely a phase-coherent
one — q_detune showed that transient concentration is not enough. Design it
behind a byte-inert knob (default 0), battery-gated, like v93's other faces.

Fallback candidates if the Hodge-star face does not yield stable binding:
- the **topological ring** (RING_DNLS A.1: holds W=+2 to t≥1000 under Strang
  in vacuum — the ONE genuinely stable object; its bath-survival is the open
  seam, A.2),
- or pivot to the **retention wall** directly (claude's path).

## NOT to do (unchanged from RECONCILIATION)
v75 three-fabric, imported charge, SI mass mapping, adoption of amp_nat>0 into
a law table, or any "particle/element" work — all gated on a confirmed stable
bound state, which does not yet exist.

## Files
- R1: `runs/r1_endurance.sh` → `runs/reckon/r1_{deep,law}.{fcs,melt,log}`.
- R2: `runs/reckon/r2_{s1.0_a0.5,s1.0_a2,tight_T1000}` (tight-seed probes + endurance).
- R3: `report/analyze_hoard_phase.py` (run on existing `runs/faceC/c_deep2.fcs`, `c_law05.fcs`).
