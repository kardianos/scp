# B1 — Winding-compensated ring seed (design note)

*2026-07-28. Produced by a parallel design agent (v89-isolated, read-only);
archived verbatim-in-substance. Status: probe fleet launched (comp12 /
unwound-matched / comp6 / naive12). MASS tech tree B1.*

## The algebra

Kernel facts (cellfab.c seeder + gates; laws_V2g: w2=2.9, q_detune=1.2,
C=1, gamma_res_m=0.10, p_gate=8, cap=2.5):

* Pitch from load: ω = w2/(1 + q·x) ⇒ x = (w2/ω − 1)/q.
* Directional gates (1:1): ψ_fwd = θ_i − ω_i·d/C − θ_j (and mirrored);
  gate = ((1+cos ψ)/2)^8; the two directions are independent.
* Both gates open together iff the pair ladder holds:
  (ω_i + ω_j)·d/C = 2πm_link.
* Equal pitches ⇒ comb det(1:1)=0 ⇒ NO internal roughness regardless of
  winding — the winding cost is purely in the gates.

**Key identity: winding IS the choice of per-link retardation.** With
forward gates open, θ_{k+1} − θ_k = −ω_k·d_k/C (mod 2π). Loop
single-valuedness gives Σ ω_k d_k / C = 2πm; with uniform links,

    φ := ω·d/C = π − 2πw/N  =  2πm/N,   m = N/2 − w   (N even).

The naive `ring_wind=w` seed already produced exactly these phase
increments — the naive and compensated rings have IDENTICAL seeded phase
patterns. The naive seed fails only because it leaves the loads on the
π-rung: every gate then sits at ψ = −2πw/N (N=6, w=1: 60°, gate ≈ 0.100
— the measured 2× faster leak). Fix: retune the LOADS so the retardation
owns the twist, and delete the phase kick.

## (a) Twist in the retardation — RECOMMENDED

    φ = π − 2πw/N          (closure integer m = N/2 − w)
    ω = φ/d                (d in the foam link window)
    x = (w2·d/(φ·C) − 1)/q_detune    (uniform load)
    θ_{k+1} = θ_k − ω·d_k/C          (pure lock recursion, ring_wind=0)

Predicted: forward gates 1.000 exactly (seeder matches actual per-link
retardations); back gates ((1+cos(4πw/N))/2)^8; seam link carries the
closure defect (accept seeds with printed closure within 0.05 of m).

Numbers (w=1):

| | N=6 | N=12 |
|---|---|---|
| φ | 2π/3 (m=2) | 5π/6 (m=5) |
| d | 1.10 (1.25 rejected: ω would sit 0.026 from w1=1.65) | 1.25 |
| ω | 1.9040 | 2.0944 |
| ring_x | 0.43593 | 0.32054 |
| fwd gate | 1.000 | 1.000 |
| back gate | 1.53e-5 (fully one-way!) | 0.1001 |
| closure/2π expected | 2.000 | 5.000 |

**Protection (the topology argument):** the seeded state shares the
frustration equally; to unwind, a link's phase increment must cross the
dead-gate desert (gate ~4e-10 at ψ≈−2.6) — full unlock and refire, not a
drift. Winding is protected by the desert, not by gate stiffness (the
gate top is flat: gate ≥ 0.9 for |ψ| ≤ 0.16).

**Risks:** (1) N=6 w=1 has dead back gates — one-way circulation is
itself the physics question; N=12 hedges (back 0.100 = no worse than
every gate of the measured naive ring). (2) Pair-temper frustration per
link — watch for slips as census mass steps. (3) Heavier voices than the
unwound ring — mass-matched controls mandatory. (4) **Diagnostic
mirage:** gg reads 0.100 (N=12) or 1.5e-5 (N=6) BY DESIGN while the ring
circulates at full forward rate — score by leak + circulation (flux
moment / DBG gate asymmetry = the chirality signature), never by gg.

**Seeder upgrade worth adding later: `ring_m`** — seed by target closure
integer (ω = 2πm·C/L_ring after cell picking): kills the seam defect
exactly and makes m the species label directly.

## (b) Twist in non-1:1 rungs — fallback / species axis, not the fix

Balanced interval multisets (Π q/p = 1) can carry winding in lock-branch
choices, but acceptance gw = Γ_m/(pq) is 6× narrower and on-tooth res
1/(pq) 6× weaker for the fifth-cycle; measured lifetime ordering says
unison ≫ fifth. Use as the B2 species axis later.

## (c) Double rings

Counter-rotating pair: sparse shear-prone axle (most axial gates dead) —
a B3 composite, not a lead. **Co-rotating tube (+w,+w): the rescue if
the single ring slips** — all axial links can sit on the exact m=1 pair
rung (both gates 1.000) with d_z = π·C/ω; restores mutuality without
giving up net winding. Needs a two-ring seeder (B1b).

## The probe fleet (zero code change) and scoring

comp12 (N=12 d=1.25 x=0.32054 wind=0; expect closure 5.00) vs
unwound-matched (N=12 d=1.5 x=0.32054 → m=6: same mass, same pitch,
differs ONLY in the closure integer — the clean topological A/B) vs
naive12 (x=auto wind=1) vs comp6 (N=6 d=1.10 x=0.43593; the one-way
extreme). T=1000, census cadence.

WIN: comp12 leak ≤ unwound-matched with sustained circulation → chase
comp6 and the w=−1 partner (m = N/2+1 — a LIGHT chiral species: the
chirality-split mass spectrum is the first charge-like observable).
PARTIAL: between unwound and naive → suspect pair-temper frustration →
go to the tube (c). KILL: comp12 ≥ naive → one-way closure worse than
detuned mutual closure; tube or (b) inherit.
