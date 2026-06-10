# v66 FINDINGS — Q-Ball Production Campaign (GPU)

**Date**: 2026-06-10. **Status**: [measured] at N=192–256, V100, f16 snapshots; diag-level
quantities at 12-digit charge precision. All runs seeded from `radial_qball` profiles
(exact η=0 solutions), standard parameters m²=2.25, m_θ²=0, μ=−41.345, κ=50, absorbing BC.
Data: `v66/results/*_diag.tsv`, SFAs in `/space/scp/v66/`. Theory: `v66/THEORY.md`.
Build verification: `v66/results/SMOKE_TEST.md`. Two-ball cluster tables:
`v66/results/twoball_clusters.md`.

## 1. The η=0 Q-ball is effectively eternal

`r1_eta0` (ω=1.39, N=192, absorbing BC, T=1000): Q = 482.18465 → 482.18394
(**retention 99.99985% over 1000 t.u.**), E=691.67 constant, r_core 3.704 frozen,
s_max 0.0688 with only the benign ~0.94 rad/t lattice breathing. No 3ω harmonic
(the v65 decay channel) at any point. Combined with the machine-floor Q-drift
(7×10⁻¹⁴/20 t.u. periodic), the bare Q-ball is stable on every accessible timescale.
**CONCEPT.md §6 items 1 (persistent core under absorbing BC) and 2 (conserved charge)
are now established for the complexified theory.**

## 2. Universal attractor at the edge of the stable branch (η=0.5)

All stable-branch seeds converge to the same dressed end state:

| run | seed ω | Q₀ | endpoint Q | endpoint E/Q | r_core |
|---|---|---|---|---|---|
| b1 | 1.36 | 1745 | 876 @ t=600 (still draining, ON branch) | 1.4175 | 4.45 |
| r2 | 1.39 | 482 | **171.7** @ t=1000 (parked) | 1.4972 | 3.187 |
| b2 | 1.42 | 210 | **173.1** @ t=600 (parked) | 1.4966 | 3.186 |

r2 and b2 agree to 0.8% in Q and 0.04% in E/Q from very different starting points.
Trajectories track the equilibrium Q(ω) branch from the radial scan throughout
(adiabatic branch-sliding; measured E/Q matches scan values at matched Q to <0.5%).
b1 (3.6× heavier) rides the same branch but drains *slower* per unit charge in the
thin-wall regime — time-to-park is strongly seed-dependent even though the destination
is not. **The marginally-bound minimal ball (Q∞ ≈ 170 at η=0.5) is the particle the
full theory selects** — a mass prediction, not an initial condition.

At η=0.75 the attractor shifts: Q∞ ≈ 156 with E/Q = 1.512 — ABOVE the bare E=mQ line.
The end state is a **dressed ball** (θ cloud carries part of E and Q); the attractor is
mildly η-dependent and the bare-branch threshold is not the dynamical boundary at
finite coupling.

## 3. Stability boundary: VK criterion + shock fragility, NOT the E=mQ line

The b3 seed (ω=1.46, E/Q=1.52 > m, "evaporation-permitted" but VK-stable) produced the
campaign's most instructive sequence:

- **b3 (η=0.5)**: core dissolved by t≈100 (s_max → 0), total dispersal, Q radiated away.
- **c1 (η=0 control)**: Q = 102.5169 frozen to 5 decimals for T=300, E/Q=1.5207 static.
  **Classically immortal despite E > mQ.**
- **c2 (η=0.1 control)**: survives the coupling switch-on; slow drain −0.014/t.u.
  (Q 102.5 → 97.1 @ t=400), no ejection.

Conclusion: the static evaporation bound E<mQ does NOT govern classical dynamics.
The bare theory is stable on the full VK branch (dQ/dω<0, ω<1.485). At finite η the
**dressing shock** (the transient when the θ cloud forms) ejects fragile edge-of-branch
balls: the ω=1.46 ball survives η≤0.1 but is destroyed at η=0.5, while mid-branch balls
(ω=1.39) survive even η=0.75. The dynamical stability window narrows with coupling
strength via shock amplitude, not via energetic bookkeeping. [Frequency-domain
diagnostics (ω_core(t), radiation spectrum) are the planned instrument to resolve the
ejection mechanism — task list.]

## 4. Drain law: perturbative η² with strong-coupling saturation

Peak |dQ/dt| in the matched window t∈[50,100], ω=1.39 seeds:

| η | 0.1 | 0.25 | 0.5 | 0.75 |
|---|---|---|---|---|
| peak −dQ/dt | 0.106 | 0.595 | 1.72 | 2.36 |
| local exponent | — | 1.88 | 1.53 | 0.78 |

The weak-coupling exponent (1.88 ≈ 2) confirms the perturbative prediction
(drain ∝ η², THEORY §5b); the exponent collapses toward ~0.8 by η=0.75 —
nonlinear saturation, consistent with back-reaction of the θ dressing on the source.
Global power-law fit over the whole range: η^1.54 (a regime average; do not use).

## 5. Two-ball interactions: the phase-force taxonomy

N=256, L=30, D=12, ω=1.39 pairs, η=0.5 (cluster analysis on phase-invariant ρ₂;
`twoball_clusters.md`):

- **Δφ=0 (co-phase): FUSION.** Single cluster by t=100; merger product relaxes
  prolate→spherical and drains down the standard branch (endpoint Q≈+147, near the
  attractor). The diag "bound oscillation" was the merger product breathing.
- **Δφ=π (anti-phase): EJECTION.** Mirror escape, D: 12 → 40.6 by t=500 at ~0.05–0.07c;
  two intact balls of Q≈+118 each. (Minor mirror-symmetric y=z drift — suspected cubic
  lattice artifact, check at higher N before interpreting.)
- **Ball/anti-ball: STEPWISE ANNIHILATION → NEUTRAL REMNANT.** Initial Q_total = 0
  held to 10⁻¹³ throughout. Weak 2ω-averaged attraction → collision at t≈100 with
  complete core dissolution (s_max ×40 down; frame captured) and an energy burst
  ~3× baseline drain; pair re-forms with only ±75 charge segregation (20% of seed),
  recollides ~every 100 t.u., and settles into a **single net-neutral two-lobed
  breathing object** (lobes Δx≈2.1, prolate) persisting ≥200 t.u. with 37% of the
  initial energy. A long-lived charge-neutral composite is a new object class for
  this theory.

The phase-dependent force confirmed in both signs is the canonical Q-ball cos(Δφ)
interaction — the clean realization of what V34's winding-dependent force measured
on unstable oscillons. Initial-charge interference (+0.8% co-phase, −2.3% anti-phase
vs 2× single) is visible at t=0 as a free consistency check.

## 6. Verification chain (one line)

Maxima 19/19 → radial solver vs theory <0.1% → CPU real-limit 2×10⁻¹⁵ → Q-drift
7×10⁻¹⁴ → GPU/CPU parity 2×10⁻¹³ → 5-lens adversarial review → smoke test → campaign.

## 7. Open items spawned

1. **Frequency-domain diagnostics** (ω_core(t), θ radiation spectrum) — explain the
   drain saturation and the shock-ejection mechanism. [task #9]
2. **Radiation-bath equilibrium** — does an ambient bath cancel the drain (detailed
   balance) and shift Q∞? The η-dependent dressed attractor says yes in kind. [task #7]
3. **Condensate fragmentation under pressure** — Affleck–Dine-style natural Q-ball
   formation; star-interior compression; 2.5D slab as later kernel extension. [task #8]
4. **Neutral remnant characterization** — is the tb3 endpoint a true bound state?
   Lifetime, spectrum, response to perturbation.
5. **Shock-threshold map** — η_crit(ω) ejection boundary; finer η scan near the
   branch edge.
6. **Anti-phase lattice drift** — verify the y=z drift is a discretization artifact.

## Infrastructure notes

- Vast.ai instances died twice mid-campaign (~4.5 h uptime each, SSH reset, total loss
  of undownloaded files). Mitigations now standard: chain runs in one self-serializing
  script; download eagerly after each run; auto_download does NOT track script-mode
  outputs (only config-mode runs) — do not rely on it for chained scripts.
- The burst snapshot mode (burst_every≈0.5 with 24 columns) is I/O-bound and stalls
  GPU compute; burst frames appeared much larger than regular f16 frames — kernel I/O
  quirk worth investigating before using burst on complex runs.
- `gen_qball_pair` (sfa/seed/) stamps two-ball seeds with arbitrary phases/charges;
  validated to Q-drift 2×10⁻¹³. The phase-invariant ρ₂ cluster tracker (built ad-hoc
  in /tmp) should be promoted into sfa/analysis/.
