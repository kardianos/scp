# v68 FINDINGS — Relativity Ladder, Charged Bath, Resonance Scan, Gauge Design

**Date**: 2026-06-11. **Status**: [measured] 13-run V100 campaign (all rc=0) + two
design documents. Data: `v68/results/*_diag.tsv`, SFAs in `/space/scp/v68/`.
Design docs: `v68/GAUGE_DESIGN.md` (59/59 Maxima), `v68/DRESSED_STATE.md`.
New tools (committed v67 set + this round): `gen_qball_bath` charged mode,
`gen_qball_boost`. Fragment spectra: `v67/results/FRAGMENT_SPECTRA.md`.

## 1. The relativity ladder — exact through v=0.7 [measured]

Boosted balls (gen_qball_boost, exact DEBROGLIE solutions, η=0, N=128):

| v | γω predicted | omega_core measured (t=0) |
|---|---|---|
| 0.1 | 1.397010 | 1.397002 |
| 0.3 | 1.457116 | 1.457116 |
| 0.5 | 1.605058 | 1.605034 |
| 0.7 | 1.946394 | 1.946389 |

5–6 digit agreement across the full range; CPU validation additionally measured
translation v=0.2967 (vs 0.30), E/E₀ = 1.0467 (vs γ=1.0483, lattice), charge
Lorentz-invariance to 2×10⁻⁶, stable shape. **The Q-ball is an exact relativistic
particle in-simulation**: E²=p²+E₀², lab-frame clock γω, de Broglie phase machinery
all confirmed. (The ℏ_eff=Q phase-tilt fingerprint analysis can be done offline from
the banked bs SFAs — snapshot phase gradient across the core; open analysis item.)

## 2. Charged bath — the drain inverts; balls accrete [measured]

Co-rotating (charged) θ bath ladder (η=0.5, periodic, T=300, Q_phi = ball-sector charge):

| e_bath | behavior | Q_phi trajectory |
|---|---|---|
| 0.005 | still drains, decelerating | 482 → 314 (vs neutral control → 268) |
| 0.02 | **transient balance** | 482 → 537 (t=50) → 438 (reservoir depletes) |
| 0.05 | strong feeding | 482 → 809 (t=50) → 661 |
| 0.15 | **runaway accretion** | 482 → 1560 (t=50) → ~1400 giant (ω≈1.37) |

- The v67 conclusion is completed: a NEUTRAL bath only erodes; a CHARGED bath
  feeds — the θ→φ charge-return channel is open and strong (CPU preview: sign of
  the flow fully reversed). Detailed balance point ≈ **e* ∈ [0.01, 0.03]**, with
  the caveat that a periodic box is a finite reservoir (all balances are transient;
  a true equilibrium needs a driven/replenished bath or much larger box).
- e=0.15 grows the ball 3× in 50 t.u. — **accretion**: a charged medium does not
  merely sustain balls, it grows them toward the thin-wall regime. Connects to the
  condensate web (FRAGMENT_SPECTRA: fragments still embedded in a charged web that
  feeds them — same physics).

## 3. Mathieu resonance scan — chirped stimulated-emission resonance [measured]

Monochromatic neutral drive (e=0.05, W ∈ {1.20, 1.31, 1.39, 1.47, 1.60}, T=200):

| W | Q_phi(100) | Q_phi(199) |
|---|---|---|
| 1.20 | 329.6 | 278.7 |
| 1.31 | 334.1 | 279.5 |
| **1.39** | **313.6** | 262.5 |
| 1.47 | 325.1 | **256.2** |
| 1.60 | 333.0 | 279.3 |

Resonance detected as **differential drain enhancement (~6–9%)**: at t=100 the
minimum sits at W=ω=1.39 (the n=1 sum-beat tongue); by t=199 the W=1.47 rung has
overtaken it — exactly the **chirp** expected as the draining ball's ω climbs
1.39→~1.44 and the instantaneous resonance sweeps up the ladder. Resonant drive
stimulates emission rather than pumping growth (WP-C's corotating-mode caveat).
The user's "frequency steps" are real, manifesting as emission resonances locked
to the internal clock. [The dramatic-pump channel remains open: requires a
corotating (charged) monochromatic drive — combine cb + md modes.]

## 4. Gauge design (v68/GAUGE_DESIGN.md, 59/59 Maxima) — decision-ready

- **"θ IS the connection" is dead as a reinterpretation** [thm]: the η coupling is
  matter-LINEAR (Proca-like two-point mixing); minimal coupling is matter-QUADRATIC;
  no local field redefinition converts one to the other. The instinct survives as
  in-medium physics: in a charged condensate, minimal coupling linearizes to an
  η_eff ~ g|Φ_bg| two-point mixing — the old phenomenology re-emerges inside matter.
- Recommended package: **new A_μ, gauge the full diagonal U(1) (Φ and Θ co-rotate),
  m_θ=1.6 retained, g=0.05, compact Kogut–Susskind links, temporal gauge.**
  Verified: Gauss law + exact static-ball Coulomb solution; ball exactly
  non-radiating; drain channel closes identically (A is neutral); AB holonomy
  2πqn restored; same-charge 1/r repulsion with contact-force watershed r*≈12;
  g=0 mode bit-identical to v66 kernel.
- Costs/risks: +25% memory; g is a new genuine input; existence-window erosion
  ~+13g²; torus requires net-neutral configurations. Predictions: Q_max≈6200 at
  g=0.05; "positronium" orbit D₁≈15.5 (T_orb≈1050) — in-box.
- **Next gates**: (i) gauged radial shooter (theory, no authorization needed);
  (ii) kernel-v3 implementation — REQUIRES EXPLICIT USER AUTHORIZATION.

## 5. Dressed-state theory (v68/DRESSED_STATE.md) — honest falsification

The isotropic neutral-Gaussian mean-field model was built and solved exactly: right
sign structure (blue→red crossover) but quantitatively dead — blue 75× short, red
unreachable at any drive, and the existence window moves UP (Gaussian smearing melts
balls; the measured sub-window e=0.8 state is impossible in this model class).
Conclusion: **the v67 dressed states are coherent, not thermal** — candidates:
(a) ±ω anti-charge admixture (0.9%/4.4% reproduces the red shifts as APPARENT
frequency; falsifiable via a 2ω beat in s_max, amplitude ≈0.42 at e=0.8); (b)
Mathieu-pumped driven states (growth rates rank-order with the shifts). Bonus
correction: exact in-core Wick slope c(X) flips FRAME's +2.27 to −0.353 (red) at
the core — blue lives only in the ball's skin. Measurement caveat: e≥0.2 rung
extractions carry ±0.02–0.07 method spread; the e=0.05 blue (+0.012) is robust.

## 6. Fragment spectra (v67/results/FRAGMENT_SPECTRA.md) [measured]

- All runs: percolating charged web + droplet spray; discrete balls only at high
  threshold. **No anti-charge fragments anywhere** (all co-rotate with the seed).
- Broad Affleck–Dine initial spectrum; attractor Q≈170 NOT yet reached at T=200–300
  (needs T≳1000); cf4's cores quantitatively match prediction (top-9 mean Q=494.5
  vs predicted 512, −3.5%).
- η accelerates individualization and drains the web but widens (not narrows) the
  early spectrum.

## 7. Next decisions / experiments

1. **Kernel-v3 (gauged U(1))** — authorization decision on GAUGE_DESIGN.md.
   Pre-work available now: gauged radial shooter.
2. **Charged monochromatic drive** (cb×md hybrid) — the true pump channel.
3. **±ω core decomposition** at e=0.8 — discriminate apparent vs driven dressed
   states (DRESSED_STATE §6; the 2ω beat fingerprint).
4. **ℏ_eff=Q phase-tilt analysis** from banked bs03/05/07 SFAs (offline).
5. **Replenished-bath equilibrium** (driven boundary or larger box) to convert the
   transient balance at e*≈0.02 into a true steady state.
6. **Quantized-orbit experiment** (deferred from this round; generator needs
   per-ball velocities — small extension of gen_qball_pair/boost).
