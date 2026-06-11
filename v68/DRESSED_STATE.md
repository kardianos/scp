# v68 — IN-CORE DRESSED-STATE THEORY (Q-ball in an isotropic theta bath)

**Date**: 2026-06-10. Basis: v66/THEORY.md (12-field complexified Cosserat, radial
Q-ball ODE), v67/FINDINGS.md §4 (bath-ladder measurements), v67/theta_dynamics/FRAME.md
(driven-condensate machinery, F6–F12). Standard parameters m² = 2.25, m_θ² = 0,
η = 0.5, μ = −41.345, κ = 50. Numerics: `v68/theory/dressed.py`
(output `dressed.out`); checks marked PASS/FAIL inside.
Claim marking: [thm] / [verified] / [estimate] / [open].

## 0. What must be explained (v67 measurements)

Matched-Q clock shifts (ω_core at the Q_core = 300 crossing, relative to the e = 0
control at ω = 1.4142), and the deep-bath state:

| e_bath | δω (FINDINGS) | re-extraction (±15 t.u. median) |
|--------|---------------|--------------------------------|
| 0.05   | **+0.012**    | +0.0120 |
| 0.20   | **−0.025**    | −0.0072 |
| 0.80   | **−0.120**    | −0.192 (matched-Q), 1.267 late-time |

plus: at e = 0.8 the ball lives at ω ≈ 1.294 **below the vacuum window edge**
ω_min = 1.3087. The re-extraction column shows the e ≥ 0.2 numbers are
extraction-method sensitive (ω_core is noisy in a strong bath); the e = 0.05 blue
shift is robust. Vacuum-perturbative theory (FRAME A.5) predicted
δω = +0.0119e − 0.0032e²: right blue sign at the bottom rung but ~20× too small,
and no red until e ≈ 3.7 — wrong crossover, no sub-window states. This document
replaces it with an in-core dressed-state (mean-field) theory.

## 1. The mean-field model

The bath is an isotropic, charge-neutral, massless-θ ensemble of energy density e.
Because the η-coupling is bilinear and V is Θ-free [THETA_DYNAMICS thm], the bath
acts on the ball ONLY through the Φ condensate it drives. Each bath mode drives
Φ_drv = (ηKB/M̃²)(…) (FRAME F6); ensemble-summed, the driven Φ field at a point is a
circular Gaussian with per-component variance

    n = n₀ e_bath,   n₀ = 0.01462   (band average K ∈ [1.1,1.7], photon branch)  [estimate]

**In-core detuning enhancement** [estimate]: inside the ball the driven response sees
the ball's own Hessian, M̃² → M̃² + W(r) with W(r) = 2Vt′(s)f⁴ = g(f²),
g(X) ≡ μX²/(1+κX³)², so the local variance is

    n(r) = n₀ e · Enh(f²),   Enh(X) = [M̃²/(M̃² + g(X))]²,   M̃² = 2.452 (band avg)

g < 0 everywhere ⟹ Enh > 1: ×1.36 at the ω=1.39 core (g(0.41) = −0.353, the W₀
flagged in FRAME), peaking at ×2.34 in the skin where f² ≈ 0.215 (g = −0.853, its
global minimum). No resonance is crossed (M̃² + g ≥ 1.60 > 0). [verified numerically]

**The dressed radial problem.** The coherent ball field Φ_a = f(r)e^{iωt} obeys the
Wick (Gaussian mean-field) average of the exact potential force
F_a = 2Vt′(s)Φ_a Π_{b≠a}w_b, w_b = |Φ_b|², giving the modified shooting problem

    f″ + (2/r) f′ = (m² − ω²) f + G_eff(f²; n) f,   f′(0) = 0, f(∞) = 0
    m_eff²(r) − m² = G_eff(f²; n) − g(f²)  →  μn² + O(n³)  as f → 0 (ambient limit)

with G_eff the Gaussian-averaged force/f (computed exactly by 6-D Gauss–Hermite
quadrature, model "M2"; first order in n below, model "M1").

### 1.1 The exact O(n) slope — FRAME A.5's g′ term was incomplete [thm, verified]

To first order in n the Wick average has FOUR pieces, not one:

    ⟨F_a⟩/Φ_a = g(X) + n·c(X),   X = f²,  t = κX³

    c(X) = g′(X)                                  (mean shift ⟨δw_b⟩ = n, all b — the
                                                   only term FRAME A.5 kept; +2.27 at core)
         − 2μtX/(1+t)³                            (δ_a×δw_a cross term, coherent projection)
         + 6μt²X/(1+t)⁴ + 2μtX(2t−4)/(1+t)⁴       (⟨δw_b²⟩ = 2Xn Hessian-trace terms)

    ⟹  **c(X) = 2μX (2t² − 6t + 1)/(1+t)⁴**      [thm; verified against the resummed
                                                   quadrature to ≤1e-2 rel., CHECK A PASS]

Sign structure [thm]: with μ < 0, c > 0 (BLUE) exactly on the band

    t ∈ ((3−√7)/2, (3+√7)/2) = (0.177, 2.823)  ⟺  f² ∈ (0.152, 0.384) ⟺ f ∈ (0.390, 0.620)

and RED outside (both the dilute tail t → 0, c ≈ 2μX, and the deep-saturated core).
Consequences:

- At the ω = 1.39 core, X = 0.410 (t = 3.45): **c = −0.353 (red)** — the FRAME blue
  core term g′ = +2.27 is almost exactly cancelled by the three missing Wick terms.
  The naive "saturated-branch g′ > 0 ⟹ blue core" mechanism is dead. [thm]
- The ball's SKIN (the shell with f between 0.39 and 0.62 — most of the profile
  volume) sits in the blue band; the window-edge amplitude f* = 0.5848 (t = 2,
  c = +1.05) is blue too.
- The net sign of any observable is a profile-weighted competition — resolved
  numerically below.

### 1.2 Matched-Q linear response [verified]

The dressed branch ω(Q; n) from the full shooting solve gives

    dδω/dn |_{Q=300} = **+0.226**  (BLUE)   [verified at n = 5e-4…2e-3, M1 analytic]

i.e. with the physical vacuum-drive variance, δω = +0.0033·e — same sign as FRAME
A.5's +0.0119e but for a different reason (the blue is profile-integrated, skin-band
dominated; the core itself is red), and 3.6× smaller.

Subtlety worth recording [verified]: the energy-envelope route — δE(Q) =
(1/2)∫n C(f²)dV on vacuum profiles with C(X) = ∫₀^X c, then δω = d(δE)/dQ — gives
+0.079/n, a factor 2.9 below the true slope. The two disagree because the
Wick-averaged EOM is NOT the variational derivative of the Wick-averaged energy:
the coherent-projection cross term (δ_a·δw_a) in c(X) has no counterpart in
⟨E⟩. omega_core measures the actual rotation rate, i.e. the eigen-parameter of
the averaged EOM — the full solve, not the energy envelope, is the model
prediction. (C(X) is sign-indefinite — C(0.152) = −0.48, C(0.41) = +0.19 — the
blue survives a near-total tail-red/skin-blue cancellation in both routes.)

## 2. Dressed branch vs measurement — the mean-field FAILS, instructively

Full dressed shooting (resummed M2, in-core enhancement), branch inverted at
Q = 300 [verified; CHECK A/B/C all PASS; vacuum control ω(Q=300) = 1.4057, the
+0.0085 offset to the measured 1.4142 is the η = 0.5 lattice dressing systematic,
differenced out]:

| e_bath | n = n₀e  | ω_min(e) (edge shift) | δω(Q=300) model | δω measured |
|--------|----------|----------------------|-----------------|-------------|
| 0.05   | 0.00073  | 1.30936 (+0.00066)   | **+0.00016**    | +0.012      |
| 0.20   | 0.00292  | 1.31126 (+0.00256)   | **+0.00060**    | −0.025      |
| 0.80   | 0.01170  | 1.31799 (+0.00929)   | **+0.00172**    | −0.120      |

(M1 first-order at e = 0.8: δω = +0.00235 — O(n²) resummation correction is −27%,
so the expansion is under control at these n.)

**Answer to the central question: NO — at the physical drive strength the dressed
branch does NOT shift down.** At matched Q the branch shifts weakly BLUE at every
rung, and the existence window moves UP, not down.

The full δω(n) map over ALL Gaussian variances (decoupling n from the vacuum-drive
value, i.e. allowing arbitrary in-core amplification) is non-monotonic [verified]:

    n      : 0.001    0.01     0.02     0.04     0.05     0.08     0.11     0.12    0.14
    δω     : +0.0002  +0.0016  +0.0021  +0.0009  −0.0006  −0.0117  −0.0240  −0.0265  branch lost
    ω_min  : 1.3096   1.3168   1.3235   1.3351   1.3402   1.3452   1.3432   1.3430  —

- **Blue maximum +0.0021 at n ≈ 0.02**: the measured +0.012 at e = 0.05 is
  unreachable by ANY Gaussian n (best case 6× short; at the physical n = 0.0007,
  75× short).
- **Red regime exists** for n > n* ≈ 0.046 (strong smearing shrinks the core,
  f₀: 0.64 → 0.35, and the matched-Q ball swells to lower ω), reaching −0.025 at
  n ≈ 0.11 — but that is ×39 the vacuum drive of the e = 0.2 rung. **Deepest
  possible red is −0.0265 (n = 0.12); the branch dissolves at n ≥ 0.14. The
  measured −0.120 is unreachable, period.**
- ω_min(n) is monotonically INCREASING everywhere (1.3096 → 1.345 at the red
  extreme). Physically: a charge-neutral Gaussian fluctuation SMEARS the sextic
  well (the Gaussian average of the bounded attractive potential is shallower
  than the potential at the mean), weakening binding — it melts Q-balls, it does
  not deepen them. Same physics as the measured bath-accelerated charge erosion
  (FINDINGS §4), now seen in the existence window. [verified within the model]

## 3. The blue anomaly and the crossover

What the model DOES get right, and its mechanism:

- **The weak-bath shift is blue** — but not for FRAME A.5's reason. The core slope
  is red (c(0.410) = −0.353); the blue is carried by the skin shell
  (f ∈ 0.39–0.62, the c > 0 band), giving dδω/dn = +0.226 at Q = 300 (§1.2).
  [verified]
- **Non-monotonicity IS produced** — blue at weak dressing (skin-envelope term,
  ∝ n), red at strong dressing (well-smearing, onset n* ≈ 0.046 in variance) —
  exactly the measured blue→red pattern, but at the wrong scale: with the physical
  vacuum-drive n(e) = n₀e the predicted crossover is e* = n*/n₀ ≈ **3.2** (the
  e-scan to 1.2 stays blue throughout), vs the measured crossover between 0.05
  and 0.2. Forcing the measured crossover requires the in-core variance to be
  amplified ≥ ×16 at e = 0.2 while ≤ ×60 at e = 0.05 — i.e. a strongly
  SUPERLINEAR n(e), which linear driving cannot give but resonant Mathieu pumping
  (growth +0.017/+0.074/+0.21 per t.u. across the rungs, THETA_DYNAMICS WP-C)
  naturally would. [verified scan + estimate]

So: right sign at e = 0.05, right qualitative blue/red competition, but the blue
magnitude is unreachable (6× short at ANY n), the crossover sits 16× too high in
e, and −0.120 is unreachable. The conclusion is not "tune n harder" — §2 bounds
the entire Gaussian model class — but that the in-core dressing is NOT a
charge-neutral Gaussian. [thm within model class]

## 4. Existence below the vacuum window — requires coherent/charged dressing

The measured e = 0.8 state at ω = 1.294 < ω_min = 1.3087 needs the well DEEPENED
(δω_min² = 1.294² − 1.3087² = −0.038). The Gaussian mean-field strictly cannot do
this (§2: edge moves up at every n until the well is smeared away [verified]).
Three coherent mechanisms can, with required scales [estimate]:

(a) **±ω anti-charged admixture (apparent shift).** A neutral bath feeds equal
    ±ω quanta; the core becomes f₊e^{+iωt} + f₋e^{−iωt}. The kernel's omega_core
    is the instantaneous phase rate θ̇ = (uv̇ − vu̇)/(u² + v²) at the s-max voxel
    (scp_sim.cu pass 3); its median over a window is exactly
    ω(f₊² − f₋²)/(f₊² + f₋²) — RED, while the true branch frequency is unshifted.
    Required contamination: f₋²/f₊² = 0.9% (e = 0.2), **4.4% (e = 0.8)** — tiny,
    entirely plausible for a neutral bath, and consistent with the same admixture
    physics as the rem1 breather (FINDINGS §5). Under this reading the e = 0.8
    "sub-window state" is partly an apparent frequency: the charge-weighted clock
    of a slightly contaminated ball, not a ball off the branch.
(b) **Driven (pumped) state.** With the in-band Mathieu n=1 sum-beat (W + ω = 2ω₀,
    THETA_DYNAMICS WP-C) the bath coherently pumps a co-rotating mode at the
    ball's own frequency: net growth +0.017/+0.074/+0.21 per t.u. at the three
    rungs — which RANK-ORDERS with the measured |red shifts| (≈0/−0.025/−0.120).
    A driven oscillator may sit outside the free spectrum: the existence window
    bounds free dressed eigenstates, not pumped states. The deficit 0.038 in ω²
    must be supplied by the pump term, i.e. a coherent in-core condensate at the
    percent level of the core amplitude. [estimate/open]
(c) **Charge exchange**: bath-mediated Q outflow keeps the interior on a
    falling-Q trajectory whose instantaneous profile lags the branch — transient,
    not an equilibrium dressed state. [open]

Discrimination: (a) predicts the 2ω beat in |Φ|² at the core (the rem1 fingerprint)
with amplitude 2f₋/f₊ ≈ 0.42 at e = 0.8 — directly visible in s_max(t); (b) predicts
the shift collapses when the bath band is moved off the sum-beat resonance
(W ≠ 2ω₀ − ω); (c) predicts the shift heals after the bath is switched off.

## 5. Honest assessment [open items marked]

What the isotropic mean-field misses, and what v68 should measure:

1. **Coherence** (the big one, §4): the model treats the driven in-core Φ as
   circular-Gaussian noise. The Mathieu sum-beat makes it partially COHERENT and
   co-rotating (charged), which evades both no-go results of §2 (a coherent
   co-rotating component adds to the mean field instead of smearing it). The
   resummed model assumed ⟨δΦ⟩ = 0 — exactly the term the pump makes nonzero.
2. **Mathieu pumping** is not in the mean-field at all (it is a Floquet
   instability, secular growth, not a static variance). The rank-correlation of
   WP-C growth rates with the measured shifts is the strongest hint.
3. **Anisotropy**: real bath modes give δm²_ab ∝ K̂_aK̂_b component structure
   (FRAME F7); isotropization averages away birefringent level splittings.
4. **Charge exchange / non-equilibrium**: matched-Q comparison assumes the ball is
   ON its (dressed) branch at the crossing; at e = 0.8 the drain time scale
   (Q: 482→300 by t = 73) is comparable to the relaxation time — the state is
   partly transient.
5. **Measurement systematics**: re-extraction of the v67 ladder (±15 t.u. medians)
   gives δω = +0.0120/−0.0072/−0.19 vs FINDINGS +0.012/−0.025/−0.120 — the e ≥ 0.2
   rungs carry ±0.02–0.07 method spread; only signs and the e = 0.05 magnitude are
   solid targets. Any future fit should use phase-slope ω over fixed Q-windows.

**Discriminating v68 measurement (one run):** rerun the e = 0.8 rung logging the
±ω decomposition of the core (project Φ(s-max voxel) onto e^{±iωt} over a 2π/ω
window — config-level, no kernel change needed if the probe columns are used) and
the s_max 2ω beat amplitude. Outcomes: beat amplitude ≈ 0.4 and f₋²/f₊² ≈ 4% ⟹
mechanism (a) (apparent shift, ball still on-branch — the "sub-window state"
dissolves); no beat but shift collapsing when the bath band moves off W = 2ω₀ − ω
⟹ mechanism (b) (genuine driven dressed state). Either outcome replaces this
mean-field with the correct coherent theory.

## 6. Summary of claims

- c(X) = 2μX(2t²−6t+1)/(1+t)⁴ exact O(n) Gaussian slope; FRAME A.5's g′ incomplete
  (core blue cancelled to red; blue band f ∈ (0.390, 0.620)). [thm, verified]
- Matched-Q linear response dδω/dn = +0.226 at Q = 300 (full solve); the energy
  envelope gives +0.079/n — the averaged EOM is non-variational (the coherent
  cross term has no energy counterpart). [verified]
- Dressed branch (resummed Gaussian + in-core detuning enhancement ×1.4–2.3):
  at physical drive, blue +0.0002…+0.0017 across the ladder; full δω(n) map is
  non-monotonic with blue max +0.0021 (n = 0.02), red onset n* = 0.046, deepest
  red −0.0265 (n = 0.12), branch lost at n = 0.14; window edge RISES at every n
  (max 1.359). **Neutral-Gaussian mean-field falsified as the explanation of the
  v67 ladder**: +0.012 blue unreachable (6× short at any n), −0.120 unreachable,
  sub-window existence impossible, crossover at e* ≈ 3.2 vs measured 0.05–0.2.
  [verified within model]
- Measured pattern requires coherent/charged dressing: ±ω contamination of
  0.9%/4.4% reproduces the red rungs as charge-weighted apparent frequency;
  Mathieu-pumped driven state supplies superlinear amplification (growth ranks
  with shifts) and can sit below the free window. [estimate]
- ω_min(e = 0.8) prediction of this model = 1.318 (vs measured state at 1.294):
  wrong direction — the strongest single falsifier. [verified]

## Files

- `v68/theory/dressed.py` — solver + checks (run: `python3 dressed.py`), output `dressed.out`
- `v68/theory/dressed_red.py` — strong-Gaussian red-regime scan, output `dressed_red.out`
- `v68/DRESSED_STATE.md` — this document
