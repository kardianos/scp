# WP-C — Frequency forcing: step-wise / resonant bands (DRAFT)

**Date**: 2026-06-10. Verification: `mathieu_forcing.mac` (17/17 PASS),
`forcing_bands.py` (Floquet 5/5 PASS, output `forcing_bands.out`).
Standard parameters m²=2.25, m_θ²=0, μ=−41.345, κ=50, η=0.5; ball ω=1.39,
core f₀=0.6405 (v66 THEORY §4); bath band |k|∈[1.1,1.7], W=|k| (massless θ).

This formalizes the user's "frequency forcing" intuition: a θ bath does NOT shift
the ball through a bare mass term — it pumps it through **discrete parametric
resonance tongues** (constructive forcing) separated by dead zones (destructive),
i.e. a step-wise, quantized-looking response in drive frequency W.

## 1. Parametric forcing model [thm structure + estimate coefficients]

The bath enters the Φ equation only via η(∇×Θ_bath)_a — an ADDITIVE force at lab
frequency W (the curl coupling is linear, so a θ background cannot directly
modulate the linearized operator). The parametric (Mathieu) channel is two-stage:

1. **Drive**: corotating force amplitude F₀ = η|k|B = ηWB. Bath normalization
   e_bath = 3W²B² per 3 complex θ components [verified C4e] ⟹
   **F₀ = η√(e_bath/3), W-independent** (massless). [thm]
2. **Response**: the driven Φ sidebands beat with the internal rotation e^{iωt};
   |Φ|² (hence s, hence the potential Hessian) is modulated at the beats
   **Ω_b = |W−ω| and W+ω** [verified C4b], with amplitude
   |X| = F₀/√((Ω₀²−Ω_b²)² + (2γΩ_b)²) [verified C4a], where Ω₀ = √(M+4ω²)
   is the k=0 corotating amplitude-mode frequency [verified C4d];
   M = M_sym(f₀) = 4μA⁴(1−2κA⁶)/(1+κA⁶)³ [verified C4c]. At f₀=0.6405:
   M=+1.862, M′=−8.55, Ω₀=3.097. [estimate: condensate Hessian applied at
   ball core amplitude]
3. **Mathieu equation** for any other ball mode q with frequency ω₀(A):

       q̈ + 2γq̇ + [ω₀² + ε cos(Ω_b t)] q = 0,
       **ε(B,η) = |M′(f₀)|·η W B/|Ω₀²−Ω_b²| = |M′| η √(e_bath/3) / |Ω₀²−Ω_b²|**

   First order in η and in bath amplitude. Damping γ ≈ 1.7×10⁻³/t.u. from the
   measured η=0.5 drain (dE/dt=−2.4, E=692 ⟹ amplitude rate ½·2.4/692).
   [estimate]

ω₀ is parameterized (ball normal modes not known analytically). Two candidates:
**ω₀ = ω = 1.39** (internal U(1) rotation / corotating shape mode) and
**ω₀ = 0.94** (the v66 lattice-breathing of s, measured on r1_eta0).

## 2. Resonance tongues — the "quantum steps" [thm, Maxima-verified]

For q̈ + [ω₀² + ε cos(Wt)]q = 0 (standard Mathieu a = 4ω₀²/W², |q_M|=2ε/W²),
instability tongues sit at **W = 2ω₀/n**:

| n | boundaries (Maxima C1–C2) | width in W | max growth | damped threshold |
|---|---|---|---|---|
| 1 | W = 2ω₀ ± ε/(2ω₀) + O(ε²) | ε/ω₀ | ε/(4ω₀) [C3a] | ε_c = 4γω₀ [C3b] |
| 2 | W ∈ [ω₀−5ε²/(24ω₀³), ω₀+ε²/(24ω₀³)] | ε²/(4ω₀³) | ≈ε²/(16ω₀³)–ε²/(8ω₀³)* | ε_c2 = 0.193 (ω₀=1.39, Floquet) |

Mathieu characteristic values verified by harmonic balance: a₁=1+q−q²/8,
b₁=1−q−q²/8, a₂=4+5q²/12, b₂=4−q²/12 [verified C1a–d]. Floquet numerics
reproduce the n=1 center growth to 10⁻⁴, the boundary locations, and the n=2
width (forcing_bands.py, 5/5 PASS). *n=2 max growth measured 9.4×10⁻³ at
ε=0.45, ω₀=1.39 — 2× the ε²/(16ω₀³) rule of thumb.

Discrete tongues + thresholds + dead zones between them = the step-wise
constructive/destructive response. Exponential growth inside a tongue,
nothing outside; tongue count increases step-wise with ε.

## 3. Application to the ball: tongue tables and in-band verdict

### Naive lab-frame tongues W_n = 2ω₀/n (instructed quick-look) [verified arithmetic]

| ω₀ | n=1 | n=2 | n=3 | n=4 | in band [1.1,1.7]? |
|------|------|------|------|------|------|
| 1.39 (U(1) internal) | 2.780 | **1.390** | 0.927 | 0.695 | **n=2 at W=1.39 IN BAND** |
| 1.44 (v66 attractor) | 2.880 | **1.440** | 0.960 | 0.720 | **n=2 at W=1.44 IN BAND** |
| 0.94 (s-breathing) | 1.880 | 0.940 | 0.627 | 0.470 | none (n=1 misses band top by 0.18) |

### Beat-channel refinement [estimate]
The honest parametric drive frequency is the beat Ω_b, so lab resonances are
W = ω ± 2ω₀/n and W = 2ω₀/n − ω. For ω₀=1.39 the ONLY in-band hit is
**W = 1.39 via the sum-beat W+ω = 2ω₀ — an n=1 (strong, O(ε)) tongue**, not
the naive n=2. Both frames single out the same drive: **W* = ω** ("drive theta
at the ball's own internal frequency"). For ω₀=0.94: no in-band hit in either
frame (nearest 1.86). [open: existence of a corotating ball mode at 1.39;
candidate finite-k phonon of the core condensate dispersion at k≈1.9.]

### Verdict for the RUNNING bath rungs (broadband [1.1,1.7], e_bath = 0.05/0.2/0.8)

Monochromatic-equivalent ε (chain of §1, upper bounds): 0.296 / 0.593 / 1.19.
Broadband dilution (only in-tongue power pumps; fixed-point estimate):

| channel | e=0.05 | e=0.2 | e=0.8 |
|---|---|---|---|
| (a) sum-beat n=1 on 1.39-mode | ε_eff=0.105, net growth **+0.017 GROW** | ε_eff=0.42, **+0.074 GROW** | ε_eff=1.19, **+0.21 GROW** |
| (c) naive n=2 (width ε²/4ω₀³) | collapses (below ε_c2) | collapses | collapses |

[estimate] n=2 broadband pumping is negligible at ALL rungs (the self-consistent
in-tongue fraction collapses for ε_mono < √(4ω₀³ΔW_band) = 2.54): the running
broadband baths test channel (a) only. If channel (a) is real (mode exists),
**every rung should show anomalous pumping** — exponential growth of an s-mode
visible as a lab-frequency-2ω ≈ 2.78 rad/t line (period 2.26 t.u.) in s_max(t),
on top of the absorption/drain balance. If no 1.39 corotating mode exists, the
broadband runs show NO parametric signature and the resonance physics is only
reachable by the monochromatic experiment below. Caveat: ε ≳ 0.3 exceeds the
perturbative window — treat growth numbers as saturated upper bounds.

## 4. Direct (non-parametric) forcing: constructive/destructive [thm + arithmetic]

- Driven response (additive channel): Lorentzian |X| = F₀/√((ω₀²−Ω_b²)²+(2γΩ_b)²)
  [verified C4a]. In-band beats |W−ω| ∈ [0, 0.31] reach ONLY soft modes
  (translation/long-wavelength shape sloshing) — constructive maximally at W→ω,
  smoothly (no steps). Off the 0.94 and 1.39 modes everywhere in band: the
  Lorentzian table (forcing_bands.out) is flat to ±12% — direct forcing predicts
  a SMOOTH, weak response across the band. Smooth-vs-stepped is the discriminator
  between direct and parametric channels.
- Two-bath-mode interference: beats |W1−W2| ∈ [0, 0.6] cannot reach either ball
  mode (0.94, 1.39 > 0.6) ⟹ **null prediction**: no intra-bath beat resonance in
  the running band. [verified arithmetic]

## 5. Sharpest experiment: monochromatic drive, n=1 tongue [proposed]

Goal: map dE/dt(W) and the s-spectrum across discrete tongues — flat-top bands of
exponential pumping with sharp edges at W = W* ± ε/(2ω₀) (step-wise signature),
vs the smooth Lorentzian of the direct channel.

- **Seed**: `gen_qball_bath N L profile omega bath_e kmin kmax nmodes rngseed out.sfa`
  with kmin = kmax = W (monochromatic shell, isotropic directions), nmodes=64,
  profile = radial_qball at ω=1.39, N=192, L=25, absorbing BC, T=200, e_bath=0.2.
- **Scan A (in-band / sum-beat n=1, strongest prediction)**: W ∈ {1.25, 1.32,
  1.36, 1.39, 1.42, 1.46, 1.53}. Predicted: pumping inside |W−1.39| ≲ ε/(2ω₀)
  ≈ 0.2·(e=0.2 chain) — growth ~0.1/t.u. (saturating), s_max spectral line at
  ≈2.78 rad/t; outside: ordinary drain. Use attractor ω=1.44 variant if the ball
  has parked: W* shifts to 1.44 — the tongue should TRACK the ball's own ω
  (clean self-consistency test).
- **Scan B (naive-frame n=1, ω₀=1.39)**: W ∈ {2.64, 2.71, 2.78, 2.85, 2.92}
  (k=2.78: 8.7 pts/λ at dx=0.26 — resolved). Predicted tongue half-width
  ε/(2ω₀) ≈ 0.04 at e_bath=0.2 (ε≈0.11 via the Ω_b=|W−ω|=1.39 chain), net
  growth ≈1.8×10⁻²/t.u. → ~3.6 amplitude e-folds in T=200.
- **Scan C (breathing, ω₀=0.94)**: W ∈ {1.80, 1.88, 1.96} — tests the
  s-breathing tongue just above the current band.
- **Controls**: e_bath=0 drain reference; W=2.2 dead-zone point (between all
  predicted tongues — must show NO pumping); one rung check at e_bath=0.05 on
  the W=1.39 peak (threshold crossing: predicted still-positive for channel (a),
  dead if only n=2 operates — discriminates the channel).
- **Observables**: interior E(t), Q(t) (pumping = dE/dt > drain baseline, or
  even >0), s_max(t) FFT (parametric smoking gun: exponential line growth at
  Ω_b, plus period-doubled subharmonic Ω_b/2 in the mode amplitude), r_core(t).

## 6. Files

- `v67/theta_dynamics/mathieu_forcing.mac` — Hill/Mathieu boundaries O(ε²),
  growth, threshold, ε-chain pieces (17/17 PASS; `maxima -b mathieu_forcing.mac`)
- `v67/theta_dynamics/forcing_bands.py` — Floquet validation (5/5), tongue
  tables, rung verdicts, dilution (`forcing_bands.out`)
- `v67/theta_dynamics/FORCING.md` — this draft
