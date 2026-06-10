# THETA DYNAMICS — How the theta sector influences field density

**Date**: 2026-06-10. Framework: v66/THEORY.md §1–2 (12-field complexified Cosserat),
standard parameters m² = 2.25, m_θ² = 0, η = 0.5, μ = −41.345, κ = 50, L = 25.
This document synthesizes workpackages A (frame/relative mass), B (winding), and
C (frequency forcing) and **supersedes the drafts** `FRAME.md`, `WINDING.md`,
`FORCING.md` (retained for provenance). All symbolic claims are machine-checked:
`theta_frame.mac` 28/28 PASS, `theta_winding.mac` 20/20 PASS, `mathieu_forcing.mac`
17/17 PASS (run with `maxima -b <file>.mac`); numerics in `theta_frame.py`,
`winding_ball_shift.py`, `forcing_bands.py` (Floquet validation 5/5 PASS).
Claim markers: [thm] exact, [verified] Maxima/numeric, [estimate], [open].

Cross-references: v66/THEORY.md §5 (theta radiative-drain honesty section) and
FUTURE.md "v67: The Theta-Boundedness Problem" (resolution routes 1–4).

---

## 1. Executive summary — answers to the three-part question

**Driving question**: theta should influence field density NOT via a bare mass term
but as (i) relative mass (bulk theta as medium/frame), (ii) winding number
(topological/quantized), (iii) frequency forcing (step-wise resonant).

**(i) RELATIVE MASS — YES, with a corrected mechanism and an explicit formula.**
The potential Vt(s) contains no Θ and is O(Φ⁵) around Φ=0, so the linear Φ mass is
exactly m² in ANY theta background — there is no bare-mass shift at any order in
the bath amplitude [verified F1]. Theta acts strictly as a MEDIUM: (a) even in
vacuum the η-curl mixing makes the theory preferred-frame (Lorentz-invariant iff
η=0 [verified F11]) with the V50 polariton speed v_ph = √(1−η²/m²) = 0.9428c
recovered exactly in the B→0 limit [verified F3], and matter-branch
c_eff² = 1+η²/m² = 1.1111 transverse; (b) a theta bath shifts apparent masses only
through the Φ condensate it DRIVES: δm²_ab = μ(ηKB/M̃²)⁴ K̂_aK̂_b, M̃² = m²+K²−W²
[thm, verified F7–F8] — an O(η⁴B⁴) tensor (anisotropic, birefringent) shift, NOT
the conjectured O(η²B²) (the η-coupling is a two-point mixing, not a vertex; the
first bath-dependent ambient self-energy needs two insertions). μ<0 ⇒ the bath
LOWERS apparent mass, only the K̂-component. Motion through the bath rescales the
shift, δm²(v) ≈ δm²(0)(1+8v²/3) [estimate] — relativity-LIKE but anchored to the
bath frame (Lorentz violation, not reciprocal relativity). The Q-ball itself DOES
get an O(η²B²) shift (its core supplies the Φ background): net prediction for the
running bath ladder δω_core(e) = +0.0119e − 0.0032e² (blue, with downward
curvature) [estimate].

**(ii) WINDING — quantized, but NOT topological for Φ.** Beltrami helices are
exact curl eigenfields ∇×Θ = ∓KΘ, dynamically realized as traveling circular
waves with Ω = K (= the bath quanta), and the pitch is box-quantized
K_n = πn/L = 0.1257n; the bath band is exactly five winding sectors n ∈ {9..13}
[verified W1,W2,W10]. The hoped-for O(ηB) chiral splitting is EXACTLY ZERO
(helicity branches degenerate, λ²=1 in the polariton determinant [thm W4]), and
there is NO Aharonov–Bohm phase from theta winding at any order — the coupling is
A·(∇×Φ) helicity-type, not A·J_Q minimal coupling [thm]. What DOES step with n:
the induced condensate amplitude c = ηK_nB/m² (linear in n at fixed B) and the box
polariton ladder ω₋(K_n) (step ≈ 0.12) [verified]; the leading spectral imprint is
O((ηB)⁴), non-chiral: δm²(Φ₂) = μc⁴/8 plus a Bragg grating at k = 2K_n [verified
W7]. The Q-ball shift is n-INDEPENDENT at fixed bath density (K cancels) — ball ω
is not stepwise under the constant-e_bath ladder [thm]. A real O(n) chiral
splitting requires a uniform Φ condensate (cholesteric mechanism, δω ~ 0.09
[estimate/open]) or gauging the U(1) (FUTURE.md route 3, which would restore true
2πqn holonomy).

**(iii) FREQUENCY FORCING — YES: Mathieu tongues, with one in-band prediction.**
The bath enters the Φ EOM only as an additive force ηW B at lab frequency W;
parametric pumping is two-stage (driven sidebands beat with e^{iωt} ⇒ Hessian
modulated at Ω_b = |W−ω| and W+ω), giving Mathieu tongues at W_n = 2ω₀/n with
widths, growth rates, and damped thresholds all derived and Floquet-verified
[verified C1–C4]. For the ball modes ω₀ = 1.39 (internal) and 0.94 (breathing),
the ONLY in-band hit of the running experiments is **W* = ω = 1.39 — drive theta
at the ball's own frequency** — and via the sum-beat W+ω = 2ω₀ it is a STRONG
n=1 tongue, not the naive weak n=2. Broadband verdict for the running rungs
(e_bath = 0.05/0.2/0.8): the n=2 channel collapses at ALL rungs; the sum-beat n=1
channel survives with net growth +0.017/+0.074/+0.21 per t.u. IF a corotating
1.39-mode exists [estimate; open] — signature: exponential s_max line at
≈ 2.78 rad/t. If absent: no parametric signature in any rung, and the resonance
physics needs the monochromatic scan (§5, Scan A).

**One unifying structural theorem** [thm, W12/F1]: the η-coupling is bilinear
(η Re[Θ̄·∇×Φ]) and Vt is Θ-free, so the linear fluctuation operator around ANY
background contains no Θ_bg. EVERY theta influence on Φ density is mediated by
the induced Φ condensate (relative mass, winding imprints) or by additive/
parametric forcing (resonance bands) — never by a bare term. This is precisely
the conjectured structure, made exact.

---

## 2. WP-A — Bulk theta as medium/frame: relative mass

Verification: `theta_frame.mac` (28/28 PASS, `theta_frame.out`); numbers
`theta_frame.py` (`theta_frame_py.out`).

### 2.1 Linearization: the potential is silent [verified F1]

Around Φ = 0 with ANY Θ background, the potential force
−2Vt′(s)u_aΠ_{b≠a}(u_b²+v_b²) is **O(Φ⁵)** (Taylor orders ε⁰–ε⁴ vanish;
Vt′(0) = μ/2 is killed by the quartic prefactor); Θ never enters Vt. Linear Φ
mass = m² exactly. Consequence [thm]: theta influences Φ only via (a) η-curl
polariton dressing, (b) the Φ condensate it drives (→ relative mass), (c)
Floquet/parametric forcing — "relative mass, not bare mass," exactly.

### 2.2 Vacuum dressing and the V50 limit [verified F2–F4]

Transverse plane waves mix (Φ_x,Θ_y) pairwise; the complexified system is two
identical copies of V50's real 2×2. Dispersion (ω²−k²−m²)(ω²−k²) = η²k².
Photon branch v_ph = √(1−η²/m²) = 0.9428c — **V50 limit check PASS**. Matter
branch ω² = m² + (1+η²/m²)k²: c_eff² = 1.1111 for the two transverse
polarizations, 1 for the longitudinal (curl decouples it); this O(η²) self-energy
is bath-independent.

**Key correction to the WP premise** [thm]: no O(η²B²) ambient self-energy
exists — the η coupling is two-point mixing, not a vertex. First bath-dependent
ambient self-energy is O(η⁴B⁴): one η²B² to build the driven condensate, one more
to scatter off it through Vt.

### 2.3 The relative-mass formula [verified F5–F8]

Charged Beltrami bath Θ = B(x̂+iŷ)e^{i(Kz−Wt)}, e_bath = (W²+K²)B², drives the
exact condensate Φ_drv = (ηKB/M̃²)(x̂+iŷ)e^{i(Kz−Wt)}, M̃² ≡ m²+K²−W². Its
potential Hessian is static, homogeneous, component-anisotropic:

    δm²_ab = μ (ηKB/M̃²)⁴ K̂_a K̂_b          [thm]
    δm² = μη⁴K⁴e_bath²/(M̃⁸(W²+K²)²);  on bare shell W=K: μη⁴e_bath²/(4m⁸), K-free.

μ<0 ⇒ bath lowers apparent mass, K̂-component only. Neutral bath (the v67
ladder): ×1/2 mean + Floquet drive at (4K,4W) [verified F12]. Band-averaged
(|K|∈[1.1,1.7], photon branch): δm²_amb ≈ μn_comp² = −0.00884 e_bath²
(isotropic Gaussian neutral; coefficient varies ±13% across band);
δm²_coh = −0.0199 e_bath² (coherent charged).

### 2.4 Frame, boosts, effective metric [verified F9, F11]

- Preferred frame even in vacuum: boost residual of D(ω,k) is −η²Δ(k²) ≠ 0 unless
  η = 0. The bath adds a second frame + tensor anisotropy.
- Effective metric per polarization: matter g^00 = −1,
  g^ij = δ^ij + (η²/m̄²)P_⊥^ij with m̄² = m² + δm²(K̂·ê)²; photon
  g^ij = (1−η²/m̄²)δ^ij. Bath slows the polariton photon:
  δ(v_ph²) = η²δm²/(m²m̄²) = −2.8×10⁻⁴ at e_bath = 0.8. Coherent-bath
  birefringence Δm_eff² = 0.0127 (k⊥K̂, e_bath = 0.8).
- Motion through the bath [estimate]: M̃² is boost-invariant; v-dependence enters
  via the lab-frame bath energy density: δm²(v) = δm²(0)γ⁴(1+v²/3)²
  ≈ ×(1+8v²/3) — ×1.28 at v = 0.3. Frame-relative, but anchored to the bath
  frame (not observer-reciprocal): Lorentz violation, not relativity proper.

### 2.5 ω_core prediction for the running bath ladder

Two terms, relative to the e_bath = 0 control at matched Q (the η = 0.5 attractor
drift ω → 1.44 must be differenced out):
1. Core mean-field, O(e): δ(ω²) = g′(f₀²)n_comp, g(X) = μX²/(1+κX³)²
   [verified F10]; at f₀² = 0.410 the saturated branch gives g′ = +2.270 > 0 ⇒
   **blue**, δω_core = +0.0119 e_bath [estimate, factor ~2 syst: in-core detuning
   smaller by W₀ = −0.353; bath band straddles ball resonances → possible
   resonant enhancement [open]].
2. Ambient window, O(e²): δω_amb = δm²_amb/(2ω) = −0.0032 e_bath² (red).

| e_bath | δω_core | δω_amb | δω_total | ω_pred |
|--------|---------|--------|----------|--------|
| 0.05 | +0.00060 | −0.00001 | **+0.00059** | 1.3906 |
| 0.20 | +0.00239 | −0.00013 | **+0.00226** | 1.3923 |
| 0.80 | +0.00955 | −0.00203 | **+0.00751** | 1.3975 |

Falsifiers: red shift at low e or pure-e² scaling kills the core-mean-field term;
no downward curvature at e = 0.8 kills the driven-condensate relative-mass term.
Needs phase-slope ω tracking (not FFT bins) for the 6×10⁻⁴ shift at e = 0.05.

---

## 3. WP-B — Winding number: quantized imprints, no topology for Φ

Verification: `theta_winding.mac` (20/20 PASS); `winding_ball_shift.py` /
`winding_ball_shift.out` (shooter cross-validated vs v66: f(0) = 0.6405,
Q = 482.2, E = 691.9 vs published 0.64/482/692).

### 3.1 Quantization [verified W1, W2, W10]

Helix B(cos Kz, sin Kz, 0) is exact Beltrami, ∇×Θ = −KΘ (mirror → +K). Static
Beltrami is NOT a solution at m_θ = 0; the dynamical realization is the traveling
circular wave ψ = Kz−Ωt with **Ω = K exactly** — i.e. the bath quanta themselves.
Box quantization K_n = πn/L = 0.1257n (2L = 50); bath band [1.1,1.7] ⟺
n ∈ {9,…,13}: five winding sectors.

### 3.2 Chiral-splitting verdict: zero at linear order [thm W4]

Helicity eigenmodes give (ω²−k²−m²)(ω²−k²−m_θ²) = η²λ²k² with λ² = 1 — the ±
branches are degenerate; by the structural theorem (§1) Θ_bg cannot enter the
fluctuation operator, so no O(ηB) splitting exists in vacuum. What steps with n:
- Induced condensate c = ηK_nB/m²: **linear in n at fixed B** (Δc = 0.00892·B per
  unit n) but n-independent at fixed e_bath = K²B² (c = η√e_bath/m² = 0.0994 at
  e_bath = 0.2) [verified W3, W6, W11].
- Box polariton ladder ω₋(K₉…K₁₃) = 1.070/1.190/1.310/1.430/1.550, step
  0.1195–0.1204 (slope → V50 polariton 0.9428) [verified W5].
- Leading spectral effect O((ηB)⁴), non-chiral [verified W7]: Hessian at the
  induced condensate has only (u₂,u₂),(v₂,v₂) = μc⁴cos²ψ sin²ψ ⇒
  δm²(Φ₂) = μc⁴/8 = −5.0×10⁻⁴ (e_bath = 0.2) + Bragg grating at k = 2K_n
  (position steps with n), gap ~5×10⁻⁴.

Route to a real O(ηBn) chiral splitting [estimate/open]: uniform charged
condensate A + helix ⇒ O(A⁴c) rotating Hessian — cholesteric mechanism,
Bragg-coupling one circular polarization; scale δω ~ |μ|A³c/(2ω) ~ 0.09 at
A = 0.4. Measurable; recommended follow-up (lives in the condensate phase).

### 3.3 Aharonov–Bohm verdict: NO [thm]

The coupling is A·(∇×Φ) helicity-type, not A·J_Q; the phase of Θ enters only
linearly in Φ, so all Φ-quadratic (propagation-dressing) terms are Θ_bg-free.
Encircling phase = 0 exactly, at any order, for any winding. Residual:
non-topological scattering phase ≈ 9.0×10⁻³ rad over ℓ = 50 (continuous in B,
not n-quantized). What WOULD give AB: (a) gauging the diagonal U(1) (standard
2πqn holonomy — FUTURE.md boundedness route 3); (b) coupling g·a_θ·J_Q with
a_θ = Σ(tu∇tv−tv∇tu): a charged-θ vortex gives ∮a_θ·dl = 2πnB² —
quasi-quantized (amplitude-weighted).

### 3.4 Q-ball in a helical bath [verified W8, W9]

First order: δω = 0 exactly (time-average, ω ≠ Ω) [thm W8]. Second order:
δ⟨Vt⟩ = c²g(f), g = (μ/2)f⁴(1−κf⁶)/(1+κf⁶)³ — sign flip at f = 0.521 (core +,
tail −, tail wins). I_g = −39.9/−29.3/−20.1 at ω = 1.39/1.41/1.44;
δE = −0.394 at ω = 1.39, e_bath = 0.2 ⇒ **δω(fixed Q) ≈ −4.9×10⁻⁴**
(δω/ω ≈ −3.5×10⁻⁴), binding deepens. Stepping: n-INDEPENDENT at fixed e_bath
[thm, K cancels]; ∝ n² at fixed B (−0.257/−0.384/−0.537 for n = 9/11/13).
The ball's ω is NOT stepwise under the constant-density bath ladder.

Note: this coherent-helix ball shift (red, −4.9×10⁻⁴) and the WP-A incoherent-bath
core shift (blue, +2.4×10⁻³ at e = 0.2) are different bath statistics/geometries —
the running isotropic-bath ladder tests WP-A's number; a single-helix seed would
test this one.

---

## 4. WP-C — Frequency forcing: Mathieu tongues, in-band prediction

Verification: `mathieu_forcing.mac` (17/17 PASS); `forcing_bands.py` /
`forcing_bands.out` (Floquet 5/5 PASS).

### 4.1 Model [thm structure / estimate coefficients]

Bath enters the Φ EOM only as additive force η(∇×Θ_bg), F₀ = ηWB = η√(e_bath/3)
(W-independent, massless θ; e_bath = 3W²B² [verified C4e]). Parametric channel is
two-stage: driven Φ sidebands beat with e^{iωt} ⇒ s (hence the Hessian) modulated
at Ω_b = |W−ω| and W+ω. Mathieu equation:

    q̈ + 2γq̇ + [ω₀² + ε cos(Ω_b t)]q = 0,
    ε = |M′(f₀)|·η√(e_bath/3)/|Ω₀²−Ω_b²|

with M_sym(0.6405) = +1.862, M′ = −8.55, Ω₀ = √(M+4ω²) = 3.097;
γ = 1.7×10⁻³ from the measured η-drain (½·2.4/692 — the drain of v66/THEORY §5).

### 4.2 Tongues [Maxima-verified C1–C3]

W_n = 2ω₀/n. n=1: boundaries 2ω₀ ± ε/(2ω₀), width ε/ω₀, growth ε/(4ω₀),
threshold ε_c = 4γω₀ = 9.6×10⁻³. n=2: width ε²/(4ω₀³), Floquet threshold
ε_c2 = 0.193 (ω₀ = 1.39). Mathieu characteristics a₁ = 1+q−q²/8, b₁ = 1−q−q²/8,
a₂ = 4+5q²/12, b₂ = 4−q²/12 by harmonic balance; Floquet integrator reproduces
all to ≤10⁻⁴. Discrete tongues + thresholds + dead zones = the step-wise
constructive/destructive response.

### 4.3 In-band verdict for the running rungs

Naive lab tongues: ω₀ = 1.39 ⇒ {2.780, **1.390 (n=2, in band)**, 0.927, 0.695};
ω₀ = 1.44 ⇒ 1.440 in band; ω₀ = 0.94 (breathing) ⇒ none (1.880 misses band top
by 0.18); intra-bath beats |W₁−W₂| ≤ 0.6 reach neither mode [verified].
**Beat-frame refinement** [estimate]: the true drive is Ω_b; the only in-band lab
hit is W = 1.39 via the SUM-beat W+ω = 2ω₀ — an **n=1 (O(ε)) tongue**, stronger
than the naive n=2. Both frames flag **W* = ω**. Requires a corotating 1.39-mode
(candidate: core-condensate phonon at k ≈ 1.9) [open].

Broadband rungs (e_bath = 0.05/0.2/0.8): monochromatic-equivalent
ε = 0.296/0.593/1.19 (upper bounds; ε ≳ 0.3 saturated). In-tongue dilution kills
the n=2 channel at ALL rungs (fixed point dies for ε_mono < √(4ω₀³ΔW) = 2.54).
Sum-beat n=1 survives: ε_eff = 0.105/0.42/1.19, net growth +0.017/+0.074/+0.21
per t.u. **If the 1.39 corotating mode exists, all three rungs pump**
(exponential s_max line at ≈2.78 rad/t, period 2.26 t.u.); if absent, no
parametric signature in any rung.

### 4.4 Direct (non-parametric) channel — the discriminator

Driven Lorentzian response is flat to ±12% across the band: direct forcing
predicts a SMOOTH weak response; parametric predicts stepped exponential tongues
with thresholds. Smooth-vs-stepped discriminates the channels [thm + arithmetic].

---

## 5. Unified experiment proposals (ranked by cost)

All are config/seed-level — **no kernel changes** (per CLAUDE.md kernel policy)
— except where flagged.

1. **FREE (already running)** — v67 bath ladder rungs e_bath = 0.05/0.2/0.8 +
   e = 0 control. Analysis additions only:
   (a) phase-slope ω_core(t) tracking (not FFT bins) → test
       δω(e) = +0.0119e − 0.0032e² (WP-A table §2.5);
   (b) s_max(t) spectrogram → look for exponential line at ≈ 2.78 rad/t
       (WP-C sum-beat n=1; presence at ALL rungs = corotating mode exists;
       absence everywhere = mode absent, monochromatic scan required);
   (c) anisotropy of ⟨|Φ_a|²⟩ per component — coherent-bath remnant
       birefringence (WP-A §2.4) [likely below noise for isotropic bath].
2. **CHEAP (config-only, existing gen_qball_bath)** — **monochromatic tongue
   scan**, N = 192, L = 25, profile = radial_qball ω = 1.39, e_bath = 0.2,
   nmodes = 64, kmin = kmax = W, T = 200, absorbing BC:
   - Scan A (strongest): W ∈ {1.25, 1.32, 1.36, 1.39, 1.42, 1.46, 1.53} —
     pumping band centered W* = ω; tracks to 1.44 if ball parked
     (self-consistency).
   - Scan B: W ∈ {2.64, 2.71, 2.78, 2.85, 2.92} (naive n=1, half-width ≈ 0.04,
     ~3.6 e-folds in T = 200; k = 2.78 resolved at 8.7 pts/λ).
   - Scan C: W ∈ {1.80, 1.88, 1.96} (breathing tongue just above band).
   - Controls: e_bath = 0; W = 2.2 dead zone (must NOT pump); e_bath = 0.05 at
     W = 1.39 (threshold discriminator n=1 vs n=2).
   - Observables: interior E(t), Q(t), s_max(t) FFT (exponential line at Ω_b +
     subharmonic Ω_b/2 = parametric smoking gun), r_core(t).
3. **CHEAP (config-only)** — **single coherent helix** seed (one mode,
   kmin = kmax = K_n, one direction) at fixed e_bath = 0.2, n ∈ {9, 11, 13}:
   tests WP-B ball shift δω ≈ −4.9×10⁻⁴ and its n-independence at fixed
   density; repeat at fixed B → δE ∝ n² check (−0.257/−0.384/−0.537).
4. **MODERATE (new seed generator in sfa/seed/, no kernel change)** —
   **cholesteric chiral splitting**: uniform charged Φ condensate A ≈ 0.4 +
   helix bath (e_bath = 0.2); measure circular-polarization-resolved Bragg gap
   at k = 2K_n; predicted δω ~ 0.09, linear in n at fixed B. The clean
   winding-quantization observable.
5. **EXPENSIVE (kernel physics change — requires explicit user authorization)** —
   gauge the diagonal U(1) (FUTURE.md boundedness route 3): restores a true
   2πqn AB holonomy AND addresses theta-boundedness; design doc first. The
   composite a_θ·J_Q coupling (quasi-quantized 2πnB²) is an intermediate
   theory-only study.

## 6. Honest limits — what was linearized away, and where the gauge route changes answers

- **All Φ-spectrum statements are linear-response around Φ = 0 or around the
  driven condensate/ball mean field.** The sextic self-interaction of ambient
  fluctuations, condensate depletion, and rung-3 nonlinearity (ε ≈ 1.19,
  beyond the perturbative window — growth numbers are saturated upper bounds)
  are not controlled.
- **Bath statistics idealized**: relative-mass coefficients assume either a
  single coherent Beltrami mode [thm] or Gaussian isotropic neutral statistics
  [estimate]; the real gen_qball_bath superposition (64 modes, band-limited) sits
  between (−0.0088 vs −0.0199 per e²; ±13% K-dependence across the band).
- **Core coefficient systematics**: δω_core has factor ~2 uncertainty (in-core
  detuning M̃²_core lowered by the Hessian W₀ = −0.353; bath band straddles ball
  resonances — possible resonant enhancement [open]). The WP-C ε-chain applies a
  condensate Hessian at the ball core amplitude [estimate].
- **Ball mode spectrum unknown**: ω₀ = 1.39 (corotating) and 0.94 (breathing)
  are candidates, not a computed normal-mode spectrum; the central WP-C
  prediction is conditional on the 1.39 corotating mode existing [open].
- **Radiative drain backdrop** (v66/THEORY §5): the massless-θ channel is
  radiative — every prediction here sits on top of a secular ω-drift toward the
  1.44 attractor; all shifts are defined relative to the e_bath = 0 control at
  matched Q. γ in the Mathieu model is this drain, treated as linear damping.
- **m_θ = 0 assumed throughout**: the boundedness routes (FUTURE.md v67 §) alter
  the answers — m_θ > ω (route 1, runs th1/th2) replaces the radiative response
  with Yukawa dressing (M̃² → m² + K² − W² + m_θ² grows; relative-mass shifts
  shrink; the polariton photon branch gaps); Cosserat binding (route 2) gaps θ
  at √(2α) and similarly detunes everything.
- **The gauge route (experiment 5) qualitatively changes WP-B**: minimal
  coupling converts the null AB verdict into a 2πqn holonomy, makes winding
  genuinely topological for Φ, and replaces the "no chirality at O(ηB)" theorem
  (the fluctuation operator then contains A_bg directly). WP-A's "no bare mass"
  theorem survives (gauge fields don't enter Vt), but the medium dispersion
  acquires the standard A²|Φ|² seagull — an honest O(A²) mass-like term. The
  no-go results here are statements about the CURRENT η-curl coupling only.

## Files

- `theta_frame.mac` / `theta_frame.out` — WP-A Maxima, 28/28 PASS
- `theta_frame.py` / `theta_frame_py.out` — WP-A band numbers, prediction table
- `theta_winding.mac` — WP-B Maxima, 20/20 PASS
- `winding_ball_shift.py` / `winding_ball_shift.out` — WP-B quantization ladder,
  shooter, δE/δω, n-stepping, AB phase estimate
- `mathieu_forcing.mac` — WP-C Mathieu/Hill boundaries, growth, thresholds,
  ε-chain, 17/17 PASS
- `forcing_bands.py` / `forcing_bands.out` — WP-C Floquet validation (5/5),
  tongue tables, rung verdicts
- `FRAME.md`, `WINDING.md`, `FORCING.md` — superseded section drafts (provenance)
