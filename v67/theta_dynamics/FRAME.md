# WP-A — Bulk theta as medium/frame: RELATIVE MASS (not a bare term)

**Date**: 2026-06-10. Basis: v66/THEORY.md §1–2 EOM (12-field complexified Cosserat).
Standard parameters m² = 2.25, m_θ² = 0, η = 0.5, μ = −41.345, κ = 50.
Symbolic verification: `theta_frame.mac` — **28/28 PASS** (`maxima -b theta_frame.mac`,
output `theta_frame.out`). Numbers: `theta_frame.py` (output `theta_frame_py.out`).
Claim marking as in v66/v67: [thm] / [verified] / [estimate] / [open].

## A.1 Linearization around a bulk theta background: the potential is silent

Around Φ = 0 with ANY Θ background, the potential force on Φ_a,

    −∂Vt/∂u_a = −2 Vt′(s) u_a Π_{b≠a}(u_b² + v_b²),

is **O(Φ⁵)**: scaling all 12 Φ-components by ε, the Taylor orders ε⁰…ε⁴ vanish
identically and the leading term is μ u_a Π_{b≠a}(u_b²+v_b²) ε⁵ (Vt′(0) = μ/2 ≠ 0 is
harmless — the quartic prefactor kills it). Θ never appears in Vt(s) at all. So the
**linear Φ mass is exactly m²** — no bare-mass shift from theta at any order in B.
[verified F1a–c]

Consequence [thm]: a theta bulk influences Φ density ONLY through (i) the η-curl
mixing (polariton dressing, §A.2), (ii) the Φ condensate it DRIVES via η∇×Θ, which
then feeds the potential (relative mass, §A.3), and (iii) Floquet/parametric forcing
when the driven condensate is time-modulated (§A.6). This is precisely the
"relative mass, not bare mass" structure conjectured.

## A.2 Vacuum polariton dressing — the B-independent O(η²) self-energy

Transverse plane waves (k ∥ ẑ) mix (Φ_x, Θ_y) pairwise; (u,tu) and (v,tv) are two
identical copies of the V50 real system. Dispersion [verified F2]:

    (ω² − k² − m²)(ω² − k²) = η² k²

- Photon branch: ω² = k²(1 − η²/m²) + O(k⁴) ⇒ v_ph = 0.9428c — **exactly the V50
  polariton speed; the complexification changes nothing** [verified F3]. (Limit check
  requested: PASS. The static-background formula is recovered because at B → 0 the
  fluctuation operator is the vacuum one; the V50 "real background" plays no role at
  linear order — same 2×2 block, two U(1) copies.)
- Matter branch: **ω² = m² + (1 + η²/m²) k² + O(k⁴)** [verified F4]. A Φ quantum's
  inertia is already k-dependent in vacuum: c_eff² = 1 + η²/m² = 1.1111 for the two
  transverse polarizations, c_eff² = 1 for the component-longitudinal one (curl
  decouples it). This O(η²) self-energy is **bath-independent**.

**Key structural fact** [thm, F1+F2]: because the Φ–Θ coupling is strictly LINEAR
(η Re[Θ̄·∇×Φ] is bilinear, a two-point mixing, not a vertex) and the only nonlinearity
is the sextic Φ potential, there is **no O(η²B²) self-energy for ambient Φ quanta**.
The orchestrating conjecture "O(η²B²) scattering off the background" is corrected:
the first bath-dependent ambient self-energy is **O(η⁴B⁴)** (one factor η²B² to build
the driven Φ condensate, one more to scatter off it through the potential). An
O(η²B²)·(core amplitude)⁴ effect DOES exist for the Q-ball itself (§A.5) because the
ball supplies the missing Φ background.

## A.3 Relative mass from the bath: the driven condensate mechanism

Background: circularly-polarized charged Beltrami theta wave (the U(1) phase and the
spatial polarization locked),

    Θ_a = B (x̂ + iŷ)_a e^{i(Kz − Wt)},   ∇×Θ = K Θ,   e_bath = (W² + K²) B²
    [verified F5a–c; W general; the self-consistent background is the photon
     polariton branch W²(K) = (2K² + m² − √(m⁴+4η²K²))/2]

It drives an exact Φ condensate [verified F6a–c]:

    Φ_drv = (ηKB/M̃²) (x̂ + iŷ) e^{i(Kz−Wt)},   M̃² ≡ m² + K² − W²  (detuning)

(at W = K this is the V50 admixture ratio ηK/m²). The potential Hessian on this
background is **static, homogeneous, and component-anisotropic** [verified F7a,b]:

    δm²_ab = μ (ηKB/M̃²)⁴ K̂_a K̂_b        ← THE RELATIVE-MASS FORMULA [thm]

i.e. only the Φ component along the bath wavevector K̂ is shifted; μ < 0 ⇒ the bath
**lowers** the apparent mass. Closed forms [verified F8a–c]:

    δm² = μ η⁴ K⁴ B⁴ / M̃⁸  =  μ η⁴ K⁴ e_bath² / (M̃⁸ (W²+K²)²)
    on bare shell W=K:  δm² = μ η⁴ e_bath² / (4 m⁸)  — K-independent!

Charge-neutral bath (the v67 ladder): a single tu-only Beltrami mode gives HALF the
mean shift plus a parametric drive at (4K, 4W) [verified F12a–c]; an incoherent
isotropic Gaussian neutral bath gives the mean-field value δm² ≈ μ n_comp², with
n_comp the per-component driven variance, isotropic in components [estimate]:

    n_comp = (2/3) η²K² e_bath / ((W²+K²) M̃⁴) = 0.01462 e_bath  (band avg)
    δm²_amb = −0.00884 e_bath²   (Gaussian neutral, isotropic)
    δm²_coh = −0.01991 e_bath²   (coherent charged, K̂-component only) [thm]

Band average over |K| ∈ [1.1, 1.7] on the photon branch; the coefficient varies only
±13% across the band (−0.0225 … −0.0173).

## A.4 Frame, boosts, effective metric

- **The theory has a preferred frame even in vacuum** [verified F11a,b]: the
  polariton polynomial D(ω,k) is Lorentz-invariant iff η = 0; the residual under a
  boost is exactly −η²[k′²_boosted − k′²] (the curl term, F11c). The theta bulk adds a
  SECOND frame (the bath rest frame) and a tensor anisotropy δm²_ab ∝ K̂_aK̂_b.
- **BLV-style effective metric, per polarization** (small k, lattice frame) [thm to
  the stated order]:
    matter branch:  −ω² + (1 + η²/m̄²) k_⊥² + k_∥² + m̄² = 0,
                    g_eff^{00} = −1, g_eff^{ij} = δ^{ij} + (η²/m̄²) P_⊥^{ij},
                    m̄² = m² + δm² (K̂·ê)²   (ê = component polarization)
    photon branch:  ω² = (1 − η²/m̄²) k², so g_eff^{ij} = (1 − η²/m̄²)δ^{ij}.
  Bath correction to the photon speed [verified F9]: δ(v_ph²) = η²δm²/(m²m̄²) < 0 —
  the bath SLOWS the polariton photon: δ(v_ph²) = −2.8×10⁻⁴ at e_bath = 0.8. A
  coherent bath is birefringent for k ⊥ K̂: Δm_eff² between the two transverse
  polarizations = |δm²_coh| = 0.0127 at e_bath = 0.8 [thm/estimate].
- **Relative motion = relative mass** [estimate]: M̃² = m² + (K²−W²) is built from the
  wave-4-vector norm, hence invariant; the entire vrel-dependence enters through the
  lattice-frame bath energy density, e → γ²(1 + v²/3) e for an isotropic bath:
        δm²(v) = δm²(0) · γ⁴ (1 + v²/3)²  ≈  δm²(0) (1 + 8v²/3)
  (×1.03 at v=0.1, ×1.28 at v=0.3, ×2.09 at v=0.5). A Φ object moving through the
  theta bulk gets HEAVIER mass-shift magnitude — a genuinely relativity-like,
  frame-relative effect, but anchored to the bath frame, not observer-reciprocal
  (the medium breaks reciprocity; this is Lorentz violation, not relativity).

## A.5 Q-ball ω_core prediction for the v67 bath ladder (measurable, runs in flight)

Two competing contributions for a ball parked at ω = 1.39 (f(0) = 0.6405):

1. **Core mean-field shift, O(η²B²) — linear in e_bath** [estimate, factor ~2 syst]:
   the bath-driven variance n_comp enters each |Φ_b|² inside the ball:
   δ(ω²) = g′(f₀²)·n_comp, g(X) = μX²/(1+κX³)² (= the v67 condensate relation,
   verified F10a,b). At f₀² = 0.410, κX³ = 3.45 > 1/2 (saturated branch) so
   g′ = +2.270 > 0: **blue shift**, δω_core = +0.01194 e_bath.
   Caveats: in-core detuning is smaller (Hessian W₀ = −0.353 lowers M̃²_core → response
   up to ~40% larger); the bath band straddles the ball's own resonances (by design) —
   resonant enhancement possible [open].
2. **Ambient window shift, O(η⁴B⁴) — quadratic in e_bath** [estimate from the [thm]
   coherent formula + Gaussian statistics]: the ambient m² → m² + δm²_amb shifts the
   whole existence window (ω_min² = m² + (μ/9)(2/κ)^{2/3}) rigidly:
   δω_amb = δm²_amb/(2ω) = −0.00318 e_bath².

**Prediction table (ω_core column of the bath runs; shifts RELATIVE to the e_bath = 0
control at matched Q/time — the η = 0.5 attractor drift ω → 1.44 must be differenced
out):**

| e_bath | n_comp   | δm²_amb   | δm_eff/m (ambient) | δω_core   | δω_amb    | δω_total  | ω_pred |
|--------|----------|-----------|--------------------|-----------|-----------|-----------|--------|
| 0.05   | 0.00073  | −2.2×10⁻⁵ | −4.9×10⁻⁶          | +0.00060  | −0.00001  | **+0.00059** | 1.3906 |
| 0.20   | 0.00292  | −3.5×10⁻⁴ | −7.9×10⁻⁵          | +0.00239  | −0.00013  | **+0.00226** | 1.3923 |
| 0.80   | 0.01170  | −5.7×10⁻³ | −1.3×10⁻³          | +0.00955  | −0.00203  | **+0.00751** | 1.3975 |

Falsifiable structure: δω(e) = +0.0119 e − 0.0032 e² — **linear blue shift at small
e_bath, with measurable downward curvature at the top rung**. A red shift at low
e_bath, or a shift ∝ e² only, falsifies the core-mean-field mechanism; absence of the
curvature term at e = 0.8 falsifies the driven-condensate (relative-mass) term.
Magnitudes (5×10⁻⁴ … 8×10⁻³) need ≳ T=100 of ω-tracking at the e=0.05 rung
(Δω ~ 2π/T resolution → use phase-slope fit, not FFT bins).

## A.6 Pointers to the other two channels (user question ii, iii)

- **Frequency forcing** (iii): the charge-neutral bath's driven condensate is
  time-modulated — H_zz = (μc₀⁴/8)(1 − cos 4(Kz−Wt)) [verified F12b]: a parametric
  drive at (4K, 4W) ⇒ Floquet resonance bands for Φ modes at ω ± 2W, k ± 2K, on top
  of the in-core resonances of §A.5. Step-wise/band structure, not smooth — matches
  the "quantum-like, constructive/destructive" intuition. [thm for single mode; band
  widths O(μc₀⁴/8) — narrow; worked out in WP-C.]
- **Winding** (ii): the Beltrami background has helicity +1 per mode (∇×Θ = +K Θ);
  the sign enters Φ_drv's polarization handedness but NOT δm² (∝ B⁴). Topological
  (winding-number) effects therefore require either chirality-odd observables or the
  κ_h helicity term (excluded from v1) — deferred to WP-B. [open]

## Files

- `v67/theta_dynamics/theta_frame.mac` — Maxima, 28/28 PASS (`theta_frame.out`)
- `v67/theta_dynamics/theta_frame.py` — band numbers + prediction table (`theta_frame_py.out`)
- `v67/theta_dynamics/FRAME.md` — this section
