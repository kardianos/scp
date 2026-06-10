# WP-B — Theta Winding Number: Quantized Impacts (draft)

**Date**: 2026-06-10. Framework: v66/THEORY.md §1–2 (12-field complexified Cosserat),
standard parameters m²=2.25, m_θ²=0, η=0.5, μ=−41.345, κ=50, L=25 (box side 2L=50).
**Verification**: `theta_winding.mac` — **20/20 PASS** (checks W1–W12);
numerics `winding_ball_shift.py` (→ `winding_ball_shift.out`; Q-ball shooter
cross-checked against v66: f(0)=0.6405, Q=482.2, E=691.9 vs published 0.64/482/692).
Markers: [thm] exact, [verified] Maxima/numeric, [estimate], [open].

---

## 0. Structural theorem that controls everything (bilinearity)

The ONLY θ-term in the Φ equation is η(∇×tu)_a (and tv-twin) — LINEAR in Θ; the only
Φ-term in the Θ equation is η(∇×u)_a — LINEAR in Φ; the potential Vt(s) contains no Θ.
Hence the linear fluctuation operator around ANY background (Φ_bg, Θ_bg) contains
**no Θ_bg at all** [thm; verified W12]: Θ_bg enters only (i) as an inhomogeneous
source for Φ, and (ii) indirectly through the Φ-condensate it induces, via the
potential Hessian. Every winding effect on the Φ spectrum is therefore mediated by
the **induced Φ condensate**, not by Θ directly.

## 1. Beltrami backgrounds and winding quantization [verified W1, W2, W10]

- Θ⁻ = B(cos Kz, sin Kz, 0) is an exact Beltrami eigenfield: ∇×Θ⁻ = **−K** Θ⁻;
  the mirror helix Θ⁺ = B(cos Kz, −sin Kz, 0) has ∇×Θ⁺ = **+K** Θ⁺ [verified W1].
- A static Beltrami Θ is NOT a solution at m_θ=0 (∇²Θ = −K²Θ ≠ 0). The dynamical
  realization is the traveling circular wave Θ = B(cos ψ, sin ψ, 0), ψ = Kz − Ωt,
  which solves the free θ-EOM **iff Ω = K** and is Beltrami at every instant
  [verified W2] — exactly the quanta of the v67 bath band.
- Periodic box, side 2L = 50: **K_n = 2πn/(2L) = πn/L = 0.125664·n** — the winding
  number n per box length is the quantized pitch. Bath band |k| ∈ [1.1, 1.7] ⟺
  **n ∈ {9,…,13}** (K₉=1.1310 … K₁₃=1.6336) — five winding sectors [verified W10].

## 2. Φ spectrum in a Beltrami background — chiral splitting verdict

**Linear-order chiral splitting is EXACTLY ZERO.** [thm; verified W4]
For helicity eigenmodes along K (curl → λk, λ=±1) the coupled (Φ,Θ) system gives the
polariton determinant

    (ω² − k² − m²)(ω² − k² − m_θ²) = η² λ² k²,

and λ² = 1: the ± helicity branches are **degenerate**. By §0 the helical Θ_bg cannot
enter the fluctuation operator, so no O(ηB) splitting exists around Φ=0. The expected
"chiral splitting ∝ ηB stepping with n" does NOT occur in vacuum — stated honestly.

What the winding background DOES do, with quantized stepping:

1. **Induced helical Φ condensate** [thm; verified W3]: the forcing η∇×Θ_bg = −ηKΘ_bg
   drives u_ind = c(cos ψ, sin ψ, 0) with

       c = −ηKB/(m² + K² − Ω²)  →(on shell Ω=K)→  c = −ηKB/m².

   At fixed B: **c_n = (ηBπ/Lm²)·n — LINEAR in the winding integer** [verified W11b],
   step Δc = 0.008920·B per unit n; at fixed bath density e_bath = K²B²
   [verified W6], c = η√e_bath/m² = 0.09938 (e_bath=0.2) — **n-independent**
   [thm; verified W11a]. This c ∝ ηB·n quantity is the real "winding charge" imprint.
2. **Polariton branch frequencies step with n** (quantized box modes on the
   helicity-degenerate dispersion): at η=0.5, lower branch ω₋(K₉…K₁₃) =
   1.0700, 1.1898, 1.3097, 1.4299, 1.5503 — step ≈ **0.1195–0.1204 per unit n**
   (→ asymptotic slope = v_phase = √(1−η²/m²) = 0.9428, matching the v50 polariton);
   upper branch ω₊ = 1.914…2.277, step ≈ 0.092 [verified W5 + numerics].
3. **Leading background effect on the spectrum is O((ηB)⁴), scalar, non-chiral**
   [verified W7]: the potential Hessian at the induced condensate (component 2 of
   Φ_bg is zero) has exactly ONE nonzero entry pattern — (u₂,u₂) and (v₂,v₂) =
   μc⁴cos²ψ sin²ψ = (μc⁴/8)(1−cos 4ψ): a uniform mass shift δm²(Φ₂) = μc⁴/8 =
   −5.0×10⁻⁴ (e_bath=0.2) plus a traveling grating at (4K_n, 4Ω_n) → Bragg/parametric
   coupling for δΦ₂ at k = 2K_n (location steps by 2π/L per n), gap width ~|μ|c⁴/8 ~
   5×10⁻⁴ [estimate]. Helicity-blind.

**Route that WOULD give an O(ηB·n) chiral splitting** [estimate/open]: in a uniform
charged Φ-condensate A (v67 §1) the Hessian at background (A+c cos ψ, A+c sin ψ, A)
acquires O(A⁴c) off-diagonal entries rotating with ψ — a chiral grating exactly
analogous to a cholesteric liquid crystal, which Bragg-couples one circular
polarization only. Scale: δω ~ |μ|A³c/(2ω) ~ 0.09 at A=0.4, e_bath=0.2 — large enough
to measure, linear in n at fixed B (c ∝ n). Chirality selectivity not yet
diagonalized [open] — this is the recommended follow-up, but it lives in the
condensate phase, not the vacuum.

## 3. Aharonov–Bohm analog — verdict: NO AB phase from θ winding [thm]

The η coupling is η Re[Θ̄·(∇×Φ)] — an A·(∇×Φ) "helicity/spin-orbit" form, NOT the
A·J_Q minimal coupling. Two independent kills:

- **Linear order**: by §0/W12 a Φ wavepacket's homogeneous propagation is unaffected
  by ANY Θ_bg (vortex Θ ~ Be^{inφ_az} or real-sector azimuthal winding alike); the
  background only adds an inhomogeneous source ⇒ scattering, never a geometric phase.
  Encircling phase = **0 exactly** [thm].
- **All orders in the U(1) phase**: the diagonal-U(1) phase of Φ couples to Θ only
  through the RELATIVE phase in Re[Θ̄·∇×Φ], which is linear in Φ — terms quadratic in
  Φ (the only ones that dress propagation) are Θ_bg-free. Integrating Θ out yields
  η²(∇×Φ̄)G(∇×Φ̄)-type self-couplings with no Θ-winding dependence.

Residual non-topological effect: the induced condensate's δm²(Φ₂) gives a path phase
δφ = −δm²ℓ/(2k) ≈ +9.0×10⁻³ rad over ℓ=50 at k=1.4, e_bath=0.2 — continuous in B,
**not n-quantized** [estimate].

**What WOULD produce a true n-quantized AB phase** [open, design options]:
(a) gauge the diagonal U(1): ∇Φ → (∇ − iqA)Φ with A_μ a new vector field — a vortex
in A gives the standard 2πqn holonomy; (b) add a coupling g·a_θ·J_Q with
a_θ = Σ_a(tu_a∇tv_a − tv_a∇tu_a) (θ's own charge current as composite gauge
potential): a charged-θ vortex Θ_a = Be^{inφ_az} then gives ∮a_θ·dl = 2πnB² —
n-stepped but amplitude-weighted ("quasi-quantized", true quantization only if |Θ|
is modulus-locked). Note the existing conserved current already contains
η(u×tv + tu×v) (v66 §2), so θ DOES transport charge — but as a field bilinear, never
as background×J_Q.

## 4. Q-ball in a helical background [verified W8, W9; numerics]

- **First order: δω = 0 exactly** — δE⁽¹⁾ = −η∫Θ_bg·(∇×u_ball) time-averages to zero
  (ball at ω, background at Ω=K ≠ ω) [thm; verified W8].
- **Second order** (via the induced condensate overlapping the ball): two-angle
  averaging gives δ⟨Vt⟩ = c²g(f), g(f) = Vt′(f⁶)f⁴ + Vt″(f⁶)f¹⁰ =
  (μ/2)f⁴(1−κf⁶)/(1+κf⁶)³ [thm; verified W9a–c]. g flips sign at κf⁶=1 (f=0.521):
  positive in the ball core, negative in the tail; the tail wins. With
  I_g = ∫4πr²g dr from the shot profiles:

  | ω | Q | E | I_g | δE = c²I_g (e_bath=0.2) | δE/E |
  |------|-------|-------|--------|--------|----------|
  | 1.39 | 482.2 | 691.9 | −39.90 | −0.394 | −5.7×10⁻⁴ |
  | 1.41 | 267.4 | 391.5 | −29.29 | −0.289 | −7.4×10⁻⁴ |
  | 1.44 | 139.5 | 209.5 | −20.05 | −0.198 | −9.5×10⁻⁴ |

  Fixed-Q frequency shift δω = d(δE)/dQ ≈ **−4.9×10⁻⁴** (δω/ω ≈ −3.5×10⁻⁴) at
  e_bath = 0.2, η = 0.5: the helical bath DEEPENS binding, ω drifts down [estimate;
  profile integrals verified against v66].
- **Winding stepping of the ball shift**: at fixed e_bath the K-dependence cancels
  exactly (c = η√e_bath/m²) — **δE is n-INDEPENDENT** [thm]. At fixed helix amplitude
  B, δE(n) ∝ n² (quadratic steps): −0.257 / −0.384 / −0.537 for n = 9/11/13
  (B=0.3194). So the ball's frequency is NOT stepwise in n under a constant-density
  bath; stepping appears only in the constant-amplitude normalization, and it is
  quadratic, not linear.

## 5. Answers to the driving question (theta → density, channel ii)

- Winding IS quantized (K_n = πn/L, Beltrami eigenvalue λ_B = ∓K_n) and imprints
  LINEAR-in-n steps on the induced Φ condensate amplitude (at fixed B) and on the
  box polariton ladder (step ≈ 0.12 in ω) [verified].
- It does NOT act as a gauge flux: zero AB phase, no topological protection of any
  Φ phase — the curl coupling is helicity-type, not minimal-coupling-type [thm].
- Spectral chirality from winding requires a Φ condensate (cholesteric mechanism,
  δω ~ 0.09 scale) — the clean falsifiable target for a follow-up run [open].
- Q-ball ω shifts at O(η²e_bath): −3.5×10⁻⁴ relative at e_bath=0.2, n-independent
  at fixed bath density [estimate].

## Files

- `v67/theta_dynamics/theta_winding.mac` — Maxima checks W1–W12 (20/20 PASS)
- `v67/theta_dynamics/winding_ball_shift.py` / `.out` — quantization table, polariton
  ladder, Q-ball shooter + δE/δω, n-stepping tables, AB scattering-phase estimate
