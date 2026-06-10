# v67 THEORY — Charged condensate, modulational instability, and bath balance

**Date**: 2026-06-10. Builds on v66 (complexified 12-field kernel, verified; `v66/THEORY.md`).
Standard parameters: m² = 2.25, m_θ² = 0, η = 0.5, μ = −41.345, κ = 50;
Vt(s) = (μ/2)s/(1+κs), s = Π_a(u_a²+v_a²) = Π_a|Φ_a|².

**Claim marking**: [thm] exact derivation, [verified] checked symbolically
(`theory/mi_dispersion.mac`, 20/20 PASS) and/or numerically (`theory/mi_band.py`,
`theory/mi_eta_scan.py`), [estimate] order-of-magnitude with stated assumptions,
[open] unresolved.

---

## 1. Uniform charged condensate

Ansatz: Φ_a = A e^{iω_c t} for all a (u_a = A cos ω_c t, v_a = A sin ω_c t), Θ = 0,
spatially uniform. Then s = A⁶ is constant, ∇×Φ = 0 (so Θ = 0 is EXACT here, unlike
the Q-ball, where ∇×Φ ≠ 0 sources Θ radiation — the condensate has no curl source). [thm]

**Self-consistent dispersion** [thm; verified D1, D2]:

    ω_c² = m² + 2 Vt′(A⁶) A⁴ = m² + μ A⁴ / (1 + κA⁶)²

(the factor 2 is the real-field convention: ∂Vt/∂u_a = 2Vt′(s) u_a Π_{b≠a}|Φ_b|²).
Since μ < 0, Vt′ < 0, so **ω_c < m for every A**: the condensate is always bound
relative to free quanta. [thm]

**Charge and energy densities** [thm; verified D3, D4]:

    ρ_Q = 3 ω_c A²
    e_0 = (3/2)(ω_c² + m²) A² + Vt(A⁶)

**Elegant special point** [thm; verified D9]: at A = f* = (2/κ)^{1/6} = 0.5848 (the
existence-window amplitude κf⁶ = 2 of v66 §3.1), ω_c² = m² + (μ/9)(2/κ)^{2/3} = ω_min²
EXACTLY — the uniform condensate at the window amplitude oscillates at the bottom edge
of the Q-ball window.

| A | κA⁶ | ω_c | ρ_Q | e_0 |
|------|--------|--------|--------|--------|
| 0.30 | 0.036 | 1.3922 | 0.3759 | 0.5509 |
| 0.40 | 0.205 | 1.2332 | 0.5919 | 0.8347 |
| 0.50 | 0.781 | 1.1982 | 0.8986 | 1.2008 |
| 0.585| 2.004 | 1.3090 | 1.3439 | 1.7588 |

## 2. Modulational instability (Affleck–Dine fragmentation)

This is the standard rotating-condensate fragmentation mechanism of Affleck–Dine
baryogenesis (Kusenko & Shaposhnikov, PLB 418 (1998) 46; Kasuya & Kawasaki, PRD 62
(2000) 023512; Enqvist & McDonald, PLB 425 (1998) 309): a U(1) condensate in a
shallower-than-quadratic potential is unstable to wavelength-scale clumping and
fragments into Q-balls. Our three-field generalization follows.

### 2.1 Linearization [thm; verified D5–D8 + 48×48 numeric cross-check to 6e-17]

Perturb Φ_a = e^{iω_c t}(A + χ_a), χ_a = (R_a + iI_a) e^{ik·x + σt}. In the corotating
frame the U(1)-invariant potential Hessian is time-independent; evaluated at the gauge
point u = (A,A,A), v = 0 [verified D5]:

    H_uu = W₀·𝟙 + C,  H_vv = W₀·𝟙,  H_uv = 0,   W₀ = 2Vt′(A⁶)A⁴
    C_ab = 4A⁴ [ A⁶ Vt″(A⁶) + Vt′(A⁶)(1 − δ_ab) ]

W₀ cancels against the dispersion relation: the I-sector (phase) is massless — the
Goldstone of the diagonal U(1) — and the R-sector (amplitude) carries C. Per eigenmode
M of C, the (R, I) pair obeys [verified D8]

    σ⁴ + σ²(2k² + M + 4ω_c²) + k²(k² + M) = 0

**Growth (σ² > 0) exists iff M < 0, with instability band 0 < k² < −M.** [thm]

Eigenvalues of C [verified D6, D7]:

    M_sym  = 8A⁴Vt′ + 12A¹⁰Vt″ = 4A⁴ μ (1 − 2κA⁶)/(1 + κA⁶)³     (mode R ∝ (1,1,1))
    M_anti = −4A⁴Vt′ = −2μA⁴/(1 + κA⁶)²  > 0 always (μ<0)         (two modes, Σ R_a = 0)

The antisymmetric (relative-amplitude) channel is ALWAYS stable. The symmetric channel
is unstable iff μ(1 − 2κA⁶) < 0, i.e. (μ < 0):

    **instability threshold:  κA⁶ < 1/2   ⟺   A < A* = (2κ)^{−1/6} = 0.46416**  [thm]

Physically: below A* the potential is shallower than quadratic (attractive, negative
pressure); above A* the sextic attraction has saturated (Vt → μ/2κ flat) and the
condensate is phi-sector STABLE. Growth rate:

    σ²(k) = ½[ −(2k²+M+4ω_c²) + √((2k²+M+4ω_c²)² − 4k²(k²+M)) ]
    σ_max ≈ |M|/(4ω_c) at k² ≈ |M|/2     (small-|M| approx; exact σ_max matches to 4 digits)

### 2.2 Full 12-field verification (theta included) [verified, `mi_band.py`]

The corotating frame also rotates the (tu, tv) doublet (the η-curl coupling commutes
with the diagonal U(1)), making the FULL linear system autonomous; with cos/sin spatial
components it is a 48×48 first-order eigenproblem, scanned over k and direction k̂
(the background picks (1,1,1) in component≡spatial space, so k̂ matters at η≠0).

- At η = 0 the 48×48 growth reproduces the analytic σ(k) to 6×10⁻¹⁷, and theta modes
  are neutrally stable (σ = ±ik). [verified]
- At A < A*, η ≤ 0.75: the dominant mode is unchanged — for k̂ ∥ (1,1,1) the curl
  decouples from the symmetric eigenvector (k̂×R = 0) and σ_max = 0.28975 (A=0.40)
  exactly as at η = 0; other directions are reduced ≤ 8%. η does NOT slow
  fragmentation below threshold. [verified]
- **NEW: above threshold (A > A*) the condensate is NOT absolutely stable at η ≠ 0.**
  A curl-mediated parametric instability appears in a narrow band k* ≈ 0.8–0.9 with
  growth FIRST ORDER in η (local exponent 0.99 at small η — resonant pumping of theta
  by the rotating condensate, distinct from the η² Q-ball drain):

      A = 0.50:  σ* = 0.0163 (η=0.125), 0.0626 (η=0.5), k* ≈ 0.80
      A = 0.585: σ* = 0.0174 (η=0.125), 0.0671 (η=0.5, worst k̂=(1,−1,0)), k* ≈ 0.87

  [verified numerically, `mi_eta_scan.py`; analytic form of this branch — open]

### 2.3 Headline testable predictions (N=192, L=25 box, side 2L=50, periodic BC)

Phi-sector MI (A < A*), η-independent rates; λ = 2π/k_max; N_frag ~ (2L/λ)³;
Q_frag = ρ_Q λ³:

| A | M_sym | band k < | k_max | σ_max | λ | N_frag | Q_box | Q/frag |
|------|---------|-------|--------|--------|-------|------|--------|------|
| 0.30 | −1.1154 | 1.056 | 0.7194 | 0.2003 | 8.73 | ~188 | 46987 | ~250 |
| 0.35 | −1.5559 | 1.247 | 0.8309 | 0.2958 | 7.56 | ~289 | 60415 | ~209 |
| **0.40** | −1.4293 | 1.196 | **0.7942** | **0.2898** | **7.91** | **~252** | 73993 | **~293** |
| 0.45 | −0.5365 | 0.637 | 0.4423 | 0.0856 | 14.2 | ~44 | 89962 | ~2060 |
| 0.50 | +1.029 | — | stable (phi); η-mode k*=0.802, σ*=0.0626 | | 7.83 | ~260 | 112325 | ~430 |
| 0.585| +2.149 | — | stable (phi); η-mode k*=0.866, σ*=0.0671 | | 7.25 | ~328 | 167988 | ~512 |

Fragmentation time from seed noise δ (relative): t_frag ≈ ln(0.1A/δ)/σ_max.
With 1% amplitude noise: **t_frag ≈ 16 t.u. at A=0.40** (visible clumping by t≈16,
formed lumps by t≈25–30); ~5/σ ≈ 17 for full nonlinear saturation. For the η-mode at
A=0.50, η=0.5: t_frag ≈ 74 t.u. — an order of magnitude slower.

**Falsifiable signatures** [predictions]:
1. **Sharp threshold at A* = 0.464**: at η = 0, A = 0.40 fragments in ~tens of t.u.;
   A = 0.50 and 0.585 do NOT fragment at all (any fragmentation there at η=0 falsifies
   the analysis). At η = 0.5 the above-threshold cases fragment ~10× slower, and the
   rate scales ∝ η (vs η-independent below threshold) — three clean discriminators.
2. **Fragment spacing λ ≈ 7.9 at A = 0.40** (structure-factor peak at k ≈ 0.79;
   lattice resolves it: dx = 0.26 → 30 pts/λ; dk = 2π/50 = 0.126 → ~9 band modes/axis).
3. **~250 fragments** in the 50³ box at A = 0.40, mean Q/fragment ≈ 290 — above
   Q_min ≈ 92 and the η=0.5 attractor Q∞ ≈ 170 (v66): fragments are born viable, then
   at η = 0.5 drain along the branch toward the Q∞ ≈ 170 attractor — a natural-formation
   route to the SAME particle, completing the Affleck–Dine analogy. [estimate: counts
   assume one lump per λ³ cell; factor ~2 from nonlinear mergers expected]
4. Weak anisotropy at η = 0.5, A < A*: fastest growth along k̂ ∥ (1,1,1) (σ 0.290 vs
   0.265–0.271 transverse) — a few-percent structure-factor anisotropy. [verified linear;
   visibility in 3D — open]
5. ω_c(0.40) = 1.233 < ω_min = 1.3087: fragments must migrate INTO the Q-ball window —
   amplitude rises locally during clumping. [thm + open: nonlinear stage]

## 3. Bath balance (Q-ball in an isotropic massless theta bath)

Measured drain at η = 0.5 (v66 FINDINGS §4, ω = 1.39, Q = 482, E = 692, r_Q = 3.71):
dQ/dt = −1.72, dE/dt ≈ −2.4. Consistency: (dE/dt)/(dQ/dt) = 1.395 ≈ ω_core = 1.39 —
the radiation carries E/Q = ω quanta, as massless theta at k = ω must. [verified
consistency]. Quality factor ωE/P = 400: weak coupling. ωr_Q = 5.2: wavelength-scale ball.

Absorption from an isotropic bath of energy density e_bath (c = 1):
P_abs = σ_abs·e_bath; **balance at e_bath* = |dE/dt|/σ_abs** [thm given σ_abs].
Cross-section estimates [estimate]:

- Geometric upper bound: σ = πr_half² = 61  → e_bath* ≥ 0.039.
- Partial-wave unitarity (ℓ ≤ 1 source, v66 §5, 2 transverse polarizations):
  σ = 2π(2ℓ+1)/ω²|_{ℓ≤1} = 13  → e_bath* ≈ 0.18.
- Mode-occupation criterion (stimulated ≈ spontaneous when n̄ ~ 1 at the emission
  frequency): e* ~ g ω³Δω/2π², g = 2–6 polarizations, Δω = (0.1–1)ω → 0.04–1.1.
- Scaling [estimate, by reciprocity with the measured η^1.88 emission law]:
  σ_abs ∝ η² R² at weak coupling, saturating toward geometric at η ≳ 0.75.

**Combined bracket: e_bath* ∈ [0.04, 1.1], center ≈ 0.2** (vs ball core e(0) = 2.25,
i.e. a balancing bath is ~2–50% of core density). [estimate]

**Recommended bath ladder: e_bath ∈ {0.05, 0.2, 0.8}** (+ e_bath = 0 control),
isotropic theta noise band-centered at ω ≈ 1.39 (width ±0.3), equal (tu, tv) power
(charge-neutral). dE/dt(e_bath) should cross zero inside the ladder; its slope near
zero measures σ_abs directly.

**Charge caveat** [estimate/open]: a charge-neutral bath restores E but not net Q
(absorbed quanta average to zero charge), so at e_bath* expect energy stasis with a
residual Q drift and a SHIFT of the attractor Q∞(e_bath) — not full stasis. Exact
detailed balance of both E and Q needs a circularly-polarized (co-rotating, charged)
bath whose quanta carry E/Q = ω. Predicted discriminator: neutral bath ⇒ E/Q of the
parked ball falls below 1.497; charged bath ⇒ stays on the v66 attractor line.
**Beware: e_bath ≥ 0.04 is within reach of the condensate MI scale (§2) — at the top
rung watch for bath-seeded nucleation of new balls [open].**

## 4. Verification summary

- `theory/mi_dispersion.mac`: 20/20 PASS (dispersion, ρ_Q, e_0, Hessian blocks,
  C eigenvalues, closed forms, threshold κA⁶ = 1/2, quartic determinant, ω_c(f*) = ω_min).
- `theory/mi_band.py`: band table above; 48×48 full-field growth equals analytic
  σ(k_max) to 5.6×10⁻¹⁷ at η = 0; η-mode discovery above threshold.
- `theory/mi_eta_scan.py`: η-exponent 0.99/0.95 (linear) for the above-threshold mode;
  below-threshold σ_max strictly η-independent along (1,1,1).
- `theory/bath_balance.py`: bracket + ladder printout.

## Files

- `v67/THEORY.md` — this document
- `v67/theory/mi_dispersion.mac` — Maxima symbolic checks (run: `maxima -b mi_dispersion.mac`)
- `v67/theory/mi_band.py` — band evaluation + 48×48 verification (`mi_band.out`)
- `v67/theory/mi_eta_scan.py` — η-scaling of the parametric mode (`mi_eta_scan.out`)
- `v67/theory/bath_balance.py` — bath-balance estimate (`bath_balance.out`)
