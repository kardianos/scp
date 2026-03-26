# Nonlinear Maxwell Corrections from the Cosserat V(P) Potential

**OQ6 Resolution** — Perturbative expansion to second order

---

## 0. Setup and Notation

The 6-field Cosserat equations:

```
d²φ_a/dt² = ∇²φ_a - m²φ_a - (∂V/∂φ_a) + η (∇×θ)_a       (a = 0,1,2)
d²θ_a/dt² = ∇²θ_a            + η (∇×φ)_a                   (a = 0,1,2)
```

The potential:

```
V(P) = (μ/2) P² / (1 + κP²),     P = φ₀ φ₁ φ₂
```

The force on φ_a from V(P) is:

```
∂V/∂φ_a = (dV/dP)(∂P/∂φ_a) = (dV/dP) · P/φ_a
```

where P/φ_a is the product of the other two fields. The potential derivative:

```
dV/dP = μP / (1 + κP²)²
```

Parameters: m² = 2.25, μ = -41.345, κ = 50, η = 0.5, A_bg = 0.1.

**Background**: The uniform standing-wave background is

```
φ_a^(bg)(z) = A cos(kz + δ_a),     θ_a^(bg) = 0
```

with A = A_bg = 0.1, δ = {0, 3.0005, 4.4325}, and k determined by the
box size (k = 2π/L for periodic BC).

---

## 1. Derivatives of V(P)

We need the first three derivatives of V with respect to P:

```
V(P) = (μ/2) P² / (1 + κP²)

V'(P) = dV/dP  = μP / (1 + κP²)²

V''(P) = d²V/dP² = μ(1 - 3κP²) / (1 + κP²)³

V'''(P) = d³V/dP³ = μ(-6κP)(1 - κP²) / ... (not needed at second order)
```

**Check of V''**: Differentiate V'(P) = μP(1+κP²)^{-2}:

```
V''(P) = μ(1+κP²)^{-2} + μP · (-2)(1+κP²)^{-3} · 2κP
       = μ(1+κP²)^{-2} - 4μκP²(1+κP²)^{-3}
       = μ(1+κP²)^{-3} [(1+κP²) - 4κP²]
       = μ(1 - 3κP²) / (1+κP²)³
```

Confirmed.

At the background, P_bg = A³ cos(kz+δ₀) cos(kz+δ₁) cos(kz+δ₂).

**Magnitude estimate**: |P_bg| ≤ A³ = 10⁻³. Then κP_bg² ≤ 50 × 10⁻⁶ = 5×10⁻⁵ ≪ 1.

This means V(P) is deep in its quadratic regime: V(P) ≈ (μ/2)P² + O(κP⁴).

Consequently:

```
V'(P_bg)  ≈ μ P_bg
V''(P_bg) ≈ μ (1 - 3κP_bg²) ≈ μ
V'''(P_bg) ≈ -6μκ P_bg
```

The nonlinear corrections from V are suppressed by powers of κP_bg² ≈ 5×10⁻⁵.

---

## 2. Perturbative Expansion

Expand:

```
φ_a = φ_a^(bg) + ε δφ_a^(1) + ε² δφ_a^(2) + ...
θ_a =             ε δθ_a^(1) + ε² δθ_a^(2) + ...
```

where ε is a formal ordering parameter (set to 1 at the end).

The triple product expands as:

```
P = (φ₀^bg + ε δφ₀^(1) + ...)(φ₁^bg + ε δφ₁^(1) + ...)(φ₂^bg + ε δφ₂^(1) + ...)
```

Define shorthand: φ_a ≡ φ_a^bg, δ_a ≡ δφ_a^(1), Δ_a ≡ δφ_a^(2).

```
P = P₀ + ε P₁ + ε² P₂ + ...

P₀ = φ₀ φ₁ φ₂

P₁ = δ₀ φ₁ φ₂ + φ₀ δ₁ φ₂ + φ₀ φ₁ δ₂
   = Σ_a δ_a (P₀/φ_a)

P₂ = Δ₀ φ₁ φ₂ + φ₀ Δ₁ φ₂ + φ₀ φ₁ Δ₂    (linear in second-order)
   + δ₀ δ₁ φ₂ + δ₀ φ₁ δ₂ + φ₀ δ₁ δ₂       (quadratic in first-order)
   = Σ_a Δ_a (P₀/φ_a) + Σ_{a<b} δ_a δ_b (P₀/(φ_a φ_b))
```

The V' force expands as:

```
∂V/∂φ_a = V'(P) · (P/φ_a)
```

**Zeroth order**: F_a^(0) = V'(P₀) · (P₀/φ_a).

**First order**: Collecting terms linear in ε:

```
F_a^(1) = V''(P₀) P₁ (P₀/φ_a) + V'(P₀) [∂/∂ε(P/φ_a)]_{ε=0}
```

The second term requires care. We have P/φ_a = ∏_{b≠a} φ_b. At first order:

```
∂(P/φ_a)/∂ε|₀ = Σ_{b≠a} δ_b ∏_{c≠a,b} φ_c
```

For a=0: ∂(P/φ₀)/∂ε|₀ = δ₁ φ₂ + φ₁ δ₂.

So:

```
F_a^(1) = V''(P₀) P₁ (P₀/φ_a) + V'(P₀) Σ_{b≠a} δ_b ∏_{c≠a,b} φ_c
```

**Second order**: Collecting terms at ε²:

```
F_a^(2) = (1/2) V'''(P₀) P₁² (P₀/φ_a)           (from expanding V'(P₀+εP₁+...))
         + V''(P₀) P₂ (P₀/φ_a)                    (from V'' × P₂)
         + V''(P₀) P₁ [∂(P/φ_a)/∂ε]               (cross term)
         + V'(P₀) [∂²(P/φ_a)/∂ε²]                  (from two first-order perturbations)
```

The last term: ∂²(P/φ_a)/∂ε² at zeroth order gives the product of two δ's from
different fields:

For a=0: ∂²(φ₁φ₂)/∂ε² = 2 δ₁ δ₂.

So ∂²(P/φ_a)/∂ε² = 2 ∏_{b≠a} δ_b (there are only two b≠a, so this is just 2δ_b δ_c
where {b,c} = {0,1,2}\{a}).

---

## 3. Zeroth Order: Background Consistency

The background must satisfy its own equation:

```
∂²φ_a^bg/∂t² = ∇²φ_a^bg - m²φ_a^bg - V'(P₀)(P₀/φ_a) + η(∇×θ^bg)_a
```

For the standing-wave background φ_a = A cos(kz + δ_a), θ^bg = 0:

```
LHS: ∂²φ_a/∂t² = 0    (static background)
RHS: -k² A cos(kz+δ_a) - m² A cos(kz+δ_a) - V'(P₀)(P₀/φ_a)
   = -(k²+m²) φ_a - V'(P₀)(P₀/φ_a)
```

This does NOT vanish — the background is not an exact solution. However, the
potential term V'(P₀)(P₀/φ_a) is of order μA³ × A² ~ μA⁵ ≈ -41 × 10⁻⁵ = 4×10⁻³,
while the mass term is m²A = 2.25 × 0.1 = 0.225, and the Laplacian term k²A depends
on the box size.

**Resolution**: The background is approximate. For a self-consistent treatment,
the background should satisfy the static equation ∇²φ - m²φ - V'(P)(P/φ) = 0,
which determines k. In practice, the φ background is a standing wave with
k² = m² (or more precisely, k is set by the simulation domain). The V(P) correction
to the background is of order μA⁴/m² ≈ 0.02, a ~2% correction that we absorb into
the "slightly modified background" without changing the perturbative structure.

For the perturbative expansion, what matters is that the PERTURBATIONS around
whatever background exists are correctly captured. We proceed with the understanding
that the zeroth-order equation is approximately satisfied.

---

## 4. First Order: Linearized System (Maxwell)

The first-order equations are:

```
∂²δφ_a^(1)/∂t² = ∇²δφ_a^(1) - m²δφ_a^(1) - F_a^(1) + η(∇×δθ^(1))_a
∂²δθ_a^(1)/∂t² = ∇²δθ_a^(1)                          + η(∇×δφ^(1))_a
```

where F_a^(1) = V''(P₀) P₁ (P₀/φ_a) + V'(P₀) Σ_{b≠a} δ_b ∏_{c≠a,b} φ_c.

**In the small-background limit** (A_bg = 0.1):

Since P₀ ~ A³ = 10⁻³ and V'(P₀) ≈ μP₀ ~ -0.04, and V''(P₀) ≈ μ ≈ -41,
the force F_a^(1) contains terms of order:

- V''(P₀) P₁ (P₀/φ_a): This is μ × (δ_a × A² product) × A² ~ μA⁴ δ ~ 4×10⁻³ δ
- V'(P₀) × δ_b × φ_c: This is μA³ × δ × A ~ μA⁴ δ ~ 4×10⁻³ δ

Both are **much smaller** than the mass term m²δ = 2.25 δ. The V(P)-derived forces
at first order are a 0.2% correction to the mass term.

**Therefore**: At first order, the linearized system is effectively:

```
∂²δφ_a/∂t² ≈ ∇²δφ_a - m²δφ_a + η(∇×δθ)_a        (massive KG + coupling)
∂²δθ_a/∂t² = ∇²δθ_a           + η(∇×δφ)_a          (massless wave + coupling)
```

The V(P) terms are negligible at first order because the background amplitude A_bg = 0.1
makes the triple product P₀ ~ 10⁻³ extremely small.

This gives the standard Maxwell identification from EM_THEORY.md §5:
- θ_a → vector potential A_a
- E = -∂θ/∂t, B = ∇×θ
- J = η∇×φ (effective current)

---

## 5. Second Order: Nonlinear Corrections

The second-order equations are:

```
∂²δφ_a^(2)/∂t² = ∇²δφ_a^(2) - m²δφ_a^(2) - F_a^(2) + η(∇×δθ^(2))_a    ...(*)
∂²δθ_a^(2)/∂t² = ∇²δθ_a^(2)                          + η(∇×δφ^(2))_a    ...(**)
```

where F_a^(2) is the second-order V(P) force from §2.

**Equation (\*\*)** is the key: the θ (electromagnetic) equation at second order
is sourced by η∇×δφ^(2), where δφ^(2) is determined by equation (\*).

The coupling chain is:

```
δθ^(1) → [via η∇×δθ in φ-eq] → δφ^(1) → [via V(P) nonlinearity] → F^(2)
→ [drives] → δφ^(2) → [via η∇×δφ in θ-eq] → δθ^(2)
```

The nonlinear Maxwell correction enters via F_a^(2) driving δφ^(2), which then
sources δθ^(2) through the curl coupling.

### 5.1 Structure of the Second-Order Source

From §2, the second-order V(P) force is:

```
F_a^(2) = (1/2) V'''(P₀) P₁² (P₀/φ_a)
         + V''(P₀) P₂ (P₀/φ_a)
         + V''(P₀) P₁ Σ_{b≠a} δ_b ∏_{c≠a,b} φ_c
         + V'(P₀) ∏_{b≠a} δ_b
```

The last term is the most important. It is:

```
V'(P₀) ∏_{b≠a} δ_b = μ P₀ × ∏_{b≠a} δφ_b^(1)     (using V'≈μP in quadratic regime)
```

This is a DIRECT coupling between first-order perturbations in different field
components. When δφ^(1) is driven by δθ^(1) (the EM wave), this generates a
nonlinear EM self-interaction.

### 5.2 Magnitude Estimates

Let us estimate the strength of each term in F_a^(2) for a perturbation of
amplitude δ propagating through the background of amplitude A = 0.1.

| Term | Order | Value (δ=0.01, A=0.1) |
|------|-------|-----------------------|
| V'''(P₀) P₁² (P₀/φ_a) | μκ A³ × (δA²)² × A² = μκ A⁹ δ² | ~10⁻⁸ δ² |
| V''(P₀) P₂^{quad} (P₀/φ_a) | μ × δ² A × A² = μ A³ δ² | ~4×10⁻² δ² |
| V''(P₀) P₁ × δ_b φ_c | μ × δA² × δA = μ A³ δ² | ~4×10⁻² δ² |
| V'(P₀) ∏_{b≠a} δ_b | μ A³ × δ² | ~4×10⁻² δ² |

The V''' term (from κ-nonlinearity in V) is suppressed by κA⁶ ~ 5×10⁻⁵ relative to
the V'' terms. **The dominant nonlinearity comes from the multilinear structure of
P = φ₀φ₁φ₂, not from the nonlinear denominator (1+κP²).**

This is a crucial insight: even if V(P) were exactly quadratic (κ=0, so V = μP²/2),
the nonlinear Maxwell corrections would still exist at the same strength. The
nonlinearity is in P = φ₀φ₁φ₂ being a PRODUCT of three fields, not in V(P) being
a nonlinear function.

### 5.3 The Dominant Second-Order Correction

The dominant source for δφ^(2) comes from the last term in F^(2):

```
F_a^(2,dom) = V'(P₀) ∏_{b≠a} δφ_b^(1) ≈ μ P₀ · δφ_b^(1) δφ_c^(1)     ({b,c}={0,1,2}\{a})
```

In the equation for δφ_a^(2):

```
□δφ_a^(2) + m² δφ_a^(2) = -μ P₀(x) · δφ_b^(1) δφ_c^(1) + η(∇×δθ^(2))_a + ...
```

Here P₀(x) = A³ cos(kz+δ₀) cos(kz+δ₁) cos(kz+δ₂) is the spatially modulated
background triple product.

**The background's periodic structure acts as a spatially-varying coupling constant.**

This is analogous to nonlinear optics in a crystalline medium, where the χ^(2) susceptibility
has the periodicity of the lattice.

### 5.4 Expressing in Terms of EM Fields

At first order, the photon-like mode is mostly δθ with a small δφ admixture.
For a transverse EM wave propagating in direction n̂ with polarization ê:

```
δθ^(1) = ε_θ ê cos(k·x - ωt)
δφ^(1) = ε_φ ê cos(k·x - ωt)     (small admixture, ε_φ ~ η ε_θ / m²)
```

The φ-admixture amplitude is ε_φ ~ η k ε_θ / m² (from the curl coupling).
For sub-gap frequencies (ω < m), the φ response is evanescent and ε_φ/ε_θ ~ η/(m²/k - k).

The second-order source involves δφ_b^(1) δφ_c^(1), which is a PRODUCT of two
φ-perturbations from (potentially different) first-order EM waves. Since ε_φ ~ η ε_θ/m²,
the nonlinear source scales as:

```
F^(2) ~ μ A³ × (η/m²)² ε_θ² = μ η² A³ ε_θ² / m⁴
```

Numerically: μ η² A³ / m⁴ = (-41.345)(0.25)(10⁻³)/(2.25²) ≈ -2.0 × 10⁻³.

This means the second-order V(P) force is about 2×10⁻³ times the square of the
first-order θ amplitude.

---

## 6. Physical Effects

### 6.1 Photon-Photon Scattering (Euler-Heisenberg Analog)

Two first-order EM waves (δθ^(1,A) and δθ^(1,B)) induce φ-perturbations
(δφ^(1,A) and δφ^(1,B)). The product δφ_b^(1,A) × δφ_c^(1,B) sources a
second-order δφ^(2), which in turn sources δθ^(2) via the curl coupling.

This is a 2→2 photon scattering process mediated by the φ field.

**Effective coupling**: The scattering amplitude scales as:

```
M ~ μ η⁴ A³ / m⁸     (four η factors: two for θ→φ conversion, two for φ→θ)
```

Wait — let us be more careful. The full chain is:

```
θ^(1,A) →[η∇×] δφ^(1,A)  (amplitude ~ η k/m² × θ_A)
θ^(1,B) →[η∇×] δφ^(1,B)  (amplitude ~ η k/m² × θ_B)
δφ^(1,A) × δφ^(1,B) →[V(P)] δφ^(2)  (amplitude ~ μA³ × δφ_A × δφ_B / m²)
δφ^(2) →[η∇×] δθ^(2)  (amplitude ~ η k/m² × ??? )
```

But δφ^(2) is NOT massive-suppressed in the same way: it is driven by a source, so
the Green's function is (□ + m²)⁻¹ ≈ 1/m² for low-frequency sources. Then:

```
δφ^(2) ~ μA³ (η k/m²)² θ_A θ_B / m²
δθ^(2) ~ η k × δφ^(2) / k² ~ η/(k m²) × μA³ (η/m²)² k² θ_A θ_B
        = μ η³ A³ k / m⁶ × θ_A θ_B
```

The effective four-photon vertex has coupling strength:

```
g_eff ~ |μ| η³ A³ k / m⁶
```

For k ~ 1 (wavelength ~ code unit ~ 0.56 fm): g_eff ~ 41 × 0.125 × 10⁻³ / (1.5⁶)
≈ 5.1 × 10⁻³ / 11.4 ≈ 4.5 × 10⁻⁴.

### 6.2 Comparison to QED Euler-Heisenberg

The QED Euler-Heisenberg effective Lagrangian for low-energy photon-photon scattering is:

```
L_EH = (2α²/45m_e⁴) [(F_μν F^μν)² + (7/4)(F_μν F̃^μν)²]
```

where α ≈ 1/137 and m_e is the electron mass. The dimensionless coupling at
energy ω is:

```
g_QED ~ α² (ω/m_e c²)⁴ / (45)
```

In our theory, the analogous coupling is:

```
g_Cos ~ |μ| η³ A³ / m⁶ × (k/m)
```

The critical difference: **in QED, the coupling is suppressed by (ω/m_e)⁴ — it
vanishes at low energy. In the Cosserat theory, the coupling is suppressed by
A³ (background amplitude) but NOT by (ω/m)⁴.**

At a given energy scale, the Cosserat nonlinearity is MUCH STRONGER than QED
because it does not require virtual pair creation — the nonlinearity is already
present in the classical background field.

**Ratio**: At photon energy ω ~ 1 keV (ω/m_e ≈ 0.002):

```
g_QED ~ (1/137)² × (0.002)⁴ / 45 ~ 10⁻¹⁶
g_Cos ~ 4.5 × 10⁻⁴
```

The Cosserat nonlinearity is ~10¹² times stronger than QED at this energy.

**However**, this comparison is misleading in an important way. The Cosserat
coupling g_Cos ∝ A³ depends on the background amplitude A_bg. In truly empty
space (A_bg = 0), the coupling vanishes EXACTLY — there is no photon-photon
scattering in vacuum without a background φ field. The nonlinearity requires
a pre-existing medium (the φ background) to mediate the interaction. This is
fundamentally different from QED, where virtual electron-positron pairs
provide the medium spontaneously.

**Physical interpretation**: The Cosserat nonlinear correction is NOT vacuum
photon-photon scattering. It is the nonlinear optical response of the φ
background — analogous to χ^(2) nonlinear optics in a crystal. The "vacuum"
of this theory (with a φ background) acts like a nonlinear optical medium.

### 6.3 Vacuum Birefringence

**Question**: Does the background φ_a = A cos(kz + δ_a) make the vacuum
birefringent (different propagation speeds for different polarizations)?

The background has a preferred direction: the z-axis (along which the standing
wave varies). Consider two orthogonal polarizations of a transverse wave
propagating along z: x-polarization and y-polarization.

For the linearized θ equation (first order), the curl coupling η∇×δφ is
isotropic — it does not distinguish polarizations. So at first order, there is
NO birefringence.

At second order, the V(P) correction has the structure:

```
F_a^(2) ~ μ P₀(z) × product of perturbations
```

The key is P₀(z) = A³ cos(kz) cos(kz+δ₁) cos(kz+δ₂). This is a function of z only.

For a wave propagating along x (perpendicular to the background modulation),
with polarization in y or z: the y-polarized wave sees a z-independent average
of P₀, while the z-polarized wave couples to P₀ differently because it creates
φ perturbations along the direction of P₀ variation.

**However**, the fields φ_0, φ_1, φ_2 are scalar fields with no intrinsic
directional coupling to the spatial x,y,z axes. The "x-polarization" of θ
means θ₀ is excited, which couples to φ₀ via the curl. The nonlinear coupling
is between the FIELD INDICES (0,1,2), not the spatial directions.

The V(P) force on δφ_a involves δφ_b × δφ_c for b,c ≠ a. For a single EM wave
with only one polarization excited (say δθ₀ ≠ 0, δθ₁ = δθ₂ = 0), we get
δφ₀ ≠ 0, δφ₁ = δφ₂ = 0. Then:

```
F₀^(2) ~ μ P₀ × δφ₁ × δφ₂ = 0     (both zero)
F₁^(2) ~ μ P₀ × δφ₀ × δφ₂ = 0     (δφ₂ = 0)
F₂^(2) ~ μ P₀ × δφ₀ × δφ₁ = 0     (δφ₁ = 0)
```

**A single-polarization wave produces NO second-order V(P) correction.**

This is because P = φ₀φ₁φ₂ requires ALL THREE field components to be nonzero.
A perturbation in only one component cannot modify P at first order (P₁ = 0 if
only one δφ_a ≠ 0 while the other two background fields are unperturbed... wait,
this is not quite right.)

Let me reconsider. P₁ = δ₀ φ₁^bg φ₂^bg + ..., so even a single δφ₀ gives P₁ ≠ 0
because φ₁^bg and φ₂^bg are the nonzero background fields. The V'' terms DO contribute.

Going back to the full F_a^(1):

```
F₀^(1) = V''(P₀) P₁ (P₀/φ₀) + V'(P₀) (δ₁ φ₂ + φ₁ δ₂)
```

With only δ₀ ≠ 0, δ₁ = δ₂ = 0:
```
P₁ = δ₀ φ₁ φ₂
F₀^(1) = V''(P₀) δ₀ φ₁ φ₂ (P₀/φ₀) + 0 = V''(P₀) δ₀ (φ₁ φ₂)² ≈ μ δ₀ A⁴
F₁^(1) = V''(P₀) δ₀ φ₁ φ₂ (P₀/φ₁) + V'(P₀) δ₀ φ₂ = μ P₀ δ₀ φ₂ (P₀/(φ₁ P₀) + 1)
```

Hmm, this is getting notationally complicated. Let me switch to a cleaner approach.

**Cleaner approach**: Work directly with the force ∂V/∂φ_a and expand.

The potential is V(P) with P = φ₀φ₁φ₂. The force on field a is:

```
∂V/∂φ_a = V'(P) ∂P/∂φ_a = V'(P) ∏_{b≠a} φ_b
```

Now expand φ_a = φ_a^bg + δφ_a. Define Π_a ≡ ∏_{b≠a} φ_b^bg.

```
∂V/∂φ_a = V'(P) ∏_{b≠a}(φ_b^bg + δφ_b)
```

Since we are in the quadratic regime (κP² ≪ 1), V'(P) ≈ μP. Then:

```
∂V/∂φ_a ≈ μ (φ_a^bg + δφ_a) ∏_{b≠a}(φ_b^bg + δφ_b) × ∏_{b≠a}(φ_b^bg + δφ_b)
```

Wait, this is circular. Let me just directly expand ∂V/∂φ_a = μP ∏_{b≠a} φ_b + O(κ).

Actually, ∂V/∂φ_a = V'(P) ∏_{b≠a} φ_b where P = φ₀φ₁φ₂ = φ_a ∏_{b≠a} φ_b, and
V'(P) = μP/(1+κP²)² ≈ μP for small P. So:

```
∂V/∂φ_a ≈ μ φ_a (∏_{b≠a} φ_b)²
```

Now expand each φ_b = φ_b^bg + δφ_b. For field 0 with {b,c} = {1,2}:

```
(∏_{b≠0} φ_b)² = (φ₁ φ₂)² = (φ₁^bg + δφ₁)² (φ₂^bg + δφ₂)²
```

Expanding to second order in δ:

```
(φ₁^bg + δφ₁)² = (φ₁^bg)² + 2φ₁^bg δφ₁ + (δφ₁)²
```

So:

```
(φ₁ φ₂)² = (φ₁^bg)²(φ₂^bg)² + 2φ₁^bg δφ₁ (φ₂^bg)² + 2(φ₁^bg)² φ₂^bg δφ₂
            + (δφ₁)²(φ₂^bg)² + 4 φ₁^bg δφ₁ φ₂^bg δφ₂ + (φ₁^bg)²(δφ₂)² + ...
```

The full force on δφ₀ at second order from V(P) involves terms like:

```
μ δφ₀ [(δφ₁)²(φ₂^bg)² + 4 φ₁^bg δφ₁ φ₂^bg δφ₂ + (φ₁^bg)²(δφ₂)²]
+ μ φ₀^bg [(δφ₁)²(φ₂^bg)² + 4 φ₁^bg δφ₁ φ₂^bg δφ₂ + (φ₁^bg)²(δφ₂)²]
```

**For a single-polarization wave (only δφ₀ ≠ 0, δφ₁ = δφ₂ = 0)**:

All the δφ₁, δφ₂ terms vanish. We only get contributions from the expansion
of φ₀^bg + δφ₀ in the prefactor:

```
∂V/∂φ₀ ≈ μ(φ₀^bg + δφ₀)(φ₁^bg φ₂^bg)²
```

This is LINEAR in δφ₀. There is NO quadratic (second-order) correction.

**For a single-polarization wave, the V(P) nonlinearity produces no
second-order correction.** The wave equation remains linear in δφ₀.

This is because V ∝ P² = (φ₀φ₁φ₂)², and if only φ₀ varies while φ₁, φ₂ are
fixed at their background values, V is quadratic in φ₀ and the force is linear.

**Result: No vacuum birefringence from V(P).**

Since a single-polarization wave gets no nonlinear correction, the dispersion
relation is the same for all polarizations. The Cosserat vacuum (with a φ
background) is NOT birefringent at second order. This is a symmetry result:
it follows from V depending on P = φ₀φ₁φ₂ as P², making the force linear in
each individual field when the others are held fixed.

### 6.4 Two-Polarization Interactions

The nonlinear correction requires AT LEAST TWO different field components to be
perturbed simultaneously. Consider two EM waves exciting δφ₀ and δφ₁:

```
δφ₀ = ε_A cos(k_A · x - ω_A t)    (from EM wave A)
δφ₁ = ε_B cos(k_B · x - ω_B t)    (from EM wave B)
```

Then the second-order force on φ₂ is:

```
F₂^(2) ≈ μ φ₂^bg × 2 φ₀^bg δφ₀ × 2 φ₁^bg δφ₁
        = 4μ φ₂^bg φ₀^bg φ₁^bg δφ₀ δφ₁
        = 4μ P₀ δφ₀ δφ₁
```

And the product δφ₀ × δφ₁ contains:

```
ε_A ε_B cos(k_A·x - ω_A t) cos(k_B·x - ω_B t)
= (ε_A ε_B / 2) [cos((k_A+k_B)·x - (ω_A+ω_B)t) + cos((k_A-k_B)·x - (ω_A-ω_B)t)]
```

This generates:

1. **Sum frequency**: ω_A + ω_B, k_A + k_B (second harmonic if A=B)
2. **Difference frequency**: ω_A - ω_B, k_A - k_B (optical rectification if A=B)

These are the standard χ^(2) nonlinear optical effects, but they require DIFFERENT
field components (different "polarizations" in the field-index sense, not the spatial
sense).

### 6.5 Harmonic Generation

For a single EM wave that excites ALL THREE field components (e.g., circularly
polarized or at an angle that mixes all three), the δφ_a δφ_b cross-terms
generate second harmonics.

The condition is that the EM wave must excite at least two of the three φ-components.
Since each θ_a couples to φ_a via the curl, a wave with wavevector along z and
polarization in the x-y plane excites θ₀ and θ₁ (for x and y polarizations),
which couple to φ₀ and φ₁. This gives the two-component interaction of §6.4.

A wave propagating along a general direction n̂ can excite all three θ-components
and hence all three φ-components. In this case, the full three-field nonlinear
coupling activates.

**Selection rules**: The coupling P = φ₀φ₁φ₂ is totally symmetric under
permutations of field indices. The second-order correction requires at least
two different field indices to be excited. There is no special enhancement or
suppression for any particular pair — the coupling is democratic.

### 6.6 Self-Focusing / Self-Defocusing

Self-focusing requires a nonlinear refractive index: n = n₀ + n₂|E|². From
the analysis above, a single-polarization wave gets NO nonlinear correction.
Therefore:

**A single-polarization EM wave does NOT self-focus or self-defocus.**

For a multi-polarization wave (e.g., unpolarized light), the cross-polarization
coupling can produce an effective n₂. The sign depends on μ:

Since μ = -41.345 < 0, the V(P) potential is ATTRACTIVE (energy decreases with P²).
The nonlinear correction makes the effective mass SMALLER where the field is
stronger, which increases the refractive index. This gives n₂ > 0:
**self-focusing** for multi-polarization light.

---

## 7. Effective Nonlinear Lagrangian

The second-order correction can be written as an effective 4-field interaction.
In the quadratic regime (V ≈ μP²/2), the V(P) energy density is:

```
V = (μ/2)(φ₀φ₁φ₂)²
```

This is already a sixth-order interaction in the fields. When we split
φ_a = φ_a^bg + δφ_a and integrate out the massive δφ, the effective
Lagrangian for the θ field picks up a quartic correction:

```
L_NL = g_NL (θ_a ∂_i θ_b - θ_b ∂_i θ_a)(θ_a ∂_i θ_b - θ_b ∂_i θ_a)
```

where the effective coupling is:

```
g_NL ~ μ η⁴ A_bg² / m⁸
```

Let us estimate this more carefully. The chain θ → φ → V → φ → θ involves:

1. Two θ→φ conversions: each contributes η k / (k² + m²) ≈ η k / m² for k ≪ m
2. One V(P) vertex: contributes μ × (background factors) ~ μ A_bg²
3. One φ→θ conversion: contributes η k / k² ≈ η/k for the outgoing θ
4. One massive propagator for δφ^(2): contributes 1/m²

Total: g_NL ~ μ × (ηk/m²)² × A_bg² × (η/k) × (1/m²) = μ η³ A_bg² k / m⁶

But we also need to account for the spatial averaging over the background
modulation. The background P₀(z) = A³ ∏ cos(kz+δ_a) oscillates with period
2π/k. For processes that average over many wavelengths, only the spatial mean
of P₀² matters:

```
⟨P₀²⟩ = A⁶ ⟨cos²(kz+δ₀) cos²(kz+δ₁) cos²(kz+δ₂)⟩
```

For generic phases δ_a, this spatial average is A⁶/8 (each cos² averages to 1/2).
With the specific phases δ = {0, 3.0005, 4.4325}, which are chosen to minimize
certain resonances, the average is still close to A⁶/8.

---

## 8. Summary of Nonlinear EM Effects

| Effect | Present? | Mechanism | Scaling |
|--------|----------|-----------|---------|
| Photon-photon scattering | YES | Cross-component V(P) coupling | ~μη⁴A²/m⁸ |
| Vacuum birefringence | **NO** | V ∝ P² is linear in each φ_a separately | — |
| Second harmonic generation | YES | δφ_a × δφ_b cross-terms | ~μη²A²/m⁴ |
| Difference frequency gen. | YES | Same mechanism, difference freq. | ~μη²A²/m⁴ |
| Self-focusing (single pol.) | **NO** | No self-nonlinearity for one component | — |
| Self-focusing (multi pol.) | YES | Cross-polarization coupling, μ<0 | n₂ > 0 |

**The key structural result**: The absence of birefringence and single-polarization
self-interaction is a direct consequence of P = φ₀φ₁φ₂ being MULTILINEAR
(degree 1 in each field). V(P²) is quadratic in each individual field component
when the others are fixed. The nonlinearity is an INTER-FIELD coupling, not a
SELF-coupling.

This is physically natural: a single "color" of the Cosserat field propagates
linearly through the background. Nonlinear effects only appear when multiple
"colors" overlap and their product P becomes dynamically relevant.

---

## 9. Comparison with QED

| Property | QED (Euler-Heisenberg) | Cosserat |
|----------|----------------------|----------|
| Mechanism | Virtual e⁺e⁻ loops | Classical φ-background |
| Scaling | α²(ω/m_e)⁴ | μη⁴A²/m⁸ (background-dependent) |
| Vacuum (no background) | Active | **Zero** |
| Birefringence | Yes (7:4 ratio) | **No** |
| Low-energy limit | ∝ ω⁴ (vanishes) | ∝ A² (constant if background exists) |
| UV behavior | Requires renormalization | Natural cutoff at m |

**Critical difference 1**: QED nonlinear effects exist in true vacuum
(they arise from virtual pairs). The Cosserat nonlinearity requires a
pre-existing φ background (A_bg ≠ 0). In regions of space far from any
matter, where A_bg → 0, the Cosserat nonlinearity vanishes.

**Critical difference 2**: QED produces vacuum birefringence (different indices
for E∥B and E⊥B orientations relative to an external field, with the famous 7:4
ratio). The Cosserat theory does NOT, because V(P²) is separately quadratic in each
field component. This is a TESTABLE DISTINCTION if the Cosserat φ background is
ever measurable near matter.

**Critical difference 3**: The Cosserat nonlinearity acts as a χ^(2) process
(three-wave mixing with the background counted as one wave), while QED is purely
χ^(3) (four-wave). This means the Cosserat theory predicts parametric frequency
conversion (SHG, DFG) near matter, while QED does not at tree level.

---

## 10. Numerical Estimates

Using m² = 2.25 (m = 1.5), μ = -41.345, η = 0.5, A_bg = 0.1:

**Characteristic nonlinear field strength** (where second-order corrections
become comparable to first-order):

The ratio of second-order to first-order force is:

```
F^(2)/F^(1) ~ μ A² δφ / m² ~ |μ| A² (η/m²) δθ / m²
             = |μ| η A² / m⁴ × δθ
             = 41 × 0.5 × 0.01 / 5.06 × δθ
             ≈ 0.04 × δθ
```

So perturbation theory breaks down when δθ ~ 25 (in code units), which corresponds
to extremely intense fields. For comparison, the θ field around a braid has
amplitude θ_rms ~ 0.03 at r = 1 (from V34 data). The nonlinear correction is:

```
F^(2)/F^(1) ~ 0.04 × 0.03 = 1.2 × 10⁻³
```

**The nonlinear correction is about 0.1% near a braid core.** This is small
but potentially measurable in high-precision simulations.

**In physical units** (using 1 code E = 9.098 MeV, 1 code L = 0.5624 fm):

The critical field for nonlinear effects: E_crit ~ m²/(|μ|ηA²) in code units,
which corresponds to an EM field energy density of ~10⁴ MeV/fm³. For comparison,
the QED critical field (Schwinger limit) corresponds to ~10⁹ MeV/fm³. The
Cosserat nonlinearity becomes important at field strengths ~10⁵ times LOWER than
the Schwinger limit.

---

## 11. Effect of Background Periodicity

The standing-wave background φ_a^bg = A cos(kz + δ_a) introduces additional
structure beyond the simple estimates above.

### 11.1 Quasi-Phase-Matching

The background triple product P₀(z) has spatial Fourier components at
wavenumbers 0, ±2k, and ±(other combinations of δ-dependent phases). The
nonlinear coupling is spatially modulated, which means:

- **Phase-matched processes** (where k_out = k_1 + k_2 + n·k_bg for integer n)
  are resonantly enhanced over a coherence length.
- **Non-phase-matched processes** average out over many background wavelengths.

This is exactly analogous to quasi-phase-matching in periodically-poled
nonlinear crystals. The background modulation can enhance certain frequency
conversion processes that would otherwise be forbidden by momentum conservation.

### 11.2 Band Structure

The periodic background creates a photonic band structure for the φ field.
The first-order φ equation:

```
□δφ_a + m²δφ_a ≈ -μ A⁴ cos²(kz+δ_b) cos²(kz+δ_c) δφ_a + η(∇×δθ)_a
```

contains a periodic potential μA⁴ cos²(...)cos²(...) with period π/k.
This creates Bragg scattering for φ waves with wavenumber near k (or
multiples), opening small band gaps of width ~μA⁴/m² ≈ 2×10⁻³ m.

These gaps are extremely small (0.1% of the mass gap m) and would only
affect φ modes very close to the Bragg condition. For the photon-like (θ)
modes, which are predominantly in the massless sector, the effect is
further suppressed by the mixing angle η²/m².

### 11.3 Spatial Averaging for Long-Wavelength Processes

For EM waves with wavelength λ ≫ 2π/k (the background period), the
relevant coupling is the spatial average of the nonlinear vertex. This
averages P₀(z) over the background period, giving:

```
⟨P₀⟩ = A³ ⟨cos(kz+δ₀) cos(kz+δ₁) cos(kz+δ₂)⟩
```

Expanding the product of three cosines:

```
cos α cos β cos γ = (1/4)[cos(α+β-γ) + cos(α-β+γ) + cos(-α+β+γ) + cos(α+β+γ)]
```

The spatial average is zero unless one of the arguments is zero (independent of z).
With α = kz+δ₀, β = kz+δ₁, γ = kz+δ₂:

- α+β-γ = kz + (2δ₀+2δ₁-2δ₂)/2... no, more carefully:
- α+β-γ = kz + δ₀ + kz + δ₁ - kz - δ₂ = kz + δ₀+δ₁-δ₂

This still depends on z (linear in z with coefficient k). Therefore ⟨P₀⟩ = 0.

**The spatial average of the background triple product vanishes.**

This means the leading nonlinear coupling for long-wavelength photons averages
to zero. The first surviving contribution comes from ⟨P₀²⟩ ≠ 0, which enters
at FOURTH order (not second). This further suppresses the nonlinear effects
for wavelengths much longer than the background modulation period.

For wavelengths COMPARABLE to the background period, the full spatially-resolved
coupling must be used, and quasi-phase-matching effects become important.

---

## 12. Conclusions

1. **The dominant nonlinearity is MULTILINEAR, not from V's functional form.**
   Even with V = μP²/2 (the κ→0 limit), the nonlinearity P = φ₀φ₁φ₂ creates
   inter-field coupling at second order. The κ-dependent corrections are
   suppressed by κA⁶ ~ 5×10⁻⁵.

2. **No vacuum birefringence.** V(P²) is quadratic in each field component
   separately, so a single-polarization wave propagates linearly. This is a
   structural result independent of parameters.

3. **No single-polarization self-interaction.** Self-focusing, self-phase
   modulation, and other single-beam effects are absent for a single field
   component. Multi-component beams DO interact.

4. **Cross-polarization coupling exists.** Two waves exciting different φ_a
   components interact via the V(P) potential, producing sum/difference
   frequency generation. This is a χ^(2)-type process (three-wave mixing
   with the background as the third wave).

5. **The background acts as a nonlinear medium.** All nonlinear effects
   vanish when A_bg → 0. The Cosserat "vacuum" is nonlinear only where
   the φ background is present (near matter).

6. **Spatial averaging kills the leading coupling.** For long-wavelength
   photons, ⟨P₀⟩ = 0 and the first surviving nonlinear vertex comes from
   ⟨P₀²⟩, entering at fourth order. Nonlinear EM effects are most important
   at wavelengths comparable to the background modulation.

7. **The nonlinear effects are ~0.1% near a braid core** and negligible in
   the far field. The critical field strength is ~10⁵ times below the QED
   Schwinger limit, but requires a pre-existing φ background.

8. **Testable prediction**: Unlike QED, which predicts vacuum birefringence,
   the Cosserat theory predicts NO birefringence but YES parametric frequency
   conversion near matter. These are qualitatively different experimental
   signatures.
