# EMF and Interaction Equations — Complete Derivation

## References

- Finster, Smoller, Yau (1999): "Particle-Like Solutions of the Einstein-Dirac Equations"
  arXiv:gr-qc/9801079, Phys. Rev. D 59:104020
- Finster, Smoller, Yau (1999): "Particle-Like Solutions of the Einstein-Dirac-Maxwell Equations"
  arXiv:gr-qc/9802012
- Liebling & Palenzuela: "Dynamical boson stars" (Living Reviews in Relativity)
- Dzhunushaliev et al.: "Dirac-boson stars" JHEP 02 (2023) 249

---

## System I: Einstein-Dirac-Soler (ED-S) — Gravitating Nonlinear Dirac

### Lagrangian

L = (1/16πG) R + ψ̄(iγ^μ ∇_μ - m)ψ + (λ/2)(ψ̄ψ)²

### Metric (Finster convention)

ds² = -T(r)⁻² dt² + A(r)⁻¹ dr² + r² dΩ²

where T(r) > 0 and A(r) > 0. Asymptotically: A → 1-2M/r, T → √(1-2M/r).

### Dirac Ansatz

Two fermions in singlet state (j = 1/2, N = 2):

ψ_{j,k} = e^{-iωt} (√T/r) [α(r) χ^k_{j-1/2}, i β(r) χ^k_{j+1/2}]^T

where α(r), β(r) are real radial functions.

**Scalar density**: S = ψ̄ψ = T(α² - β²)/r² (includes curved-space factor)

### Coupled Radial ODE System

**Dirac equations** (from arXiv:2503.15995, eq. 2.6–2.7, with N=2 → κ/2 = 1):

  √A · α' = +(1/r) α - (ωT + m) β + λS·β       ... (D1)
  √A · β' = +(ωT - m) α - (1/r) β - λS·α        ... (D2)

where S = T(α² - β²)/r² is the curved-space scalar density.

The self-interaction enters as m_eff = m - λS, so:
  √A · α' = +(1/r) α - (ωT + m_eff) β            ... (D1')
  √A · β' = +(ωT - m_eff) α - (1/r) β             ... (D2')

**Einstein equations** (eq. 2.8–2.9, with κ=2):

  rA' = 1 - A - 16πG ωT²(α² + β²) - 4πGr²λS²   ... (E1)
  2rA(T'/T) = A - 1 - 16πGT[ωT(α² + β²)
              - (2/r)αβ - m(α² - β²)] + 4πGr²λS² ... (E2)

### Boundary Conditions

**At r = 0:**
  α(r) = α₁ r + O(r³)      (shooting parameter: α₁)
  β(r) = β₂ r² + O(r⁴)     where β₂ = (ωT₀ - m + λT₀α₁²)α₁/3
  A(0) = 1                   (regularity)
  T(0) = T₀                  (gauge freedom — set to 1, rescale ω afterward)

**At r → ∞:**
  α, β → 0                   (exponentially: ~ e^{-κr}, κ = √(m²-ω²/T_∞²))
  A → 1 - 2GM/r
  T → √(1 - 2GM/r) → 1 (for Newtonian limit)

### Normalization

‖ψ‖² = 4π ∫₀^∞ (α² + β²) T/√A dr = 1

### Shooting Strategy

- Fix m = 1, λ = 1 (same as flat-space Soler)
- Set T₀ = 1 (gauge)
- Free parameter: α₁ (amplitude at origin)
- Eigenvalue: ω (frequency)
- Method: Bisect on α₁ for each ω, or use SCF iteration

### Newtonian Limit

For weak gravity (GM << r), set A = 1, T = 1 + Φ(r)/c² ≈ 1:
Reduces exactly to the flat-space Soler system with gravitational perturbation.

The leading correction is O(Gm²/ℏc) ≈ 10⁻³⁸ — negligible for nuclear-scale solitons.

**Implication**: Gravity is a perturbative correction to the Soler soliton, NOT the
primary confinement mechanism (unlike the boson star where gravity IS the confinement).

---

## System II: Einstein-Dirac-Maxwell (EDM) — Finster-Smoller-Yau (1999)

### References

- Finster, Smoller, Yau, Phys. Lett. A 259, 431-436 (1999), arXiv:gr-qc/9802012

### Lagrangian

L = (1/16πG) R + ψ̄(iγ^μ D_μ - m)ψ - (1/4)F_{μν}F^{μν}

where D_μ = ∇_μ + ieA_μ (Finster convention: +ie).

### Fields

Same metric ds² = -T⁻²dt² + A⁻¹dr² + r²dΩ² and Dirac ansatz as System I.
Electromagnetic: A_μ = (-φ(r), 0, 0, 0) — electrostatic potential.

### The Complete 5-ODE Radial System

**Dirac equations** (gauge coupling shifts ω → ω - eφ):

  √A · α' = +(1/r) α - ((ω - eφ)T + m) β              ... (D1-EM)
  √A · β' = +((ω - eφ)T - m) α - (1/r) β               ... (D2-EM)

**Einstein equations** (including Maxwell stress-energy):

  rA' = 1 - A - 16π(ω - eφ)T²(α² + β²)
        - r²AT²(φ')²                                     ... (E1-EM)

  2rA(T'/T) = A - 1 - 16π(ω - eφ)T²(α² + β²)
              + 32π(1/r)Tαβ + 16πmT(α² - β²)
              + r²AT²(φ')²                                ... (E2-EM)

**Maxwell equation:**

  r²Aφ'' = -8πe(α² + β²)
           - (2rA + r²A(T'/T) + (r²/2)A') · φ'          ... (M1)

This is **5 coupled ODEs** (4 first-order + 1 second-order = 6 first-order total).

### Boundary Conditions

**At r = 0** (gauge: T(0)=1, φ(0)=0, m=1):
  α(r) = α₁r + O(r²), β(r) = O(r²), A(0)=1, T(0)=1, φ(0)=0, φ'(0)=0

**At r → ∞:**
  α, β → 0 (exponential decay)
  A ~ T⁻² ~ 1 - 2ρ/r + (2e)²/r²  (Reissner-Nordström)
  φ ~ 2e/r

### Shooting Parameters

3 parameters: e (charge), ω (energy), α₁ (amplitude). Shoot on ω for each (e, α₁).
Solutions exist for ALL values of (e/m)², including (e/m)² > 1 where Coulomb > gravity.

### Key Results (FSY)

- (e/m)² < 1: infinite number of ground-state solutions on a spiral in (m, m-ω) plane.
  Lower branch is linearly stable.
- (e/m)² > 1: solutions only for large enough m. Binding energy > 0 (unstable).
- Binding: individual fermions always bound (m - ω > 0).
  Total binding: ρ - 2m < 0 for stable solutions when (e/m)² < 1.

### EDMS (with Soler self-interaction)

Add (λ/2)(ψ̄ψ)² to the EDM Lagrangian. The modifications are:
- Dirac eqs: m → m_eff = m - λS, where S = T(α²-β²)/r²
- Einstein eqs: add -4πr²λS² to (E1-EM), +4πr²λS² to (E2-EM)
- Maxwell eq: unchanged

### Flat-Space Charged Soler (A=1, T=1) — IMPLEMENTED

  f' = (m + (ω + eV) - λ(f²-g²)) · g                   ... (D1-flat)
  g' + 2g/r = (m - (ω + eV) - λ(f²-g²)) · f            ... (D2-flat)
  V'' + 2V'/r = -4πe · (f² + g²)                        ... (M1-flat)

Code: `src/charged_soler.c`. SCF iteration with under-relaxation α=0.2.
See RESULTS.md for numerical findings.

---

## System III: Charged Boson Star (Einstein-KG-Maxwell)

### Lagrangian

L = R/(16πG) + |D_μφ|² - μ²|φ|² - (1/4)F_μν F^{μν}

### Ansatz

Metric: ds² = -α²dt² + a²dr² + r²dΩ²
Scalar: φ = φ₀(r) e^{iωt}
EM: A_μ = (A₀(r), 0, 0, 0)

### Radial System (from Liebling & Palenzuela review)

  da/dr = a²[(1-a²)/r + 4πG(Φ² + a²((ω+eA₀)²/α² + μ²)φ₀²)]
  dα/dr = (α/2a²)[(a²-1)/r + 4πG(Φ² + a²((ω+eA₀)²/α² - μ²)φ₀²)]

  dΦ/dr = ...  (from KG equation)
  dA₀/dr = ... (from Maxwell equation)

where Φ = dφ₀/dr, and the gauge field satisfies Gauss's law.

### Maximum charge bound

The critical dimensionless charge: ẽ² = e²M_Pl²/(8πm²) < 1/2 for stable solutions.
Due to binding energy, solutions exist slightly beyond ẽ² = 1/2.

### Newtonian limit (extends boson_star.c)

  u'' + 2u'/r + 2(E - Φ_grav)u = 0    (Schrödinger, as before)
  Φ'' + 2Φ'/r = 4πG u²                 (Poisson for gravity)
  V'' + 2V'/r = -4πe · 2ω u²           (Poisson for Coulomb)

Energy eigenvalue shifts: E → E + eV(r), or equivalently ω_eff = ω + eV(r).

---

## Implementation Priority

### Step 1: Charged Soler (flat space) — FIRST TARGET
- Extend soler.c with Poisson solver for V(r)
- SCF iteration: Soler + Coulomb
- Scan over charge coupling e
- Key question: does charge unbind the soliton? At what critical e?

### Step 2: Charged Boson Star (Newtonian)
- Extend boson_star.c with Coulomb potential
- Three-way SCF: Schrödinger + Poisson(grav) + Poisson(EM)
- Maximum charge-to-mass ratio

### Step 3: Einstein-Dirac-Soler (full GR)
- New code: coupled (α, β, A, T) system
- Shooting on (α₁, ω) with A=1, T=1 at origin
- Compare with flat-space Soler: quantify gravitational correction
- For nuclear-scale solitons: correction is O(10⁻³⁸), negligible

### Step 4: Full EDMS
- Combine Steps 1 and 3: charged + gravitating Dirac soliton
- Five functions: (α, β, A, T, V)
- This is the Finster-Smoller-Yau system + Soler self-interaction

---

## System IV: Einstein-Born-Infeld-Dilaton (EBId) — Clement & Gal'tsov (2000)

### References

- Clement & Gal'tsov, PRD 62, 124013 (2000), arXiv:hep-th/0012059
- Garfinkle, Horowitz & Strominger, PRD 43 (1991) 3140
- Gibbons & Maeda, Nucl. Phys. B 298 (1988) 741

### Motivation

This is the ONLY model surveyed that produces BOTH solitonic solutions AND a massless
mediator giving 1/r long-range force. The dilaton φ is massless and couples universally
to energy density, acting as a scalar graviton.

### Action (Einstein frame)

S = ∫ d⁴x √(-g) [ R/(16πG) - (1/2)(∇φ)² + e^{-2aφ} b²(1 - √(1 + F²/(2b²))) ]

where:
- φ = dilaton (real scalar, massless)
- a = dilaton coupling (a=1 for heterotic string theory)
- b = Born-Infeld maximal field strength
- F² = F_μν F^μν

### String frame action (a=1)

S_string = ∫ d⁴x √(-g) e^{-2φ} [ R + 4(∇φ)² - F²/4 ]

### Radial System (schematic)

Metric: ds² = -f(r)dt² + f(r)⁻¹dr² + R(r)²dΩ²

- Einstein equations for f(r) and R(r)
- Dilaton: φ'' + (f'/f + 2R'/R)φ' = a·e^{-2aφ}E²/f
- BI Gauss law: E = Q·e^{2aφ}/√(R⁴ + Q²e^{4aφ}/b²)  (algebraic)

### Soliton BCs

At r=0: f(0)=1, R(0)=0, R'(0)=1, φ(0)=φ₀ (free parameter)
At r→∞: f→1, R→r, φ→φ_∞ + D/r (asymptotic flatness)

### Key Results (Clement-Gal'tsov)

- **Globally regular solitons exist** for any effective coupling ratio
- One-parameter family labeled by φ₀ (central dilaton value)
- Mass M and dilaton charge D both diverge at limiting parameter value
- For a=1: SL(2,R) symmetry

### Numerical Results — CONFIRMED (src/ebid.c)

Uses the Tamaki-Torii field equations (PRD 62, 061501) with metric
ds² = -f·e^{-2δ}dt² + dr²/f + r²dΩ² and κ²=2.

Solver: log-grid (s=ln r) RK4, auto-bisection on CG family parameter c for φ_inf=0.

**Q-scan** (b=1, γ=1):

| Q    | c_best | M      | D       | M/|D| | (M²+D²)/Q² | D/M    |
|------|--------|--------|---------|-------|-------------|--------|
| 0.01 | 18.05  | 0.0012 | -0.0023 | 0.514 | 0.069       | -1.950 |
| 0.1  | 3.80   | 0.032  | -0.055  | 0.586 | 0.40        | -1.706 |
| 0.2  | 2.70   | 0.080  | -0.127  | 0.631 | 0.56        | -1.583 |
| 0.5  | 2.19   | 0.251  | -0.356  | 0.706 | 0.76        | -1.416 |
| 1.0  | 2.17   | 0.568  | -0.738  | 0.770 | 0.87        | -1.299 |
| 2.0  | 2.32   | 1.235  | -1.486  | 0.831 | 0.93        | -1.204 |
| 5.0  | 2.65   | 3.284  | -3.669  | 0.895 | 0.97        | -1.117 |
| 10.0 | 2.95   | 6.685  | -7.237  | 0.924 | 0.97        | -1.083 |
| 50.0 | 3.70   | 34.17  | -35.51  | 0.962 | 0.97        | -1.039 |

**Key results**:
- Dilaton 1/r falloff: r·φ → D = -0.739 (Q=1) at large r (verified to r=2000)
- BPS limit: M²+D² → Q², D/M → -1 as Q → ∞
- D < 0 → attractive scalar force
- Dilaton force / gravitational force = (D/M)² ≈ 1.1–3.8
- BI field bounded: E(0) = √b (finite self-energy, no singularity)

### Dilaton 1/r Force

The asymptotic dilaton field φ(r) ~ φ_∞ + D/r mediates a TRUE 1/r force
between any two solitons:

  V_dilaton(R) = -D₁ D₂ / (4π R)

This is MASSLESS (no Yukawa screening!) and UNIVERSAL (couples to all energy).
This is structurally identical to Newtonian gravity.

### The Hierarchy Problem

The dilaton coupling strength D/M is determined by the parameter a and the
particular solution, NOT by the model. For the string value a=1, the coupling
is O(1) — far too strong compared to G_N ~ 10⁻³⁸. The hierarchy G_N/g_dilaton
is a free parameter, just as in v2/v3.

### Connection to v4 Program

- The BI sector provides self-trapping (Phase 3 mechanism)
- The dilaton provides massless 1/r force (what Phase 5 was missing)
- Together: a self-trapped, charged, gravitating soliton
- BUT: dilaton coupling strength is unconstrained (hierarchy problem)
- This is the most complete model, but it does NOT explain gravity's weakness

---

## Key Physics Questions

### Q1: What determines the charge?
- In the Soler model: Q_particle = e × (particle number Q)
- Q = 4π∫(f²+g²)r²dr is continuous (not quantized classically)
- Quantization requires quantum mechanics or topology
- The charge coupling e is a FREE PARAMETER — not determined by the model

### Q2: What's the maximum charge?
- Coulomb self-energy E_C ~ e²Q²/(R) scales as Q²
- Mass energy E_M ~ mQ scales as Q
- For E_C > E_M: the soliton is unstable to fission
- Critical charge: Q_max ~ mR/e² ~ (m/e)²R

### Q3: Can these form atoms?
- Two oppositely-charged solitons: V(R) = -e²/(4πR) - GM₁M₂/R
- For Soler solitons with m=1, R~2, λ=1: mass M ~ ωQ ~ 50
- Bohr radius: a_B = 1/(α_EM · M_reduced) where α_EM = e²/(4π)
- If e is the physical electron charge: a_B >> soliton size → well-separated atoms

### Q4: What mediates the force?
- **EM force**: mediated by photon (massless → 1/r) ✓
- **Scalar exchange**: mediated by (ψ̄ψ) fluctuation → Yukawa (massive → e^{-mr}/r) ✓
- **Gravity**: NO mechanism in this model — must be added by hand via Einstein equations
- **Lesson**: The Soler model + EM gives electromagnetism and short-range nuclear forces,
  but NOT gravity. Gravity requires additional structure (dilaton, or accepting GR as given).

### Q5: Connection to v2/v3 program?
- The Soler soliton is a DIFFERENT object from the Skyrmion:
  - Skyrmion: static topological soliton, SU(2) → S³ winding
  - Soler: oscillating dynamical soliton, U(1) charge conservation
- But both have spin-½ structure (Dirac for Soler, quantized isospin for Skyrmion)
- The nuclear force (Yukawa) emerges from scalar density exchange in BOTH models
- Gravity remains unsolved in BOTH
