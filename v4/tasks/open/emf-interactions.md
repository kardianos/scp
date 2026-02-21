# Task: EMF and General Interaction Equations for Self-Confined Solitons

## Motivation

Two models from the v4 program successfully produce self-confined localized solutions:
- **Model D (Boson Star)**: Complex scalar + gravity. Gravity IS the self-interaction.
- **Model E (Soler)**: Nonlinear Dirac with (ψ̄ψ)² coupling. Spin-½ from Dirac structure.

Neither model includes electromagnetism. For a complete particle model, we need:
1. **U(1) gauge coupling** — how do these solitons interact with EM fields?
2. **Multi-soliton interactions** — do they scatter, bind, form atoms?
3. **Gravitational coupling for Soler** — Model E has no gravity; Model D has no spin.

The goal: write down the coupled field equations for a **spinning, charged, gravitating soliton** and identify which pieces are tractable.

---

## Phase 1: Einstein-Dirac-Maxwell (EDM) System

The most natural unification: a Dirac field with self-interaction, minimally coupled to both gravity and electromagnetism.

**Lagrangian**:
L = R/(16πG) + ψ̄(iγ^μ D_μ - m)ψ + λ(ψ̄ψ)² - (1/4)F_μν F^{μν}

where D_μ = ∂_μ - ieA_μ (U(1) gauge covariant derivative).

### Step 1.1: Write the coupled field equations
- Einstein: G_μν = 8πG T_μν (Dirac + EM stress-energy)
- Dirac: (iγ^μ D_μ - m + 2λ(ψ̄ψ))ψ = 0
- Maxwell: ∇_μ F^{μν} = j^ν = eψ̄γ^ν ψ

### Step 1.2: Spherically symmetric reduction
- Metric: ds² = -A(r)dt² + B(r)dr² + r²dΩ²
- Dirac: same (f,g) ansatz but with gauge-covariant frequency ω → ω - eΦ(r)
- EM: only A₀ = Φ(r) for static electric field, A_i = 0 (no magnetic monopole)

### Step 1.3: Literature survey
- Finster, Smoller, Yau (1999): Einstein-Dirac-Maxwell particlelike solutions
- Herdeiro & Radu: charged boson stars
- Jetzer (1992): review of boson stars

### Step 1.4: Radial ODE system
Write the coupled (f, g, Φ_grav, M, Φ_em) system and identify if it's tractable
with the same shooting/SCF approach used for boson_star.c and soler.c.

---

## Phase 2: Charged Boson Star (Einstein-Klein-Gordon-Maxwell)

Simpler test case: complex scalar + gravity + EM, no spin.

**Lagrangian**:
L = R/(16πG) + |D_μ φ|² - μ²|φ|² - (1/4)F_μν F^{μν}

### Step 2.1: Radial equations
- Metric: same static spherically symmetric
- Scalar: φ = σ(r)e^{-iωt}, D_t φ = (-iω - ieΦ)φ → ω_eff = ω + eΦ(r)
- EM: Φ'' + (2/r - (A'/A+B'/B)/2)Φ' = -e·j⁰·B (sourced by Noether charge)

### Step 2.2: Implementation
Extend boson_star.c to include Coulomb self-energy. The EM field adds a repulsive
1/r contribution competing with gravitational attraction.

### Step 2.3: Key question
Does the charged boson star have a maximum charge-to-mass ratio? Is there a BPS-like
bound Q/M ≤ 1 (in Planck units)?

---

## Phase 3: Soliton-Soliton Interactions

### Step 3.1: Long-range forces
For two well-separated solitons, what are the interaction channels?
- Gravity: always attractive, 1/r
- EM: repulsive for same-sign charge, 1/r
- Scalar/Dirac exchange: Yukawa e^{-μr}/r (from mass gap)

### Step 3.2: Binding condition
Stable atom = gravity + Yukawa attraction > Coulomb repulsion at some radius.
This is the classical hydrogen analog. Derive the effective potential.

### Step 3.3: What determines charge?
In boson star: charge = e × particle number Q. Quantized? Or continuous?
In Soler: charge = e × fermion number. Quantized by Dirac structure.

---

## Phase 4: DBI Generalization

Replace standard kinetic term with Born-Infeld/DBI:
- Scalar DBI: P(X) = b²(√(1+2X/b²) - 1) instead of X
- Dirac DBI: more subtle — need to define DBI extension of γ^μ∂_μ

### Step 4.1: Does DBI + gravity give something new?
The DBI refractive index creates an ADDITIONAL self-trapping mechanism beyond gravity.
Two length scales: gravitational (R_grav ~ 1/√(Gμ²)) and DBI (R_DBI ~ 1/b).

### Step 4.2: Can DBI generate effective gravity?
From Phase 5 results: NO. The BLV effective metric is Yukawa, not 1/r.
But: DBI + dilaton might work (string theory low-energy action).

---

## Success Criteria

- [x] Write down complete EDM radial ODE system (Phase 1.4) — see emf-derivation.md
- [x] Identify number of shooting parameters and boundary conditions — 3 params (e, ω, α₁)
- [x] Determine if system is tractable with existing numerical methods — YES (SCF + shooting)
- [x] Literature comparison with Finster-Smoller-Yau results — solutions for ALL (e/m)²
- [x] Implement flat-space charged Soler — DONE (src/charged_soler.c)
- [x] Implement full GR Einstein-Dirac-Soler solver (4 ODEs) — DONE (src/ed_soler.c)
- [x] Implement full EDMS solver — DONE (src/edm_soler.c, SCF iteration)
- [x] Implement Einstein-BI-Dilaton soliton solver — DONE (src/ebid.c)
- [ ] If tractable: implement and find first charged spinning soliton

## Priority

1. **Flat-space charged Soler** — DONE. Critical charge e_crit≈0.23 at ω=0.9.
2. **Einstein-Dirac-Soler** — DONE. Gravitational corrections: contracts soliton,
   E_grav/M = 2-5% at GM/R~0.03-0.06, O(10⁻³⁸) for physical nucleons.
3. **Full EDMS** — DONE. SCF iteration (NOT simultaneous integration).
   Charge scan: e_crit≈0.12-0.13 (N=2) at ω=0.9. Combined G+e verified.
4. **Einstein-BI-Dilaton** — DONE. Globally regular solitons found for Q=0.01–50.
   Massless dilaton 1/r force CONFIRMED numerically (φ~D/r, D<0 → attractive).
   Q=1,b=1,γ=1: M=0.569, D=-0.739, (M²+D²)/Q²=0.87. BPS limit approached at large Q.
   Code: `src/ebid.c`. See RESULTS.md for full numerical results.
