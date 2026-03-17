# SCP Discoveries — Complete Record

## V28: Bimodal Synergy — Complete Equations and Replication Guide (March 2026)

### 1. Lagrangian

Three real scalar fields φ_a(x,t), a ∈ {0,1,2}, with NO mass term:

    L = Σ_a [½(∂_t φ_a)² - ½(∂_i φ_a)²] - V(P)

    V(P) = (μ/2) P² / (1 + κ P²)

    P = φ₀ φ₁ φ₂    (triple product)

Parameters at bimodal sweet spot (t=0.85):
- μ = -41.3 (negative = attractive binding)
- κ = 50.0  (regularization, prevents V → -∞)
- NO mass term (m=0 in Lagrangian; effective mass from dynamics)

Equations of motion:

    ∂²φ_a/∂t² = ∇²φ_a - ∂V/∂φ_a

    ∂V/∂φ_a = μ P (∂P/∂φ_a) / (1 + κP²)²

    ∂P/∂φ₀ = φ₁φ₂,  ∂P/∂φ₁ = φ₀φ₂,  ∂P/∂φ₂ = φ₀φ₁

### 2. Initial Condition — Propagating Helical Braid

Each field is a helical wave along z with per-field elliptical envelope:

    φ_a(x,y,z,t=0) = E_a(x,y) · cos(k·z + δ_a)

    ∂_t φ_a(x,y,z,t=0) = ω · E_a(x,y) · sin(k·z + δ_a)

where ω = √(k² + m²), and m here is the initialization mass (1.50)
used to set the dispersion relation (NOT in the Lagrangian).

**Per-field envelope** (the key symmetry-breaking structure):

    E_a(x,y) = A_a · G_a(x,y) + A_bg

    G_a(x,y) = exp(-r²_{eff,a} / (2R²_tube))

    r²_{eff,a} = (x'_a)²/(1+ε)² + (y'_a)²/(1-ε)²

    x'_a = (x - c^x_a)cos(θ_a) + (y - c^y_a)sin(θ_a)
    y'_a = -(x - c^x_a)sin(θ_a) + (y - c^y_a)cos(θ_a)

Per-field centers (strand displacement, R_disp = 0 at sweet spot):

    c^x_a = R_disp · cos(2πa/3)
    c^y_a = R_disp · sin(2πa/3)

Per-field ellipse angle (ell_rot = 0 at sweet spot → all same angle):

    θ_a = θ_0 + [2πa/3 if ell_rot > 0.5, else 0]

**Bimodal sweet spot parameters** (t=0.85 interpolation of Path A and B):

| Parameter | Value | Meaning |
|-----------|-------|---------|
| A₀ = A₁ = A₂ | 0.80 | Field amplitudes (equal) |
| δ₀ | 0.00 | Phase of field 0 (reference) |
| δ₁ | 3.00 | Phase of field 1 (= 0.15×0.00 + 0.85×3.53) |
| δ₂ | 4.43 | Phase of field 2 (= 0.15×1.67 + 0.85×4.92) |
| R_tube | 3.00 | Gaussian envelope width |
| ε (ellip) | 0.333 | Ellipticity (0=circular, 1=flat) |
| θ₀ (ell_ang) | 0.00 | Ellipse orientation |
| k | π/L | Axial wavenumber (k_fac=1.0, one half-twist) |
| A_bg | 0.00 | Background amplitude |
| R_disp | 0.00 | Strand displacement (not used at sweet spot) |
| ell_rot | 0.00 | Per-field rotation (not used at sweet spot) |
| m_init | 1.50 | Initialization mass (dispersion only) |
| μ | -41.3 | Triple-product coupling |
| κ | 50.0 | Regularization |
| λ_pw | 0.00 | Pairwise coupling (not used) |

The dispersion mass m_init=1.50 enters ONLY the initial velocity:
ω = √(k² + m²_init). It is NOT in the Lagrangian. The soliton mass
is entirely emergent from the field dynamics.

### 3. Boundary Conditions

- **z-direction**: Periodic (braid closes on itself)
- **x,y-directions**: Absorbing damping layer

    damp(r_⊥) = 1 - 0.98·f²    where f = (r_⊥ - 0.70L) / (0.95L - 0.70L)

    φ_a → damp · φ_a,  ∂_t φ_a → damp · ∂_t φ_a    for r_⊥ > 0.70L

### 4. Numerical Method

- **Grid**: N³ points on [-L, L]³, dx = 2L/(N-1), L=20
- **Time integration**: Velocity Verlet, dt = 0.20·dx
- **Laplacian**: 7-point stencil (6 neighbors + center), periodic in z
- **Boundary margin**: i,j ∈ [1, N-2] for force computation; acc=0 outside

Verlet step:
1. v += (dt/2)·a
2. φ += dt·v
3. Recompute forces a
4. v += (dt/2)·a
5. Apply damping

### 5. Diagnostic Observables

**Field concentration** (localization):

    fc = (Σ_{r_⊥ < R_core} |φ|²) / (Σ_all |φ|²),    R_core = 8

**Transverse quadrupole** (spin-2 gravity proxy):

    I_xx = Σ ρ(x)·x²,  I_yy = Σ ρ(x)·y²,  ρ = Σ_a φ²_a

    trans_l2 = |I_xx - I_yy| / (I_xx + I_yy)

This measures xy-plane anisotropy (NOT cylindrical elongation along z).

**Torsion flux** (EM proxy):

    ω_{xy} = ½(∂_x φ₁ - ∂_y φ₀)    (antisymmetric gradient tensor)

    Φ_T = ∫∫ ω_{xy} dx dy    (integrated over z=0 slice)

Also computed: ω_{xz} = ½(∂_x φ₂ - ∂_z φ₀), ω_{yz} = ½(∂_y φ₂ - ∂_z φ₁).
Report max(|Φ_xy|, |Φ_xz|, |Φ_yz|).

**Phase winding number** (topological invariant):

    W = (1/2π) Σ_{k=0}^{N-1} Δθ_k

    θ_k = atan2(φ₁(x_c, y_c, z_k), φ₀(x_c, y_c, z_k))

    Δθ_k = atan2(sin(θ_{k+1} - θ_k), cos(θ_{k+1} - θ_k))

measured along z-axis at domain center (x_c, y_c).

**Triple product peak**: |P|_max = max over grid of |φ₀ φ₁ φ₂|

### 6. The Bimodal Discovery

The bimodal configuration was found by interpolating between two
configurations identified by CMA-ES automated search:

**Path A** (gravity-optimized):
- δ = (0.00, 0.00, 1.67), ε = 0.80, μ = -29.7
- High ellipticity → strong xy-quadrupole, weak torsion

**Path B** (EM-optimized):
- δ = (0.00, 3.53, 4.92), ε = 0.25, μ = -43.4
- B-like phase structure → strong torsion, weak quadrupole

**Parameter interpolation** at mixing fraction t:

    p(t) = (1-t)·p_A + t·p_B    for each parameter

The optimum is at **t = 0.85** (85% Path B, 15% Path A).

The interpolation landscape is nonlinear: t=0.50 is a valley where
both metrics DROP, and t=0.85 is a ridge where both PEAK above
either pure control. This is a genuine synergy from the triple-product
coupling, not a linear averaging effect.

### 7. Validated Results (N=128, T=500)

| Config | trans_l2 | |Φ_T| | fc | |P| | winding |
|--------|----------|-------|------|------|---------|
| Pure A (t=0) | 0.083 | 0.092 | 0.519 | 0.000 | +1.0 |
| Pure B (t=1) | 0.100 | 0.512 | 0.958 | 0.477 | -1.0 |
| **Bimodal (t=0.85)** | **0.207** | **1.021** | **0.932** | **0.730** | **-1.0** |

The bimodal config exceeds Pure A in trans_l2 by 2.5× and Pure B in
torsion by 2.0×. Both pure controls degrade more at high resolution
than the bimodal config — the synergy strengthens with resolution.

### 8. Physical Interpretation

The three fields φ₀, φ₁, φ₂ form a propagating helical braid along z.
The field index a maps to spatial direction a (elastic interpretation):

    ε_{ij} = ½(∂_i φ_j + ∂_j φ_i)    (strain tensor → gravity)
    ω_{ij} = ½(∂_i φ_j - ∂_j φ_i)    (torsion tensor → EM)

The moderate ellipticity (ε=0.33) breaks xy-symmetry, creating a
quadrupolar strain pattern (spin-2). The B-like phase offsets
(δ ≈ 3.0, 4.4) create torsion circulation (spin-1). The triple-product
V(P) binds them into a single coherent structure.

The effective mass M = |E|/c² is entirely from field dynamics — no
mass parameter in the Lagrangian. The topological winding W=-1 is
conserved by the continuous dynamics with periodic z-boundary.

---

## V26-Dyn: Pulsating Dynamic Braid — Self-Concentrating Through Rotation (March 2026)

**Discovery**: A rotating+propagating braided soliton (DynAB, m=1, Ω=combined)
exhibits PULSATING self-concentration: fc oscillates between 0.07 (diffuse)
and 0.89 (concentrated) with period ~100 time units. The braid survives
all 500 time units with a persistent core at fc≈0.30 (2.5% of initial energy).

**Pulsation peaks**: fc = 0.70 (1st) → 0.89 (2nd) → 0.78 (3rd) → 0.35 (4th).
The pulsation is DAMPED but persistent. Energy stabilizes during concentrated
phases (self-sustaining periods).

**Key findings from the 7-variant scan**:
- m=0 always collapses (DynC, DynBC, DynAC, DynABC — all fatal)
- m=1 + rotation creates pulsating braids (DynB, DynAB)
- Propagation + rotation together give the strongest pulsation (DynAB)
- The static braid (V26 Phase 4) slowly unwinds; the dynamic braid pulsates

**The pulsation mechanism**: angular momentum from rotation creates a
gyroscopic effect. During the diffuse phase, the braid spreads. During
the concentrated phase, the angular momentum pulls it back together.
Each cycle radiates some energy, damping the pulsation. The remaining
core persists as a slowly-decaying dynamical soliton.

**Mass term is essential**: m² provides the confinement that prevents
runaway growth. No dynamical mechanism (propagation, rotation) can
substitute for it. The mass term is a NECESSARY PARAMETER, not emergent.

## V27-M4: MASSLESS Propagating Braid SURVIVES — Mass is Emergent (March 2026)

**BREAKTHROUGH**: A propagating helical braid with m=0 (NO mass term)
survives t=500 with |P|=1.742 — the HIGHEST triple-product coupling in
the entire project. The mass term was not just unnecessary — it was
HOLDING THE BRAID BACK.

**The two-regime structure**:
- Regime 1 (m≈1.0): positive energy, moderate |P|=0.574
- Dead zone (m=0.6-0.8): neither regime works
- Regime 2 (m≤0.4): negative energy, |P| INCREASES as m decreases

| m   | |P|(500) | E       | Survived? |
|-----|---------|---------|-----------|
| 1.0 | 0.574   | +281    | YES       |
| 0.8 | 0.003   | +137    | NO        |
| 0.6 | 0.000   | +31     | NO        |
| 0.4 | 0.589   | -2,871  | YES       |
| 0.2 | 1.105   | -6,381  | YES       |
| 0.1 | 1.591   | -7,518  | YES       |
| 0.0 | 1.742   | -7,915  | YES       |

**The Lagrangian with NO mass term**:

    L = Σ_a [½(∂_t φ_a)² - ½(∂_i φ_a)²] - (μ/2)P²/(1+κP²)

Just kinetic energy + triple product. The "mass" of the soliton is
ENTIRELY from the field dynamics: M = |E|/c² = 7915 code units.

**Why m=0 works with propagation but not static**: V26 showed m=0 static
braids collapse in ~6 time units. The propagating helical wave maintains
the three-field overlap (|P|) through continuous flow — the braid is
REMADE at each step. Without propagation, there's nothing to maintain
the overlap and the triple product drives the fields to infinity.

**V24 "self-reconstruction" hypothesis PARTIALLY confirmed**: propagation
maintains |P| at m=1 (Regime 1), but at m=0 (Regime 2) the topology
self-conserves regardless of motion. The V26 "unwinding" was amplitude
death from weak binding (μ=-20), not topological loss. At μ=-50, both
static and propagating braids maintain winding=-1.000 exactly.

## V27-M5: Winding Number Exactly Conserved — Intrinsic Topological Protection

**Discovery**: The phase winding number of the helical braid is EXACTLY
conserved (to machine precision) by the field dynamics. No Chern-Simons
term or additional topological coupling is needed.

Verified for: propagating braid, static braid, and CS-augmented braid.
All give winding=-1.000 at every timestep from t=0 to t=500.

**Mechanism**: periodic BC + smooth dynamics make the winding a topological
invariant that cannot change under continuous evolution. The Jacobian
density is a total divergence (d_m[cofactor]=0), so its integral is
exactly conserved.

## V27-M3c: Two Braids ATTRACT (Both Same and Opposite Twist)

**Discovery**: Two m=0 propagating braids at D=30 attract each other:
- Same twist: Δsep = -7.2 (stronger attraction)
- Opposite twist: Δsep = -5.0 (weaker attraction)

Both attract, with same-twist attracting MORE. This is gravity-like
(universal attraction regardless of "charge"), not EM-like (same charges
repel). The force is from triple-product tail overlap, not from torsion.

---

**DynA (propagation only, m=1) maintains braid structure (|P|=0.574 vs
static's 0.006 — 96× better!)**: Axial propagation at v_g=0.30c sustains
the three-field overlap that defines the braided soliton. The static braid
unwinds (|P|→0.006); the propagating braid HOLDS (|P|→0.574). This is the
"self-reconstruction through motion" hypothesis CONFIRMED: the helical
wave continuously flows through the braid, remaking the triple-product
coupling at each step. l=2 = 26.7% (vs 41.5% static, 0.07% oscillon).

**DynB (rotation only, m=1) is the most stable braid**: fc=0.71 at t=200,
rotation-speed-independent. The azimuthal phase structure (θ in cos(θ+kz))
provides stability beyond what propagation (DynA: fc=0.27) or the
combination (DynAB: fc=0.28) achieves.

**Complete V26 dynamic braid table (all 3D, N=128)**:

| Variant | m | Config | fc | l=2 | Status |
|---------|---|--------|-----|-----|--------|
| DynA | 1 | Propagate | 0.27 | 6.8% | Survives |
| **DynB** | **1** | **Rotate** | **0.71** | **7-12%** | **Most stable** |
| DynAB | 1 | Both | 0.28 | **18%** | Most aspherical |
| DynC | 0 | Propagate | 0.00 | — | Collapse |
| DynAC | 0 | Prop | 0.00 | — | Collapse |
| DynBC | 0 | Rotate | 0.10 | 15% | Barely |
| DynABC | 0 | All | ~0 | 11-29% | Collapse |

---

## V26 Phase 4: 41.5% Quadrupole from Braided Soliton — Highest Spin-2 Ever (March 2026)

**Discovery**: The twisted-tube braided soliton with mass=1 and the simplest
Lagrangian (NO torsion or strain coupling, just kinetic + mass + triple product)
produces **l=2 = 41.5%** — the highest quadrupole fraction in any SCP test.

**The key insight**: extra coupling terms (κ_T, κ_S) DESTROY the braid by
penalizing the very gradients that create the twist. The torsion term
κ_T ω² fights the braid because ω_{ij} IS the braid's defining structure.
Penalizing it unwinds the braid.

**The correct picture**: the torsion ω_{ij} is the OBSERVABLE (EM field),
not an energy cost. The strain ε_{ij} is the METRIC (gravity). Both emerge
from the same gradient tensor ∂_iφ_j without needing explicit coupling terms.

**Self-consistent metric**: STABLE at the baseline. The framework g^{ij} =
δ^{ij} - 2ε_{ij} works through t=500 without instability.

**Comparison of l=2 content across the project**:

| Configuration | l=2 fraction | Notes |
|--------------|-------------|-------|
| Oscillon (V25 Phase 2) | 0.07% | Spherical, breathing |
| Oscillon + tidal ε=0.01 (V25 Phase 6) | 26% | But during destruction |
| Binary oscillon (V25 Phase 5) | 2.7% | Merged, monopole dominates |
| **Braided soliton baseline (V26 Phase 4)** | **41.5%** | **Stable, self-consistent** |
| Braided + torsion κ_T=2 | <5% | Torsion unwinds the braid |

---

## V26: Braided Solitons — First Non-Breathing Aspherical Configuration (March 2026)

**Discovery 1: Borromean rings are the most stable braided topology.**
Three interlocking rings (one per field, in orthogonal planes) survive
t=500 with fc=0.92 and |P| growing 70× (from 0.006 to 0.39). The triple
product coupling draws the rings together, tightening the braid.

**Discovery 2: Twisted tube without mass term is non-breathing.**
Mode 2 (triple product only, no m²φ²) does NOT breathe — the DFT of
ρ(center,t) shows ω=0.07 (essentially DC) vs ω=1.37 for the oscillon.
This is the FIRST non-breathing configuration in the SCP project.

**Discovery 3: Non-breathing configuration has 97× more l=2 content.**
The twisted tube (mode 2) has l=2 = 6.76% vs the oscillon's 0.07%.
The aspherical braid geometry naturally produces quadrupolar strain.

**Discovery 4: Mass term m²φ² destroys braided configurations.**
Mode 0 (twisted tube + full couplings including m²) disperses immediately.
Mode 2 (same tube, no mass term) survives. The mass term pushes fields
toward zero, unwinding any spatial structure. This confirms: braided
solitons require a DIFFERENT binding mechanism than the mass gap.

**Key negative: braid topology alone doesn't prevent breathing.**
The Borromean rings (mode 1, most stable) breathe at ω=1.37, the SAME
frequency as the oscillon. The triple product coupling creates breathing
regardless of the topological structure. Non-breathing requires removing
the mass term, which weakens localization.

**The Lagrangian for the non-breathing mode (mode 2):**

    L = Σ_a [½(∂_t φ_a)² - ½(∂_i φ_a)²] - (μ/2)P²/(1+κP²)

No mass term, no pairwise, no elastic. Pure kinetic + triple product.
The braided initial condition provides the spatial structure; the triple
product provides the binding. This is the SIMPLEST Lagrangian that
supports a non-breathing, aspherical, topological field configuration.

---

## V25: 3D Elastic Gravity — Tidal Quadrupole Resonance (March 2026)

**Achievement**: The 3D elastic oscillon with self-consistent metric works.
The strain tensor is computable, causal, and has spin-2 structure. A tidal
quadrupolar resonance exists at Ω_T ≈ 0.1 with 26% l=2 content and
pole/equator suppression ratio 0.65.

**The V25 Lagrangian**:

    L = Σ_a [½(∂_t φ_a)² - ½(∂_i φ_a)²] - ½m²Σφ_a²
        - (μ/2)P²/(1+κP²)
        - λ_pw(φ₁φ₂+φ₂φ₃+φ₃φ₁)
        - ½η(∂_i φ_j)(∂_j φ_i)
        - ½λ_L(∂_i φ_i)²

    Self-consistent metric: h_{ij} = α_g(∂_i φ_j + ∂_j φ_i)
    Backreaction: (δ_{ij}+h_{ij})∂_i∂_j replaces ∇²; m²(1+½tr(h)) replaces m²

**Key V25 Results:**

| Finding | Value | Significance |
|---------|-------|-------------|
| 3D elastic oscillon | Survives at η=0.1, λ_L=0.1 | First 3D elastic soliton |
| Self-consistent metric | Converges at α_g=0.001 | ds=c·dt built in |
| Tidal Love number | k₂ ≈ 100-1000 | Very soft (no shear rigidity) |
| Quadrupole resonance | 26% l=2 at Ω_T=0.1 | First l=2 dominance > 25% |
| Pole suppression | h₊_pole/h₊_eq = 0.65 | First angular spin-2 signature |
| Binary GW | 97% l=0 (monopole) | Breathing overwhelms orbital |

**The fundamental obstacle**: The oscillon is a breathing scalar object.
Its monopolar (l=0) radiation from spherical breathing overwhelms any
quadrupolar (l=2) signal from orbital or tidal dynamics. The 26% l=2
occurs only during oscillon destruction (ε_T=0.01), not in steady state.

**What spin-2 gravity requires**: A non-breathing ground state (static
soliton like a Skyrmion), or a mechanism that forbids l=0 radiation
while allowing l=2, or a rigid object (high stiffness, low k₂) that
deforms quadrupolarly instead of dispersing.

---

## V24 Complete Theory: Three-Field Oscillon with Forces and Causal Gravity

### The Full Lagrangian (as established by V24)

    L = Σ_a [½(∂_t φ_a)² - ½(∂_x φ_a)²] - ½m²Σ_a φ_a²
        - V(P) - λ(φ₁φ₂ + φ₂φ₃ + φ₃φ₁)
        + ½(∂_t Φ)² - ½c²(∂_x Φ)²

    V(P) = (μ/2)P²/(1+κP²),  P = φ₁φ₂φ₃

    Source:   □Φ = ∂²Φ/∂t² - c²∂²Φ/∂x² = -α·ρ
    where ρ = Σ_a [½(∂_t φ_a)² + ½(∂_x φ_a)² + ½m²φ_a²] + V + V_pw

    Backreaction:  m²_eff(x) = m²(1 + 2Φ),  c²_eff(x) = 1 + 4Φ

### Components and Their Roles

| Term | Role | Parameters | Discovery |
|------|------|-----------|-----------|
| m²φ_a² | Mass gap → oscillon protection | m=1.0 | v21 |
| V(P) = (μ/2)P²/(1+κP²) | Triple-product binding → matter | μ=-20, κ=20 | v21 |
| λ(φ₁φ₂+φ₂φ₃+φ₃φ₁) | Pairwise coupling → mass splitting, forces | λ=0.5 (lattice) | V24-ME |
| □Φ = -αρ | Causal gravity → 1/r potential at c | α~1e-4 | V24-L4 |
| m²(1+2Φ), c²(1+4Φ) | Metric backreaction | from α | V24-F |

### Mass Spectrum

The pairwise coupling splits the three-field spectrum:

    Symmetric mode:      m²_S = m² + 2λ     (oscillon lives here)
    Antisymmetric modes: m²_A = m² - λ       (Proca mediator, 2-fold)

At λ=0.5: m_S = 1.414, m_A = 0.707 (range 1/m_A = 1.41)
At λ=0.99: m_S = 1.726, m_A = 0.100 (range 10)
At λ=0.9999: m_S = 1.732, m_A = 0.010 (range 100)

The oscillon frequency tracks m_S: ω ≈ 0.95·m_S (3-5% gap margin).

---

## Discovery 1: Causal Gravity via Wave Equation (V24-L1, L4)

**Equation**: □Φ = ∂²Φ/∂t² - c²∂²Φ/∂x² = -αρ

**Verified**: Signal at x=50 arrives at t=50.5±0.5 (L4) and t=51.5 (L1).
Speed = c to 1-3% accuracy. Exact causality with no free parameters.

**Static limit**: ⟨Φ(0)⟩ = -4.9e-3 to -9.1e-3 (78-100% of Poisson,
depending on implementation). Correct spatial profile ∝ |x|.

**Gravitational radiation**: Φ oscillates at 2ω ≈ 1.75 (from breathing ρ).
This IS physical gravitational radiation. The DC (static) component gives
the Newtonian potential; the AC component is gravitational waves.

**Implementation**: Evolve Φ as a 4th field array alongside φ₁,φ₂,φ₃
using Velocity Verlet. Same O(Nx) cost as the matter fields. Absorbing
boundary conditions for Φ.

---

## Discovery 2: Stable Oscillon Lattice from Pairwise Coupling (V24-E)

**Equation**: V_pw = λ(φ₁φ₂ + φ₂φ₃ + φ₃φ₁),  force: -λ(φ_b + φ_c)

At λ=0.5: lattice stability improves from max_u/d=4.0 (melting) to 0.7
(marginally stable). This is a 10⁴× improvement in effective well depth.

**Phonon spectrum**: Real linear dispersion with c_s ≈ 0.9.
The symmetric (acoustic) branch: ω = c_s·k (gapless, compression).
The antisymmetric (optical) branch: ω ≈ 0.51 (gapped, Proca-like).

**Antisymmetric mode propagation**: A localized antisymmetric perturbation
(φ₁-φ₂) at one oscillon PROPAGATES to all others AND is AMPLIFIED
(1.1-1.5× source). Parametric pumping from the breathing mode.

---

## Discovery 3: Unlimited Proca Range (V24-ME, P1, P5)

The antisymmetric Proca mass m_A = √(m²-λ) can be tuned to near zero
by pushing λ → m². Tested to λ=0.99999 (m_A=0.003, range=316).
**No stability wall exists.** The oscillon survives even at λ=1.05
(past the tachyonic boundary), stabilized by the triple-product Hessian.

**Multi-scale forces (V24-S3)**: Asymmetric pairwise couplings
(λ₁₂≠λ₂₃) produce THREE distinct eigenvalues → three ranges → two
simultaneous force scales with 13× spread.

---

## Discovery 4: Newton's Law and Inertial Deformation (V24-B)

**Verified**: a = F/M to 4% accuracy for the oscillon, using momentum
injection at small F (3e-6 to 1e-5).

**Deformation**: QUADRATIC in F (tidal, second-order). The trailing
side stretches ~20%; the leading side is unchanged. This is the
spin-2 signature of inertial deformation (in 1D: dipolar asymmetry).

**Critical force**: F_crit ≈ 2e-5 (v_max ≈ 0.05c). Above this the
oscillon is destroyed.

---

## Discovery 5: Self-Consistent Gravitational Potential (V24-F)

At α=1e-4: Φ(0) = -0.006 (0.6% well depth), oscillon perfectly stable
(fc=0.999, E retained 98.2%, ω=0.864).

At α=1e-3: Φ(0) = -0.038 (3.8% well depth), oscillon survives but
degraded (fc=0.80, ω shifted to 0.954 from gravitational blueshift).

Critical α ≈ 0.005: above this, c²_eff = 1+4Φ goes negative → tachyonic.

The self-consistency loop converges: ρ → Φ → modified dynamics → new ρ
→ same Φ (within numerical precision).

---

## Discovery 6: Z₂ Internal Symmetry (V24-180)

The transformation φ₃ → -φ₃ maps 0° state (φ,φ,φ) to 180° state
(φ,φ,-φ). Since V depends on P², and P → -P under the flip,
the Lagrangian is exactly invariant.

**Verified**: dE/dt identical to 6 significant figures. Spectra
bit-for-bit identical. Both states stable against perturbation.
This is the oscillon's simplest internal quantum number — a "charge."

---

## Discovery 7: UDD > UUD Mass Ordering (V24-UUD/UDD)

With per-field masses m_U ≠ m_D and pairwise coupling λ=0.30:
UDD (neutron analog) retains 27-187% MORE energy than UUD (proton).
Matches the real M_n > M_p ordering.

**Mechanism**: DD-antisymmetric eigenvalue m_D²-λ is lighter than
UU-antisymmetric m_U²-λ → slower oscillation → less radiation → heavier.

In the strongly-bound 180° regime: ordering REVERSES (bare mass dominates).

---

## Discovery 8: Vacuum Equation of State (V24-A)

The oscillon's time-averaged stress tensor:
- Core (|x|<5): P/ρ ≈ 0 (pressureless dust, virial balance)
- Vacuum (|x|>30): P/ρ ≈ +0.8 (radiation, from leaked harmonics)
- Integrated: ∫P/∫ρ ≈ 5e-5 (total pressure vanishes — virial theorem)

The oscillon IS conventional matter (pressureless dust at macroscopic
level), not exotic vacuum energy.

---

## Discovery 9: Spin-0 Mediator (V24-S6)

The Proca mediator from pairwise coupling is definitively SPIN-0.
2D multipole analysis: l=0 (monopole) = 100%, l=2 (quadrupole) < 0.01%.
Adding cross-gradient η=0.5: l=2 rises to 10.9% but as anisotropic
scattering, not spin-2 radiation.

**Representation-theory obstruction**: scalar fields transform trivially
under rotations. No coupling can change spin-0 to spin-2. Spin-2
requires a fundamental rank-2 tensor field.

---

## Discovery 10: Explicit Gauge Coupling Kills Oscillons (V24-MA, MB)

U(1) gauge coupling (covariant derivative D_μ = ∂_μ - ieA_μ) destroys
the oscillon at e ≥ 0.1. The e²A²|Ψ|² term shifts ω above the mass gap.

**Root cause**: the triple-product phase-locking (which makes the oscillon
work) is fundamentally incompatible with U(1) charge (which requires a
phase DIFFERENCE between fields).

---

## Negative Results and Lessons

| Approach | Result | Lesson |
|----------|--------|--------|
| 120° phase separation | P≈0, no binding | Triple product vanishes at 120° |
| Boosted oscillon | v=0 most stable | Forward beaming > backward Doppler |
| Thermal bath | Noise destroys at T=0.001 | Binding too weak per mode |
| Confining cusp V=-σ|P| | Near-gapless mode is INSTABILITY | Cusp = saddle point, not minimum |
| Vortex (2D) | Mass kills 1/r tail in ~50 t.u. | Mass gap prevents extended structures |
| Dissipative soliton | No stable balance point | Polynomial gain always blows up |
| Local field metric | Causal but LOCAL | Self-interaction, not gravity |
| Telegraph equation | Φ too small at large γ | Overdamped = slow convergence |

---

## V23 Results Referenced

- V23-D: Inter-oscillon equilibrium at D≈15.7, phonon speed c_s=0.40
  (matches spring constant prediction 0.44 within 9%)
- V23-A: Hessian splitting — antisymmetric mode HARDENED by μ<0
  (λ_anti = +|μ|f⁴/(1+κf⁶)² > 0 always)
- V23 complete: all single-oscillon gravity paths closed

---

## Open Questions for V25

1. **Spin-2**: The Proca mediator is spin-0. How to get spin-2 for gravity?
   The lattice deformation (Combo 2+5 optical branch at ω≈0.51) is the
   closest, but it's gapped. Need gapless spin-2 → rank-2 tensor field.

2. **Coupling constant α**: What determines the gravitational coupling?
   Test F showed α < 0.005 for stability. S5 showed no self-consistency
   mechanism fixes α. It remains a free parameter.

3. **3D**: All V24 results are in 1D. The 3D oscillon survives pairwise
   coupling (V24-S4), but the lattice, phonons, and causal gravity have
   not been tested in 3D.

4. **Gravitational radiation vs static potential**: The wave equation
   □Φ = -αρ produces BOTH. Can the ratio be controlled? In GR, the
   static component dominates for slow sources (v << c).

5. **Confinement**: The oscillons still leak (dE/dt ~ -10⁻⁵). True
   stability requires topological protection or a confining mechanism
   that doesn't create saddle-point instabilities.
