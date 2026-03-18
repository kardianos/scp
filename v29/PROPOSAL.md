# V29: Validation and Stress-Testing the Bimodal Braid

Starting point: the V28 bimodal sweet spot (t=0.85 interpolation).
Goal: answer the open questions raised by external review and our own analysis.

## Directory Structure

```
v29/
├── PROPOSAL.md              ← this file
├── RESULTS.md               ← aggregate results (written after tests complete)
├── src/
│   └── braid_core.h         ← shared PDE solver, braid init, diagnostics
│
├── T1_thermal/              ← thermal bath equilibrium
│   ├── PROPOSAL.md
│   ├── src/t1.c             ← #include "../src/braid_core.h"
│   ├── data/
│   └── RESULTS.md
│
├── T2_minit/                ← initialization mass independence
│   ├── PROPOSAL.md
│   ├── src/t2.c
│   ├── data/
│   └── RESULTS.md
│
├── T3_radiation/            ← late-time radiation rate measurement
│   ├── PROPOSAL.md
│   ├── src/t3.c
│   ├── data/
│   └── RESULTS.md
│
├── T4_fragility/            ← topological fragility under perturbation
│   ├── PROPOSAL.md
│   ├── src/t4.c
│   ├── data/
│   └── RESULTS.md
│
├── T5_metric/               ← self-consistent metric coupling (own solver)
│   ├── PROPOSAL.md
│   ├── src/t5.c
│   ├── data/
│   └── RESULTS.md
│
├── T6_universality/         ← bimodal synergy across (mu, kappa) pairs
│   ├── PROPOSAL.md
│   ├── src/t6.c
│   ├── data/
│   └── RESULTS.md
│
├── T7_twobraid/             ← two-braid interaction and force law
│   ├── PROPOSAL.md
│   ├── src/t7.c
│   ├── data/
│   └── RESULTS.md
│
├── T8_boost/                ← Lorentz boost validation
│   ├── PROPOSAL.md
│   ├── src/t8.c
│   ├── data/
│   └── RESULTS.md
│
└── T9_substrate/            ← the model behind the model
    ├── PROPOSAL.md
    ├── T9A_lattice/         ← lattice reverse-engineering
    │   ├── src/t9a.c
    │   └── data/
    ├── T9C_expansion/       ← expansion cosmology simulation
    │   ├── src/t9c.c
    │   └── data/
    └── RESULTS.md
```

### Shared Code: `src/braid_core.h`

Single header containing:
- Grid type, alloc, free, setup
- Bimodal sweet spot parameters (PATH_A, PATH_B, BIMODAL_085)
- Braid initialization (init_braid, with configurable m_init override)
- Force computation (compute_forces, OMP-parallelized)
- Time stepping (verlet_step, verlet_finish)
- Damping (apply_damping — called optionally by each test)
- Diagnostics (compute_diagnostics, compute_winding, check_blowup)
- Result type and printing

Each test's .c file includes this header and implements its own main()
with test-specific logic (BC choice, perturbations, scans, etc.).

### Build Convention

Each test builds from its own directory:
```
cd v29/T2_minit && gcc -O3 -fopenmp -o t2 src/t2.c -lm
```
The #include path "../src/braid_core.h" resolves relative to the test's src/.
Alternatively from v29 root: gcc -O3 -fopenmp -I src -o T2_minit/t2 T2_minit/src/t2.c -lm

## Test List

### T1: Thermal Bath Equilibrium
**Question**: Does the braid reach dynamic equilibrium in a radiation bath,
absorbing as much as it radiates?

**Motivation**: The criticism that m=0 solitons "bleed energy and die" assumes
empty space. In a universe of these fields, radiation is not loss — it's the
interaction mechanism. Real particles exist in thermal equilibrium, not vacuum.

**Method**:
- Initialize bimodal braid (t=0.85 params) in a box with low-amplitude
  thermal noise: δφ_a ~ A_noise × gaussian_random at each grid point
- Periodic BC in ALL directions (no absorbing layer)
- Scan A_noise ∈ {0.001, 0.01, 0.05, 0.1}
- Run to T=1000+
- Measure: does the braid stabilize? What is the equilibrium energy?
  Does the radiation spectrum approach Planck (T⁴)?
- Diagnostics: E(t), fc(t), |P|(t), radiation spectrum via FFT of
  φ at far-field points, spectral energy density ε(ω)

**Key observable**: If dE/dt → 0 at late times with noise present, the
soliton is in dynamic equilibrium. If dE/dt < 0 always, it truly bleeds.

---

### T2: Initialization Mass Independence
**Question**: Does the bimodal synergy survive with m_init = 0?

**Motivation**: Gemini's strongest critique — m_init=1.5 injects excess
energy. If the braid works with ω = k (m_init=0), the emergent mass
claim is solid.

**Method**:
- Run bimodal t=0.85 with m_init = 0, 0.5, 1.0, 1.5 (current), 2.0
- Same grid (N=128, T=500), absorbing BC
- Compare all metrics: trans_l2, torsion, fc, |P|, winding

---

### T3: Late-Time Radiation Rate
**Question**: How fast does the braid radiate in open space?

**Motivation**: Is it oscillon-like (power-law decay, quasi-stable for
astronomical times) or genuinely unstable (rapid energy loss)?

**Method**:
- Large domain L=60, N=192 (dx ≈ 0.63, adequate resolution)
- No damping layer at all, absorbing at r > 0.90L only
- Run to T=2000
- Track E_core(t), E_far(t), dE_core/dt at late times (T > 500)
- Fit: dE/dt ~ -C/t^α (oscillon) or -C·exp(-γt) (resonance decay)?
- Also: is radiated power ~ T⁴ for some effective temperature?

---

### T4: Topological Fragility — Critical Perturbation
**Question**: How strong a perturbation can the braid survive without
losing its winding number?

**Motivation**: Gemini correctly notes π₁(R³)=0 — the protection is
dynamical, not fundamental. How robust is it?

**Method**:
- Initialize bimodal braid, let it settle for T=200
- Apply localized Gaussian perturbation at the core:
  φ_a → φ_a + ε × exp(-r²/σ²) × random_direction
- Scan ε ∈ {0.1, 0.3, 0.5, 1.0, 2.0, 5.0} with σ = R_tube
- Continue evolution to T=500
- Track: winding(t), fc(t), |P|(t)
- Find ε_crit where winding changes
- Compare ε_crit to the field amplitude (~0.8) — if ε_crit ≫ A₀,
  the protection is effectively topological

---

### T5: Self-Consistent Metric Coupling
**Question**: Does the bimodal braid produce spin-2 radiation when
coupled to the metric it creates?

**Motivation**: The strain/torsion are currently "proxies." V25 showed
self-consistent coupling works for oscillons. The bimodal braid (non-
breathing, aspherical) is the best candidate for true spin-2 radiation.

**Method**:
- Implement V25-style backreaction on the V28 braid:
  g_ij = δ_ij + α_g × (∂_i φ_j + ∂_j φ_i)
  Replace ∇² with g^{ij}∂_i∂_j in the EOM
- Start with small α_g (0.0001, 0.001) and increase
- Measure: l=2 content of the metric perturbation h_ij at far field
- Compare to l=0 (monopole) — does the braid produce MORE l=2 than
  the oscillon did (V25: l=2 was 0.07% in steady state)?

---

### T6: Universality — Different (μ, κ) Pairs
**Question**: Is the bimodal synergy specific to μ=-41.3, κ=50, or does
it persist across parameter space?

**Motivation**: If it's universal, it's a generic feature of triple-product
braids. If it's fine-tuned, it's a numerical curiosity.

**Method**:
- Fix the bimodal GEOMETRY (phases, ellip from t=0.85) but vary dynamics:
  (μ, κ) ∈ {(-20,20), (-30,30), (-41,50), (-60,60), (-80,80), (-100,100)}
- Also vary the ratio μ/κ independently
- Run at N=80, T=200
- Check: does the synergy (beating both pure A and B) persist?

---

### T7: Two-Braid Interaction at Bimodal Sweet Spot
**Question**: Do two bimodal braids attract? What is the force law?

**Motivation**: V27-M3c showed attraction between simple braids. The bimodal
braid has both quadrupolar strain AND torsion — the interaction should be
richer (both "gravitational" and "electromagnetic" components).

**Method**:
- Two bimodal braids separated by D=30 along x-axis
- Same twist, opposite twist, and perpendicular orientations
- N=256 or N=192 at L=40 for adequate resolution
- Track separation vs time, decompose force into multipoles
- Is the force 1/r² (Coulomb/Newton) or shorter-range (Yukawa)?

---

### T8: Boost Test — Lorentz Contraction
**Question**: Does a boosted braid Lorentz-contract correctly?

**Motivation**: The Lagrangian is Lorentz-invariant, so a boosted soliton
should contract. This validates the relativistic structure.

**Method**:
- Initialize bimodal braid with a Lorentz boost: φ(x) → φ(γ(x-vt))
  with corresponding velocity field
- Test v = 0.1c, 0.3c, 0.5c
- Measure the transverse and longitudinal profiles
- Check: does the aspect ratio match γ?

---

### T9: The Model Behind the Model — Substrate Derivation
**Question**: What underlying physical structure, when coarse-grained or
rapidly expanded from a dense state, naturally produces L = ½(∂φ)² - V(P)?

**Motivation**: The Lagrangian was postulated. But its symmetries are highly
specific: S₃ permutation of three fields, Z₂ per field, and the triple
product P = φ₀φ₁φ₂ is a VOLUME FORM (det of the displacement gradient).
This means V(P) penalizes local volume change. These are the symmetries
of a 3D elastic medium. The question is: what discrete or continuous
microstructure has these symmetries as its natural low-energy limit?

**Approach — Three parallel tracks**:

**Track A: Lattice reverse-engineering (computational)**
- Define candidate discrete structures:
  (a) Cubic lattice with nearest-neighbor springs
  (b) FCC/HCP close-packed lattice
  (c) Random network (Voronoi/Delaunay)
  (d) Oriented graph with 3-body interactions
  (e) Spin-ice / vertex model with 3 species
- For each, define the most general nearest-neighbor interaction Hamiltonian
  compatible with the lattice symmetry
- Coarse-grain analytically (Taylor expand displacements u_i = φ_i(x))
  or numerically (block-spin RG)
- Check: does the continuum limit contain a term proportional to
  (u₁u₂u₃)² or det(∂u)²? If so, that lattice IS the substrate.

The key constraint: the triple product arises from THREE-BODY interactions
(not pairwise). So the substrate must have irreducible 3-body couplings.
In lattice terms: the energy depends on the VOLUME of a simplex (tetrahedron),
not just the lengths of its edges.

**Track B: Symmetry / algebraic classification (analytical)**
- The Lagrangian symmetry group is S₃ × (Z₂)³ acting on (φ₀, φ₁, φ₂)
- Classify all 3D lattice symmetry groups that reduce to this in the
  continuum limit
- The triple product P = det(F) where F_ij = δ_ij + ∂_j φ_i is the
  deformation gradient. This is a standard object in finite-strain elasticity.
- The potential V(det F) = (μ/2)(det F - 1)²/(1 + κ(det F - 1)²) is a
  specific volumetric strain energy function
- Question: what MICROSCOPIC interaction (pair potential, many-body potential,
  bond-angle potential) produces this macroscopic strain energy?
- Known results from materials science: Cauchy-Born rule gives the
  connection between lattice and continuum. Apply it in reverse.

**Track C: Expansion cosmology (physical scenario)**
- Start with a dense, isotropic "pre-field" state: uniform energy density ρ₀
  in a small volume, with 3 internal degrees of freedom per point
- Rapidly expand (inflate) the volume by factor α³
- As the system cools, the 3 DOF per point "freeze" into displacement fields
- The frozen-in structure contains topological defects from the cooling:
  braids = dislocations, winding = Burgers vector
- Simulate this process:
  1. Dense random initial state on small lattice
  2. Adiabatic expansion (increase lattice spacing, rescale interactions)
  3. Monitor: do braided structures spontaneously form?
  4. Coarse-grain the final state: does it match the V28 braid?

This is analogous to the Kibble-Zurek mechanism for cosmic string formation
during symmetry-breaking phase transitions, but for a 3D elastic medium
instead of a gauge field.

**Track D: Information-theoretic minimum (theoretical)**
- What is the SIMPLEST discrete dynamical system that:
  (i) Has 3 equivalent, interchangeable DOF per site
  (ii) Has a conserved winding number under smooth evolution
  (iii) Supports localized propagating structures
  (iv) Has both symmetric (strain) and antisymmetric (torsion) gradients
- This is a cellular automaton / lattice gas approach: enumerate rules,
  check which produce soliton-like behavior
- The answer tells us the MINIMUM STRUCTURE needed — the irreducible
  substrate

**Expected outcome**: Track A is the most likely to produce concrete results.
If a specific lattice + 3-body interaction reproduces V(P) under coarse-
graining, that lattice is the candidate for "what space is made of." Track C
tests whether braids form spontaneously during expansion (they should, by
Kibble-Zurek). Track B provides the mathematical classification. Track D
gives the information-theoretic lower bound.

**Computational approach for Track A**:
- For each lattice type, parameterize the 3-body interaction:
  H = Σ_{triangles} f(V_tri) where V_tri is the signed volume of the
  triangle formed by 3 displaced nodes
- Taylor expand: V_tri ≈ V₀(1 + ∂_i u_j ε_{ijk} + ...) → det(δ+∂u)
- The continuum limit of Σ f(V_tri) naturally produces ∫ V(det(∂u)) d³x
- The specific form V = (μ/2)P²/(1+κP²) constrains f to be a saturating
  quadratic in volume deviation

**Computational approach for Track C**:
- 3D lattice, N=64³, start at high density (small lattice spacing a₀)
- Each node has position r_i and 3 internal DOF (spin/color/displacement)
- Nearest-neighbor + triangle interactions
- Slowly increase a (expansion): a(t) = a₀ × exp(H×t), H = expansion rate
- Cool: reduce kinetic energy (thermostat at decreasing T)
- At final state: identify topological defects, measure field profiles
- Compare defect structure to V28 braid

---

## Priority Order

**Tier 1 — Answer the critics (fast, targeted)**:
1. **T2** (m_init independence) — quick, directly answers strongest critique
2. **T4** (topological fragility) — quantifies the protection strength
3. **T3** (radiation rate) — characterizes the quasi-stability

**Tier 2 — Test the dynamic equilibrium hypothesis (the big idea)**:
4. **T1** (thermal bath) — does the braid reach steady state in a radiation bath?
5. **T6** (universality) — is bimodal synergy generic or fine-tuned?

**Tier 3 — Push toward real physics**:
6. **T5** (self-consistent metric) — the path to real spin-2
7. **T7** (two-braid interaction) — force law characterization
8. **T8** (Lorentz boost) — relativistic validation

**Tier 4 — The deep question**:
9. **T9** (substrate derivation) — what is the model a model of?
   Track A (lattice reverse-engineering) first, Track C (expansion) second

## Bimodal Sweet Spot Parameters (from V28)

For all tests, the baseline configuration is the t=0.85 interpolation:

    L = (Σ_a [½(∂_t φ_a)² - ½(∂_i φ_a)²]) - (μ/2)P²/(1+κP²)
    P = φ₀φ₁φ₂
    μ = -41.3, κ = 50.0

    Init: A = 0.8, δ = (0, 3.00, 4.43), R_tube = 3.0, ε = 0.333
    k = π/L, m_init = 1.50 (or 0 for T2), periodic z, absorbing xy

    Validated at N=128, T=500:
    trans_l2 = 0.207, torsion = 1.021, fc = 0.932, |P| = 0.730, W = -1.000
