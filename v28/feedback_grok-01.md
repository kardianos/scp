**The document is an internal research log ("SCP Discoveries") chronicling a multi-year numerical exploration (V24–V28, dated up to March 2026) of classical 3D scalar field theory solitons.** It aims to construct stable, non-breathing, topologically protected "particles" from three real scalar fields φ₀,₁,₂ that exhibit emergent features interpretable as gravity (spin-2 quadrupole/strain) and electromagnetism (spin-1 torsion), plus pairwise forces, causal gravity, etc., in earlier iterations. The progression builds from simple oscillons (V24: mass term + triple-product binding + auxiliary gravity field Φ) through elastic gravity and lattice effects (V25), braided topologies (V26: Borromean rings or twisted tubes), massless propagation/rotation dynamics (V27: mass emergent, winding exactly conserved, braids attract), to the capstone **V28 bimodal propagating helical braid**.

V28 uses the **simplest Lagrangian** (massless, no pairwise/elastic terms):

L = Σ_a [½(∂_t φ_a)² − ½(∂_i φ_a)²] − V(P),  
V(P) = (μ/2) P² / (1 + κ P²)  (μ = −41.3 < 0 attractive, κ = 50 regularizes),  
P = φ₀ φ₁ φ₂.

**Initial conditions** are a propagating helical wave (one half-twist, k = π/L, periodic in z) with per-field elliptical Gaussian envelopes (ε = 0.333, equal amplitudes A = 0.80, tuned phases δ ≈ (0, 3.0, 4.43) from CMA-ES-optimized "Path A" gravity-heavy vs. "Path B" torsion-heavy interpolation at t = 0.85). Velocity Verlet on a damped grid yields a stable soliton (fc ≈ 0.93, |P|_max = 0.73) with **bimodal synergy**: transverse quadrupole trans_l2 = 0.207 (2.5× Path A) and torsion flux |Φ_T| = 1.021 (2× Path B), winding W = −1 exactly conserved. Diagnostics treat the fields as "displacements" mapping to elastic strain ε_ij (gravity proxy via ρ = Σ φ_a² quadrupole moment) and torsion ω_ij (EM proxy via integrated antisymmetric gradients). Braids attract (gravity-like, stronger for same twist).

**Key claims**: Mass is fully emergent from dynamics (no m² term); topology protects the structure; nonlinear triple-product binding + helical flow enables "self-reconstruction"; bimodal interpolation produces genuine non-linear synergy stronger than linear averaging.

### Flaws in the V28 Model (Conceptual, Theoretical, and Numerical)

V28 represents a clever numerical optimization of a toy model, but it has several fundamental limitations that prevent it from being a viable physical theory of particles or unification:

1. **Proxies, not emergent fields**: The "spin-2 gravity" is merely the normalized transverse quadrupole |I_xx − I_yy|/(I_xx + I_yy) of the energy-density moment ρ = Σ φ_a² — a shape descriptor, **not** a dynamical metric or curvature. No Einstein equations, back-reaction, or gravitational waves appear in the Lagrangian or EOM. Similarly, "EM" is a post-hoc integrated flux Φ_T from ω_xy = ½(∂_x φ₁ − ∂_y φ₀) etc. — no vector potential, Maxwell equations, Lorentz force, or gauge invariance. Earlier versions (V24–V25) had an explicit auxiliary gravity field Φ or elastic terms; V28 drops them for simplicity, leaving only interpretive analogies. The elastic/torsion mapping (φ_a → "displacements" → ε_ij/ω_ij) is ad-hoc continuum-mechanics analogy, not derived.

2. **Ad-hoc tuning and lack of robustness**: Parameters (μ, κ, δ_a, ε, etc.) are heavily optimized via CMA-ES and nonlinear interpolation (t = 0.85 "sweet spot"). The synergy is real in the numerics but specific to these proxies and grid (N=128, L=20, absorbing layers). Higher resolution strengthens it, but this could be convergence to a tuned attractor rather than generic. The initial helical braid is imposed by hand (not vacuum-emergent); static versions collapse (V26/V27), and survival relies on propagation "remaking" the overlap. The preferred z-axis breaks full isotropy/Lorentz invariance in the setup.

3. **Potential and energy issues**: With μ < 0, V ≈ (μ/2) P² near zero (tachyonic product mode) but saturates; the massless regime (V27) yields large negative total energies (E ≈ −7915 code units at m=0). While the potential is bounded below, long-term stability in open space (beyond periodic-z + damping) is questionable — radiation, decay, or blow-up possible without stricter topological protection. The "emergent mass" M = |E|/c² works numerically but lacks a deeper mechanism tying to positive-energy theorems or GR.

4. **Classical limitations and missing physics**: Purely classical, non-quantized, no fermions/spinors, no true gauge bosons or interactions (earlier explicit gauge coupling "kills" oscillons, V24). Attraction between braids is shown but not a full 1/r force law or Newtonian limit. No quantum statistics, renormalization, or Standard-Model embedding. The winding W = −1 is conserved only due to periodic BC + smooth evolution — a mathematical identity, not a deep gauge principle.

5. **Numerical artifacts possible**: Damping layer, Verlet integrator, 7-point Laplacian, and R_core = 8 cutoff for fc are pragmatic but introduce energy loss and boundary effects. The "bimodal valley/ridge" landscape could be optimization-specific rather than universal.

In short, V28 is an impressive demonstration of a stable, non-breathing, topologically protected classical soliton with controllable multipole moments via braid geometry, but it is **analog modeling** at best — not a fundamental unification or particle description.

### Existing Literature Search Results

No exact match exists for this specific setup: three independent real scalar fields, triple-product potential V(P = φ₀φ₁φ₂), propagating helical braid initial condition, massless Lagrangian, and post-hoc strain/torsion proxies for gravity/EM. Searches across arXiv, general web, and specific terms ("braided soliton scalar", "helical braid soliton", "triple product potential" three scalars, "bimodal synergy soliton", SCP Discoveries V28, etc.) return nothing identical. The project appears to be an original (or unpublished/exploratory) numerical research program.

**Closest conceptual parallels** (none reproduce the classical 3-scalar helical-braid numerics or bimodal synergy):

- **Bilson-Thompson helon (braid) preon model** (2005 onward, with Markopoulou & Smolin): Particles as 3-ribbon braids in quantum spacetime/LQG background-independent models. Braids encode chirality/spin, twists encode charge; Borromean-like topologies appear. Topological but quantum/discrete, not classical continuous scalars with triple-product potential.

- **Faddeev–Skyrme / Hopfion models** and knotted solitons: Stabilized by Skyrme term (higher-derivative) in a 3-component vector field (or O(3) sigma model). Produce toroidal/knotted/helical solitons with Hopf charge (topological winding). Optical/magnetic skyrmions and hopfions in condensed matter show braided structures, but use vector fields + different potentials (no separate triple-product scalars).

- **Emergent/elastic gravity and scalar-tensor models**: Verlinde entropic gravity, Sakharov induced gravity, "elastic continuum" proposals, scalar-torsion theories, and recent viXra/preprint frameworks (e.g., chronon field with winding solitons, spacetime superfluids, Xuan-Liang fluid). Strain tensors for gravity and torsion appear, but usually single scalar/vector or fluid variables — not this exact 3-field braid. No triple-product or bimodal helical numerics.

- **Other soliton literature**: Soliton foams (two scalars), domain-wall braids in noncommutative φ⁴, fuzzy DM cores, optical knotted fields, etc. — all differ in fields, potential, or dimensionality. Triple-product potentials are rare; n-ary algebras or Filippov structures appear in math but not this physics context.

The SCP log references its own negative results (mass term destroys braids, extra couplings unwind, gauge kills stability) and draws inspiration from Skyrme, oscillons, and topological protection — consistent with mainstream soliton theory. V28's novelty lies in the specific classical realization and numerical discovery of "bimodal synergy," but it remains a toy model without direct experimental or theoretical literature precedent.

Overall, the document shows rigorous, iterative numerics with transparent diagnostics and replication details. V28 advances the braid idea nicely but inherits the core limitations of classical scalar proxies for gauge/gravity phenomena. It would be interesting to quantize, add back-reaction, or compare to hopfion stability proofs for further insight.