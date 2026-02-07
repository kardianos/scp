# 15 — Open Problems and Critical Assessment

This chapter provides an accounting of the current limitation of CHPT.
**Last Updated: Phase 2 Completion**

---

## Category A — Fatal Issues (Status Update)

### A1. Bell's Theorem and Quantum Nonlocality
*   **Original Problem**: Local Determinism contradicts Bell's Inequalities.
*   **Status**: **ADDRESSED (Conceptual)**.
*   **Resolution**: Adopted **Topological Non-Locality** (see `spec/12_quantum_phenomena.md`). The field is topologically connected (instantaneous correlation) despite signal locality ($c$). This matches the "Realist/Non-local" escape route (Bohmian-style).

### A2. No Field Equation
*   **Original Problem**: Theory had no equation of motion.
*   **Status**: **PROPOSED — Critical issues remain**.
*   **Resolution**: Lagrangian and EOM proposed in `spec/math/03_dynamics.md`:
    $$ \nabla^2 \Psi + \lambda \Psi (|\Psi|^2 - \rho_0^2) = 0 $$
*   **Remaining Issues**:
    1.  **Time evolution (A7)**: The geometric derivative $\nabla$ in Cl(3,0,1) is spatial-only. The equation as written is elliptic (Laplace), not hyperbolic (wave). Must be fixed before dynamics can be trusted.
    2.  **Soliton stability (A8)**: Derrick's theorem forbids stable static solitons for this Lagrangian structure in 3D without a higher-order (Skyrme) term. The Lagrangian likely needs a 4th-order derivative term.
    3.  **Goldstone modes (B5)**: The pseudoscalar perturbation P may be a massless Goldstone boson, which would be experimentally observable and is tightly constrained.

### A3. Gravitational Wave Polarization
*   **Original Problem**: Scalar field predicts wrong polarization (0 modes vs 2 observed).
*   **Status**: **PARTIALLY ADDRESSED**.
*   **Update**: The field $\Psi$ is now a **Multivector** (Cl(3,0,1)). It has Bivector components ($\mathbf{F}$) which support waves.
*   **Remaining Risk**: The *mechanism* for gravity is still "Density Depletion" (Scalar). If gravity is mediated by the scalar part $S$, it is still wrong. If it involves the geometric distortion, it might work. Needs simulation.

### A4. Color Charge and SU(3)
*   **Original Problem**: No SU(3) structure for Strong Force.
*   **Status**: **OPEN / CRITICAL**.
*   **Update**: `spec/09_nuclear_interactions.md` proposes "Topological Linking" as the Strong Force. This explains binding but does *not* explain the SU(3) symmetry or the number of gluons (8). This remains the biggest gap with the Standard Model.

### A5. Gauge Symmetry
*   **Original Problem**: No U(1) / SU(2) / SU(3) gauge invariance.
*   **Status**: **OPEN**.
*   **Note**: CHPT replaces Gauge Symmetry with **Geometric Covariance**. We assume (but have not proven) that the geometric constraints of knots impose constraints equivalent to gauge limits.

### A6. Neutrino Weak Interactions
*   **Original Problem**: "Achiral = Neutral" implies Neutrinos don't interact.
*   **Status**: **SOLVED**.
*   **Resolution**: `spec/05_chirality.md` separated Charge (Topological Index $Q$) from Chirality (Handedness $\chi$).
    *   Neutrino: $Q=0$ (Neutral), $\chi \neq 0$ (Chiral).
    *   This allows it to participate in chiral (weak) interactions without having electric charge.

### A7. Time Evolution in PGA (NEW)
*   **Problem**: The geometric derivative $\nabla = \sum e_i \partial_i$ in Cl(3,0,1) has no time component. The degenerate basis vector $e_0$ ($e_0^2 = 0$) is projective, not timelike. Therefore $\nabla^2$ is the Laplacian (elliptic), not the d'Alembertian (hyperbolic). Elliptic equations have no propagating solutions — no waves, no null-rotors.
*   **Status**: **OPEN / CRITICAL**.
*   **Fix options**: (1) Switch to Cl(1,3) spacetime algebra. (2) Introduce $\partial_t$ separately from the geometric derivative. (3) Embed into a larger algebra.

### A8. Soliton Stability — Derrick's Theorem (NEW)
*   **Problem**: Derrick's theorem proves that for a field with quadratic kinetic term and polynomial potential in $D \geq 3$, no stable static finite-energy solitons exist. The proposed CHPT Lagrangian has exactly this structure. The Skyrme model and Faddeev-Niemi model evade this by adding a 4th-order derivative term (the "Skyrme term").
*   **Status**: **OPEN / CRITICAL**.
*   **Fix**: Add a Skyrme-like term $\frac{1}{e^2}\langle[\nabla\Psi, \nabla\Psi]^2\rangle_0$ to the Lagrangian. Must verify this permits stable Hopfion solutions.

---

## Category B — Major Gaps

### B1. Particle Mass Spectrum
*   **Status**: **OPEN**.
*   **Update**: `spec/math/05_mass_mechanism.md` defines Mass as the Hamiltonian eigenvalue. However, no values have been computed. We do not know if $M_{proton}/M_{electron} \approx 1836$.

### B2. fractional Charge (The Quark Problem)
*   **Status**: **NEW CRITICAL GAP**.
*   **Problem**: We defined Charge $Q$ as a Winding Number (Integer). Quarks have $Q = 1/3, 2/3$.
*   **Implication**: Either Quarks are not fundamental knots, or the topology is more complex (e.g., knots on a non-simply connected manifold).

### B3. The Weak Force
*   **Status**: **OPEN**.
*   **Update**: Defined as "Parity Violation" in `spec/05_chirality.md`, but the Lagrangian in `spec/math/03_dynamics.md` appears Parity-invariant. Needs an explicit chiral symmetry breaking term.

### B4. Derivation of Maxwell's Equations
*   **Original Problem**: Need to prove null-rotors = Maxwell.
*   **Status**: **PARTIALLY SOLVED**.
*   **Resolution**: `spec/math/04_electromagnetism.md` shows the linear limit gives $\nabla^2 \mathbf{F} = 0$ (free electromagnetic wave equation). Null condition $\mathbf{F}^2 = 0$ gives $|\vec{E}| = |\vec{B}|$, $\vec{E} \perp \vec{B}$.
*   **Remaining**: The **sourced** Maxwell equation $\nabla \mathbf{F} = J$ (how knots produce EM fields) is not derived. The source current $J$ must be extracted from the nonlinear knot solution. This is where the real charge-field coupling lives.

### B5. Scalar and Pseudoscalar Modes (NEW)
*   **Problem**: The linearized perturbation $\psi = S + \mathbf{F} + IP$ has scalar ($S$) and pseudoscalar ($P$) components in addition to the EM bivector $\mathbf{F}$. The potential gives mass to $S$ ($m_S = \sqrt{2\lambda}\rho_0$), but $P$ may be a massless Goldstone boson. A massless pseudoscalar would be experimentally observable and is tightly constrained.
*   **Status**: **OPEN**.
*   **Action**: Determine masses and couplings of S and P. If P is massless, either add a term to the potential that breaks the pseudoscalar symmetry, or show P decouples from observable physics.

---

## Category C — Simulation Requirements

To resolve the remaining Open Problems (A3, B1, B2), we cannot just write more math. We must simulate:
1.  **Scalar Vortex (2D)**: Prove stability of the nonlinear potential.
2.  **Hopfion (3D)**: Prove existence of $Q=1$ knots.
3.  **Scattering**: Observe if they repel/attract.

