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
*   **Status**: **SOLVED**.
*   **Resolution**: Derived in `spec/math/03_dynamics.md`:
    $$ \nabla^2 \Psi + \lambda \Psi (|\Psi|^2 - \rho_0^2) = 0 $$

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
*   **Status**: **SOLVED**.
*   **Resolution**: `spec/math/04_electromagnetism.md` shows the linear limit of the vector wave equation is $\nabla^2 \mathbf{F} = 0$, which is the free Maxwell equation.

---

## Category C — Simulation Requirements

To resolve the remaining Open Problems (A3, B1, B2), we cannot just write more math. We must simulate:
1.  **Scalar Vortex (2D)**: Prove stability of the nonlinear potential.
2.  **Hopfion (3D)**: Prove existence of $Q=1$ knots.
3.  **Scattering**: Observe if they repel/attract.

