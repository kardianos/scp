# Math Spec 02 — Topology and Invariants

## 1. The Field Manifold

Our field $\Psi(x)$ takes values in the even subalgebra of $Cl(3,0,1)$.
To support topological stability, the "vacuum manifold" (the space of lowest-energy values) must have non-trivial homotopy.

### The Vacuum Condition
Assume the potential $V(\Psi)$ imposes a constraint on the magnitude of the field at infinity (and in the vacuum).
$$ \Psi \widetilde{\Psi} = \rho_0^2 $$
This defines a manifold of allowed vacuum states. For the spinor representation of Cl(3,0) (and by extension P(Cl(3,0,1))), this manifold is isomorphic to $S^3$ (the 3-sphere), since unit spinors in 3D are quaternions ($\cong S^3$).

## 2. Topological Charge (The Hopf Invariant)

Since the field $\Psi$ at any point maps to $S^3$, and physical space (compactified) is $S^3$, we have a map:
$$ \Psi : S^3 \to S^3 $$
The homotopy group $\pi_3(S^3) = \mathbb{Z}$.
This integer invariant is the **Winding Number** or **Hopf Index**.

### Definition
We identify this topological index $Q$ with **Electric Charge**.
$$ Q = \int_{\mathbb{R}^3} J_{top} \, d^3x $$

Where $J_{top}$ is the topological current density. In terms of the bivector field $\mathbf{F}$ (from `01_algebra.md`):
$$ Q \propto \int \mathbf{A} \wedge \mathbf{F} $$
(This is the Chern-Simons form).

### Physical Implications
1.  **Quantization**: $Q$ is strictly integer-valued for smooth fields.
2.  **Conservation**: $Q$ is conserved because a continuous change in the field cannot jump from one integer homotopy class to another.
3.  **Particles**:
    *   $Q = +1$: Proton/Positron? (Fundamental knot)
    *   $Q = -1$: Electron/Antiproton? (Mirror knot)
    *   $Q = 0$: Neutrino / Photon.

## 3. Chirality (The Handedness)

Chirality $\chi$ is independent of $Q$. It refers to the **orientation** of the mapping.
Given a knot with $Q \neq 0$, it can still have an internal "twist" direction.
*   We define $\chi$ via the Pseudoscalar projection of the field gradient.
*   $\chi = \text{sgn}(\langle \Psi^\dagger \nabla \Psi \rangle_4)$.

This allows us to distinguish:
*   **Electron**: $Q=-1, \chi=L$
*   **Positron**: $Q=+1, \chi=R$ (CPT conjugate)
*   **Sterile/Mirror Matter?**: $Q=-1, \chi=R$ (Does Nature use this?)

## 4. Stability Mechanism

The knot is stable because it wraps around the vacuum manifold. To dissipate, it would have to unwrap, which requires passing through a singularity (infinite energy magnitude) or breaking the continuity of the field.
This is **Topological Stability** (Type A).

### Knotted Vortices (Higher Order)
Higher order knots (trefoils, etc.) may correspond to unstable resonances or higher generations. The simplest stable soliton is the **Hopfion** (a twisted torus).
