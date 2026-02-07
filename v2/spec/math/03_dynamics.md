# Math Spec 03 — Dynamics and Equation of Motion

## 1. The Action Principle

The dynamics of the field $\Psi(x)$ are governed by the principle of stationary action:
$$ S = \int \mathcal{L}(\Psi, \nabla \Psi) \, d^4x $$

### The Lagrangian Density
We propose a Lagrangian density consisting of a **Kinetic Term** (geometric covariance) and a **Potential Term** (vacuum stability).

$$ \mathcal{L} = \underbrace{\frac{1}{2} \langle \nabla \Psi \widetilde{\nabla \Psi} \rangle_0}_{\text{Kinetic}} - \underbrace{\frac{\lambda}{4} (|\Psi|^2 - \rho_0^2)^2}_{\text{Potential}} $$

*   $\nabla = \sum e_i \partial_i$: Geometric gradient.
*   $\widetilde{\nabla \Psi}$: Reverse of the gradient field (needed for scalar product).
*   $\langle \dots \rangle_0$: Scalar projection.
*   $|\Psi|^2 = \Psi \widetilde{\Psi}$: Squared magnitude.
*   $\rho_0$: Vacuum density parameter.
*   $\lambda$: Self-interaction coupling constant.

## 2. The Equation of Motion

Applying the Euler-Lagrange equations to this Lagrangian yields the **Nonlinear Vector Wave Equation**:

$$ \nabla^2 \Psi + \lambda \Psi (|\Psi|^2 - \rho_0^2) = 0 $$
*(Note: Time evolution is implicit in the d'Alembert operator $\square$ or explicitly separated depending on metric signature. In Cl(3,0,1) with degenerate time, this is formally $\nabla^2$ but interpreted as evolution along the process parameter).*

### Linear Limit (Weak Field)
Far from a knot, $|\Psi| \approx \rho_0$. Let $\Psi = \rho_0 + \psi$ where $\psi$ is small.
The potential term vanishes (to first order), leaving:
$$ \nabla^2 \psi \approx 0 $$
This is the **Wave Equation** for free radiation (null-rotors). Propagating at speed $c$.

### Nonlinear Regime (The Knot)
Near the core, $|\Psi|$ deviates significantly from $\rho_0$. The nonlinear term $\lambda \Psi (|\Psi|^2 - \rho_0^2)$ provides the "restoring force" or "mass term" that allows localized solutions (solitons) to exist.

## 3. Stability and Chirality

### Topological Stability
The potential $V(\Psi) = (|\Psi|^2 - \rho_0^2)^2$ enforces the boundary condition $|\Psi| \to \rho_0$ at infinity. This compactifies the domain to $S^3$, allowing non-trivial Hopf invariants ($Q \neq 0$) as defined in [02_topology.md](02_topology.md).

### Chirality in Dynamics
The operator $\nabla$ is sensitive to orientation.
*   A solution $\Psi_L$ (Left-twisted) and $\Psi_R$ (Right-twisted) are distinct.
*   If the coupling $\lambda$ or the vacuum $\rho_0$ has any internal structure (e.g., if $\rho_0$ is a spinor condensate rather than a scalar), the equation can spontaneously break parity, giving different masses/dynamics to $L$ and $R$ solutions (Weak interaction mechanism).

## 4. Origin of Mass

Mass is the energy of the field configuration.
$$ E = \int T_{00} \, d^3x = \int \left( \frac{1}{2} |\nabla \Psi|^2 + V(\Psi) \right) d^3x $$
*   **Radiation**: Has Kinetic energy but zero Potential energy (it stays near either vacuum or linear regime).
*   **Knots**: Have both Kinetic (internal twist) and Potential (displacement from vacuum) energy. This total integrated energy behaves as inertial mass.
