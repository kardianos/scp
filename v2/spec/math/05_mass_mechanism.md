# Math Spec 05 — General Mass Mechanism

## 1. Mass Definition

In CHPT, the **Rest Mass** $M$ of a particle is defined as the total static energy of its field configuration (the soliton solution).
From the Lagrangian density $\mathcal{L}$ defined in `03_dynamics.md`:

$$ \mathcal{H} = \frac{1}{2} |\nabla \Psi|^2 + V(\Psi) $$

The mass is the integral of the Hamiltonian density over all space:
$$ M c^2 = \int_{\mathbb{R}^3} \left[ \frac{1}{2} |\nabla \Psi|^2 + \frac{\lambda}{4} (|\Psi|^2 - \rho_0^2)^2 \right] d^3x $$

*   **Kinetic Term** ($\nabla \Psi$): Encodes the "Twist" or "Winding" energy.
*   **Potential Term** ($V$): Encodes the "Displacement" from vacuum.

## 2. Stability and Scale (Derrick's Theorem)

Hobart and Derrick's Theorem states that for a scalar field in 3D, static solitons are unstable against rescaling without a potential term.
If we scale $x \to \alpha x$:
*   Kinetic energy $E_K \to \alpha E_K$ (Decreases as knot expands).
*   Potential energy $E_V \to \alpha^3 E_V$ (Increases as knot expands).
*   Stable equilibrium exists at $\frac{dE}{d\alpha} = 0$, implying $E_K = 3 E_V$.

This **Virial Theorem** for CHPT knots gives us a robust scaling law:
$$ M \propto \sqrt{\lambda} \rho_0 \times (\text{Topological Factor}) $$

## 3. The Mass Spectrum (Hierarchy)

Why do different particles have different masses?
The integral $M$ depends on the specific geometry of the solution $\Psi(x)$.

### Topological Multipliers
The "Topological Factor" is roughly proportional to the complexity of the knot.
1.  **Fundamental Knot** ($Q=1$, Hopfion):
    *   Minimum energy configuration.
    *   Mass = $M_0$. identified with the **Electron** (if light) or **Proton** (if heavy).
    *   *Note: Standard Model mass ratios (proton/electron ~ 1836) suggest the Electron might be the "Fundamental" and Proton a composite, OR Electron is Q=1 and Proton is a complex Q=1 state.*

2.  **Higher Winding**:
    *   Energy scales non-linearly with charge $Q$.
    *   Typically $E(Q) \approx |Q| \cdot E_1$ to first order, but knot geometry adds interaction terms.

### Resonance Modes (Excited States)
Just as a guitar string has harmonics, a stable knot can have internal vibrational modes that add to its total energy.
$$ M_{n} = M_{ground} + n \cdot \hbar \omega $$
These correspond to unstable, heavier particles (Muon? Tau?) which eventually decay ($n \to 0$) by emitting null-rotors.

## 4. The Origin of Inertia ("Drag")

Why does this energy $M$ resist acceleration?
When a knot accelerates, the field configuration $\Psi(x,t)$ must update.
$$ \vec{F} = \frac{d\vec{p}}{dt} = \frac{d}{dt} \int \vec{\mathcal{P}}_{field} d^3x $$
The field momentum density $\vec{\mathcal{P}}_{field}$ is proportional to $\mathcal{H} \vec{v}$.
Thus, the resistance to change in motion (Inertia) is exactly the resistance of the field energy density to redistribution.
**$M_{inertial} \equiv M_{field}$**.

## Conclusion
CHPT Mass is **not** an arbitrary parameter. It is a **calculated eigenvalue** of the field equation.
*   **Input**: Coupling $\lambda$, Vacuum density $\rho_0$.
*   **Output**: Discrete Mass Spectrum $M_1, M_2, \dots$ corresponding to stable topological solutions.
