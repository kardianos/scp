# Math Spec 05 — General Mass Mechanism

## 1. Mass Definition

In CHPT, the **rest mass** $M$ of a particle is defined as the total static energy of its field configuration (the soliton solution).

From the Lagrangian density defined in `03_dynamics.md`, the static Hamiltonian has **four terms**:

$$ \mathcal{H} = \underbrace{\frac{1}{2}\langle \nabla\Psi\, \widetilde{\nabla\Psi}\rangle_0}_{E_2\text{ density}} + \underbrace{\frac{1}{4e^2}\sum_{i<j}\langle [R_i, R_j]^2 \rangle_0}_{E_4\text{ density}} + \underbrace{\frac{\lambda}{4}(\langle\Psi\tilde{\Psi}\rangle_0 - \rho_0^2)^2}_{E_V\text{ density}} + \underbrace{\frac{\mu^2}{2}(|\vec{J}|^2 + \tau^2)}_{E_D\text{ density}} $$

where $R_i = \tilde{\Psi}\partial_i\Psi$ are spatial right-currents.

The mass is the integral over all space:
$$ Mc^2 = \int_{\mathbb{R}^3} \mathcal{H}\, d^3x = E_2 + E_4 + E_V + E_D $$

*   **$E_2$ (Kinetic/gradient term)**: Encodes the "twist" or "winding" energy from spatial gradients.
*   **$E_4$ (Skyrme term)**: Encodes higher-order derivative stabilization. Essential for soliton stability (Derrick's theorem evasion).
*   **$E_V$ (Bulk potential)**: Encodes the "displacement" from the vacuum in the rotational sector ($|q|^2 = \rho^2 + |\mathbf{F}|^2$).
*   **$E_D$ (Degenerate potential)**: Constrains the translational sector ($|p|^2 = |\vec{J}|^2 + \tau^2$). Gives mass $\mu$ to the pseudoscalar and flux modes.

## 2. Stability and Scale (Derrick's Theorem)

### The Problem

Hobart–Derrick's theorem states that for a field with only a quadratic kinetic term and polynomial potential in $D \geq 3$ spatial dimensions, no stable static finite-energy solitons exist. The original CHPT Lagrangian (without the Skyrme term) had exactly this structure.

### The Fix: Four-Term Scaling

Under the rescaling $x \to \alpha x$ (with $\alpha = 1$ at equilibrium), the energy contributions scale as:

| Term | Scaling | Direction |
|:---|:---|:---|
| $E_2$ (gradient) | $\alpha E_2$ | Decreases as knot expands |
| $E_4$ (Skyrme) | $\alpha^{-1} E_4$ | Increases as knot expands |
| $E_V + E_D$ (potentials) | $\alpha^3 (E_V + E_D)$ | Increases as knot expands |

Both potential terms ($E_V$ and $E_D$) scale as $\alpha^3$ — they have no derivatives. Define $E_\text{pot} = E_V + E_D$.

The total energy under rescaling:
$$ E(\alpha) = \alpha E_2 + \alpha^{-1} E_4 + \alpha^3 E_\text{pot} $$

### Equilibrium Condition

Setting $dE/d\alpha = 0$ at $\alpha = 1$:
$$ E_2 - E_4 + 3E_\text{pot} = 0 $$

This is the **virial theorem** for CHPT solitons. It has nontrivial solutions because $E_4$ enters with a **negative** sign (opposing expansion).

### Stability Condition

The second derivative at equilibrium:
$$ \frac{d^2 E}{d\alpha^2}\bigg|_{\alpha=1} = 2E_4 + 6E_\text{pot} > 0 $$

This is automatically satisfied since both $E_4 > 0$ and $E_\text{pot} > 0$.

### Mass Formula

Combining the total energy $Mc^2 = E_2 + E_4 + E_\text{pot}$ with the virial relation $E_2 = E_4 - 3E_\text{pot}$:
$$ Mc^2 = 2E_4 - 2E_\text{pot} = 2E_4 - 2E_V - 2E_D $$

Soliton mass is determined by the balance of the Skyrme term against both potentials.

## 3. The Mass Spectrum (Hierarchy)

Why do different particles have different masses?
The integral $M$ depends on the specific geometry of the solution $\Psi(x)$.

### Topological Multipliers
The mass depends on the topological charge $Q$ (Hopf invariant) and the internal geometry:

1.  **Fundamental knot** ($Q=1$, Hopfion):
    *   Minimum energy configuration in the $Q=1$ sector.
    *   Mass = $M_0$. Identified with the **electron** (if light) or **proton** (if heavy/composite).
    *   *Note: Standard Model mass ratios (proton/electron ~ 1836) suggest the electron might be the fundamental $Q=1$ knot and the proton a composite, OR both are $Q=1$ with different internal structure.*

2.  **Higher winding**:
    *   Energy scales non-linearly with charge $Q$.
    *   Typically $E(Q) \approx |Q| \cdot E_1$ to first order (Bogomolny bound), but knot geometry adds interaction terms.

### Resonance Modes (Excited States)
A stable knot can have internal vibrational modes that add to its total energy:
$$ M_n = M_{\text{ground}} + n \cdot \hbar\omega $$
These correspond to heavier particles (muon? tau?) which eventually decay ($n \to 0$) by emitting null-rotors.

### Scaling Law
From the virial theorem and dimensional analysis with the five parameters $(\rho_0, \lambda, e, \mu, c)$:
$$ M \propto \frac{\rho_0}{e} \times f(Q, \lambda e^2) $$
where $f$ is a dimensionless function of the topological charge and the coupling ratio. The function $f$ must be computed numerically for each knot type.

## 4. The Origin of Inertia ("Drag")

Why does this energy $M$ resist acceleration?
When a knot accelerates, the field configuration $\Psi(x,t)$ must update.
$$ \vec{F} = \frac{d\vec{p}}{dt} = \frac{d}{dt} \int \vec{\mathcal{P}}_{\text{field}}\, d^3x $$
The field momentum density $\vec{\mathcal{P}}_{\text{field}}$ is proportional to $\mathcal{H} \vec{v}$.
Thus, the resistance to change in motion (inertia) is exactly the resistance of the field energy density to redistribution.
**$M_{\text{inertial}} \equiv M_{\text{field}}$**.

## Conclusion
CHPT mass is **not** an arbitrary parameter. It is a **calculated eigenvalue** of the field equation.
*   **Input**: Five parameters $(\rho_0, \lambda, e, \mu, c)$.
*   **Output**: Discrete mass spectrum $M_1, M_2, \dots$ corresponding to stable topological soliton solutions.
*   **Mechanism**: Mass arises from the balance of gradient energy ($E_2$), Skyrme stabilization energy ($E_4$), bulk potential ($E_V$), and degenerate potential ($E_D$), with the constraint $E_2 = E_4 - 3(E_V + E_D)$.
