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

## 5. Sigma Model Limit and the Faddeev-Bogomolny Bound

### The Sigma Model ($\lambda \to \infty$)

In the strong-coupling limit $\lambda \to \infty$, the potential $E_V$ enforces $|q| = \rho_0$ exactly. The degenerate sector decouples and relaxes to zero ($E_D = 0$). Only $E_2$ and $E_4$ survive, and the virial theorem simplifies to:
$$ E_2 = E_4, \qquad Mc^2 = E_2 + E_4 = 2E_4 = E_{\text{total}} $$

### The Faddeev-Bogomolny Bound

For the CHPT Lagrangian in the sigma model limit, the Faddeev-Bogomolny (FB) topological energy bound is:
$$ E \geq E_{FB} = \frac{6\sqrt{2}\,\pi^2\,\rho_0^3}{e}\,|B| $$

where $B$ is the baryon number (topological charge). The $\sqrt{2}$ factor arises from the specific ratio of $E_2$ and $E_4$ prefactors in the CHPT normalization. Under the rescaling $r = R\sqrt{c_4}$ with $c_4 = 2\rho_0^2/e^2$, the hedgehog ODE maps to the standard Skyrmion ODE, giving $E = (\sqrt{2}\rho_0^3)/(2e) \times E_{\text{std}}$ where $E_{\text{std},FB} = 12\pi^2$.

### Numerical Result ($B=1$)

The $B=1$ hedgehog soliton has been computed numerically via a shooting method (RK4 integration + bisection). The hedgehog ansatz $q = \rho_0(\cos f(r) + \sin f(r)\,\hat{r}\cdot\boldsymbol{\sigma})$ reduces the 3D problem to the radial ODE:
$$ f''(r^2 + 2c_4\sin^2 f) + 2rf' + c_4 f'^2\sin 2f - \sin 2f - c_4\frac{\sin 2f\,\sin^2 f}{r^2} = 0 $$
with $f(0) = \pi$, $f(\infty) = 0$. Results for $\rho_0 = 1$, $e = 4$:

| Quantity | Value |
|----------|-------|
| $E_{\text{total}}$ | $25.782$ |
| $E_2$ | $12.891$ |
| $E_4$ | $12.891$ |
| $E_2/E_4$ | $0.99999$ |
| $Q$ | $1.000000$ |
| $E/E_{FB}$ | **$1.2314$** |
| $-f'(0)$ | $5.678$ |

The ratio $E/E_{FB} = 1.232$ matches the standard result from the Skyrme model literature (the $B=1$ Skyrmion is 23.2% above the Bogomolny bound), confirming that the CHPT bulk sector is numerically equivalent to the classical Skyrme model.

The energy scales as $M \propto \rho_0^3/e$, with the precise coefficient:
$$ Mc^2 = 1.232 \times \frac{6\sqrt{2}\,\pi^2\,\rho_0^3}{e} $$

### Higher-B Results (Rational Map Ansatz)

For $B \geq 2$, the hedgehog ansatz generalizes via the rational map approximation (Houghton, Manton, Sutcliffe 1998): $\hat{q} = \cos f(r) + \sin f(r)\,\hat{n}(\theta,\phi)\cdot\boldsymbol{\sigma}$ where $\hat{n}$ is determined by a degree-$B$ rational map $R: S^2 \to S^2$. The energy factorizes into:
$$ E_2 = 2\pi\rho_0^2 \int_0^\infty \left(f'^2 r^2 + 2B\sin^2 f\right) dr, \qquad E_4 = \frac{4\pi\rho_0^4}{e^2} \int_0^\infty \left(2Bf'^2\sin^2 f + I\frac{\sin^4 f}{r^2}\right) dr $$
where $I = (4\pi)^{-1}\int_{S^2} b^2\,d\Omega$ is the quartic angular integral of the map (with $b$ the local baryon density). The generalized ODE is:
$$ f'' = \frac{\sin 2f\left(B + Ic_4\sin^2 f/r^2\right) - 2rf' - Bc_4 f'^2\sin 2f}{r^2 + 2Bc_4\sin^2 f} $$
with $f(0)=\pi$, $f(\infty)=0$. Near $r=0$, $f(r) \sim \pi - ar^\alpha$ where $\alpha = (-1+\sqrt{1+8B})/2$.

Results for $\rho_0 = 1$, $e = 4$ (sigma model limit):

| $B$ | $I$ | $\alpha$ | $E/E_{FB}$ | $E_2/E_4$ | $E(B)/(B\cdot E_1)$ | Shape |
|-----|-----|----------|------------|-----------|---------------------|-------|
| 1 | 1.000 | 1.000 | 1.231 | 1.000 | 1.000 | Spherical |
| 2 | 5.808 | 1.562 | 1.208 | 1.000 | 0.981 | Toroidal |
| 3 | 13.577 | 2.000 | 1.184 | 1.000 | 0.962 | Tetrahedral |
| 4 | 20.650 | 2.372 | 1.137 | 1.000 | 0.923 | Cubic |

The binding energy per baryon $E(B)/(B\cdot E_1) < 1$ for all $B > 1$, confirming that multi-Skyrmions are bound. The binding fraction increases with $B$: 1.9% for $B=2$, 3.8% for $B=3$, 7.7% for $B=4$. These values are consistent with the standard Skyrme model literature.

## 6. Static Decoupling Theorem

A key structural result confirmed numerically: **the static energy functional decouples**.

The 8-component field $\Psi = (s, f_1, f_2, f_3, j_1, j_2, j_3, p)$ decomposes into bulk quaternion $q = (s, f_1, f_2, f_3)$ and degenerate part $(j_1, j_2, j_3, p)$. The energy terms depend on these sectors as:

| Term | Depends on |
|------|-----------|
| $E_2 = \frac{1}{2}\int\|\nabla q\|^2 d^3x$ | Bulk only |
| $E_4 = \frac{1}{4e^2}\int\sum_{i<j}\|[A_i, A_j]\|^2 d^3x$ | Bulk only (via $A_i = \tilde{q}\,\partial_i q$) |
| $E_V = \frac{\lambda}{4}\int(\|q\|^2 - \rho_0^2)^2 d^3x$ | Bulk only |
| $E_D = \frac{\mu^2}{2}\int(\|\vec{J}\|^2 + \tau^2) d^3x$ | Degenerate only |

The scalar extraction $\langle M \rangle_0$ in $Cl^+(3,0,1)$ kills all terms containing $e_0$, so $E_2$, $E_4$, and $E_V$ depend only on the quaternion part. The degenerate sector sees only $E_D$ (a simple quadratic potential), which has its unique minimum at $\vec{J} = 0$, $\tau = 0$.

**Consequence**: For static solitons, the bulk sector evolves under $E_2 + E_4 + E_V$ — which is exactly the **standard Skyrme model** on $S^3$. The degenerate sector trivially relaxes to zero. This was verified numerically: with all 8 components active, the degenerate components remain $< 10^{-10}$ while the bulk converges to the Skyrmion profile.

## Conclusion
CHPT mass is **not** an arbitrary parameter. It is a **calculated eigenvalue** of the field equation.
*   **Input**: Five parameters $(\rho_0, \lambda, e, \mu, c)$.
*   **Output**: Discrete mass spectrum $M_1, M_2, \dots$ corresponding to stable topological soliton solutions.
*   **Mechanism**: Mass arises from the balance of gradient energy ($E_2$), Skyrme stabilization energy ($E_4$), bulk potential ($E_V$), and degenerate potential ($E_D$), with the constraint $E_2 = E_4 - 3(E_V + E_D)$.
*   **First numerical result**: The $B=1$ soliton mass is $Mc^2 = 1.232 \times 6\sqrt{2}\pi^2\rho_0^3/e$, confirming equivalence with the standard Skyrme model in the sigma model limit.
