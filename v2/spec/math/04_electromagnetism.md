# Math Spec 04 — Electromagnetism (Linear Limit)

## 1. Linearization of the Dynamics

From `03_dynamics.md`, the full equation of motion is:
$$ \frac{1}{c^2}\partial_t^2 \Psi - \nabla^2 \Psi + \frac{\delta \mathcal{L}_4}{\delta \Psi} + \lambda \Psi(\langle\Psi\tilde{\Psi}\rangle_0 - \rho_0^2) = 0 $$

In the vacuum, $\Psi \approx \rho_0$ (scalar constant). Writing $\Psi = \rho_0 + \psi$ where $\psi$ is a small perturbation, the **linearized** equation is:
$$ \frac{1}{c^2}\partial_t^2 \psi - \nabla^2 \psi + m_\psi^2\, \psi = 0 $$
where the mass term $m_\psi^2$ depends on the grade of the perturbation (see Section 2).

**Note on operators**: The spatial derivative $\nabla = \sum e_i \partial_i$ is the PGA geometric derivative (spatial only). Time enters via the explicit $\partial_t^2$ term. The combination $\frac{1}{c^2}\partial_t^2 - \nabla^2$ gives the d'Alembertian $\square$ — this hyperbolic structure is NOT built into Cl(3,0,1) itself but arises from the Lagrangian's explicit time treatment.

## 2. Decomposition of the Perturbation

The perturbation field $\psi$ is an even multivector in $Cl^+(3,0,1)$:
$$ \psi = S + \vec{J}e_0 + \mathbf{F} + IP $$
where:
*   $S$: Scalar fluctuation (grade 0)
*   $\vec{J}e_0$: Degenerate bivector fluctuation (grade 2, degenerate)
*   $\mathbf{F}$: Spatial bivector field (grade 2, spatial)
*   $P$: Pseudoscalar fluctuation (grade 4)

### Mass Spectrum (from `03_dynamics.md`, Section 3)

| Mode | Equation | Mass | Identification |
|:---|:---|:---|:---|
| $\mathbf{F}$ (spatial bivector) | $\square \mathbf{F} = 0$ | $m_F = 0$ | **Electromagnetic field** |
| $S$ (scalar) | $\square S + 2\lambda\rho_0^2 S = 0$ | $m_S = \sqrt{2\lambda}\rho_0$ | Massive scalar (Higgs-like?) |
| $P$ (pseudoscalar) | $\square P + \mu^2 P = 0$ | $m_P = \mu$ | **Massive pseudoscalar** |
| $\vec{J}e_0$ (degenerate bivector) | $\square \vec{J} + \mu^2 \vec{J} = 0$ | $m_J = \mu$ | **Massive flux** (3 components) |

Only the spatial bivector $\mathbf{F}$ is massless — this is electromagnetism. The scalar $S$ is massive from the bulk potential $V$. The pseudoscalar $P$ and flux $\vec{J}$ are massive from the degenerate mass term $V_D = (\mu^2/2)(|\vec{J}|^2 + \tau^2)$, which uses the dual quaternion weight norm (see `03_dynamics.md`, Section 5). Without $V_D$, $P$ and $\vec{J}$ would be massless Goldstone-like modes because $e_0^2 = 0$ makes the standard Clifford norm blind to the degenerate sector.

## 3. Identification of the Electromagnetic Field

We identify the **spatial bivector component** $\mathbf{F}$ with the **Faraday tensor** of electromagnetism.
In Geometric Algebra, the electromagnetic field is a bivector:
$$ \mathbf{F} = \vec{E} + I \vec{B} $$
*   $\vec{E}$: Electric field (vector part of bivector in spacetime split).
*   $\vec{B}$: Magnetic field (bivector part spatial, or pseudovector).

### The Wave Equation for F
From the mass table above, the linearized equation for $\mathbf{F}$ is:
$$ \square \mathbf{F} = \frac{1}{c^2}\partial_t^2 \mathbf{F} - \nabla^2 \mathbf{F} = 0 $$
This decouples into:
$$ \nabla^2 \vec{E} - \frac{1}{c^2} \frac{\partial^2 \vec{E}}{\partial t^2} = 0 $$
$$ \nabla^2 \vec{B} - \frac{1}{c^2} \frac{\partial^2 \vec{B}}{\partial t^2} = 0 $$
**Result**: The bivector perturbations propagate as standard electromagnetic waves at speed $c$.

## 4. Maxwell's Equations vs Wave Equation

The wave equation $\square \mathbf{F} = 0$ is a consequence of Maxwell's free-field equations:
$$ \nabla \mathbf{F} = 0 $$
(Vacuum Maxwell equation in GA form).

Does CHPT imply $\nabla \mathbf{F} = 0$?
*   The second-order equation $\square \mathbf{F} = 0$ is the primary dynamic from linearization.
*   Maxwell's first-order equation $\nabla \mathbf{F} = 0$ is likely a **constraint** or **integrability condition** derived from the source definition.
*   If we define the source current $J = \nabla \mathbf{F}$, then in vacuum ($J=0$) we recover Maxwell.

**Open Problem (B4)**: The **sourced** Maxwell equation $\nabla \mathbf{F} = J$ — how knots produce EM fields — is not derived. The source current $J$ must be extracted from the nonlinear knot solution.

## 5. Polarization and Chirality

A plane wave solution for $\mathbf{F}$ has the form:
$$ \mathbf{F}(x,t) = F_0 e^{i(k \cdot x - \omega t)} $$
The geometric algebra naturally handles polarization.
*   The "null rotors" (§06) are solutions where $\mathbf{F}^2 = 0$ (light-like).
*   This corresponds to $|\vec{E}| = |\vec{B}|$ and $\vec{E} \perp \vec{B}$.
*   This matches the properties of physical photons.

## Conclusion
The linear limit of the CHPT equation of motion naturally contains the **free electromagnetic wave equation** for the spatial bivector mode $\mathbf{F}$. Identifying $\mathbf{F}$ with the Faraday tensor recovers the correct propagation speed ($c$), polarization structure, and wave mechanics. The theory also predicts: a massive scalar $S$ (mass $m_S = \sqrt{2\lambda}\rho_0$), and a massive degenerate sector consisting of the pseudoscalar $P$ and flux $\vec{J}$ (both with mass $\mu$, short-range).
