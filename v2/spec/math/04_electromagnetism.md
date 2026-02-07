# Math Spec 04 — Electromagnetism (Linear Limit)

## 1. Linearization of the Dynamics

From `03_dynamics.md`, the field equation in the vacuum (where $|\Psi| \approx \rho_0$) is the **Vector Wave Equation**:
$$ \nabla^2 \psi = 0 $$
where $\psi = \Psi - \rho_0$ is the perturbation field.

In $Cl(3,0,1)$ (or spacetime algebra $Cl(1,3)$), the operator $\nabla^2$ is the d'Alembertian $\square = \partial_t^2 - \nabla^2$.
Thus, every component of the multivector field $\psi$ propagates at speed $c$.

## 2. Decomposition of the Perturbation

The field $\psi$ is an even multivector. We decompose it into components:
$$ \psi = S + \mathbf{F} + I P $$
*   $S$: Scalar fluctuation (Higgs-like? Dilaton? Gravity?)
*   $\mathbf{F}$: Bivector field (6 degrees of freedom).
*   $P$: Pseudoscalar fluctuation (Axion?).

## 3. Identification of the Electromagnetic Field

We identify the **Bivector Component** $\mathbf{F}$ with the **Faraday Tensor** of electromagnetism.
In Geometric Algebra, the electromagnetic field is a bivector:
$$ \mathbf{F} = \vec{E} + I \vec{B} $$
*   $\vec{E}$: Electric Field (Vector part of bivector in space-time split).
*   $\vec{B}$: Magnetic Field (Bivector part spatial, or pseudo-vector).

### The Wave Equation for F
Since $\nabla^2 \psi = 0$ holds for the whole multivector, it holds linearly for the bivector part:
$$ \nabla^2 \mathbf{F} = 0 $$
$$ \implies \square (\vec{E} + I \vec{B}) = 0 $$
This decouples into:
$$ \nabla^2 \vec{E} - \frac{1}{c^2} \frac{\partial^2 \vec{E}}{\partial t^2} = 0 $$
$$ \nabla^2 \vec{B} - \frac{1}{c^2} \frac{\partial^2 \vec{B}}{\partial t^2} = 0 $$
**Result**: The bivector perturbations propagate as standard electromagnetic waves at speed $c$.

## 4. Maxwell's Equations vs Wave Equation

The Wave Equation $\nabla^2 \mathbf{F} = 0$ is a consequence of Maxwell's Free Field Equations:
$$ \nabla \mathbf{F} = 0 $$
(Vacuum Maxwell Equation in GA).

Does CHPT imply $\nabla \mathbf{F} = 0$?
*   The second-order equation $\nabla^2 \psi = 0$ is the primary dynamic.
*   Maxwell's first-order equation $\nabla \mathbf{F} = 0$ is likely a **Constraint** or **Integrability Condition** derived from the source definition.
*   If we define the Source Current $J = \nabla \mathbf{F}$, then in vacuum ($J=0$) we recover Maxwell.

## 5. Polarization and Chirality

A plane wave solution for $\mathbf{F}$ has the form:
$$ \mathbf{F}(x,t) = F_0 e^{i(k \cdot x - \omega t)} $$
The geometric algebra naturally handles polarization.
*   The "Null Rotors" (§06) are solutions where $\mathbf{F}^2 = 0$ (Light-like).
*   This corresponds to $|\vec{E}| = |\vec{B}|$ and $\vec{E} \perp \vec{B}$.
*   This matches the properties of physical photons.

## Conclusion
The linear limit of the CHPT Equation of Motion naturally contains the **Free Electromagnetic Wave Equation**.
Identifying the bivector component $\mathbf{F}$ with the EM field recovers the correct propagation dynamics ($c$), polarization structure, and wave mechanics.
