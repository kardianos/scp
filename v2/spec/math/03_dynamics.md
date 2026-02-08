# Math Spec 03 — Dynamics and Equation of Motion

## 1. Spacetime Structure

The field $\Psi(x, t)$ is defined over 3+1D spacetime. We separate the treatment of space and time:

- **Space**: Modeled by the PGA geometric derivative $\nabla = \sum_{i=1}^3 e_i \partial_i$ from Cl(3,0,1).
- **Time**: An external evolution parameter $t$, entering via the ordinary derivative $\partial_t$.

This separation is consistent with **Process Ontology**: time is not a geometric dimension alongside space, but the parameter along which the spatial configuration evolves. It is also standard practice in topological soliton models (Skyrme, Faddeev-Niemi), where the target space algebra handles spatial structure and time enters separately.

**Lorentz Covariance**: Despite the separate treatment, the Lagrangian is constructed to be Lorentz invariant. The spacetime derivatives $\partial_\mu = (\frac{1}{c}\partial_t, \partial_1, \partial_2, \partial_3)$ contract with the Minkowski metric $\eta^{\mu\nu} = \text{diag}(+1, -1, -1, -1)$. The PGA algebra Cl(3,0,1) handles the *internal* (target space) structure of the field, not the spacetime signature.

---

## 2. The Lagrangian Density

The dynamics are governed by the action $S = \int \mathcal{L} \, d^4x$ with Lagrangian density:

$$ \mathcal{L} = \mathcal{L}_2 + \mathcal{L}_4 - V(\Psi) - V_D(\Psi) $$

### Term 1 — Quadratic Kinetic Term ($\mathcal{L}_2$)

$$ \mathcal{L}_2 = \frac{1}{2c^2} \langle \partial_t \Psi \, \widetilde{\partial_t \Psi} \rangle_0 - \frac{1}{2} \langle \nabla \Psi \, \widetilde{\nabla \Psi} \rangle_0 $$

*   The **temporal** part $\langle \partial_t \Psi \, \widetilde{\partial_t \Psi} \rangle_0$ is the time-kinetic energy density.
*   The **spatial** part $\langle \nabla \Psi \, \widetilde{\nabla \Psi} \rangle_0$ is the gradient/elastic energy density.
*   The relative **minus sign** between temporal and spatial terms gives **hyperbolic** (wave-like) character. This ensures propagating solutions exist.
*   In covariant notation: $\mathcal{L}_2 = \frac{1}{2} \eta^{\mu\nu} \langle \partial_\mu \Psi \, \widetilde{\partial_\nu \Psi} \rangle_0$

### Term 2 — Skyrme Term ($\mathcal{L}_4$) — Topological Stabilizer

$$ \mathcal{L}_4 = \frac{1}{4e^2} \sum_{\mu < \nu} \langle [R_\mu, R_\nu]^2 \rangle_0 $$

where $R_\mu = \widetilde{\Psi} \partial_\mu \Psi$ is the **left-invariant current** and $[A, B] = AB - BA$ is the commutator in $Cl^+(3,0,1)$.

*   This is a **4th-order derivative** term (each $R_\mu$ contains one derivative; the commutator squared has four factors).
*   It **vanishes** for configurations where all gradients commute (abelian/scalar perturbations).
*   It is **nontrivial** for configurations with internal twist — exactly the topologically non-trivial solitons (knots).
*   $e$ is the **Skyrme coupling constant** (dimensionless in natural units), controlling the soliton size.
*   **Precedent**: This is the standard stabilization mechanism used in the Skyrme model (baryons) and the Faddeev-Niemi model (Hopf solitons). Without this term, Derrick's theorem forbids stable static solitons in 3D (see §4).

### Term 3 — Potential ($V$)

$$ V(\Psi) = \frac{\lambda}{4} \left( \langle \Psi \widetilde{\Psi} \rangle_0 - \rho_0^2 \right)^2 $$

*   $\langle \Psi \widetilde{\Psi} \rangle_0$ is the **scalar part** of the norm (ensuring $V$ is a real number; see §5 for why this matters).
*   Enforces $|\Psi| \to \rho_0$ at spatial infinity (vacuum state).
*   The vacuum manifold is $\langle\Psi\widetilde{\Psi}\rangle_0 = \rho_0^2$, which is $S^3$ in the even subalgebra.
*   $\lambda$: self-coupling constant. $\rho_0$: vacuum density.

### Term 4 — Degenerate Mass Term ($V_D$) — Weight Sector Constraint

$$ V_D(\Psi) = \frac{\mu^2}{2}\left(|\vec{J}|^2 + \tau^2\right) $$

where $\vec{J}$ and $\tau$ are the degenerate bivector and pseudoscalar coefficients of $\Psi$ (the "weight" part in the dual quaternion decomposition $\Psi = q + e_0 p$; see §5).

*   This gives **mass $\mu$** to the pseudoscalar ($P$) and flux ($\vec{J}$) modes, which would otherwise be massless Goldstone-like modes (the standard potential $V$ cannot see them because $e_0^2 = 0$).
*   It breaks the **pseudoscalar shift symmetry** $\tau \to \tau + \alpha$.
*   It uses the **dual quaternion norm** $|p|^2 = |\vec{J}|^2 + \tau^2$, the natural invariant of the weight sector.
*   **Precedent**: Analogous to the pion mass term in the Skyrme model.

### Full Lagrangian

$$ \boxed{ \mathcal{L} = \mathcal{L}_2 + \mathcal{L}_4 - V - V_D } $$

Expanded:
$$ \mathcal{L} = \frac{1}{2c^2} \langle \partial_t \Psi \, \widetilde{\partial_t \Psi} \rangle_0 - \frac{1}{2} \langle \nabla \Psi \, \widetilde{\nabla \Psi} \rangle_0 + \frac{1}{4e^2} \sum_{\mu < \nu} \langle [R_\mu, R_\nu]^2 \rangle_0 - \frac{\lambda}{4} \left( \langle \Psi \widetilde{\Psi} \rangle_0 - \rho_0^2 \right)^2 - \frac{\mu^2}{2}(|\vec{J}|^2 + \tau^2) $$

### Parameters

| Parameter | Symbol | Role |
|-----------|--------|------|
| Vacuum density | $\rho_0$ | Sets the vacuum state |
| Bulk self-coupling | $\lambda$ | Controls potential stiffness / scalar mass |
| Skyrme coupling | $e$ | Controls soliton size |
| Degenerate mass | $\mu$ | Controls pseudoscalar and flux masses |
| Process speed | $c$ | Sets temporal/spatial ratio (speed of light) |

The theory has **five** free parameters: ($\rho_0$, $\lambda$, $e$, $\mu$, $c$). All particle masses, charges, and coupling constants must emerge from these five inputs plus the topology of the solutions.

---

## 3. The Equation of Motion

Applying the Euler-Lagrange equations to the full Lagrangian yields:

$$ \boxed{ \frac{1}{c^2}\partial_t^2 \Psi - \nabla^2 \Psi + \frac{\delta \mathcal{L}_4}{\delta \Psi} + \lambda \Psi \left( \langle \Psi \widetilde{\Psi} \rangle_0 - \rho_0^2 \right) = 0 } $$

Where $\frac{\delta \mathcal{L}_4}{\delta \Psi}$ is the Euler-Lagrange contribution from the Skyrme term (involves third-order spatial derivatives; explicit form to be derived for simulation).

Key properties:
*   **Hyperbolic**: The $\frac{1}{c^2}\partial_t^2 - \nabla^2$ structure is the **d'Alembertian** $\square$. This supports propagating wave solutions.
*   **Nonlinear**: The potential term and the Skyrme term provide the nonlinearity needed for solitons.
*   **Lorentz covariant**: The equation is invariant under Lorentz boosts (the Lagrangian is a Lorentz scalar by construction).

### Linear Limit (Weak Field)

Far from any knot, $|\Psi| \approx \rho_0$. Let $\Psi = \rho_0 + \psi$ where $|\psi| \ll \rho_0$.

The Skyrme term is $O(\psi^4)$ and negligible. The potential linearizes to a mass term. The perturbation $\psi = S + \vec{J}e_0 + \mathbf{F} + IP$ decomposes into four sectors:

| Component | Equation | Mass | Physical Identity |
|-----------|----------|------|-------------------|
| $\mathbf{F}$ (spatial bivector) | $\square \mathbf{F} = 0$ | $0$ | **Electromagnetic waves** ([04_electromagnetism.md](04_electromagnetism.md)) |
| $S$ (scalar) | $\square S + m_S^2 S = 0$ | $m_S = \sqrt{2\lambda}\,\rho_0$ | Massive scalar (Higgs-like?) |
| $P$ (pseudoscalar) | $\mu^2 P = 0$ | — | **Non-dynamical** (see §5.6) |
| $\vec{J}e_0$ (degenerate bivector) | $\mu^2 \vec{J} = 0$ | — | **Non-dynamical** (see §5.6) |

*   $\mathbf{F}$ is massless — this IS electromagnetism.
*   $S$ is massive from the bulk potential $V$ — standard Mexican-hat result.
*   **$P$ and $\vec{J}$ are non-dynamical**: The kinetic term $\langle\partial_\mu\Psi\,\widetilde{\partial_\mu\Psi}\rangle_0 = |\partial_\mu q|^2$ does not include the degenerate sector at all (because $e_0^2 = 0$ kills all weight terms in the scalar extraction). Their only contribution to the Lagrangian is $-V_D = -(\mu^2/2)(|\vec{J}|^2 + P^2)$, which is algebraic — no derivatives. The Euler-Lagrange equations give $\mu^2 J_i = 0$, $\mu^2 P = 0$: they are forced to zero. They are **not** propagating Klein-Gordon fields. See §5.6 for the full analysis and implications.

### Nonlinear Regime (The Knot)

Near the core of a knot, $|\Psi|$ deviates significantly from $\rho_0$. The Skyrme term becomes important and provides the crucial **topological stabilization** against collapse.

---

## 4. Derrick's Theorem — Now Evaded

For a static configuration $\Psi(x)$ (no time dependence), the total energy is:

$$ E[\Psi] = E_2 + E_4 + E_V + E_D $$

where:
*   $E_2 = \frac{1}{2}\int \langle\nabla\Psi\,\widetilde{\nabla\Psi}\rangle_0 \, d^3x$ — quadratic gradient energy
*   $E_4 = \frac{1}{4e^2}\int \sum_{i<j}\langle[R_i, R_j]^2\rangle_0 \, d^3x$ — Skyrme energy (spatial indices only for static case)
*   $E_V = \frac{\lambda}{4}\int(\langle\Psi\widetilde{\Psi}\rangle_0 - \rho_0^2)^2 \, d^3x$ — bulk potential energy
*   $E_D = \frac{\mu^2}{2}\int(|\vec{J}|^2 + \tau^2) \, d^3x$ — degenerate potential energy

Under **Derrick rescaling** $\Psi_\alpha(x) = \Psi(x/\alpha)$ (stretching the soliton by factor $\alpha$):

| Term | Scaling | Behavior |
|------|---------|----------|
| $E_2$ | $\to \alpha \, E_2$ | Grows as soliton expands |
| $E_4$ | $\to \alpha^{-1} E_4$ | Grows as soliton **shrinks** |
| $E_V + E_D$ | $\to \alpha^3 (E_V + E_D)$ | Grows as soliton expands |

(Both $E_V$ and $E_D$ scale as $\alpha^3$ — they are both potential-type terms with no derivatives.)

Define the **total potential energy** $E_\text{pot} = E_V + E_D$.

### Equilibrium ($dE/d\alpha = 0$ at $\alpha = 1$):

$$ E_2 - E_4 + 3 E_\text{pot} = 0 $$

### Stability ($d^2E/d\alpha^2 > 0$ at $\alpha = 1$):

$$ \frac{d^2E}{d\alpha^2}\bigg|_{\alpha=1} = 2E_4 + 6E_\text{pot} > 0 $$

This is always satisfied for non-trivial solutions ($E_4 > 0$ or $E_\text{pot} > 0$).

### Why the Skyrme Term Is Essential

**Without** $\mathcal{L}_4$ ($E_4 = 0$): The equilibrium condition becomes $E_2 + 3E_\text{pot} = 0$. Since both $E_2 \geq 0$ and $E_\text{pot} \geq 0$, the only solution is $E_2 = E_\text{pot} = 0$ — the vacuum. **No soliton exists.** This is Derrick's theorem.

**With** $\mathcal{L}_4$ ($E_4 > 0$): The equilibrium condition $E_2 = E_4 - 3E_\text{pot}$ can be satisfied non-trivially. The Skyrme energy resists collapse ($E_4$ grows as $\alpha^{-1}$), while the potential and gradient energies resist expansion. The soliton is trapped at a definite size.

---

## 5. The Degenerate Sector — Pseudoscalar and Flux Modes

### 5.1 The Problem: The Standard Potential Misses Half the Field

The pseudoscalar $I = e_{0123}$ in Cl(3,0,1) satisfies $I^2 = 0$ (because $e_0^2 = 0$). More broadly, ALL products involving $e_0^2$ vanish. This has a direct structural consequence for the potential.

**The dual quaternion decomposition**: Any even multivector $\Psi \in Cl^+(3,0,1)$ decomposes as:
$$ \Psi = \underbrace{(\rho + \mathbf{F})}_{q \,\in\, Cl^+(3,0)} + \; e_0 \underbrace{(\vec{J} + \tau I_3)}_{p \,\in\, Cl^-(3,0)} $$
where $q$ is the **bulk** (rotational) quaternion, $p$ is the **weight** (translational) part, and $I_3 = e_{123}$ is the spatial pseudoscalar ($I_3^2 = -1$).

This is the standard PGA bulk/weight decomposition. The even subalgebra $Cl^+(3,0,1)$ is isomorphic to the dual quaternions $\mathbb{H} + \varepsilon\mathbb{H}$ with $\varepsilon = e_0$, $\varepsilon^2 = 0$.

**Why the standard norm only sees $q$**: In any Clifford product involving $\Psi\widetilde{\Psi}$, the weight part $e_0 p$ contributes terms proportional to $e_0^2 = 0$. Explicitly:

$$\langle\Psi\widetilde{\Psi}\rangle_0 = \langle q\tilde{q}\rangle_0 + \underbrace{\langle(e_0 p)\widetilde{(e_0 p)}\rangle_0}_{=\;0\;\text{(contains }e_0^2)} = \rho^2 + |\mathbf{F}|^2$$

The potential $V = (\lambda/4)(\langle\Psi\widetilde{\Psi}\rangle_0 - \rho_0^2)^2 = (\lambda/4)(\rho^2 + |\mathbf{F}|^2 - \rho_0^2)^2$ is **completely flat** in both $\tau$ (pseudoscalar) and $\vec{J}$ (flux). These components have no mass term from $V$.

### 5.2 Coupling Analysis: P Is Not Decoupled

Could the pseudoscalar and flux modes simply be ignored (decoupled from observable physics)? **No.** The Skyrme term couples them.

The right-current $R_\mu = \widetilde{\Psi}\partial_\mu\Psi$ contains cross-terms between all components of $\Psi$. Specifically, for a soliton background with nonzero $\rho$ and $\mathbf{F}$, a pseudoscalar perturbation $\delta\Psi = I\delta\tau$ produces:

$$\delta R_\mu = \widetilde{\Psi}_0 \cdot I\partial_\mu(\delta\tau) + \widetilde{(I\delta\tau)} \cdot \partial_\mu\Psi_0 $$

Since $I$ maps spatial bivectors to degenerate bivectors ($I \cdot e_{ij} = e_0 \cdot (\text{spatial vector})$), these cross-terms are generically nonzero when $\Psi_0$ has bivector structure. The commutator $[R_\mu, R_\nu]$ then mixes the pseudoscalar perturbation with the soliton's internal structure.

**Conclusion**: The Skyrme term provides a nonlinear coupling between $P$ and solitons. A massless pseudoscalar propagating at $c$ and coupled to matter is experimentally excluded (analogous to axion bounds). The degenerate sector **must** be given a mass.

### 5.3 The Degenerate Mass Term ($V_D$)

Since $e_0^2 = 0$ prevents any polynomial in $\Psi\widetilde{\Psi}$ from seeing the weight components, we use the **dual quaternion norm** directly.

The weight component $p = \vec{J} + \tau I_3$ is an element of $Cl^-(3,0)$ (odd part of the spatial algebra). Its norm in $Cl(3,0)$ is:

$$|p|^2 = \langle p\tilde{p}\rangle_0^{(3)} = |\vec{J}|^2 + \tau^2$$

where the tilde denotes $Cl(3,0)$ reversion (vectors unchanged, $I_3 \to -I_3$), and $\langle\cdot\rangle_0^{(3)}$ is the $Cl(3,0)$ scalar part. The cross-terms $\vec{J}I_3$ cancel because $I_3$ commutes with all elements of $Cl(3,0)$ (it is central in odd-dimensional Clifford algebras).

We add the **degenerate mass term** to the potential:

$$ \boxed{V_D(\Psi) = \frac{\mu^2}{2}\left(|\vec{J}|^2 + \tau^2\right) = \frac{\mu^2}{2}|p|^2} $$

Properties:
- **Gives mass $\mu$ to both $P$ and $\vec{J}$** (the entire degenerate sector).
- **Preserves the vacuum**: At $\Psi = \rho_0$ (scalar), $p = 0$ so $V_D = 0$.
- **Preserves spatial rotation invariance**: $|p|^2$ is invariant under $p \to RpR^{-1}$ for spatial rotors $R$.
- **Breaks the pseudoscalar shift symmetry** $\tau \to \tau + \alpha$, which is the symmetry whose spontaneous breaking produced the Goldstone mode.
- **Uses the dual quaternion structure** intrinsic to $Cl^+(3,0,1) \cong \mathbb{H} + \varepsilon\mathbb{H}$. The bulk norm $|q|^2$ constrains the rotational sector; the weight norm $|p|^2$ constrains the translational sector.

### 5.4 Physical Interpretation of the Massive Degenerate Modes

With $V_D$, the degenerate modes acquire mass $\mu$ and become **short-range**:

- **$P$ (massive pseudoscalar, 1 component)**: Parity-odd. A massive pseudoscalar field coupled to solitons mediates a short-range, parity-violating interaction. This is suggestive of aspects of the **weak interaction** (which also violates parity and is short-range).
- **$\vec{J}$ (massive flux, 3 components)**: A massive vector-like mode with three spatial components. Three massive vector degrees of freedom are reminiscent of the $W^+$, $W^-$, $Z^0$ bosons of the weak force.

**Caution**: These identifications are highly speculative. The weak force has specific SU(2) gauge structure, chiral coupling, and distinct $W/Z$ masses that are not reproduced here. The connection is suggestive at the level of counting degrees of freedom and parity properties, not at the level of detailed dynamics. It should be treated as a direction for investigation, not a result.

### 5.5 Why the Dual Quaternion Norm Is Natural

The degenerate mass term is not ad hoc — it follows from the internal structure of $Cl^+(3,0,1)$:

1. **Algebraic motivation**: The dual quaternion algebra has two independent norms: the "study norm" $|q|^2$ and the "dual norm" $|p|^2$. The existing potential $V$ uses only $|q|^2$. Adding $V_D$ completes the picture by constraining $|p|^2$. A general potential should depend on both independent invariants of the algebra.

2. **Geometric motivation**: In PGA, the bulk $q$ encodes **rotational** structure (orientation, spin), while the weight $p$ encodes **translational** structure (displacement, momentum). A theory that constrains only rotations but leaves translations unconstrained is geometrically incomplete.

3. **Precedent**: In the Skyrme model, the "pion mass term" $V_\pi = (f_\pi^2 m_\pi^2/4)\text{Tr}(U - 1)$ explicitly breaks the chiral symmetry that would otherwise produce massless Goldstone bosons (pions). The CHPT degenerate mass term $V_D$ plays an analogous role.

### 5.6 Critical Issue: The Degenerate Sector Is Non-Dynamical

**Problem**: The kinetic term $\mathcal{L}_2 = \frac{1}{2}\eta^{\mu\nu}\langle\partial_\mu\Psi\,\widetilde{\partial_\nu\Psi}\rangle_0$ does not include the degenerate sector. This is a direct algebraic consequence of $e_0^2 = 0$:

**Proof**: For $\Psi = q + e_0 p$ with $q \in Cl^+(3,0)$ and $p \in Cl^-(3,0)$:

$$\partial_\mu\Psi = \partial_\mu q + e_0\,\partial_\mu p$$
$$\widetilde{\partial_\mu\Psi} = \widetilde{\partial_\mu q} + \widetilde{(e_0\,\partial_\mu p)}$$

The product $(\partial_\mu\Psi)(\widetilde{\partial_\mu\Psi})$ has four terms:
1. $(\partial_\mu q)(\widetilde{\partial_\mu q})$ — pure bulk, scalar part = $|\partial_\mu q|^2$
2. $(\partial_\mu q)\widetilde{(e_0\,\partial_\mu p)}$ — contains $e_0$, scalar part = 0
3. $(e_0\,\partial_\mu p)(\widetilde{\partial_\mu q})$ — contains $e_0$, scalar part = 0
4. $(e_0\,\partial_\mu p)\widetilde{(e_0\,\partial_\mu p)}$ — contains $e_0^2 = 0$

Therefore: $\langle(\partial_\mu\Psi)(\widetilde{\partial_\mu\Psi})\rangle_0 = |\partial_\mu q|^2$

The same argument applies to the Skyrme term: the right-current $R_\mu = \widetilde{\Psi}\,\partial_\mu\Psi$ has quaternion part $\tilde{q}\,\partial_\mu q$ (all $e_0$-containing terms drop out of $\langle[R_\mu, R_\nu]^2\rangle_0$). And $E_V$ depends only on $|q|^2$.

**Consequence**: The only term in the Lagrangian that involves $P$ or $\vec{J}$ is $V_D = (\mu^2/2)(|\vec{J}|^2 + P^2)$. This is algebraic (no derivatives). The Euler-Lagrange equations give:

$$\frac{\partial V_D}{\partial J_i} = \mu^2 J_i = 0, \qquad \frac{\partial V_D}{\partial P} = \mu^2 P = 0$$

The degenerate sector is **forced to zero** — it has no dynamics, no propagation, and no coupling to solitons. This is verified numerically: the 3D gradient flow code confirms $|J|_{\max} < 10^{-10}$, $|P|_{\max} < 10^{-10}$ (see `proposal/hopfion_search/src/main.c`).

**Resolution options**: To make the degenerate sector dynamical, the Lagrangian must be extended:

1. **Degenerate kinetic term** (minimal fix):
$$\mathcal{L}_{2,D} = \frac{1}{2c^2}|\dot{p}|^2 - \frac{1}{2}|\nabla p|^2$$
where $|p|^2 = |\vec{J}|^2 + \tau^2$. This gives $P$ and $\vec{J}$ proper Klein-Gordon dynamics ($\square P + \mu^2 P = 0$), but they remain **decoupled** from solitons (free massive fields in the soliton background).

2. **Bulk-degenerate coupling** (needed for interactions):
An explicit coupling term is required, e.g.:
$$\mathcal{L}_{\text{int}} = \frac{g^2}{2}|q|^2|\nabla p|^2$$
This creates a position-dependent effective mass for the degenerate modes near the soliton core (where $|q| \neq \rho_0$), enabling soliton-dependent scattering. The coupling constant $g$ becomes a sixth parameter.

3. **Component-wise Skyrme term**: Replace $\langle[R_\mu, R_\nu]^2\rangle_0$ with a norm that sees all 8 components. This would create nonlinear coupling between bulk and degenerate sectors through the topological stabilization mechanism.

**Status**: This issue is **open**. The choice of coupling mechanism has physical implications (parity violation, interaction range, coupling strength) that should be guided by the target phenomenology (weak force properties). See `15_open_problems.md`, B3.

---

## 6. Chirality in Dynamics

The operator $\nabla$ is sensitive to orientation.
*   A solution $\Psi_L$ (Left-twisted) and $\Psi_R$ (Right-twisted) are distinct.
*   If the coupling $\lambda$ or the vacuum $\rho_0$ has any internal structure (e.g., if $\rho_0$ is a spinor condensate rather than a scalar), the equation can spontaneously break parity, giving different masses/dynamics to $L$ and $R$ solutions (weak interaction mechanism).

---

## 7. Origin of Mass

Mass is the total static energy of the field configuration:

$$ M c^2 = E_2 + E_4 + E_\text{pot} = \int \left[ \frac{1}{2}\langle\nabla\Psi\,\widetilde{\nabla\Psi}\rangle_0 + \frac{1}{4e^2}\sum_{i<j}\langle[R_i, R_j]^2\rangle_0 + V(\Psi) + V_D(\Psi) \right] d^3x $$

Using the Derrick equilibrium condition $E_2 = E_4 - 3E_\text{pot}$ (where $E_\text{pot} = E_V + E_D$):

$$ M c^2 = 2E_4 - 2E_\text{pot} = 2E_4 - 2E_V - 2E_D $$

The mass depends on:
*   The Skyrme coupling $e$ (sets the soliton size scale)
*   The bulk potential parameters $\lambda$, $\rho_0$ (set the rotational energy scale)
*   The degenerate mass $\mu$ (sets the translational energy scale)
*   The topological sector $Q$ (determines the solution geometry)

**Radiation** (null-rotors, $\mathbf{F}$ mode): $E_4 = 0$, $E_\text{pot} = 0$ → massless, propagating at $c$.
**Knots** (solitons): $E_4 > 0$, $E_\text{pot} > 0$ → massive, localized.
**Degenerate radiation** ($P$, $\vec{J}$ modes): massive at $\mu$, short-range.

---

## 8. Comparison with Established Models

| Model | Field Target | Topological Charge | Skyrme Term | Soliton Type |
|-------|-------------|-------------------|-------------|-------------|
| Skyrme (1961) | SU(2) $\cong S^3$ | $\pi_3(S^3) = \mathbb{Z}$ (Baryon #) | Yes | Skyrmion |
| Faddeev-Niemi | $S^2$ | $\pi_3(S^2) = \mathbb{Z}$ (Hopf) | Yes (modified) | Hopfion |
| **CHPT** | **$Cl^+(3,0,1)$ / $S^3$** | **$\pi_3(S^3) = \mathbb{Z}$ (Charge $Q$)** | **Yes** | **Hopfion** |

CHPT's Lagrangian is structurally analogous to the Skyrme model, with two key differences:
1.  The field is valued in $Cl^+(3,0,1)$ (dual quaternions) rather than SU(2), giving additional degrees of freedom (scalar density fluctuations, pseudoscalar helicity, degenerate bivector flux).
2.  The potential term $V(\Psi)$ is essential (it fixes the vacuum norm to $\rho_0$), whereas the Skyrme model's potential (pion mass term) is optional.

---

## 9. Status of Former Critical Issues

### M1. Time Evolution — RESOLVED

Time enters as an explicit external parameter $\partial_t$, separate from the PGA spatial derivative $\nabla$. The Lagrangian $\mathcal{L}_2$ has the relative minus sign between temporal and spatial kinetic terms, making the EOM hyperbolic ($\square\Psi + \ldots = 0$). Propagating wave solutions now exist.

The PGA algebra Cl(3,0,1) is retained for **spatial/internal** structure (projective geometry, dual quaternion kinematics, Hopf fibration). Time is the **evolution parameter**, consistent with Process Ontology.

### M2. Derrick's Theorem — RESOLVED

The Skyrme term $\mathcal{L}_4$ evades Derrick's theorem by introducing a 4th-order derivative term that penalizes collapse. The Derrick equilibrium condition now has non-trivial solutions with $E_4 > 0$, permitting stable solitons at a definite size (§4). This is the same mechanism used in all successful topological soliton models.

### M3. Pseudoscalar Goldstone Mode — RESOLVED

The pseudoscalar $P$ and flux $\vec{J}$ modes were massless because $e_0^2 = 0$ makes the standard Clifford norm $\langle\Psi\widetilde{\Psi}\rangle_0$ blind to the degenerate sector. Resolution: the degenerate mass term $V_D = (\mu^2/2)(|\vec{J}|^2 + \tau^2)$ uses the dual quaternion weight norm to constrain these modes. Both $P$ and $\vec{J}$ now have mass $\mu$. See §5.

### Remaining Issues

*   **Skyrme term EOM**: **RESOLVED**. The Euler-Lagrange contribution from $\mathcal{L}_4$ has been derived analytically. It involves right-currents $A_d = \tilde{q}\,\partial_d q$, commutators $C_{d,d'} = [A_d, A_{d'}]$, and the Skyrme G-tensor $G_d = \sum_{d' \neq d}[A_{d'}, C_{d,d'}]$. The force on bulk component $a$ is:
    $$ F_{4,a}(x) = \frac{1}{2e^2}\sum_d\left[\sigma_a\langle\epsilon_a\,\partial_d q,\, G_d\rangle_0 - D_d\langle\tilde{q}\,\epsilon_a,\, G_d\rangle_0\right] $$
    where $\epsilon_a$ are the quaternion basis elements and $\sigma_a$ are reversion signs. This has been verified against finite differences for all 8 field components to $\sim 10^{-8}$ relative error. Implementation: `proposal/hopfion_search/src/field.c`.
*   **Sourced Maxwell equations**: The free EM wave equation $\square\mathbf{F} = 0$ follows from the linear limit. The sourced equation $\nabla\mathbf{F} = J$ must be derived from the nonlinear knot solution. See [../15_open_problems.md](../15_open_problems.md), B4.
*   **Weak force connection**: The massive degenerate sector (parity-odd $P$ + 3-component $\vec{J}$) has suggestive parallels with the weak interaction. This must be investigated but is currently speculative.
