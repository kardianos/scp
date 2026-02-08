# 15 — Open Problems and Critical Assessment

This chapter provides an accounting of the current limitation of CHPT.
**Last Updated: Phase 5 — Numerical Soliton Verification**

---

## Category A — Fatal Issues (Status Update)

### A1. Bell's Theorem and Quantum Nonlocality
*   **Original Problem**: Local Determinism contradicts Bell's Inequalities.
*   **Status**: **ADDRESSED (Conceptual)**.
*   **Resolution**: Adopted **Topological Non-Locality** (see `spec/12_quantum_phenomena.md`). The field is topologically connected (instantaneous correlation) despite signal locality ($c$). This matches the "Realist/Non-local" escape route (Bohmian-style).

### A2. No Field Equation
*   **Original Problem**: Theory had no equation of motion.
*   **Status**: **RESOLVED — Equation established**.
*   **Resolution**: Full Lagrangian and hyperbolic EOM established in `spec/math/03_dynamics.md`:
    $$ \frac{1}{c^2}\partial_t^2 \Psi - \nabla^2 \Psi + \frac{\delta \mathcal{L}_4}{\delta \Psi} + \lambda\Psi(\langle\Psi\tilde{\Psi}\rangle_0 - \rho_0^2) = 0 $$
    Four-term Lagrangian with kinetic ($\mathcal{L}_2$), Skyrme ($\mathcal{L}_4$), bulk potential ($V$), and degenerate mass ($V_D$) terms. Five parameters: $(\rho_0, \lambda, e, \mu, c)$.
*   **Sub-issues resolved**: A7 (time evolution), A8 (Derrick stability), B5 (pseudoscalar Goldstone) — see below.

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

### A7. Time Evolution in PGA
*   **Problem**: The geometric derivative $\nabla = \sum e_i \partial_i$ in Cl(3,0,1) has no time component. $\nabla^2$ is the Laplacian (elliptic), not the d'Alembertian (hyperbolic).
*   **Status**: **RESOLVED**.
*   **Resolution**: Time derivative $\partial_t$ introduced as an explicit external parameter, separate from the spatial PGA derivative $\nabla$. The EOM is now hyperbolic: $\frac{1}{c^2}\partial_t^2\Psi - \nabla^2\Psi + \cdots = 0$. This is consistent with Process Ontology (time is evolution, not geometry) and follows the precedent of the Skyrme and Faddeev-Niemi models. See `spec/math/03_dynamics.md`, Section 1.

### A8. Soliton Stability — Derrick's Theorem
*   **Problem**: Derrick's theorem forbids stable static solitons in 3D for fields with only quadratic kinetic + polynomial potential.
*   **Status**: **RESOLVED — Numerically Verified**.
*   **Resolution**: Skyrme term $\mathcal{L}_4 = \frac{1}{4e^2}\sum_{\mu<\nu}\langle[R_\mu, R_\nu]^2\rangle_0$ added to the Lagrangian. Under spatial rescaling $x \to \alpha x$: $E_2 \to \alpha E_2$, $E_4 \to \alpha^{-1}E_4$, $E_V \to \alpha^3 E_V$. Equilibrium: $E_2 - E_4 + 3E_V = 0$. Stability: $d^2E/d\alpha^2 = 2E_4 + 6E_V > 0$ (automatic). See `spec/math/03_dynamics.md`, Section 4, and `spec/math/05_mass_mechanism.md`, Section 2.
*   **Numerical verification**: The $B=1$ hedgehog Skyrmion has been found numerically in the sigma model limit ($\lambda \to \infty$, $|q| = \rho_0$). The shooting method solver confirms $E/E_{FB} = 1.232$, matching the standard Skyrme model literature, with the virial theorem $E_2 = E_4$ satisfied to 6 significant figures. See `proposal/hopfion_search/`.

---

## Category B — Major Gaps

### B1. Particle Mass Spectrum
*   **Status**: **PARTIALLY ADDRESSED — First numerical result obtained**.
*   **Update**: `spec/math/05_mass_mechanism.md` defines Mass as the Hamiltonian eigenvalue. The $B=1$ soliton mass has been computed in the sigma model limit:
    $$ Mc^2 = E_{\text{total}} = 1.232 \times \frac{6\sqrt{2}\,\pi^2\,\rho_0^3}{e} $$
    For $\rho_0=1$, $e=4$: $Mc^2 = 25.78$ (dimensionless units). The mass depends on $\rho_0$ and $e$ as $M \propto \rho_0^3/e$.
*   **Remaining**: Higher-charge ($B=2,3,\ldots$) soliton masses, and the mapping to actual particles (proton, electron), have not been computed. We do not yet know if $M_{proton}/M_{electron} \approx 1836$.

### B2. fractional Charge (The Quark Problem)
*   **Status**: **NEW CRITICAL GAP**.
*   **Problem**: We defined Charge $Q$ as a Winding Number (Integer). Quarks have $Q = 1/3, 2/3$.
*   **Implication**: Either Quarks are not fundamental knots, or the topology is more complex (e.g., knots on a non-simply connected manifold).

### B3. The Weak Force
*   **Status**: **OPEN — New direction from B5 resolution**.
*   **Update**: Defined as "Parity Violation" in `spec/05_chirality.md`. The bulk Lagrangian is parity-invariant, but the newly-identified massive degenerate sector (B5 resolution) provides a candidate: the pseudoscalar $P$ is parity-odd and massive, and the 3-component flux $\vec{J}$ is also massive. Together, these 4 massive modes with parity-odd coupling are structurally suggestive of weak bosons. **Investigation needed**: does the degenerate sector couple asymmetrically to left-handed vs right-handed solitons?

### B4. Derivation of Maxwell's Equations
*   **Original Problem**: Need to prove null-rotors = Maxwell.
*   **Status**: **PARTIALLY SOLVED**.
*   **Resolution**: `spec/math/04_electromagnetism.md` shows the linear limit gives $\nabla^2 \mathbf{F} = 0$ (free electromagnetic wave equation). Null condition $\mathbf{F}^2 = 0$ gives $|\vec{E}| = |\vec{B}|$, $\vec{E} \perp \vec{B}$.
*   **Remaining**: The **sourced** Maxwell equation $\nabla \mathbf{F} = J$ (how knots produce EM fields) is not derived. The source current $J$ must be extracted from the nonlinear knot solution. This is where the real charge-field coupling lives.

### B5. Pseudoscalar Goldstone Mode
*   **Problem**: $I^2 = 0$ in Cl(3,0,1) makes the bulk potential $V$ flat in the pseudoscalar and flux directions. Both $P$ and $\vec{J}$ were massless.
*   **Status**: **RESOLVED**.
*   **Resolution**: The **degenerate mass term** $V_D = (\mu^2/2)(|\vec{J}|^2 + \tau^2)$ uses the dual quaternion weight norm $|p|^2$ to give mass $\mu$ to both the pseudoscalar $P$ and flux $\vec{J}$ modes. This exploits the $Cl^+(3,0,1) \cong \mathbb{H} + \varepsilon\mathbb{H}$ structure: the bulk potential $V$ constrains $|q|^2$, while $V_D$ constrains $|p|^2$. See `spec/math/03_dynamics.md`, Section 5.
*   **New parameter**: $\mu$ (degenerate mass scale). Theory now has 5 parameters: $(\rho_0, \lambda, e, \mu, c)$.
*   **Speculative direction**: The massive degenerate sector (parity-odd $P$ + 3-component $\vec{J}$) has suggestive parallels with the weak interaction (parity violation, short-range, 4 massive modes). This is speculative and requires investigation.

---

## Category C — Simulation Requirements

With the Lagrangian now complete (A2, A7, A8, B5 all resolved), the remaining open problems (A3, A4, B1, B2, B3) require **numerical simulation**. The equation of motion is well-defined and can be discretized.

### Completed

1.  **Static decoupling theorem (verified)**: The static energy functional was implemented for the full 8-component $Cl^+(3,0,1)$ field on a 3D lattice. Numerical experiments confirm the theoretical prediction: the static problem **decouples** — the bulk quaternion sector $(s, f_1, f_2, f_3)$ evolves under $E_2 + E_4 + E_V$ (the standard Skyrme model), while the degenerate sector $(j_1, j_2, j_3, p)$ relaxes trivially to zero under $E_D$. This was verified by gradient computation and finite-difference testing (all 8 components agree to $\sim 10^{-8}$ relative error).

2.  **$B=1$ Skyrmion profile (solved)**: In the sigma model limit ($\lambda \to \infty$, enforcing $|q| = \rho_0$), the hedgehog ansatz $q = \rho_0(\cos f(r) + \sin f(r)\,\hat{r}\cdot\boldsymbol{\sigma})$ reduces the 3D problem to a 1D ODE. The Euler-Lagrange equation is:
    $$ f''(r^2 + 2c_4\sin^2 f) + 2rf' + c_4 f'^2\sin 2f - \sin 2f - c_4\frac{\sin 2f \sin^2 f}{r^2} = 0 $$
    where $c_4 = 2\rho_0^2/e^2$. This was solved via a shooting method (RK4 + bisection), yielding:
    - $E/E_{FB} = 1.232$ (matching the standard Skyrmion result from the literature)
    - $E_2/E_4 = 1.000$ (Derrick virial theorem verified to 6 significant figures)
    - $Q = 1.000$ (topological charge exactly unity)
    - Faddeev-Bogomolny bound: $E_{FB} = 6\sqrt{2}\,\pi^2\rho_0^3/e$
    - Shooting parameter: $a = -f'(0) \propto e/\sqrt{2}$ (scales with Skyrme coupling)
    - Code: `proposal/hopfion_search/src/radial.c`

3.  **Skyrme force derivation (verified)**: The analytical gradient of the full energy functional (including the Skyrme term $E_4$) was derived and implemented. It involves right-currents $A_d = \tilde{q}\,\partial_d q$, commutators $[A_d, A_{d'}]$, and the G-tensor $G_d = \sum_{d' \neq d}[A_{d'}, C_{d,d'}]$. Verified against finite differences for all 8 field components (bulk and degenerate) at multiple test points, with relative errors $\lesssim 10^{-8}$. Code: `proposal/hopfion_search/src/field.c`, `src/verify.c`.

### Remaining

4.  **Full 3D soliton relaxation**: The 3D gradient flow driver is implemented but encounters topology loss: the soliton shrinks below grid resolution during energy minimization (equilibrium size $\sim \sqrt{c_4} \ll h$ for typical parameters). Solutions: (a) use the 1D radial profile to initialize a well-resolved 3D grid, (b) implement rational map ansatz for higher-$B$ sectors, (c) use adaptive mesh refinement.
5.  **Higher-charge solitons ($B=2, 3, \ldots$)**: Compute masses and shapes. The $B=2$ Skyrmion is known to be toroidal in the standard Skyrme model; verify this carries over.
6.  **Mass spectrum**: Compare $B$-dependence of soliton masses with known particle mass ratios. The standard Skyrme model gives $E(B) \approx 1.232 \times B \times E_{FB}$ for small $B$ (near the Bogomolny bound).
7.  **Scattering**: Collide two solitons. Observe if they repel, attract, scatter elastically, or produce new solitons.
8.  **Degenerate sector dynamics**: Study how the massive $P$ and $\vec{J}$ modes interact with solitons. Do they mediate parity-violating interactions (weak force connection)?
9.  **Finite $\lambda$ effects**: Move beyond the sigma model limit. Solve with finite $\lambda$ and study how $E_V$ modifies the soliton profile and mass formula ($Mc^2 = 2E_4 - 2E_V - 2E_D$).

