# 15 — Open Problems and Critical Assessment

This chapter provides an accounting of the current limitation of CHPT.
**Last Updated: Phase 9.3 — Degenerate Sector Coupling Implementation**

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
*   **Update (Phase 6)**: Higher-charge solitons ($B=1$–$4$) have been computed via the rational map ansatz. Binding energies: $E(B)/(B\cdot E_1)$ = 0.981 ($B=2$), 0.962 ($B=3$), 0.923 ($B=4$). Multi-Skyrmions are bound. See `spec/math/05_mass_mechanism.md`, §5.
*   **Remaining**: The mapping to actual particles (proton, electron) has not been established. We do not yet know if $M_{proton}/M_{electron} \approx 1836$.

### B2. fractional Charge (The Quark Problem)
*   **Status**: **NEW CRITICAL GAP**.
*   **Problem**: We defined Charge $Q$ as a Winding Number (Integer). Quarks have $Q = 1/3, 2/3$.
*   **Implication**: Either Quarks are not fundamental knots, or the topology is more complex (e.g., knots on a non-simply connected manifold).

### B3. The Weak Force
*   **Status**: **OPEN — Coupling mechanism implemented, bound states not yet computed**.
*   **Update**: Defined as "Parity Violation" in `spec/05_chirality.md`. The massive degenerate sector (B5 resolution) was a candidate: 4 massive modes (parity-odd $P$ + 3-component $\vec{J}$) suggestive of weak bosons.
*   **Critical Problem (Phase 7)**: The degenerate sector is **completely decoupled** from solitons in the original Lagrangian (scalar extraction $\langle\cdot\rangle_0$ kills all $e_0$ terms). See B6.
*   **Resolution (Phase 8–9.3)**: Three coupling terms implemented and numerically verified:
    *   $\mathcal{L}_{2,D}$: degenerate kinetic energy (propagation)
    *   $E_{4,C} = (1/4e^2)\sum|F^w_{ij}|^2$: Skyrme cross-coupling (short-range repulsion)
    *   $E_{\text{int}} = (g^2/2)|q|^2|\nabla p|^2$: attractive binding at finite $\lambda$
    *   Gradient verified to $< 10^{-7}$ relative error. Integration tests confirm Lennard-Jones-like potential (repulsion + attraction). Energy conservation $\Delta E/E < 2 \times 10^{-4}$.
*   **Remaining**: Whether this potential supports discrete bound states (trapped degenerate modes = weak bosons?) requires solving the radial Schrödinger equation with the computed effective potential. See `proposal/hopfion_search/results/README.md`.

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
*   **Speculative direction**: Superseded by B6; the degenerate sector must first be made dynamical before its phenomenological role can be assessed.

### B6. Degenerate Sector Decoupling — Non-Dynamical Weight Modes
*   **Status**: **PARTIALLY RESOLVED (Phase 9.3) — Coupling implemented and verified; bound states not yet computed**.
*   **Problem**: The kinetic term $\mathcal{L}_2 = \frac{1}{2}\eta^{\mu\nu}\langle\partial_\mu\Psi\,\widetilde{\partial_\nu\Psi}\rangle_0$ evaluates to $\frac{1}{2}|\partial_\mu q|^2$ — it does not include the degenerate sector ($\vec{J}$, $P$) at all. This is because $e_0^2 = 0$ in Cl(3,0,1), so the scalar extraction $\langle\cdot\rangle_0$ is blind to the weight part $e_0 p$. The same holds for $\mathcal{L}_4$ and $V$. The only term involving $P$ and $\vec{J}$ is $V_D = (\mu^2/2)(|\vec{J}|^2 + P^2)$, which is algebraic (no derivatives). The EOM is $\mu^2 J_i = 0$, $\mu^2 P = 0$ — the degenerate sector is forced to zero. It has no dynamics, no propagation, and no coupling to solitons.
*   **Verified**: Numerically confirmed — 3D gradient flow gives $|J|_{\max} < 10^{-10}$, $|P|_{\max} < 10^{-10}$.
*   **Implication**: The linear-limit table in `spec/math/03_dynamics.md` §3 previously claimed $\square P + \mu^2 P = 0$ (Klein-Gordon). This was **incorrect** — corrected to $\mu^2 P = 0$ (algebraic).
*   **Resolution (Phase 8 analysis + Phase 9.3 implementation)**:
    Three coupling terms derived from the Clifford product structure (cross-terms discarded by scalar extraction) and implemented in `src/coupling.c` + `src/coupling.h`:
    1. $\mathcal{L}_{2,D} = (1/2)|\nabla p|^2$ — degenerate kinetic energy (propagation).
    2. $E_{4,C} = (1/4e^2)\sum_{i<j}|F^w_{ij}|^2$ — weight-sector Skyrme cross-coupling (repulsive, no new parameters).
    3. $E_{\text{int}} = (g^2/2)|q|^2|\nabla p|^2$ — attractive coupling at finite $\lambda$ (one new parameter $g$).
    **Gradient verification**: All 8 field components match finite differences to $< 10^{-7}$ relative error (`src/verify_coupling.c`).
    **Integration tests** (e=1, λ=10000, N=128, leapfrog):
    - Energy conservation $\Delta E/E < 2 \times 10^{-4}$ over 2–3 time units.
    - $E_{4,C}$ repulsion confirmed: degenerate field expelled from soliton core ($E_{4,C}$ drops 94–99%).
    - $E_{\text{int}}$ attraction ($g=1$): retains 14% more coupling energy than $g=0$ (weak trapping at $\lambda=10000$ where $\rho(0)=0.997$; deeper trapping expected at lower $\lambda$).
    - Lennard-Jones-like potential structure confirmed dynamically: short-range repulsion from topology + long-range attraction from $\rho(r)$ variation.
    **Parameters**: Theory now has 6 parameters $(\rho_0, \lambda, e, \mu, g, c)$.
*   **Remaining**: Bound state computation — solve the radial Schrödinger equation with the effective potential to determine if discrete trapped modes exist. See `proposal/hopfion_search/results/README.md` and `results/extended_lagrangian_analysis.md`.

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

4.  **Full 3D soliton verification** ✓: 3D initialization from 1D profiles works to $<0.1\%$ for all $B=1$–$4$ (`src/verify3d.c`). Gradient flow **always** loses topology at every tested resolution ($N=128$–$192$). The sigma-model Skyrmion is a **lattice saddle point** — the discrete topology cannot prevent unwinding once the core shrinks below grid resolution. Resolution: abandon gradient flow; use **time evolution** with finite $\lambda$ instead. See item 7 below.
5.  **Higher-charge solitons ($B=2, 3, 4$)** ✓: Solved via rational map ansatz (`src/rational_map.c`). Angular integrals $I$ match literature. Binding energies: 1.9% ($B=2$), 3.8% ($B=3$), 7.7% ($B=4$) per baryon. Shapes: toroidal ($B=2$), tetrahedral ($B=3$), cubic ($B=4$).
6.  **Mass spectrum**: The $B$-dependence is now known: $E/E_{FB}$ decreases from 1.231 ($B=1$) to 1.137 ($B=4$). Mapping to actual particle mass ratios remains open.
7.  **Soliton scattering** ✓ (COLLISION OBSERVED): Full 3D time-dependent code implemented (`src/scatter.c`). Uses leapfrog (Störmer-Verlet) integration with product ansatz initialization and 3D charge-weighted centroid tracking. **Key result**: Head-on collision of two $B=1$ Skyrmions in the repulsive channel (π-isorotation) at $v = 0.5c$ with **perfect topology preservation** ($Q = 1.9999$) for $2.35$ time units — through the entire approach phase and deep core interpenetration. The solitons decelerate continuously ($18\times$ slowdown) under the repulsive interaction, reaching closest approach at $r \approx 1.23$ ($0.87$ core radii). Lattice topology loss at $t \approx 2.4$ prevents observing the bounce. **Critical parameter**: grid points across soliton core — using $e=1$ (core radius $\sqrt{2}$) at $N=192$, $L=10$ gives $13.6$ grid points and $\sim 2.4t$ stability. The $\sigma$-model profile with product ansatz eliminates the breathing mode ($|q_1 \cdot q_2/\rho_0| = \rho_0$ exactly). See `results/README.md` for full numerical data.
8.  **Degenerate sector dynamics** ✓ (IMPLEMENTED AND VERIFIED): All three coupling terms implemented in `src/coupling.c` + `src/coupling.h`, gradient verified to $< 10^{-7}$ relative error (`src/verify_coupling.c`), and integrated into the time-dependent code (`src/scatter.c` with `-degenerate -g <val>` flags):
    - **Option 1** (kinetic $\mathcal{L}_{2,D}$): degenerate gradient energy, gives propagation.
    - **Option 2** ($g^2|q|^2|\nabla p|^2$): **attractive** well at finite $\lambda$. Integration tests confirm: retains 14% more coupling energy with $g=1$ vs $g=0$.
    - **Option 3** (full 8-component Skyrme norm cross-terms $|F^w_{ij}|^2$): **repulsive** near soliton. Integration tests confirm: $E_{4,C}$ drops 94–99% as degenerate field is expelled from core.
    - Combined: Lennard-Jones-like potential confirmed dynamically. Energy conservation $\Delta E/E < 2 \times 10^{-4}$.
    - **Remaining**: Bound state computation (solve radial Schrödinger equation with effective potential to find trapped modes).
9.  **Finite $\lambda$ effects** ✓: Solved via coupled f-shooting + $\rho$-BVP Newton iteration with under-relaxation (`proposal/hopfion_search/src/finite_lambda.c`). Self-consistent solutions (virial $\approx 0$, $Q = 1.000$) obtained for $\lambda$ from $10^8$ down to $8000$. Key results:
    - $\rho(0)$ decreases from $\rho_0$ as $\lambda$ decreases, with $\delta\rho(0) \sim -\text{source}/(2\lambda)$ for large $\lambda$.
    - Finite-$\lambda$ effects **lower** the soliton mass: $E/E_{FB}$ decreases from 1.232 ($\lambda=\infty$) to 0.891 ($\lambda=8000$).
    - The Faddeev-Bogomolny bound is violated ($E/E_{FB} < 1$) at $\lambda \approx 9000$, because the FB bound assumes $\rho \equiv \rho_0$.
    - The virial theorem $E_2 - E_4 + 3E_V = 0$ holds for all converged solutions. Mass formula: $Mc^2 = 2E_4 - 2E_V$.
    - Below $\lambda \approx 7000$–$8000$, the solver collapses — the soliton core hollows out ($\rho \to 0$). The $\sigma$-model soliton ($\rho \equiv \rho_0$) is a **local** energy minimum; the **global** minimum has $\rho \to 0$ (topological collapse). The constraint $\rho \equiv \rho_0$ provides topological protection.
    - The $\rho$ BVP is exponentially stiff (growth rate $\sim e^{\sqrt{2\lambda}\,R}$), requiring BVP methods rather than shooting.
    - See `spec/math/05_mass_mechanism.md`, §7 for the full numerical table.

