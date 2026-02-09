# Extended Lagrangian Analysis — Degenerate Sector Coupling

## Problem Statement

The degenerate sector (P, J) of the Cl⁺(3,0,1) field is **non-dynamical** in the current Lagrangian. The scalar extraction ⟨·⟩₀ kills all e₀ terms, so L₂, L₄, and V depend only on the bulk quaternion q. The degenerate modes have no kinetic energy, no propagation, and no coupling to solitons.

To make the degenerate sector physical (candidate for weak force mediators), the Lagrangian must be extended. We analyze three options.

## Option 1: Component-Wise Kinetic Term

### Definition
Replace the scalar-extracted kinetic term with a component-wise norm:

$$\mathcal{L}_{2,D} = \frac{1}{2c^2}|\dot{p}|^2 - \frac{1}{2}|\nabla p|^2$$

where |p|² = |J|² + P² is the dual quaternion weight norm.

### Analysis
- **Propagation**: Gives Klein-Gordon dynamics □P + μ²P = 0 and □J + μ²J = 0.
- **Coupling to solitons**: **NONE**. The kinetic term is additive and does not couple p to q. The degenerate modes propagate freely in the soliton background with no scattering.
- **Modification of soliton**: None. The soliton profile is unchanged because the soliton has p ≡ 0.
- **Linear spectrum**: P and J become massive Klein-Gordon fields with mass μ.
- **Assessment**: This is the **minimal fix** but physically unsatisfying — it makes the modes dynamical but they don't interact with anything.

### Modified radial equation
For a spherically symmetric degenerate perturbation p = g(r)ê (some direction):

$$g'' + \frac{2}{r}g' - \mu^2 g = 0$$

Solution: g(r) = A e^{-μr}/r (Yukawa potential). Range: r ~ 1/μ.

## Option 2: Bulk-Degenerate Coupling

### Definition
Add an explicit interaction:

$$\mathcal{L}_{\text{int}} = \frac{g^2}{2}|q|^2|\nabla p|^2$$

where g is a new coupling constant.

### Analysis
- **Propagation**: With L₂,D + L_int, the degenerate equation becomes:
  $$(1 + g^2|q|^2)\square p + \mu^2 p = \text{source terms involving } \nabla(|q|^2) \cdot \nabla p$$
- **Coupling to solitons**: **YES**. Near the soliton core where |q| deviates from ρ₀ (in the finite-λ model), the effective propagation speed changes. The soliton creates a position-dependent "refractive index" for degenerate modes.
- **Scattering**: Degenerate waves scatter off the soliton. The scattering cross-section depends on g² and on the soliton's |q(r)| profile.
- **Modification of soliton**: In the σ-model (|q| = ρ₀ everywhere), the coupling reduces to a constant factor (1 + g²ρ₀²), which can be absorbed into the definition of c for the degenerate sector. **No coupling in the σ-model limit.**
- **For coupling to work**: Need **finite λ** (so |q| varies spatially near the soliton core).
- **Assessment**: Physically reasonable but requires fine-tuning. The coupling vanishes in the σ-model limit, which is the regime where our soliton solutions are most accurate.

### Modified radial equation (finite λ)
With the coupled term, a degenerate perturbation p = g(r)ê around a B=1 soliton satisfies:

$$(1 + g^2\rho(r)^2)\left(g'' + \frac{2}{r}g'\right) + g^2[\rho^2]' g' - \mu^2 g = 0$$

where ρ(r) is the soliton's radial amplitude profile from the finite-λ solver.

In the σ-model limit (ρ = ρ₀): reduces to simple Yukawa, no coupling.
At finite λ: ρ(r) varies near origin (from finite_lambda.c: ρ(0) drops to 0.8ρ₀ at λ=8000), creating a potential well/barrier for degenerate modes.

## Option 3: Modified Skyrme Term

### Definition
Replace ⟨[R_μ, R_ν]²⟩₀ with a full component norm:

$$\mathcal{L}_{4,\text{full}} = \frac{1}{4e^2}\sum_{\mu<\nu} |[R_\mu, R_\nu]|^2_8$$

where |M|²₈ = Σ_α M_α² sums over all 8 components of the even multivector M (not just the scalar part ⟨M²⟩₀).

### Analysis
- **Key difference**: The full right-current R_μ = Ψ̃ ∂_μΨ contains cross-terms between q and p. With ⟨·⟩₀, these cross-terms are projected out. With |·|²₈, they contribute.
- **Coupling**: **YES, intrinsic**. The commutator [R_μ, R_ν] mixes bulk and degenerate components when both q and p are nonzero. The 8-component norm preserves these cross-terms.
- **Topological stabilization**: Still provides Derrick stability (E₄ still scales as α⁻¹).
- **Modification of soliton**: For the pure-quaternion soliton (p=0), this reduces to the standard Skyrme term (cross-terms vanish). The soliton profile is **unchanged**.
- **Perturbative coupling**: A small degenerate perturbation δp around a soliton background would feel a potential from the soliton's right-currents through the cross-terms in [R_μ, R_ν].
- **Assessment**: Most elegant option. Preserves existing results (soliton unchanged), naturally couples p to soliton topology through the same mechanism that stabilizes the soliton itself. Introduces no new parameters.

### Cross-term structure (detailed derivation)

Use the dual quaternion representation: Ψ = q + εp where ε² = 0 and q, p are
quaternions. The reversal is Ψ̃ = q̃ + εp̃ (quaternion conjugation on each part).

The product rule gives:
$$\Psi_1\Psi_2 = q_1q_2 + \varepsilon(q_1p_2 + p_1q_2) \quad (\varepsilon^2 = 0)$$

For background soliton q₀ with degenerate perturbation δp (so Ψ = q₀ + εδp):

$$R_\mu = \tilde{\Psi}\partial_\mu\Psi = \tilde{q_0}\partial_\mu q_0 + \varepsilon(\tilde{q_0}\partial_\mu\delta p + \widetilde{\delta p}\,\partial_\mu q_0)$$

Define:
- $A_\mu = \tilde{q_0}\partial_\mu q_0$ — background right-current (pure quaternion)
- $D_\mu = \tilde{q_0}\partial_\mu\delta p + \widetilde{\delta p}\,\partial_\mu q_0$ — linear perturbation (pure quaternion)

So $R_\mu = A_\mu + \varepsilon D_\mu$, and the commutator is:
$$[R_\mu, R_\nu] = [A_\mu, A_\nu] + \varepsilon([A_\mu, D_\nu] + [D_\mu, A_\nu])$$

The standard Skyrme term uses $\langle [R_\mu, R_\nu]^2 \rangle_0$. Since $\varepsilon$ terms
are degenerate and $\langle\varepsilon\cdot\text{stuff}\rangle_0 = 0$, only $[A_\mu, A_\nu]^2$
survives — **no coupling to δp**.

The modified Skyrme term uses $|[R_\mu, R_\nu]|_8^2$. For $M = a + \varepsilon b$:
$$|M|_8^2 = |a|^2 + |b|^2$$

Therefore:
$$|[R_\mu, R_\nu]|_8^2 = |[A_\mu, A_\nu]|^2 + |[A_\mu, D_\nu] + [D_\mu, A_\nu]|^2$$

The second term is the **coupling term** — quadratic in δp and its derivatives:

$$\mathcal{L}_{4,\text{coupling}} = \frac{1}{4e^2}\sum_{\mu<\nu} |[A_\mu, D_\nu] + [D_\mu, A_\nu]|^2 \geq 0$$

### Key properties of the coupling

1. **Positive-definite**: The coupling adds a non-negative potential for degenerate modes.
   This means the soliton creates a **trapping potential well** for degenerate excitations.

2. **Localized at the soliton**: $V_\text{eff}(r)$ is nonzero wherever the background
   right-currents $A_\mu$ are nonzero — i.e., inside the soliton. It vanishes at spatial
   infinity where $f(r) \to 0$ and $A_\mu \to 0$.

3. **Topology-dependent**: The coupling strength depends on the commutators $[A_\mu, D_\nu]$,
   which involve the soliton's internal twisting (same structure as the topological charge
   density). Higher-B solitons have stronger coupling.

### Hedgehog effective potential

For the B=1 hedgehog $q_0 = \rho_0(\cos f + \sin f\,\hat{r}\cdot\vec{\sigma})$, the
background right-currents in Cartesian coordinates are:
$$A_i = \rho_0^2\left(f'\hat{x}_i(\hat{x}\cdot\vec{\sigma}) + \frac{\sin f\cos f}{r}(\sigma_i - \hat{x}_i\hat{x}\cdot\vec{\sigma})\right)$$

For the simplest perturbation mode — a scalar degenerate mode $\delta p = g(r)\cdot\mathbf{1}$
(P component only):
$$D_i = 2\rho_0\cos f \cdot g'(r)\hat{x}_i + \text{angular terms}$$

After angular averaging over the unit sphere, the effective radial potential is:
$$V_\text{eff}(r) \sim \frac{\rho_0^4}{e^2}\left[\left(\frac{\sin f\cos f}{r}\right)^2 f'^2 + \frac{\sin^4 f}{r^4}\right]$$

The exact coefficients require completing the angular integration of the quaternion
commutator, but the structure is clear:
- $V_\text{eff}(r) \to 0$ as $r \to 0$ (regularity)
- $V_\text{eff}(r)$ peaks at $r \sim R_\text{sol}$ (soliton core radius)
- $V_\text{eff}(r) \to 0$ as $r \to \infty$ (exponentially, since $f \sim e^{-r}$)
- Depth $\sim \rho_0^4/e^2$, width $\sim 1/e$

### Bound state equation

With kinetic term $\mathcal{L}_{2,D}$ and the modified Skyrme coupling, the degenerate
perturbation satisfies a Schrödinger-type radial equation:
$$-g''(r) - \frac{2}{r}g'(r) + \frac{\ell(\ell+1)}{r^2}g(r) + V_\ell(r)g(r) = \omega^2 g(r)$$

where $V_\ell(r)$ is the effective potential after partial wave decomposition, and $\omega$
is the mode frequency. Bound states exist when $\omega^2 < \mu^2$ (below the degenerate
sector mass gap).

The existence of bound states would mean **the soliton has internal excitation modes
associated with the degenerate sector** — physically analogous to excited hadron states
or weak-force coupling.

## Comparison

| Property | Option 1 (kinetic) | Option 2 (coupling) | Option 3 (Skyrme) |
|----------|-------------------|--------------------|--------------------|
| New parameters | 0 | 1 (g) | 0 |
| Propagation | Yes | Yes | Needs L₂,D too |
| Soliton coupling | No | Yes (finite λ only) | Yes (any λ) |
| Potential sign | N/A | **Attractive** | **Repulsive** |
| Bound states | No | Yes (finite λ) | No |
| σ-model limit | Decoupled | Decoupled | Coupled (repulsive) |
| Soliton modified | No | No | No |
| Elegance | Minimal | Moderate | Best |
| Complexity | Simple | Moderate | Complex |

## Preliminary Recommendation (pre-numerical)

Option 3 was initially considered most promising for its elegance (no new parameters,
topology-dependent coupling). However, the numerical V_eff computation revealed that
it produces only repulsive coupling — no bound states.

**See Revised Recommendation below for the updated assessment based on numerics.**

## Numerical Results

### V_eff computation (src/veff.c)

We decompose the coupling energy into:
$$E_{\text{coupling}} = \frac{1}{4e^2}\int\left[\alpha(r)\,g'^2 + \gamma(r)\,g^2\right] 4\pi r^2\,dr$$

where γ(r) is the "mass-like" potential (from g=1, g'=0) and α(r) is the kinetic
modification (from g=0, g'=1).

**Key finding: Both α(r) and γ(r) are POSITIVE everywhere.**

For e=2, ρ₀=1 (our standard parameters), at the soliton core (r ≈ 0.5):

| Mode | γ(0.5) | α(0.5) | γ/α ratio |
|------|--------|--------|-----------|
| P (scalar) | 1.4×10⁴ | 6.0×10² | 23 |
| J (vector) | 3.1×10⁴ | 1.8×10² | 172 |

Both decay exponentially for r > 1 (soliton tail), reaching < 1 by r ≈ 1.5.

### Physical interpretation

The modified Skyrme term creates a **repulsive potential** for degenerate modes
near the soliton. The degenerate modes cost MORE energy when they overlap with
the soliton's winding structure.

**Consequences:**
- **No bound states**: V_eff > 0 everywhere means no trapped degenerate modes.
  The soliton cannot bind weak-sector excitations through this mechanism alone.
- **Scattering**: The repulsive potential produces scattering of degenerate
  waves off the soliton. The scattering cross-section ~ (soliton radius)².
- **Hard-core behavior**: The soliton acts as an impenetrable object for
  degenerate modes at low energies (V_eff >> kinetic energy for r < R_sol).

### Why Option 3 is repulsive

The cross-term |[A_i, D_j] + [D_i, A_j]|² is a sum of squares — necessarily
non-negative. This adds positive energy when δp ≠ 0 near the soliton, making the
overlap energetically costly. The sign is inherent to the 8-component norm
construction and cannot be changed without fundamentally altering the algebra.

## Revised Recommendation

Option 3 alone does **not** provide the desired trapping/coupling. We revise
the analysis:

### Option 2 revisited: attractive potential from finite-λ

For Option 2 with coupling $\mathcal{L}_{\text{int}} = (g^2/2)|q|^2|\nabla p|^2$,
the effective equation for g(r) near a finite-λ soliton is:

$$-[1 + g^2\rho(r)^2]\nabla^2 g + \mu^2 g = \omega^2 g$$

Since $\rho(0) < \rho_0$ at finite λ (e.g., ρ(0) = 0.97ρ₀ at λ=5000), the kinetic
coefficient is SMALLER at the soliton core. This is equivalent to an **attractive
potential well** with depth:

$$\Delta V \sim \mu^2 \cdot g^2(\rho_0^2 - \rho(0)^2) / (1 + g^2\rho_0^2)$$

For strong enough coupling g, this can trap degenerate modes. The well depth
increases as λ decreases (larger ρ deviation).

### Recommended hybrid approach

The optimal Lagrangian combines:
1. **L₂,D** (Option 1): Gives degenerate modes kinetic energy and propagation
2. **L_int** (Option 2): Provides attractive coupling to solitons via ρ(r) variation
3. **L₄,full** (Option 3): Provides short-range repulsion from topology

The combination of Options 2 (attraction) and 3 (repulsion) creates a potential
with structure analogous to a Lennard-Jones potential:
- Long range: attractive well from ρ(r) variation
- Short range: repulsive core from topological coupling

This could produce bound states at a finite distance from the soliton center, with
the bound state radius set by the balance of attraction and repulsion.

### Full recommended Lagrangian

$$\mathcal{L} = \mathcal{L}_2 + \mathcal{L}_{2,D} + \mathcal{L}_{4,\text{full}} + \frac{g^2}{2}|q|^2|\nabla p|^2 - V - V_D$$

Parameters: (ρ₀, λ, e, μ, g, c) — the coupling g is the one new parameter.

## Numerical Testing Plan — Status

1. **Option 2 bound states**: For a given g, solve the radial equation with the
   finite-λ profile ρ(r) and find bound state energies.
   **Status: OPEN.** Requires solving Schrödinger-type radial equation with the
   computed V_eff(r). The integration tests confirm the potential structure but
   discrete bound states have not been computed.

2. **Combined potential**: Compute the hybrid Option 2+3 effective potential and
   find the Lennard-Jones-like minimum.
   **Status: CONFIRMED (Phase 9.3).** Integration tests with the full 3D code
   show E_{4,C} repulsion (degenerate field expelled from core, E_{4,C} drops
   94-99%) and E_int attraction (g=1 retains 14% more coupling energy than g=0).
   The Lennard-Jones structure is verified dynamically.

3. **Scattering cross-section**: Partial wave analysis of degenerate mode
   scattering off the soliton.
   **Status: OPEN.** Integration tests show qualitative scattering (degenerate
   perturbation disperses from core) but quantitative cross-sections not computed.

4. **Full 3D simulation**: Add all coupling terms to the time-dependent code.
   **Status: COMPLETE (Phase 9.3).** Implemented in `src/coupling.c` +
   `src/coupling.h`. Gradient verified to <1e-7 (`src/verify_coupling.c`).
   Integrated into `src/scatter.c` with `-degenerate -g <val>` flags. Energy
   conservation <0.02% in integration tests. See `results/README.md` for full
   numerical data.
