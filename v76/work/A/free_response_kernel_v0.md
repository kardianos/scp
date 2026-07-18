# Free-Response Kernel \(K\) as a Monist Single-Sector Object — v0

**Approach:** A (forward theoretical)  
**Date:** 2026-07-18  
**Round:** 2  
**Theme (O-001):** Make the free-response monist or kill it.  
**Depends on:** `axioms_v0.md` (v0.1), C path-cost / no-go, B kernel residue, D degeneracy  
**Status:** definition + kill criteria + candidate dynamics — not a closed PDE theory

---

## 0. The critical gap (why this note exists)

Round 1 preferred shape:

\[
\text{compact lock }E_\star
\;+\;
\text{nonlocal free-response}
\;\Rightarrow\;
\ell(r)-\ell_0 \propto \frac{E_\star}{c^2 r}.
\]

B’s sandbox implements

\[
\Phi(x)=\alpha\int\frac{\rho_{\mathrm{bound}}(x')}{R+\varepsilon}\,dA',
\quad n=1+\Phi,
\]

which **hits** the monopole target but **postulates** a Poisson-like integral.  
Algebraically that is the same map dualist gravity uses. D shows ray isomorphism: fit alone cannot certify monism.

**Question for Round 2:**  
Is there a **single-sector free continuum law** whose stationary linear response *is* that Green kernel (or an observationally equivalent path-cost tail), without a second “gravity” field?

If yes → monist \(K\). If every workable \(K\) requires an independent metric/potential sector → kill the kernel branch as dualism in monist clothing.

---

## 1. Vocabulary

| Symbol | Meaning |
|--------|---------|
| \(\rho_b,\rho_f\) | Bound / free **budget** densities (ledger of one continuum) |
| \(E_\star\) | Unlock = rest ledger of a lock \(=\int\mathcal{E}[\rho_b]\) (Ax6) |
| \(c\) | Free-field locality (Ax3) |
| \(\psi\) | **Free-medium state field(s)** beyond bare budget (capacity potential, pressure, free-phase stiffness, graph voltages, …) |
| \(K\) | Linear free-response operator / Green kernel of free dynamics about uniform free vacuum |
| \(\ell\) | Free-signal path-cost density in a global chart (Ax8) |
| \(\alpha,G_{\mathrm{eff}}\) | Emergent amplitudes from free constitutive constants — **not** primitive Newton constants of a second sector |

**Split to avoid smuggling:**

- **Ledger sector (bookkeeping):** \(\rho_f,\rho_b,E_\star\) — what is locked vs free.  
- **Free dynamics sector (same continuum):** \(\psi\) and its update law — *how free medium rearranges and carries signals*.  
- **Not allowed:** an independent potential \(\Phi_{\mathrm{grav}}\) with no free-medium interpretation whose only job is to bend rays.

\(K\) lives in the free dynamics sector. Path cost \(\ell\) is a functional of free state \((\rho_f,\psi)\), not of a foreign field.

---

## 2. What \(K\) **is** (monist definition)

### 2.1 Definition

**Free-response kernel \(K\)** = the linear map from an infinitesimal lock / bound source to the first-order change in free-medium observables that control free-signal path cost, **when that map is induced by the free continuum’s own constitutive dynamics**.

Formally, about a uniform free vacuum \(\psi=\psi_0\), \(\rho_b=0\), \(\rho_f=\rho_0\):

\[
\delta\psi(x)
\;=\;
\int K_\psi(x,x')\, J_{\mathrm{lock}}(x')\,dV'
\;+\;
O(J^2),
\]

where \(J_{\mathrm{lock}}\) is a **source built only from bound formation / lock occupancy** (e.g. \(J\propto \partial_t\rho_b\), or \(\rho_b\) as a frozen constraint, or free-capacity sink rate). Then

\[
\delta\ell(x)
\;=\;
\mathcal{L}'[\psi_0;\rho_0]\,\delta\psi(x)
\;+\;
\mathcal{L}_\rho'\,\delta\rho_f(x)
\;+\;\cdots
\]

with \(\mathcal{L}'\) the Fréchet derivative of path-cost functional w.r.t. free state.

**\(K\) is monist only if:**

1. \(\psi\) is a degree of freedom of the **same** continuum as \(\rho_f/\rho_b\) (or is \(\rho_f\) itself under nonlocal dynamics);  
2. the evolution equation for \(\psi\) is **one** medium law (no second pass);  
3. \(J_{\mathrm{lock}}\) is not “stress-energy of matter on a stage” but **budget conversion free↔bound** inside that law;  
4. \(c\) remains the free-update bound (Ax3), and local frames keep free speed \(c\) (Ax7).

### 2.2 What \(K\) is **not**

| Impostor | Why dualist |
|----------|-------------|
| \(\Phi\) solved from \(\nabla^2\Phi=4\pi G\rho_b\) with \(G\) foreign and \(\Phi\) unused by free budget dynamics | S2 / S9: second sector |
| Hand-set \(\Phi=\alpha\int\rho_b/R\) with no free evolution that produces it | **Postulated dualist residue** (B-003 honesty; O-001 gap) even if labeled “monist” |
| Path cost from flat \(\int T_{00}\) of a soliton on fixed grid | S6 / S8 |
| Independent metric \(g_{\mu\nu}\) evolved by Einstein with \(T[\phi]\) | S1+S2 complete dualism |
| Soft-coded Einstein residuals labeled free response | D softE cheat |

**Algebraic look-alike is not automatic dualism:** the free Laplacian Green function *is* \(\sim 1/r\) in 3D. Monism requires that Green function to be **of free medium**, not of an extra field.

### 2.3 Slogan

> **\(K\) is how free fabric rearranges when some fabric locks — not how matter tells spacetime to curve.**

---

## 3. Kill criteria (when \(K\) is dualist smuggling)

Apply these to any candidate (theory or B sandbox). Fail **any** hard kill ⇒ dualist / ineligible for PROBLEM §7.

### 3.1 Hard kills (instant S2)

| ID | Kill if… |
|----|----------|
| **HK1** | Ray warp requires a solver that does **not** update free-medium state \(\psi\) or \(\rho_f\) — only writes a potential for rays. |
| **HK2** | Free budget / free dynamics can be frozen and rays still bend via an independent \(\Phi\). |
| **HK3** | Coupling \(G\) or \(\alpha\) is a free parameter of a gravity sector, **not** fixed by free constitutive constants \((D,\kappa,c,\rho_0,\ldots)\). (Calibration of emergent \(\alpha\) from free law is OK; floating Newton \(G\) is not.) |
| **HK4** | Two-pass schedule: (i) place locks, (ii) solve gravity, (iii) integrate rays — with (ii) not a free-medium equation. |
| **HK5** | Ledger link broken: \(E_\star\) not tied to free capacity removed; \(M_{\mathrm{ray}}\) free-fits independently of \(E_\star\) by design. |

### 3.2 Soft kills (needs redesign)

| ID | Soft kill if… |
|----|----------------|
| **SK1** | Stationary \(\delta\ell\) matches dualist Plummer/Poisson **and** no measurable free-state field \(\psi\) exists that D/B can export as 1-sector data. |
| **SK2** | Only ray isomorphism (D-003) is offered as monist proof — Occam/ledger ignored. |
| **SK3** | Local \(n(\rho_f)\) only, compact \(\rho_b\), claims long-range Einstein-class (already Round-1 dead). |
| **SK4** | Kernel works only in chart units that violate Ax7 (local free speed not \(c\)). |

### 3.3 Dualist-isomorphism rule (from D)

If monist \(\delta\ell=\beta M/R_\varepsilon\) and dualist \(\Phi=-GM/R_a\) are ray-isomorphic at matched params:

- **L_fit alone never certifies monism.**  
- Monism requires **exportable free-state evolution** + **sector count 1** + **free–bound ledger link**.  
- Congruence package (§ congruence note) encodes this for B/D.

### 3.4 Strip-list map (C S2)

| Strip ID | How it applies to \(K\) |
|----------|-------------------------|
| S2 | \(K\) must not be \(T\to G\); may be free Green of one medium |
| S9 | Bolt-on Poisson after soliton energy = hard fail |
| S12 | One evolution of free+bound; not matter then gravity |
| S1 | No independent metric; \(g\) only from free null structure |
| S6 | Mass not flat \(\int T_{00}\) alone |

---

## 4. Candidate dynamical principles (single sector)

Each candidate: **law → linear response → path cost → exterior monopole?**  
Verdict: monist-eligible / dualist risk / kill.

### Candidate F1 — Free capacity potential with lock sink (constrained relaxation)

**Law (schematic):** Free medium carries a capacity potential \(\psi\) (units: free reconfiguration activity). In the quasistatic limit,

\[
-\nabla\cdot\bigl(\sigma(\rho_f)\,\nabla\psi\bigr)
\;=\;
s\,\partial_t\rho_b
\quad\text{or, frozen lock:}\quad
-\nabla\cdot(\sigma\nabla\psi)= s\,\rho_b^{\mathrm{eff}},
\]

with conductivity \(\sigma>0\) a free constitutive function, \(s\) a conversion factor. Free path cost

\[
\ell = \ell_0 + \gamma\,\psi
\quad\text{(or }\ell=\ell_0+\gamma|\nabla\psi|^2\text{, model choice)}.
\]

**Linear vacuum:** \(\sigma=\sigma_0\), \(-\sigma_0\nabla^2\psi = s\rho_b\) ⇒

\[
\psi = \frac{s}{\sigma_0}\int G(x,x')\rho_b(x')\,dV',
\quad
G_{3D}\sim\frac{1}{4\pi R},\quad
G_{2D}\sim-\frac{1}{2\pi}\log R.
\]

Exterior monopole in 3D: \(\psi\sim (s E_\star/\sigma_0)/(4\pi r)\) if \(E_\star\sim\int\rho_b\).

**Monist reading:** \(\psi\) is free continuum state; RHS is budget conversion, not foreign \(T_{00}\).  
**Dualist risk:** if code **sets** \(\psi=\alpha\int\rho_b/R\) without evolving \(\psi\) from \(\sigma,\rho_f\), F1 collapses to postulated kernel (B residue).  
**Verdict:** **Monist-eligible** iff \(\psi\) is evolved/solved as free law with \(\sigma,\rho_f\) free DOF; **dualist** if \(\psi\) is only a ray bookkeeping field.

**FOR_B:** Implement relaxation/diffusion for \(\psi\) with source \(\propto\rho_b\) (or \(\partial_t\rho_b\)); set \(\ell=\ell_0+\gamma\psi\); tag `sector=1`, `phi_origin=free_relaxation`.

### Candidate F2 — Diffusive free-budget deficit with nonlocal optics from strain

**Law:**

\[
\partial_t \rho_f = D\nabla^2\rho_f - \Gamma_{\mathrm{lock}}[\rho_b,\rho_f],
\quad
\rho_f+\rho_b=\rho_{\mathrm{tot}}\ \text{(or integral form)}.
\]

Steady exterior may be \(\rho_f\to\rho_0\) (compact deficit only) — then **local** \(n(\rho_f)\) dies by C no-go.

**Escape inside F2:** path cost from **strain / flux**, not \(\rho_f\):

\[
\ell-\ell_0 = \gamma\,|\mathbf{j}|,\quad \mathbf{j}=-D\nabla\rho_f,
\]

or cumulative potential of past lock formation (history kernel). Steady \(\nabla^2\rho_f=0\) outside with multipole sources can give slow exterior gradients; need careful 3D integrability.

**Verdict:** **Conditional.** Pure local \(n(\rho_f)\) dead. Strain/flux path cost **monist-eligible** if \(\mathbf{j}\) is free DOF; long-range \(1/r\) for \(\ell\) not automatic — B must measure.

### Candidate F3 — Graph / causal free Laplacian Green (B1-aligned)

**Law:** Free medium on graph; free-update weights \(w_{ij}>0\). Locks remove free capacity on vertices or freeze edges. Stationary free “voltage” \(V\) for unit free-current injection at locks:

\[
L_{\mathrm{free}} V = J_{\mathrm{lock}},
\quad
K = L_{\mathrm{free}}^{+}.
\]

Path cost between nodes = resistance distance / free commute time (locality \(c\) = max free hop rate).

**Exterior continuum limit:** graph Green → continuum Green of free Laplacian → \(1/r\) (3D) / log (2D).

**Monist reading:** Distance **is** free-path cost; no chart truth (S10 strip).  
**Verdict:** **Strongest monist eligibility** philosophically. Numerics heavier (B Round 2+).

### Candidate F4 — Wave / Helmholtz free envelope

**Law:** Free field \(u\) with

\[
\partial_{tt}u = c^2\nabla^2 u - m_{\mathrm{eff}}^2[\rho_b]\,u,
\]

or index from free stiffness. Stationary / time-averaged intensity defines \(\ell\). Yukawa-like kernels if \(m_{\mathrm{eff}}>0\); massless limit → long-range.

**Verdict:** **Monist-eligible** if \(u\) is free continuum; risk of reintroducing fixed background wave operator as stage (flag S1 carefully: operator coefficients must come from free state).

### Candidate F5 — Hand Poisson \(\Phi=\alpha\int\rho_b/R\) (B Round-1 kernel)

**Law:** none for free medium; definition only.

**Verdict:** **Dualist residue / provisional phenomenology.** Allowed as **target shape** for monist candidates to match, **not** as monist proof. Tag `sector=1_claimed` but `phi_origin=postulated` → D must not treat as theory win (O-001).

### Summary table

| ID | Principle | Exterior \(\sim 1/r\)? | Monist if… | Kill if… |
|----|-----------|----------------------|------------|----------|
| **F1** | Free potential relaxation / conductivity | Yes (3D Green) | \(\psi\) free DOF; \(\sigma(\rho_f)\); one solve | \(\psi\) only for rays; foreign \(G\) |
| **F2** | Diffuse \(\rho_f\) + strain path cost | Maybe | \(\ell\) from free fluxes | local \(n(\rho_f)\) long-range claim |
| **F3** | Graph free Laplacian Green | Yes (continuum limit) | Path cost = free resistance | Fixed edge length as true distance |
| **F4** | Free wave envelope | Yukawa→Coulomb limit | \(c\), mass from free state | Fixed Minkowski stage ontology |
| **F5** | Hand \(\int\rho_b/R\) | Yes by fiat | **Never as theory** | Always provisional |

**Preferred Round-2 theory target:** **F1** (analytic clarity) + **F3** (relational purity). F5 only as D/B fit target.

### Round-3 update (O-003 / B-M2)

- **B-M2** (2D free Laplace) **endorsed** as monist F3/F1 mechanism — exterior is **log**, not a monism fail.  
- Einstein-class \(\ell\sim 1/r\) requires **3D free Green** — see `free_response_3d_v0.md`.  
- F5 hand \(1/r\) remains dead as monist theory (`monist_kernel_failed`).

---

## 5. How exterior path cost \(\propto E_\star/(c^2 r)\) can follow **without** a second gravity sector

### 5.1 Derivation sketch (F1, 3D)

1. Free vacuum: \(\psi=0\), \(\rho_b=0\), \(\rho_f=\rho_0\), \(\sigma=\sigma_0\).  
2. Compact lock: \(\int\rho_b\,dV = E_\star\) (ledger units \(\mathcal{E}=\mathrm{id}\)).  
3. Quasistatic free law: \(-\sigma_0\nabla^2\psi = s\rho_b\).  
4. Multipole: \(\psi(r)=\dfrac{s E_\star}{4\pi\sigma_0 r}+O(1/r^2)\).  
5. Path cost: \(\ell-\ell_0=\gamma\psi\) (linear response of free signalling cost to free potential).  
6. Match C target \(\ell-\ell_0 \approx 2 G_{\mathrm{eff}} M/(c^2 r)\) with \(M=E_\star/c^2\):

\[
\gamma\frac{s E_\star}{4\pi\sigma_0}
\;=\;
\frac{2 G_{\mathrm{eff}} E_\star}{c^4}
\quad\Rightarrow\quad
G_{\mathrm{eff}}
\;=\;
\frac{\gamma s\, c^4}{8\pi\sigma_0}.
\]

**All of \(\gamma,s,\sigma_0,c\) are free-medium constitutive parameters.**  
\(G_{\mathrm{eff}}\) is **emergent** — not a second-sector Newton constant.

### 5.2 What was *not* used

- No Einstein equation.  
- No independent metric dynamics.  
- No \(T_{\mu\nu}\) of a matter field on fixed \(g\).  
- Local free speed remains \(c\) in free frames (Ax7): chart slowing is \(\ell\), not local \(c\) change.

### 5.3 2D sandbox note (B)

In 2D, free Laplacian Green is logarithmic. Exterior “monopole” for \(\psi\) is \(\sim\log r\), while B’s Round-1 kernel used \(1/R\) (3D Coulomb in a 2D integral).  

**Congruence honesty:**  
- Either B works in 3D / quasi-3D,  
- or theory states 2D path-cost target \(\propto\log r\) (or \(\nabla\psi\sim 1/r\) deflection laws carefully),  
- or F1 uses a **screened / Yukawa / slab** free operator whose 2D Green is chosen as free law (must be justified as free medium, not gravity).

Do not silently paste 3D Newton into 2D and call it derived monism.

---

## 6. Kinetic inertia: \(\tfrac12(E_\star/c^2)v^2\) from free response (partial)

### 6.1 Setup

Lock \(L\) moves at slow speed \(v\ll c\) through free medium. Free response \(\psi\) must track the moving source \(J_{\mathrm{lock}}(x-vt)\).

### 6.2 Quasi-static comoving free energy

Suppose free constitutive energy (F1-like)

\[
U_{\mathrm{free}}[\psi]
\;=\;
\frac{\sigma_0}{2}\int|\nabla\psi|^2\,dV
\]

(positive definite free “field” energy of the capacity potential). For a static lock,

\[
U_{\mathrm{free}}^{\mathrm{static}}
\;=\;
\frac{s}{2}\int\psi\rho_b\,dV
\;=\;
\frac{s^2}{8\pi\sigma_0}\frac{E_\star^2}{R_{\mathrm{eff}}}
\quad\text{(order-of-magnitude; self-energy)}.
\]

Self-energy is sensitive to lock size (renormalization). For **inertia**, use the **excess free energy of a slowly moving lock** relative to comoving static configuration.

### 6.3 Liénard–Wiechert / multipole sketch

In the free medium, signals rearrange at most at \(c\). A moving source’s free potential lags; the free energy to order \(v^2\) typically has the structure

\[
U_{\mathrm{free}}(v)
\;=\;
U_{\mathrm{free}}(0)
\;+\;
\frac{1}{2}\mu_{\mathrm{em}} v^2
\;+\;O(v^4/c^2),
\]

where \(\mu_{\mathrm{em}}\) is an **electromagnetic-mass-like** coefficient built from free field energy / \(c^2\).

Classic EM analogy (not dualist matter): electromagnetic mass \(\sim U_{\mathrm{field}}/c^2\).  
Here the “field” is **free continuum \(\psi\)**, and the source energy scale is \(E_\star\).

**Target identification:**

\[
\mu_{\mathrm{em}} \;\stackrel{?}{=}\; \frac{E_\star}{c^2}.
\]

### 6.4 What would force equality (still open)

Full match needs either:

1. **Equivalence postulate (weak):** free medium response is universal in \(E_\star\) only (C EQ5) → \(\mu_{\mathrm{em}}=k E_\star/c^2\) with \(k=1\) fixed by free constitutive normalization (same \(\gamma,s,\sigma_0\) that fix \(G_{\mathrm{eff}}\)); or  
2. **Explicit boost calculation:** expand free action for lock velocity; show coefficient of \(v^2/2\) equals unlock ledger \(E_\star/c^2\).

**Partial result (Round 2):**  
- Inertia **must** be a property of free-response cost of accelerating a lock (relational).  
- Dimensional skeleton \(\mu\sim E_\star/c^2\) is forced once free energy scales with \(E_\star\) and free speed \(c\).  
- **Coefficient 1** (exact \(m=E_\star/c^2\)) is **not yet derived** — it is a killable claim (K1 in `mass_from_locality.md`) for F1 with fixed constitutive normalization.

**Kill if B finds:** push-derived \(m_{\mathrm{inertial}}\) systematically \(\neq E_\star/c^2\) while free law is held fixed and unlock energy is \(E_\star\).

### 6.5 Minimal operational definition for B

\[
m_{\mathrm{inertial}}
\;=\;
\lim_{a\to 0}\frac{F_{\mathrm{ext}}}{a},
\]

where \(F_{\mathrm{ext}}\) is free-medium force needed to hold acceleration \(a\) of the lock center (drag from free \(\psi\) wake). Compare to \(E_\star/c^2\) and to \(M_{\mathrm{ray}}\) from path-cost monopole amplitude.

---

## 7. Relation to B’s \(\Phi\sim\alpha\int\rho_b/R\)

| Aspect | B Round-1 kernel | Monist F1 target |
|--------|------------------|------------------|
| Map | \(\Phi=\alpha\int\rho_b/R\) | \(\psi\) solves free law \(-\nabla\cdot(\sigma\nabla\psi)=s\rho_b\) |
| \(\alpha\) | free fit param | \(\alpha_{\mathrm{eff}}=s\gamma/(4\pi\sigma_0)\) from free constants |
| Sector tag | claimed 1, origin postulated | 1, origin free_relaxation |
| D score | ray-OK; theory incomplete | ray-OK + dynamical export of \(\psi\) |
| O-001 | dualist residue risk | closes residue **if implemented** |

**Bridge:** treat B’s kernel as the **Green’s function evaluation of F1**, not as ontology. Round 2 B success = evolve/solve \(\psi\) from free law and recover \(\alpha_{\mathrm{eff}}\) without hand Poisson gravity.

---

## 8. Kill / keep decision tree

```text
Does path cost come only from free-medium state (ρ_f, ψ, fluxes)?
  NO → dualist (HK1–HK2)
  YES → Does a free dynamical law determine ψ (or equivalent)?
           NO → provisional phenomenology only (F5); not monist theory
           YES → Is the law independent of lock ledger (no free–bound link)?
                    YES → dualist / incomplete (HK5)
                    NO → Is G_eff foreign or emergent from free constants?
                           foreign → dualist (HK3)
                           emergent → MONIST-ELIGIBLE K
                                      then test: exterior ℓ ∝ E_star/(c²r)?
                                      inertia m = E_star/c²?
                                      local c preserved?
```

---

## 9. Open items

1. Exact F1 action principle (variational form of free continuum).  
2. 2D vs 3D Green honesty for B sandboxes.  
3. Close coefficient in \(\tfrac12(E_\star/c^2)v^2\).  
4. Rod/clock C1 interaction with \(\psi\) (does \(\psi\) affect free oscillators?).  
5. Multi-lock linearity and nonlinear free saturation.

---

## 10. Bottom line

- **\(K\) monist** = Green / linear response of **free continuum dynamics** to lock budget conversion.  
- **\(K\) dualist** = independent potential / second-pass Poisson / hand \(\int\rho/R\) without free law.  
- **F1 (free capacity potential)** is the leading analytic candidate: yields \(\ell\sim E_\star/r\) with \(G_{\mathrm{eff}}\) from \((\gamma,s,\sigma_0,c)\).  
- **F5 (B hand kernel)** remains a **shape target**, not monism.  
- Inertia sketch: free-response energy of moving lock \(\Rightarrow\mu\sim E_\star/c^2\); exact factor 1 still open.

Cross-links: `congruence_package_v0.md`, `axioms_v0.md`, `mass_from_locality.md`, C strip S2, B `sandbox_b2_kernel.py`, D degeneracy D-003.
