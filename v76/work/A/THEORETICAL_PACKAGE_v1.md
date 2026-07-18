# v76 Theoretical Package v1 — Free-Capacity Monism (F1-3D)

**Approach:** A (forward theoretical)  
**Date:** 2026-07-18  
**Round:** 4  
**Status:** **complete workable theoretical idea** for weak-field free-response monism  
**Congruence tier:** `goal2_PC3D` theory side **MET** (with documented residuals)  
**Numeric twin:** B Round 3 — `monist_3d_1r_pass=True`, parent N=32 SOR \(R^2_{1/r}\approx 0.998\)  
**Depends on:** PROBLEM.md §7, APPROACHES.md §1–2, axioms_v0, free_response_kernel/3d, C NCs/checklist v1, O-005

---

## 0. One-sentence claim

**There is a single continuum whose free capacity potential \(\psi\) responds to bound (mass-form) ledger \(\rho_b\); free-signal path cost is set by \(\psi\); in 3D the exterior multipole is \(\propto 1/r\); mass is bound ledger over \(c^2\); no second gravity sector is required as ontology.**

---

## 1. Primitives (shared v76 language)

| Primitive | Definition in this package |
|-----------|----------------------------|
| **Field / continuum** | Only physical medium: free + bound states of one fabric |
| **Free field** | Medium free to update / carry signals at local locality bound \(c\) |
| **Bound / mass-form** | Medium locked in stable particle-form (lock); not free to leave without unlock |
| **Energy** | Ledger of continuum reconfiguration / causal capacity (frame-dependent integral) |
| **Mass** | \(m = E_\star/c^2\) with \(E_\star =\) rest/unlock bound ledger |
| **\(c\)** | Free-field locality: sup free-influence speed in free frames |
| **Warp** | Nontrivial free-signal geometry in a global chart when local free speed is held \(=c\) around locks |
| **Path cost \(\ell\)** | Free-signal optical cost density derived from free state (not an independent metric field) |

**Not primitives:** independent spacetime metric as distance carrier; foreign Newton \(G\); matter stress-energy *on* a stage as a second substance.

---

## 2. Dynamical core — F1-3D free capacity

### 2.1 State

On a 3D continuum (or large box approximating \(\mathbb{R}^3\)):

\[
\rho_b(\mathbf{x})\ge 0,\quad
\rho_f(\mathbf{x})\ge 0,\quad
\psi(\mathbf{x})\in\mathbb{R},\quad
\sigma(\rho_f)>0.
\]

| Symbol | Meaning |
|--------|---------|
| \(\rho_b\) | Bound budget density (lock occupancy) |
| \(\rho_f\) | Free budget density |
| \(\psi\) | Free capacity potential — **free continuum DOF** |
| \(\sigma\) | Free conductivity / stiffness (constitutive) |
| \(s,\gamma\) | Free constitutive couplings (source strength; path-cost gain) |
| \(c\) | Free locality scale |

### 2.2 Budget identity

**Integral / strong form used in sandboxes:**

\[
\rho_f + \rho_b = \rho_{\mathrm{tot}}
\quad\text{(often }\rho_{\mathrm{tot}}=\rho_0=\mathrm{const}\text{)}.
\tag{B}
\]

**Monist content:** forming bound ledger **uses** free budget. Free fabric is less (as free capacity) where mass-form has formed.

**Weaker admissible form:** integral conservation only, with compact support of \(\rho_b\) and free capacity removed \(\sim E_\star\).

**Kill (local optics):** pointwise (B) + local \(n=n(\rho_f)\) + compact \(\rho_b\) **cannot** give long-range Einstein-class path cost (C no-go; B/D confirmed). Path cost must come from free-response \(\psi\), not local \(\rho_f\) alone.

### 2.3 Free law (quasistatic)

\[
\boxed{
-\nabla\cdot\bigl(\sigma(\rho_f)\,\nabla\psi\bigr)
\;=\;
s\,\rho_b
}
\tag{F1-3D}
\]

Vacuum linearization \(\sigma=\sigma_0\):

\[
-\sigma_0\nabla^2\psi = s\rho_b.
\]

**Equivalent monist form (M2 lift):** \(\nabla^2\psi=0\) on free set; lock = free hole (Dirichlet/flux BC); free capacity removed inside lock.

**Dynamics (optional relaxational):**

\[
\partial_t\psi = \kappa\nabla\cdot\bigl(\sigma\nabla\psi\bigr) + s\,\mathcal{S}[\rho_b].
\]

### 2.4 Bound ledger and mass

\[
E_\star = \int \mathcal{E}[\rho_b]\,dV,
\quad
\mathcal{E}[\rho]=\rho\ \text{(v1 units)},
\quad
M_{\mathrm{ledger}} = \frac{E_\star}{c^2}.
\tag{M}
\]

Unlock energy target: \(E_{\mathrm{unlock}}=E_\star\) (single ledger; C EQ checklist).

### 2.5 Exterior multipole (3D theorem)

For compact \(\rho_b\) with \(E_\star=\int\rho_b\,dV\):

\[
\boxed{
\psi(\mathbf{r})
= \frac{s E_\star}{4\pi\sigma_0\, r}
+ O(r^{-2})
}
\tag{ψ-mono}
\]

**Dimension law:** free Laplace Green is \(\log r\) in 2D and \(1/r\) in 3D.  
B M2 (2D) is monist mechanism with log multipole.  
B R3 (3D) realizes \(\psi\sim 1/r\) dynamically (`monist_3d_1r_pass`).

Diagnostic:

\[
\mathcal{R}_\psi(r) = \frac{4\pi r\,\psi(r)\,\sigma_0}{s E_\star}
\;\xrightarrow{\text{large }r}\; 1.
\]

Parent verify: N=32 SOR multipole prefer \(1/r\), \(R^2_{1/r}\approx 0.998\).

### 2.6 Path cost and rays

\[
\ell = \ell_0 + \gamma\,\psi,
\qquad
n_{\mathrm{chart}}-1 \propto (\ell-\ell_0)
\quad\text{(Born / eikonal; local free speed still }c\text{ in free frames)}.
\]

Exterior:

\[
\ell(r)-\ell_0 = \frac{\gamma s E_\star}{4\pi\sigma_0\, r} + O(r^{-2}).
\]

Match C weak-field path-cost target \(\ell-\ell_0 \approx 2 G_{\mathrm{eff}} M/(c^2 r)\) with \(M=E_\star/c^2\):

\[
\boxed{
G_{\mathrm{eff}}
= \frac{\gamma s\, c^4}{8\pi\,\sigma_0}
}
\tag{Geff}
\]

**All of \(\gamma,s,\sigma_0,c\) are free-medium constitutive parameters.**  
\(G_{\mathrm{eff}}\) is **emergent**, not a second-sector Newton constant.

Rays: integrate free eikonal / Born series on \(\ell[\psi]\) with `gravity_solver=none`.  
Vacuum \(\rho_b\equiv 0\Rightarrow\psi\equiv 0\Rightarrow\) null deflection/delay.

### 2.7 Locality-\(c\) and warp (hinge)

1. Free updates cannot outrun \(c\) (locality).  
2. Observers build rods/clocks from free medium so local free speed \(=c\) isotropic.  
3. Locks rearrange free capacity \(\Rightarrow\psi\neq\mathrm{const}\).  
4. Global charts then see warped free paths — **warp = constant local \(c\) around locks**, not \(T\to G\).

---

## 3. Ontology discipline (monist vs dualist)

### 3.1 Single sector

| Allowed | Forbidden as ontology |
|---------|----------------------|
| \(\psi\) free continuum state | Independent \(\Phi_{\mathrm{grav}}\) unused by free dynamics |
| \(\rho_b\) bound ledger of same medium | Matter \(T_{00}\) on fixed stage as sole mass def |
| One evolution / one solve of free law | Two-pass: place mass → Poisson gravity → rays |
| Emergent \(G_{\mathrm{eff}}\) from free constants | Foreign \(G\) fit with free DOF idle |
| Tags: `phi_origin=free_relaxation`, `sector=1` | Soft-label dualist solve as monist (softE) |

### 3.2 Poisson-form residual (honest)

The **linear vacuum equation** \(-\sigma_0\nabla^2\psi = s\rho_b\) is **mathematically Poisson-shaped**.

**Monism is not “a different PDE shape.”** Monism is:

1. \(\psi\) **is** free fabric capacity (same continuum as \(\rho_f,\rho_b\));  
2. RHS is **budget conversion / lock occupancy**, not a second substance’s \(T_{00}\);  
3. Path cost and rods/clocks are free-medium operational;  
4. No independent gravity sector in the evolution schedule;  
5. Free–bound ledger link (T1 deficit) is part of the theory.

Dualist 3D Poisson \(\Phi\propto M/r\) is a **ray-isomorphic twin**. Separation requires:

- sector tags / Occam (D),  
- free-origin of \(\psi\),  
- free–bound link,  
- optional \(\sigma=\sigma(\rho_f)\) so free budget affects conductivity.

**Residual R1:** ray fit alone never proves monism (D-003; O-005).

---

## 4. What is proven vs residual

### 4.1 Theory claims at v1 (proven as closed workable package)

| Claim | Status |
|-------|--------|
| One continuum free/bound language | Stated; monist-eligible |
| Budget identity + depletion content | Stated; Lean ledger fragment |
| \(c\) free locality + warp schema | Stated |
| F1-3D free law | Written |
| Exterior \(\psi\sim 1/r\) in 3D from free Green | **Theorem** (linear PDE + multipole) |
| \(\ell\sim 1/r\) and \(G_{\mathrm{eff}}(\gamma,s,\sigma_0,c)\) | Derived under \(\ell=\ell_0+\gamma\psi\) |
| Local GRIN long-range dead | C lemma + B/D kill |
| Hand F5 \(\int\rho_b/R\) not monist theory | O/B/D/A consensus |
| 2D free response monist but log | Dimension law |

### 4.2 Numerically congruent (B/D R3 — not theory alone)

| Item | Result |
|------|--------|
| 3D free-capacity SOR/Jacobi | `monist_3d_1r_pass=True` |
| Multipole | prefer \(1/r\), \(R^2\approx 0.998\) (N=32 parent) |
| Free deficit at lock | positive |
| Rays without gravity solver | nonzero defl/delay; vacuum 0 |
| D score | monist_3d wins fit+Occam; softE killed; dualist twin iso on rays |

### 4.3 Residuals (documented — not hidden)

| ID | Residual | Blocks full goal2? |
|----|----------|--------------------|
| **R1** | Poisson-form isomorphism with dualist 3D twin | Ontology: tags/Occam required; not a PDE disproof |
| **R2** | Inertia triad \(m_{\mathrm{inertial}}=E_\star/c^2=M_{\mathrm{ray}}\) not closed | **Yes for goal2_FULL / J5** |
| **R3** | Exact Einstein light factor \(4GM/(c^2b)\) needs space+time path-cost split (PPN-lite) | No for PC-3D multipole tier |
| **R4** | Rod/clock C1 operational construction partial | Soft; Ax7 stated |
| **R5** | Box BC / finite-\(N\) amplitude calibration (ratio to infinite-space analytic) | Numeric systematics |
| **R6** | Full nonlinear \(\sigma(\rho_f)\), multi-lock, cosmology | Out of v1 scope |

### 4.4 Dead branches (do not reopen as monist proof)

| Branch | Verdict |
|--------|---------|
| Local \(n(\rho_f)\) long-range gravity | **DEAD** |
| Hand F5 / postulated \(\alpha\int\rho_b/R\) as monist derivation | **DEAD** (`monist_kernel_failed`) |
| 2D free Laplace as Einstein-class \(1/r\) | **DEAD** (wrong \(d\); mechanism still valid) |
| Soft Einstein / dualist labeled monist | **DEAD** (D softE) |
| Q-ball on fixed grid as monist gravity | **DEAD** (PROBLEM §3) |

---

## 5. Kill conditions (theory package fails if…)

| ID | Kill if… |
|----|----------|
| KP1 | Warp requires independent \(\Phi\) with free \(\psi\) idle |
| KP2 | \(G_{\mathrm{eff}}\) foreign forever and free constitutive law unused |
| KP3 | Free budget never reduced by locks (ledger link broken) |
| KP4 | Only F5 hand \(1/r\) works; free 3D dynamics never prefer \(1/r\) — **refuted by B R3** |
| KP5 | Inertia systematically \(\neq E_\star/c^2\) under honest free-drag (if J5 tested) |
| KP6 | Rod/clock require external non-medium geometry as ontology |

---

## 6. Congruence map (A ↔ B/C/D)

```text
THEORY (this package)          NUMERICS
  F1-3D  −σ∇²ψ = s ρ_b    ←→   B free-capacity 3D SOR
  ψ ~ 1/r exterior        ←→   multipole prefer 1/r, R²~0.998
  ℓ = ℓ0 + γψ             ←→   Born rays / delays
  budget free+bound       ←→   free_deficit_core > 0
  gravity_solver=none     ←→   B/D tags
  G_eff from free const   ←→   amplitude calibration (ongoing)
  sector 1 free origin    ←→   D Occam vs dualist twin
```

C checklist: `goal2_minimal` PASS; `goal2_PC3D` theory T12–T14 **PASS** at v1; joint PC3D **PARTIAL→workable** with R1–R2 residuals (O-005).

---

## 7. Theory-side declaration for goal condition (2)

### 7.1 Orchestrator wording

Goal (2): ≥1 complete theoretical workable idea **AND** ≥1 numerically workable approach, mutually congruent.

### 7.2 A declaration (Round 4)

| Tier | Theory | Numeric (B/D) | Joint |
|------|--------|---------------|-------|
| **goal2_minimal** | **MET** | **MET** (M2+) | **MET** |
| **goal2_PC3D_workable** | **MET** (this package) | **MET** (`monist_3d_1r_pass`) | **MET with residuals R1,R2** |
| **goal2_FULL** (incl. J5 inertia triad) | Theory target stated; coefficient open | **NOT MET** until B inertia | **NOT MET** |

**Theory side of goal (2) as “complete workable idea” for weak-field free-response monism with Einstein-class exterior multipole:**  
**MET** — F1-3D package is closed, killable, dimension-correct, and congruent with B’s 3D free-capacity numerics.

**Not claimed:** full GR derivation; inertia coefficient theorem; that Poisson shape alone is dualism-free without tags; cosmology/DM.

### 7.3 Recommended orchestrator language

> **goal2_PC3D_workable = True** (theory+numerics congruent free-capacity monism in 3D),  
> with documented residuals: dualist Green isomorphism (Occam/tags), inertia triad open (J5).

---

## 8. Pointers

| Doc | Role |
|-----|------|
| `axioms_v0.md` | Ax1–Ax9 locality-first |
| `free_response_kernel_v0.md` | K monist def; F1–F5; kills |
| `free_response_3d_v0.md` | 3D specialization; M2 endorse |
| `congruence_package_v0.md` | T1–T8 tests; export tags |
| `mass_from_locality.md` / `kinetic_inertia_v0.md` | \(m=E/c^2\) sketches |
| `inertia_target_v0.md` | **R4:** non-tautological B inertia protocol |
| C `congruence_checklist_v1.md` | Tier A/B scoring |
| B `outputs/round3_*` | Numeric twin |

---

## 9. Bottom line

**Workable monist package:** one continuum; budget free+bound; free capacity \(\psi\) from F1-3D; path cost \(\ell=\ell_0+\gamma\psi\); exterior \(1/r\) in 3D; \(m=E_\star/c^2\); \(G_{\mathrm{eff}}=\gamma s c^4/(8\pi\sigma_0)\); warp from local \(c\) around locks.

**Dead:** local GRIN long-range; hand F5 as monism; 2D log as Einstein exterior.

**Open for full completeness:** non-tautological inertia (J5); deeper rod/clock; PPN light factor.

**Theory goal-(2) PC3D package: MET (workable complete idea, residuals documented).**
