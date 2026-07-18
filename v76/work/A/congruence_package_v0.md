# Theory–Numerics Congruence Package v0

**Approach:** A (forward theoretical)  
**Date:** 2026-07-18  
**Round:** 2–3 (v0.2: 3D multipole / M2 endorsement)  
**Audience:** B (implement/evolve), C (stress-test), D (score)  
**Depends on:** `free_response_kernel_v0.md`, `free_response_3d_v0.md`, `axioms_v0.md`, C NCs, O-001, O-003

---

## 0. Purpose

State the **minimal theoretical package** that counts as a congruent monist free-response solution for v76 goal condition (2):

> ≥1 complete theoretical workable idea **and** ≥1 numerically workable approach, **mutually congruent**.

Without this package, B’s ray success and D’s fit success can all be dualist-isomorphic (O-001, D-003).

---

## 1. Ontological package (must all hold)

| # | Claim | Formal / operational |
|---|--------|----------------------|
| O1 | One continuum | State = \((\rho_f,\rho_b,\psi,\ldots)\) on event/graph medium; no independent metric actor |
| O2 | Free / bound ledger | \(E_\star=\int\mathcal{E}[\rho_b]\); free capacity reduced by lock formation (integral sense) |
| O3 | \(c\) = free locality | Free updates cannot outrun \(c\); local frames keep free speed \(c\) |
| O4 | Free dynamics | \(\psi\) (or free graph voltages) obey a **single** medium law; path cost \(\ell=\ell[\rho_f,\psi]\) |
| O5 | No second gravity pass | Rays read free state only; no Poisson/Einstein solver as ontology |
| O6 | Emergent \(G_{\mathrm{eff}}\) | Amplitude of exterior \(\ell\) from free constitutive constants, not foreign \(G\) |

**Forbidden as “congruent monism”:** F5-only hand \(\Phi=\alpha\int\rho_b/R\) without free law, even if rays look right.

### 1.1 Round-3 embedding rule (O-003 / B-M2)

| Mechanism | Embedding | Exterior multipole | Congruence role |
|-----------|-----------|--------------------|-----------------|
| **M2 free Laplace / F1** | **2D** | \(\log r\) | Monist mechanism **PASS**; Einstein-class T2 **FAIL** (wrong \(d\)) |
| **M2 free Laplace / F1** | **3D** | \(1/r\) | Primary path to full T2 + C path-cost class |
| Hand kernel / dualist Poisson | any | can fake \(1/r\) | **Not** monist theory |

**A endorsement:** B-M2 is the correct monist free-response **mechanism class**.  
**Requirement:** GR-class exterior \(\ell\propto M/(c^2 r)\) demands **3D free Green**, not postulated 3D Coulomb.

---

## 2. Dynamical package (equations — F1 reference form)

### 2.1 Ledger (minimal)

\[
\rho_f + \rho_b = \rho_{\mathrm{tot}}
\quad\text{or integral conservation with compact support of }\rho_b.
\]

Lock formation: \(\partial_t\rho_b = +\Gamma\), \(\partial_t\rho_f = -\Gamma\) (schematic).

### 2.2 Free response (F1 reference — replaceable by F3 graph)

\[
\partial_t\psi
\;=\;
\kappa\,\nabla\cdot\bigl(\sigma(\rho_f)\nabla\psi\bigr)
\;+\;
s\,\mathcal{S}[\rho_b]
\quad\text{or quasistatic }\;
-\nabla\cdot(\sigma\nabla\psi)=s\rho_b.
\]

Constitutive: \(\sigma(\rho_f)>0\), \(s,\kappa\) free constants; \(\mathcal{S}[\rho_b]=\rho_b\) (frozen) or \(\partial_t\rho_b\) (causal).

**3D exterior (required for Einstein-class T2):**

\[
\psi(\mathbf{r})=\frac{s E_\star}{4\pi\sigma_0\, r}+O(r^{-2}),
\qquad
G_{\mathrm{eff}}=\frac{\gamma s\, c^4}{8\pi\sigma_0}
\quad\text{(with }\ell-\ell_0=\gamma\psi\text{)}.
\]

See `free_response_3d_v0.md`. Alternative monist form: \(\nabla^2\psi=0\) on free set with lock as Dirichlet/flux hole (M2 lift).

### 2.3 Path cost and rays

\[
\ell = \ell_0 + \gamma\,\psi
\quad\text{(linear; alternatives allowed if exported)},
\]

\[
\frac{d\mathbf{x}}{dt} = \frac{c}{n_{\mathrm{chart}}}\,\hat{\mathbf{k}},
\quad
n_{\mathrm{chart}}-1 \propto (\ell-\ell_0)\ \text{or eikonal from free null structure}.
\]

Local frames: free proper speed \(=c\) (Ax7).

### 2.4 Mass ledger

\[
M_{\mathrm{ledger}} = \frac{E_\star}{c^2},\qquad E_\star = \int\mathcal{E}[\rho_b]\,dV.
\]

---

## 3. Observable package (what B must export)

Export tags (machine-readable JSON recommended):

```text
sector_count: 1 | 2
sector_tag: monist_1sector | dualist_2sector
phi_origin: free_relaxation | free_diffusion | graph_green | free_laplace_3d
            | postulated_kernel | poisson_solve
gravity_solver: none | poisson | einstein
budget_identity: local | integral | none
embedding_dim: 2 | 3
multipole_class: compact | log_2d | 1/r | other   # MEASURED fit, not aspiration
```

### 3.1 Required fields

| Field | Meaning |
|-------|---------|
| \(\rho_b(x)\), \(\rho_f(x)\) | Ledger densities |
| \(\psi(x)\) or free voltages | Free-state field (not only n) |
| \(\ell(x)\) or \(n(x)\) | Path cost / chart index |
| \(E_\star\), \(M_{\mathrm{ledger}}\) | Lock ledger |
| ray table \((b,\Delta\theta,\Delta t)\) | Observables |
| vacuum control rays | Null warp check |
| free constitutive constants | \(\sigma_0,s,\gamma,c,D,\ldots\) used |

### 3.2 Optional but high-value

| Field | Meaning |
|-------|---------|
| \(m_{\mathrm{inertial}}\) | Push test (F_ext/a) |
| \(E_{\mathrm{unlock}}\) | If unlock dynamics exist |
| exterior profile \(\ell(r)\cdot r / E_\star\) | 3D monopole ratio \(\to\gamma s/(4\pi\sigma_0)\) |
| \(\mathcal{R}_\psi(r)=4\pi r\psi\sigma_0/(s E_\star)\) | Should \(\to 1\) in 3D free Laplace |
| formation history | \(\rho_b(t)\) |

### 3.3 Tag rules (for D)

| phi_origin | N_sec for D Occam | Notes |
|------------|-------------------|--------|
| free_relaxation / free_diffusion / graph_green / free_laplace_3d | 1 | Monist candidate (M2 class) |
| postulated_kernel | 1_claimed | Phenomenology; **not** theory congruence win (`monist_kernel_failed`) |
| poisson_solve | 2 | Dualist adversary (M3 class; multipole twin possible) |

---

## 4. Quantitative congruence tests

### T1 — Free deficit link

\[
\int(\rho_0-\rho_f)\,dV \;\ge\; f\,E_\star
\quad(f\sim O(1)\ \text{model-dependent; not optional zero}).
\]

**Fail:** lock with \(E_\star>0\) and zero free-capacity change (HK5).

### T2 — Exterior path-cost multipole (dimension-aware)

**T2-3D (Einstein-class / C path-cost target):**

\[
\lim_{r\to\infty} r\,(\ell(r)-\ell_0)
\;=\;
\frac{2 G_{\mathrm{eff}} M_{\mathrm{ledger}}}{c^2}
\;=\;
\frac{\gamma s E_\star}{4\pi\sigma_0},
\qquad
M_{\mathrm{ledger}}=\frac{E_\star}{c^2},
\quad
G_{\mathrm{eff}}=\frac{\gamma s c^4}{8\pi\sigma_0}.
\]

Equivalent free-state form: \(\mathcal{R}_\psi(r)\to 1\) with \(\psi\sim s E_\star/(4\pi\sigma_0 r)\).

**T2-2D (M2 class — monist mechanism only):**

\[
\ell(r)-\ell_0 \;\sim\; A\log(R_{\mathrm{box}}/r)
\quad\text{(document }A\text{; do NOT claim C }1/r\text{ monopole)}.
\]

**Pass rules:**

| Claim | Need |
|-------|------|
| “Dynamical free-response monism works” | T2-2D **or** T2-3D with free origin |
| “Congruent with C Einstein-class exterior” | **T2-3D only** + free origin + O1–O6 |
| Hand \(1/r\) without free solve | **Always fail** theory congruence |

**Fail:** exterior flat while claiming long-range; \(M_{\mathrm{ray}}\) free-fit with \(M\neq M_{\mathrm{ledger}}\) by design; 2D log sold as \(1/r\); F5 hand Coulomb.

### T3 — No second solver

`gravity_solver=none` and rays use only exported free state.

**Fail:** HK1/HK4.

### T4 — Vacuum / unlock controls

- \(\rho_b\equiv 0\Rightarrow\Delta\theta=\Delta t=0\).  
- Unlocking lock \(\Rightarrow\) exterior \(\ell\) collapses toward \(\ell_0\) (when free dynamics relax).

### T5 — Inertia triad (when available)

\[
m_{\mathrm{inertial}} \;\stackrel{?}{=}\; \frac{E_{\mathrm{unlock}}}{c^2} \;\stackrel{?}{=}\; M_{\mathrm{ray}}.
\]

**Pass band:** relative error within sandbox tolerance (declare \(\varepsilon_{\mathrm{triad}}\)).  
**Fail:** systematic split (C EQ checklist).

### T6 — Dualist separation (D)

Against dualist adversary:

- **2D:** M3 Poisson log twin of M2 — \(L_{\mathrm{fit}}\) may tie → Occam + `sector_tag` (D-008).  
- **3D:** dualist \(\Phi\propto M/r\) twin of monist free Laplace \(1/r\) — same rule.  
- Non-isomorphic dualist (wrong multipole family) separates on pure \(L_{\mathrm{fit}}\).

Congruent monism requires T1–T4 on **forward free-state data**, not rays alone.

### T7 — Local-optics no-go compliance

Local \(n(\rho_f)\) compact channel may pass T1 and compact rays but **must not** claim T2-3D. Tag `channel=compact_grin` vs `channel=free_response`.

### T8 — Dimension honesty (Round 3)

| embedding_dim | multipole_class allowed for C long-range claim |
|---------------|--------------------------------------------------|
| 2 | log_2d only — monist mechanism OK |
| 3 | 1/r for Einstein-class claim |

**Fail T8:** `embedding_dim=2` + claim multipole_class=1/r without 3D measure; or force-fit 1/r on log data (D R2).

---

## 5. Kill conditions for the package (when monism fails)

| Kill | Condition |
|------|-----------|
| K-pack-1 | Every free law that hits T2 requires an independent \(\Phi\) unused by free dynamics (HK1–2) |
| K-pack-2 | \(G_{\mathrm{eff}}\) cannot be reduced to free constitutive constants (HK3) |
| K-pack-3 | Inertia triad fails under honest free dynamics (K1 mass_from_locality) |
| K-pack-4 | Only postulated_kernel (F5) ever works; free evolution never recovers multipole |
| K-pack-5 | Rod/clock construction forces second substance (C1 fail) |
| K-pack-6 | 3D free Laplace never yields measured \(1/r\) (then Einstein-class monist free-response branch fails) |

Any one sustained kill → report goal condition (3) partial for kernel branch.

---

## 6. Dimensional / unit conventions (shared)

Prefer code units \(c=1\), \(\rho_0=1\). Report:

- how \(E_\star\) relates to \(\int\rho_b\);  
- 2D vs 3D Green used;  
- definition of \(n_{\mathrm{chart}}\) vs local free speed \(c\).

---

## 7. Minimal success story (what “congruent” looks like)

### 7.1 Partial success already (Round 2) — mechanism monism

1. **A:** F1/F3 written; 2D log honesty.  
2. **B-M2:** free Laplace dynamical; `monist_1sector`; exterior **log**.  
3. **D:** M2 beats dualist_log on Occam; non-iso separates multipole families.  
4. **C:** schema free-response conditionally monist (NC-K).

⇒ Dynamical free-response monism is **alive** (goal (2) still incomplete for GR-class).

### 7.2 Full free-response congruence (Round 3 target)

1. **A:** F1-3D + \(G_{\mathrm{eff}}=\gamma s c^4/(8\pi\sigma_0)\) (`free_response_3d_v0.md`).  
2. **B:** 3D free Laplace / free capacity \(\psi\); \(\mathcal{R}_\psi\to 1\); rays without second solver; tags `embedding_dim=3`, `multipole_class=1/r`, `phi_origin=free_laplace_3d`.  
3. **C:** NC checklist pass for 3D monist maps; dimension note in NCs.  
4. **D:** monist 3D free Laplace wins fit+Occam vs dualist 3D Poisson twin; \(M_{\mathrm{ledger}}\) recovered.

Then goal condition (2) can be argued for the free-response branch at Einstein-class exterior order.

---

## 8. Explicit non-success / partial status

| Item | Status |
|------|--------|
| Local GRIN B2-lite / M1 | Compact monist only — not T2-3D |
| B hand kernel R1 / F5 | `monist_kernel_failed` — not theory monism |
| D monist_kernel R1 win | Synthetic F5 — not dynamical free law |
| **B-M2 free Laplace 2D** | **Monist mechanism PASS**; T2-3D **FAIL** (log) |
| B-M3 dualist 2D Poisson | Multipole twin of M2; Occam lose |
| 3D free Laplace monist | **OPEN** — Round 3 critical path |

---

## 9. FOR handoffs (Round 3)

**FOR_B:** Implement **3D** free Laplace / F1; locks as sources or free holes; measure \(\psi\sim 1/r\); export tags above; dualist 3D Poisson control optional; avoid F5 hand \(M/r\); inertia only if non-tautological free-drag.

**FOR_C:** Dimension note 2D log vs 3D 1/r in NC; update congruence checklist: T2-3D + free origin for pass; keep S2 stress on F1-3D (A: still monist if free capacity).

**FOR_D:** Score B 3D maps; dualist 3D Poisson adversary with sector_tag; do not award monism to postulated 1/r; expect M2-class 3D monist to tie multipole with dualist and win Occam.

---

## 10. Bottom line

Congruence = **same free law** on paper (A) and in sandbox (B), in the **right embedding dimension** for the claimed multipole, scored with **ledger + sector + path cost** (D), reverse-checked (C).

**M2 = monist mechanism. 3D = GR-class exterior. Hand 1/r = dead as monism.**
