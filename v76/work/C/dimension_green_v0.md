# Dimension and Free Green Functions (2D log vs 3D \(1/r\)) — v0

**Approach:** C — reverse theoretical  
**Date:** 2026-07-18 (Round 3)  
**Trigger:** O-003; B M2 free Laplace (2D log); D multipole separation; A F1  
**Purpose:** Update path-cost NC for embedding dimension; place B M2 correctly  
  (mechanism-right, dimension-wrong for GR-class exterior).

---

## 0. Core fact (not optional)

For free linear response of **Laplace type** on \(\mathbb{R}^d\),

\[
-\Delta \psi \;=\; s\, J_{\mathrm{lock}}
\quad\text{(or discrete free-graph Laplace)},
\]

the fundamental solution (Green function) \(G_d\) satisfies \(-\Delta G_d=\delta\) and has exterior radial form:

| \(d\) | Free Green \(G_d(r)\) (leading, infinite space) | Exterior multipole of compact source |
|------|--------------------------------------------------|--------------------------------------|
| **2** | \(\displaystyle G_2(r)=\frac{1}{2\pi}\log\frac{1}{r}+\mathrm{const}\) | **logarithmic** growth/decay class |
| **3** | \(\displaystyle G_3(r)=\frac{1}{4\pi r}\) | **\(1/r\) monopole** |
| \(d>3\) | \(\propto r^{2-d}\) | power \(2-d\) |

**Consequence for monist free Laplace (M2 / A F1 class):**  
path-cost contrast inherited from \(\psi\sim G_d * J_{\mathrm{lock}}\) **must** follow the table.  
You cannot get 3D GR-class exterior \(1/r\) from honest 2D free Laplace without smuggling a 3D kernel by hand.

---

## 1. Path-cost NC update (C3 dimension-aware)

### 1.1 Round 1 forced profile (recalled)

Weak-field **Einstein-class / Shapiro-class phenomenological targets** (3D lab world) used:

\[
\ell(r)-\ell_0 \;\sim\; \frac{2 G_{\mathrm{eff}} M}{c^2 r}.
\tag{PC-3D}
\]

That remains the **target for physical 3D monism** aiming at those tests.

### 1.2 Dimension-correct free-response NC

**NC-dim.1 (Green–dimension match).**  
If free response is in the free-Laplace / free-capacity class

\[
-\nabla\cdot\bigl(\sigma\nabla\psi\bigr) = s\,\rho_b
\quad\text{(A F1 / B M2)},
\qquad
\delta\ell = \gamma\,\psi + \cdots,
\]

then the forced exterior path-cost class is:

| Embedding | Forced \(\delta\ell\) class | Matches PC-3D? |
|-----------|----------------------------|----------------|
| 2D continuum / 2D grid | \(\delta\ell \sim A\log(r_0/r)\) (or log growth with IR cutoff) | **No** |
| 3D continuum / 3D grid | \(\delta\ell \sim A/r\) | **Yes** (leading monopole) |

**NC-dim.2 (do not fake dimension).**  
Postulating \(\Phi=\alpha\int\rho_b/R\) on a **2D** domain to force \(1/r\) while free dynamics are 2D is **dimension smuggling** (F5 / soft kernel) — fails NC-K.2 even if rays look “GR-like.”

**NC-dim.3 (score multipole honestly).**  
D/B must report multipole class (`log_2d` vs `inv_r_3d`) from data (e.g. \(\psi(r_2)/\psi(r_1)\) ratios), not only from intention tags.

**NC-dim.4 (IR in 2D).**  
2D log Green is IR-sensitive (grows without bound). Finite-box Neumann/Dirichlet sets the zero of \(\psi\). Congruence tests in 2D must fix gauge/cutoff; **do not** treat 2D log as failed monism — treat as **wrong dimension for PC-3D**.

### 1.3 Revised slogan

> **Forced path cost = free Green of the medium’s dimension, not a universal \(1/r\) pasted from GR.**  
> GR-class \(1/r\) is the **3D** free-Laplace exterior — a dimensional theorem, not an extra gravity law.

---

## 2. Why B M2 is mechanism-right, dimension-wrong

### 2.1 What M2 is (from B design + D scores)

**M2 — free Laplace on free graph / free neighbors:**
- One sector: free capacity / free-site potential relaxes among free neighbors.
- Locks = sinks / Dirichlet holes / bound cells that remove free edges.
- Path cost / chart index from free \(\psi\) (or free travel on remaining graph).
- Tags: `monist_1sector`, dynamical, not hand \(\int\rho_b/R\).

**D R2:** M2 multipole \(\sim\) log; beats dualist_log on Occam; non-iso 1/r dualist loses pure fit on M2 data.

### 2.2 Mechanism-right (PASS monist class)

| Test | M2 |
|------|-----|
| NC-K.1 single sector | **PASS** (if \(\psi\) is free medium) |
| NC-K.2 free-origin of \(K\) | **PASS** (Laplace from free relaxation, not postulated 1/r kernel) |
| NC-K.3 bound = lock of same medium | **PASS** (design) |
| NC-K.4 free capacity reduced at locks | **PASS** (graph / budget) |
| D-K3 postulated soft Poisson | **PASS kill avoided** (unlike R1 hand kernel) |
| O-001 “derived free response” | **PASS class** |

So M2 **clears the dualism-suspect status of Round-1 hand kernel** for the *mechanism*.

### 2.3 Dimension-wrong for Einstein-class exterior (FAIL PC-3D)

| Test | M2 on 2D grid |
|------|----------------|
| Exterior multipole | **log**, not \(1/r\) |
| Match path_cost_profile PC-3D | **FAIL** |
| Congruence J3 as written in R1 (\(1/r\) only) | **FAIL** if J3 requires PC-3D |
| “Wrong theory” vs “wrong dimension” | **Wrong dimension** |

**Verdict (C Round 3):**  
\[
\boxed{\text{M2 = right monist mechanism class; wrong embedding dimension for GR-class }1/r.}
\]
Do **not** kill free-Laplace monism because 2D M2 is log.  
Do **not** claim goal-(2) Einstein-class exterior from 2D M2.

### 2.4 What would make M2 dimension-right

Implement the **same mechanism in 3D**:
- free Laplace / free-capacity \(\psi\) on 3D grid or graph;
- compact lock as free sink / hole;
- measure \(\psi\sim A/r\) for \(r\) in intermediate zone (between lock size and box size);
- rays / Born delay \(\propto M/b\) class from \(\delta\ell=\gamma\psi\);
- tag `monist_1sector`, `K_origin=derived`, `multipole=inv_r_3d`.

That is O-003 critical path for B.

---

## 3. Dualist twin still exists in every dimension

| Dimension | Monist free Laplace | Dualist multipole twin |
|-----------|---------------------|------------------------|
| 2D | M2 log | M3 / Poisson on \(\rho_b\) → same log class |
| 3D | M2-3D \(1/r\) | 3D Poisson on matter → same \(1/r\) class |

**NC-obs (unchanged):** fit alone cannot separate monist free Green from dualist Poisson when multipoles match. Need:
- `sector_tag` / free–bound link / free DOF essential (KB2),
- Occam \(N_{\mathrm{sec}}\),
- optional dynamics (how \(\psi\) is produced).

**3D does not remove dualism risk** — it only aligns multipole with PC-3D. SoftE and sector tags remain mandatory (D R2 false-positive M3).

---

## 4. Monist \(K\) re-affirmation for free Laplace class (NC-K-L)

Free Laplace / free-capacity class is **monist-eligible** when all hold:

| ID | Condition |
|----|-----------|
| **NC-K-L.1** | Equation is free-medium law: \(-\nabla\cdot(\sigma[\rho_f]\nabla\psi)=s\,J_{\mathrm{lock}}\) with \(\sigma,\psi\) free; not \(\nabla^2\Phi_{\mathrm{grav}}=4\pi G\rho_m\) as second sector |
| **NC-K-L.2** | \(J_{\mathrm{lock}}\) from free↔bound conversion / free holes, not foreign dust |
| **NC-K-L.3** | Path cost \(\delta\ell=\mathcal{L}[\psi,\rho_f]\) of free signals; local \(c\) preserved in free frames |
| **NC-K-L.4** | Multipole = Green of **actual** dimension (NC-dim.1) |
| **NC-K-L.5** | Amplitude \(A\propto E_\star\) (or \(\int J\)) with medium constants \(\sigma,s,\gamma\) only — \(G_{\mathrm{eff}}\) emergent |
| **NC-K-L.6** | Dualist rewrite with idle free medium changes physics (free \(\sigma\), free graph matter) |
| **NC-K-L.7** | Tags: `K_origin=derived`, `gravity_solver=none`, `sector=monist_1sector`, `multipole=log_2d\|inv_r_3d` |

**Fails monist Laplace class:**
- Hand \(\Phi=\alpha\int\rho_b/R\) (R1) — NC-K-L.1/2 fail (no free evolution).
- M3 Poisson on \(\rho_b\) labeled monist — dualist twin (softE if mis-tagged).
- 2D code claiming GR \(1/r\) without 3D Green — NC-dim.2 fail.

**Passes (provisional):** M2 2D for **log-class monist demo**; M2 3D for **PC-3D monist candidate**.

---

## 5. Amplitude sketch (3D free Laplace → \(G_{\mathrm{eff}}\))

Stationary free capacity (constant \(\sigma\)):

\[
-\sigma\Delta\psi = s\,\rho_b
\quad\Rightarrow\quad
\psi(\mathbf{x})
=
\frac{s}{4\pi\sigma}
\int\frac{\rho_b(\mathbf{y})}{|\mathbf{x}-\mathbf{y}|}\,dV_y.
\]

For compact lock \(\int\rho_b=E_\star\) (units where ledger density = energy density), exterior:

\[
\psi(r)
\;\approx\;
\frac{s E_\star}{4\pi\sigma\, r}.
\]

If \(\delta\ell=\gamma\psi\) and PC-3D wants \(\delta\ell\approx 2 G_{\mathrm{eff}} M/(c^2 r)\) with \(M=E_\star/c^2\):

\[
\frac{\gamma s E_\star}{4\pi\sigma\, r}
\;=\;
\frac{2 G_{\mathrm{eff}} E_\star}{c^4\, r}
\quad\Rightarrow\quad
\boxed{
G_{\mathrm{eff}}
=
\frac{\gamma s\, c^4}{8\pi\sigma}
}
\quad\text{(schematic; normalize to eikonal convention)}.
\]

**Reverse demand:** \(G_{\mathrm{eff}}\) is **predicted** by free constants \((\sigma,s,\gamma,c)\), not fitted as foreign Newton \(G\) on matter.  
A Round-3 job: lock this normalization to their F1 conventions.

---

## 6. Deflection scaling (Born, isotropic \(\delta\ell\))

For \(\delta\ell\sim A/r\) (3D), straight-line Born deflection scales as

\[
\Delta\theta \sim \frac{\mathrm{const}\cdot A}{b}
\;\propto\;
\frac{M}{b}
\]

(with const depending on whether only “time” or time+space pieces of path cost enter — PPN-lite still open).  

For \(\delta\ell\sim A\log r\) (2D), Born scaling differs (derivative \(\sim A/r\) still gives \(\Delta\theta\sim A/b\) in 2D transverse geometry, but **delay integrals and multipole diagnostics** are log-class; D’s \(\alpha\)-ratio tests separate log vs \(1/r\) profiles of \(\psi\) itself).

**NC-dim.5:** Congruence for PC-3D requires **profile of \(\psi\) or \(\delta\ell\)** \(\sim 1/r\), not only “some deflection \(\propto 1/b\).”

---

## 7. Implications for other approaches

### FOR_A
- Specialize F1 to **3D**: Green \(1/(4\pi r)\); publish \(G_{\mathrm{eff}}(\sigma,s,\gamma,c)\).
- Ax8 target \(\ell-\ell_0\propto E_\star/r\) is **3D free-Laplace theorem**, not independent postulate.
- 2D theory toys should state log exterior explicitly (do not hide).

### FOR_B
- M2 2D: keep as monist mechanism proof; tag `multipole=log_2d`.
- Critical path: **3D M2** — measure \(\psi r \to \mathrm{const}\) in shell; rays; free deficit; no gravity_solver.
- Do not reintroduce hand 1/r on 2D grid as monist win.

### FOR_D
- Separate losses: multipole class vs sector Occam.
- On 3D monist maps: dualist adversary = 3D Poisson on \(\rho_b\) with `sector_tag=dualist`; expect L_fit tie + Occam + tag discipline (M3 lesson).
- Require multipole diagnostic PASS `inv_r_3d` before claiming PC-3D congruence.

---

## 8. Bottom line

| Claim | Status |
|-------|--------|
| Free Green is log in 2D, \(1/r\) in 3D | **Theorem** (Laplace) |
| PC-3D \(\propto 1/r\) needs 3D free response (or non-Laplace law with 3D \(1/r\) Green) | **NC-dim** |
| B M2 | **Mechanism monist; dimension wrong for PC-3D** |
| Hand 2D \(1/r\) kernel | **Dead as monist theory** (R1+R2) |
| Goal (2) Einstein-class | Needs **3D monist free Green** + checklist 3D criteria |

Cross-links: `path_cost_profile_v0.md`, `kernel_dualism_stress_v0.md`, `congruence_checklist_v1.md`, O-003.
