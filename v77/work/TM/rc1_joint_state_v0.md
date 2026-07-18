# RC1 Joint State v0 — Co-field \(\psi\) + Dynamical Maxwell

**Agent:** TM (Theory — Matter)  
**Date:** 2026-07-19  
**Phase / round:** Phase 2 · **P2-R3**  
**Checkpoint:** **CP-RC1-SPEC**  
**Status:** frozen for NM implementation (co-agree: TE, NE, NM)  
**Depends on:**  
- **CP-M1-NUM ADOPTED** — `work/NE/outputs/m1_result.json` (`m1_claim=true`; parent O-011)  
- `dual_channel_v0.md` (STATE-DC / DUAL-1 baseline)  
- `lock_ontology_v0.md` (S0–S3; Supp)  
- TE `full_maxwell_monist_v0.md` FROZEN; `dual_channel_joint_v0.md` JC1–JC7  
- NE `Maxwell2D` API (`sandbox_m1_2d.py`)  
- CAMPAIGN_MAP §1 RC1: *one grid; dynamical Maxwell + \(\psi\); fixed multi-locks*

**Hard rule (map):** RC1 **REJECT** if only Maxwell-lite \(\Phi\) without dynamical \((\mathbf{E},\mathbf{B})\).

---

## 0. Claim (RC1 bar)

**On a single continuum grid, fixed multi-locks carry \((\rho_b,\rho_Q)\) that source free-capacity \(\psi\) (F1) and free-gauge \((\mathbf{E},\mathbf{B})\) via dynamical Maxwell (Yee M1-class), with shared \(c\), budget identity, sibling independence, and force diagnostics \(F^\psi + Q\mathbf{E}\) (locks do not move).**

When all RC1 gates PASS: `rc1_claim = true` → unlocks **CP-RC1-NUM** co-agreement and Phase-2 STOP eligibility (with M1).

---

## 1. Joint state vector (one grid)

### 1.1 Continuum fields

\[
\boxed{
\mathcal{S}_{\mathrm{RC1}}
=
\bigl\{\,
\rho_f,\;
\rho_b,\;
\rho_Q,\;
\mathbf{J}_Q,\;
\psi,\;
\mathbf{E},\;
\mathbf{B},\;
\sigma_0,\; s,\; \gamma,\;
\varepsilon,\; \mu,\;
c
\,\bigr\}
}
\tag{STATE-RC1}
\]

| Symbol | Type | Role |
|--------|------|------|
| \(\rho_b\ge 0\) | scalar | Bound / mass-form ledger density |
| \(\rho_Q\in\mathbb{R}\) | scalar | Gauge-charge ledger density |
| \(\mathbf{J}_Q\) | vector | Gauge current (RC1 default \(\mathbf{0}\) for fixed locks) |
| \(\rho_f\) | scalar | Free budget; \(\rho_f+\rho_b=\rho_{\mathrm{tot}}\) (strong or integral) |
| \(\psi\) | scalar | Free-capacity potential (C-ψ) |
| \(\mathbf{E},\mathbf{B}\) | vectors | Free-gauge dynamical fields (C-A) — **not** optional diagnostics of \(\Phi\) alone |
| \(\sigma_0,s,\gamma\) | const | C-ψ constitutive (v76 / dual_channel) |
| \(\varepsilon,\mu\) | const | C-A constitutive; \(c=1/\sqrt{\varepsilon\mu}\) |

**Required export tags:**

```text
sector=1
dual_channel=1
phi_origin=free_capacity_f1
E_origin=free_maxwell_full        # NOT free_maxwell_lite / free_gauge_lite for RC1 bar
em_solver=free_maxwell_yee_m1     # or import Maxwell2D
gravity_solver=none
rc1_locks_fixed=true
embedding_dim=2                   # default RC1 (M1 is 2D); 3D deferred to RC3
C_LOCAL = c = 1/sqrt(eps*mu)
```

### 1.2 Fixed locks (no motion)

\[
L_i=\bigl(\mathbf{X}_i,\; E_{\star,i},\; Q_i,\; \mathcal{C}_i\bigr),
\qquad
\mathbf{V}_i\equiv\mathbf{0},\quad
\partial_t\mathbf{X}_i=\mathbf{0}.
\tag{LOCK-RC1}
\]

Reconstruction on the grid:

\[
\rho_b(\mathbf{x})=\sum_i E_{\star,i}\,f_i(\mathbf{x}-\mathbf{X}_i),\qquad
\rho_Q(\mathbf{x})=\sum_i Q_i\,g_i(\mathbf{x}-\mathbf{X}_i),
\]

with compact kernels \(f_i,g_i\ge 0\), \(\int f_i=\int g_i=1\), and **same centers** \(\mathbf{X}_i\).

**Support law:**

\[
\boxed{\mathrm{supp}(|\rho_Q|)\subseteq\mathrm{supp}(\rho_b)}
\tag{Supp}
\]

| Forbidden | Why |
|-----------|-----|
| \(\rho_Q\equiv\rho_b\) always | JC3 / sibling collapse |
| \(\psi\equiv\Phi\) or \(\psi\) fitted to \(|\mathbf{E}|\) as identity | JC6 / K-RC1 |
| Charge blob with \(\rho_b=0\) | Supp violation (Matter default) |
| **Φ-only Coulomb** as sole EM path for `rc1_claim` | CAMPAIGN_MAP hard rule |

### 1.3 Recommended units / grid

| Choice | Default | Note |
|--------|---------|------|
| \(c,\varepsilon,\mu\) | \(1,1,1\) | Plus optional off-unit regression |
| \(\sigma_0,s,\gamma\) | as R2 dual sandbox | Do not retune to fake EM |
| Grid | Same Cartesian mesh for \(\psi\) and Yee EM | **One medium** |
| Dim | **2D** (align M1 `Maxwell2D`) | 3D multipole aesthetic deferred |
| Locks | \(N_L\ge 2\) Gaussians | Cases: neutral, same \(Q\), opposite \(Q\) |
| \(N\) | \(\ge 64\) smoke; \(\ge 128\) preferred | Document |

---

## 2. Evolution laws (RC1)

### 2.1 C-ψ — free capacity (fixed sources)

Quasistatic F1 (default RC1; TD T0 acceptable for fixed locks):

\[
\boxed{
-\nabla\cdot\bigl(\sigma_0\nabla\psi\bigr)=s\,\rho_b
}
\tag{F1}
\]

Optional TD T1 relaxational \(\partial_t\psi=\kappa\nabla\cdot(\sigma\nabla\psi)+s\rho_b\) **only if** static limit recovers F1 and does not block gates — not required for first `rc1_claim`.

Path-cost diagnostic (optional rays): \(\ell=\ell_0+\gamma\psi\).

### 2.2 C-A — dynamical Maxwell (**required for RC1 bar**)

\[
\boxed{
\begin{aligned}
\nabla\cdot\mathbf{B}&=0,\\
\nabla\times\mathbf{E}+\partial_t\mathbf{B}&=\mathbf{0},\\
\nabla\cdot(\varepsilon\mathbf{E})&=\rho_Q,\\
\nabla\times(\mathbf{B}/\mu)-\partial_t(\varepsilon\mathbf{E})&=\mathbf{J}_Q,\\
\partial_t\rho_Q+\nabla\cdot\mathbf{J}_Q&=0,\\
c&=1/\sqrt{\varepsilon\mu}.
\end{aligned}
}
\tag{MX-dyn}
\]

**RC1 fixed locks:** \(\mathbf{J}_Q=\mathbf{0}\), \(\partial_t\rho_Q=0\), \(\rho_Q=\rho_Q(\mathbf{x})\) prescribed.  
Dynamical content = free Maxwell **update** of \((\mathbf{E},\mathbf{B})\) consistent with Gauss (initial projection / Cont-safe setup as in M1-G5), **not** a single static Poisson solve offered as the only EM field.

**Allowed as diagnostics / warm-start (not substitute for MX-dyn):**

- Electrostatic \(\Phi\) with \(\mathbf{E}_0=-\nabla\Phi\) as **initial** \(\mathbf{E}\); then Yee steps run.  
- Snapshot \(\lvert\mathbf{E}\rvert\) Coulomb multipole after relax / many steps.

**Forbidden for `rc1_claim=true`:**

- EM fields exist only as \(\Phi\) from \(-\varepsilon\nabla^2\Phi=\rho_Q\) with **no** dynamical \((\mathbf{E},\mathbf{B})\) array evolved.  
- Tag `E_origin=free_maxwell_lite` as sole channel (R2 dual lite is **regression only**, not RC1 pass).

### 2.3 Shared locality (JC1)

\[
c_{\mathrm{path}} \equiv c_{\mathrm{EM}} \equiv \frac{1}{\sqrt{\varepsilon\mu}} \equiv C_{\mathrm{LOCAL}}.
\tag{c-share}
\]

If NE/ND free-wave module co-runs: \(c_{\mathrm{wave}}\approx c\) within M1 tolerance; else constitutive equality + M1 inheritance is enough for RC1-SPEC numeric (document `c_wave=inherited_m1`).

### 2.4 Joint static multipole targets (after EM settled / averaged)

| Channel | Exterior (2D note) | 3D reference (stretch) |
|---------|--------------------|-------------------------|
| \(\psi\) | Free Green (2D log or large-box F1); mass monopole from \(E_\star\) | \(\sim 1/r\) |
| \(\mathbf{E}\) | Coulomb-class from \(Q\); dynamical consistency with Gauss | \(\sim 1/r^2\) |

**Sibling tests (JC6–JC7):** neutral \(Q=0\Rightarrow\mathbf{E}\approx 0\), \(\psi\neq 0\); opposite \(Q\) does not flip \(\psi\) monopole sign structure.

---

## 3. NE Maxwell module API (consume, do not re-invent)

From M1 export (`m1_result.json` `api`, `sandbox_m1_2d.py`):

```text
class Maxwell2D:
    step(rho_Q=None, Jx=None, Jy=None, Jz=None) -> None
    fields() -> Ex, Ey, Ez, Hx, Hy, Hz   # B = μ H
```

| NM duty | Spec |
|---------|------|
| Own | \(\rho_b,\rho_f,\psi\), lock placement, F1 solve, force diagnostics, RC1 gates export |
| Call | `Maxwell2D.step(rho_Q=..., J*=0)` each EM substep (or documented batch) |
| Read | \(\mathbf{E}\) (and \(\mathbf{B}/\mathbf{H}\)) from `fields()` for forces and multipoles |
| Must not | Replace `step` with pure \(\Phi\) Poisson for claim path |
| May | Import NE module by path; copy-minimal if packaging requires — **same Yee physics** |

**Initial \(\mathbf{E}\):** Poisson projection of \(\rho_Q\) onto TE\(_z\) or TM\(_z\) sector **or** M1 Cont-safe projector; document method. After \(N_{\mathrm{eq}}\) steps, report residual Gauss and use \(\mathbf{E}\) for \(F^{\mathrm{EM}}\).

**ψ grid alignment:** same \(N_x,N_y,L_x,L_y\) as Maxwell2D (cell-centered \(\psi\) OK if interpolation to force points documented).

---

## 4. Forces on **fixed** locks (diagnostics only)

RC1 does **not** update \(\mathbf{X}_i\). Forces are **measured**, not integrated.

### 4.1 Path-cost channel

\[
U_\psi=\frac{s}{2}\int\psi\,\rho_b\,dV,
\qquad
\mathbf{F}_i^{\psi}=-\nabla_{\mathbf{X}_i}U_\psi\Big|_{\mathrm{no\,self}}
\tag{F-ψ}
\]

(Equivalent finite-difference \(-\partial U/\partial R\) for two-lock separation \(R\).)  
Expect: always **attractive** for \(E_{\star,i}>0\).

Tag: `force_closure_psi=virtual_work_energy_gradient_v0`.

### 4.2 EM channel (dynamical \(\mathbf{E}\))

\[
\boxed{
\mathbf{F}_i^{\mathrm{EM}}
=
Q_i\,\mathbf{E}_{-i}(\mathbf{X}_i)
\quad\text{or}\quad
\int \rho_Q^{(i)}\,\mathbf{E}_{-i}\,dV
}
\tag{F-E}
\]

Self-field excluded / soft-kernel regularized.  
With \(\mathbf{V}_i=\mathbf{0}\), Lorentz magnetic term vanishes; full \(Q(\mathbf{E}+\mathbf{V}\times\mathbf{B})\) is **RC2**.

Tag: `force_closure_em=qE_from_dynamical_E_v0`  
**Reject tag for RC1 claim:** `force_closure_em=coulomb_phi_only` as sole method.

### 4.3 Total diagnostic

\[
\mathbf{F}_i^{\mathrm{diag}}=\mathbf{F}_i^{\psi}+\mathbf{F}_i^{\mathrm{EM}}.
\tag{F-tot}
\]

| Config | \(F^\psi\) | \(F^{\mathrm{EM}}\) |
|--------|-----------|---------------------|
| Neutral–neutral \(Q=0\) | attract | \(\approx 0\) |
| Like charge | attract | repel |
| Opposite charge | attract | attract |
| Vacuum (no locks) | 0 | 0 |

### 4.4 Energy split (diagnostic)

Report separately (no double-count claim required beyond tags):

| Ledger | Form |
|--------|------|
| Bound | \(E_\star=\int\rho_b\) |
| Free-capacity stress | \(U_\psi\) (or \(\frac{\sigma_0}{2}\int|\nabla\psi|^2\) equivalent) |
| Free-gauge | \(U_{\mathrm{EM}}=\int\bigl(\varepsilon|\mathbf{E}|^2/2+|\mathbf{B}|^2/(2\mu)\bigr)\) |

Vacuum control: all small without sources.

---

## 5. What RC1 is / is not

| RC1 **is** | RC1 **is not** |
|------------|----------------|
| Co-field on one grid | Lock motion / orbits (→ RC2) |
| Dynamical Maxwell + F1 ψ | Φ-only dual lite rebranded |
| Force **diagnostics** | \(F=ma\) integration |
| Closes R-compose residual at co-field tier | Full 3D co-evolution (RC3) |
| Inherits M1 EM quality | Re-proving all M1-G* inside NM |

---

## 6. FOR_NM — RC1 kill-gates

**Sandbox target:** `work/NM/sandbox_rc1_cofield.py`  
**Export:** `work/NM/outputs/rc1_result.json`, `rc1_tm_gates.json`  
**Claim:** `rc1_claim = true` iff **RG0 ∧ RG1 ∧ RG2 ∧ RG3 ∧ RG4 ∧ RG5 ∧ RG6** (mandatory).

### RG0 — Joint medium / no Φ-only (P0)

| ID | Check | PASS |
|----|-------|------|
| RG0.1 | Single run owns \(\psi\) and dynamical \(\mathbf{E},\mathbf{B}\) (or H) arrays | Both present; EM from Yee/`Maxwell2D.step` |
| RG0.2 | `E_origin=free_maxwell_full` (or `free_maxwell_yee_m1`) | Not lite-only |
| RG0.3 | Budget: free deficit or \(\rho_f+\rho_b=\rho_{\mathrm{tot}}\) | Positive free deficit at cores |
| RG0.4 | Supp | \(\max|\rho_Q|\) only where \(\rho_b>\epsilon\) |
| RG0.5 | Regression optional | R2 dual lite may run as `regression_dual_lite=PASS` but **does not** set `rc1_claim` |

**FAIL hard:** only \(\Phi\) Poisson for EM.

### RG1 — Dual sources / multipoles (P0)

| ID | Check | PASS |
|----|-------|------|
| RG1.1 | \(\psi\) responds to \(\rho_b\); vacuum \(\rho_b=0\Rightarrow\psi\approx 0\) | Floor |
| RG1.2 | \(\mathbf{E}\) responds to \(\rho_Q\); vacuum \(\rho_Q=0\Rightarrow\mathbf{E}\approx 0\) after settle | Floor |
| RG1.3 | Sibling independence | Scale \(E_\star\)@fixed \(Q\) changes \(\psi\) not \(Q\); scale \(Q\)@fixed \(E_\star\) changes \(\mathbf{E}\) not \(E_\star\) |
| RG1.4 | \(\ge 2\) locks | \(N_L\ge 2\) |

### RG2 — Shared \(c\) (P0)

| ID | Check | PASS |
|----|-------|------|
| RG2.1 | Config \(c_{\mathrm{path}}=1/\sqrt{\varepsilon\mu}\) | Exact match in params |
| RG2.2 | M1 inheritance or joint wave spot-check | `shared_c` block in JSON; \(C_{\mathrm{LOCAL}}\) reported |

### RG3 — Force taxonomy diagnostics (P0)

| ID | Config | PASS |
|----|--------|------|
| RG3.1 | Neutral–neutral | \(F^{\mathrm{EM}}\approx 0\); \(F^\psi\) attractive (signed) |
| RG3.2 | Like \(Q\) | \(F^{\mathrm{EM}}\) repel; \(F^\psi\) attract |
| RG3.3 | Opposite \(Q\) | both attract (static \(\mathbf{E}\)) |
| RG3.4 | Vacuum | forces 0 |
| RG3.5 | \(F^{\mathrm{EM}}\) from **dynamical** \(\mathbf{E}\) | tag `qE_from_dynamical_E_v0` |

### RG4 — Sibling non-collapse (P0 kill)

| ID | Check | PASS / FAIL |
|----|-------|-------------|
| RG4.1 | Neutral mass does not source Coulomb monopole \(\mathbf{E}\) | Pass if \(\max|\mathbf{E}|\) at floor |
| RG4.2 | Opposite \(Q\) does not invert \(\psi\) mass monopole | Pass |
| RG4.3 | Forced global \(\psi\approx a|\mathbf{E}|\) or \(\psi\approx a\Phi\) fails when \(Q/M\) varies | Pass if residual large |

**Fail RG4 → report K-RC1 stress to TU** (not auto V77-K alone).

### RG5 — Dynamical Maxwell honesty (P0)

| ID | Check | PASS |
|----|-------|------|
| RG5.1 | \(\ge N_{\mathrm{step}}\) Maxwell `step` calls with fixed \(\rho_Q\) (recommend \(\ge 20\)) | Count in JSON |
| RG5.2 | Gauss residual documented (init and after steps); non-catastrophic | Same order as M1-G5 spirit |
| RG5.3 | \(\max|\mathbf{B}|\) reported (may be small for electrostatic settle) | Field present in state |
| RG5.4 | Optional: incomplete-Ampère adversary **not** required inside NM if M1-G8 inherited | Tag `m1_claim_inherited=true` |

### RG6 — Vacuum + energy split (P1 but mandatory for claim)

| ID | Check | PASS |
|----|-------|------|
| RG6.1 | No locks: \(\psi,\mathbf{E},\mathbf{B}\) at floor after steps | Vacuum |
| RG6.2 | With locks: report \(U_\psi\), \(U_{\mathrm{EM}}\), \(E_\star\) | Finite; no NaN |
| RG6.3 | Hierarchy optional soft | If run: \(|F^E/F^\psi|\gg 1\) without \(s=0\) — soft for first green |

### Optional (not blocking first `rc1_claim`)

| ID | Content |
|----|---------|
| RG7 | Path-cost rays / delay nonzero with \(\rho_b\) |
| RG8 | Off-unit \(c=1/2\) joint smoke |
| RG9 | Import NE `m1_claim` regression subprocess |

### Export schema (minimum)

```json
{
  "checkpoint": "CP-RC1-NUM",
  "rc1_claim": true,
  "tags": {
    "E_origin": "free_maxwell_full",
    "phi_origin": "free_capacity_f1",
    "rc1_locks_fixed": true,
    "gravity_solver": "none"
  },
  "gates": {
    "RG0": true, "RG1": true, "RG2": true,
    "RG3": true, "RG4": true, "RG5": true, "RG6": true
  },
  "forces": { "neutral": {}, "like": {}, "opposite": {} },
  "shared_c": { "C_LOCAL": 1.0, "eps": 1.0, "mu": 1.0 },
  "maxwell": {
    "api": "Maxwell2D",
    "n_steps": 0,
    "gauss_rel": null
  },
  "m1_claim_inherited": true
}
```

```text
FOR_TM: RG0=... RG1=... ... RG6=... rc1_claim=true|false
FOR_TU: D-RC1-cofield status=...
FOR_NE: maxwell_api_ok=...
```

---

## 7. Demo IDs

| Demo ID | Meaning | Owner |
|---------|---------|-------|
| **D-RC1-cofield** | Primary co-field \(\psi\)+dynamical Maxwell | NM (+NE API) |
| **D-MAT-rc1** | Alias Matter | NM |
| D-DUAL-channel (R2) | Static lite joint — **regression**, not RC1 pass | NM |
| D-EM-M1-suite | M1 inherited | NE |

---

## 8. Co-agreement roles (CP-RC1-SPEC)

| Agent | Spec duty | Stamp meaning |
|-------|-----------|---------------|
| **TM** | This document: state, Supp, F1, forces, RG* | **ADOPT** when implementable & monist |
| **TE** | JC1–JC7 / no \(\psi\equiv\Phi\); MX-dyn = FROZEN M1–M4 | ADOPT if no ontology fork |
| **NE** | `Maxwell2D` API stable; `rc1_ready` | ADOPT if API matches §3 |
| **NM** | Implementability of §1–6 | ADOPT if can code; else DEFER with gaps |
| **TD** | Optional note: T0 F1 OK for fixed locks | note |
| **TU** | Board flip when dual stamps land | board |

**Downstream:** CP-RC1-NUM after NM `rc1_claim=true` + TM/TE/NE/NM stamps.

---

## 9. Kill conditions (RC1)

| ID | Kill if… |
|----|----------|
| **K-RC1-a** | Co-field only works as two non-shared grids with idle free budget |
| **K-RC1-b** | Shared \(c\) fails (config or measured wave) on joint run |
| **K-RC1-c** | Forces require \(\psi\equiv\Phi\) identification |
| **K-RC1-d** | Only Φ-lite path ever passes gates (dynamical Maxwell unusable with \(\rho_b\)) |
| **K-RC1-e** | Neutral locks must source \(\mathbf{E}\) monopole under Supp+F1+MX |

---

## 10. Residuals (honest, non-blocking SPEC)

| ID | Residual |
|----|----------|
| R-RC1-1 | 2D Green multipole shape ≠ 3D \(1/r\) (cosmetic; RC3/M2) |
| R-RC1-2 | Fixed locks: no Cont from motion; \(\mathbf{J}_Q=0\) |
| R-RC1-3 | \(\alpha_\psi\) / \(\gamma\) optical vs \(F^\psi\) (MR2 / J5-β) |
| R-RC1-4 | Poisson isomorphism residual per channel (tags/Occam) |
| R-RC1-5 | Gauss cleanse docs inherited from M1-G5 spirit |

---

## 11. Implementation sketch for NM (non-normative)

```text
1. Import Maxwell2D from work/NE/sandbox_m1_2d.py
2. Place ≥2 locks → ρ_b, ρ_Q on grid (Supp)
3. Solve F1 for ψ (SOR or analytic Green in 2D with same BC honesty)
4. Init E from ρ_Q (Poisson TE_z or NE helper); B=0
5. for n in range(N_step): maxwell.step(rho_Q=ρ_Q, J*=0)
6. Measure F_ψ, F_EM; run RG0–RG6 configs
7. Export rc1_result.json
```

Do **not** edit `scp_sim` / `sfa.h`.

---

## 12. Bottom line

| Item | Status in this SPEC |
|------|---------------------|
| State \((\rho_b,\rho_Q,\psi,\mathbf{E},\mathbf{B})\) | §1 |
| Supp + F1 + **dynamical** Maxwell | §2 |
| Fixed-lock force diagnostics \(F^\psi+Q\mathbf{E}\) | §4 |
| FOR_NM RG0–RG6 | §6 |
| Φ-only RC1 | **Forbidden** for `rc1_claim` |
| CP-M1-NUM | Prerequisite **ADOPTED** (`m1_claim=true`) |

**TM freezes CP-RC1-SPEC content with this file.** Partners stamp ADOPT/DEFER/REJECT in their logs; NM implements toward CP-RC1-NUM.
