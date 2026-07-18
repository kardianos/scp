# Constitutive table — free-capacity vs free-gauge (v0)

**Agent:** TE  
**Date:** 2026-07-18  
**Round:** 1  
**Parent:** [`maxwell_monist_v0.md`](maxwell_monist_v0.md)  
**Seed:** v76 F1-3D constitutive \((\sigma,s,\gamma,c)\) + Maxwell-lite \((\varepsilon,\mu)\)

---

## 1. Side-by-side constitutive map

| Item | **Free capacity** (path-cost / F1-3D) | **Free gauge** (Maxwell-lite) |
|------|--------------------------------------|-------------------------------|
| **Channel name** | Free capacity | Free gauge |
| **Primary DOF** | Scalar potential \(\psi\) | \(\mathbf{E},\mathbf{B}\) (or \(\Phi,\mathbf{A}\)) |
| **Linear vacuum operator** | \(-\sigma_0\nabla^2\) (elliptic quasistatic) | \(\varepsilon\partial_t^2 - \nabla\times(\mu^{-1}\nabla\times)\) (hyperbolic) + Gauss constraint |
| **Quasistatic law** | \(-\nabla\cdot(\sigma\nabla\psi)=s\rho_b\) | \(\nabla\cdot(\varepsilon\mathbf{E})=\rho_Q\), \(-\nabla\cdot(\varepsilon\nabla\Phi)=\rho_Q\) |
| **Constitutive params** | \(\sigma(\rho_f)>0\), \(s\), \(\gamma\) | \(\varepsilon>0\), \(\mu>0\) |
| **Vacuum defaults** | \(\sigma=\sigma_0\) const | \(\varepsilon=\varepsilon_0\), \(\mu=\mu_0\) const |
| **Locality / signal speed** | Free locality \(c\) (path cost holds local free speed) | Wave speed \(c=1/\sqrt{\varepsilon\mu}\) |
| **Shared \(c\) rule** | Same free-medium \(c\) | **Must match** capacity-channel \(c\) in joint world |
| **Source ledger** | \(\rho_b\ge 0\) bound / mass-form | \(\rho_Q\in\mathbb{R}\) gauge charge |
| **Current / flux** | Free capacity flux \(\boldsymbol{\Pi}_\psi=-\sigma\nabla\psi\) | \(\mathbf{J}_Q\); Poynting \(\mathbf{E}\times\mathbf{H}\) |
| **Exterior Green (3D)** | \(\psi\sim s E_\star/(4\pi\sigma_0 r)\) | \(\Phi\sim Q/(4\pi\varepsilon_0 r)\), \(E\sim 1/r^2\) |
| **Observable map** | \(\ell=\ell_0+\gamma\psi\); rays / delay | Forces on \(Q\); radiation; Gauss flux |
| **Emergent coupling** | \(G_{\mathrm{eff}}=\gamma s c^4/(8\pi\sigma_0)\) | Coulomb strength \(\propto 1/\varepsilon\); \(\alpha_{\mathrm{eff}}\) later |
| **Energy density (linear)** | Free-capacity stress (v76 \(U[\psi]\); Dynamics) | \(u_{\mathrm{EM}}=\varepsilon|\mathbf{E}|^2/2 + |\mathbf{B}|^2/(2\mu)\) |
| **Budget role** | Forming \(\rho_b\) depletes free capacity | \(u_{\mathrm{EM}}\) is free-channel stress; \(\rho_Q\) does **not** replace \(\rho_b\) |
| **Sign / force structure** | Path-cost: same-sign mass attraction (weak-field) | \(Q_i Q_j\): attract or repel |
| **Null control** | \(\rho_b=0\Rightarrow\psi=0\Rightarrow\) null deflection | \(\rho_Q=0,\mathbf{J}_Q=0\Rightarrow\) static vacuum \(\mathbf{E}=\mathbf{B}=0\) (gauge) |
| **Nonlinear extension (later)** | \(\sigma=\sigma(\rho_f)\); multi-lock | Media \(\varepsilon(\mathbf{x})\); lock form-factors — **not** R1 required |
| **Numeric sandbox family** | v76 B free-capacity 3D | NE Gauss/Coulomb/wave |
| **Demo IDs** | D-ψ-3D (seed LIVE) | D-EM-* (this focus) |

---

## 2. Joint constitutive constraints (dual-channel)

When both channels live in one medium (V77-2):

| Constraint ID | Statement | Owner if fails |
|---------------|-----------|----------------|
| **JC1** | Single \(c\): path-cost locality \(\equiv 1/\sqrt{\varepsilon\mu}\) in free vacuum | TE + TD |
| **JC2** | Budget: \(\rho_f+\rho_b=\rho_{\mathrm{tot}}\) still holds; free-gauge energy counted in free ledger | TE + TD + TM |
| **JC3** | Sources: \(\rho_b\) does not automatically equal \(\rho_Q\) | TM design; TE forbids collapse |
| **JC4** | No second-sector constants: \(G_{\mathrm{eff}}\) and Coulomb strength both from free constitutive | TU Occam |
| **JC5** | Sector tags: free-origin for \(\psi\) and for \(A/\mathbf{E}\) | NE + TU |

---

## 3. Minimal sandbox unit choices (recommendation to NE)

| Choice | Value | Rationale |
|--------|-------|-----------|
| \(\varepsilon\) | \(1\) | Coulomb \(\Phi=Q/(4\pi r)\) simple |
| \(\mu\) | \(1\) | \(c=1\) |
| \(c\) | \(1\) | Match Dynamics natural units |
| Grid | 3D Cartesian, open/absorbing or large Dirichlet \(\Phi\to 0\) | Exterior multipole |
| Charge model | Compact blob or single-cell \(Q\) with smooth kernel | Avoid pure 1-point singularity if discrete Gauss noisy |
| Diagnostics | Gauss residual, \(\mathcal{R}_E(r)\), \(c_{\mathrm{meas}}\), vacuum floor | Kill-gates KG1–KG4 |

Path-cost constitutive \((\sigma_0,s,\gamma)\) **unchanged** from v76 when joint demos appear; do not retune to fake EM.

---

## 4. What constitutive table forbids

1. Setting \(c_{\mathrm{EM}}\neq c_{\mathrm{locality}}\) without a documented medium-frame story (default: **forbidden** in R1).  
2. Using \(\varepsilon(\rho_b)\) as a **hidden gravity** channel to replace \(\psi\) (re-opens local-GRIN-class dualism under another name).  
3. Absorbing \(\gamma\) into \(\varepsilon\) to force \(\psi\equiv\Phi\).

---

## 5. Bottom line

Free capacity and free gauge share **one continuum and one \(c\)**, but use **independent constitutive tuples** \((\sigma,s,\gamma)\) vs \((\varepsilon,\mu)\) and **independent sources** \(\rho_b\) vs \(\rho_Q\). That split is the monist dual-channel design for V77-2.
