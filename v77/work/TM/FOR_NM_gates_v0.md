# FOR_NM — Numeric Gates (Matter Round 1)

**From:** TM  
**To:** NM  
**Date:** 2026-07-18  
**Theory:** `lock_ontology_v0.md`  
**Constraint:** sandboxes under `v77/work/NM/` only; no `scp_sim` / `sfa.h`; may adapt v76 free-capacity solvers.

---

## 0. Priority for Round 1

| Priority | Demo ID | Goal |
|----------|---------|------|
| **P0** | `D-MAT-lock-S0` | ≥1 compact lock marker on free medium |
| **P0** | `D-MAT-dual0` | Dual static sources: \(\psi\) from \(\rho_b\), \(\Phi\) from \(\rho_q\) |
| **P1** | `D-MAT-force-tax` | Multi-lock force sign structure (path-cost vs Coulomb) |
| **P2** | `D-MAT-hier` | Constitutive hierarchy ratio \(F^C/F^\psi \gg 1\) for charged pair |

If TE/NE Coulomb solve is not ready, implement \(\Phi\) as linear Poisson twin **on the same grid and budget tags** with `sector_tag=free_gauge_lite`, not as foreign gravity. Path-cost \(\psi\) must remain free-capacity / F1-style (`phi_origin=free_capacity` or equivalent).

---

## 1. Shared medium setup (recommended)

Minimal 3D (or 2D provisional) grid:

```text
state:
  rho_b  >= 0     # bound unlock density (mass ledger)
  rho_q  signed   # charge density; supp(|rho_q|) ⊆ supp(rho_b)
  rho_f  >= 0     # free budget; rho_f + rho_b = rho_tot (strong form OK)
  psi             # free capacity potential
  Phi             # free-gauge electrostatic potential (lite)
  sigma0, s, gamma, eps0, c   # constitutive
```

**Solves (static):**

\[
-\sigma_0\nabla^2\psi = s\,\rho_b,
\qquad
-\varepsilon_0\nabla^2\Phi = \rho_q
\quad(\mathbf{E}=-\nabla\Phi).
\]

**BC:** Dirichlet or mixed on large box; report multipole diagnostics outside cores.

**Locks:** place \(N_L\ge 2\) Gaussian (or compact) cores; parameters \((E_{\star,i}, q_i, \mathbf{X}_i, R_{\mathrm{core}})\).

---

## 2. Gate G0 — Lock S0 (`D-MAT-lock-S0`)

**Setup:** one or more compact \(\rho_b\) cores; \(\rho_q\) optional.

| Check | Pass criterion | Fail |
|-------|----------------|------|
| G0.1 Compact mass | \(E_\star=\int\rho_b > 0\); \(R_{\mathrm{rms}} < L_{\mathrm{box}}/4\) | Diffuse fill |
| G0.2 Free deficit | Core free density below vacuum floor (or integral free capacity drop \(\sim E_\star\)) | No free–bound link |
| G0.3 Not hand mass | Mass used in reports is \(E_\star/c^2\), not independent | Dualist mass primitive |

**Export tags:** `demo=D-MAT-lock-S0`, `lock_class=S0`, `budget=strong|integral`.

---

## 3. Gate G1 — Dual source statics (`D-MAT-dual0`)

**Setup:** ≥2 locks; mix of neutral and charged (recommended: one neutral \(q=0\), one charged \(q\neq 0\), or two opposite charges).

| Check | Pass criterion | Fail |
|-------|----------------|------|
| G1.1 Support constraint | \(\max|\rho_q|\) only where \(\rho_b>\epsilon\) (Supp) | Free-floating charge |
| G1.2 Mass multipole | Far-field \(\psi\) prefers \(1/r\) (3D) or documented 2D log; vacuum \(\rho_b=0\Rightarrow\psi=0\) | Wrong multipole / nonzero vacuum |
| G1.3 Charge multipole | Far-field \(\Phi\) prefers \(1/r\); total \(Q=\int\rho_q\); vacuum \(\rho_q=0\Rightarrow\Phi=0\) | Fail Gauss / multipole |
| G1.4 Sibling independence | Scaling \(E_\star\) at fixed \(q\) changes \(\psi\) amplitude, not \(Q\); scaling \(q\) at fixed \(E_\star\) changes \(\Phi\), not \(E_\star\) | Forced \(\psi\propto\Phi\) always |
| G1.5 Shared grid | Same mesh, same \(c\) parameter, both solves free-sector tagged | Two unrelated foreign solvers with idle medium |

**Optional G1.6:** free conductivity \(\sigma=\sigma(\rho_f)\) still yields exterior \(1/r\) (weak nonlinear).

**Export:** multipole \(R^2\) or fit residual; \(E_{\star,i}\), \(q_i\); `phi_origin`, `gauge_origin`.

---

## 4. Gate G2 — Force taxonomy (`D-MAT-force-tax`)

**Method v0 (virtual work):** exclude self-fields; evaluate

\[
\mathbf{F}_i^{\psi} = -\alpha_\psi E_{\star,i}\nabla\psi_{-i}(\mathbf{X}_i),
\qquad
\mathbf{F}_i^{C} = q_i\mathbf{E}_{-i}(\mathbf{X}_i)
\]

with \(\alpha_\psi\) chosen so two neutral masses recover \(G_{\mathrm{eff}}M_iM_j/r^2\) within tolerance **or** report raw \(-\nabla\psi\) signs only if coefficient open (tag `force_closure=sign_only`).

| Config | Expect |
|--------|--------|
| G2.1 Neutral–neutral \(q=0\) | \(F^C\approx 0\); \(F^\psi\) attractive (inward) |
| G2.2 Like charge, equal \(M\) | \(F^C\) repulsive; \(F^\psi\) attractive; report which wins at chosen params |
| G2.3 Opposite charge | \(F^C\) attractive; \(F^\psi\) attractive |
| G2.4 Control vacuum | No locks → forces 0 |

**Pass:** all four configs match sign table; magnitudes finite; no NaNs.  
**Soft fail OK Round 1:** coefficient of \(F^\psi\) not calibrated to \(G_{\mathrm{eff}}\) if signs + multipoles pass.

---

## 5. Gate G3 — Hierarchy (`D-MAT-hier`)

**Setup:** two elementary charged locks with \(\lambda_q=|q|/E_\star=\lambda_0\); choose constitutive \(\Xi=8\pi\sigma_0 k_C/(\gamma s)\) with \(k_C=1/(4\pi\varepsilon_0)\).

| Check | Pass |
|-------|------|
| G3.1 | Measured \(F^C/F^\psi \approx \Xi\lambda_0^2\) within 20% (or theory formula at same \(r\)) |
| G3.2 | Parameter choice with \(\Xi\lambda_0^2 > 10^2\) (EM-dominant pair) **and** a neutral pair with \(F^C\approx 0\), \(F^\psi>0\) |
| G3.3 | Document that hierarchy used constitutive constants, not channel deletion |

**Fail:** only achieve EM dominance by setting \(s=0\) or \(\gamma=0\) (kills path-cost channel).

---

## 6. Explicit non-goals (do not block Round 1)

- Orbital integration S3 (nice-to-have).
- Full Maxwell waves (NE may own waves; Matter needs static Gauss at minimum).
- Matching real \(G\) or \(\alpha\).
- Loading `.sfa` Q-balls as monism proof.
- scp_sim coupling.

---

## 7. Suggested file layout for NM

```text
v77/work/NM/
  README.md
  medium_dual_v0.py      # or split modules
  outputs/
    dual0_result.json
    force_tax.tsv
    hierarchy.json
```

Inherit patterns from `v76/work/B/sandbox_r3_3d_free.py` for \(\psi\); add second Poisson for \(\Phi\).

---

## 8. Report back tags (for TU / TM)

Please log:

```text
FOR_TM: G0=PASS|FAIL|SKIP ...
FOR_TM: G1=...
FOR_TU: demo D-MAT-* status
```

If dual Poisson is blocked on dimensional grounds (2D log), prefer **3D small-N** over wrong multipole claims.
