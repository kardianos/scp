# Approach B Round 3 — 3D free-capacity medium (F1)

**Mandate (O-003):** 3D free Green ≃ 1/r monist dynamical free response.  
**Prior:** R2 M2 free Laplace = monist but **2D log**; R1 hand 1/r killed as monist-derived.

---

## Dynamics (A F1)

\[
-\sigma_0\nabla^2\psi = s\,\rho_b,\qquad
\ell=\ell_0+\gamma\psi,\qquad
\rho_f+\rho_b=\rho_0.
\]

| Symbol | Role |
|--------|------|
| \(\psi\) | Free continuum capacity potential (single sector DOF) |
| \(\rho_b\) | Bound ledger of **same** continuum |
| \(\sigma_0,s,\gamma\) | Free constitutive constants |
| \(c\) | Free locality (rays use local free speed in frame) |

Infinite-space exterior monopole:

\[
\psi\sim\frac{s}{\sigma_0}\frac{M}{4\pi r},\quad M=\int\rho_b\,dV/c^2.
\]

Path-cost coupling:

\[
n-1=\gamma\psi \implies \alpha_{\mathrm{eff}}=\frac{\gamma s}{4\pi\sigma_0}
\quad(n-1=\alpha_{\mathrm{eff}} M/r).
\]

**Not monist if:** code hard-sets \(\psi=\alpha\int\rho/R\) without free law (F5), or treats \(\psi\) as independent gravity with foreign \(G\).

**Dualist control:** same field tagged `dualist_2sector` for D Occam (math isomorphic in linear vacuum).

---

## Codes

| File | Role |
|------|------|
| `sandbox_r3_3d_free.py` | N³ SOR + Born + exports |
| `offline_round3_3d.py` | Analytic Green + mini SOR |
| `outputs/round3_*` | Maps for D |

---

## Success

- Exterior fit prefers **1/r** over log (R² high).
- Free deficit > 0 where bound rises.
- Ray defl/delay nonzero; vacuum ≈ 0.
- `gravity_solver=None`, `phi_origin=free_capacity_3d_green` / free_relaxation.
- `sector_tag=monist_1sector`.
