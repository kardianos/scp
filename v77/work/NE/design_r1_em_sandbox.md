# NE Round 1 — Monist free-gauge EM sandbox design

**Date:** 2026-07-18  
**Agent:** NE  
**Status:** TE-aligned Maxwell-lite (KG1–4); `full_maxwell_claim=false`  
**Seed:** v76 F1-3D free-capacity monism; FOCI §2 Focus EM  
**TE twin:** `work/TE/maxwell_monist_v0.md`, `for_ne_kill_gates_r1.md`

---

## 0. Claim (what this demo is / is not)

**Is:** A killable numeric package for a **free-gauge channel** sibling to path-cost \(\psi\):

1. Quasistatic Coulomb from a compact **charge lock** with discrete Gauss check.
2. Linear free-gauge wave propagating at shared free locality \(c\).
3. Explicit monist sector tags + optional dualist twin for Occam honesty.

**Is not:** Full Maxwell theory; claim that Nature is monist; claim that Poisson shape alone proves monism; bridge to `scp_sim` kernel.

---

## 1. Shared monist primitives (from v76 seed + FOCI)

| Primitive | EM-channel use |
|-----------|----------------|
| One continuum | Free gauge lives in the same monist medium language as free capacity |
| Free vs bound | Charge lock = bound ledger coupling to free gauge |
| \(c\) free locality | \(c = 1/\sqrt{\varepsilon\mu}\) **and** path-cost free-signal bound — **same number in frame** |
| Energy ledger | EM free energy + lock charge ledger (bookkeeping only this round) |
| Sibling channels | \(\psi\) (path cost, F1-3D) ≠ free gauge \((\Phi,\mathbf{E})\) — different constitutive laws |

**Provisional TE interface (NE assumes until TE supersedes):**

\[
c \;=\; \frac{1}{\sqrt{\varepsilon_0\mu_0}}
\quad\text{(free-gauge constitutive locality)}
\tag{c-shared}
\]

\[
-\nabla\cdot\bigl(\varepsilon_0\nabla\Phi\bigr)
\;=\;
\rho_Q
\quad\text{(quasistatic free-gauge lock coupling)}
\tag{FG-Q}
\]

\[
\mathbf{E} = -\nabla\Phi,
\qquad
\oint_{\partial V}\mathbf{E}\cdot d\mathbf{A}
= \frac{Q_{\mathrm{encl}}}{\varepsilon_0}
\tag{Gauss}
\]

\[
\partial_t^2 A = c^2\,\partial_x^2 A
\quad\text{(1D free-gauge wave; Maxwell-lite)}
\tag{FG-W}
\]

**Charge lock:** compact \(\rho_Q\) (Gaussian in R1) tagged as lock–gauge coupling, not foreign stage charge. Budget identity with path-cost \(\rho_b\) is **not** enforced numerically this round (Matter/V77-2 joint demo).

---

## 2. Demo A — D-EM-gauss-coulomb (3D)  [alias: D-EM-coulomb]

### Setup

- Box \([-L/2,L/2]^3\), \(N^3\) grid, \(\Phi=0\) Dirichlet on boundary.
- Charge lock: \(\rho_Q(\mathbf{x}) = A\exp(-r^2/(2\sigma^2))\), \(Q=\int\rho_Q\,dV\).
- SOR solve of \(-\nabla^2\Phi = \rho_Q/\varepsilon_0\).
- \(\mathbf{E}=-\nabla\Phi\) by central differences.

### Kill gates

| Gate | Pass criterion (R1 defaults) |
|------|------------------------------|
| **G-Gauss** | On a spherical shell outside lock: \(\lvert 4\pi r^2\langle E_r\rangle - Q/\varepsilon_0\rvert / (|Q|/\varepsilon_0) < 0.08\) for mid-range \(r\) |
| **G-1/r** | Exterior \(\Phi\) fit prefer \(1/r\) over log/\(1/r^2\) (\(R^2_{1/r}\) highest among \(\{1/r,\log,1/r^2\}\)) |
| **G-1/r²** | \(\lvert E_r\rvert\) multipole prefer \(1/r^2\) |
| **G-vacuum** | \(\rho_Q\equiv 0\Rightarrow \max\lvert\Phi\rvert,\max\lvert E\rvert\) at floor |

### Tags

```text
demo_id: D-EM-gauss-coulomb
sector: 1
sector_tag: monist_free_gauge_channel
E_origin: free_maxwell_lite
em_solver: free_maxwell_lite
channel: free_gauge_quasistatic
phi_origin: free_gauge_poisson_3d
embedding_dim: 3
c_shared: true
gravity_solver: none
full_maxwell_claim: false
te_gates: KG1, KG2, KG3
```

### Dualist twin (optional Occam)

Same discrete PDE with:

```text
sector_tag: dualist_2sector_poisson
phi_origin: dualist_stage_charge
note: multipole-isomorphic; monism not proven by fit alone
```

---

## 3. Demo B — D-EM-wave-c (1D)

### Setup

- 1D line length \(L_w\), \(N_x\) points, periodic or absorbing soft edges.
- Initial Gaussian pulse on \(A(x,0)\), \(\partial_t A = 0\) (split left/right) **or** pure right-going.
- Leapfrog / CTCS: \(A^{n+1}_i = 2A^n_i - A^{n-1}_i + (c\Delta t/\Delta x)^2(A^n_{i+1}-2A^n_i+A^n_{i-1})\).
- Courant \(c\Delta t/\Delta x \le 1\).

### Kill gates

| Gate | Pass criterion |
|------|----------------|
| **G-v=c** | Peak-track or cross-correlation travel speed \(\lvert v_{\mathrm{meas}}/c - 1\rvert < 0.03\) |
| **G-CFL** | Stable under CFL; energy not exploding |
| **G-vac-wave** | Zero IC → remains zero |

### Tags

```text
demo_id: D-EM-wave-c
sector_tag: monist_free_gauge_channel
channel: free_gauge_wave
c_shared: true
c_def: 1/sqrt(eps0*mu0)
full_maxwell_claim: false
```

---

## 4. Shared \(c\) language (path-cost link)

| Symbol | Path-cost (v76) | Free gauge (this sandbox) |
|--------|-----------------|---------------------------|
| \(c\) | Free-field locality / interaction-rate bound | Same numerical constant; \(c=1/\sqrt{\varepsilon\mu}\) |
| Free DOF | \(\psi\) capacity potential | \(\Phi\) (static) / \(A\) (wave) |
| Bound source | \(\rho_b\) mass-form | \(\rho_Q\) charge lock (R1 independent of \(\rho_b\)) |
| Path cost | \(\ell=\ell_0+\gamma\psi\) | Not used in EM channel R1 |
| Long-range | Free response \(\psi\sim 1/r\) | Free-gauge Coulomb \(\Phi\sim 1/r\) |

**V77-2 target (later):** dual-channel medium with shared \(c\) and budget identity; R1 only shows **independent** free-gauge kill-gates + \(c\) identity in constitutive constants.

---

## 5. Residuals / kill conditions for the *program*

1. **Poisson isomorphism** — fit cannot separate monist free-gauge from dualist electrostatics; needs ontology + TE.
2. If TE requires free-gauge dynamics **incompatible** with shared free locality \(c\), tag FOR_TU kill candidate.
3. If charge lock **must** be a second substance with no free–bound ledger path, monist reduction fails (Matter track).
4. Do **not** promote D-EM-* to full Maxwell until TE writes Faraday + Ampère + continuity in monist language and NE implements matching gates.

---

## 6. Exports

| File | Content |
|------|---------|
| `outputs/r1_result.json` | Gates, fits, tags, parameters |
| `outputs/r1_coulomb_radial.tsv` | \(r,\Phi,E_r\), analytic refs |
| `outputs/r1_gauss_shells.tsv` | shell \(r\), flux, \(Q_{\mathrm{encl}}/\varepsilon_0\), residual |
| `outputs/r1_wave_track.tsv` | \(t\), peak_x, \(v\) estimate |
| `outputs/r1_summary.txt` | human one-pager |

---

## 7. FOR_TE / FOR_TU

**FOR_TE:** Confirm or replace FG-Q / FG-W / (c-shared); specify \(\varepsilon(\rho_f),\mu\); charge ledger vs \(\rho_b\); whether \(\mathbf{B}\) is required before V77-2.

**FOR_TU:** Promote D-EM-gauss-coulomb, D-EM-vacuum-control, D-EM-wave-c to LIVE_PASS; dualist DOCUMENTED; full_maxwell_claim=false; no V77-2 claim (sibling-ψ deferred).
