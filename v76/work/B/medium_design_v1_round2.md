# Approach B — Round 2 medium design (dynamical free-response)

**Date:** 2026-07-18  
**Mandate (O-001):** Make free-response monist or kill it.  
**Prior:** Round 1 B2-lite local GRIN (compact only) + postulated kernel Φ=α∫ρ_b/R (Poisson-like).

---

## 0. Verdict target

| Tag | Meaning |
|-----|---------|
| `monist_1sector` | Single continuum; free updates only; bound is freeze of same continuum; path cost from free dynamics; **no** independent Φ solve from ρ as matter |
| `dualist_2sector` | Independent gravity/Φ sector sourced by ρ_bound (Poisson/T→G style) |
| `monist_kernel_failed` | Long-range Einstein-class effect only appears when 2-sector Φ is added |

---

## 1. Three mechanisms (all pure Python sandboxes)

### M1 — Free-graph geodesic monism (`monist_1sector`)

- Nodes carry free budget; locks freeze budget into bound.
- Edge travel cost = function of free budgets on endpoints only.
- Rays = Dijkstra / multi-source free travel times (no Φ field).
- **Expected:** compact path-cost excess (C no-go for 1/r from local free only).

### M2 — Free-medium Laplace relaxation (`monist_1sector` candidate)

- Single free field \(u\) on free cells only.
- Dynamics: Jacobi relaxation of discrete Laplace **among free neighbors**:
  \[
  u_i \leftarrow \mathrm{mean}(u_{\mathrm{free\ nbrs}})\quad\text{or fixed BC}.
  \]
- Bound cells are **holes** (excluded from free update) — same continuum frozen.
- Interface BC: \(u = -\kappa\,\rho_{\mathrm{bound}}\) on free cells adjacent to lock
  (lock depth sets free-envelope stress — still one medium state).
- Far boundary: \(u=0\).
- Chart optical factor \(n = 1 + \beta\,|u|\) (or path cost from free eikonal on \(u\)).
- **Expected in 2D:** exterior multipole ~ **log** (2D Green), not 1/r.
- **Eligibility:** one sector if \(u\) is free continuum equilibrium, not a second
  gravity field. Honest risk: Dirichlet BC + Laplace looks like electrostatics;
  monist claim is “free medium equilibration,” not T→G. D scores N_sec.

### M3 — Dualist Poisson baseline (`dualist_2sector`)

- Solve \(\nabla^2\Phi = \rho_{\mathrm{bound}}\) (Jacobi) on full grid.
- Rays from \(n=1+\beta|\Phi|\).
- **Explicitly dualist** for D3 Occam; not monist success.

### M4 — Inertia triad (smoke)

- Ledger \(m_L = E_{\mathrm{bound}}/c^2\).
- Push: free-medium directed bias (pressure gradient) on lock envelope.
- Measure force proxy and acceleration of lock centroid under overdamped
  dynamics; \(m_{\mathrm{inert}} \sim F/a\).
- Ray-inferred \(M_{\mathrm{ray}}\) from weak Born slope of deflection vs 1/b.
- Report triad \((m_L, m_{\mathrm{inert}}, M_{\mathrm{ray}})\).

---

## 2. Kill criterion (orchestrator)

If **only** M3 produces long-range monopole-class rays while M1–M2 stay compact
(or log-only without Einstein 1/r target), log:

> `monist_kernel_failed` for Einstein-class long-range from local free dynamics.

Keep M3 as dualist baseline for D. Keep M1–M2 as monist-eligible **compact**
or **2D-log** free dynamics — do not relabel M3 as monist.

---

## 3. Export for D

TSV maps under `outputs/round2_*.tsv`:

- `sector_tag`, `mechanism`, `b`, `defl`, `delay`, `m_ledger`, …
- path-cost radial samples
- free deficit

---

## 4. Code

| File | Role |
|------|------|
| `sandbox_r2_dynamical.py` | M1–M4 pure stdlib; writes all outputs |
| `outputs/round2_result.json` | full numeric package |
| `outputs/round2_summary.txt` | human summary |
