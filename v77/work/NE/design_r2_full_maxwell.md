# NE Round 2 — Full Maxwell Yee sandbox design

**Date:** 2026-07-18  
**Agent:** NE  
**O-003:** Full Maxwell progress (not lite only)  
**TE:** `maxwell_monist_v0.md` M1–M4 + Cont + Wave

---

## 1. Scheme

**Name:** `Yee_staggered_2D_TE+TM`

2D Cartesian grid, \(\partial_z=0\), linear isotropic \(\varepsilon,\mu\), \(c=1/\sqrt{\varepsilon\mu}\).

| Sector | Fields | Advances |
|--------|--------|----------|
| **TM_z** | \(E_z, H_x, H_y\) | Faraday (H from E) + Ampère (E from H) |
| **TE_z** | \(E_x, E_y, H_z\) | Faraday (H from E) + Ampère (E from H) |

Linear isotropic media **decouple** TE/TM; running both covers the full 6-component 2D Maxwell content.

**Staggering (Yee):**
- TM: \(E_z\) at nodes \((i,j)\); \(H_x\) at \((i,j+\tfrac12)\); \(H_y\) at \((i+\tfrac12,j)\)
- TE: \(E_x\) at \((i+\tfrac12,j)\); \(E_y\) at \((i,j+\tfrac12)\); \(H_z\) at \((i+\tfrac12,j+\tfrac12)\)

**CFL:**
- Plane-wave (\(y\)-uniform TM → 1D): \(c\Delta t/\Delta x \le 1\); **CFL=1 exact shift** for pure right-going.
- Isotropic 2D pulse: \(c\Delta t \le \mathrm{cfl}\,\min(\Delta x,\Delta y)/\sqrt{2}\).

---

## 2. Equation map (TE monist)

| Law | Numeric realization |
|-----|---------------------|
| **M1** \(\nabla\cdot\mathbf{B}=0\) | TM discrete div B; conserved if initially 0 (Yee) |
| **M2** Faraday | H update from curl E |
| **M3** Gauss | TE static Poisson projection; R1 3D Coulomb recovery |
| **M4** Ampère–Maxwell | E update from curl H − J |
| **Cont** | Prescribed \(\rho(x-vt)\), \(J=v\rho\) |
| **Wave** | Peak track \(v\approx c_{\mathrm{th}}=1/\sqrt{\varepsilon\mu}\) unit + off-unit |

---

## 3. Gates (FM1–FM7)

| Gate | Pass criterion |
|------|----------------|
| FM1 vacuum | max\|E\|, max\|H\| floor |
| FM2 wave unit | \|v/c_th−1\|<5%, c_th=1 |
| FM2 off-unit | same with ε=4,μ=1 → c_th=1/2 |
| FM3 div B | max\|div B\| ~0 (conservation) |
| FM4 Faraday | discrete update residual ~0 |
| FM5 Gauss static | relative Gauss residual small on TE |
| FM6 continuity | \|ΔQ\|/Q small; Cont residual small |
| FM7 Coulomb 3D | R1 KG2∧KG3 recovery |

**full_maxwell_claim** = all FM1–FM7 PASS.

---

## 4. Ontology tags

```text
sector=1
E_origin=free_maxwell
em_solver=yee_2d_full
gravity_solver=none
monist_channel=free_gauge
```

Dualist twin: `sector=2`, same solver, stage ontology / idle free DOF.

---

## 5. What is NOT claimed

- 3D Yee multipole radiation suite  
- Dual-channel joint \(\psi\)+Maxwell (V77-2)  
- Real \(\alpha_{\mathrm{EM}}\)  
- That CFL=1 theorem replaces all discrete grid tests — dynamical sandbox exists for regression

---

## 6. Files

| File | Role |
|------|------|
| `sandbox_full_maxwell_r2.py` | Full dynamical Yee TE+TM gates |
| `offline_r2_full_maxwell.py` | `--quick` runner |
| `outputs/r2_result.json` | Gate export |
| `outputs/r2_wave_track*.tsv` | Peak tracks |
| `outputs/r2_summary.txt` | Human summary |
