# v81 path2_ell — Path-measure \(\ell\) + locks (P2)

**OP:** `v81/OP.md` §2 P2  
**Status:** see `STATUS.md`  
**Scope:** CPU sandbox only. No `scp_sim` / `sfa.h` edits. No multi-fab. No pairwise Coulomb force in the push.

## Thesis

```text
REAL  = FreeSubstrate (Yee E,B; fixed free c)
      + Locks (typed structs, not field humps)
      + scalar path-measure ℓ(x,t)
```

- Locks **source** \(\ell\) (thicken).
- Free-frame locality \(c\) is **fixed** (C-chart engine law).
- \(\ell\) obeys a **hyperbolic** sourced wave equation — never Poisson-per-tick.
- Non-decorative coupling: lock force includes \(-\kappa\nabla\ell\) in addition to \(q(\mathbf{E}+\mathbf{v}\times\mathbf{B})\).
- Optical \(c_{\mathrm{eff}}\sim c/\ell\) is an **M-chart readout only** — not free-frame gravity (v76 GRIN kill).

## Hyperbolic law (implemented)

\[
\ddot\ell - c^2\nabla^2\ell + \gamma\,\dot\ell + m^2(\ell-\ell_0) = S_{\mathrm{lock}}(x,t)
\]

- Leapfrog on \(\ell,\pi=\dot\ell\).
- \(S_{\mathrm{lock}}\) = Gaussian footprint weighted by lock \(E_\star\).
- Mild \(m^2\) → Yukawa IR control (2D); \(\gamma\) → stationary non-runaway.

**Kill if** elliptic Poisson-per-tick required for sanity.  
**Kill if** dynamics unchanged when \(\kappa=0\) / source off (decorative \(\ell\)).  
**Kill if** free-frame variable-\(c\) GRIN gravity claimed (v76).

## State

| Piece | Representation |
|-------|----------------|
| Free medium | 2D collocated TE: `Ex, Ey, Bz` + deposit `ρ, Jx, Jy` |
| Locks | `Lock {id, type, q, m, E_star, x,y, vx,vy, pinned, footprint}` array |
| Path measure | `ell[NX][NY]`, `pi[NX][NY]` |
| Fixed free \(c\) | `C_FREE = 1`, CFL-limited `dt` |

## Step

1. Maxwell update (fixed \(c\)) + charge/current deposit from locks  
2. Hyperbolic \(\ell\) step with lock source  
3. Push locks: \(q(E+v\times B) - \kappa\nabla\ell\) (toggles for kill checks)

## Build & run

```bash
make -C v81/path2_ell/src
./v81/path2_ell/src/path2_ell all v81/path2_ell/out
# or: e0 | e1 | e2
```

Requires: `gcc`, `libm`. OpenMP not required.

## Experiments

| ID | Setup | Pass criteria |
|----|--------|----------------|
| **E0** | One heavy lock pinned; \(\ell\) relaxes | Stationary non-runaway \(\ell(r)\); finite light-travel time for impulse |
| **E1** | Light test flyby | Deflection with \(\kappa\) on; **null** with \(\kappa\) off / source off; GRIN kill holds |
| **E2** | Two-lock \(q=\pm1\) ± \(\ell\) | Improvement vs pure Maxwell **or** honest null; locks intact by construction |

Outputs under `out/`: `e*_summary.txt`, trajectories, radial \(\ell\).

## Relation to P1

Vendored **minimal** medium+lock numerics inside `src/path2_ell.c` (no dependency on `path1_locks/`).  
If P1 greens and kernel port happens later, \(\ell\) can ride as a co-field after P1 medium is stable — not in this MVP.

## Explicit non-claims

- Not 3D Einstein multipole (2D Yukawa/\(\log\)-class exterior).  
- Not durable multi-rev positronium park (E2 is collapse/attract MVP).  
- Not free-frame GRIN gravity.  
- Not multi-fab L / Cosserat light.
