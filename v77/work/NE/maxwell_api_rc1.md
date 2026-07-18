# Maxwell2D API freeze for RC1 co-field (CP-RC1-SPEC)

**Agent:** NE  
**Date:** 2026-07-19  
**Round:** P2-R3  
**Checkpoint:** **CP-RC1-SPEC** (API freeze)  
**Depends on:** CP-M1-NUM green (`m1_claim=true`, parent-verified)  
**Implementation:** [`sandbox_m1_2d.py`](sandbox_m1_2d.py) class `Maxwell2D`  
**Theory:** TE `full_maxwell_monist_v0.md` (FROZEN); dual-channel TE-IA1 \(\psi\neq\Phi\)  
**Consumers:** **NM** (lead RC1 co-field); TE/TM interface stamps

---

## 0. One-line freeze

**NM must evolve free-gauge fields by calling `Maxwell2D.step(...)` on the shared grid — not by reimplementing quasistatic \(\Phi\)-only Poisson as “dynamical Maxwell.”**

Reject RC1 if EM channel is only `−ε∇²Φ=ρ_Q` without Yee `step`.

---

## 1. Import path

```python
import sys
sys.path.insert(0, "/home/d/code/scp/v77/work/NE")  # or relative
from sandbox_m1_2d import Maxwell2D, C_LOCAL, EPS0, MU0, TAGS
```

**Module:** `v77/work/NE/sandbox_m1_2d.py`  
**Class:** `Maxwell2D`

---

## 2. Construction

```python
m = Maxwell2D(
    nx=64,          # grid points in x (same as ψ grid recommended)
    ny=64,          # grid points in y
    Lx=16.0,        # box length x
    Ly=16.0,        # box length y
    eps=1.0,        # permittivity (ε)
    mu=1.0,         # permeability (μ)
    cfl=0.5,        # Courant fraction of 2D limit 1/√2
    periodic=False, # True → wrap; False → PEC-ish (no update on boundary E)
    incomplete_ampere=False,  # True only for adversary tests — NOT for RC1 production
)
```

| Attribute | Meaning |
|-----------|---------|
| `m.dx`, `m.dy`, `m.dt` | Spacing and timestep (set in `__init__`) |
| `m.c` | \(c=1/\sqrt{\varepsilon\mu}\) — **must equal** path-cost `C_LOCAL` when constitutive shared |
| `m.eps`, `m.mu` | Constitutive |
| `m.t`, `m.nstep` | Simulation time / step count |
| `m.periodic` | BC mode |

**Shared-\(c\) rule (JC1):** For RC1 joint medium, set `eps`, `mu` so `m.c == C_LOCAL` (default both 1.0). Off-unit media allowed only if ψ channel uses the same free locality story.

---

## 3. Advance one Maxwell step

```python
m.step(
    rho_Q=None,   # List[List[float]] density on nodes; reserved / diagnostics
    Jx=None,      # current J_x on TE staggering (drives Ex Ampère)
    Jy=None,      # current J_y
    Jz=None,      # current J_z (drives Ez TM Ampère)
)
```

### Semantics

1. Advances **TM** \((E_z,H_x,H_y)\) and **TE** \((E_x,E_y,H_z)\) by one Yee leapfrog step.  
2. Implements free Maxwell M2 (Faraday) + M4 (Ampère–Maxwell) with optional \(\mathbf{J}\).  
3. `rho_Q` is **not** used inside the update (caller owns Cont: \(\partial_t\rho_Q+\nabla\cdot\mathbf{J}=0\)).  
   Provide Cont-safe \((\rho_Q,\mathbf{J})\) pairs; Gauss is diagnostic / projection on NM side if needed.  
4. `None` currents ⇒ vacuum free step.  
5. Arrays: `nx × ny` nested lists (`List[List[float]]`), same shape as `m.Ex`, etc.  
6. Mutates fields **in place**; increments `m.t` by `m.dt`.

### Cont-safe source recipe (NM locks)

For a prescribed charge density of lock motion \(\mathbf{v}\) with \(|\mathbf{v}|<c\):

\[
\rho_Q(\mathbf{x},t)=\rho_0(\mathbf{x}-\mathbf{v}t),\qquad
\mathbf{J}_Q=\rho_Q\,\mathbf{v}.
\]

For **fixed** multi-locks (RC1 default): \(\mathbf{v}=\mathbf{0}\), \(\mathbf{J}=\mathbf{0}\); hold \(\rho_Q\) on lock support; optionally re-project \(\mathbf{E}\) quasistatically between steps **or** let free Maxwell radiate (document choice).

---

## 4. Read fields

```python
f = m.fields()
# f["Ex"], f["Ey"], f["Ez"], f["Hx"], f["Hy"], f["Hz"]
# each List[List[float]], shape (nx, ny)
```

| Key | Sector | Role |
|-----|--------|------|
| `Ex`, `Ey`, `Hz` | TE\(_z\) | In-plane E; useful for planar Coulomb / lock forces \(q\mathbf{E}\) |
| `Ez`, `Hx`, `Hy` | TM\(_z\) | Out-of-plane E; radiation / waves |

**Lorentz / force on lock \(i\) (NM Tier B when available):**

\[
\mathbf{F}_i^{\mathrm{EM}} \approx q_i\,\mathbf{E}(\mathbf{x}_i)
\quad\text{(RC1 fixed locks: diagnostics; \(\mathbf{v}\times\mathbf{B}\) deferred to RC2)}
\]

Sample \(\mathbf{E}\) by bilinear interpolate of `Ex`,`Ey` (and `Ez` if needed) at lock centroid.

**Energy (TM sector helper):** `m.energy_tm()` → \(U_{\mathrm{TM}}=\int(\varepsilon E_z^2 + \mu(H_x^2+H_y^2))/2\,dA\).

---

## 5. Grid alignment with ψ (NM co-field)

| Item | Requirement |
|------|-------------|
| Same `nx, ny, Lx, Ly` | Preferred: identical Cartesian box as F1 \(\psi\) solve |
| Same origin | Node \((i,j)\) at \((-L/2+i\,dx,\ldots)\) **or** document NM convention `(i*dx, j*dy)` from corner — **NE uses** \(x=i\,dx\in[0,L_x]\), \(y=j\,dy\in[0,L_y]\) |
| \(\rho_b\) vs \(\rho_Q\) | Independent ledgers; `Supp(\rho_Q)\subseteq Supp(\rho_b)` per TM DUAL-0 when claiming monist dual-source |
| Solver order | Each co-field step: update \(\psi\) (F1/T1) **and** `m.step(...)` once (or documented subcycle); **not** Φ-only instead of `step` |

**TE-IA1:** Never set \(\psi\equiv\Phi\) or identify \(\psi\) with any of `Ex,Ey,Ez`.

---

## 6. Tags (export on RC1 JSON)

```text
sector=1
E_origin=free_maxwell_full
em_solver=free_maxwell_yee_m1
gravity_solver=none
embedding_dim_dynamics=2
dual_channel=1              # when ψ co-present
c_shared=true
C_LOCAL = m.c
```

Forbidden for RC1 production: `incomplete_ampere=True`, `free_gauge_lite` as sole EM.

---

## 7. Minimal NM co-field loop (sketch)

```python
from sandbox_m1_2d import Maxwell2D

# ψ solver is NM-owned (F1-3D / SOR); same nx, L
m = Maxwell2D(nx=N, ny=N, Lx=L, Ly=L, eps=1.0, mu=1.0, cfl=0.5)

# fixed locks: rho_Q from lock charges, J=0
for n in range(nsteps):
    # 1) assemble rho_b, rho_Q on grid from locks
    # 2) solve / relax ψ from rho_b  (NM)
    # 3) dynamical Maxwell
    m.step(rho_Q=rho_Q, Jx=Jx, Jy=Jy, Jz=None)
    # 4) forces: F_psi from ∇ψ; F_EM from q E (sample m.fields())
    # 5) locks fixed in RC1 — diagnostics only
```

---

## 8. Kill / reject conditions (API contract)

| ID | Reject RC1 if… |
|----|----------------|
| **K-API-1** | EM advanced only by Poisson \(\Phi\) with no `Maxwell2D.step` |
| **K-API-2** | Different \(c\) for \(\psi\) path and \(1/\sqrt{\varepsilon\mu}\) without frame story |
| **K-API-3** | \(\psi\) identified with electrostatic potential |
| **K-API-4** | Production run uses `incomplete_ampere=True` |

---

## 9. Freeze stamp

| Item | Status |
|------|--------|
| Class name | **`Maxwell2D`** frozen |
| `step` signature | **`step(rho_Q=None, Jx=None, Jy=None, Jz=None)`** frozen |
| `fields()` keys | **`Ex,Ey,Ez,Hx,Hy,Hz`** frozen |
| Grid list type | `List[List[float]]` nx×ny |
| File | `sandbox_m1_2d.py` |

Amend only via NE+NM+TE co-stamp (new version `maxwell_api_rc1_v1.md` if breaking).

---

## 10. Relation to M1 green

Parent-verified `outputs/m1_result.json`:

| Metric | Value |
|--------|------:|
| `m1_claim` | **true** |
| G2 \(v/c\) (2D) | ≈ 1.044 |
| G3 off-unit \(c_{\mathrm{th}}\) | 0.5, \(v/c\) ≈ 1.044 |
| G4 energy drift | ≈ 0.0049 |
| G8 adversary \(v\) | 0.0 (fails free wave) |

M1 LIVE unlocks this API for **CP-RC1-NUM** (NM implements co-field using this freeze).
