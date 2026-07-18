# RC1 EM interface — dual-source dynamical Maxwell (v0)

**Agent:** TE  
**Date:** 2026-07-19  
**Round:** P2-R3  
**Checkpoint:** **CP-RC1-SPEC**  
**Depends:** CP-M1-NUM green (`m1_claim=true`); `full_maxwell_monist_v0.md` FROZEN  
**API source:** `v77/work/NE/sandbox_m1_2d.py` — class **`Maxwell2D`**  
**Partners:** TM (joint state), NM (co-field lead), NE (module), TD (ψ time-dep note)

---

## 0. One-line claim

**RC1 co-field uses one continuum: fixed multi-locks source \(\psi\) via \(\rho_b\) (F1) and dynamical free-gauge \((\mathbf{E},\mathbf{B})\) via lock charge/current fed to `Maxwell2D.step` — not Φ-only Poisson, not \(\psi\equiv\Phi\).**

---

## 1. Dual-source attachment (theory)

### 1.1 Ledgers on fixed locks \(L_i\) (RC1: \(\mathbf{V}_i=\mathbf{0}\))

| Ledger | Continuum field | Free channel |
|--------|-----------------|--------------|
| Mass-form | \(\rho_b=\sum_i E_{\star,i} f_i(\mathbf{x}-\mathbf{X}_i)\) | C-ψ: \(-\nabla\cdot(\sigma\nabla\psi)=s\rho_b\) |
| Gauge charge | \(\rho_Q=\sum_i Q_i f_i(\mathbf{x}-\mathbf{X}_i)\) | C-A: M3 \(\nabla\cdot(\varepsilon\mathbf{E})=\rho_Q\) |
| Gauge current | \(\mathbf{J}_Q\) (RC1 default \(\mathbf{0}\) if static charge; else Cont-safe) | C-A: M4 Ampère–Maxwell |

**Supp:** \(\mathrm{supp}(|\rho_Q|)\subseteq\mathrm{supp}(\rho_b)\) (TM DS-B lite).  
**TE-IA1:** \(\psi \not\equiv \Phi_{\mathrm{EM}}\); neutral \(Q=0\Rightarrow\) no Coulomb \(\mathbf{E}\) from charge; mass still sources \(\psi\).

### 1.2 RC1 dynamical system (co-field, fixed locks)

\[
\boxed{
\begin{aligned}
\rho_f+\rho_b&=\rho_{\mathrm{tot}},\\
-\nabla\cdot(\sigma\nabla\psi)&=s\,\rho_b
\quad\text{(quasistatic F1; or TD T1 if chosen)},\\
&\text{C-A: Yee advance of }(\mathbf{E},\mathbf{B})\text{ with sources }(\rho_Q,\mathbf{J}_Q),\\
c&=1/\sqrt{\varepsilon\mu}=c_{\mathrm{path}}.
\end{aligned}
}
\tag{RC1}
\]

**Hard rule (campaign):** RC1 **REJECT** if EM half is only Maxwell-lite \(\Phi\) without a dynamical \((\mathbf{E},\mathbf{B})\) path.

### 1.3 Force diagnostics (fixed locks — no motion)

\[
\mathbf{F}_i^{\psi} \text{ from }\nabla\psi,\qquad
\mathbf{F}_i^{\mathrm{EM}} \approx Q_i\,\mathbf{E}(\mathbf{X}_i)
\quad(+\;Q_i\mathbf{v}\times\mathbf{B}\text{ only if }v\neq 0\text{; RC1 }v=0).
\]

Total diagnostic: \(\mathbf{F}_i=\mathbf{F}_i^{\psi}+\mathbf{F}_i^{\mathrm{EM}}\).  
Sign structure: path-cost attract; Coulomb by \(Q_i Q_j\). Hierarchy constitutive (TM Ξ).

---

## 2. NE module API (binding for NM)

**Module:** `work/NE/sandbox_m1_2d.py`  
**Class:** `Maxwell2D`  
**Provenance:** M1 suite green; `api.rc1_ready=true` in `outputs/m1_result.json`.

### 2.1 Construct

```text
Maxwell2D(
  nx, ny,           # grid
  Lx, Ly,           # box size → dx=Lx/(nx-1), dy=Ly/(ny-1)
  eps=1, mu=1,      # c = 1/sqrt(eps*mu)
  cfl=0.5,          # dt = cfl * min(dx,dy) / (c*sqrt(2))
  periodic=False,
  incomplete_ampere=False   # adversary only; RC1 must use False
)
```

Shared-\(c\): use same `eps,mu` (and thus `c`) as path-cost locality / ND `C_LOCAL` (default 1).

### 2.2 Step (drive free gauge)

```text
Maxwell2D.step(
  rho_Q=None,   # 2D array (nx,ny) cell-centered charge density (diagnostics / Gauss)
  Jx=None,      # TE edge current → updates Ex
  Jy=None,      # TE edge current → updates Ey
  Jz=None,      # TM face current → updates Ez
)
```

| Arg | Role in RC1 |
|-----|-------------|
| `Jx, Jy, Jz` | **Dynamical sources** for Ampère–Maxwell (M4). Fixed static locks: all `None` or zeros after IC. |
| `rho_Q` | Charge density for Gauss diagnostics / optional projection (see §2.4). Cont must hold if \(\rho_Q\) time-varies: \(\partial_t\rho_Q+\nabla\cdot\mathbf{J}=0\). |

Advance: one Yee leapfrog substep of **TM** \((E_z,H_x,H_y)\) + **TE** \((E_x,E_y,H_z)\) = full 2D Maxwell.

### 2.3 Read fields

```text
f = m.fields()   # dict: Ex, Ey, Ez, Hx, Hy, Hz  (lists of lists)
# B = μ H in linear constitutive (mu stored on m)
m.clear()        # zero fields; reset t, nstep
m.t, m.dt, m.c, m.eps, m.mu, m.dx, m.dy
```

**NM sampling at lock center:** interpolate `Ex,Ey` (and `Ez` if used) to \(\mathbf{X}_i\) for \(Q_i\mathbf{E}\).

### 2.4 Gauss / \(\rho_Q\) honesty (RC1)

Yee `step` advances via **currents** (and free curls). For static charged locks:

1. **Preferred RC1-A (static EM):** project electrostatic TE fields once from \(\rho_Q\) (Poisson / Gauss cleanse as in M1-G5), then `step()` with \(\mathbf{J}=\mathbf{0}\) for free evolution / radiation tests; **or** hold quasi-static by re-projecting every \(N\) steps.  
2. **RC1-B (dynamic sources):** prescribe Cont-safe \((\rho_Q(t),\mathbf{J}_Q(t))\) and call `step(rho_Q=…, Jx=…, …)` each tick (M1-G5 pattern).

Do **not** claim monist Coulomb from \(\mathbf{J}=0\) with zero IC and never projecting \(\rho_Q\).

### 2.5 Tags (export with co-field JSON)

```text
sector=1
E_origin=free_maxwell_full
em_solver=free_maxwell_yee_m1
gravity_solver=none
dual_channel=1
psi_origin=free_capacity_*     # NM
embedding_dim_dynamics=2
```

---

## 3. Co-field schedule (recommended for NM)

```text
1. Place fixed locks → build ρ_b, ρ_Q on shared (nx,ny) grid [or ψ on coarser grid with documented interp]
2. Solve / relax ψ from ρ_b (F1; TD default quasistatic OK for RC1)
3. Init Maxwell2D(same Lx,Ly,eps,mu,cfl); project E from ρ_Q if static charge
4. Loop:
     - optional: refresh ψ if ρ_b fixed once is enough (skip)
     - m.step(rho_Q=ρ_Q, Jx=…, Jy=…, Jz=…)   # J=0 if static
     - sample E at X_i; F_i = Fψ_i + Q_i E_i
5. Gates: vacuum; neutral (E~0, ψ>0); same/opposite Q; shared c; energy split; sibling
```

**Same grid preferred.** If ψ solver differs in N, document interpolation; do not use different \(c\).

---

## 4. Kill / reject (RC1 EM half)

| ID | REJECT if… |
|----|------------|
| **RC1-EM1** | EM is only \(-\varepsilon\nabla^2\Phi=\rho_Q\) with no `Maxwell2D` / no dynamical \(\mathbf{B}\) path available |
| **RC1-EM2** | \(\psi\) field array is copied into \(\Phi\) or \(\mathbf{E}=-\nabla\psi\) by identity |
| **RC1-EM3** | Neutral locks produce Coulomb-scale \(\mathbf{E}\) from mass-only \(\rho_b\) |
| **RC1-EM4** | `incomplete_ampere=True` used as production physics |
| **RC1-EM5** | Shared \(c\) split: path-cost \(c\) ≠ `m.c` without frame story |

---

## 5. Minimal gate hints for CP-RC1-NUM (EM channel)

| Gate | PASS sketch |
|------|-------------|
| Vacuum | No locks → \(\psi\approx 0\), \(\|\mathbf{E}\|,\|\mathbf{B}\|\) floor |
| Neutral | \(\rho_b>0,\rho_Q=0\) → \(\psi\sim 1/r\) (or 2D Green); \(Q_i E\) ~ 0 |
| Same-sign \(Q\) | \(F^{\mathrm{EM}}\) repel; \(F^\psi\) attract |
| Opposite \(Q\) | \(F^{\mathrm{EM}}\) attract; \(\psi\) still same-sign mass |
| Shared \(c\) | `m.c == C_LOCAL` (unit constitutive) |
| Not Φ-only | Export proves `Maxwell2D.step` called or dynamical fields present (`Hz`/`Hx` path live) |

Full NM gate list: TM/NM ownership; TE stamps EM half only.

---

## 6. What TE freezes / does not own

| Freeze | Owner |
|--------|-------|
| C-A equations M1–M4, TE-IA1, JC1–JC7 | TE (`full_maxwell_monist_v0`, dual_channel_joint) |
| This interface + `Maxwell2D.step` contract | TE + NE API |
| Lock ontology, \(F^\psi+F^{\mathrm{EM}}\), hierarchy | TM |
| Co-field sandbox + numeric gates | **NM lead** |
| ψ time-dep vs quasistatic | TD note (default: quasistatic F1 OK for RC1) |

---

## 7. Bottom line

**CP-RC1-SPEC (TE):** dual-source RC1 = F1 \(\psi(\rho_b)\) + **dynamical** Maxwell via  
`Maxwell2D.step(rho_Q, Jx, Jy, Jz)` / `fields()` from `sandbox_m1_2d.py`.  
Fixed locks; force diagnostics; no channel collapse. NM may implement co-field against this note.
