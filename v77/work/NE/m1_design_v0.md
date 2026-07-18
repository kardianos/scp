# NE M1 design v0 — True 2D dynamic full Maxwell

**Agent:** NE  
**Date:** 2026-07-19  
**Phase:** 2 · Round P2-R1  
**Checkpoint:** **CP-M1-SPEC**  
**Depends on:** M0 = `sandbox_full_maxwell_r2.py` (`full_maxwell_claim=true`);  
`CAMPAIGN_MAP.md` §4.2; TE `full_maxwell_monist_v0.md` (FROZEN);  
TE `for_ne_kill_gates_r2.md` (KG-F1..F6 baseline)  
**TE twin (authoritative gates):** [`../TE/m1_gates_v0.md`](../TE/m1_gates_v0.md)  
**Status:** **ALIGNED with TE** — CP-M1-SPEC co-agree ADOPT; **not** M1 numeric green

---

## 0. Why M1 (beyond M0)

| M0 (done) | M1 (this design) |
|-----------|------------------|
| 2D Yee skeleton TE+TM | Same solver family, **true 2D packets** |
| Wave proof often **1D TEM reduction** (y-uniform) | Wave proof must **not** rely only on 1D TEM |
| Static Gauss / Cont prescribed | **Dynamic** Gauss with Cont-safe sources |
| KG-F6 structural / optional | Explicit **incomplete-Ampère adversary** numeric kill |
| Soft box / Dirichlet | Absorbing **or** large-box honesty + radiation check |
| No energy ledger gate | **Energy / Poynting** vacuum gate |

**Hard rule (CAMPAIGN_MAP / NE CAMPAIGN.md):**  
M1 wave PASS must **not** be claimed solely from 1D TEM; 1D remains **regression only**.

---

## 1. Solver plan — `sandbox_m1_2d.py`

### 1.1 Scheme (inherit + extend M0)

Reuse Yee staggered 2D from `sandbox_full_maxwell_r2.py`:

| Sector | Fields | Role in M1 |
|--------|--------|------------|
| **TM_z** | \(E_z, H_x, H_y\) | 2D cylindrical / Gaussian beam in \(E_z\); div B monitor |
| **TE_z** | \(E_x, E_y, H_z\) | Dynamic Gauss with \(\rho_Q\); in-plane radiation |

**Tags (carry forward):**

```text
sector=1
E_origin=free_maxwell_full
em_solver=free_maxwell_yee_m1
gravity_solver=none
embedding_dim_dynamics=2
m1_wave_is_true_2d=true   # kill if only 1D TEM used for M1-G2
```

**Constitutive:** \(\varepsilon,\mu\) const; \(c=1/\sqrt{\varepsilon\mu}\); unit + one off-unit (\(\varepsilon=4,\mu=1\Rightarrow c=\tfrac12\)).

### 1.2 Boundary options

| Mode | Use |
|------|-----|
| **PEC box** (default smoke) | Simple; document reflection contamination window |
| **Large box** | Packet free-flight before wall hit; primary wave/energy window |
| **PML / polynomial absorber** (preferred M1) | Radiation honesty; optional if pure-Python cost high — large-box + early stop is acceptable **if labeled** |

### 1.3 Module layout (implementation outline)

```text
sandbox_m1_2d.py
  ├── Yee2D state (TE + TM grids, dx, dy, dt, eps, mu)
  ├── step_tm / step_te  (from M0, keep byte-stable physics)
  ├── diagnostics:
  │     divB, divE_gauss, energy_U, poynting_flux, peak_track_2d
  ├── BC: pec | large_box_window | soft_absorb (simple σ layer)
  ├── gates M1-G0 .. M1-G9
  ├── adversary: step_tm_incomplete_ampere (drop ∂t D or set ε→∞ fake)
  └── export outputs/m1_result.json + tsv tracks

API foreshadow (CP-RC1-SPEC later, not required this round):
  class Maxwell2D:
      def step(self, rho_Q=None, J_Q=None) -> None
      @property
      def E, B  # or Ex,Ey,Ez, Hx,Hy,Hz
```

**Regression:** call / subprocess `sandbox_full_maxwell_r2.py` or import M0 gates as M1-G0.

---

## 2. Gate list (M1) — **ADOPTS TE `m1_gates_v0.md` as authoritative**

NE implements **exactly** TE IDs and thresholds. Short implementer checklist
(full PASS/FAIL tables → TE file):

| TE ID | Name | NE implement note |
|-------|------|-------------------|
| **M1-R0** | M0 regression | Re-run/load `sandbox_full_maxwell_r2` / r2_result |
| **M1-G1** | Vacuum | max\|E\|,\|B\| ≤ 1e-12; ≥50 steps |
| **M1-G2** | **True 2D beam** \(v\approx c\) | 2D Gaussian both \(x,y\); coupled E,B; \|v/c−1\|≤**0.08**; `reduction=true_2d_packet`; **not** 1D TEM |
| **M1-G3** | Off-unit \(c\) | ε=4,μ=1 → c_th=0.5; same 8% bar |
| **M1-G4** | Energy / Poynting | periodic: \|ΔU\|/U0≤**0.02** or flux balance 5%; sub-volume if open BC |
| **M1-G5** | Dynamic Gauss+Cont | ρ(x−vt), J=vρ; \|ΔQ\|/Q≤**1e-3**; Gauss rel init≤0.05, non-growing |
| **M1-G6** | div B | norm ≤ **1e-6** double; no secular growth |
| **M1-G7** | Faraday | Yee identity floor or continuum loop ≤5% |
| **M1-G8** | Incomplete-Ampère adversary | full G2-class PASS; adversary fails v/c or freezes; tags |
| **M1-G9** | BC honesty | document `bc` per gate + top-level `bc_honesty` |
| M1-G10 | Radiation flux (opt) | recommended |
| M1-G11 | Shared-\(c\) note (opt) | for ND |

**`m1_claim = true` iff**  
M1-R0 ∧ G1 ∧ G2 ∧ G3 ∧ G4 ∧ G5 ∧ G6 ∧ G7 ∧ G8 ∧ G9  
(TE §1 — NE adopts without modification.)

---

## 3. Pass bar for **CP-M1-NUM** (future)

```text
m1_claim = M1-R0 ∧ G1 ∧ G2 ∧ G3 ∧ G4 ∧ G5 ∧ G6 ∧ G7 ∧ G8 ∧ G9
```

(G10–G11 recommended.)  
**Unlocks:** CP-RC1-SPEC, M2, M3 (per CAMPAIGN_MAP).

**Export:** `work/NE/outputs/m1_result.json` per TE §5 schema  
(`gates.M1_G2_beam2d`, `bc_honesty`, `shared_c`, `m1_claim`, …).

---

## 4. Crosswalk TE R2 / M0 → M1

| Prior | M1 |
|-------|-----|
| KG-F1 / FM3 divB | M1-G6 |
| KG-F2 / FM4 Faraday | M1-G7 |
| KG-F3 wave (1D TEM ok for M0) | M1-G2 **true 2D only** + G3 |
| KG-F4 Cont | M1-G5 dynamic |
| KG-F5 Coulomb | M1-R0 / optional G10 path |
| KG-F6 displacement | M1-G8 adversary (numeric twin) |
| (new) energy | M1-G4 |
| (new) BC honesty | M1-G9 |

---

## 5. Alignment stamp vs TE

| Item | Status |
|------|--------|
| TE `m1_gates_v0.md` present | **YES** (P2-R1) |
| Gate IDs R0, G1–G9 match | **YES** |
| True-2D hard rule | **YES** (shared) |
| Ampère adversary | **YES** (shared) |
| Energy + BC honesty | **YES** (shared) |
| Thresholds | **NE adopts TE numbers** (no DEFER) |
| Conflicts | **none** |

**Non-negotiable (map + TE):** true 2D G2; G8 adversary; G4 energy; monist tags; no scp_sim.

---

## 6. RC1 foreshadow (not this checkpoint)

After M1 green, freeze API for NM:

```text
step(rho_Q, J_Q) → update E,B
read E, B on same grid as ψ
```

NM **must** call this for RC1 (CAMPAIGN: REJECT if NM reimplements Φ-only as “dynamical Maxwell”).

---

## 7. Files / next

| Path | Role |
|------|------|
| `m1_design_v0.md` | This SPEC |
| `sandbox_m1_2d.py` | Implement at CP-M1-NUM (skeleton outline OK at SPEC) |
| `outputs/m1_result.json` | CP-M1-NUM export |
| `sandbox_full_maxwell_r2.py` | M0 regression (M1-G0) |

**P2-R1 done when:** design written + CP-M1-SPEC stamp with TE.  
**P2-R2:** implement gates → m1_result.json.
