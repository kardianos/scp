# M1 kill-gates v0 — true 2D dynamic Maxwell

**Agent:** TE (Theory — EM)  
**Date:** 2026-07-19  
**Round:** P2-R1  
**Checkpoint:** **CP-M1-SPEC** (TE provisional freeze; NE co-agree required)  
**Theory authority:** [`full_maxwell_monist_v0.md`](full_maxwell_monist_v0.md) **FROZEN** (no ontology fork)  
**Campaign:** [`../../CAMPAIGN_MAP.md`](../../CAMPAIGN_MAP.md) layer **M1**  
**Implements:** NE under `work/NE/sandbox_m1_*.py` → `outputs/m1_result.json`  
**M0 baseline (do not regress):** `sandbox_full_maxwell_r2.py` FM1–FM7  

---

## 0. Purpose and bar

### 0.1 Why M1 exists

M0 (Phase 1) established monist full Maxwell **skeleton**: Yee TE+TM, 1D TEM wave reduction, div B, Faraday identity, Cont ledger, 3D static Coulomb recovery.  

**M1 closes the “thin wave” residual:** prove **true 2D** dynamical free-gauge Maxwell with energy/Poynting, dynamic Gauss+sources, radiation honesty, and incomplete-Ampère adversary — still on **one free continuum** (`sector=1`, `E_origin=free_maxwell_full`).

### 0.2 Layer claim when green

| Claim | When all **mandatory** gates PASS |
|-------|-----------------------------------|
| `m1_claim` | `true` |
| `full_maxwell_2d_dynamic` | `true` |
| Unlocks | **CP-M1-NUM** → RC1, M2, M3 (per campaign map) |

### 0.3 Ontology tags (required on every export)

```text
sector=1
E_origin=free_maxwell_full
em_solver=free_maxwell_yee   # or free_maxwell_yee_m1
gravity_solver=none
embedding_dim_dynamics=2
monist_channel=free_gauge
c = 1/sqrt(eps*mu)          # C_LOCAL language shared with ND
```

Dualist adversary runs: `sector=2` (or explicit `adversary=incomplete_ampere`).

### 0.4 Units and solver defaults

| Item | Default | Note |
|------|---------|------|
| \(\varepsilon,\mu\) | 1, 1 → \(c=1\) | Plus **one** off-unit case (M1-G3) |
| Scheme | 2D Yee staggered TE\(_z\) and/or TM\(_z\) | Full 2D Maxwell (\(\partial_z=0\)) |
| CFL | \(\le 1/\sqrt{2}\) for isotropic 2D; document | Prefer stable ≤ 0.95 of limit |
| Grid | \(N_x,N_y \ge 128\) recommended for beam | State \(dx,dy,dt\) in JSON |
| Precision | float64 preferred | Relax residual floors ×10 for float32 if declared |

### 0.5 M0 regression (pre-flight)

| ID | Check | PASS |
|----|-------|------|
| **M1-R0** | Re-run or load M0 suite (`sandbox_full_maxwell_r2.py` / r2_result) | All M0 FM* still PASS; no claim M1 if M0 broken |

---

## 1. Mandatory gate list (CP-M1-NUM bar)

**`m1_claim = true` iff**  

\[
\text{M1-R0} \land \text{M1-G1} \land \text{M1-G2} \land \text{M1-G3} \land \text{M1-G4} \land \text{M1-G5} \land \text{M1-G6} \land \text{M1-G7} \land \text{M1-G8} \land \text{M1-G9}
\]

Optional: M1-G10 radiation outward flux (recommended, not blocking alone if G2+G4 strong).

---

### M1-G1 — Vacuum control

| | |
|--|--|
| **Maps** | Campaign #1; R2 FM1 |
| **Setup** | \(\rho_Q=\mathbf{J}_Q=\mathbf{0}\); \(\mathbf{E}=\mathbf{B}=\mathbf{0}\) (or pure gauge zero fields) |
| **Evolve** | ≥ 50 steps free Maxwell update |
| **PASS** | \(\max|\mathbf{E}| \le 10^{-12}\) and \(\max|\mathbf{B}| \le 10^{-12}\) (unitless; or \(\le 10^{-10}\) if float32 declared) |
| **FAIL** | Spontaneous fields without source |
| **Export** | `gates.M1_G1_vacuum` |

---

### M1-G2 — True 2D wave packet / beam (\(v \approx c\))

| | |
|--|--|
| **Maps** | Campaign #2 — **not** 1D TEM-only |
| **Setup** | Localized **2D** packet or finite-width beam in the plane: amplitude varies in **both** \(x\) and \(y\) (e.g. Gaussian envelope \(e^{-(x^2+y^2)/(2\sigma^2)}\) times transverse polarization; or beam with \(\sigma_y\) finite and \(k\) primarily in \(x\)). **Forbidden as sole evidence:** y-uniform / 1D TEM reduction used for M0 wave speed. |
| **Polarization** | TM\(_z\) (\(E_z,H_x,H_y\)) **or** TE\(_z\) (\(E_x,E_y,H_z\)); fields must be **coupled** \((\mathbf{E},\mathbf{B})\) dynamical |
| **Evolve** | Propagate long enough that peak travels \(\ge 8\sigma\) (or \(\ge 0.25 L_{\mathrm{box}}\)) |
| **Metric** | \(c_{\mathrm{th}}=1/\sqrt{\varepsilon\mu}\); \(v_{\mathrm{meas}}\) from peak track or envelope centroid |
| **PASS** | (i) \(|v_{\mathrm{meas}}/c_{\mathrm{th}} - 1| \le 0.08\) (8%; looser than 1D CFL=1 exact shift); (ii) packet remains 2D-localized (second moment in transverse direction stays \(O(\sigma^2)\), no collapse to y-uniform only); (iii) report `reduction=none` or `true_2d_packet` |
| **FAIL** | Only 1D TEM pass offered as M1-G2; \(v\) wrong by >8%; packet dies immediately from instability |
| **Export** | `gates.M1_G2_beam2d` with `v_meas`, `c_th`, `v_ratio`, `sigma`, track TSV optional |
| **Honesty** | M0 1D TEM may still appear as **regression** elsewhere; it does **not** satisfy M1-G2 |

---

### M1-G3 — Off-unit constitutive \(c\)

| | |
|--|--|
| **Maps** | Campaign #3 |
| **Setup** | Same 2D packet family as G2 (or shorter twin) with \((\varepsilon,\mu)\neq(1,1)\), e.g. \(\varepsilon=4,\mu=1\Rightarrow c_{\mathrm{th}}=0.5\) (preferred, matches M0 off-unit) |
| **PASS** | \(|v_{\mathrm{meas}}/c_{\mathrm{th}}-1|\le 0.08\); \(c_{\mathrm{th}}=1/\sqrt{\varepsilon\mu}\) must be used in CFL |
| **FAIL** | Speed stuck at \(c=1\) independent of \(\varepsilon,\mu\); or Courant number alone sets \(v\) without constitutive \(c\) |
| **Export** | `gates.M1_G3_offunit` |

---

### M1-G4 — Energy / Poynting (vacuum free stress)

| | |
|--|--|
| **Maps** | Campaign #4; theory §3.6 Poynting |
| **Setup** | Source-free 2D packet (G2 IC class) in **periodic** domain **or** measurement region well inside absorbing layers |
| **Definitions** | \(u=\frac{\varepsilon}{2}|\mathbf{E}|^2+\frac{1}{2\mu}|\mathbf{B}|^2\), \(\mathbf{S}=\mathbf{E}\times\mathbf{H}\), \(U=\int u\,dA\), discrete consistent with Yee staggering |
| **PASS** | (i) **Periodic vacuum:** \(|U(t)-U(0)|/U(0) \le 0.02\) (2%) over travel time of G2, **or** relative drift per light-crossing \(\le 0.005\); (ii) **Poynting identity residual** (discrete): domain integral of \(\partial_t u + \nabla\cdot\mathbf{S}\) consistent with \(O(h^2)\) / CFL truncation — mean absolute residual \(\le 0.05\,u_{\mathrm{peak}}/\Delta t\) when normalized, **or** show flux through closed curve accounts for \(\Delta U\) within 5% of \(|\Delta U|\) on a sub-box tracking the packet; (iii) \(U>0\) while packet present |
| **FAIL** | Energy grows secularly without source; \(U\) collapses while fields nonzero; no Poynting diagnostic at all |
| **Export** | `gates.M1_G4_energy_poynting` with `U0`, `U_end`, `rel_drift`, optional residual |
| **Note** | With PML/open BC, global \(U\) **need not** conserve — then PASS via **sub-volume** energy balance + outward flux (see G9/G10) |

---

### M1-G5 — Dynamic Gauss + Cont-safe sources

| | |
|--|--|
| **Maps** | Campaign #5; M3 + Cont |
| **Setup** | Prescribed compact \(\rho_Q(\mathbf{x},t)\), \(\mathbf{J}_Q\) with **analytic Cont:** \(\partial_t\rho_Q+\nabla\cdot\mathbf{J}_Q=0\) (e.g. drifting Gaussian \(\rho=\rho_0(\mathbf{x}-\mathbf{v}t)\), \(\mathbf{J}=\rho\mathbf{v}\), \(|\mathbf{v}|<c\)) |
| **Evolve** | ≥ 40 steps with full Maxwell driven by \(\mathbf{J}_Q\) (Ampère–Maxwell) |
| **PASS** | (i) Charge ledger: \(|Q(t)-Q(0)|/|Q_0| \le 10^{-3}\) interior (no BC leak), or flux-corrected \(Q\) conserved to same tol; (ii) Gauss residual \(r_G=\|\nabla\cdot(\varepsilon\mathbf{E})-\rho_Q\|_\infty / (\|\rho_Q\|_\infty+\epsilon)\): initial after consistent projection \(\le 0.05\); **does not grow** by more than factor 2 relative to initial over the run (or stays \(\le 0.08\)); (iii) Cont residual of prescribed fields at discrete floor or \(\le 0.02\) relative |
| **FAIL** | \(Q\) drifts \(O(1)\); Gauss residual secular \(O(1)\) growth with Cont-safe sources |
| **Export** | `gates.M1_G5_gauss_dynamic` |
| **Cleanse** | Optional electrostatic projection each N steps — if used, document; residual after cleanse still reported |

---

### M1-G6 — Divergence-free \(\mathbf{B}\)

| | |
|--|--|
| **Maps** | Campaign #6; M1; R2 KG-F1 |
| **Setup** | 2D dynamical run (G2 packet and/or G5 sources); IC with \(\mathrm{div}_h\mathbf{B}=0\) |
| **PASS** | \(\max|\mathrm{div}_h\mathbf{B}| / (\|\mathbf{B}\|_\infty/h + \epsilon) \le 10^{-6}\) (double) **and** no secular growth (end/start ≤ 10× if both near floor) |
| **FAIL** | divB grows as free dynamical error |
| **Export** | `gates.M1_G6_divB` |

---

### M1-G7 — Faraday (discrete)

| | |
|--|--|
| **Maps** | Campaign #7; M2; R2 KG-F2 |
| **Setup** | Yee update; optional continuum loop diagnostic on smooth pulse |
| **PASS** | (i) **Yee identity:** magnetic update is discrete Faraday (residual ≤ \(10^{-12}\) relative to Stokes scale), **or** (ii) continuum loop: \(|\oint\mathbf{E}\cdot d\boldsymbol{\ell} + \partial_t\Phi_B| / (|\oint\mathbf{E}|+|\partial_t\Phi_B|+\epsilon) \le 0.05\) time-averaged on smooth 2D pulse |
| **FAIL** | Faraday violated at \(O(1)\) under free Maxwell claim |
| **Export** | `gates.M1_G7_faraday` |

---

### M1-G8 — Incomplete-Ampère adversary

| | |
|--|--|
| **Maps** | Campaign #8; theory K-EM9 / R2 KG-F6 |
| **Setup** | **Control:** full M4 (with \(\partial_t(\varepsilon\mathbf{E})\)). **Adversary:** same code path with **displacement current removed** (magnetostatic Ampère \(\nabla\times\mathbf{H}=\mathbf{J}\) only, or \(\partial_t\mathbf{E}\) term zeroed in vacuum \(\mathbf{J}=0\)) |
| **Probe** | Same 2D free packet IC as G2 (vacuum) |
| **PASS** | (i) Full Maxwell: G2-class \(v\approx c\) **PASS**; (ii) Incomplete: either fails to propagate at \(c\) (\(|v/c_{\mathrm{th}}-1|>0.25\) or no traveling EM wave), **or** energy/causality diagnostic fails (fields freeze / non-radiating), **documented**; (iii) tags: full `sector=1`, adversary `adversary=incomplete_ampere` |
| **FAIL as M1** | Incomplete Ampère **also** “passes” free 2D wave as well as full Maxwell — then displacement current is untested |
| **FAIL as monism kill K-M1** | Only full Maxwell fails while dualist stage hand-wave is required (escalate to TU) |
| **Export** | `gates.M1_G8_ampere_adversary` with `full_pass`, `adversary_fails` |

---

### M1-G9 — Box / PML / BC honesty

| | |
|--|--|
| **Maps** | Campaign #9 |
| **Setup** | Document BC for **each** gate: `periodic` | `pec` | `pmc` | `pml` | `absorbing` | `large_box` |
| **PASS** | (i) Written note in `m1_result.json` → `bc_honesty` string **and** per-gate `bc` field; (ii) if PML/absorbing: show reflection figure-of-merit — e.g. late-time energy in launch half-domain after packet exits \(\le 0.1\,U_{\mathrm{peak}}\) **or** reflected amplitude ≤ 15% of incident for normal packet; (iii) if pure periodic: state that wrap-around can re-enter (G2 travel time < box light-crossing preferred); (iv) if large-box Dirichlet: measurement shell far from boundary |
| **FAIL** | Silent open BC with uncontrolled reflections sold as vacuum physics; no BC documentation |
| **Export** | `gates.M1_G9_bc_honesty` + top-level `bc_honesty` |
| **Not a physics fail** | Using periodic for energy (G4) and PML for radiation (G10) in **different** runs — preferred |

---

## 2. Recommended (not mandatory for `m1_claim`)

### M1-G10 — Radiation / outward flux

| | |
|--|--|
| **Setup** | Compact source turn-on or scattering proxy; absorbing/PML BC |
| **PASS** | Outward Poynting flux through large loop positive when source radiates; vacuum control still holds with sources off |
| **Export** | `gates.M1_G10_radiation` optional |

### M1-G11 — Shared-\(c\) cross-note for ND

| | |
|--|--|
| **PASS** | Report `C_LOCAL = c_th` unit case equals 1.0 (or declared shared value) for ND regression observation at CP-M1-NUM |
| **Export** | `shared_c` block |

---

## 3. Explicit non-goals (M1)

- 3D Yee (→ **M2**)  
- Self-consistent lock motion / Lorentz orbit (→ **M3** / RC2)  
- Joint \(\psi+(\mathbf{E},\mathbf{B})\) co-field (→ **RC1**; needs M1 green first)  
- Real \(\alpha_{\mathrm{EM}}\), SI units, `scp_sim` edits  
- Replacing \(\psi\) path-cost with \(\varepsilon(\rho_b)\) GRIN  

---

## 4. Kill conditions (M1 layer)

| ID | Kill if… | Escalation |
|----|----------|------------|
| **K-M1-num** | Mandatory gates systematically FAIL after honest redesign | DEFER CP-M1-NUM; redesign NE |
| **K-M1** (campaign) | 2D dynamical energy+wave **requires** dualist stage that breaks free/bound monism | TU DISPROVE path |
| **K-EM8/9** | Faraday/divB/displacement cannot live as free continuum | TE amend package only if proven; else kill proposal |

M0 success **does not** auto-pass M1 (1D TEM ≠ 2D beam).

---

## 5. Export schema (minimum `m1_result.json`)

```json
{
  "round": "P2-R2",
  "layer": "M1",
  "te_ref": "v77/work/TE/m1_gates_v0.md",
  "theory_ref": "v77/work/TE/full_maxwell_monist_v0.md",
  "tags": {
    "sector": 1,
    "E_origin": "free_maxwell_full",
    "em_solver": "free_maxwell_yee",
    "gravity_solver": "none",
    "embedding_dim_dynamics": 2
  },
  "m0_regression": { "M1_R0": true },
  "gates": {
    "M1_G1_vacuum": { "pass": true },
    "M1_G2_beam2d": { "pass": true, "v_ratio": 1.0, "reduction": "true_2d_packet" },
    "M1_G3_offunit": { "pass": true, "c_th": 0.5, "v_ratio": 1.0 },
    "M1_G4_energy_poynting": { "pass": true, "rel_drift": 0.0 },
    "M1_G5_gauss_dynamic": { "pass": true, "dQ_rel": 0.0 },
    "M1_G6_divB": { "pass": true },
    "M1_G7_faraday": { "pass": true },
    "M1_G8_ampere_adversary": { "pass": true, "full_pass": true, "adversary_fails": true },
    "M1_G9_bc_honesty": { "pass": true }
  },
  "bc_honesty": "periodic energy runs; ...",
  "m1_claim": true,
  "shared_c": { "C_LOCAL": 1.0, "c_from_eps_mu": 1.0 }
}
```

---

## 6. Demo IDs (for TU board)

| Demo ID | Gates | Status after SPEC |
|---------|-------|-------------------|
| **D-EM-M1-suite** | all mandatory | **SPEC** (await NUM) |
| **D-EM-M1-beam2d** | G2, G3 | SPEC |
| **D-EM-M1-energy** | G4 | SPEC |
| **D-EM-M1-gauss-dyn** | G5 | SPEC |
| **D-EM-M1-ampere-adv** | G8 | SPEC |
| D-EM-maxwell-full (M0) | M0 | LIVE (do not drop) |

---

## 7. Co-agreement (CP-M1-SPEC)

| Agent | Role |
|-------|------|
| **TE** | Author of this list; **STAMP CP-M1-SPEC: ADOPT** (provisional freeze for NE implement) |
| **NE** | Co-agree ADOPT / DEFER / REJECT after design review; then implement |
| **TU** | Note on CHECKPOINT_BOARD |

**TE rule:** SPEC is **provisionally frozen** at TE-ADOPT so NE can code; if NE DEFER with listed fixes, TE amends `m1_gates_v0.md` → v0.1 and re-stamps.  
**NUM stamp** only after `m1_result.json` exists (CP-M1-NUM).

---

## 8. Relation to prior gate sheets

| Prior | Relation to M1 |
|-------|----------------|
| `for_ne_kill_gates_r1.md` | Lite static + scalar wave — **M0-adjacent**; keep as regression |
| `for_ne_kill_gates_r2.md` | Full Maxwell skeleton (often 1D TEM wave) — **M0 bar** |
| **This file** | **M1 bar** — strict 2D beam + energy + dynamic Gauss + Ampère adversary + BC honesty |

---

## 9. Bottom line

**M1 = true 2D dynamic monist Maxwell**, not a re-label of M0 1D TEM.  
Mandatory: vacuum, **2D beam**, off-unit \(c\), energy/Poynting, dynamic Gauss+Cont, divB, Faraday, incomplete-Ampère fail, BC honesty.  
Theory remains `full_maxwell_monist_v0` freeze. NE implements; TE stamps NUM later.
