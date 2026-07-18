# Dual-Channel Matter Package v0 — Locks on Joint Free Medium

**Agent:** TM (Theory — Matter)  
**Date:** 2026-07-18  
**Round:** 2  
**Status:** V77-2 theory package — dual-channel congruence for Matter  
**Depends on:**  
- R1 `lock_ontology_v0.md` (lock S0–S3; force taxonomy; hierarchy skeleton)  
- R1 `FOR_NM_gates_v0.md` + NM G0–G3 **PASS** (`work/NM/outputs/r1_tm_gates.json`)  
- TE `maxwell_monist_v0.md` + `constitutive_table_v0.md` (M1–M4, Cont, JC1–JC5)  
- v76 F1-3D seed; O-002/O-003 (V77-2 dual-channel theme; full Maxwell progress)  
**Does not claim:** real \(G\)/\(\alpha\); full nonlinear media; SCP kernel monism; closed \(\alpha_\psi\) dynamics (MR2/J5 residual).

---

## 0. One-sentence claim

**Locks are bound continuum regions carrying independent ledgers \(\rho_b\) (mass-form) and \(\rho_Q\) (gauge charge); on one free medium they jointly source path-cost \(\psi\) and free-gauge \((\mathbf{E},\mathbf{B})\); multi-lock forces are virtual-work / Lorentz responses on those free channels with constitutive hierarchy — not two dualist stages.**

---

## 1. Round-1 → Round-2 upgrade map

| R1 (V77-1 Matter) | R2 (V77-2 dual-channel) |
|-------------------|-------------------------|
| DUAL-0: \(\psi\) + scalar \(\Phi\) Poisson lite | Joint medium: \(\psi\) + **Maxwell** (lite→full when NE ready) |
| Notation \(\rho_q\) | Align TE: **\(\rho_Q\)** (same object) |
| `free_gauge_lite` tag | Prefer `free_gauge` / M3–M4 origin; lite allowed as **subset** of static M3 |
| Separate conceptual siblings | **One sandbox state** \(\{\rho_f,\rho_b,\rho_Q,\psi,\mathbf{E},\mathbf{B}\}\) |
| G0–G3 PASS (analytic Green) | Joint gates **JG0–JG6** (this note) for V77-2 credit |
| Forces: \(F^\psi + F^C\) Coulomb only | Forces: \(F^\psi + F^{\mathrm{EM}}\) with Lorentz when \(\mathbf{B},\mathbf{J}\) live |

**R1 is not regressing.** Dual0 multipoles + force signs remain LIVE inputs. R2 adds joint constitutive identity, Maxwell upgrade path, and V77-2 kill criteria.

---

## 2. Lock ledgers on the joint medium

### 2.1 State (one continuum, one box)

\[
\boxed{
\bigl\{\,
\rho_f,\;
\rho_b,\;
\rho_Q,\;
\mathbf{J}_Q,\;
\psi,\;
\mathbf{E},\;
\mathbf{B},\;
\sigma,\; s,\; \gamma,\;
\varepsilon,\; \mu,\;
c
\,\bigr\}
}
\tag{STATE-DC}
\]

| Field | Role | Bounds / notes |
|-------|------|----------------|
| \(\rho_b\ge 0\) | Bound unlock density (mass-form ledger) | \(E_\star=\int\rho_b\) (v1 units); \(M=E_\star/c^2\) |
| \(\rho_Q\in\mathbb{R}\) | Gauge-charge density (lock–gauge ledger) | \(Q=\int\rho_Q\); TE M3 source |
| \(\mathbf{J}_Q\) | Gauge current | Continuity Cont; static demos may set \(\mathbf{J}_Q=\mathbf{0}\) |
| \(\rho_f\) | Free budget | Budget \(\rho_f+\rho_b=\rho_{\mathrm{tot}}\) (strong or integral) |
| \(\psi\) | Free-capacity / path-cost channel | F1-3D; not \(\Phi_{\mathrm{EM}}\) |
| \(\mathbf{E},\mathbf{B}\) | Free-gauge channel | TE M1–M4; quasistatic \(\mathbf{B}=\mathbf{0}\) OK for early gates |
| \(c\) | Shared free locality | \(c=1/\sqrt{\varepsilon\mu}\) **and** path-cost locality (JC1) |

### 2.2 Per-lock effective state

Lock \(L_i\):

\[
L_i = \bigl(\mathbf{X}_i,\; E_{\star,i},\; Q_i,\; \mathbf{V}_i,\; \mathcal{C}_i\bigr)
\]

with continuum reconstruction:

\[
\rho_b = \sum_i E_{\star,i}\, f_i(\mathbf{x}-\mathbf{X}_i),\qquad
\rho_Q = \sum_i Q_i\, g_i(\mathbf{x}-\mathbf{X}_i),
\]

\(f_i,g_i\) compact positive kernels, \(\int f_i=\int g_i=1\).

**Support law (monist charge on mass-form):**

\[
\boxed{
\mathrm{supp}(|\rho_Q|) \subseteq \mathrm{supp}(\rho_b)
}
\tag{Supp}
\]

**Forbidden collapses (TE + TM):**

| Collapse | Why forbidden |
|----------|----------------|
| \(\rho_Q \equiv \rho_b\) always | Mass = charge; kills neutral locks and sibling structure (JC3, K-EM5) |
| \(\rho_Q\) with \(\rho_b=0\) support | Free-floating charge without mass-form (unless TE defines free charge mode — out of Matter default) |
| \(\psi \equiv \Phi\) by parameter tuning | Sibling identification kill (TE-IA1, DS4) |

### 2.3 Neutral vs charged locks

| Species | \(E_\star\) | \(Q\) | Sources |
|---------|------------|-------|---------|
| Neutral lock | \(>0\) | \(0\) | \(\psi\) only (leading multipole) |
| Charged lock | \(>0\) | \(\neq 0\) | Both \(\psi\) and gauge |
| Opposite pair | \(>0,>0\) | \(+Q,-Q\) | Net \(Q_{\mathrm{tot}}=0\); dipole \(\Phi\); mass monopole \(\psi\) |

---

## 3. Free laws (dual channel)

### 3.1 Path-cost channel (unchanged seed)

\[
-\nabla\cdot\bigl(\sigma(\rho_f)\nabla\psi\bigr) = s\,\rho_b
\tag{F1}
\]

Path cost: \(\ell=\ell_0+\gamma\psi\). Exterior (3D linear): \(\psi\sim s E_\star/(4\pi\sigma_0 r)\).

### 3.2 Free-gauge channel — tiered Maxwell

**Tier A — Quasistatic Maxwell-lite (already congruent with R1 dual0 when \(\Phi\) from M3):**

\[
\nabla\cdot(\varepsilon\mathbf{E})=\rho_Q,\qquad
\mathbf{E}=-\nabla\Phi,\qquad
\mathbf{B}=\mathbf{0},\quad \mathbf{J}_Q=\mathbf{0}.
\tag{MX-A}
\]

**Tier B — Full free Maxwell (O-003 / campaign goal; TE M1–M4):**

\[
\begin{aligned}
\nabla\cdot\mathbf{B}&=0,\\
\nabla\times\mathbf{E}+\partial_t\mathbf{B}&=\mathbf{0},\\
\nabla\cdot(\varepsilon\mathbf{E})&=\rho_Q,\\
\nabla\times(\mathbf{B}/\mu)-\partial_t(\varepsilon\mathbf{E})&=\mathbf{J}_Q,\\
\partial_t\rho_Q+\nabla\cdot\mathbf{J}_Q&=0,\\
c&=1/\sqrt{\varepsilon\mu}.
\end{aligned}
\tag{MX-B}
\]

**Matter rule:** forces and gates must be written so **Tier A is a consistent limit of Tier B**, and Tier B is used when NE provides time-dep \((\mathbf{E},\mathbf{B})\).

### 3.3 Joint constitutive constraints (from TE JC1–JC5 + Matter)

| ID | Statement | Matter ownership |
|----|-----------|------------------|
| **JC1** | Single \(c\): path-cost locality \(\equiv 1/\sqrt{\varepsilon\mu}\) | Verify in joint sandbox; fail = V77-2 stress |
| **JC2** | Budget \(\rho_f+\rho_b\) holds; \(u_{\mathrm{EM}}\) is free-channel stress, not a second mass substance | Locks deplete \(\rho_f\) via \(\rho_b\) only |
| **JC3** | \(\rho_b\not\equiv\rho_Q\) | Enforce Supp + independent \(Q_i,E_{\star,i}\) |
| **JC4** | \(G_{\mathrm{eff}}\) and Coulomb/\(k_C\) from free constitutive tuples only | Hierarchy §5 |
| **JC5** | Sector tags: `phi_origin=free_capacity_*`, `gauge_origin=free_maxwell_*` | No silent dualist gravity solver |

**Joint static system (V77-2 minimal theory close):**

\[
\boxed{
\begin{aligned}
-\nabla\cdot(\sigma\nabla\psi) &= s\rho_b,\\
\nabla\cdot(\varepsilon\mathbf{E}) &= \rho_Q,\\
\rho_f+\rho_b&=\rho_{\mathrm{tot}},\\
\mathrm{supp}(|\rho_Q|)&\subseteq\mathrm{supp}(\rho_b),\\
c_{\mathrm{path}} &= c_{\mathrm{EM}} = 1/\sqrt{\varepsilon\mu}.
\end{aligned}
}
\tag{DUAL-1}
\]

DUAL-1 ⊃ R1 DUAL-0 + explicit shared-\(c\) identity.

---

## 4. Multi-lock forces (ψ + EM)

### 4.1 Decomposition

\[
\mathbf{F}_i = \mathbf{F}_i^{\psi} + \mathbf{F}_i^{\mathrm{EM}} + \mathbf{F}_i^{\mathrm{core}}
\tag{F-tot}
\]

| Term | Meaning | Round-2 default |
|------|---------|-----------------|
| \(\mathbf{F}^\psi\) | Path-cost / free-capacity virtual work | Virtual-work energy gradient (R1 PASS) |
| \(\mathbf{F}^{\mathrm{EM}}\) | Free-gauge force on lock charge/current | Coulomb (Tier A) or full Lorentz (Tier B) |
| \(\mathbf{F}^{\mathrm{core}}\) | Overlap / nonlinear \(\sigma(\rho_f)\) / contact | Deferred; tag if used |

### 4.2 Path-cost force (Tier-independent)

With free-capacity energy functional consistent with F1 linear vacuum (R1 NM form):

\[
U_\psi = \frac{s}{2}\int \psi\,\rho_b\,dV
\quad\Rightarrow\quad
\mathbf{F}_i^{\psi} = -\nabla_{\mathbf{X}_i} U_\psi
\quad\text{(exclude self)}.
\tag{F-ψ}
\]

Far-field two-body (3D):

\[
\mathbf{F}_{ij}^{\psi} \approx -\,\frac{s^2}{4\pi\sigma_0}\frac{E_{\star,i}E_{\star,j}}{r^2}\,\hat{\mathbf{r}}_{ij}
\quad(\text{attractive for }E_\star>0).
\]

**Optical path-cost** \(\ell=\ell_0+\gamma\psi\) remains the ray observable; relating \(\gamma\) to inertial/virtual-work coefficient is **MR2** (TD/J5), not required for V77-2 multipole+force-sign pass.

Tag: `force_closure_psi=virtual_work_energy_gradient_v0`.

### 4.3 EM force — Tier A (Coulomb)

\[
\mathbf{F}_i^{C} = Q_i\,\mathbf{E}_{-i}(\mathbf{X}_i)
\quad\text{or}\quad
U_C=\frac{\varepsilon}{2}\int|\mathbf{E}|^2
\Rightarrow
F_i=-\partial_{X_i}U_C.
\tag{F-C}
\]

Sign: \(Q_i Q_j>0\) repel; \(<0\) attract. Neutral: \(Q=0\Rightarrow F^C=0\).

Tag: `force_closure_em=coulomb_virtual_work_v0`.

### 4.4 EM force — Tier B (full Maxwell / Lorentz)

When \(\mathbf{B}\) and \(\mathbf{J}_Q\) (or lock velocity \(\mathbf{V}_i\)) are live:

\[
\boxed{
\mathbf{F}_i^{\mathrm{EM}}
=
\int\bigl(\rho_Q\,\mathbf{E} + \mathbf{J}_Q\times\mathbf{B}\bigr)_i\,dV
\;\approx\;
Q_i\bigl(\mathbf{E}_{-i}+\mathbf{V}_i\times\mathbf{B}_{-i}\bigr)
}
\tag{F-L}
\]

**Requirements for claiming “full Maxwell force” in Matter demos:**

1. \(\mathbf{E},\mathbf{B}\) satisfy M1–M4 residual floors (NE ownership).  
2. Continuity Cont held for moving locks / prescribed \(\mathbf{J}_Q\).  
3. Self-field regularization documented (exclude self or soft kernel).  
4. Static limit \(\mathbf{V}=\mathbf{0},\mathbf{B}=\mathbf{0}\) recovers F-C.

**Optional radiation reaction / Poynting momentum:** out of V77-2 minimum; may appear as residual.

Tag: `force_closure_em=lorentz_v0` when (F-L) used.

### 4.5 Taxonomy (joint, unchanged signs)

| Config | \(F^\psi\) | \(F^{\mathrm{EM}}\) (static) |
|--------|-----------|-------------------------------|
| Neutral–neutral | attract | \(0\) |
| Like charge | attract | repel |
| Opposite charge | attract | attract |
| Vacuum | \(0\) | \(0\) |
| Moving charges (Tier B) | attract (quasistatic \(\psi\)) | Lorentz; magnetic if \(\mathbf{B}\neq 0\) |

---

## 5. Hierarchy (constitutive, joint medium)

### 5.1 Ratio (static leading multipole)

\[
\frac{|F^{C}|}{|F^{\psi}|}
=
\frac{k_C |Q_i Q_j|}{\kappa_\psi E_{\star,i} E_{\star,j}},
\qquad
k_C=\frac{1}{4\pi\varepsilon_0},\quad
\kappa_\psi=\frac{s^2}{4\pi\sigma_0}
\quad\text{(energy-gradient convention)}.
\tag{H-ratio}
\]

TE-style optical factor (path-cost rays) still uses \(G_{\mathrm{eff}}=\gamma s c^4/(8\pi\sigma_0)\); **do not silently mix** \(\gamma\) into \(F^\psi\) without stating convention (R1 NM residual).

### 5.2 Design knobs (keep both channels alive)

| Knob | Effect |
|------|--------|
| \(\varepsilon_0\downarrow\) or \(k_C\uparrow\) | EM stronger |
| \(s\downarrow\) or \(\sigma_0\uparrow\) | Path-cost force weaker |
| \(\lambda_Q=\|Q\|/E_\star\) | Species charge-to-mass-ledger ratio |
| Neutral composite \(Q_{\mathrm{tot}}=0\) | External EM multipoles higher-order; \(\psi\) monopole remains |

**Pass hierarchy without killing channels:** achieve \(|F^C|/|F^\psi|\gg 1\) for elementary charged pairs **with** \(s\neq 0\), \(\sigma_0<\infty\), \(\gamma\) free for rays.

### 5.3 Phenomenological stacking (world story, not fitted SI)

1. **Micro / charged locks:** EM dominates pair force (constitutive).  
2. **Neutral bulk / composites:** external force ≈ path-cost only.  
3. **Same \(c\):** EM waves and free locality share \(c\) (JC1) — hierarchy is **coupling**, not speed split.

---

## 6. V77-2 success criteria (Matter share)

Parent V77-2 bar (PROBLEM §5): *EM free gauge **and** path-cost \(\psi\) coexist with shared \(c\) and budget identity (theory + ≥1 joint numeric).*

**Matter theory side MET when this package is accepted and interfaces TE JC1–JC5.**  
**Matter numeric side MET when NM (with NE Maxwell fields as available) passes joint gates §7.**

| Claim for V77-2 | Evidence |
|-----------------|----------|
| One medium | Single STATE-DC sandbox |
| Dual sources | \(\rho_b\to\psi\), \(\rho_Q\to\) gauge; JC3 |
| Shared \(c\) | Parameter + wave/locality check |
| Budget | Free deficit with \(\rho_b\) |
| Forces joint | \(F^\psi+F^{\mathrm{EM}}\) on ≥2 locks |
| Not lite-only forever | Path to MX-B / Lorentz documented; Tier A acceptable for first joint PASS if tagged |

**Not required for first V77-2 joint PASS:** orbits S3; full Poynting; real constants; scp_sim.

---

## 7. FOR_NM — joint gates (V77-2)

**Priority:** implement under `work/NM/`; reuse R1 multilock; **replace or dual-tag** `free_gauge_lite` with TE/NE Maxwell field source when available.  
**Coordinate FOR_NE:** import \(\mathbf{E}\) (and \(\mathbf{B}\) if Tier B) from joint state, not a foreign second grid.

### 7.1 Demo IDs

| Demo ID | Gate set | V77-2 role |
|---------|----------|------------|
| `D-DUAL-channel` | JG0–JG4 | Primary joint medium demo (shared with TU/NE) |
| `D-MAT-dual1` | JG0–JG3 | Matter alias of dual-channel statics |
| `D-MAT-force-maxwell` | JG5 | EM force from Maxwell fields (A or B) |
| `D-MAT-hier-joint` | JG6 | Hierarchy on joint medium |
| R1 `D-MAT-*` | — | Remain LIVE; do not re-kill |

### 7.2 Gate table

#### JG0 — Joint state / budget (P0)

| Check | Pass |
|-------|------|
| JG0.1 | Single run exports \(\psi\) and \(\mathbf{E}\) (or \(\Phi\)) from **one** medium object |
| JG0.2 | \(\rho_f+\rho_b=\rho_{\mathrm{tot}}\) (or integral free deficit \(\sim E_\star\)) |
| JG0.3 | \(\mathrm{supp}(|\rho_Q|)\subseteq\mathrm{supp}(\rho_b)\) |
| JG0.4 | Tags: `phi_origin=free_capacity_*`, `gauge_origin=free_maxwell_*` or `free_gauge_lite` (lite only if NE full not ready; document) |

#### JG1 — Dual multipoles (P0)

| Check | Pass |
|-------|------|
| JG1.1 | \(\psi\) exterior prefers \(1/r\) (3D); vacuum \(\rho_b=0\Rightarrow\psi=0\) |
| JG1.2 | \(\|\mathbf{E}\|\) or \(\Phi\) exterior Coulomb multipole; vacuum \(\rho_Q=0\Rightarrow\mathbf{E}=0\) |
| JG1.3 | Sibling independence: scale \(E_\star\) at fixed \(Q\) moves \(\psi\) not \(Q\); scale \(Q\) at fixed \(E_\star\) moves \(\mathbf{E}\) not \(E_\star\) |
| JG1.4 | Opposite-charge pair: \(\psi\) monopole from total \(M\); \(\Phi\) monopole cancelled (dipole) |

#### JG2 — Shared \(c\) (P0 for V77-2)

| Check | Pass |
|-------|------|
| JG2.1 | Constitutive \(c_{\mathrm{EM}}=1/\sqrt{\varepsilon\mu}\) equals path-cost locality parameter \(c_{\mathrm{path}}\) (same number in config) |
| JG2.2 | If free wave available (NE/ND): measured \(c_{\mathrm{wave}}\approx c\) within tol; **else** tag `c_wave=deferred` and pass JG2.1 only with `shared_c=constitutive` |

#### JG3 — Force taxonomy on joint fields (P0)

| Check | Pass |
|-------|------|
| JG3.1 | Neutral–neutral: \(F^{\mathrm{EM}}\approx 0\), \(F^\psi\) attract |
| JG3.2 | Like charge: \(F^{\mathrm{EM}}\) repel, \(F^\psi\) attract |
| JG3.3 | Opposite: both attract (static) |
| JG3.4 | Vacuum forces 0 |
| JG3.5 | \(F^{\mathrm{EM}}\) computed from **joint** \(\mathbf{E}\) (not a separately fitted \(\kappa_c\) only) — at least one config |

#### JG4 — No sibling collapse (P0 kill)

| Check | Pass / Fail |
|-------|-------------|
| JG4.1 | Neutral mass does **not** source static \(\mathbf{E}\) (beyond numeric floor) |
| JG4.2 | Opposite charges do **not** source opposite \(\psi\) |
| JG4.3 | Fit \(\psi\approx a\Phi\) forced globally fails when \(Q/M\) varies across locks |

**Fail any JG4 = dual-channel ontology kill stress (report to TU; not automatic V77-K alone).**

#### JG5 — Maxwell upgrade (P1; full Maxwell path)

| Check | Pass |
|-------|------|
| JG5.A | Tier A: \(\mathbf{E}=-\nabla\Phi\) from M3 on joint \(\rho_Q\); Gauss residual at floor |
| JG5.B | Tier B (when NE ready): time-dep M1–M4; Cont; optional moving lock or \(\mathbf{J}_Q\); Lorentz \(Q(\mathbf{E}+\mathbf{V}\times\mathbf{B})\) sign check |
| JG5.C | Static limit of Tier B matches Tier A forces within tol |

**V77-2 Matter minimum:** JG5.A PASS.  
**Campaign full-Maxwell Matter credit:** JG5.B+C PASS.

#### JG6 — Hierarchy constitutive (P1)

| Check | Pass |
|-------|------|
| JG6.1 | Parameter set with \(|F^C|/|F^\psi|>10^2\) for charged pair **without** \(s=0\) or deleting \(\psi\) solve |
| JG6.2 | Neutral pair still \(F^\psi>0\), \(F^{\mathrm{EM}}\approx 0\) at same constitutive \(\varepsilon,\sigma,s\) |
| JG6.3 | Report convention: energy-gradient \(\kappa_\psi\) vs optical \(\gamma\) (do not claim equality without TD) |

### 7.3 Suggested NM implementation path

```text
1. Fork R1 offline_r1_multilock / sandbox → sandbox_r2_dual_channel.py
2. STATE-DC: keep ρ_b, ρ_Q, ψ; set ε,μ,c with c=1/√(εμ)
3. Phase R2a: joint static DUAL-1 (analytic Green or SOR) → JG0–JG4, JG5.A, JG6
4. Phase R2b: plug NE Maxwell E (and B) fields → recompute F_EM; JG5.B if available
5. Export: outputs/r2_dual_result.json, r2_tm_joint_gates.json
```

Adapt NE Coulomb/wave modules by **import/copy into joint runner**, not by claiming two separate Occam wins as V77-2.

### 7.4 Report tags

```text
FOR_TM: JG0=... JG1=... JG2=... JG3=... JG4=... JG5=... JG6=...
FOR_TU: D-DUAL-channel status=...
FOR_NE: maxwell_tier=A|B; gauss_residual=...; shared_c=...
```

---

## 8. Interfaces

| Agent | Message |
|-------|---------|
| **NM** | Joint gates §7; R1 PASS inherited; V77-2 needs JG0–JG4 minimum |
| **NE** | Provide Maxwell \(\mathbf{E}\) (Tier A) / \((\mathbf{E},\mathbf{B})\) (Tier B) on shared grid; Cont; do not require \(\psi\) inside pure EM gates |
| **TE** | Notation \(\rho_Q\); DUAL-1 = TE joint §4.4 + Matter Supp + forces; freeze M1–M4 |
| **TD/ND** | MR2/\(\alpha_\psi\); moving-lock Lorentz ↔ free drag later; J5-β default does not block dual-channel statics |
| **TU** | Register `D-DUAL-channel` / `D-MAT-dual1`; promote when NM joint PASS; score V77-2 |
| **O** | Matter R2 theory ready for dual-channel congruence; full Maxwell force = JG5.B path |

---

## 9. Residuals and kills

### 9.1 Residuals (honest)

| ID | Residual |
|----|----------|
| MR1 | Poisson-form isomorphism each channel (Occam/tags) |
| MR2 | \(\gamma\) optical vs virtual-work \(F^\psi\) coefficient (TD/J5) |
| MR3 | R1 analytic Green vs full SOR systematics |
| MR4 | Tier B Lorentz + continuity on moving locks not yet numeric |
| MR5 | SCP Noether \(Q\) ↔ monist \(\rho_Q\) still vocabulary |

### 9.2 Kill conditions (Matter dual-channel)

| ID | Kill if… |
|----|----------|
| DC1 | Joint demo requires two non-communicating continua |
| DC2 | Only works if \(\rho_Q\equiv\rho_b\) or \(\psi\equiv\Phi\) |
| DC3 | Shared \(c\) impossible without destroying multipole of one channel |
| DC4 | Hierarchy only by deleting a free channel |
| DC5 | Neutral locks inevitably source Coulomb monopole under Supp+F1+M3 |

---

## 10. SCP vocabulary (unchanged non-claim)

Ball / Q-ball ↔ lock; kernel \(Q\) ↔ monist \(\rho_Q\) **analogy only**.  
**Fixed-grid Q-balls do not prove dual-channel monism.** No `scp_sim` edits.

---

## 11. Bottom line

- **R1:** lock ontology + dual Poisson forces + hierarchy — **numeric PASS**.  
- **R2:** **DUAL-1** joint medium with \((\rho_b,\rho_Q)\), shared \(c\), forces from \(\psi\) + Maxwell (Coulomb now, Lorentz when ready), constitutive hierarchy, **JG0–JG6** for NM.  
- **V77-2 Matter theory:** this package.  
- **V77-2 Matter numeric:** NM joint PASS (minimum JG0–JG4 + JG5.A).
