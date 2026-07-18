# ψ Time Dependence for RC1 — Fixed Locks (T0 / T1 / T2)

**Agent:** TD (Theory — Dynamics)  
**Round:** P2-R3  
**Checkpoint:** **CP-RC1-SPEC**  
**Status:** **FROZEN recommendation** for RC1 co-field  
**Stamp:** `STAMP CP-RC1-SPEC: ADOPT` (TD vote on ψ-sector law)  
**Date:** 2026-07-19  

**Depends on:**  
- `dynamics_package_v0.md` §2 (T0–T4 option table; TD-S stationary recovery)  
- `j5_beta_default_v0.md` (J5-β locked; dual-channel time-dep sketch)  
- v76 F1-3D seed (`goal2_PC3D_workable`)  
- `CAMPAIGN_MAP.md` §4.3 D1; RC1 = fixed multi-locks, dynamical Maxwell + ψ  

**Partners:** TM/NM (joint state + force diagnostics); NE (Maxwell API); ND (shared \(c\), optional T1 probe)  
**Does not own:** Maxwell Yee (TE/NE); lock Supp ontology (TM); RC2 moving locks (later `lock_motion_rc2_v0.md`)

---

## 0. RC1 scope (what is fixed)

From campaign map:

| Item | RC1 | Not RC1 |
|------|-----|---------|
| Locks | **Fixed** multi-lock centers | Moving locks (RC2) |
| Maxwell | **Dynamical** \(\mathbf{E},\mathbf{B}\) via NE API | Φ-only lite as sole EM |
| ψ | Free-capacity / path-cost channel | Identified with Φ (TE-IA1 kill) |
| Forces | Diagnostics \(F^\psi + q\mathbf{E}\) on fixed sites | \(F=ma\) evolution of \(X_i(t)\) |
| Shared \(c\) | \(c_{\mathrm{EM}}=1/\sqrt{\varepsilon\mu}=C_{\mathrm{LOCAL}}\) | New free speed |

**Implication for ψ law:** with \(\rho_b(\mathbf{x})\) **static** (fixed locks), the free-capacity sector only needs a law whose **end-state is F1-3D**. True retardation / wake is **optional** in RC1; it becomes mandatory for RC2 inertia / motion.

---

## 1. Option recap (from dynamics_package_v0)

| ID | Law (schematic) | Static limit | Free signals on ψ | RC1 fitness |
|----|-----------------|--------------|-------------------|-------------|
| **T0 Quasistatic** | Re-solve \(-\nabla\cdot(\sigma\nabla\psi)=s\rho_b\) each snapshot | Exact F1 | None (infinite \(c_{\mathrm{eff}}\) for ψ) | **Primary** for fixed locks |
| **T1 Relaxational** | \(\partial_t\psi=\kappa\nabla\cdot(\sigma\nabla\psi)+s\rho_b\) (sign so \(t\to\infty\Rightarrow\) F1) | F1 | Diffusive | **Allowed upgrade** |
| **T2 Hyperbolic** | \(c^{-2}\partial_{tt}\psi=\nabla\cdot(\sigma\nabla\psi/\sigma_0)+\cdots\) | F1 if \(\partial_{tt}\to0\) | Waves at \(c\) | **Not required** for RC1 fixed locks |
| **T3 Hybrid** | T0/T1 for ψ + Maxwell (or free \(u\)) for radiation | F1 + EM vacuum | EM carries \(c\) | **Default dual-channel picture** |
| **T4 Graph** | Discrete free hop | Graph Green | Causal hop | Deferred |

**Stationary recovery (hard, all options):**

\[
\partial_t\to 0,\quad \rho_b=\rho_b(\mathbf{x})
\quad\Rightarrow\quad
-\nabla\cdot(\sigma\nabla\psi)=s\rho_b
\quad\text{(F1-3D)}.
\tag{TD-S}
\]

**Kill:** any ψ law whose static exterior is not \(\sim 1/r\) in 3D (or fails free deficit at lock) — same as v76 congruence.

---

## 2. RC1 recommendation (locked)

### 2.1 Primary: **T0 quasistatic F1**

**For RC1 fixed multi-locks, TD recommends T0 as the default ψ solver.**

| Reason | Detail |
|--------|--------|
| Scope match | \(\rho_b\) fixed \(\Rightarrow\) one F1 solve (or warm-start re-solve) is exact for the path-cost sector |
| Congruence | Identical to v76 F1-3D seed already congruent with numerics |
| Occam | No free \(\kappa\) or Courant number for ψ; NM co-field complexity stays on Maxwell + dual sources |
| Forces | \(F^\psi=\int(-\rho_b\nabla\psi)\,dV\) (or multipole form) is well-defined on the static solution |
| Shared \(c\) | **Not carried by ψ dynamics** under T0; **carried by Maxwell** (NE M1) — T3 hybrid default |

**Implementation note for NM:**

```text
each RC1 frame / diagnostic step:
  given fixed ρ_b (and ρ_f budget tags):
    solve F1-3D for ψ   # T0: elliptic SOR / multigrid / whatever
  given ρ_Q, J_Q (often J=0 for fixed locks):
    step or relax Maxwell via NE API → E, B
  report F^ψ, F^EM, multipoles, shared-c tags
```

With fixed locks and \(J_Q=0\), Maxwell may itself be near-static Coulomb; the **campaign hard rule** still requires a **dynamical Maxwell path** available (not Φ-only ontology), even if the steady state is Coulomb-like.

### 2.2 Allowed: **T1 relaxational** (optional)

T1 is **ADOPT-compatible** for RC1 if partners want a time-dep free-capacity probe without moving locks:

\[
\partial_t\psi
=
\kappa\,\nabla\cdot(\sigma\nabla\psi)
+
\kappa s\,\rho_b
\quad\text{(convention: cold start } \psi=0 \to \text{ F1 as }t\to\infty\text{)}.
\tag{T1-RC1}
\]

| Requirement | Bar |
|-------------|-----|
| TD-S | After long \(t\), exterior \(\psi\sim 1/r\) in 3D; free deficit at core |
| Cold-start gate | From \(\psi=0\), residual of F1 PDE decreases monotonically within tolerance |
| No fake \(c\) | Do **not** claim T1 wavefront speed \(=c\); \(c\) remains Maxwell / ND free-wave sibling |
| κ documentation | Report \(\kappa\) and relaxation time vs light-crossing time of box |

**When to use T1 in RC1:** ND Dyn1 probe; or NM “turn on locks” transient demos. **Not** required for co-field force taxonomy PASS.

### 2.3 Deferred: **T2 hyperbolic ψ**

T2 remains **valid free-capacity radiation law** (ND R1 free-wave PASS on vacuum T2) but is **not recommended as RC1 co-field default**:

- Fixed \(\rho_b\) \(\Rightarrow\) no need for ψ wave dynamics for force diagnostics.  
- Shared \(c\) already demonstrated on **Maxwell** (NE M1) + ND dual-channel.  
- Hyperbolic ψ + Yee Maxwell doubles Courant constraints and confuses energy split.  

**Use T2 for:** ND free-capacity wave regression; RC2 / inertia wake studies; optional T3 if TE/TU demand capacity radiation distinct from EM.

### 2.4 Default dual-channel picture (T3 hybrid)

\[
\boxed{
\text{RC1 free medium}
=
\underbrace{\text{T0 (or T1) }\psi}_{\text{path-cost / mass multipole}}
\;+\;
\underbrace{\text{dynamical Maxwell }(\mathbf{E},\mathbf{B})}_{\text{gauge / shared }c}
}
\tag{RC1-ψ}
\]

| Forbidden | Why |
|-----------|-----|
| \(\psi\equiv\Phi\) or single Poisson for both | TE-IA1 / K-RC1 |
| T0 ψ + **only** dualist foreign gravity | Ontology kill |
| Claiming T0 proves full dynamics / inertia coefficient | J5-β residual; RC2 |

---

## 3. Interface to joint state (for TM/NM SPEC)

Minimal ψ-sector fields (alongside Maxwell):

```text
rho_b, rho_f          # budget
psi                   # free capacity
sigma0, s, gamma      # F1 constitutive (path cost)
psi_law               # "T0_quasistatic_F1" | "T1_relaxational"
kappa                 # only if T1
```

**Coupling to locks (RC1 fixed):**  
Sources \(\rho_b,\rho_Q\) held fixed; forces evaluated, centers not integrated.

**Energy tags (J5-β consistent):**  
- Path-cost / ray mass from \(\int\rho_b\).  
- Free-capacity field energy \(U_\psi=\frac{\sigma_0}{2}\int|\nabla\psi|^2\) available for diagnostics; not required equal to \(\int\rho_b\) (form factor).  
- EM free energy \(U_{\mathrm{EM}}\) from NE; sibling free ledger.

---

## 4. Kill-gates (theory-facing; NM/ND implement)

| ID | Gate | Pass |
|----|------|------|
| **RC1-ψ-0** | Law tag is T0 or T1 (not silent) | `psi_law` field present |
| **RC1-ψ-1** | TD-S: static end-state is F1 | Exterior multipole prefer \(1/r\) (3D); vacuum \(\rho_b=0\Rightarrow\psi\approx 0\) |
| **RC1-ψ-2** | Sibling: \(\psi\) responds to \(\rho_b\), not \(\rho_Q\) | Scale \(Q\) at fixed \(E_\star\): \(\psi\) amplitude stable |
| **RC1-ψ-3** | Shared \(c\) not assigned to T0 ψ | `c_shared` from Maxwell / constitutive; T0 does not invent \(c_\psi\neq c_{\mathrm{EM}}\) |
| **RC1-ψ-4** (if T1) | Cold start → F1 | Residual drop documented; no wave-speed \(=c\) claim |
| **RC1-ψ-K** | Kill | \(\psi:=\Phi\); static exterior wrong multipole; T2 required falsely for fixed locks |

**ND Dyn1 (after SPEC):** optional T1 cold-start or vacuum T2 regression — does **not** block NM T0 co-field.

---

## 5. Decision table (quick)

| Scenario | Choose |
|----------|--------|
| RC1 default co-field (fixed locks) | **T0** |
| RC1 + free-capacity transient demo | **T1** |
| Shared-\(c\) wave proof | **Maxwell (NE)** + ND dual-channel; not T0 |
| RC2 moving locks / wake inertia | T1 or T2 (open in `lock_motion_rc2_v0`) |
| J5-β A2 boost energy | T0 comoving re-solve (v76 R4) still valid control |

---

## 6. Stamp and co-agreement

**TD vote on CP-RC1-SPEC (ψ sector):**

```text
STAMP CP-RC1-SPEC: ADOPT
```

**Scope of this ADOPT:**  
ψ time-dep law for RC1 fixed locks = **T0 primary, T1 allowed, T2 deferred**; TD-S mandatory; T3 hybrid with dynamical Maxwell.

**Does not alone adopt full CP-RC1-SPEC:**  
Board co-stamps still required from **TM, TE, NM, NE** (joint state, Maxwell API, force law). TD ADOPT is the Dynamics half of ψ constitutive choice.

**Evidence path:** `work/TD/psi_time_dep_rc1_v0.md` (this file).

---

## 7. Residuals

| Residual | Owner |
|----------|--------|
| T1 κ scale vs \(c\) | TD/ND if T1 used |
| T2 capacity radiation vs EM | TE/TD dual-channel (optional) |
| Which \(m\) in \(F=ma\) | **RC2** — `lock_motion_rc2_v0.md` |
| Form factor / J5-β | Locked; not reopened for RC1 |

---

## 8. Bottom line

> **RC1 fixed locks: solve quasistatic F1 for ψ (T0). Optionally relax T1. Do not require hyperbolic ψ. Dynamical Maxwell carries shared \(c\). Static recovery TD-S is non-negotiable.**

**V77 / Phase-2:** Unblocks NM `sandbox_rc1_cofield.py` ψ sector without waiting for T2 production.
