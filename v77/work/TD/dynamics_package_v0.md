# Dynamics Package v0 — Time-Dependent Free Response, Inertia, Frames

**Agent:** TD (Theory — Dynamics)  
**Round:** 1  
**Date:** 2026-07-18  
**Status:** draft laws + kill-gates for ND; not a closed theorem  
**Seed:** v76 F1-3D (`goal2_PC3D_workable`); residuals R2 (J5 partial), R4 (rods/clocks)  
**Depends on:**  
- `v76/work/A/THEORETICAL_PACKAGE_v1.md`  
- `v76/work/A/inertia_target_v0.md`, `kinetic_inertia_v0.md`  
- `v76/work/B/outputs/round4_result.json` (J5 raw FAIL / renorm PASS)  
- `v76/work/C/rods_clocks_C1_v0.md`  
- `v77/PROBLEM.md`, `v77/FOCI.md` §3  

**Demo IDs (proposed):** `D-DYN-timedep`, `D-DYN-J5`, `D-DYN-rods`, `D-DYN-cwave`  
**Partner:** ND — implement kill-gates under `work/ND/` (adapt `v76/work/B/sandbox_r4_inertia.py`, free-wave probes)

---

## 0. Mandate and claim

**Claim (killable):**  
The v76 monist free continuum admits a **time-dependent free law** whose quasistatic limit is F1-3D; **inertial mass** is free-medium response to lock motion (not a hand-set \(F=ma\)); **rods/clocks** are free-medium operations so local \(c\) is medium-made; **Lorentz structure** is *not* assumed proved from locality-\(c\) alone in this package.

**Owns for v77 Dynamics focus:**

| Item | TD role | ND role |
|------|---------|---------|
| Time-dep free \(\psi\) | Law options + selection criteria | Implement ≥1; wave/relax tests |
| J5 inertia | Non-tautology targets; renorm necessity statement | Measure triad; kill naïve equality |
| Rods/clocks C1 | Operational sketch (from C) | Minimal free-signal \(c_{\mathrm{op}}\) probe |
| Lorentz-from-\(c\) | Caution + minimal claims | Do **not** claim full SR from \(c\) alone |

**Does not own:** Maxwell channel (TE/NE); lock microphysics (TM/NM); program kill (TU).

---

## 1. Shared primitives (seed — no re-proof)

Inherited without re-litigation (v76 CONVERGENCE):

| Symbol | Meaning |
|--------|---------|
| \(\rho_b,\rho_f\) | Bound / free budget of **one** continuum |
| \(\psi\) | Free capacity potential (free DOF) |
| \(c\) | Free-field locality bound (local frames) |
| \(E_\star\) | Unlock / bound ledger \(\int\rho_b\,dV\) (v1 units) |
| \(M_{\mathrm{ledger}}=E_\star/c^2\) | Ledger mass |
| F1-3D | \(-\nabla\cdot(\sigma\nabla\psi)=s\rho_b\); exterior \(\psi\sim 1/r\) in 3D |
| \(\ell=\ell_0+\gamma\psi\) | Path cost; \(G_{\mathrm{eff}}=\gamma s c^4/(8\pi\sigma_0)\) |

**Dead (do not reopen as monist proof):** local GRIN long-range; hand Poisson F5 as theory; 2D log as Einstein exterior.

---

## 2. Time-dependent free response options

Quasistatic F1-3D is the **stationary slice**. Dynamics needs a law for \(\partial_t\psi\) (and free signals) so that:

1. Stationary locks \(\Rightarrow\) F1-3D (or M2 free-Laplace with lock BC).  
2. Free disturbances propagate / relax with **locality bound \(c\)**.  
3. Moving locks produce **wake / retardation** that can source kinetic free energy (inertia).  
4. No second gravity sector; \(\psi\) remains free continuum.

### 2.1 Option table

| ID | Law (schematic) | Stationary limit | Free signals | Inertia path | Risk |
|----|-----------------|------------------|--------------|--------------|------|
| **T0 Quasistatic** | Re-solve F1 each time for \(\rho_b(t)\) | Exact F1 | None (infinite \(c_{\mathrm{eff}}\) for \(\psi\)) | Only geometric comoving \(U[\psi_v]\) if boosted by hand | **No true dynamics**; A2 boost energy may miss radiation/drag |
| **T1 Relaxational** | \(\partial_t\psi=\kappa\nabla\cdot(\sigma\nabla\psi)+s\rho_b\) | F1 as \(t\to\infty\) | Diffusive, not ballistic | Wake at finite \(\kappa\); low-\(\omega\) drag | No sharp light-cone; \(c\) must enter via free-signal **sibling** channel or \(\kappa\sim c^2/\Lambda\) |
| **T2 Hyperbolic free capacity** | \(\frac{1}{c^2}\partial_{tt}\psi=\nabla\cdot(\sigma\nabla\psi/\sigma_0)+\frac{s}{\sigma_0}\rho_b\) (or KG-like) | F1 if \(\partial_{tt}\to0\) and \(\rho_b\) static | Waves at \(c\) (or \(c\sqrt{\sigma/\sigma_0}\)) | Field momentum / boost energy \(\propto U/c^2\) natural | Stage Minkowski residual if operator coefficients not free-state |
| **T3 Hybrid** | Hyperbolic free waves \(u\) **sibling** of elliptic/relaxational \(\psi\) | \(\psi\) → F1; \(u\) vacuum | \(u\) carries free radiation at \(c\) | Inertia from \(\psi\) energy + optional \(u\)-drag | Two free DOFs — must share budget/ontology (TE dual-channel) |
| **T4 Graph hop** | Free voltage / hop rates with max hop \(c\) | Graph Green → continuum Green | Causal free commute | Resistance / memory cost of moving sink | Heavier numerics; lattice ≠ rod (C1) |

### 2.2 Round-1 theory preference

| Priority | Choice | Why |
|----------|--------|-----|
| **P0 for J5** | **T0 + analytic/boost model** (v76 R4) then **T1 or T2 on grid** | Reproduces residual; then upgrades non-tautology |
| **P0 for free radiation / \(c\)** | **T2** (or T3 if TE owns sibling EM waves) | Direct locality-\(c\) test: free wavefront speed |
| **Avoid as sole dynamics** | **T0 alone** for claiming “closed dynamics” | Infinite signal speed for \(\psi\) |

**Recommended ND stack (Round 1):**

1. **T0 control:** comoving F1 solve for \(\rho_b(x-vt)\); measure \(\Delta U(v)\) (inherit R4).  
2. **T1 or T2 lite:** 1D or 3D small box; source pulse or moving Gaussian; measure wake and/or wavefront.  
3. Do **not** require full 3D hyperbolic production in Round 1 if T0+T1 already kill-gated.

### 2.3 Stationary recovery (hard requirement)

Any accepted time-dep law must satisfy:

\[
\partial_t\to 0,\ \rho_b=\rho_b(\mathbf{x})
\quad\Rightarrow\quad
-\nabla\cdot(\sigma\nabla\psi)=s\rho_b
\quad\text{(F1-3D)}.
\tag{TD-S}
\]

**Kill:** time-dep law whose static end-state exterior is not \(\sim 1/r\) in 3D (or fails free deficit at lock).

### 2.4 Free radiation vs path-cost statics

| Channel | Role |
|---------|------|
| Static / slow \(\psi\) | Path-cost monopole, \(M_{\mathrm{ray}}\), weak “gravity-class” optics |
| Dynamic free waves (T2/T4 or free \(u\)) | Radiation, causality, operational rods/clocks, possible EM sibling (TE) |

**Slogan:**  
> Path cost is the **slow free capacity** bookkeeping; free waves are **fast free transport**. Same continuum budget; different constitutive response.

Interface to TE: Maxwell-lite may **be** the free-wave sibling (T3), not a replacement for \(\psi\). TD does not absorb TE’s constitutive \(\varepsilon,\mu\); TD only requires a free channel with speed \(c\).

---

## 3. Non-tautological inertia targets (J5)

### 3.1 What v76 left open (facts, not re-run)

From `round4_result.json` / `round4_summary.txt` (B R4):

| Tier | Comparison | Result |
|------|------------|--------|
| A | \(M_{\mathrm{ray}}\) vs \(M_{\mathrm{ledger}}\) | **PASS** (by F1 construction) |
| B | \(m_{\mathrm{boost}}\) vs \(M_{\mathrm{ledger}}\) (raw \(s=1\)) | **FAIL** (\(U/Mc^2\approx 0.141\) form factor) |
| C | \(m_{\mathrm{boost}}\) vs \(m_{\mathrm{field}}=U/c^2\) | **PASS** |
| Raw three-way J5 | | **FAIL** |
| R1 renorm \(s\to s_*\) so \(U=Mc^2\) | | **PASS** (calibration, not dynamics-free proof) |

**Honest residual:** free-field self-energy \(U\) tracks **form factor**, not raw \(\int\rho_b\), unless free constants are renormalized or the ledger is redefined as free unlock energy \(E_\star:=U\).

### 3.2 Target equality (triad — restated for v77)

\[
\boxed{
m_{\mathrm{inertial}}
\;\stackrel{?}{=}\;
\frac{E_\star}{c^2}
\;\stackrel{?}{=}\;
M_{\mathrm{ray}}
}
\tag{J5}
\]

| Symbol | Operational (non-tautological) meaning |
|--------|----------------------------------------|
| \(E_\star\) | Bound ledger **or** measured unlock work (tag which) |
| \(M_{\mathrm{ray}}\) | Exterior path-cost amplitude with **fixed** \(G_{\mathrm{eff}}\) from free constants (no free-fit \(M\)) |
| \(m_{\mathrm{inertial}}\) | From free-field energy expansion / free force integrals — **never** from \(a=F/m_{\mathrm{ledger}}\) |

### 3.3 Forbidden protocols (tautology / dualism)

| Forbidden | Why |
|-----------|-----|
| Set \(a=F/m_{\mathrm{ledger}}\) and report \(m=m_{\mathrm{ledger}}\) | Tautology |
| Fit \(m\) from dualist \(\Phi\) with foreign \(G\) | Dualist sector |
| Prescribe Newtonian trajectory and back out \(m\) | Circular |
| Call ray slope alone “inertia” | That is \(M_{\mathrm{ray}}\), not independent \(m_{\mathrm{inertial}}\) |
| Insert \(m=E_\star/c^2\) into the free integrator by hand | Smuggled coefficient |

### 3.4 Allowed protocols (inherit A R4; TD endorses)

| ID | Protocol | Extracts |
|----|----------|----------|
| **A2** | Slow boost: translate \(\rho_b(x-vt)\); re-solve/relax \(\psi\); \(\Delta U(v)=\tfrac12 m_{\mathrm{inertial}}v^2+\cdots\) | Preferred Round-1 inertia mass |
| **A3** | Free drag under external free drive \(\psi_{\mathrm{ext}}\); \(m=F_\psi/\ddot X\) with \(F_\psi=\int(-\rho_b\nabla\psi)\,dV\) | Easy to code wrong; secondary |
| **B unlock** | Work to dissolve lock vs \(\int\rho_b\) | Ledger cross-check |
| **C ray** | Fixed \(G_{\mathrm{eff}}\) multipole mass | Already closed at F1 |

### 3.5 Theory fork: close J5 vs renorm necessity

Two **honest** closures (either is V77-3 success if numeric confirms):

#### Fork **J5-α** — Raw coefficient 1 (hard)

Require, without retuning \(s\) after measuring \(U\):

\[
m_{\mathrm{inertial}} = \frac{E_\star}{c^2} = M_{\mathrm{ray}}
\quad\text{within band }\varepsilon_{\mathrm{triad}}.
\]

**Theory status:** **not expected** for bare Gaussian + default free constants (v76 form factor \(\sim 0.14\)). Would need:

- lock form / constitutive law that forces \(U=E_\star\), or  
- hyperbolic free action + budget identity that equates free self-energy to unlock ledger (missing step in `kinetic_inertia_v0`).

**If ND finds stable raw FAIL under A2 across form factors:** do **not** soft-pass; go to J5-β.

#### Fork **J5-β** — Renorm necessity (documented)

State as theory:

1. **Ledger mass** \(M_{\mathrm{ledger}}=E_\star/c^2\) sources **path cost** (\(M_{\mathrm{ray}}\)) via F1 multipole — closed.  
2. **Inertial mass** from free continuum is  
   \[
   m_{\mathrm{inertial}} = \xi\,\frac{U[\psi]}{c^2},
   \]
   with \(\xi=O(1)\) geometric (EM-mass analog; \(\xi=1\) or \(4/3\) class).  
3. **Naïve equality** \(m_{\mathrm{inertial}}=\int\rho_b/c^2\) **fails** when \(U\neq E_\star\) (form factor).  
4. **Closure options (mutually exclusive tags):**  
   - **β1 R1:** choose free constitutive \(s,\sigma_0\) so \(U=E_\star\) for the lock class under study (calibration of free constants — **not** a free-fit of \(m\) per object).  
   - **β2 Ledger-as-U:** set rest energy for inertia to free unlock \(E_\star^{\mathrm{(in)}}=U\); path-cost multipole may still use bound integral if dual bookkeeping is monist-justified (must convince TU).  
   - **β3 Bare + field mass:** total \(m = m_{\mathrm{bare}}+U/c^2\) with \(m_{\mathrm{bare}}\) from bound internal ledger — dual-looking; **high dualism risk** unless bare is same continuum.

**V77-3 bar (PROBLEM):**  
> J5 pass under non-tautology **or** documented renorm necessity + **kill of naïve** \(m=\int\rho_b\).

**TD Round-1 prediction:** **J5-β** is the working theory; **J5-α** remains aspirational. ND must:

- **Kill naïve** \(m_{\mathrm{inertial}}\stackrel{?}{=}\int\rho_b/c^2\) under A2 (expect FAIL at default constants).  
- **Pass** \(m_{\mathrm{inertial}}\stackrel{?}{=}U/c^2\) (boost vs field).  
- **Document** β1 renorm that restores three-way if free constants allowed to fix \(U=Mc^2\) once per constitutive family (not per particle).

### 3.6 Non-tautology criterion (export)

```text
tautology_flag = false   # mandatory for any J5 claim
protocol       = A2_boost_energy | A3_free_drag | ...
E_star_def     = integral_rho_b | unlock_work | free_U
m_inertial     = from ΔU(v) or F_ψ/a only
M_ray          = fixed G_eff, not free M fit
renorm_tag     = none | R1_s | ledger_as_U | bare_plus_field
```

Score:

\[
L_{\mathrm{triad}}
=
\left|\frac{m_{\mathrm{inertial}}}{E_{\mathrm{ref}}/c^2}-1\right|
+
\left|\frac{M_{\mathrm{ray}}}{E_{\mathrm{ref}}/c^2}-1\right|,
\]

with \(E_{\mathrm{ref}}\) tagged (\(E_\star\) raw vs \(U\) vs renormed).

---

## 4. Rods and clocks (C1 sketch for Dynamics)

Full reverse constraints: `v76/work/C/rods_clocks_C1_v0.md`. TD restates the **dynamics-facing** operational core.

### 4.1 Circularity problem

If rods/clocks are external stage instruments, “local \(c=\mathrm{const}\)” is dualist measurement.  
Monism requires: **rods and clocks are free-medium states** (or weakly bound free oscillators).

### 4.2 Operational definitions

| Object | Definition |
|--------|------------|
| **Local free clock** | Free oscillator \(O_U\) whose period \(\tau_U\) is a free-medium process |
| **Local free rod** | Segment calibrated so free round-trip time \(2L/c\) with isotropic free speed in a free patch |
| **Local \(c\)** | \(c\equiv 2L_{\mathrm{rod}}/\tau_{\mathrm{round}}\) after Einstein-style free-signal synchrony in the patch |

**Normalize** \(c\) to a constant **by construction** in each free local frame (Ax7 / NC-1.2). Warp = failure of a single global chart once locks rearrange free capacity.

### 4.3 Necessary conditions (carry forward)

| ID | Condition |
|----|-----------|
| **RC-C1** | Clock is free/weakly bound medium; kill free budget → rate change or death |
| **RC-C2** | Free-patch transport: identical free clocks agree (isotropy) |
| **RC-C3** | Rate = functional of free state, not independent metric field |
| **RC-R1** | Length = free-signal calibrated, not lattice \(\Delta x\) as ontology |
| **RC-R2** | Rotating rod in free patch leaves \(c\) unchanged |
| **RC-R3** | Nearby locks can change operational length vs distant standard |

### 4.4 Minimal dynamics claim (Round 1)

TD does **not** claim unique free-oscillator microphysics or full Lorentz group from free updates.

**Round-1 numeric bar (soft demo `D-DYN-rods`):**

1. Free-signal proxy path cost between two points (from \(\ell[\psi]\) or free hop).  
2. Define \(c_{\mathrm{op}}=2L_{\mathrm{op}}/\tau_{\mathrm{rt}}\) after calibration in free region.  
3. **Pass:** \(c_{\mathrm{op}}\) constant in free patches after calibration.  
4. **Near lock:** excess round-trip tracks path-cost integral (consistency with F1 rays).

Lattice \(\Delta x\) remains **storage** (scaffolding tag), not true monist length.

### 4.5 Link to inertia

\(m=E/c^2\) uses **this** operational \(c\). If \(c\) is a foreign chart speed, mass smuggles stage structure (C1 §6). J5 protocols should use the same \(c\) as free-wave / free-signal calibration when T2 is active.

---

## 5. Lorentz-from-\(c\) — caution (anti-overclaim)

### 5.1 What locality-\(c\) **does** force

| Forced / natural | Statement |
|------------------|-----------|
| Local causal cone | Free updates cannot outrun local free locality bound |
| Operational isotropy procedure | Free frames can be built so free speed is isotropic (C1) |
| Chart warp with fixed local \(c\) | Path cost \(\ell[\psi]\) around locks (v76 seed hinge) |

### 5.2 What locality-\(c\) does **not** automatically force

| Not forced by “\(c\) is free locality” alone | Why |
|-----------------------------------------------|-----|
| Full Lorentz group on matter fields | Needs transformation law for \(\rho_b,\psi\) and free oscillators |
| Minkowski metric as ontology | Dualist stage risk (S1) |
| Exact \(E=\gamma mc^2\) for locks | Needs hyperbolic free action + renorm (J5 open) |
| Thomas precession, relativistic beaming, full SR kinematics | Not Round-1 deliverables |
| GR / PPN light factor 4 | Residual R3; space+time path-cost split |

### 5.3 TD policy (Round 1)

1. **May claim:** free medium admits a **local causal structure** with isotropic free speed \(c\) after rod/clock calibration.  
2. **May use:** hyperbolic free capacity (T2) as **analog** of EM / relativistic field energy for inertia sketches.  
3. **Must not claim:** “therefore special relativity is derived” or “therefore Lorentz invariance of the monist world is proved.”  
4. **Kill overclaim (FOR_TU):** any demo that equates (local \(c\)) \(\Rightarrow\) (full SR) without transformation tests.

**Minimal positive next step (later rounds):** boost a free-wave packet and a lock; check whether free-signal rods/clocks reconstruct invariant interval **approximately** — killable, not assumed.

---

## 6. Cross-focus interfaces

| Focus | Shared | Dynamics-specific | Do not smuggle |
|-------|--------|-------------------|----------------|
| **EM (TE/NE)** | One continuum; free locality \(c\); free vs bound budget | \(\psi\) path-cost vs free gauge sibling | Maxwell ≠ \(\psi\); same \(c\) OK |
| **Matter (TM/NM)** | Locks source \(\rho_b\); dual source later | Inertia of locks; multi-lock forces not required R1 | Particle masses from SCP kernel |
| **Unification (TU)** | Seed primitives; G1–G3 | J5 status; Lorentz caution | Private dynamics ontology |

---

## 7. Kill conditions (Dynamics package fails if…)

| ID | Kill if… |
|----|----------|
| **KD1** | Time-dep free law does not recover F1-3D exterior \(\sim 1/r\) when static |
| **KD2** | Inertia only via \(a=F/m_{\mathrm{ledger}}\) (tautology) claimed as J5 pass |
| **KD3** | Naïve \(m=\int\rho_b/c^2\) asserted while A2 systematically measures \(m\sim U/c^2\neq\int\rho_b/c^2\) without renorm/ledger tag |
| **KD4** | Rods/clocks require non-medium external geometry as **ontology** (scaffolding OK if tagged) |
| **KD5** | Full Lorentz/SR claimed from locality-\(c\) slogan alone |
| **KD6** | Free waves require independent dualist gravity sector; free \(\psi\) idle |
| **KD7** | \(G_{\mathrm{eff}}\) or \(m_{\mathrm{inertial}}\) is a floating second-sector parameter with free constitutive law unused |

---

## 8. FOR_ND — Round-1 kill-gates (implement under `work/ND/`)

Numeric partner implements pass/fail tests. Prefer local Python; may copy/adapt `v76/work/B/sandbox_r4_inertia.py`. No `scp_sim`.

### Gate ND-G0 — Seed integrity (smoke)

| | |
|--|--|
| **Setup** | Static compact \(\rho_b\) (Gaussian OK); solve F1-3D (analytic Green or SOR). |
| **Pass** | Exterior multipole prefers \(1/r\); free deficit at core \(>0\); vacuum \(\rho_b=0\Rightarrow\psi=0\). |
| **Fail** | Static package broken — stop; do not interpret J5. |
| **Export** | `psi_multipole_prefer`, `free_deficit_core`, `vacuum_psi_max` |

### Gate ND-G1 — Naïve inertia kill (required for V77-3 path)

| | |
|--|--|
| **Setup** | Protocol **A2**: \(\rho_b(x-vt)\), \(v\ll c\); quasistatic or T1 re-solve \(\psi\); \(U=\frac{\sigma_0}{2}\int|\nabla\psi|^2\). Fit \(m_{\mathrm{inertial}}\) from \(\Delta U(v)\approx\tfrac12 m v^2\). Default free constants (e.g. \(s=\sigma_0=1\)). |
| **Pass (kill naïve)** | \(\lvert m_{\mathrm{inertial}}/(E_\star/c^2) - 1\rvert > \varepsilon_{\mathrm{kill}}\) (recommend \(\varepsilon_{\mathrm{kill}}=0.5\)) **and** \(\lvert m_{\mathrm{inertial}}/(U/c^2)-1\rvert < \varepsilon_{\mathrm{field}}\) (recommend \(0.3\)). |
| **Meaning** | Documents that naïve \(m=\int\rho_b/c^2\) **fails** as free-field inertia; inertia tracks \(U\). |
| **Fail gate** | If \(m_{\mathrm{inertial}}\approx E_\star/c^2\) **without** renorm at default constants, report **J5-α surprise** (theory must revise form-factor story). |
| **Tautology** | `tautology_flag` must be false (no \(F=m_{\mathrm{ledger}}a\)). |
| **Export** | `m_inertial`, `m_ledger`, `m_field`, `U`, `ratios`, `protocol=A2` |

### Gate ND-G2 — Renorm necessity (optional same run)

| | |
|--|--|
| **Setup** | From G1, compute \(s_*\) (or \(\sigma_{0*}\)) such that \(U(s_*)=E_\star\) for the **same lock shape family**. Re-eval triad at \(s_*\). |
| **Pass** | At single constitutive calibration \(s_*\), \(m_{\mathrm{inertial}}\approx E_\star/c^2\approx M_{\mathrm{ray}}\) within \(\varepsilon_{\mathrm{triad}}\) (e.g. 0.3). Tag `renorm_tag=R1_s`. |
| **Fail** | No single \(s_*\) closes boost+ray+ledger for **two** different lock widths — renorm is shape-by-shape fudge (dualism smell). |
| **Export** | `s_renorm`, `triad_after_renorm`, `shape_invariance_ok` |

### Gate ND-G3 — Time-dependent free response (first dynamics)

| | |
|--|--|
| **Setup** | Implement **T1** or **T2** in ≥1D; static source → recover F1; then (a) free pulse with \(\rho_b=0\), or (b) suddenly placed / moved lock. |
| **Pass T1** | \(\psi(t\to\infty)\to\) F1 solution (L2 error band); moving lock shows lag/wake. |
| **Pass T2** | Free pulse wavefront speed \(\approx c\) (within grid error); static limit F1. |
| **Fail** | Static limit wrong (KD1); or superluminal free update beyond numerical artifact. |
| **Export** | `law=T1|T2`, `static_recovery_err`, `wave_speed_meas` or `wake_lag` |

### Gate ND-G4 — Free-signal \(c_{\mathrm{op}}\) (C1 lite)

| | |
|--|--|
| **Setup** | Two free-region probes; round-trip free-signal time via path cost or T2 hops; calibrate \(L_{\mathrm{op}}\). |
| **Pass** | After calibration, \(c_{\mathrm{op}}\) constant in free patch (variation \(<\varepsilon_c\)); near-lock excess delay correlates with \(\int\gamma\psi\,dl\). |
| **Fail** | Only fixed \(\Delta x/dt\) “works” with no free-state dependence and no calibration story — dualist rods. |
| **Export** | `c_op_free`, `c_op_near_lock`, `delay_vs_pathcost_corr` |
| **Priority** | P2 if G1–G3 consume Round 1 |

### Gate ND-G5 — Lorentz overclaim control (meta)

| | |
|--|--|
| **Setup** | Documentation / result tags only. |
| **Pass** | No result file claims “SR derived” or “Lorentz group proved” from Round-1 gates. |
| **Fail** | Overclaim tags present → TU/O reject demo as V77-1 dynamics credit. |

### Priority order for Round 1

| Priority | Gate | Maps to |
|----------|------|---------|
| **P0** | ND-G0, ND-G1 | Seed + naïve kill (V77-3 path) |
| **P1** | ND-G2, ND-G3 | Renorm necessity + time-dep option |
| **P2** | ND-G4, ND-G5 | C1 lite + anti-overclaim |

### Suggested filenames (ND)

```text
work/ND/kill_gates_r1.md          # optional mirror of this §8
work/ND/sandbox_j5_a2.py          # adapt v76 sandbox_r4_inertia
work/ND/sandbox_timedep_T1T2.py   # relaxational / hyperbolic lite
work/ND/outputs/round1_*.json
```

---

## 9. Demo registry stubs (for TU)

| Demo ID | Theory | Numeric gate | Success | Residual |
|---------|--------|--------------|---------|----------|
| `D-DYN-J5` | §3 J5-β | ND-G1, ND-G2 | Naïve kill + renorm or raw triad | Form factor; ledger choice |
| `D-DYN-timedep` | §2 T1/T2 | ND-G3 | Static F1 recovery + wake or \(c\)-wave | Lattice \(c\) |
| `D-DYN-rods` | §4 C1 | ND-G4 | \(c_{\mathrm{op}}\) free-patch constancy | Oscillator microphysics open |
| `D-DYN-cwave` | §2 T2 | ND-G3 wave | Wavefront \(\approx c\) | Not full SR |

---

## 10. Residuals and non-goals (Round 1)

| Residual | Owner later |
|----------|-------------|
| Exact coefficient 1 without renorm (J5-α) | TD+ND deeper free action |
| Full Lorentz group / SR tests | TD+ND R2+ |
| PPN light factor 4 | Path-cost time+space split |
| Machian inertia from distant free medium | Speculative; not required |
| Multi-lock binding dynamics | TM/NM |
| Real \(G\), particle masses | Out of scope (PROBLEM §6) |

---

## 11. Bottom line

1. **Time-dep:** Prefer T1/T2 with F1 static recovery; T0 only as control for A2 energy.  
2. **Inertia:** Non-tautological \(m\) from free-field \(\Delta U(v)\) (A2). v76: tracks \(U\), not raw \(\int\rho_b\). **J5-β renorm necessity** is the Round-1 theory default; **kill naïve** \(m=\int\rho_b/c^2\).  
3. **Rods/clocks:** Free-signal calibration; lattice edges ≠ ontology.  
4. **Lorentz:** Local causal \(c\) ≠ derived SR — do not overclaim.  
5. **ND:** Execute gates G0–G3 (P0–P1); G4–G5 as capacity allows.

**V77-3 readiness:** G1 kill + G2 renorm doc (or unexpected G1 raw pass) ⇒ dynamics focus can claim closed-or-killed inertia residual.
