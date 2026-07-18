# Full Maxwell monist free-gauge package (v0)

**Agent:** TE (Theory — EM)  
**Date:** 2026-07-18  
**Round:** 2 (written) · **3 (FROZEN authoritative)**  
**Status:** **FROZEN AUTHORITATIVE** EM free-gauge law set for v77 C-A  
**Freeze stamp:** TE-007 (2026-07-19) — no new ontology forks (O-005)  
**Supersedes (scope):** Maxwell-lite draft `maxwell_monist_v0.md` (historical R1 subset; still valid for quasistatic DUAL-0)  
**Seed:** v76 F1-3D  
**Numeric twins:** NE R1 lite LIVE_PASS; NE R2 full Maxwell Yee LIVE_PASS (`full_maxwell_claim=true`)  
**Dual-channel:** NM DUAL-0 (ψ + Maxwell-lite Φ) V77-2 PASS; composition note `dual_channel_composition_r3.md`  
**Unification:** TU C-A channel; V77-2 dual-channel MET (parent O-004)

---

## 0. One-sentence claim

**Full classical Maxwell (M1–M4 + continuity) is the free linear gauge channel of the same continuum that hosts free-capacity \(\psi\): \((\mathbf{E},\mathbf{B})\) are free-medium stress; \(\rho_Q,\mathbf{J}_Q\) are lock–gauge ledgers; \(\varepsilon,\mu\) fix \(c=1/\sqrt{\varepsilon\mu}\) equal to free locality; \(\psi\) is a sibling free channel, never identified with \(\Phi_{\mathrm{EM}}\).**

---

## 1. Completeness + congruence (R3 freeze)

| Question | Answer |
|----------|--------|
| Is full Maxwell monist-**complete on theory side**? | **YES** — this doc FROZEN |
| Is full Maxwell monist-**complete as campaign EM demo**? | **YES** under labeled scope: 2D Yee TE+TM dynamics + 3D Coulomb recovery (NE R2) |
| `te_equation_match` (lite R1) | **YES** (TE-004) |
| `te_equation_match` (full R2) | **YES** (TE-007) — see gate map §8 |
| Dual-channel with \(\psi\) | **Composes** under WORLD: DUAL-0 joint LIVE (NM); full Yee is same C-A channel (composition residual: single-sandbox dynamical dual optional) |

**What “full Maxwell” means here:** all four macroscopic Maxwell equations in free vacuum constitutive form, plus charge continuity, plus energy/Poynting bookkeeping of free-gauge stress — **not** QED, materials microstructure, or real \(\alpha_{\mathrm{EM}}\).

### 1.1 Freeze rules (Round 3)

1. **This file is the single authoritative C-A equation source** for TU UNIFIED package.  
2. Do **not** fork new Maxwell ontologies; amend only via TE log + version bump if NE kill.  
3. Lite quasistatic M3 remains the **static sector** of the same package (§3.5), used by NM dual-channel.  
4. Yee 2D is an allowed discrete realization (§4.6), not a different theory.

---

## 2. Ontology (monist free gauge)

### 2.1 Shared WORLD primitives

Same as `maxwell_monist_v0` §1 and TU `WORLD_v0`:

| Primitive | Full-Maxwell reading |
|-----------|----------------------|
| One continuum | Free + bound only; no Minkowski stage ontology |
| Free field | Carries \(\mathbf{E},\mathbf{B}\) updates and free waves ≤ \(c\) |
| Bound / lock | Compact mass-form; may carry charge ledger \(Q\) and current \(\mathbf{J}_Q\) |
| Energy | Includes free-gauge stress \(u_{\mathrm{EM}}\) + bound \(E_\star\) |
| Mass | \(m=E_\star/c^2\) from \(\rho_b\), **not** from \(Q\) |
| \(c\) | Free locality \(=1/\sqrt{\varepsilon\mu}\) in free vacuum |
| Free channels | **C-ψ** path-cost; **C-A** free gauge (this package); **C-W** radiation packets of either |

### 2.2 What \(\mathbf{E},\mathbf{B}\) *are*

| Allowed | Forbidden |
|---------|-----------|
| Free continuum **stress / orientation** of free-gauge channel | Fields painted on independent empty space as second substance |
| Linear free DOFs with constitutive \(\varepsilon,\mu\) | EM sector that never depletes/participates in free ledger |
| Gauge potentials as free-channel redundancy | Primitive particle \(q\) with no lock ledger |

### 2.3 Charge and current as lock–gauge ledgers

\[
Q=\int\rho_Q\,dV,\qquad
\partial_t\rho_Q+\nabla\cdot\mathbf{J}_Q=0.
\tag{Cont}
\]

| | |
|--|--|
| \(\rho_Q\) | Signed **gauge-charge density** on continuum (typically supported on locks; TM Supp) |
| \(\mathbf{J}_Q\) | Gauge **current** — lock motion / bound flow: \(\mathbf{J}_Q=\rho_Q\mathbf{v}_{\mathrm{lock}}+\mathbf{J}_{\mathrm{int}}\) |
| Free charge mode | Optional later; R2 default: \(\mathrm{supp}(\|\rho_Q\|)\subseteq\mathrm{supp}(\rho_b)\) (TM DS-B lite) |

**Not mass:** \(\rho_b\ge 0\) sources C-ψ; \(\rho_Q\) sources C-A. Neutral lock: \(\rho_Q=0\), \(\rho_b>0\).

---

## 3. Complete free Maxwell set (continuum)

### 3.1 State

\[
\rho_Q(\mathbf{x},t),\;
\mathbf{J}_Q(\mathbf{x},t),\;
\mathbf{E}(\mathbf{x},t),\;
\mathbf{B}(\mathbf{x},t),\;
\varepsilon>0,\;\mu>0.
\]

Linear free constitutive (Maxwell vacuum monist default):

\[
\mathbf{D}=\varepsilon\mathbf{E},\qquad
\mathbf{H}=\mathbf{B}/\mu,
\qquad
c=\frac{1}{\sqrt{\varepsilon\mu}},\qquad
Z=\sqrt{\mu/\varepsilon}.
\tag{Const}
\]

### 3.2 The four Maxwell equations (boxed)

**Gauss for \(\mathbf{B}\) (no magnetic monopoles):**

\[
\boxed{\nabla\cdot\mathbf{B}=0}
\tag{M1}
\]

**Faraday:**

\[
\boxed{\nabla\times\mathbf{E}+\partial_t\mathbf{B}=\mathbf{0}}
\tag{M2}
\]

**Gauss for \(\mathbf{E}\) (lock–gauge source):**

\[
\boxed{\nabla\cdot(\varepsilon\mathbf{E})=\rho_Q}
\tag{M3}
\]

**Ampère–Maxwell:**

\[
\boxed{\nabla\times\bigl(\mathbf{B}/\mu\bigr)-\partial_t(\varepsilon\mathbf{E})=\mathbf{J}_Q}
\tag{M4}
\]

### 3.3 Continuity (necessary consistency)

Taking \(\partial_t\) of M3 and \(\nabla\cdot\) of M4 yields **(Cont)**.  
**Monist content:** charge ledger conservation is not an extra dualist law — it is free-channel consistency of C-A under lock currents.

If (Cont) is violated by prescribed \((\rho_Q,\mathbf{J}_Q)\), discrete Gauss will **drift** — NE must either enforce Cont or project.

### 3.4 Wave equations (derived, free vacuum)

Source-free \(\rho_Q=\mathbf{J}_Q=\mathbf{0}\), constant \(\varepsilon,\mu\):

\[
\bigl(\partial_t^2-c^2\nabla^2\bigr)\mathbf{E}=\mathbf{0},
\qquad
\bigl(\partial_t^2-c^2\nabla^2\bigr)\mathbf{B}=\mathbf{0},
\tag{Wave}
\]

with transversality \(\nabla\cdot\mathbf{E}=0\), \(\nabla\cdot\mathbf{B}=0\) preserved if initial data satisfy M1, M3.

**Shared-\(c\) hinge:** \(c\) here **is** free locality (same number as path-cost / Dynamics free waves).

### 3.5 Quasistatic Coulomb reduction (NE R1 subset)

\(\partial_t=0\), \(\mathbf{J}_Q=\mathbf{0}\), \(\mathbf{B}=\mathbf{0}\):

\[
-\nabla\cdot(\varepsilon\nabla\Phi)=\rho_Q,\quad\mathbf{E}=-\nabla\Phi.
\tag{Coulomb}
\]

Exterior: \(\Phi\sim Q/(4\pi\varepsilon r)\), \(E\sim Q/(4\pi\varepsilon r^2)\).  
**NE R1 LIVE_PASS** is congruent with this reduction + scalar wave proxy for \(c\).

### 3.6 Energy and Poynting (free ledger)

\[
u_{\mathrm{EM}}=\frac{\varepsilon}{2}|\mathbf{E}|^2+\frac{1}{2\mu}|\mathbf{B}|^2,
\qquad
\mathbf{S}=\mathbf{E}\times\mathbf{H}=\frac{1}{\mu}\mathbf{E}\times\mathbf{B}.
\tag{Poynting}
\]

Local balance (standard Maxwell identity):

\[
\partial_t u_{\mathrm{EM}}+\nabla\cdot\mathbf{S}=-\mathbf{J}_Q\cdot\mathbf{E}.
\tag{Work}
\]

**Monist placement:** \(u_{\mathrm{EM}}\) is **free continuum stress** of C-A; \(\mathbf{J}_Q\cdot\mathbf{E}\) is work on lock-gauge ledger (bound/kinetic of locks — Dynamics/Matter). Radiation = free-in-flight C-W packets of C-A.

**Residual R-EM1:** Coulomb self-energy vs \(Mc^2\) renorm (parallels J5) — not required for theory completeness of Maxwell form.

### 3.7 Potential form (optional, same physics)

\[
\mathbf{B}=\nabla\times\mathbf{A},\quad
\mathbf{E}=-\nabla\Phi-\partial_t\mathbf{A}.
\]

Gauge freedom is **redundancy of free-channel description**, not a second sector. Lorenz gauge \(\nabla\cdot\mathbf{A}+c^{-2}\partial_t\Phi=0\) yields wave equations for \(\Phi,\mathbf{A}\) with sources \(\rho_Q/\varepsilon\), \(\mu\mathbf{J}_Q\).

M1 is automatic if \(\mathbf{B}=\nabla\times\mathbf{A}\). Discrete: prefer Yee \(\mathbf{E},\mathbf{B}\) primary (below) so M1 is preserved by construction.

---

## 4. Discrete / Yee-friendly form (for NE)

Goal: NE can implement **structure-preserving** free Maxwell without inventing dualist stage ontology. Tags: `em_solver=free_maxwell_yee`, `E_origin=free_maxwell_full`, `sector=1`.

### 4.1 Staggered layout (standard Yee)

On Cartesian grid spacing \(h\), time step \(\Delta t\):

| Field | Location |
|-------|----------|
| \(E_x\) | edges \((i+\tfrac12,j,k)\) at integer times \(n\) |
| \(E_y,E_z\) | analogous edge centers |
| \(B_x\) | faces \((i,j+\tfrac12,k+\tfrac12)\) at half times \(n+\tfrac12\) |
| \(B_y,B_z\) | analogous |
| \(\rho_Q\) | cell centers (or vertices consistently) |
| \(\mathbf{J}_Q\) | colocated with \(\mathbf{E}\) components (edge currents) |

### 4.2 Discrete operators

Let \(\mathrm{curl}_h\), \(\mathrm{div}_h\) be standard staggered differences. Identities to preserve:

\[
\mathrm{div}_h(\mathrm{curl}_h\,\bullet)=0
\quad\text{(exact on closed grid / periodic; careful BC)}.
\tag{Id}
\]

### 4.3 Leapfrog update (vacuum Const)

\[
\begin{aligned}
\mathbf{B}^{n+1/2}
&=\mathbf{B}^{n-1/2}
-\Delta t\,\mathrm{curl}_h\,\mathbf{E}^{n},
\\
\mathbf{E}^{n+1}
&=\mathbf{E}^{n}
+\frac{\Delta t}{\varepsilon}
\Bigl(
\frac{1}{\mu}\mathrm{curl}_h\,\mathbf{B}^{n+1/2}
-\mathbf{J}_Q^{n+1/2}
\Bigr).
\end{aligned}
\tag{Yee}
\]

**CFL:** \(c\Delta t/h \le 1/\sqrt{d}\) (\(d=2\) or \(3\)). Prefer \(c=1/\sqrt{\varepsilon\mu}\) from constitutive, not a foreign chart speed.

### 4.4 Constraint preservation

| Constraint | Discrete status |
|------------|-----------------|
| M1 \(\mathrm{div}_h\mathbf{B}=0\) | **Preserved** if initial \(\mathrm{div}_h\mathbf{B}=0\) and (Id) holds |
| M3 \(\mathrm{div}_h(\varepsilon\mathbf{E})=\rho_Q\) | Preserved if (Cont) holds discretely and IC satisfy Gauss; else **project** or cleanse |
| Cont | Update \(\rho_Q\) by discrete continuity with edge \(\mathbf{J}_Q\), or co-evolve locks |

**Gauss cleanse (optional):** after update, \(\mathbf{E}\leftarrow\mathbf{E}-\nabla_h\chi\) with \(-\varepsilon\nabla_h^2\chi=\mathrm{div}_h(\varepsilon\mathbf{E})-\rho_Q\).

### 4.5 Static Coulomb from same continuum

For quasistatic tests, either:

1. Relax \(\mathbf{E}\) under M3 with \(\mathbf{B}=\mathbf{0}\), or  
2. Solve \(-\varepsilon\nabla_h^2\Phi=\rho_Q\), \(\mathbf{E}=-\nabla_h\Phi\) (NE R1 path).

Both are free-channel C-A; tag consistently. Full Maxwell claim requires path (Yee) or equivalent time-dep \(\mathbf{E},\mathbf{B}\).

### 4.6 2D TE/TM reduction (allowed R2)

For cheap Faraday + wave-with-B demos:

- **2D TM\(_z\):** \(E_z, B_x, B_y\) — wave + M1 in-plane  
- **2D TE\(_z\):** \(B_z, E_x, E_y\) — Faraday loop tests  

3D Coulomb multipole may stay on separate static solver (same tags) until joint 3D Yee affordable. **Honest:** mark `embedding_dim` in exports.

### 4.7 Dualist adversary (discrete)

Twin: hand Biot–Savart / Coulomb on fixed stage with free DOF idle (`sector=2`). Fit may match; monism needs free-origin tags + Occam (R-EM2).

---

## 5. Interface to \(\psi\) (sibling, not identical)

### 5.1 TE-IA1 (frozen)

> \(\psi\) and free-gauge \((\Phi,\mathbf{A})\) or \((\mathbf{E},\mathbf{B})\) are **independent free DOFs** of one medium, coupled through locks and shared \(c\)/budget — **not** by field identification.

### 5.2 Dual free system (static dual-channel; TM DUAL-0)

\[
\boxed{
\begin{aligned}
-\nabla\cdot\bigl(\sigma(\rho_f)\nabla\psi\bigr) &= s\,\rho_b,\\
\nabla\cdot(\varepsilon\mathbf{E}) &= \rho_Q,\\
\rho_f+\rho_b&=\rho_{\mathrm{tot}},\\
\mathrm{supp}(|\rho_Q|)&\subseteq\mathrm{supp}(\rho_b).
\end{aligned}
}
\tag{DUAL-0}
\]

**Time-dep dual-channel (R2 target sketch):**

\[
\begin{aligned}
&\text{C-A: full M1–M4 + Cont on }(\mathbf{E},\mathbf{B},\rho_Q,\mathbf{J}_Q),\\
&\text{C-ψ: F1-3D quasistatic or TD time-dep free response on }\psi,\\
&\mathbf{J}_Q \text{ from lock velocities (NM/TD)}.
\end{aligned}
\tag{DUAL-1}
\]

Same geometry, same \(c\), two constitutive tuples. Full joint note: [`dual_channel_joint_v0.md`](dual_channel_joint_v0.md).

### 5.3 Sibling stress tests (kill if broken)

| Case | \(\psi\) | \(\mathbf{E}\) |
|------|----------|----------------|
| Neutral mass \(\rho_b>0,\rho_Q=0\) | monopole \(1/r\) | **≈0** (static) |
| Opposite charges, equal \(|Q|\), both \(\rho_b>0\) | same-sign mass monopole | dipole / opposite \(\Phi\) |
| \(\rho_Q\neq 0\) free of all \(\rho_b\) | (forbidden by Supp default) | — |

---

## 6. Constitutive summary

| Param | Role |
|-------|------|
| \(\varepsilon\) | Gauss/Coulomb strength; \(\mathbf{D}=\varepsilon\mathbf{E}\) |
| \(\mu\) | Induction; with \(\varepsilon\) sets \(c\) |
| \(c=1/\sqrt{\varepsilon\mu}\) | Free locality = EM wave speed |
| \(\sigma,s,\gamma\) | C-ψ only (v76); **independent** of \(\varepsilon,\mu\) |
| \(Z=\sqrt{\mu/\varepsilon}\) | Impedance diagnostics |

Sandbox default: \(\varepsilon=\mu=1\Rightarrow c=1\). Off-unit: vary \(\varepsilon\) or \(\mu\) and demand \(c_{\mathrm{meas}}\) tracks (NE R1 already for scalar wave).

---

## 7. Kill conditions (full Maxwell monism)

### 7.1 Package kills (EM focus)

| ID | Kill if… |
|----|----------|
| **K-EM1** | Full Maxwell requires independent EM sector that cannot be free continuum stress |
| **K-EM2** | \(c_{\mathrm{wave}}\) cannot equal free locality from \(\varepsilon,\mu\) without foreign speed forever |
| **K-EM3** | Continuity / Gauss systematically fail under free dynamics when Cont is enforced |
| **K-EM4** | Free-medium M3 never yields \(E\sim 1/r^2\) (lite already **refuted** by NE R1) |
| **K-EM5** | Sibling collapse forced: neutral mass must source \(\mathbf{E}\), or opposite \(Q\) must reverse \(\psi\) monopole |
| **K-EM6** | Multi-focus unity: C-A contradicts budget identity with no repair (TU K2) |
| **K-EM7** | Occam: monist free gauge never preferred over dualist twin on any joint criterion |
| **K-EM8** | **Full-Maxwell-specific:** Faraday (M2) or M1 cannot be preserved as free-channel geometry without dualist metric stage as ontology |
| **K-EM9** | **Full-Maxwell-specific:** Ampère–Maxwell displacement current cannot live in free continuum (only magnetostatics on stage works) |
| **K-EM10** | Time-dep \(\mathbf{B}\) dynamics **require** second substance current foreign to lock ledger under all redesigns |

### 7.2 When full Maxwell *cannot* stay monist → program path

If after NE R2 + one redesign round:

1. Yee/free Maxwell fails K-EM8–10 with monist tags, **and**  
2. Only dualist stage-Maxwell (foreign \(J\), idle free budget) passes Faraday/divB/wave-B gates,

then TE proposes **focus kill** of full C-A dynamic Maxwell as monist (may **keep** quasistatic Coulomb lite as monist-eligible residual) and TU evaluates **K1** / partial retreat vs **V77-K**.

**Theory stance now:** no evidence of K-EM8–10; full Maxwell is monist-**eligible** and equation-complete. Numeric must try.

### 7.3 Not automatic kills

- Poisson isomorphism R-EM2 (static Coulomb twin).  
- R-EM1 self-energy renorm.  
- Using 2D TE/TM for Faraday tests while Coulomb stays 3D.  
- Gauge choice.  
- SI \(\alpha_{\mathrm{EM}}\).

---

## 8. Congruence map

### 8.1 NE R1 (lite) — **STAMPED congruent**

```text
TE Coulomb + Wave(c)     NE R1
  M3 + E=-∇Φ        ←→  continuum Gauss/Coulomb LIVE_PASS
  E ~ 1/r²          ←→  KG3 PASS
  vacuum            ←→  KG1 PASS
  c=1/√(εμ) wave    ←→  KG4 unit+off-unit PASS
```

`te_equation_match_r1 = yes` (TE-004).

### 8.2 NE R2 (full) — **STAMPED congruent** (TE-007)

```text
TE M1 div B=0           ←→  FM3 / KG-F1 PASS (divB ~ 1e-15)
TE M2 Faraday           ←→  FM4 / KG-F2 PASS (Yee identity)
TE M4 + Wave E,B        ←→  FM2 unit+off / KG-F3 PASS (v/c=1; ε=4→c=0.5)
TE Cont                 ←→  FM6 / KG-F4 PASS
TE M3 Coulomb           ←→  FM7+FM5 / KG-F5 PASS
TE full_maxwell_claim   ←→  TRUE (NE + TE endorse)
embedding_dim_dyn=2     ←→  §4.6 allowed; labeled honest
```

Detail: [`congruence_stamp_r3.md`](congruence_stamp_r3.md).  
Gates: [`for_ne_kill_gates_r2.md`](for_ne_kill_gates_r2.md).

### 8.3 Dual-channel with ψ — **composes** (TE-007)

```text
TE DUAL-0 + TE-IA1      ←→  NM R2 joint KG7/8/9 PASS (ψ + lite Φ)
TE full C-A             ←→  NE R2 Yee = dynamical sector of same C-A
COMP-1                  ←→  dual_channel_composition_r3.md
```

---

## 9. Demo IDs (update for TU)

| Demo ID | Scope | Theory | Numeric target | Status after TE R2 |
|---------|-------|--------|----------------|--------------------|
| **D-EM-maxwell-lite** | M3+Coulomb+scalar wave | v0 + this §3.5 | NE R1 | **theory PASS; numeric LIVE_PASS** |
| **D-EM-maxwell-full** | M1–M4+Cont+Yee | **this doc** | NE R2 | **theory DRAFT→READY; numeric PLANNED** |
| **D-EM-faraday** | M2 | §3.2 M2 | NE R2 | **PROPOSED** |
| **D-EM-divB** | M1 | §3.2 M1 | NE R2 | **PROPOSED** |
| **D-EM-ampere** | M4 | §3.2 M4 | NE R2 | **PROPOSED** |
| **D-EM-continuity** | Cont | §3.3 | NE R2 | **PROPOSED** |
| **D-EM-wave-EB** | Wave with B | §3.4 | NE R2 | **PROPOSED** |
| **D-EM-sibling-psi** | TE-IA1 | §5 | NE+NM | **V77-2 critical** |
| **D-DUAL-channel** | DUAL-0/1 | dual_channel_joint | NE+NM | **V77-2** |

---

## 10. Residuals

| ID | Residual |
|----|----------|
| R-EM1 | Coulomb self-energy vs mass ledger |
| R-EM2 | Static Poisson isomorphism |
| R-EM3 | Radiation reaction / lock back-reaction |
| R-EM4 | Micro origin of \(\varepsilon,\mu\) |
| R-EM5 | SCP kernel U(1) vocabulary bridge only |
| R-EM6 | Magnetostatics of bound currents / multipoles |
| R-EM7 | Full 3D Yee + dual-channel cost (may stage 2D first) |

---

## 11. Bottom line

**Full classical Maxwell is monist-complete on theory and (scoped) numeric sides:** M1–M4, Cont, Const, Wave, Poynting, lock ledgers, sibling \(\psi\), Yee form; NE R2 LIVE_PASS; `te_equation_match_full=yes` (TE-007).

**Package FROZEN** as authoritative C-A for v77. Dual-channel with \(\psi\) **composes** (DUAL-0 LIVE + full Yee same channel). Residual: optional single-box DUAL-1 (RC1), not ontology.

**TU may cite this file** for UNIFIED C-A; no TE ontology forks further unless gates reverse.
