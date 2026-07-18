# J5-β Default Theory — Locked (v0)

**Agent:** TD (Theory — Dynamics)  
**Round:** 2  
**Date:** 2026-07-18  
**Status:** **LOCKED as Dynamics default** for residual R2 / V77-3 path  
**Parent:** [`dynamics_package_v0.md`](dynamics_package_v0.md)  
**Numeric twin:** ND R1 — `work/ND/outputs/j5_formfactor_result.json`, `D-DYN-j5-formfactor`, `D-DYN-free-wave-c`  
**Orchestrator:** O-002 / O-003 — lock J5-β; dual-channel + full Maxwell progress  

---

## 0. Lock declaration

\[
\boxed{\text{J5-}\boldsymbol{\beta}\ \text{is the default Dynamics theory of inertia.}}
\]

| Claim | Status |
|-------|--------|
| Naïve \(m_{\mathrm{inertial}}=\displaystyle\int\rho_b/c^2\) at generic free constants | **KILLED** (ND-G1 PASS; v76 R4) |
| Free-field inertia \(m_{\mathrm{inertial}}=\xi\,U[\psi]/c^2\) | **DEFAULT** (tier_C structural + A2 protocol) |
| Path-cost / ray mass \(M_{\mathrm{ray}}=M_{\mathrm{ledger}}=\int\rho_b/c^2\) via F1 | **HELD** (tier_A; seed) |
| Raw three-way triad without renorm/ledger split | **NOT required** for V77-3 |
| J5-α (dynamics-free coefficient 1, all shapes, fixed free constants) | **DEFERRED** (not default; open research) |

**V77-3 reading under this lock:**  
> Documented renorm / ledger-split necessity + kill of naïve \(m=\int\rho_b\)  
> **is** the closed Dynamics residual — not a failure mode.

Parent may declare V77-3 **MET** once TU accepts this package + ND G1 kill + G2 policy below.  
Raw triad PASS is **bonus**, not the bar.

---

## 1. Why β (evidence, not preference)

### 1.1 Numeric facts (ND R1 / v76 R4 baseline)

Defaults: \(A=0.4\), \(\sigma=1\), \(s=\kappa=\sigma_0=1\), \(\gamma=0.5\), \(c=1\).

| Quantity | Value |
|----------|------:|
| \(M_{\mathrm{ledger}}=\int\rho_b/c^2\) | 6.299844 |
| \(U[\psi]\) free self-energy | 0.890953 |
| form factor \(ff=U/(Mc^2)\) | **0.141421** |
| \(m_{\mathrm{field}}=U/c^2\) | 0.890953 |
| \(m_{\mathrm{boost}}\) (\(\xi=1\)) | 0.890953 |
| \(M_{\mathrm{ray}}\) (F1 Born) | 6.299844 |
| \(s_*\) for \(U=Mc^2\) (this shape) | 7.071068 |
| rel split \(\lvert m_{\mathrm{field}}-M\rvert/M\) | 0.8586 |
| `tautology_flag` | **false** |

Gates: tier_A **PASS**; tier_B raw **FAIL**; tier_C **PASS**; naïve kill **PASS**.

### 1.2 Form-factor law (Gaussian free capacity)

For F1 vacuum \(-\sigma_0\nabla^2\psi=s\rho_b\), \(\kappa=s/\sigma_0\), infinite-space Gaussian \(\rho_b\):

\[
\boxed{
ff
\equiv
\frac{U}{M c^2}
=
\frac{s\,\kappa\, M}{8\pi\,\sigma\,\sqrt{\pi}}
=
\frac{s^2 M}{8\pi\sigma_0\,\sigma\sqrt{\pi}}
}
\tag{FF}
\]

ND scan (fixed \(s=\kappa=1\)): \(ff = 0.141421356237\times(A/0.4)\times\sigma^2\) — analytic identity verified on 18 points.

**Theory content of (FF):**  
Free self-energy scales as \(M^2/R_{\mathrm{eff}}\) (here \(R_{\mathrm{eff}}\propto\sigma\)).  
Bound ledger scales as \(M\).  
Ratio \(ff\propto M/\sigma\) is **shape- and constitutive-dependent** — not a bug, a structural EM-mass-class fact.

### 1.3 Universal \(s_*\) fails (ND-G2)

| Fact | Implication |
|------|-------------|
| \(s_*=1/ff\) ranges ~1.18–39 across \((A,\sigma)\) scan | No **single** free-constant \(s\) closes \(U=Mc^2\) for all lock shapes |
| Per-shape \(s_*\) always works algebraically | Calibration of a **lock class**, not a universal medium law |
| Accidental raw PASS (e.g. \(A=0.6,\sigma=2\), \(ff\approx0.85\)) | Parameter tuning ≈ hidden renorm — **not** J5-α proof |

**Therefore:** pure “universal R1 renorm of free vacuum constants” is **not** the default closure. Shape-class calibration remains allowed as a **demo tag**, not as a theorem of coefficient 1.

---

## 2. Locked ontology: two bookkeepings, one continuum

J5 residual is **not** dualism. It is a split between two **monist ledgers** of the same fabric:

| Ledger | Symbol | What it counts | What it sources |
|--------|--------|----------------|-----------------|
| **Bound occupancy** | \(E_\star^{\mathrm{path}}=\int\rho_b\,dV\) | How much medium is locked (mass-form budget) | F1: \(\psi\) multipole, path cost, \(M_{\mathrm{ray}}\) |
| **Free capacity stress** | \(E_\star^{\mathrm{in}}=U[\psi]=\frac{\sigma_0}{2}\int\lvert\nabla\psi\rvert^2\) | Free continuum energy of the capacity channel around the lock | Inertial response \(m_{\mathrm{inertial}}=\xi E_\star^{\mathrm{in}}/c^2\) |

\[
\boxed{
\begin{aligned}
M_{\mathrm{ray}}
&=
\frac{E_\star^{\mathrm{path}}}{c^2}
\quad\text{(F1 + fixed }G_{\mathrm{eff}}\text{)},
\\
m_{\mathrm{inertial}}
&=
\xi\,\frac{E_\star^{\mathrm{in}}}{c^2}
=
\xi\,\frac{U[\psi]}{c^2}
\quad\text{(A2 free-field; }\xi=1\text{ default)}.
\end{aligned}
}
\tag{J5-β}
\]

**Naïve collapse** \(E_\star^{\mathrm{path}}\equiv E_\star^{\mathrm{in}}\) (hence \(m=M_{\mathrm{ray}}\) from free dynamics alone) is **false** at generic \((\sigma,s,\text{shape})\).

### 2.1 Sub-options under β (priority)

| ID | Rule | When to use | Status under lock |
|----|------|-------------|-------------------|
| **β-ledger (preferred)** | Always report \(m_{\mathrm{inertial}}\) vs \(U/c^2\) and \(M_{\mathrm{ray}}\) vs \(\int\rho_b/c^2\) **separately**; equality only if \(ff\to 1\) | Default theory language; V77-3 | **LOCKED DEFAULT** |
| **β-class** | For a **named lock family**, choose \(s_*\) (or \(\sigma_{0*}\)) once so \(U=E_\star^{\mathrm{path}}\) for that family; tag `renorm_tag=R1_class` | Demo convenience; single-shape textbooks | **Allowed, not universal** |
| **β-U-rest** | Call monist rest energy for **inertia** \(E_\star^{\mathrm{in}}=U\); path multipole still uses \(\int\rho_b\) | Cleanest non-tautological inertia story | **Equivalent to β-ledger** |
| **β-bare+field** | \(m=m_{\mathrm{bare}}+U/c^2\) with independent bare mass | Avoid unless bare is same-continuum internal lock stress with equation | **Discouraged** (dualism smell) |
| **J5-α** | Prove \(ff=1\) from free action + formation | Future free-formation research | **Not default** |

### 2.2 Forbidden under the lock

| Forbidden | Why |
|-----------|-----|
| Claim J5 raw three-way PASS without \(ff\approx 1\) or tagged renorm | Dishonest |
| \(a=F/m_{\mathrm{ledger}}\) as inertia proof | Tautology |
| Universal free vacuum \(s\) that fits all shapes silently | ND-G2 FAIL |
| Foreign inertial mass unrelated to free or bound continuum ledgers | Dualist; V77-K risk |
| Collapsing path-cost \(\psi\) into Coulomb \(\Phi\) to “fix” energy split | Sibling kill (TE K-ID) |

### 2.3 Export tags (mandatory for ND/TU)

```text
j5_theory_default     = beta
E_star_path           = integral_rho_b
E_star_inertial       = U_free_capacity   # or unlock_work when available
m_inertial_protocol   = A2_boost_energy
m_ray_protocol        = F1_fixed_Geff
form_factor           = U / (E_star_path)
renorm_tag            = none | R1_class | (not R1_universal)
tautology_flag        = false
naive_integral_killed = true
```

---

## 3. Time-dependent free response in a dual-channel world

Round-2 theme: **one medium**, two free channels (V77-2), dynamics time-dep if free.

### 3.1 Dual-channel state (shared with TE/TM)

\[
\bigl(
\rho_b,\;\rho_f,\;\psi,\;
\rho_Q,\;\mathbf{J}_Q,\;\mathbf{E},\;\mathbf{B}
\bigr)
\quad\text{on one continuum}.
\]

| Channel | Quasistatic | Time-dependent (Dynamics + EM) |
|---------|-------------|-------------------------------|
| **C-ψ** free capacity | F1-3D: \(-\nabla\cdot(\sigma\nabla\psi)=s\rho_b\) | T1 relaxational and/or T2 hyperbolic (below) |
| **C-A** free gauge | Coulomb: \(\nabla\cdot(\varepsilon\mathbf{E})=\rho_Q\) | Full Maxwell M1–M4 (TE R2 target) |

**Hard sibling rules (unchanged):**

1. \(\psi\not\equiv\Phi_{\mathrm{EM}}\) as ontology.  
2. \(\rho_b\not\equiv\rho_Q\) (TM dual-source may place both on same lock support).  
3. Shared free locality \(c\): \(c_{\psi\text{-wave}}=c_{\mathrm{EM}}=1/\sqrt{\varepsilon\mu}\) (JC1).  
4. Budget: \(\rho_f+\rho_b=\rho_{\mathrm{tot}}\); free energies \(U_\psi+U_{\mathrm{EM}}\) live in free ledger.

### 3.2 Capacity-channel time laws (TD owns)

**T1 — Relaxational free capacity (wake / approach to F1):**

\[
\partial_t\psi
=
\kappa_\psi\,\nabla\cdot\bigl(\sigma\nabla\psi\bigr)
+
s\,\rho_b.
\tag{T1}
\]

- Static fixed point: F1-3D (TD-S).  
- Moving lock \(\rho_b(\mathbf{x}-\mathbf{v}t)\): lag/wake at finite \(\kappa_\psi\).  
- No sharp light cone for \(\psi\) alone — free radiation for rods/clocks and EM is **C-A** (and optional T2).

**T2 — Hyperbolic free capacity (capacity free waves):**

\[
\frac{1}{c^2}\partial_{tt}\psi
=
\nabla\cdot\bigl((\sigma/\sigma_0)\nabla\psi\bigr)
+
\frac{s}{\sigma_0}\rho_b.
\tag{T2}
\]

- Static: F1.  
- Source-free: \(\partial_{tt}\psi=c^2\nabla^2\psi\) (ND `D-DYN-free-wave-c` **PASS**).  
- Same \(c\) as Maxwell vacuum waves under JC1.

**T0 — Quasistatic comoving (control only):**  
Re-solve F1 for \(\rho_b(t)\) each step — valid for A2 \(\Delta U(v)\) when \(v\ll c\) and radiation neglected; **not** dual-channel dynamics claim alone.

### 3.3 Dual-channel dynamics package (joint world)

Minimal closed system for V77-2 + Dynamics:

```text
(1) Bound / charge placement (TM locks):  ρ_b(x,t), ρ_Q(x,t), J_Q
(2) Capacity channel:                    T1 or T2 for ψ
(3) Gauge channel:                       full Maxwell M1–M4 for (E,B)
(4) Shared c:                            c = 1/√(εμ) = T2 wave speed
(5) Forces on locks (virtual work):
      F_ψ = ∫ (−ρ_b ∇ψ) dV
      F_EM = ∫ (ρ_Q E + J_Q × B) dV     # Lorentz; monist free-gauge stress
(6) Path cost (slow):                    ℓ = ℓ0 + γ ψ
(7) Inertia (J5-β):                      m_inertial from U_ψ (and U_EM if charged moving)
```

**Inertia note in dual channel:**  
Charged moving locks also store free-gauge field energy \(U_{\mathrm{EM}}\).  
Then

\[
m_{\mathrm{inertial}}^{\mathrm{(tot)}}
=
\xi_\psi\frac{U_\psi}{c^2}
+
\xi_{\mathrm{EM}}\frac{U_{\mathrm{EM}}}{c^2}
+
\cdots
\]

Same β logic: free stress energies, not bare \(\int\rho_b\). Neutral locks (\(Q=0\)) reduce to pure \(U_\psi\) story (ND R1).

### 3.4 What “time-dep free response” must recover

| Limit | Must recover |
|-------|----------------|
| \(\partial_t\to 0\), static \(\rho_b\) | F1-3D exterior \(\psi\sim 1/r\) |
| \(\partial_t\to 0\), static \(\rho_Q\), \(\mathbf{J}_Q=0\) | Coulomb / Gauss (TE) |
| \(\rho_b=\rho_Q=0\), free waves | \(v=c\) on **both** channels (shared \(c\)) |
| Slow boost, \(Q=0\) | J5-β: \(\Delta U_\psi\sim\tfrac12(U_\psi/c^2)v^2\) |

**Kill (dual-channel dynamics):** static limits wrong; \(c_{\psi}\neq c_{\mathrm{EM}}\) without frame story; forces from foreign \(G,m\) inserted by hand.

---

## 4. Interface to **full** Maxwell free waves (TE R2 / NE)

Round-2 human goal includes **full Maxwell**, not only lite quasistatic + wave equation.

### 4.1 What TD needs from full Maxwell

| TE/NE full-Maxwell piece | Dynamics use |
|--------------------------|--------------|
| M1 \(\nabla\cdot\mathbf{B}=0\) | Free-gauge identity; no monopole dualism |
| M2 Faraday \(\nabla\times\mathbf{E}+\partial_t\mathbf{B}=0\) | Induction; time-dep \(\mathbf{B}\) for rods/clocks & radiation |
| M3 Gauss \(\nabla\cdot(\varepsilon\mathbf{E})=\rho_Q\) | Dual-source Coulomb; sibling to F1 (not identical) |
| M4 Ampère–Maxwell \(\nabla\times\mathbf{B}/\mu-\partial_t(\varepsilon\mathbf{E})=\mathbf{J}_Q\) | Displacement current; free EM waves; continuity |
| Cont \(\partial_t\rho_Q+\nabla\cdot\mathbf{J}_Q=0\) | Ledger consistency for charge |
| Vacuum waves \(v=1/\sqrt{\varepsilon\mu}\) | JC1 shared-\(c\) with T2 / free packets |
| Poynting / \(u_{\mathrm{EM}}\) | Free-ledger radiation; EM contribution to inertia (β) |

### 4.2 Sibling vs identity (interface table)

| Item | Capacity (TD) | Gauge (TE full Maxwell) | Joint rule |
|------|---------------|-------------------------|------------|
| DOF | \(\psi\) | \(\mathbf{E},\mathbf{B}\) (or \(A,\Phi\)) | Independent free channels |
| Source | \(\rho_b\ge 0\) | \(\rho_Q\in\mathbb{R}\), \(\mathbf{J}_Q\) | Dual-source locks OK; no forced identity |
| Static Green | \(1/r\) for \(\psi\) | \(1/r\) for \(\Phi\), \(1/r^2\) for \(E\) | Same multipole *shape class*, different ontology |
| Waves | T2: \(\partial_{tt}\psi=c^2\nabla^2\psi\) | \(\partial_{tt}\mathbf{E}=c^2\nabla^2\mathbf{E}\) (transverse) | **Same \(c\)** (JC1) |
| Energy | \(U_\psi=\frac{\sigma_0}{2}\int\lvert\nabla\psi\rvert^2\) (+ kinetic if T2) | \(U_{\mathrm{EM}}=\int(\varepsilon E^2/2+B^2/(2\mu))\) | Both free ledger; sum for charged inertia |
| Path cost | \(\ell=\ell_0+\gamma\psi\) | **Not** from \(\Phi\) | Gravity-class optics stays on \(\psi\) |
| Inertia default | J5-β on \(U_\psi\) | EM-mass analog on \(U_{\mathrm{EM}}\) | Same β logic; no naïve \(\int\rho_b\) |

### 4.3 Coupling that is **allowed** vs **forbidden**

**Allowed (weak dual-channel):**

- Same grid, same \(c\), simultaneous evolution of T1/T2 and M1–M4.  
- Locks with \((\rho_b,\rho_Q)\) support overlap (TM).  
- Forces \(F_\psi+F_{\mathrm{EM}}\) on lock COM (virtual work).  
- Optional mild constitutive \(\sigma=\sigma(\rho_f)\), \(\varepsilon=\varepsilon(\rho_f)\) **without** replacing \(\psi\) by \(\varepsilon(\rho_b)\) gravity (GRIN-kill).

**Forbidden:**

- Setting \(\psi:=\Phi\) or \(\gamma\psi\) as Coulomb potential.  
- Using Maxwell only as dualist stage fields with free \(\psi\) idle (HK1-class).  
- Claiming full SR from shared \(c\) alone (dynamics_package §5).  
- Absorbing \(G_{\mathrm{eff}}\) into \(\varepsilon\) to force monism-by-fit.

### 4.4 Radiation split (C-W channel)

| Radiation type | Carrier | Numeric owner |
|----------------|---------|---------------|
| Capacity free waves | T2 \(\psi\) packets | ND |
| EM free waves | Full Maxwell \(\mathbf{E},\mathbf{B}\) | NE |
| Shared-\(c\) cross-check | Compare \(v_\psi\) vs \(v_E\) | ND+NE joint (V77-2) |

TD theory: both are **free continuum in flight**; neither is a second substance.

### 4.5 Rods/clocks note (C1 residual)

Full Maxwell free waves supply the natural free-signal for C1 rods (round-trip light).  
Capacity T2 waves are a valid **sibling** free-signal for the same calibration if JC1 holds.  
Round-2: prefer **EM free signals for C1** once NE full leapfrog exists; ND free-wave-\(c\) remains capacity-side control.

---

## 5. V77-3 closure criteria (under locked β)

| Criterion | Evidence | Owner |
|-----------|----------|-------|
| Naïve \(m=\int\rho_b/c^2\) killed non-tautologically | ND-G1 PASS | ND (done R1) |
| Form-factor / ledger-split theory written | This doc + (FF) | TD (**this lock**) |
| Renorm policy: no false universal \(s_*\) | ND-G2 FAIL documented; β-class allowed with tag | TD+ND |
| Free wave \(v=c\) on capacity channel | ND G3 PASS | ND (done R1) |
| No SR overclaim | G5 | TD+ND |
| Optional: dual-channel shared \(c\) with full Maxwell | Joint NE+ND | R2+ |

**Not required for V77-3:** J5-α; universal \(s_*\); full Lorentz group; PPN-4.

---

## 6. FOR_ND — Round 2 gates

Partner implements under `work/ND/`. Prefer local Python; reuse R1 sandboxes. No `scp_sim`.

### ND-R2-G0 — Accept J5-β in exports (docs + tags)

| | |
|--|--|
| **Action** | Update result JSON / README: `j5_theory_default=beta`; `naive_integral_killed=true`; report \(m_{\mathrm{inertial}}\) vs \(U/c^2\) and \(M_{\mathrm{ray}}\) vs \(\int\rho_b/c^2\) as **two** equalities, not one collapsed triad without \(ff\). |
| **Pass** | Tags present; no claim of raw three-way PASS at baseline. |
| **Fail** | Re-assert naïve equality or universal \(s_*\) as theorem. |

### ND-R2-G1 — Form-factor identity (theory–numeric lock)

| | |
|--|--|
| **Action** | Keep/export (FF) scan: \(ff=s\kappa M/(8\pi\sigma\sqrt{\pi})\); already largely done R1. |
| **Pass** | Analytic vs numeric agreement to \(\sim 10^{-6}\) relative on ≥3 shapes. |
| **Priority** | P1 if R1 outputs already suffice — may **cite R1** without re-run. |

### ND-R2-G2 — T1 or T2 with **sourced** F1 recovery (time-dep free response)

| | |
|--|--|
| **Setup** | Compact \(\rho_b\); evolve T1 or T2 from cold \(\psi=0\); measure approach to static F1 solution. Optional: translate lock, measure wake lag. |
| **Pass** | Static recovery L2 error within band (e.g. \(<5\%\) energy or residual); exterior multipole still \(\sim 1/r\) class if 3D/1D-analog documented. |
| **Fail** | End state not F1; or superluminal \(\psi\) update beyond scheme artifact (T2 CFL). |
| **Priority** | **P0** Round 2 Dynamics numeric for dual-channel readiness. |

### ND-R2-G3 — Shared-\(c\) with NE full Maxwell (joint V77-2)

| | |
|--|--|
| **Setup** | Same \(c\) units as NE Maxwell vacuum wave. Capacity free wave (existing sandbox) vs NE \(v_E\) (or Faraday/Ampère leapfrog). |
| **Pass** | \(\lvert v_\psi/c-1\rvert\) and \(\lvert v_E/c-1\rvert\) both in band; report joint table. |
| **Fail** | Systematic \(c_{\psi}\neq c_{\mathrm{EM}}\) with both claiming free vacuum. |
| **Priority** | **P0** when NE full Maxwell sandbox exists; else **defer** and tag `shared_c=pending_NE`. |
| **Coordinate** | FOR_NE — same \(c=1\) (or common off-unit \(c=0.5\) control). |

### ND-R2-G4 — Charged inertia sketch (optional)

| | |
|--|--|
| **Setup** | Neutral vs charged compact lock; quasistatic \(U_\psi+U_{\mathrm{EM}}\); compare \(m\sim (U_\psi+U_{\mathrm{EM}})/c^2\) story (analytic or coarse). |
| **Pass** | Neutral reduces to R1; charged \(m\) increases with \(U_{\mathrm{EM}}\) (directionally). |
| **Priority** | P2 — needs NE Coulomb energy export. |

### ND-R2-G5 — Anti-overclaim (meta)

| | |
|--|--|
| **Pass** | No “SR derived,” no “universal s* proves m=E/c²,” no \(\psi\equiv\Phi\). |
| **Priority** | Standing. |

### Round-2 priority stack

| Priority | Gate | Goal |
|----------|------|------|
| **P0** | R2-G0, R2-G2 | Lock β in data; time-dep F1 recovery |
| **P0*** | R2-G3 | Dual-channel shared \(c\) (depends on NE full Maxwell) |
| **P1** | R2-G1 | Cite form-factor law |
| **P2** | R2-G4, G5 | Charged inertia lite; meta |

### Suggested ND files

```text
work/ND/sandbox_timedep_F1.py      # T1/T2 + static recovery
work/ND/sandbox_shared_c_joint.py  # optional; or joint note with NE outputs
work/ND/outputs/round2_*.json
```

---

## 7. FOR_TE / FOR_NE / FOR_TU / FOR_TM (theory cross)

**FOR_TE:** TD locks J5-β (free stress ≠ bound integral generically). Full Maxwell free energy \(U_{\mathrm{EM}}\) is the same β-class object as \(U_\psi\). Interface §4 is TD’s dual-channel demand: M1–M4 + shared \(c\), no \(\psi\equiv\Phi\). Capacity T2 is **sibling** free wave, not Maxwell substitute.

**FOR_NE:** When full leapfrog Maxwell is up, export \(v_E/c\) and vacuum \(U_{\mathrm{EM}}\) for ND R2-G3/G4. Static Coulomb alone insufficient for joint radiation-\(c\) lock.

**FOR_TU:**  
- Set Dynamics inertia residual to **J5-β LOCKED**.  
- V77-3 path = naïve kill (done) + this theory lock + no universal-\(s_*\) lie.  
- Register `D-DYN-j5-formfactor` as LIVE under β (not raw triad).  
- Dual-channel V77-2 still needs joint ψ+Maxwell numeric (R2-G3).

**FOR_TM:** Dual-source locks feed \(\rho_b\to\psi\), \(\rho_Q\to\) Maxwell. Inertia of multi-lock systems sums free-channel field energies (β), not \(\sum\int\rho_b\) alone. Force closure \(F_\psi+F_{\mathrm{EM}}\) matches §3.3.

**FOR_O:** J5-β locked as requested. Time-dep free response specified for dual-channel world; Maxwell interface for full M1–M4 written. Blocking numeric: ND R2-G2; joint shared-\(c\) with NE full Maxwell for V77-2.

---

## 8. Bottom line

1. **J5-β is locked:** free inertial mass tracks free capacity (and EM) **stress energy** \(U/c^2\); path-cost mass tracks **bound ledger** \(\int\rho_b/c^2\); naïve collapse killed; universal \(s_*\) rejected.  
2. **Form factor (FF)** is the quantitative residual — theory + ND scan.  
3. **Dual-channel time-dep:** T1/T2 for \(\psi\) + full Maxwell for \((\mathbf{E},\mathbf{B})\); shared \(c\); forces from both free stresses.  
4. **Full Maxwell interface:** M1–M4, continuity, Poynting, same free-ledger β logic for \(U_{\mathrm{EM}}\).  
5. **ND R2:** tag β; sourced T1/T2 F1 recovery; shared-\(c\) joint when NE ready.

**Parent package path:**  
[`dynamics_package_v0.md`](dynamics_package_v0.md) (R1 options) + **this lock** (R2 default) → Dynamics theory ready for V77-3 under β.
