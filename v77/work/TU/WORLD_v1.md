# WORLD v1 — Cohesive Monist Free-Channel World

**Agent:** TU · **Round:** 2 · **Status:** coherent draft (not V77-4 / not terminal done)  
**Supersedes:** `WORLD_v0.md` (kept as history)  
**Inputs:** v76 F1-3D seed; TE `maxwell_monist_v0` + constitutive table; TD `dynamics_package_v0`; TM `lock_ontology_v0`; NE/ND/NM R1 PASS outputs; O-002/O-003  
**Human terminal bar:** cohesive unified theory **and** model **or** NULL/DISPROVE; **full Maxwell** required before campaign done (O-003).

---

## 0. One-page claim

There is **one continuum**. Its budget splits into **free** and **bound** (locks).  
**Mass** is bound ledger over locality: \(m = E_\star/c^2\). **Energy** is the continuum ledger (bound + free stress, including free-gauge \(u_{\mathrm{EM}}\) and free-capacity \(U[\psi]\)). **\(c\)** is free-field locality — one number in free vacuum for all free channels.

Locks **source free sibling channels** with **independent ledgers and constitutive laws**:

| Channel | Free DOF | Source | Constitutive | Long-range role |
|---------|----------|--------|--------------|-----------------|
| **C-ψ** free capacity | \(\psi,\sigma\) | \(\rho_b\ge 0\) | \(-\nabla\cdot(\sigma\nabla\psi)=s\rho_b\); \(\ell=\ell_0+\gamma\psi\) | Path cost / gravity-class exterior \(\sim 1/r\) |
| **C-A** free gauge | \(\mathbf{E},\mathbf{B}\) (or \(\Phi,\mathbf{A}\)) | \(\rho_Q\in\mathbb{R}\), \(\mathbf{J}_Q\) | Maxwell (target full M1–M4+Cont); \(c=1/\sqrt{\varepsilon\mu}\) | Coulomb \(E\sim 1/r^2\); radiation at \(c\) |
| **C-W** free waves | hyperbolic free packets of C-ψ and/or C-A | lock motion / multipoles | wave ops at same \(c\) | Causality, rods/clocks feed |

**Forbidden:** \(\psi\equiv\Phi_{\mathrm{EM}}\); local free-density GRIN as monist gravity; hand Poisson as monist *proof*; foreign dualist stage; \(c_{\mathrm{EM}}\neq c_{\mathrm{locality}}\) without medium-frame story.

**Warp:** what constant *local* free \(c\) looks like around locks in a global chart (path cost from \(\psi\), not from electrostatic \(\Phi\)).

---

## 1. Shared primitives (frozen)

| Primitive | Definition | Owner residual |
|-----------|------------|----------------|
| Continuum | Only medium; free + bound states | all |
| Free / bound budget | \(\rho_f+\rho_b=\rho_{\mathrm{tot}}\) (or integral) | seed LIVE |
| \(c\) | Free locality = free-signal bound | JC1 shared |
| Mass / ledger | \(M=E_\star/c^2\) path-cost sources \(\rho_b\) | TD: inertia ≠ raw \(\int\rho_b\) |
| Lock | Compact free-depleting mass-form; may carry \(q\) | TM S0–S3 |
| Free channel | Linear free DOF + constitutive on same continuum | TE/TD |
| Energy ledger | Bound unlock + free-channel stresses | TD J5-β; TE \(u_{\mathrm{EM}}\) |

---

## 2. Dynamics of free channels (laws stack)

### 2.1 Path-cost (v76 seed — LIVE)

\[
-\nabla\cdot(\sigma\nabla\psi)=s\rho_b,\quad
\ell=\ell_0+\gamma\psi,\quad
\psi\sim\frac{s E_\star}{4\pi\sigma_0 r}\ \text{(3D exterior)}.
\]

Time-dep options (TD): T1 relaxational / T2 hyperbolic free capacity; static recovery **TD-S** → F1-3D. Free hyperbolic medium: \(\partial_t^2\psi=c^2\nabla^2\psi\) (ND **PASS**).

### 2.2 Maxwell (terminal target = **full**)

**Written (TE):** M1 \(\nabla\cdot B=0\); M2 Faraday; M3 Gauss \(E\); M4 Ampère–Maxwell; Cont; Wave; Coulomb quasistatic.

**Numeric today (NE R1):** Maxwell-**lite** only — M3 quasistatic Coulomb + 1D free wave at \(c=1/\sqrt{\varepsilon\mu}\) (unit + off-unit).  
`full_maxwell_claim = false`.

**R2+ requirement (O-003 / human):** time-dependent \(\mathbf{E},\mathbf{B}\); Faraday + Ampère–Maxwell dynamics; \(\nabla\cdot B=0\); continuity; wave+Coulomb consistency in ≥2D — before terminal **done**.

### 2.3 Inertia (J5 — Dynamics) — **LOCKED J5-β** (TD R2)

Source: TD `j5_beta_default_v0.md` + ND R1 `j5_formfactor_result.json`.

| Claim | Status |
|-------|--------|
| \(M_{\mathrm{ray}}=M_{\mathrm{ledger}}=\int\rho_b/c^2\) via F1 | **HELD** |
| \(m_{\mathrm{inertial}}=\xi U[\psi]/c^2\) | **DEFAULT** |
| Naïve \(m=\int\rho_b/c^2\) | **KILLED** (ND-G1) |
| Universal \(s^*\) all shapes | **REJECTED** (ND-G2) |
| Form factor (Gaussian default) | \(ff\approx 0.141421\) |
| Preferred ledger split | \(E_\star^{\mathrm{path}}=\int\rho_b\); \(E_\star^{\mathrm{in}}=U\) (β-ledger) |
| J5-α raw triad | **DEFERRED** (not required for V77-3) |

**TU acceptance:** J5-β lock + ND naïve kill **satisfies V77-3 bar** as written in PROBLEM.md (renorm necessity + kill naïve), not raw three-way.

### 2.4 Locks and dual source (Matter)

Lock \(L_i=(X_i,E_{\star,i},q_i,\mathcal{C}_i)\).  
**DUAL-0 (R1 LIVE):** \(\rho_b\to\psi\), \(\rho_Q\to\) gauge lite, \(\mathrm{supp}(|\rho_Q|)\subseteq\mathrm{supp}(\rho_b)\).  
**DUAL-1 (R2 theory — TM `dual_channel_v0`):** joint state  
\(\{\rho_f,\rho_b,\rho_Q,\mathbf{J}_Q,\psi,\mathbf{E},\mathbf{B},\sigma,s,\gamma,\varepsilon,\mu,c\}\)  
with \(c_{\mathrm{path}}=c_{\mathrm{EM}}=1/\sqrt{\varepsilon\mu}\) and JC1–JC5.  
Forces: \(F=F^\psi+F^{\mathrm{EM}}\) (Tier A Coulomb now; Tier B Lorentz when full Maxwell).  
Hierarchy constitutive (R1 NM G3 κ scan PASS).

---

## 3. Joint constitutive constraints (dual-channel)

From TE JC1–JC5 — **required for V77-2 / G2**:

| ID | Constraint | R1 status |
|----|------------|-----------|
| **JC1** | One \(c\): path-cost locality \(\equiv 1/\sqrt{\varepsilon\mu}\) | Separate demos both use \(c=1\); **joint not stress-tested** |
| **JC2** | Budget identity; free-gauge energy in free ledger | Stated; not joint-accounted |
| **JC3** | \(\rho_b\not\equiv\rho_Q\) | NM dual0 PASS (independence by construction) |
| **JC4** | No second-sector \(G\) / \(\alpha\) foreign constants | Occam tags; residual R1 Poisson iso |
| **JC5** | Free-origin tags both channels | Separate sector tags; joint sandbox missing |

**V77-2 bar:** one medium, \(\rho_b\) and \(\rho_Q\) together, shared \(c\), theory + **joint numeric**.

---

## 4. What is proven vs residual (honest)

### Proven workable (theory ↔ numeric congruence at stated tier)

| Claim | Evidence |
|-------|----------|
| F1-3D path-cost monism | v76 `goal2_PC3D_workable`; D-ψ-3D LIVE |
| Free-gauge Coulomb + wave lite | NE KG1–4 PASS; D-EM-* LIVE_PASS |
| Free medium wave at \(c\) | ND D-DYN-free-wave-c PASS |
| Naïve inertial \(m=\int\rho_b\) dead | ND D-DYN-j5-formfactor PARTIAL/G1 PASS |
| Multi-lock ψ + dual-source lite | NM G0–G3 PASS |
| Shared monist language across 3 foci | TE/TD/TM packages + WORLD |

### Residuals (do not hide)

| ID | Residual | Blocks |
|----|----------|--------|
| **R-iso** | Poisson-form isomorphism (ψ and Φ) | Monism needs Occam/tags, not fit alone |
| **R-J5** | Form factor; non-universal \(s^*\) | Full V77-3 / operational inertia |
| **R-joint** | No single joint ψ+Maxwell sandbox | **V77-2** |
| **R-Max** | Full time-dep Maxwell not numeric | **Human terminal done** |
| **R-C1** | Rods/clocks operational incomplete | Dynamics soft |
| **R-force** | \(\alpha_\psi\) not dynamics-derived (virtual_work_v0) | Matter force depth |
| **R-SCP** | Q-ball kernel does not prove monism | Vocabulary only |

---

## 5. Terminal success / failure (campaign)

| Outcome | Criteria |
|---------|----------|
| **(A) Cohesive agreement** | Unified WORLD accepted; dual-channel joint LIVE; full Maxwell theory+numeric; G1–G2 declared; no fatal residual contradiction; demos stacked under one model |
| **(B) NULL / DISPROVE** | KILL_CRITERIA K1–K4 met (irreducible dualism, shared-\(c\) contradiction, ledger war, seed regression) — honest stop |

**Not enough for (A):** Maxwell-lite alone; separate single-channel PASSes; language unity without joint medium.

---

## 6. Focus interfaces (no private worlds)

```text
                    one continuum · free/bound · c · budget
         ┌──────────────────┼──────────────────┐
         ▼                  ▼                  ▼
      EM (C-A)         Dynamics            Matter
   full Maxwell      J5-β, timedep ψ      locks, dual source
   NE time-dep         ND triad             NM multi-lock
         └──────────────────┼──────────────────┘
                            ▼
              TU: WORLD_v1 → V77-4 or V77-K
```

| Focus | Must share | Specialize | Terminal block if… |
|-------|------------|------------|-------------------|
| EM | \(c\), monist tags, Cont | \(\varepsilon,\mu\), full M1–M4 | Full Maxwell unreachable monistically |
| Dynamics | F1 + free \(c\) | J5-β, T1/T2, C1 | Inertia forces foreign mass with no renorm story |
| Matter | Lock = bound continuum | Dual source, hierarchy | Locks require second substance |

---

## 7. Document control

| Version | Change |
|---------|--------|
| v0 | R1 map; free-channel table sketch |
| **v1** | Absorb TE/TD/TM; R1 LIVE numerics; full Maxwell terminal track; JC table; (A)/(B) stop |

Related: `DEMOS.md`, `KILL_CRITERIA.md`, `UNIFICATION_SCORE_r2.md`.
