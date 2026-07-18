# Lock Ontology v0 — Matter Focus (Theory)

**Agent:** TM (Theory — Matter)  
**Date:** 2026-07-18  
**Round:** 1  
**Status:** draft package for V77-1 Matter track  
**Seed:** v76 F1-3D free-capacity monism (`v76/CONVERGENCE.md`, `v76/work/A/THEORETICAL_PACKAGE_v1.md`)  
**Does not claim:** real \(G\), \(\alpha_{\mathrm{EM}}\), particle masses; that fixed-grid SCP Q-balls prove monism; full multi-body dynamics closed.

---

## 0. Mandate and scope

FOCI Matter open problems this note addresses:

1. What a **lock** is (stability, formation criteria) in monist language.
2. Multi-lock forces: **path-cost** vs **Coulomb**.
3. **Dual source**: \(\rho_b\) (or related ledgers) → \(\psi\) and gauge.
4. Hierarchy: EM strong at small scale, path-cost weak at large — **constitutive**, not hand-waved.
5. Vocabulary bridge to SCP Q-balls / charge **without** monism-from-kernel claims.

**Interface:** TE owns free-gauge Maxwell-lite; TM owns lock as mass-form object + how it couples to both free channels. TD owns inertia/J5 of a lock in free medium. Shared primitives only (one continuum, free/bound budget, local \(c\)).

---

## 1. What a lock *is*

### 1.1 Definition (monist)

A **lock** is a compact region of the single continuum in which a positive amount of continuum budget is held in a **bound / mass-form** state that is:

| Property | Meaning |
|----------|---------|
| **Localized** | Support of bound density \(\rho_b\) is compact (or exponentially localized) on free-medium scales |
| **Ledger-positive** | Unlock energy \(E_\star = \int \mathcal{E}[\rho_b]\,dV > 0\) |
| **Mass-form** | Rest mass \(M = E_\star / c^2\) with \(c\) = free locality (v76 ledger) |
| **Free-depleting** | Bound occupancy reduces free capacity available for free response (budget identity, integral or strong) |
| **Stable (or metastable)** | Against free-medium fluctuations, radiation channels, and weak external free-response gradients, for a timescale \(\tau_{\mathrm{lock}} \gg\) free light-crossing of the core |

**One-line:** a lock is **field locked into particle-form** — not a second substance placed on a stage.

**Not a lock:**

- Pure free radiation (free packets in flight; \(E_\star=0\) as rest unlock ledger).
- A painted “matter density” that does not deplete free budget (dualist mass primitive — NC-2.5).
- A hand Poisson source \(\rho\) with idle free continuum (dead monist claim; v76).

### 1.2 Ontology relation to free channels

The continuum admits (at least) two **free sibling channels** that locks couple to (TE + TM interface):

| Channel | Free DOF | Role near locks |
|---------|----------|-----------------|
| **Path-cost / free capacity** | \(\psi\), conductivity \(\sigma(\rho_f)\) | Long-range free-signal optical cost; exterior \(\psi\sim 1/r\) (F1-3D) |
| **Free gauge** | \(A_\mu\) or \((\mathbf{E},\mathbf{B})\) + constitutive \(\varepsilon,\mu\) | Coulomb + EM waves at same local \(c\) |

Locks are **not** these free DOFs. Locks are bound states that **source** free DOFs through ledgers (dual source, §3) and **feel** free-channel gradients as forces (§2).

### 1.3 Minimal state for one lock \(L_i\)

In Round-1 theory language (effective, continuum-level):

\[
L_i \;=\;
\bigl(
  \mathbf{X}_i,\;
  E_{\star,i},\;
  q_i,\;
  \mathcal{C}_i
\bigr)
\]

| Symbol | Meaning |
|--------|---------|
| \(\mathbf{X}_i\) | Effective center (centroid of \(\rho_b^{(i)}\)) |
| \(E_{\star,i}\) | Bound unlock ledger (rest) |
| \(q_i\) | Gauge charge ledger (signed; see §3.2) — may be 0 |
| \(\mathcal{C}_i\) | Core / constitution: shape profile class, stability class, optional multipoles |

Continuum densities recovered as superpositions of cores:

\[
\rho_b(\mathbf{x}) = \sum_i \rho_b^{(i)}(\mathbf{x}-\mathbf{X}_i),\qquad
\rho_q(\mathbf{x}) = \sum_i q_i\, f_i(\mathbf{x}-\mathbf{X}_i)
\]

with \(\int \rho_b^{(i)}=E_{\star,i}\) (v1 units), \(\int f_i=1\).

**Budget (seed):** \(\rho_f + \rho_b = \rho_{\mathrm{tot}}\) (strong form in sandboxes; integral form admissible).

---

## 2. Stability criteria (sketch)

Full dynamical stability is TD/ND + full continuum evolution territory. Matter focus needs **killable operational criteria** so NM can score “this is a lock” vs “this is a smear.”

### 2.1 Classes

| Class | Name | Sketch |
|-------|------|--------|
| **S0** | **Ledger lock** | Compact \(\rho_b>0\), \(E_\star>0\), free deficit positive at core; no dynamics required. Minimal monist *mass-form marker*. |
| **S1** | **Quasistatic lock** | S0 + stationary under free-response solve (F1-3D, gauge statics); no spontaneous delocalization of core on free grid. |
| **S2** | **Metastable lock** | S1 + survives free waves / weak external \(\psi\) or \(A\) gradients for \(T \ge N_{\mathrm{cross}}\, L_{\mathrm{core}}/c\). |
| **S3** | **Bound multi-lock** | ≥2 locks with net force channel closed into orbit or bound pair for \(T \ge N_{\mathrm{orb}}\) orbital times (optional Round-1+). |

**Round-1 default target:** S0–S1 for dual-source demos; S2 if dynamics sandboxes exist; S3 soft.

### 2.2 Operational checks (theory → NM gates)

For an isolated candidate region \(\Omega\):

1. **Compact support / localization:** \(\int_{\Omega}\rho_b = E_\star\), exterior \(\rho_b\) below threshold; second moment \(R_{\mathrm{rms}}^2 = E_\star^{-1}\int |\mathbf{x}-\mathbf{X}|^2\rho_b\) finite and \(\ll L_{\mathrm{box}}\).
2. **Free deficit:** free budget reduced at core relative to vacuum floor (T1 monist link) — \(\delta\rho_f|_{\mathrm{core}} < 0\) if strong budget, or integral free capacity removed \(\sim E_\star\).
3. **Unlock ledger:** \(E_\star\) is the work to return region to free floor (single ledger with mass; J5 residual owned by Dynamics).
4. **No free-only mass:** if \(\rho_b\equiv 0\) but “mass” is assigned by hand, **not a lock**.
5. **Radiation floor (S2):** energy leaving \(\Omega\) as free packets over test window \(\ll E_\star\) (or locked fraction stable).

### 2.3 Formation sketch (not required numeric Round 1)

Formation = free → bound conversion:

\[
\partial_t E_\star \ge 0
\quad\text{when free capacity concentrates past a lock threshold},
\]

with free budget depleted by \(\Delta E_{\mathrm{free}} \simeq -E_\star\) (plus radiation).  
**Formation mechanism is not unique** at v0: any monist continuum that can create compact \(\rho_b\) with free deficit is eligible. SCP-style soliton formation is a **vocabulary analogy**, not a monism proof (§6).

### 2.4 Kill (lock ontology)

| ID | Kill if… |
|----|----------|
| LK1 | “Lock” requires a second substance that never participates in free/bound budget |
| LK2 | Stability defined only by external stage geometry, not free-medium response |
| LK3 | Free deficit never appears when \(E_\star>0\) (ledger link broken — v76 KP3) |
| LK4 | Multi-lock forces require foreign potentials unused by free \(\psi\)/gauge dynamics |

---

## 3. Dual-source design: ledgers → \(\psi\) and gauge

### 3.1 Problem

v76 single-source free law:

\[
-\nabla\cdot(\sigma\nabla\psi) = s\,\rho_b.
\tag{F1-3D}
\]

EM free channel needs a **charge** source for Gauss law (TE), e.g.

\[
\nabla\cdot\mathbf{D} = \rho_q,\qquad \mathbf{D}=\varepsilon\mathbf{E}.
\tag{G}
\]

**Dual-source question:** how do lock ledgers feed **both** free channels without dualist “matter on a stage”?

### 3.2 Two ledger options (both monist-eligible)

#### Option A — Single bound density, two couplings (simple dual)

Same \(\rho_b\) appears in both free laws with constitutive weights:

\[
-\nabla\cdot(\sigma\nabla\psi) = s\,\rho_b,
\qquad
\nabla\cdot(\varepsilon\mathbf{E}) = g_q\,\chi\,\rho_b.
\tag{DS-A}
\]

Here \(\chi(\mathbf{x})\in\{-1,0,+1\}\) (or continuous) is a **lock charge-sign field** living on the bound sector — still continuum state, not foreign charge fluid. Total charge of lock \(i\):

\[
q_i = g_q\int \chi\,\rho_b^{(i)}\,dV.
\]

**Pros:** one mass-form density; easy multi-lock numerics.  
**Cons:** charge and mass forced proportional unless \(\chi\) structure is rich; same-\(\rho_b\) isomorphism residual for both channels.

#### Option B — Split bound ledgers (recommended for hierarchy)

Partition bound budget into **neutral mass-form** and **charged mass-form** (still one continuum):

\[
\rho_b = \rho_m + \rho_{b,q},\qquad
\rho_q = g_q\,\rho_{b,q}^{\mathrm{(signed)}},
\tag{DS-B}
\]

with

\[
\begin{aligned}
-\nabla\cdot(\sigma\nabla\psi) &= s\,\rho_b = s(\rho_m + |\rho_{b,q}|_{\mathrm{mass}}),\\
\nabla\cdot(\varepsilon\mathbf{E}) &= \rho_q.
\end{aligned}
\tag{DS-B'}
\]

Interpretation:

- \(\rho_m\): bound occupancy that sources path-cost only (neutral lock mass).
- \(\rho_{b,q}\): bound occupancy that carries gauge charge **and** still counts in unlock/mass ledger (charged lock has both \(q\) and \(E_\star\)).

**Pros:** \(q/M\) free per lock species; hierarchy (§4) natural; TE charge ledger clean.  
**Cons:** one more continuum field; NM must track two bound densities (or \(\rho_b\) + sign).

**v0 recommendation:** implement **DS-B lite** — store \(\rho_b\ge 0\) and \(\rho_q\) (signed) with constraint that charged support \(\subset\) bound support:

\[
\mathrm{supp}(|\rho_q|) \subseteq \mathrm{supp}(\rho_b).
\tag{Supp}
\]

Mass source for \(\psi\): \(\rho_b\). Gauge source: \(\rho_q\). Neutral locks: \(\rho_q=0\), \(\rho_b>0\). Pure free charge without mass forbidden by (Supp) unless TE defines a free charge mode (out of Matter Round-1).

### 3.3 Shared free-medium identity

Both free solves live on the **same** continuum geometry and **same** local \(c\):

\[
c = \frac{1}{\sqrt{\varepsilon\mu}}
\quad\text{(EM free locality; TE)}
\quad=\quad
c_{\mathrm{path}}
\quad\text{(path-cost locality; v76)}.
\tag{c-share}
\]

**Not required:** \(\psi \equiv |\mathbf{A}|\) or identifying path-cost with Coulomb potential. Sibling channels; shared medium; shared \(c\); shared budget.

### 3.4 Static dual free system (Round-1 package)

On a fixed multi-lock \(\rho_b,\rho_q\):

\[
\boxed{
\begin{aligned}
-\nabla\cdot\bigl(\sigma(\rho_f)\nabla\psi\bigr) &= s\,\rho_b,\\
\nabla\cdot\bigl(\varepsilon\mathbf{E}\bigr) &= \rho_q,\\
\rho_f + \rho_b &= \rho_{\mathrm{tot}},\\
\mathrm{supp}(|\rho_q|) &\subseteq \mathrm{supp}(\rho_b).
\end{aligned}
}
\tag{DUAL-0}
\]

Vacuum linearization: \(\sigma=\sigma_0\), \(\varepsilon=\varepsilon_0\):

\[
-\sigma_0\nabla^2\psi = s\rho_b,\qquad
\varepsilon_0\nabla\cdot\mathbf{E} = \rho_q.
\]

Exterior multipoles (isolated locks, 3D):

\[
\psi \sim \frac{s E_\star}{4\pi\sigma_0 r},\qquad
\Phi_{\mathrm{EM}} \sim \frac{Q}{4\pi\varepsilon_0 r}
\quad(\mathbf{E}=-\nabla\Phi_{\mathrm{EM}}).
\]

**Demo ID proposal:** `D-MAT-dual0` (static dual free response of ≥2 locks).

### 3.5 Kill (dual source)

| ID | Kill if… |
|----|----------|
| DS1 | Gauge source requires foreign charge not tied to bound support (and TE forbids free-only charge) |
| DS2 | \(\psi\) and EM must live on different continua to get both \(1/r\) multipoles |
| DS3 | Enforcing (Supp) makes either multipole fail systematically |
| DS4 | Shared-\(c\) identity forces \(\psi=\Phi_{\mathrm{EM}}\) identically (sibling collapse) under all constitutive choices |

---

## 4. Multi-lock force taxonomy

### 4.1 Two force channels

For locks \(i,j\) with separation \(\mathbf{r}_{ij}=\mathbf{X}_i-\mathbf{X}_j\), \(r=|\mathbf{r}_{ij}|\gg R_{\mathrm{core}}\):

| Channel | Origin | Force sketch (weak, static) | Sign structure |
|---------|--------|------------------------------|----------------|
| **Path-cost / free-capacity** | Gradients of \(\psi\) (and \(\sigma\)) sourced by \(\rho_b\) | \(\mathbf{F}_{ij}^{\psi} \propto - G_{\mathrm{eff}} M_i M_j\,\hat{\mathbf{r}}_{ij}/r^2\) | **Always attractive** for \(M>0\), \(G_{\mathrm{eff}}>0\) (v76 \(G_{\mathrm{eff}}=\gamma s c^4/(8\pi\sigma_0)\)) |
| **Coulomb / free gauge** | \(\mathbf{E}\) sourced by \(\rho_q\) | \(\mathbf{F}_{ij}^{C} = k_C q_i q_j\,\hat{\mathbf{r}}_{ij}/r^2\) | **+/− by charge product** |

Total quasistatic force on lock \(i\) from free channels:

\[
\mathbf{F}_i = \mathbf{F}_i^{\psi} + \mathbf{F}_i^{C} + \mathbf{F}_i^{\mathrm{higher}}
\tag{F-tot}
\]

where \(\mathbf{F}^{\mathrm{higher}}\) includes induction, radiation reaction, core overlap, nonlinear \(\sigma(\rho_f)\), magnetostatics — **deferred**.

### 4.2 Path-cost force (constitutive derivation sketch)

From v76: free signals see path cost \(\ell=\ell_0+\gamma\psi\). Slow locks as quasistatic free-capacity “test masses” feel an effective potential proportional to \(\psi\) if their bound ledger couples back to free capacity (action principle / virtual work on \(E_\star\) in a background \(\psi\)):

\[
V_i^{\psi} \approx \alpha_\psi E_{\star,i}\,\psi(\mathbf{X}_i)
\quad\Rightarrow\quad
\mathbf{F}_i^{\psi} = -\nabla_{\mathbf{X}_i} V_i^{\psi}.
\]

With exterior \(\psi_j \approx s E_{\star,j}/(4\pi\sigma_0 r)\) and \(\alpha_\psi\) fixed by free constitutive law so that

\[
G_{\mathrm{eff}} = \frac{\gamma s c^4}{8\pi\sigma_0}
\]

matches two-body Newtonian form with \(M=E_\star/c^2\).

**Honest residual:** full two-body derivation of \(\alpha_\psi\) from free dynamics is **not** closed at Round 1 (related to J5 / free drag — TD). For Matter demos, use **virtual-work / prescribed force law consistent with exterior multipoles**, tagged `force_closure=virtual_work_v0`, not “proven geodesic.”

### 4.3 Coulomb force

Standard free-gauge virtual work:

\[
V_i^{C} = q_i\,\Phi_{\mathrm{EM}}(\mathbf{X}_i),\qquad
\mathbf{F}_i^{C} = q_i\mathbf{E}(\mathbf{X}_i)
\]

with \(\Phi_{\mathrm{EM}}\) from (G) excluding self-field (or regularized).

### 4.4 Taxonomy table (multi-lock phenomenology)

| Config | \(q_i\) | Dominant far force | Near / core |
|--------|--------|--------------------|-------------|
| Two neutral locks | 0,0 | Attractive path-cost \(\propto M_i M_j/r^2\) | Core merge / free-deficit overlap |
| Opposite charges, light \(M\) | \(+q,-q\) | Attractive Coulomb \(\propto |q_i q_j|/r^2\) | Possible bound pair if dynamics exist |
| Like charges, light \(M\) | \(+q,+q\) | Repulsive Coulomb | Path-cost still attracts; usually Coulomb wins at small \(r\) if hierarchy §4 holds |
| Neutral + charged | \(0,q\) | Path-cost only (to leading multipole) | Weak polarization residuals possible |
| Heavy neutral + light charged pair | mixed | Path-cost sets “gravity”; Coulomb sets internal EM | Atom analog — Stage language only |

### 4.5 Path-cost vs Coulomb — not the same potential

| | Path-cost \(\psi\) | Coulomb \(\Phi_{\mathrm{EM}}\) |
|--|-------------------|-------------------------------|
| Source | \(\rho_b\) (mass ledger) | \(\rho_q\) (charge ledger) |
| Exterior | \(\sim M/r\) | \(\sim Q/r\) |
| Force sign | \(M_i M_j > 0\) → attract | \(q_i q_j\) → attract/repel |
| Free wave sibling | free capacity dynamics / eikonal | EM waves at \(c\) |
| Dualist twin residual | Poisson isomorphism (v76 R1) | Standard Maxwell vacuum residual |

**Kill confusion:** treating path-cost attraction as “just Coulomb with \(q=M\)” **collapses siblings** and fails neutral multi-lock demos.

---

## 5. Hierarchy: EM strong (small scale) vs path-cost weak (large) — constitutive

### 5.1 Requirement

Target structural pattern (analogy, not fitted constants):

- At **atomic / lock-core scales**, free-gauge (Coulomb) forces dominate inter-lock binding when charges are nonzero.
- At **macro / multi-lock assembly scales**, path-cost (mass-mass) forces dominate for neutral composites.

This must follow from **constitutive parameters** of the free medium + ledgers, not “because we want atoms then planets.”

### 5.2 Coupling constants from free medium

Define dimensionless and dimensionful free constants:

| Constant | Role |
|----------|------|
| \(\sigma_0,s,\gamma\) | Free-capacity channel → \(G_{\mathrm{eff}}\) |
| \(\varepsilon_0\) (and \(\mu_0\) with \(c\)) | Free-gauge channel → \(k_C = 1/(4\pi\varepsilon_0)\) |
| \(c\) | Shared locality |

Two-body force magnitudes at separation \(r\):

\[
F^{\psi} \sim \frac{G_{\mathrm{eff}} M_i M_j}{r^2},\qquad
F^{C} \sim \frac{k_C |q_i q_j|}{r^2}.
\]

Ratio (same \(r\)):

\[
\boxed{
\frac{F^{C}}{F^{\psi}}
=
\frac{k_C |q_i q_j|}{G_{\mathrm{eff}} M_i M_j}
=
\frac{c^4\,|q_i q_j|}{G_{\mathrm{eff}}\, E_{\star,i} E_{\star,j}}
\cdot k_C
}
\tag{H-ratio}
\]

with \(M=E_\star/c^2\). Using \(G_{\mathrm{eff}}=\gamma s c^4/(8\pi\sigma_0)\):

\[
\frac{F^{C}}{F^{\psi}}
=
\frac{8\pi\sigma_0\, k_C\, |q_i q_j|}{\gamma s\, E_{\star,i} E_{\star,j}}.
\tag{H-ratio'}
\]

**All quantities on the RHS are free-medium or lock-ledger constitutive data.** No foreign “gravity sector vs EM sector” ontology.

### 5.3 Scale hierarchy without new ontology

Let elementary charged locks have characteristic

\[
\lambda_q := \frac{|q|}{E_\star}
\quad\text{(charge per unlock energy)},
\qquad
\Xi := \frac{8\pi\sigma_0\, k_C}{\gamma s}
\quad\text{(free-medium hierarchy factor)}.
\]

Then

\[
\frac{F^{C}}{F^{\psi}} = \Xi\,\lambda_{q,i}\,\lambda_{q,j}.
\tag{H-λ}
\]

**Constitutive hierarchy design (Round-1):**

1. Choose free constants so \(\Xi \lambda_0^2 \gg 1\) for elementary charged species with \(\lambda_q\sim\lambda_0\) (EM-dominated pair force).
2. Neutral composites: internal \(q_{\mathrm{tot}}=0\), residual multipoles, but net \(M_{\mathrm{comp}}=\sum E_\star/c^2\) large → external force is path-cost only:

\[
F^{\psi}_{\mathrm{comp}} \sim \frac{G_{\mathrm{eff}} M_A M_B}{r^2},\qquad F^{C}_{\mathrm{comp}}\approx 0.
\]

3. **Crossover radius** is not fundamental: for two charged locks the ratio (H-λ) is **\(r\)-independent** at leading \(1/r^2\). Hierarchy of *phenomena* comes from **which objects exist** (neutral bulk vs charged micro), not from EM “running” with \(r\) at this multipole order.

**Optional higher-order hierarchy (not required Round 1):** screening, multipole falloff, or \(\sigma(\rho_f)\) nonlinearities can introduce scale dependence; document if used.

### 5.4 What is *not* claimed

- Numerical match to real \(\alpha_{\mathrm{EM}}\) vs \(G\).
- That path-cost **is** GR beyond weak multipole tier.
- That choosing \(\Xi\lambda_0^2\gg 1\) is derived from deeper theory — it is a **constitutive existence claim**: there exist free-medium parameters realizing the pattern.

### 5.5 Kill (hierarchy)

| ID | Kill if… |
|----|----------|
| H1 | Only way to get EM≫path-cost is to disable \(\psi\) channel (not constitutive coexistence) |
| H2 | Neutral composites inevitably retain Coulomb-class net charge by dual-source construction |
| H3 | Shared \(c\) forces \(G_{\mathrm{eff}}\) and \(k_C\) into a unique ratio forbidding \(\Xi\lambda_0^2\gg 1\) |

---

## 6. Vocabulary bridge to SCP particles (without monism proof)

### 6.1 Purpose

SCP U(1)/Q-ball language is a **mature particle lab** (balls, charge \(Q\), flavor, nuclei). v77 monism is a **continuum ontology + free response**. Bridging **words** helps agents and demos; it must **not** smuggle “Q-balls on fixed grid prove free-capacity monism” (v76 dead branch).

### 6.2 Bridge table (language only)

| SCP (U(1) era) term | Monist lock language (v77 Matter) | Shared? | Not shared / residual |
|---------------------|-----------------------------------|---------|------------------------|
| Ball / Q-ball | Lock (S0–S2 candidate) | Localized bound energy | Fixed-grid kernel; not free-capacity monism proof |
| \(Q\) (global U(1) charge) | Gauge charge ledger \(q\) **or** separate Noether charge | Conserved scalar | Identity of \(Q\) with monist \(\rho_q\) is **analogy** until TE maps Noether↔Gauss |
| \(\omega\), \(Q_a\), flavor | Core constitution \(\mathcal{C}\) / internal clocks | Internal structure label | Not required for dual-source statics |
| \(\eta\)-soliton / ledger-closed | Unlock ledger stability / process-form | Stability bookkeeping spirit | Different equations |
| Nucleus / multi-ball droplet | Multi-lock bound (S3) / fused \(\rho_b\) | Multi-center bound form | Fusion dynamics not in DUAL-0 |
| \(A\) gauge (kernel) | Free gauge channel | Coulomb/Gauss sibling | Kernel \(A\) is not automatically monist free medium |
| \(m=E/c^2\) particle mass | \(M=E_\star/c^2\) | Ledger slogan | SCP mass from field energy on fixed stage residual |
| Carbon Stage 2A liquid drop | Single heavy lock or fused \(\rho_b\) blob | Compact multi-nucleon analog | Not path-cost gravity demo |
| Absorbing BC / box | Free-medium truncation | Numerics | Dualist stage risk if free DOF idle |

### 6.3 Explicit non-claims

1. **Existence of gauged Q-balls does not prove F1-3D monism.**  
   They prove a field theory can host stable charged lumps on a prescribed lattice/stage.
2. **Conserved \(Q\) in scp_sim does not by itself identify monist dual-source \(\rho_q\).**  
   Mapping requires TE constitutive + Gauss law monist package.
3. **Fusion / carbon parks are nuclear-style bound ledgers**, useful metaphors for multi-lock \(\rho_b\), not weak-field path-cost congruence.
4. **No scp_sim edits** for this bridge; NM sandboxes are independent Python media under `work/NM/`.

### 6.4 Allowed positive use of SCP vocabulary

- Call a compact positive \(\rho_b\) region a “ball-like lock” in prose.
- Use “charge” for \(q\) when TE Gauss law is in play.
- Use “binding” for S3 multi-lock when forces close.
- Cite SCP results as **existence intuitions** for stable locks in *some* continuum theory.

---

## 7. Demo IDs and interfaces

| Demo ID | Theory element | Numeric owner | Pass sketch |
|---------|----------------|---------------|-------------|
| `D-MAT-lock-S0` | §1–2 S0 lock | NM | Compact \(\rho_b\), \(E_\star>0\), free deficit |
| `D-MAT-dual0` | §3 DUAL-0 | NM (+ TE/NE for gauge) | \(\psi\sim 1/r\) from \(M\); \(\Phi\sim 1/r\) from \(Q\); (Supp) |
| `D-MAT-force-tax` | §4 taxonomy | NM | Neutral attract (path-cost); like-charge repel (Coulomb); opposite attract |
| `D-MAT-hier` | §5 \(\Xi\lambda^2\) | NM | Parameter set with \(F^C/F^\psi \gg 1\) for charged pair; neutral pair Coulomb-null |

**Interfaces:**

| To | Message |
|----|---------|
| **TE** | Charge = lock–gauge ledger \(\rho_q\) with (Supp); shared \(c\); sibling not identical to \(\psi\) |
| **NE** | Need static Coulomb solve coexisting with F1 \(\psi\) on same grid for `D-MAT-dual0` |
| **TD** | Force on lock from \(\nabla\psi\) needs dynamics closure; J5 residual acknowledged |
| **ND** | Optional: time-dep free response with moving locks later |
| **NM** | Gates in `FOR_NM_gates_v0.md` |
| **TU** | Register demos; residual R1 still applies per channel |

---

## 8. Residuals (honest)

| ID | Residual |
|----|----------|
| MR1 | Poisson isomorphism for \(\psi\) (v76 R1) and Maxwell vacuum for \(\Phi\) — monism is ontology + budget, not PDE uniqueness |
| MR2 | Two-body path-cost force coefficient \(\alpha_\psi\) not derived from free dynamics (virtual-work_v0) |
| MR3 | Full S2/S3 stability needs time-dependent free medium |
| MR4 | Map SCP Noether \(Q\) ↔ monist \(\rho_q\) incomplete without TE |
| MR5 | Hierarchy is constitutive existence, not derived unique constants |

---

## 9. Round-1 deliverable checklist

| Item | Status |
|------|--------|
| Lock definition | §1 |
| Stability sketch S0–S3 | §2 |
| Dual-source DUAL-0 + DS-B lite | §3 |
| Force taxonomy path-cost vs Coulomb | §4 |
| Constitutive hierarchy | §5 |
| SCP vocabulary bridge + non-claims | §6 |
| FOR_NM gates | `FOR_NM_gates_v0.md` |

**Package claim:** Matter focus has a **workable lock + dual-force theory draft** aligned with v76 seed and V77 dual-channel milestone — killable, not yet numerically closed.
