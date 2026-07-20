# v78 PHYSICS_RELATIONS — R1–R10 Relation Matrix

**Date:** 2026-07-19  
**Owner:** W (World) · co-agree CP-W with U  
**Status:** **FREEZE v0** — Phase 0 complete  
**Inputs:** `GOALS.md` §1–2; v76 `CONVERGENCE.md` + F1-3D package; v77 `PHASE2_CONVERGENCE.md` + WORLD/UNIFIED packages; CONCEPT.md; v75 `STATE.md` / multi-fabric; v74 carbon map.

**Scope rule:** monist sandboxes (v76/v77) supply **language and free-channel laws** for R1, R2, R4, R9 (and Maxwell sibling for R3/R6). SCP kernel particles supply **engineered stable lumps** and measured nuclear/atomic ladders. Neither replaces the other; this document maps **how each relation is expressed in both stacks**.

---

## 0. Two stacks, one campaign

| Stack | What it is | What it is not |
|-------|------------|----------------|
| **Monist free-capacity + free gauge** (v76 F1-3D; v77 M1∧RC1) | One continuum; free/bound budget; \(\psi\) path-cost; Maxwell free channel; dual-source locks | Not a production particle engine; not Q-ball fusion numerics |
| **SCP kernel** (complex Cosserat + gauged U(1); multi-fabric C/Q/L) | Stable Q-balls, nuclei, P/N analogs, atom path on lattice | Does not prove monism; \(\psi\) free-capacity not co-evolved in `scp_sim` |

**Shared ontology claim for v78:** particles are **stable localized field configurations** (locks / Q-balls / multi-fabric composites). Mass, charge, binding, and stability are **ledger relations on field content**, not point masses on empty stage.

```text
  monist free channels          SCP particle ladder
  ───────────────────          ───────────────────
  c locality (free)     ←→     lattice wave/gauge speed
  m = E_★/c²            ←→     E = ∫T₀₀ soliton mass
  ∇·E ∝ ρ_em            ←→     discrete Gauss + Q_em
  ψ path-cost (opt.)    ←→     (not required nuclear T)
  locks / free budget   ←→     bound lumps / continuum rest
  dual ρ_b ⟂ ρ_Q        ←→     bag vs EM (C vs Q fabric)
```

---

## 1. Relation matrix (summary)

| ID | Relation | Monist (v76/v77) | Kernel particles (SCP) | Particle consequence |
|----|----------|------------------|------------------------|----------------------|
| **R1** | Locality \(c\) | Free update bound; \(c=1/\sqrt{\varepsilon\mu}\) | Hyperbolic matter + gauge CFL / link transport | Shared causal structure |
| **R2** | Mass ↔ energy | \(m=E_\star/c^2\); J5-β inertia \(m\sim U/c^2\) | \(E=\int T_{00}\); \(E/Q\); mass defect | Binding / rest mass of lumps |
| **R3** | Charge ↔ Gauss | \(\nabla\cdot(\varepsilon E)=\rho_Q\); Cont | Discrete Gauss \(\nabla\cdot E - g\rho_{\mathrm{em}}\approx0\); `gauss_max`~1e-13 | P EM vs N⁰ ~0 EM |
| **R4** | Path-cost / gravity-like | F1-3D: \(-\nabla\cdot(\sigma\nabla\psi)=s\rho_b\); \(\ell=\ell_0+\gamma\psi\) | **Optional bookkeeping** — not in nuclear recipe | Large-scale geometry |
| **R5** | Short-range binding | Lock cohesion ontology only (not nuclear \(V_t\)) | \(V_t(s)\), phase coherence, multi-ball fusion | Nuclear force / fusion |
| **R6** | Coulomb | Free-gauge Maxwell; same/opposite sign forces | Shared \(A\); fabric \(q_f\); measured \(1/r^2\) | Atom binding; P–P repulsion |
| **R7** | Flavor / multiplet | Lock multiplet labels optional | \((\omega_a,Q_a)\); flavored profiles | p/n candidates (not EM-neutral alone) |
| **R8** | Multi-fabric isolation | Dual channel: \(\rho_b\not\equiv\rho_Q\); channels independent | Private \(s_f\); C/Q/L; no C–L bag merge | Bag vs light; atom isolation |
| **R9** | Free/bound budget | \(\rho_f+\rho_b=\rho_{\mathrm{tot}}\); free deficit at locks | Continuum vacuum vs compact lumps; charge/energy ledgers | Mass-formation language |
| **R10** | Stability | Fixed locks RC1; no free evaporation of lock labels | VK \(\mathrm{d}Q/\mathrm{d}\omega<0\); park; long-T drift | P/N/C₁₂/atom survival |

---

## 2. R1 — Locality \(c\)

### 2.1 Statement

There is a single free-signal locality bound \(c\) in free vacuum: free influence (waves, free capacity radiation, free gauge waves) cannot outrun \(c\). In a preferred monist frame, local free \(c\) is held fixed; warp is path cost around locks, not a second light-speed sector.

### 2.2 Monist map (v76/v77)

| Element | Content | Evidence |
|---------|---------|----------|
| Axiom | \(c\) = free-field locality | v76 seed hinge; WORLD_v1 JC1 |
| EM constitutive | \(c=1/\sqrt{\varepsilon\mu}\) | v77 TE full Maxwell; NE M1 gates G2–G3 |
| Free capacity | Hyperbolic free medium \(\partial_t^2\psi=c^2\nabla^2\psi\) (ND) | dual-channel shared \(c\) PASS |
| Path cost | \(\ell=\ell_0+\gamma\psi\) **at fixed local \(c\)** | F1-3D; not local GRIN |

**Kill (agreed dead):** local free-density GRIN as long-range “gravity”; \(c_{\mathrm{EM}}\neq c_{\mathrm{path}}\) without medium-frame story.

### 2.3 Kernel map (SCP)

| Element | Content |
|---------|---------|
| Matter | Second-order hyperbolic PDE for \(\Phi_a,\Theta_a\) (and multi-fabric copies) |
| Gauge | Lattice links + \(E\); massless photon analog \(A\) propagates at lattice light speed (units \(c=1\) standard) |
| Numerics | CFL / `dt` bound; absorbing BC damp outgoing waves without reflecting acausal bulk |

**Particle consequence:** all nucleon / nuclear / atom dynamics share one causal cone. Collision, radiation, and Coulomb adjustment are not superluminal bookkeeping tricks.

### 2.4 Assignment for ladder

| Stage | R1 role |
|-------|---------|
| N | Profile + evolution at standard `dt`; no superluminal seed tricks |
| P/N⁰, C | Multi-ball fusion inside light-cone of contact |
| L, A | Light sector and Coulomb orbit/park use same \(c\) as nuclear bag |

---

## 3. R2 — Mass ↔ energy

### 3.1 Statement

Rest mass of a localized object is a **ledger of continuum energy** associated with the bound configuration, divided by \(c^2\). Binding is a **mass defect** relative to free constituents.

### 3.2 Monist map (v76/v77)

| Ledger | Definition | Status |
|--------|------------|--------|
| Path-cost mass | \(M=E_\star^{\mathrm{path}}/c^2\), \(E_\star^{\mathrm{path}}=\int\rho_b\) | F1 exterior \(\psi\sim 1/r\); \(m_{\mathrm{ray}}=m_{\mathrm{ledger}}\) |
| Inertial (J5-β) | \(m_{\mathrm{inertial}}=\xi U[\psi]/c^2\) (free stress) | Naïve \(m=\int\rho_b\) **killed** as universal inertia |
| Free energy | Free-channel stresses \(U[\psi]\), \(u_{\mathrm{EM}}\) | Part of continuum ledger |

**Residual (honest):** Poisson-form isomorphism of \(\psi\); form-factor non-universality for inertia — does **not** block nuclear particle work.

### 3.3 Kernel map (SCP)

| Quantity | Definition |
|----------|------------|
| Soliton energy | \(E=\int T_{00}\) (kinetic + gradient + potential + gauge) over compact support |
| Specific energy | \(E/Q\) on the Q-ball branch |
| Mass defect | \(\Delta E = E_{\mathrm{composite}} - \sum E_{\mathrm{free}}\) after fusion/park |
| Charge pressure | \(Q^2\) term (fixed-\(Q\) / \(\omega\)) stabilizes against Derrick collapse |

Standard light nucleon template: \(\omega=1.46\), \(Q_N\approx 114\) (`f_w146_g005`); branch quantities from radial/gauged shooter.

### 3.4 Particle consequence

| Object | Mass-energy content |
|--------|---------------------|
| Nucleon | On-branch \(E(Q)\); VK stability from \(\mathrm{d}Q/\mathrm{d}\omega\) |
| Nucleus | Parked droplet \(E\), defect vs free \(N\) balls; super-critical evaporates toward \(Q_{\max}\) |
| Atom | Nuclear \(E\) + L-sector \(E\) + Coulomb field energy; large-box bookkeeping |

---

## 4. R3 — Charge ↔ Gauss

### 4.1 Statement

Electromagnetic charge is a **signed ledger** on matter that **sources** the free gauge field. Gauss law is an identity (or exact discrete constraint), not a soft diagnostic.

### 4.2 Monist map (v76/v77)

\[
\nabla\cdot(\varepsilon\mathbf{E})=\rho_Q,\qquad
\partial_t\rho_Q+\nabla\cdot\mathbf{J}_Q=0.
\]

| Item | Content |
|------|---------|
| Dual source | \(\rho_b\) → \(\psi\); \(\rho_Q\) → Maxwell; \(\mathrm{supp}(|\rho_Q|)\subseteq\mathrm{supp}(\rho_b)\) preferred |
| Independence | Neutral lock: \(\psi\neq0\), \(\|E\|=0\) (RC1 TE-IA1) |
| Numeric | M1 G5 dynamic Gauss+Cont; RC1 co-field on same grid |

**Forbidden collapse:** \(\rho_Q\equiv\rho_b\) always; \(\psi\equiv\Phi_{\mathrm{EM}}\).

### 4.3 Kernel map (SCP)

**Single-fabric (standard Q-ball):** diagonal U(1) Noether

\[
Q=\int\sum_a\big(u_a\dot v_a-v_a\dot u_a+\cdots\big)\,dV,
\]

gauged with \(D_\mu=\partial_\mu+igA_\mu\), \(g=0.05\). Discrete Gauss law flat at \(\sim10^{-13}\) (`gauss_max` floor). Tripwire: Gauss drift ⇒ implementation bug.

**Multi-fabric (v75):**

\[
\rho_{\mathrm{em}}=q_C\rho_C+q_Q\rho_Q+q_L\rho_L,\qquad
\boxed{q_C=0,\;q_Q=+1,\;q_L=-1}.
\]

| Species | EM bookkeeping |
|---------|----------------|
| Proton-analog | C bag + Q co-located → \(Q_{\mathrm{em}}\propto Q_Q\neq0\) |
| Neutron-analog | C bag only → \(Q_{\mathrm{em}}\approx0\), \(E_{\mathrm{em}}\approx0\) |
| Light / e-analog | L with \(q_L=-1\) |

Proxies: `Q_flux` \(\propto\int(\nabla\cdot E)/g\); `Q_phi` bag; per-fabric `Q_C,Q_Q,Q_L` when enabled.

### 4.4 Particle consequence

R3 is the **definitional** split P vs N⁰ and the atom neutrality condition \(Z\) L-charges cancel nuclear \(Q_{\mathrm{em}}\). Flavor multiplet alone (R7) does **not** replace R3 for EM-neutral neutron.

---

## 5. R4 — Path-cost / gravity-like

### 5.1 Statement

Around bound mass-form, free capacity \(\psi\) develops a long-range exterior (3D \(\sim1/r\)); free-signal **path cost** is \(\ell=\ell_0+\gamma\psi\) with local \(c\) held fixed. This is the monist gravity-class relation.

### 5.2 Monist map (v76/v77)

\[
-\nabla\cdot(\sigma\nabla\psi)=s\rho_b,\qquad
\psi\sim\frac{s E_\star}{4\pi\sigma_0 r}\ \text{(3D)},\qquad
G_{\mathrm{eff}}=\frac{\gamma s c^4}{8\pi\sigma_0}.
\]

| Status | Note |
|--------|------|
| v76 | `goal2_PC3D_workable` MET |
| v77 RC1 | \(\psi\) + dynamical Maxwell co-field (2D free Laplace multipole class) |
| Residual | 3D co-evolution with moving locks optional (RC2) |

### 5.3 Kernel map (SCP)

| Item | Content |
|------|---------|
| In `scp_sim` | **No** free-capacity \(\psi\) field; nuclear/atomic T does not integrate F1-3D |
| Role in v78 | **Parallel bookkeeping / ontology** — what “mass as lock” means |
| Optional later | Post-process \(\rho_b\sim T_{00}\) → \(\psi\) diagnostic; not a kill-gate for N→A |

### 5.4 Particle consequence

R4 is **optional for nuclear timescales** (GOALS §6). Nuclear binding is R5; atomic binding is R6. R4 frames why soliton energy can source large-scale geometry without a second dualist stage — for campaign completeness, not for Stage-2A carbon park.

---

## 6. R5 — Short-range binding (nuclear)

### 6.1 Statement

Compact multi-component matter coheres via a **short-range saturating potential** and **phase coherence**, producing fusion, mass defect, and liquid-drop nuclei — not Coulomb.

### 6.2 Monist map (v76/v77)

Monist packages define **locks** as compact free-depleting mass-form but **do not** implement Cosserat \(V_t(s)\) or three-component bag physics. Nuclear force is **SCP-primary**. Monist language: “lock cohesion” / multi-lock \(\psi\) attraction is **not** the nuclear residual force of CONCEPT § nuclear channel.

### 6.3 Kernel map (SCP)

\[
s=\prod_{a=0}^{2}|\Phi_a|^2,\qquad
V_t(s)=\frac{\mu}{2}\frac{s}{1+\kappa s},\quad \mu<0,\;\kappa=50.
\]

| Mechanism | Role |
|-----------|------|
| Triple-product \(s\) | No single-component particle; “color” binding |
| Phase alignment | Co-phase multi-ball → contact fusion (v71) |
| Coulomb self-repulsion | Caps branch \(Q_{\max}=921\) at \(g=0.05\) |
| η, Θ | Optional dressing; standard particle package \(\eta=0\), \(m_\theta=1.6\) |

**Carbon (v74):**

| Run | Fate |
|-----|------|
| **c6_light** (Z=6 map) | Parks mid-branch \(Q\to650\) — primary stable carbon nucleus |
| **c12_light** (A=12 free) | Super-critical \(Q\to1411>921\) — evaporates charge; not parked isotope at \(g=0.05\) |

### 6.4 Particle consequence

| Stage | R5 requirement |
|-------|----------------|
| N | Single-ball on VK branch; \(V_t\) balances charge pressure |
| C | Multi-nucleon fusion → parked droplet (Z-carbon) or characterized super-critical A=12 |
| A | Nuclear core uses R5; **must not** fuse L into bag (R8) |

---

## 7. R6 — Coulomb

### 7.1 Statement

Long-range EM force between charges is carried by the **shared free gauge field** \(A\) (kernel) / \((\mathbf{E},\mathbf{B})\) (monist Maxwell). Same-sign repel; opposite attract. Atoms bind by R6 with scale separation from R5.

### 7.2 Monist map (v76/v77)

| Law | Content |
|-----|---------|
| Gauss / Coulomb | \(E\sim 1/r^2\) exterior for pointlike \(Q\) |
| Full Maxwell (M1) | Dynamic E,B; Faraday; Ampère–Maxwell; waves at \(c\) |
| Dual-channel forces | \(F=F^\psi+F^{\mathrm{EM}}\); opposite charges → EM attract |
| RC1 | Co-field \(\psi\) ≠ EM (neutral lock test) |

### 7.3 Kernel map (SCP)

| Element | Content |
|---------|---------|
| Field | Temporal-gauge \(A_i\), conjugate \(E_i\); link-covariant derivatives |
| Coupling | \(g=0.05\); multi-fabric \(U^{q_f}=\mathrm{e}^{i q_f\theta}\) |
| Measured | Massless Coulomb analog; P–P repulsion; opposite L attraction (v75 B1–B4) |
| Atom | 6 L at \(q_L=-1\) Coulomb-bound to Z≈6 nuclear \(Q_{\mathrm{em}}\) |

**Seed rule (opposite EM):** same sign of Noether \(\omega\) on C and L when fabric \(q\) supplies EM sign flip (v75 STATE).

### 7.4 Particle consequence

| Interaction | Dominant relation |
|-------------|-------------------|
| N–N contact fusion | R5 (short) + soft R6 self-energy |
| P–P separation | R6 repel; R5 if co-phase contact |
| Nucleus–L cloud | **R6 only** for binding if R8 holds |
| Positronium analog | L⁺L⁻ via R6 (prep for atom) |

---

## 8. R7 — Flavor / multiplet

### 8.1 Statement

Internal **frequency / charge partition** \((\omega_a,Q_a)\) labels multiplet structure of a multi-component ball. Flavor is real conserved bookkeeping; it is **not** by itself EM proton vs neutron.

### 8.2 Monist map (v76/v77)

Optional lock multiplet labels \(\mathcal{C}_i\); no required Cosserat flavor dynamics. R7 is **SCP-primary**.

### 8.3 Kernel map (SCP)

| Structure | Content |
|-----------|---------|
| Symmetry | Exact U(1)³; only diagonal U(1) gauged |
| Conserved | Each \(Q_a\) separately; total \(Q=\sum Q_a\) |
| Equal-\(f\) ball | \(Q_a=Q/3\), common \(\omega\) — standard nucleon |
| Flavored | Per-component \(\omega_a\), profiles (`gen_qball_flavored`) |
| v71 multiplet | Stable flavor family; **same-sign EM** — not EM-neutral n |
| Cancel multiplet | \(\omega=(+w,+w,-w)\) leaves residual \(Q\sim Q/3\), not zero |

### 8.4 Particle consequence

| Use | Status |
|-----|--------|
| Internal quark-like bookkeeping | Required (components exist) |
| p/n-analog **EM** definition | **R3 multi-fabric** (Q on/off), not R7 alone |
| Nuclear modes / GDR-like | Optional research (Stage 2B) |

---

## 9. R8 — Multi-fabric isolation

### 9.1 Statement

Distinct matter sectors may share space and a gauge field while keeping **private short-range bags**. Light sector must not enter nuclear product \(s_C\); charge fabric can be unlocked from bag (B2) for P/N.

### 9.2 Monist map (v76/v77)

| Dual-channel idea | SCP multi-fabric analog |
|-------------------|-------------------------|
| \(\rho_b\not\equiv\rho_Q\) | Bag (C) vs EM charge (Q) |
| Independent free channels \(\psi\), Maxwell | Private bags \(s_f\) + shared \(A\) |
| Neutral lock with mass | Neutron: bag without EM source |

Monist dual-channel is the **ontological twin** of selective engagement; implementation is multi-fabric Option B (v75).

### 9.3 Kernel map (SCP)

Three fabrics C / Q / L (when multi-fabric path enabled):

```text
   C ←ε→ Q ←g→ A ←g→ L
   (bag)  (heavy EM)     (light, q_L=-1)
```

| Rule | Content |
|------|---------|
| Private \(s_f=\prod_a|\Phi_f^a|^2\) | No C–L merge at head-on (B1 proven) |
| B1 lock | \(\Phi_Q\equiv\Phi_C\) after steps — every nuclear lump charged |
| B2 unlock | `mf_lock_CQ=0` — Q independent → true P/N |
| L engagement | Shared \(A\) only; never \(s_C\) |

### 9.4 Particle consequence

| Target | Isolation need |
|--------|----------------|
| P | C+Q co-located; L absent or distant |
| N⁰ | C only at site |
| Atom | Nuclear C(+Q) + 6L; **R8 + R6**; block today = multi-fabric kernel path + L stability |

---

## 10. R9 — Free / bound budget

### 10.1 Statement

Continuum capacity splits into **free** (propagating / available) and **bound** (locked in mass-form). Forming particles **uses** free budget; free deficit near locks is monist signature.

### 10.2 Monist map (v76/v77)

\[
\rho_f+\rho_b=\rho_{\mathrm{tot}}
\quad\text{(strong or integral)}.
\]

| Evidence | Content |
|----------|---------|
| Free deficit | Core free deficit ~0.13–0.17 (v76 B) |
| Multi-lock | NM free deficit + path-cost gates |
| RC1 | Locks source \(\psi\) without equating to EM |

### 10.3 Kernel map (SCP)

| Bound | Free / continuum |
|-------|------------------|
| Q-ball cores, nuclear droplets, L lumps | Vacuum continuum; radiation; θ sector; outgoing waves to sponge |
| Conserved \(Q,Q_a\) locked in lumps | Radiated charge/energy absorbed at BC |
| Gauge energy in/near lumps | Propagating \(A,E\) radiation |

SCP does not store \(\rho_f,\rho_b\) arrays; the **language** maps: particle = bound ledger; vacuum + radiation = free sector.

### 10.4 Particle consequence

Mass formation (R2) is locking continuum energy into R10-stable lumps. Super-critical nuclei **return** charge to free/radiation channels until on-branch (A=12 story). Atoms park when bound nuclear + bound L + Coulomb field settle without fatal free radiation.

---

## 11. R10 — Stability

### 11.1 Statement

A particle target **survives**: charge/energy retention, on-branch diagnostics, park criteria, no fatal radiation, Gauss floor held.

### 11.2 Monist map (v76/v77)

| Level | Content |
|-------|---------|
| RC1 | Fixed multi-locks — composition stability of dual sources |
| Not claimed | Moving-lock RC2; dynamical Q-ball formation in monist sandbox |

Monist stability ≠ VK Q-ball theorem; it only shows dual-channel media do not self-destruct under static locks.

### 11.3 Kernel map (SCP)

| Gate | Content |
|------|---------|
| **VK** | \(\mathrm{d}Q/\mathrm{d}\omega<0\) classical stability branch |
| **Long-T** | Charge retention (e.g. 99.99985% over 1000 t.u. bare ball) |
| **Park** | \(c_{Q,\mathrm{park}}\) / rolling death checks; droplet not evaporating |
| **Gauss** | `gauss_max` ~ 1e-13 floor |
| **Fragmentation** | Connected-component / multi-cluster tracking |
| **L survival** | massL / Ql drift bounds (atom) |

**Defaults:** \(\eta=0\), \(m_\theta=1.6\), \(g=0.05\), absorbing BC; for \(\eta>0\) must use `eta_qflow` stationary seed.

### 11.4 Particle consequence (success bars)

| Target | Stability bar |
|--------|---------------|
| **N** | On VK branch; long-T \(Q\) retention |
| **P / N⁰** | Parked / long-T; EM ledger matches species (F17–F19 class) |
| **C nucleus** | Z-carbon park (c6_light class); A=12 characterized if super-critical |
| **Atom** | Nuclear park + L cloud without merge; long-T; Gauss floor |

---

## 12. Who owns which relation (campaign)

| Relation | Primary evidence owner | Phase consumers |
|----------|------------------------|-----------------|
| R1 Locality | W (monist+kernel map); N runtime | all |
| R2 Mass-energy | W + N (E/Q); C (defect) | N,P,C,A |
| R3 Charge-Gauss | P (P/N); W map; U board | P,C,A |
| R4 Path-cost | W (monist); optional post | U completeness |
| R5 Nuclear binding | N, C | N,C,A core |
| R6 Coulomb | L, A, P (P–P); v75 | P,L,A |
| R7 Flavor | N (bookkeeping); P note | N,P |
| R8 Multi-fabric | P, L, A; v75 | P,L,A |
| R9 Budget | W language; N/C radiation | all |
| R10 Stability | N,P,C,L,A kill-gates | all |

---

## 13. Kill / honesty list (relation-level)

| Kill | Reason |
|------|--------|
| Local GRIN as R4 | Dead (v76) |
| \(\psi\equiv\Phi_{\mathrm{EM}}\) | Dead (v77 dual independence) |
| Flavor alone as neutron | Dead (v75; residual Q) |
| Same-fabric atom (no opposite L) | Dead (v74 rings ≠ shells) |
| A=12 free park at g=0.05 | Blocked by \(Q_{\max}\); use Z-carbon or lower g |
| Monist sandbox as particle proof | Vocabulary only (WORLD residual R-SCP) |
| Kernel edit without auth | Policy — ask human |

---

## 14. Freeze declaration

**CP-W Relations freeze v0:** R1–R10 assigned as above. Downstream phases (N→P→C→L→A→U) **update evidence**, not the relation identities, without a new W revision stamp.

**World primitive stack:** see [`work/W/world_freeze_v0.md`](work/W/world_freeze_v0.md).

**Handoff:** [`work/W/FOR_U.md`](work/W/FOR_U.md).

---

## 15. References (absolute paths)

| Doc | Path |
|-----|------|
| Goals | `/home/d/code/scp/v78/GOALS.md` |
| Campaign map | `/home/d/code/scp/v78/CAMPAIGN_MAP.md` |
| v76 convergence | `/home/d/code/scp/v76/CONVERGENCE.md` |
| v76 theory | `/home/d/code/scp/v76/work/A/THEORETICAL_PACKAGE_v1.md` |
| v77 phase2 | `/home/d/code/scp/v77/PHASE2_CONVERGENCE.md` |
| v77 WORLD | `/home/d/code/scp/v77/work/TU/WORLD_v1.md` |
| v77 unified | `/home/d/code/scp/v77/work/TU/UNIFIED_PACKAGE_v1.md` |
| CONCEPT | `/home/d/code/scp/CONCEPT.md` |
| v75 state | `/home/d/code/scp/v75/STATE.md` |
| v74 carbon | `/home/d/code/scp/v74/CARBON_MAP.md`, `RESULTS.md` |
