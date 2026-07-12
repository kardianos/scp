# v75 Proposal — Multi-Fabric Architecture toward Atoms

**Date**: 2026-07-12  
**Status**: proposal (design only — no kernel change without explicit authorization)  
**Prior**: v74 Stage 2A — Z-carbon liquid-drop parks; free A=12 evaporates; Q-ball
≠ atom (missing light opposite-charge sector and multi-scale binding)  
**Standing goal** (CLAUDE.md): carbon atom from fabric — Stages 3–4 blocked on
electron sector

---

## 0. Why multi-fabric

v74 closed a structural fact: **one same-sign Q-ball fabric produces nuclei
(liquid-drop charge droplets), not atoms.** Multi-nucleon seeds fuse; late ρ²
rings are radiation, not K/L shells.

An atom needs:

1. Heavy compact core (nuclear fabric — largely in hand),
2. Light opposite-charge matter (leptonic fabric — missing),
3. Long-range binding that does **not** force short-range fusion of (1)+(2),
4. Scale separation \(m_e \ll M_N\), Bohr radius ≫ nuclear radius.

**Multi-fabric** = several full 3D field systems on the **same** space, with
**selective engagement** (sparse cross-couplings). Primary design path for v75:
**three fabrics first.**

---

## 1. Mathematical skeleton (all options)

\(S\) fabrics over shared \(\mathbb{R}^{1,3}\):

\[
\Phi^{(s)}(x),\quad s = 1,\ldots,S
\]

each \(\Phi^{(s)}\) may be multi-component (Cosserat-like). Lagrangian:

\[
\mathcal{L}
=
\sum_{s}
\mathcal{L}_s\!\big[\Phi^{(s)}, A^{(s)}\big]
\;+\;
\sum_{s < t}
\mathcal{L}_{st}\!\big[\Phi^{(s)},\Phi^{(t)}, A^{\mathrm{shared}}\big].
\]

- \(\mathcal{L}_s\): internal dynamics (own masses, potentials, optional gauges).
- \(\mathcal{L}_{st}\): **only** allowed engagement (the design content).
- Absent \(\mathcal{L}_{st}\) ⇒ complete decoupling.

Engagement graph \(K_{st}\in\{0,1\}\) or continuous strengths \(\varepsilon_{st}\):

\[
\mathcal{L}_{\mathrm{engage}}
=
\sum_{s,t}
K_{st}\,
\Big(
\eta_{st}\,\Theta^{(s)}\!\cdot\!(\nabla\!\times\!\Phi^{(t)})
+
\lambda_{st}\,V_{\mathrm{mix}}(|\Phi^{(s)}|^2,|\Phi^{(t)}|^2)
+
g_{st}\,J^{(s)}\!\cdot\! A^{(t)}
\Big).
\]

---

## 2. Option catalog (all listed paths)

### Option A — Two matter fabrics + shared U(1)  [minimal atom]

| Piece | Content | Role |
|-------|---------|------|
| Fabric **N** | Current complex Cosserat + short-range binding | Heavy Q-balls / nuclei |
| Fabric **e** | Light complex scalar(s) | Electron analogs |
| Shared **A** | Massless abelian gauge | Coulomb binding |

- Charges: \(g_N g_e < 0\).
- **No** strong \(V_{Ne}\) fusion potential.
- Smallest model that can host positronium → hydrogenoid atoms.

**Pros:** minimal fields, maps to nucleus+electron+photon.  
**Cons:** nuclear binding and “color” not separated from charge-carrying matter.

---

### Option B — Three fabrics: Bind / Charge / Light  **[v75 PRIMARY]**

| Fabric | Name | Couples to | Role |
|--------|------|------------|------|
| **C** | Color / binding | Itself (triple-product \(s_C\)); weak portal to Q | Quark-like binding medium |
| **Q** | Charge / baryonic | Shared \(A\); optional mix with C | Charged nuclear matter |
| **L** | Light / lepton | Shared \(A\) only (opposite \(g\)); **not** \(s_C\) | Electrons |

\[
\mathcal{L}
=
\mathcal{L}_C[\Phi_C]
+
\mathcal{L}_Q[\Phi_Q,A]
+
\mathcal{L}_L[\Phi_L,A]
+
\varepsilon_{CQ}\,\mathcal{L}_{\mathrm{mix}}[\Phi_C,\Phi_Q].
\]

**Selective rule:** \(\Phi_L\) never enters \(s_C\) ⇒ light lumps cannot become
nuclear “quarks.” Engagement graph:

```
   C ←ε→ Q ←g→ A ←g→ L
         (heavy)         (light, opposite charge)
```

No direct C–L or fusion Q–L short-range channel.

**Pros:** clean separation of binding vs EM vs leptons; matches three-media
poetry; room to keep v66–v74 nuclear results inside C+Q.  
**Cons:** more parameters; co-design of three sectors required.

**v75 works this option first.**

---

### Option C — Three parallel Cosserat copies (symmetric triple)

\[
\Phi^{(s)}_a,\ \Theta^{(s)}_a,\quad s=1,2,3,\ a=0,1,2
\]

each a full 12-real-field Cosserat block; engagement via sparse \(K_{st}\).

Atom-friendly assignment still needed:

- \(s=1\): nuclear (heavy \(m\), self-binding \(V\), charge \(+\)),
- \(s=2\): mediator / torsion–EM hybrid,
- \(s=3\): leptonic (light \(m\), charge \(-\)),
- \(K_{13}=0\), \(K_{12},K_{23}\neq 0\).

**Pros:** uniform code path (three instances of one kernel block).  
**Cons:** if masses/charges/potentials stay equal, **no atoms** — only richer
same-sign mergers. Symmetry must be **broken by design**, not assumed.

---

### Option D — Nuclear fabric + dark/neutral fabric + EM

Third fabric is **neutral** (no \(A\) charge): torsion bath, dark sector, or
gravity-analog clock field. Atom uses A+B only; third is optional environment.

**Pros:** exploratory; gravity/clock program (FUTURE).  
**Cons:** does not by itself supply electrons.

---

### Option E — Single fabric, two representations (not multi-fabric)

One spacetime field content, two charge representations (e.g. \(\Phi\) charge
\(+1\), \(\psi\) charge \(-1\)) under one \(A\). Mathematically a limit of Option
A with shared kinetic structure.

**Pros:** smallest code delta.  
**Cons:** harder to enforce “L never feels nuclear \(V\)” without careful
potential design; less modular co-design.

---

### Option F — Extra dimensions / branelike separation

Fabrics on different leaves of a product space with small mixing. Deferred —
overkill before Options A/B fail.

---

### Comparison (atom capability)

| Option | Can host atom? | v75 priority |
|--------|----------------|--------------|
| **B (C/Q/L)** | Yes, cleanest modular | **1 — design + toy math** |
| A (N+e+A) | Yes, minimal | 2 — implement after B spec freezes |
| C (triple Cosserat) | Yes only if asymmetric | 3 — implementation variant of B |
| E (two reps) | Yes | parallel simplification of A |
| D (neutral third) | Partial | later |
| F (extra dim) | Maybe | not now |

---

## 3. v75 primary: three-fabric (C, Q, L) target model

### 3.1 Field content (proposal)

**Fabric C — binding (color-analog)**  
- Complex triple \(\Phi^C_a\), \(a=0,1,2\)  
- Potential depends only on \(s_C = \prod_a |\Phi^C_a|^2\) (confine-like)  
- **Ungauged** or gauged only under a *global* symmetry not shared with L  
- Optional massive \(\Theta^C\) with \(\eta_C\) curl (internal only)

**Fabric Q — nuclear charge**  
- Complex \(\Phi^Q\) (one or three components)  
- Couples to shared \(A_\mu\) with charge \(g_Q > 0\)  
- Mixes to C via \(\varepsilon_{CQ}\) portal so “baryons” are C+Q composites  
- Limit \(\varepsilon_{CQ}\to\infty\) or algebraic constraint can recover today’s
  single-ball identification “components = quarks + charge”

**Fabric L — light leptons**  
- Complex \(\psi\) (or \(\psi_a\))  
- Couples to **same** \(A_\mu\) with \(g_L < 0\), \(|g_L|\) comparable to \(|g_Q|\)  
- Own mass \(m_L\) and potential \(V_L\) with **window for light Q-balls**  
- **Forbidden:** any term placing \(\psi\) inside \(s_C\) or a strong \(V(s_C,|\psi|^2)\)

**Shared gauge**  
- \(A_\mu\), \(E_i\) as today (kernel-v3 abelian sector)  
- Gauss law: \(\nabla\cdot E = g_Q\rho_Q + g_L\rho_L\)

### 3.2 Engagement map (what is on / off)

| Channel | On? | Form | Purpose |
|---------|-----|------|---------|
| C self-binding | **ON** | \(V_C(s_C)\) | nuclear cohesion / confinement analog |
| Q self-mass + kinetic | **ON** | standard | charged nuclear kinetic |
| L self-mass + kinetic | **ON** | light parameters | electron Q-balls |
| C–Q portal | **ON (weak–moderate)** | \(\varepsilon_{CQ}\) mix or shared components | baryon = bound C+Q |
| Q–A minimal | **ON** | \(D=\partial+ig_Q A\) | nuclear Coulomb |
| L–A minimal | **ON** | \(D=\partial+ig_L A\) | electronic Coulomb |
| C–A | **OFF or tiny** | — | color not long-range EM |
| C–L direct | **OFF** | — | no nuclear–lepton fusion bag |
| Q–L direct potential | **OFF** | — | only EM binds N to e |
| L–L short-range | optional weak | — | positronium structure |

### 3.3 Success criteria (structural atom)

1. Isolated N-ball (C+Q composite) stable as in v72/v74.  
2. Isolated light L-ball stable, \(M_L \ll M_N\).  
3. Opposite-charge pair (positronium analog): long-lived orbit or bound ground
   state, \(E < M_N+M_L\) if applicable / \(E < 2M_L\) for e⁺e⁻.  
4. **Hydrogenoid:** one N + one L, separation \(D \gg r_N\), no merger for
   \(T \ge 10^3\) nuclear clocks.  
5. Optional: multi-L “shell” density peaks in \(|\psi|^2\) (true K/L analog).  
6. Bookkeeping: separate Noether charges \(Q_C\) (if any), \(Q_Q\), \(Q_L\);
   total Gauss flux \(= g_Q Q_Q + g_L Q_L\); fabric ledgers nearly closed per
   sector (v73).

### 3.4 Explicit non-goals for v75

- Quantitative chemistry / real MeV–eV matching  
- Fermions / Pauli (classical shells first)  
- Kernel implementation until user authorizes  
- Gravity fabric (Option D later)

---

## 4. First-principle concepts for selective engagement

How to **derive or constrain** who engages whom — not arbitrary knobs.

### P1 — Symmetry and representations (primary)

Engagement is allowed only when fields transform so that \(\mathcal{L}_{st}\) is
a **scalar under the full symmetry group**.

- Shared continuous symmetry \(G\) (e.g. diagonal U(1)\(_\mathrm{EM}\)): only
  charged operators couple to \(A\).
- Fabric-private symmetries \(G_s\): if \(G_C\) is not shared with L, **no
  invariant** local C–L coupling exists ⇒ engagement forbidden by symmetry
  (strongest selection rule).
- **Charge lattice:** assign \((q_C, q_Q, q_L)\); allowed operators must have
  total charge 0 under every gauged U(1).

**Design rule:** put nuclear binding under a **private** symmetry of C (or C+Q);
put EM under a **shared** U(1) with \(q_Q = -q_L\).

### P2 — Conservation laws / Noether sectors

Each fabric’s continuous phase symmetry yields a charge \(Q_s\). Engagement
that **violates** a desired conservation law is disallowed.

- Want \(Q_L\) conserved separately from \(Q_Q\) ⇒ no \(\Phi_Q^\dagger\psi\) mass
  mixing (or only extremely weak portal).
- Want baryon-like conservation ⇒ C+Q portal may conserve \(B = Q_Q\) or
  \(B = Q_C+Q_Q\) but not convert to L.

**Design rule:** list conserved charges first; write only operators that respect
them. Engagement graph = dual of conservation list.

### P3 — Scale separation and decoupling (EFT)

If fabric \(s\) has characteristic mass \(m_s\) and \(m_L \ll m_N\), then at
electronic energies the nuclear fabric is **integrated out** up to multipoles:

\[
\mathcal{L}_{\mathrm{eff}}[\psi,A]
=
\mathcal{L}_L
+
\mathcal{L}_A
+
\sum_k
\frac{c_k}{M_N^{k}}\,O_k[\psi,A,J_N^{\mathrm{ext}}].
\]

Short-range N–e fusion operators are **irrelevant** (suppressed by powers of
\(1/M_N\)) if the UV completion doesn’t include a light shared bag field.

**Design rule:** choose \(m_L/m_N\) and omit light shared scalars that mediate
fusion; remaining engagement is Coulomb (+ weak higher multipoles).

### P4 — Stability and energy positivity

Cross terms must not destroy:

- positive energy (or bounded-below energy) of each sector,
- VK / fixed-Q existence for N and for L separately,
- Derrick balance (complex phases / charge pressure still required per fabric).

Forbidden examples: wrong-sign portals that run away; mixings that open a
radiation channel with \(\omega > m\) in the light sector.

**Design rule:** check quadratic fluctuation spectrum about vacuum and about
each isolated soliton after enabling each \(\varepsilon_{st}\).

### P5 — Process ledgers (v73)

Each fabric has energy flux \(S^{(s)}\) and current \(J^{(s)}\):

\[
\partial_t e^{(s)}+\nabla\!\cdot S^{(s)}
=
\sum_t \mathcal{P}_{st},\qquad
\mathcal{P}_{st}=-\mathcal{P}_{ts}.
\]

Selective engagement ⇔ sparse \(\mathcal{P}_{st}\). An atom is a **composite
process**: nuclear ledger + leptonic ledger exchange power only through the EM
Poynting channel, each nearly closed on its own timescale.

**Design rule:** engagement channels are exactly the allowed \(\mathcal{P}_{st}\)
ports; measure cross-flux in simulation as an engagement diagnostic.

### P6 — Locality and causality

Only local (or exponential Yukawa) operators; shared massless \(A\) is the
unique long-range engagement unless a second massless field is introduced.

**Design rule:** long-range selective engagement = shared gauge; short-range =
contact/portal with mass \(\ge \min(m_s,m_t)\).

### P7 — Minimal coupling uniqueness

Once charges under \(A\) are assigned, the covariant derivative is fixed
(\(D=\partial+iqA\)). No free “engagement shape” for EM — only \(q_s\).

**Design rule:** settle charge assignment table first; EM engagement follows.

### P8 — Discrete selection rules

Parity, charge conjugation, fabric-exchange \(\mathbb{Z}_2\) (if any) can forbid
portals. Example: if L is odd under a dark parity that C is even under, C–L
scalars are forbidden.

### P9 — Topological / winding selection (optional)

If fabric C supports windings (v73 Q-ring) and L does not, engagement can
couple only to charge density, not to winding, preventing light fabric from
unwinding nuclear spin carriers.

---

## 5. Co-design of each fabric “in phase”

Multi-fabric atoms fail if sectors are designed in isolation with incompatible
clocks, lengths, or ansätze. **Co-design in phase** means simultaneous
constraints on \((\omega_s, Q_s, g_s, m_s)\) so a joint stationary (or
Born–Oppenheimer) state exists.

### C1 — Shared spacetime and common rest frame

All fabrics use the same coordinates and a common rest frame for the atom.
Boosts act diagonally on the whole multiplet.

### C2 — Commensurate clocks (phase locking)

Rotating ansätze:

\[
\Phi^{(s)}(x,t)
=
f^{(s)}(x)\,e^{i\omega_s t}
\quad\text{(or }\omega_s t + n_s\phi\text{)}.
\]

**In phase** does **not** require \(\omega_C=\omega_Q=\omega_L\). It requires:

- each \(\omega_s\) inside its fabric’s existence window,
- gauge-invariant combinations consistent with static \(|\Phi^{(s)}|\) and static
  \(A_0\) when seeking stationary atoms,
- for bound orbits, **commensurability** or adiabatic separation (C4).

Gauge-covariant local frequency:

\[
\omega_s^{\mathrm{eff}}(x)
=
\omega_s + q_s A_0(x)
\]

must solve the coupled radial system together with Gauss’s law.

### C3 — Joint variational principle (fixed multi-charge)

Generalize v72 fixed-Q flow to several charges:

\[
E[\{f^{(s)}\},A]
+
\sum_s
\frac{Q_s^2}{2N_s[f^{(s)}]}
\;\rightarrow\;
\min
\quad\text{at fixed }\{Q_s\}.
\]

Lagrange multipliers \(\omega_s = Q_s/N_s\). Engagement enters only through
shared \(A\) and allowed portals. **Co-design** = one optimization over all
fabrics, not sequential independent Q-balls glued by hand.

### C4 — Born–Oppenheimer (nuclear frozen, light relaxes)

Because \(M_N \gg m_e\):

1. Solve nuclear fabric (C+Q) → fixed \(\rho_N(x)\), multipoles.  
2. Solve light fabric + \(A\) in that background.  
3. Optional back-reaction iterate.

Phase relation: nuclear clock may be averaged; electronic cloud sees time-averaged
Coulomb field. This is the first practical **co-design algorithm** for
simulation without a full joint 3-fabric relaxer.

### C5 — Length co-design

| Scale | Set by | Target hierarchy |
|-------|--------|------------------|
| Nuclear radius \(r_N\) | \(m_N\), \(V_C\), \(Q_Q\) | \(\sim\) few lattice units |
| Compton / light radius \(r_L\) | \(m_L\), \(V_L\), \(Q_L\) | \(r_L \ll r_N\) for pointlike e, or moderate |
| Bohr radius \(a_0\) | \(\alpha = g^2/4\pi\), reduced mass | \(a_0 \gg r_N\) |
| Box \(L\) | numerics | \(L \gtrsim 4a_0\) |

**In phase** includes choosing \((g, m_L, m_N, Q_L, Q_Q)\) so these scales fit
one grid.

### C6 — Charge co-design

- Neutrality for true atoms: \(g_Q Q_Q + g_L Q_L = 0\) (classical hydrogenoid:
  \(Q_L = - (g_Q/g_L) Q_Q\)).
- For positronium: \(Q\) and \(-Q\) in L only, N absent.
- Integer-like quantization later via \(\hbar_{\mathrm{eff}}=Q\) per fabric
  (CONCEPT §7) — optional v75+ story.

### C7 — Adiabatic engagement ramp

Turn on \(\varepsilon_{st}\) and \(g\) slowly from decoupled fabrics (each
pre-solved) so the composite tracks the joint branch — multi-fabric analog of
η-dressing (v72).

### C8 — Diagnostics for “in phase”

- Relative phase drift \(\dot\delta_{st} = \omega_s-\omega_t\) in overlap region  
- Cross-ledger power \(\mathcal{P}_{st}\) time-averaged to ~0 for stationary atom  
- Gauss law residual with **both** charge densities  
- No secular growth of fusion contact density \(\int |\Phi_Q|^2|\psi|^2\)

---

## 6. Simulation fidelity ladder (L0–L3) — CRITICAL

Hardware (larger than V100) and multi-adaptive / nested grids change **cost**,
not **when** reduced models lie. Geometry-sensitive claims need 3D.

| Layer | What | Authoritative for | Not authoritative for |
|-------|------|-------------------|------------------------|
| **L0** | Full **3D** multi-field PDE (gauged kernel; later multi-fabric) | Engagement discovery; fusion vs orbit; winding/spin; multi-center | — (truth layer) |
| **L1** | **Nested / multi-adaptive 3D** grids (see `sfa/MULTI_RESOLUTION.md`) | Atom hierarchy \(r_N \ll a_0 \ll L\); long-range EM + fine cores | New physics not in L0 |
| **L2** | Radial / 1D / BO **shooters** | Spherical existence windows; order-of-magnitude \(E_b(g,m)\) | Knots, spin, impact parameter, multi-center standoff, multipole radiation |
| **L3** | **Kinematic n-bar** (point/rigid bodies + allowed force graph) | Many-body population, “chemistry-scale” constructs **after** L0 forces known | Inventing engagement; replacing topology/ledgers |

**Rules**

1. Selective engagement is defined in the **Lagrangian (L0 design)**, not as an
   n-bar UI switch.
2. L2 shooters are **design thermometers** — stamp every result
   “not topology-safe.” Knots can slip/bind wrong when angular structure is
   averaged away.
3. L3 n-bar is **not a distraction** only after L0 shows: stable N, stable L
   (or ±Q), Coulomb \(\sim 1/r^2\), **no** N–L fusion. Then n-bar *inherits*
   the engagement graph for multi-electron / abundance runs.
4. Larger GPUs + adaptive grids → invest in **L0/L1 3D**, not more L2 alone.

**Effective L3 Hamiltonian (only after L0 validation):**

\[
H_{\mathrm{eff}}
=
\sum_i \sqrt{p_i^2+M_i^2}
+
\sum_{i<j}
V_{ij}^{\mathrm{allowed}}(r_{ij},\mathrm{phase}_{ij},\mathrm{fabric}_{ij})
\]

with e.g. \(V_{NL}^{\mathrm{fusion}}\equiv 0\) *because PDE said so*, not by hope.

---

## 7. Suggested v75 work plan (fidelity-aware)

Detail and checklists: **`v75/FIRST_STEPS.md`**. Setup/findings log:
**`v75/FINDINGS.md`**.

### Phase 0 — Design freeze (docs)

- [x] Option catalog + three-fabric primary  
- [x] Fidelity ladder L0–L3  
- [x] First principles for engagement + co-design in phase  
- [ ] Freeze Option B charge/symmetry table as non-negotiable  
- [ ] Parameter sketch \(m_C,m_Q,m_L,g_Q,g_L,\varepsilon_{CQ}\)

### Phase 1 — L0 3D with **existing kernel** (no kernel edit)

Opposite-charge already seedable (`gen_qball_multi`: negative \(\omega\) = opposite
charge). Same fabric, two signs = Option **E-lite** probe of EM engagement.

1. **±Q pair at rest** (D large): measure force sign/exponent (attraction).  
2. **±Q with tangential velocity**: orbit vs annihilation (positronium analog).  
3. **Document nulls** (slow annihilation known CONCEPT §5) as bounds on “atom
   without multi-fabric.”

### Phase 2 — L2 radial scouts (parallel, labeled)

Spherical BO / light-ball windows only — **not** atom proof. Deliverable
`v75/RADIAL_SCOUT.md` if pursued.

### Phase 3 — L1 nested 3D (large GPU)

Once ± binding is clean: fine core + coarse cloud, \(L \gtrsim 4a_0\). Use
nested-grid design; rent hardware beyond V100 as needed.

### Phase 4 — Multi-fabric kernel (needs **explicit user authorization**)

True C/Q/L field split only if Phase 1 shows EM binding is real but same-fabric
annihilation/fusion blocks atoms.

### Phase 5 — L3 n-bar (optional)

Multi-electron / carbon-atom packing only after L0 engagement graph is measured.

### Phase 6 — Carbon atom goal

v74 Z-carbon (or A=12 at lower \(g\)) + 6 light opposite charges under validated
engagement.

---

## 8. Risks and null results that would kill Option B

| Null | Meaning |
|------|---------|
| No light Q-ball window at any \(m_L,V_L\) | Electron sector impossible in this potential class |
| Opposite charges always annihilate / merge | Need stronger selection or fermionic exclusion later |
| Coulomb binds but always captures into one hybrid blob | \(V\) contact too strong; must kill Q–L potential |
| Scales won’t fit on V100 grids | Need coarser nuclear / finer electronic multi-grid |
| Portal \(\varepsilon_{CQ}\) destabilizes nuclear branch | Recover v74 limit only at \(\varepsilon\to\) special value |

Any of these is publishable as a boundary on multi-fabric atoms.

---

## 8. Relation to standing goal and docs

| Doc | Role after v75 starts |
|-----|------------------------|
| CLAUDE.md standing goal | Still carbon atom; Stage 3 = multi-fabric light sector |
| v74/RESULTS.md | Nuclear-only baseline (Z-carbon, A=12 control) |
| v75/PROPOSAL.md | This file — architecture |
| FUTURE.md | Positronium / electron items become v75 Phase 3 |
| CONCEPT.md | Update only when a fabric split is **measured**, not at proposal |

---

## 9. Decision requested from user

1. **Freeze Option B (C/Q/L)** as the design target?  
2. Proceed Phase 1 **radial Born–Oppenheimer atom** (no kernel change)?  
3. Or prefer **Option A minimal** (N+e+A) for faster path to 3D?

Default if uninstructed: **Option B design freeze + Phase 1 radial math**, no
kernel edits.

---

## 10. One-page summary

- **Problem:** one fabric ⇒ nuclei, not atoms (v74).  
- **Proposal:** multi-fabric with selective engagement.  
- **All options:** A (N+e+A), **B (C/Q/L primary)**, C (triple Cosserat), D
  (neutral third), E (two reps), F (extra dim).  
- **Engagement from principles:** symmetry/reps, conservation, EFT decoupling,
  stability, ledgers, locality, minimal coupling, discrete/topology rules.  
- **Co-design in phase:** joint fixed-multi-Q variational principle, commensurate
  / BO clocks, scale and charge tables, adiabatic ramp, cross-ledger diagnostics.  
- **First work:** radial multi-sector atom ODEs; kernel only with authorization.
