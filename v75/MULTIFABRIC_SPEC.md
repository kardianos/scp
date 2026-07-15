# v75 MULTIFABRIC SPEC — Option B staged integration

**Date**: 2026-07-13  
**Status**: **B(1) Shape β CPU implemented** (2026-07-13); GPU port pending  
**Parent**: `PROPOSAL.md` (Option B primary), `FINDINGS.md` F8–F10  
**Authorization**: multi-fabric kernel authorized 2026-07-12; user re-confirmed
2026-07-13. `sfa.h` unchanged (L columns use names `lphi_*` / `lth_*`).

---

## 0. Decision (locked)

| Item | Choice |
|------|--------|
| Architecture | **Option B** — fabrics C / Q / L + shared massless \(A\) |
| First ship | **B(1) Shape β** — full three-array layout; C and Q **locked** |
| Follow-on | **B(2)** — unlock independent Q + \(\varepsilon_{CQ}\) |
| Experimentation | **Outcome-driven sequence** (gates below); branch if a gate fails or
  a cheaper probe is enough |

**Not chosen for first ship:** pure Option A (two fabrics only), Shape α
(aliases only), Shape γ (Q absent). Those remain fallbacks if β is too heavy
on memory or schedule (see §6).

---

## 1. Physics skeleton (Option B)

### 1.1 Fabrics

| Fabric | Role | Bag \(V_t(s)\) | Couples to \(A\) | Notes |
|--------|------|---------------|------------------|--------|
| **C** | Binding / “color” medium | **Yes** — \(s_C=\prod_a\|\Phi_C^a\|^2\) | optional \(q_C\) (B1: often 0) | Private nuclear bag |
| **Q** | Charge / baryonic matter | **No** (B2 ideal); B1: locked to C | **Yes** \(q_Q=+g\) | Charge carrier |
| **L** | Light / lepton | **Never** | **Yes** \(q_L=-g\) | Opposite charge |

Shared: one massless abelian **\(A\)** (existing gauged sector).

Engagement graph (target):

```
   C ←ε_CQ→ Q ←g→ A ←g→ L
         (heavy)         (light)
```

**Hard off-list (all stages):**

- L never enters any nuclear bag product  
- No C–L or Q–L short-range fusion / mass-mix portal  
- No second massless gauge in v1  

### 1.2 Charge table (freeze for B1; refine in B2)

| Stage | \(q_C\) | \(q_Q\) | \(q_L\) | Comment |
|-------|--------|--------|--------|---------|
| **B(1) β** | 0 | \(+g\) | \(-g\) | All EM charge on Q; C neutral binder. With lock \(\Phi_Q=\Phi_C\), effective heavy charge is carried by the Q copy of the same profile. |
| **B(2) option** | 0 | \(+g\) | \(-g\) | Preferred long-term (bag ≠ charge) |
| **B(2) alt** | \(+g\) | 0 or small | \(-g\) | If Q alone cannot hold lumps |

Default coupling magnitude: `g_gauge = 0.05` (same package as v69+).

### 1.3 B(1) Shape β lock

**Allocate** full independent storage for \(\Phi_C,\Phi_Q,\Phi_L\) (+ θ sectors per
existing complex Cosserat) and one \(A\).

**Dynamics constraint (B1):**

\[
\Phi_Q(x,t) \equiv \Phi_C(x,t)
\quad\text{(and likewise velocities / θ if present)}
\]

Implementation options (pick in code review; prefer the one that keeps G0 easy):

| Mode | Mechanism | Use |
|------|-----------|-----|
| **β-copy** | Evolve C fully; each substep set Q ← C (and \(\dot Q\leftarrow\dot C\)) | Simple; Q arrays “live” for I/O |
| **β-dual-step** | Evolve both with identical force from \(s_C\) only; charge current from Q only | Stresses dual paths; must stay locked |
| **β-init-lock** | Identical seeds only; free evolution (will drift) | **Not for B1 science** — only for drift diagnostics |

**B1 science uses β-copy or β-dual-step.** Bag force from \(s_C\) only; gauge
current from Q (and L); L has its own bag **off** or own light potential later
(v1: same \(V_t\) form on L **alone** so L can be a Q-ball, never mixed into \(s_C\)).

### 1.4 L potential in B1

For the first atom gate, L is an **E-lite-class light Q-ball**:

- Own product \(s_L=\prod_a\|\Phi_L^a\|^2\) and \(V_t(s_L)\) — self-binding only  
- **Not** mixed into \(s_C\)  
- Parameters: start with nuclear-class light ω=1.46 package; hierarchy later  

This matches F8–F10 seeds moved onto two matter sectors (heavy locked C=Q, light L).

---

## 2. Config surface (proposed)

All new; defaults preserve today’s single-fabric path.

```
# Multi-fabric (proposed names — finalize at implement)
n_fabrics        = 1     # 1 = legacy single (G0); 3 = C,Q,L
mf_stage         = 0     # 0=off/legacy, 1=B1 lock, 2=B2 free Q
eps_CQ           = 0.0   # portal; B1 ignored or forced 0
q_C              = 0.0
q_Q              = 1.0   # × g_gauge for physical charge
q_L              = -1.0
mf_lock_CQ       = 1     # B1: force Φ_Q=Φ_C; B2: 0
# existing package still applies:
complex_phi=1, complex_gauge=1, g_gauge=0.05, m_theta=1.6, eta=0, …
```

**G0 regression:** `n_fabrics=1` (or `mf_stage=0`) must match current kernel
behavior (same arrays, same forces) within existing bit/ULP policy (v69 style).

---

## 3. Gate sequence (outcome-driven)

```
                    ┌─ G0 PASS ──► B1 β implementation complete
                    │
B1 implement ───────┤
                    └─ G0 FAIL ──► fix plumbing; do not proceed

G0  n_fabrics=1 regression (byte/ULP vs current)

G1  Heavy+L, g=0 or A inert: no Coulomb, no cross-bag force
    FAIL → engagement leak (L seeing s_C or shared bag bug)

G2  ± large D (F8 redo on two fabrics): attract / 2 clusters
    FAIL → gauge/charge assignment bug

G3  Head-on contact (F10 redo): MUST NOT merge ★
    PASS → atom path open; proceed to hydrogenoid / lighter L
    FAIL → bag isolation broken; do not unlock B2

G4  B2 prep: mf_lock_CQ=0, eps_CQ=0 — Q alone can/can’t hold a lump?
    outcome branches §6

G5  Ramp eps_CQ (C7 / outer tune) — C+Q nuclear composite

G6  Repeat G3 under true C≠Q — still no L merge
```

### 3.1 Experimental flexibility (user intent)

Sequence is **not** rigid beyond G0→G3 for B1:

| If … | Then … |
|------|--------|
| Memory tight on GPU for 3 full copies | Fall back Shape **α** (alias Q=C) for G1–G3 only; keep β on CPU |
| G3 PASS early on small N | Skip long orbit; do light hierarchy next |
| G3 FAIL | Debug isolation; optional probe: zero L bag and pure Coulomb cloud |
| G4 shows Q alone unstable | Keep charge on C for B2-alt; ε portal binds Q into C |
| G4 shows Q alone stable | Prefer q_C=0, bag on C only, ε weak–moderate |
| Self-tuning desired | Only after G3; float `{eps_CQ, m_L/m_N}` under cost (§5) — not before |

---

## 4. Validation seeds (B1)

Reuse E-lite tooling, dual-fabric write:

| Gate | Seed sketch |
|------|-------------|
| G2 | Heavy (C=Q lock) ω=+1.46 @ +D/2; L **ω=+1.46** (same sign) @ −D/2; rest |
| G3 | Same + radial closing boosts (head-on). **Not** E-lite ±ω — fabric q already opposite |
| Optional | Tangential orbit at vt≈0.0154 after G3 |

Generators: extend `gen_qball_pair_boost` / multi to **per-fabric** columns once
SFA layout for multi-fabric is defined (prefer: fabric tag or column name
prefix `c_`, `q_`, `l_` — **no sfa.h change** if names fit 11-char rule via
existing component scheme; else ask before format change).

---

## 5. Self-tuning (follow-on only)

Not required for B1.

After G3 PASS, optional outer loop on soft θ:

\[
\theta \in \{\varepsilon_{CQ},\, m_L/m_N,\, g_L/g_Q\}
\]

Cost ingredients: merge penalty, F8 \(a_C\) ratio, nuclear branch vs v74,
stability, non-EM cross-ledger power.  
**Hard zeros stay hard** (L∉bag). Tuner must not reintroduce off-list operators.

---

## 6. Branch table (Shape / stage)

| Outcome | Action |
|---------|--------|
| β too heavy (RAM) | α for science gates; β retained as “full layout” CI on small N |
| G1–G3 all PASS | **Default next:** lighter L co-design **or** B2 G4 — choose by carbon priority (hierarchy vs nuclear split) |
| Need carbon packing ASAP after G3 | Multi-L + Z-carbon heavy template; delay B2 |
| Need correct nuclear substructure | B2 G4–G6 before multi-electron |
| ε_CQ → 0 under tuning | Report “B collapsed toward A physics”; keep β layout |
| ε_CQ finite stable | Full B narrative validated |

---

## 7. Implementation phases (kernel)

### Phase B1 (this spec)

1. Config + allocation for 3 matter copies + 1 A when `n_fabrics=3`  
2. Force: bag from C only; L self-bag only; currents Q+L → A  
3. Lock Q←C each step (β-copy)  
4. Diag: per-fabric Q, E pieces; `gauss_max`; optional cross-ledger  
5. G0–G3 campaigns  

### Phase B2 (separate PR)

1. `mf_lock_CQ=0`  
2. Portal \(\varepsilon_{CQ}\) operator (specify continuum + lattice in sub-spec)  
3. G4–G6 + optional self-tune  

### Explicit non-goals for B1

- Nested grids  
- Full C/Q/L three-way η  
- Multi-electron shells  
- sfa.h redesign  
- Dropping single-fabric path  

---

## 8. Success criteria

**B1 done when:**

- [x] CPU kernel: `n_fabrics=3`, `mf_lock_CQ`, private bags, shared A, Q←C  
- [x] Seed tool: `bin/gen_mf_pair_boost`  
- [x] G0 smoke: `n_fabrics=1` runs; MF rest T=5 stable gauss floor  
- [x] GPU port of B1 (`scp_sim.cu`); smoke N=48 COMPLETE, gauss floor held  
- [ ] G0 formal ULP regression vs pre-MF binary (optional)  
- [ ] G1 PASS  
- [ ] G2 PASS (signs + 2 clusters)  
- [ ] **G3 PASS** (head-on no merge; E drift / s_max no fusion signature)  
- [ ] FINDINGS F11 documents B1 gates  

**B2 start when:** B1 G3 done **and** user prioritizes nuclear split over multi-L carbon packing.

### 8.1 Implementation notes (2026-07-13)

| Item | Detail |
|------|--------|
| Config | `n_fabrics`, `mf_lock_CQ`, `mf_stage`, `q_C/q_Q/q_L`, `eps_CQ`, `init_sfa_L` |
| CPU | `scp_sim.c`: C primary arrays; Q+L 36 blocks each; `compute_forces_complex_gauge_mf` |
| Energy | C + L matter; E_em once; Q not double-counted when locked |
| Gauss | `em_rho = q_C ρ_C + q_Q ρ_Q + q_L ρ_L` |
| Output | 54 cols when MF+gauge (30 + 24 L) |
| GPU | **Ported** in `scp_sim.cu` (same β semantics; needs 24+ GB for N=192 MF) |
| Seeds | `gen_mf_pair_boost N L D profile ωC ωL vC… vL… out_C.sfa out_L.sfa` |

---

## 9. Relation to prior E-lite

| E-lite result | B1 use |
|---------------|--------|
| F8 Coulomb signs | G2 acceptance |
| F10 orbit kinematics | Optional after G3; not blocking |
| F10 head-on merge | **G3 must invert** — primary multi-fabric justification |

Skip Step 2b multi-rev until G3 PASS (or run only as non-blocking parallel).

---

## 10. Open items before first kernel commit

1. Exact continuum portal \(\mathcal{L}_{\mathrm{mix}}\) for B2 (can stub as “TBD” in B1)  
2. SFA column naming for 3× matter without format version bump  
3. GPU: N max with 3× arrays on V100-16/32 / 4090 (document table)  
4. Confirm β-copy vs β-dual-step in code review  
5. Charge on C vs Q only — default **q_C=0, q_Q=+g** unless G4 forces alt  

---

## 11. One-line summary

> **Ship Option B layout with C and Q locked equal (Shape β); prove private-bag
> isolation against L (G3); then unlock Q and \(\varepsilon_{CQ}\) (B2) or branch
> to lighter L / carbon by outcome.**
