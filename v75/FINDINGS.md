# v75 FINDINGS — Setup, baseline, and multi-fabric design log

**Updated**: 2026-07-12  
**Docs**: `PROPOSAL.md` (architecture), `FIRST_STEPS.md` (actions), this file (record)

---

## 1. Setup (where things live)

### Code / theory baseline (pre-v75)

| Item | Location / value |
|------|------------------|
| Kernel | `sfa/sim/scp_sim.c`, `scp_sim.cu` — **do not edit without explicit auth** |
| SFA | `sfa/format/sfa.h` |
| Standard U(1) package | `complex_phi=1`, `complex_gauge=1`, `g_gauge=0.05`, `m_theta=1.6`, `eta=0` |
| Potential | \(V_t(s)=(\mu/2)s/(1+\kappa s)\), \(s=\prod\|\Phi_a\|^2\), \(m^2=2.25\), \(\mu=-41.345\), \(\kappa=50\) |
| Gauged Q-ball window | \(\omega\in(1.406,1.5)\), **\(Q_{\max}=921\)** at \(g=0.05\) |
| Branch data | `v69/theory/gscan.tsv` |
| Multi-ball seeds | `bin/gen_qball_multi` (neg. \(\omega\) = opposite charge) |
| Nested grids (design) | `sfa/MULTI_RESOLUTION.md` (fixed nesting; not full AMR yet) |
| Runner | `scp-runner` MCP — local CPU / remote GPU |
| Large data | `/space/scp/` (not git) |

### Standing goal

CLAUDE.md: **carbon atom from fabric**. Stages 0–2A (nuclear) done in v74;
Stages 3–4 need light opposite-charge + multi-scale binding → **v75 multi-fabric**.

### v74 nuclear baseline (measured)

Data: `/space/scp/v74/results/`, writeup `v74/RESULTS.md`, renders
`v74/results/{c6,c12}_light.png`.

| Run | Recipe | Seed Q | Final Q (T=300) | vs \(Q_{\max}=921\) | Morphology |
|-----|--------|--------|-----------------|---------------------|------------|
| **c6_light** | 6× light (ω=1.46), octahedron D=10 | 907.8 | **649.8** | 0.71× subcritical | **One** droplet; mid-branch E/Q +0.5%; mass defect 4.85% |
| **c12_light** | 12× light, icosahedron D=10 | 1959.9 | **1410.5** | 1.53× **super** | **One** droplet; still evaporating (dQ/dt≈−0.30) |

**Setup notes (v74 production)**

- N=192, L=36, T=300, g=0.05, η=0, m_θ=1.6, absorbing BC  
- Profiles: `v74/profiles/f_w146_g005.txt`  
- GPU: Vast V100, ~38 min/run, ~72 ms/step  
- **Cfg must set N,L matching seed** (cfg path does not auto-read SFA grid)  
- Co-phase interference inflates seed Q (+33% c6, +43% c12 vs free inventory)

### Computational environment decisions (2026-07-12)

| Constraint | Decision |
|------------|----------|
| V100 not a hard ceiling | Prefer **larger GPU** for L0/L1 atom hierarchy when needed |
| Multi-adaptive / nested grids | **In scope** for L1 (`MULTI_RESOLUTION.md`); use for \(r_N \ll a_0\) |
| Non-3D shooters | **Allowed only as L2 scouts**; never claim knot/orbit/multi-center |
| Kinematic n-bar | **L3 only after L0 engagement validated** — not a substitute fabric |

---

## 2. Findings (physics + design)

### F1 — One same-sign Q-ball fabric makes nuclei, not atoms [v74 measured]

- Multi-nucleon seeds **fuse** into a single charge droplet (liquid-drop).  
- Late concentric ρ² rings are **radiation / multipoles**, not K/L electron shells.  
- Peak diagnostics: multi-center → **one core** by t≈50 (c12: \(Q_{\mathrm{core}}/Q\) 0.18→0.93).  
- Free **A=12 always super-critical** at g=0.05 (\(12 Q_N^{\min}>Q_{\max}\)).  
- Stable parked nuclear analog at g=0.05: **Z-carbon (A=6 light)**.

### F2 — Q-ball is the wrong object to *be* the whole atom [design]

An atom needs heavy core + **light opposite-charge** cloud + long-range bind
**without** short-range fusion. Current kernel: all matter same U(1) charge
sign convention on co-rotating balls; opposite via \(\omega\to-\omega\) exists
as a seed but is not a separate fabric with private binding potential.

CONCEPT: opposite-charge pairs **annihilate slowly** (charge segregation −4×
over ~600 t.u.) — bound orbit not yet demonstrated as immortal.

### F3 — Multi-fabric is the architectural response [v75 proposal]

Primary: **three fabrics C (bind) / Q (charge) / L (light)** with sparse
engagement (PROPOSAL §2 Option B, §3).  
Minimal alternative: N+e+shared A (Option A) or two reps in one fabric (E).

### F4 — Selective engagement must come from principles, not knobs [design]

Documented in PROPOSAL §4: symmetry/reps, Noether sectors, EFT decoupling,
stability, process ledgers, locality, minimal coupling, discrete/topology.

**Strongest atom-enabling rule:** L never enters nuclear bag potential \(s_C\);
N–L long-range only through shared \(A\) with \(g_Q g_L<0\).

### F5 — Co-design must be joint “in phase” [design]

PROPOSAL §5: fixed multi-Q variational principle, Born–Oppenheimer, scale table
\(r_N \ll a_0 \ll L\), charge neutrality, adiabatic engagement ramp, cross-ledger
diagnostics.

### F6 — Fidelity ladder prevents self-deception [design 2026-07-12]

| Layer | Role |
|-------|------|
| L0 3D PDE | Truth / engagement discovery |
| L1 nested 3D | Scale separation for atoms |
| L2 radial | Parameter windows only (“not topology-safe”) |
| L3 n-bar | Higher constructs after L0 force graph known |

Kinematic selective engage/disengage of “relationships” among particles is the
**L3** picture of the same graph L0 defines in the Lagrangian. Building L3
first is a **distraction**; building it after L0 is how n-body physics always
works.

### F7 — Concrete seed capability already present [tooling]

`gen_qball_multi` documents **negative omega = opposite charge**.  
⇒ Phase 1 L0 positronium / ± force tests need **no kernel change**, only
seeds + configs + analysis. That is Option **E-lite** (same fabric, two signs),
not full C/Q/L — but it tests whether Coulomb multi-center binding can exist
before investing in multi-fabric kernel work.

---

## 3. Open nulls / risks (to record when hit)

| Risk | Status |
|------|--------|
| ±Q always annihilate / merge; no long-lived orbit | Open — Phase 1 target |
| No light-ball window under any \(V_L\) | Open — needs multi-fabric or new pot |
| Nested grid not implemented in production kernel | Design only (`MULTI_RESOLUTION.md`) |
| True three-fabric needs kernel authorization | Blocked until user says so |
| A=12 parked isotope needs lower \(g\) or lighter \(Q_N\) | Known from v74 Q-scale survey |

---

## 4. Chronology (v75 session arc)

| When | Event |
|------|-------|
| 2026-07-11 | v74 Stage 0 map + c6/c12 production; standing goal in CLAUDE |
| 2026-07-12 | c6/c12 renders; merge-vs-shell analysis; atom gap discussion |
| 2026-07-12 | Multi-fabric options; first principles; co-design in phase |
| 2026-07-12 | Fidelity ladder L0–L3; large GPU + adaptive grids; n-bar deferred |
| 2026-07-12 | `PROPOSAL.md`, `FIRST_STEPS.md`, this findings log |

---

## 5. Bottom line

**Setup:** U(1) gauged Cosserat kernel + v74 nuclear campaign + multi-fabric
proposal; big-GPU and nested grids welcomed for L0/L1.

**Findings:** Nuclei yes, atoms no, on one fabric; multi-fabric (3 first) is
the path; engagement from symmetry/conservation/EFT/ledgers; 3D mandatory for
geometry; n-bar only after PDE validates forces; **first concrete step** is
3D ±Q with existing seeds (FIRST_STEPS.md).
