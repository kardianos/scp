# v72–v73 Review: Actionable Findings, Retroactive Notes, Future Direction

**Date**: 2026-07-09  
**Scope**: Unstaged root docs + untracked `v72/`, `v73/`, `eta_qflow`, `gen_sfa_pair` (and related seed fixes).  
**Companion**: earlier accuracy pass in `tmp_root_docs_review.md` (largely superseded by later doc/tool fixes).

---

## 1. Executive judgment

| Layer | Judgment |
|---|---|
| Classical soliton numerics (v72 η-stationarity, v73 ring/yrast/symmetry) | **Valid and strong** |
| Process/fabric ledger language | **Useful reinterpretation**, not new dynamical laws |
| QFI “multipartite entanglement” claim | **Numbers real; quantum interpretation not established** |
| Classical ↔ quantum unification of real physics | **Odds slim to none**; still a useful *analog / forced-structure* program |

**Keep:** engineering honesty, negative results, symmetry audits.  
**Tighten:** quantum-language claims (QFI, ℏ_eff, “spin-1 particle”).  
**Prioritize next:** magneton, η+gauge window, skeptical QFI controls — not more process poetry.

---

## 2. Actionable findings (do these)

### 2.1 Documentation / bookkeeping

| # | Action | Why | Where |
|---|---|---|---|
| D1 | Prefer **PROCESS.md** as source of truth for v73 chronology; fix any DISCOVERIES line that still says fixJ is “designed not yet implemented” if it conflicts with § continued-3 | Self-consistency | `DISCOVERIES.md` V73 continued-2 item 18 vs continued-3 |
| D2 | In CONCEPT §9 “Measured status”, state explicitly: witness is **nQFI[ρ_Q] at hrange=1**; single-ball **field** observable is *not* null | Avoid overclaim; field nQFI ~ O(1) on singles | `CONCEPT.md` §9 |
| D3 | Keep gauged vs ungauged η-drift numbers **split** (do not re-merge under “gauged or ungauged”) | Already fixed once; easy to re-break | `CONCEPT.md` §3 |
| D4 | Bank large SFAs under `/space/scp/v72/`, `/space/scp/v73/` and archive; keep logs/TSVs in-repo | Disk policy; reproducibility | runner + rclone |
| D5 | Commit or consciously hold: seed tools (`eta_qflow`, `gen_sfa_pair`), Makefile SEED_TOOLS expansion, root doc bumps | Currently untracked/unstaged | git |

### 2.2 Code / tooling

| # | Action | Why | Where |
|---|---|---|---|
| T1 | Treat **`eta_qflow` stationary seed as mandatory** for any η>0 particle experiment | Θ=0 sheds 1.4–5.6%; false “drain floor” | `CLAUDE.md` already; enforce in run checklists |
| T2 | **SFA semantics** must stay `POSITION=0, ANGLE=1, VELOCITY=2` on every new writer | Scrambled-seed bug cost a full false-instability campaign | seed generators |
| T3 | Profile convention: **f(r) = per-component amplitude** (never f/√3) | 27× off in s | `radial_eta_soliton`, `eta_relaxer`, docs |
| T4 | Ring diagnostics: **do not trust winding@peak alone** on standing waves (two π nodes → fake winding 1); use `azi_profile` mode decomposition | Documented trap in PROCESS §4.3 | analysis + future ring papers |
| T5 | Add `gen_sfa_pair` to Makefile `SEED_TOOLS` if not already after commit | Pair QFI / force seeds | `sfa/Makefile` |
| T6 | Implement or queue **magneton tool**: μ from gauge columns of gring SFA vs μ\* ~ Q·J/(2M) | Highest-value open measurement | new analysis + `FUTURE.md` v73 list |
| T7 | **Gauged fixed-(Q,J) flow** (PROCESS §5.2 with `complex_gauge=1`) for exact spinning ring | Kernel gring_k used ungauged seed + projection | `eta_qflow` extension |
| T8 | Optional: twist-constrained relaxer only if (n,m) taxonomy is pursued | Unconstrained flow unwinds m | low priority unless chirality is a goal |

### 2.3 Physics experiments (priority order)

| # | Experiment | Success criterion | Falsifier |
|---|---|---|---|
| E1 | **Magneton / g-factor** of gauged n=1 ring | Stable μ; compare to Dirac-like μ\* | μ noise / no sustained B through hole |
| E2 | **η + gauge window** on spinning ring | Lifetime ≫ 10² t.u. at some (η,g) with L_z/Q flat and A_v/A_u≈1 | All η>0 still fission on e-fold ~30 t.u. scale |
| E3 | **Skeptical QFI suite** | Pair nQFI[ρ_Q]≫1 and single ≪1 **stable under** hrange scan, phase Δφ taxonomy, time-window, scrambled/surrogate controls | Bound crossing only at hrange=1 or only for one operator |
| E4 | Operator suite: ρ_Q, s, |Φ|², per-component ρ, Θ | Map which observables null on single soliton | Field-like ops always “entangled” → diagnostic is operator-artifact |
| E5 | n=2 gauged ring + ring–ball / ring–ring collisions | Integer spin transfer? | Continuous J washout only |
| E6 | Flavored stationary dressing (multi-ω `eta_qflow`) + GDR witness | Stationary flavored multiplet; clean mode spectrum | No multi-ω fixed point |
| E7 | Fusion tracking of f_Q across merger | Witness evolves continuously through fusion | Numerical garbage at contact |

### 2.4 Claim hygiene (how to speak about results)

| Claim style | Prefer | Avoid |
|---|---|---|
| ℏ_eff = Q | “U(1)-soliton identity; action per cycle ∝ charge (measured)” | “derives quantum ħ / quantizes the field” |
| QFI | “classical dynamical-structure diagnostic using soliton ℏ_eff; pair modes cross nQFI>1 at hrange=1 on ρ_Q” | “measured multipartite entanglement of the classical field” without caveats |
| Process form | “continuity / stationarity / VK in throughput language” | “new stability theorem beyond classical criteria” |
| Spin | “topological orbital J = nQ; free polarization spin on U(1)³; ring is K-isomer” | “SM spin-1 particle” without representation/statistics work |
| Unification | “analog program: which quantum-looking structures are forced by stability + symmetry” | “classical field that is quantum mechanics” |

---

## 3. Retroactive notes (what the last campaigns taught)

### 3.1 False physics from tooling (June-26 → v72)

Three independent defects produced the “intrinsic η-drain floor”:

1. **Wrong variational principle** — `eta_relaxer` fixed both ω and norm → charge leak → pressure cage theater.  
2. **√3 seed mis-normalization** — Φ_a = f/√3 vs per-component f → s off by 27×; baselines −2.93/−2.66/−2.27% were off-shell.  
3. **SFA semantic scramble** — `{1,3,2}` loaded Φ into Θ; “−8.1% and dies” was header corruption.

**Retroactive rule:** before declaring a radiation channel “intrinsic,” require (a) same stencil as the kernel, (b) fixed-*Q* (or full dynamics) stationarity, (c) byte-level load audit of one frame, (d) comparison to an η=0 *same-grid* floor.

### 3.2 What v72 actually established (keep)

- Exact stationary (Φ, Θ, ω) at η>0; kernel drift ≤ η=0 floor (ungauged); gauged matches floor.  
- VK branch; dressing binds ~0.3% M; **extends** Q_min slightly.  
- Axisymmetric û-prolate core + oblate Θ belt (η² deformation).  
- T=10³ endurance at η=0.5.  
- Single-ball ρ_Q QFI null on clean seeds; June single-ball η-scan was **transient-dominated**.  
- Pair ρ_Q nQFI ~ 3.5–3.9 at inter-ball wavevector (hrange=1).

### 3.3 What v73 actually established (keep)

- **P1–P3** as measured bookkeeping (2E_kin=ωQ; closed vs open flux; pass-through layment ~99.75%).  
- **n=1 Q-ring** stable at η=0; L_z/Q ≈ n; ring branch holds ~2.5× ball charge at same ω.  
- Thin ring does not bind; fat smoke-ring seed required.  
- **(1,1) twist** not a fixed-Q minimum; dynamically metastable at best.  
- **η fissions rings** with no threshold; energy matches two free balls at Q/2.  
- **Symmetry audit:** η and V(s) destroy exact continuous rotation; η=0 ring protection is equivariance, not barrier height alone.  
- **Gauged ring** holds spin (first charge + integer-ish J + EM loop energy).  
- **Polarization spin is free** (flat E(J) to saturation); ring is **K-isomer** (~+3% M); orbital ↔ spin interconversion under fixJ+η.  
- Kernel exactness cut: flat-band spin balls and rings±spin exact; Ω≠0 magnitude-traded rotors are transients.

### 3.4 Interpretive overreach to roll back mentally

| Overreach | Better reading |
|---|---|
| “Radiation-free soliton is 2D” as June-26 *physics* limit | Geometry of *dressing* is axisymmetric; stationarity is 3D fixed-Q |
| “η*≈0.25 feasibility wall” for QFI | Wall was Θ=0 / wrong seeds; dissolved for stationary seeds |
| “nQFI crossing = multipartite entanglement” | Protocol output under chosen ħ_eff, operator, hrange |
| “Particle *is* fabric throughput” as new law | Continuity + translating envelope; standard continuum |
| “Spin-1 vector particle” | Three complex components + free elliptical polarization (analogy) |

### 3.5 Prior accuracy review items — status

| Earlier issue | Status after later edits |
|---|---|
| maxF 5e-15 / 14k unbanked | **Fixed** — `qflow_w145_eta025_20k.log` |
| “(b),(c) fixed” while √3 still broken | **Fixed** in `eta_relaxer` + `radial_eta_soliton` |
| CONCEPT “gauged or ungauged” lumping −0.04% | **Fixed** — split wording |
| README/CLAUDE era end at v71 | **Fixed** — v72/v73 mentioned |
| Makefile missing U(1) seeds | **Mostly fixed** — gen_qball_* + eta_* listed; confirm `gen_sfa_pair` |

---

## 4. Future direction

### 4.1 Recommended program framing (1 sentence)

> A classical gauged multi-scalar continuum used as a laboratory for **which quantum-looking structures (stable particles, forces, spin carriers, band structure, correlation witnesses) are forced by stability, continuity, topology, and symmetry** — not a claim to replace quantum mechanics until hard gates fire.

### 4.2 Near-term track (highest scientific weight)

1. **Gauged spinning ring characterization**  
   - Magneton / g-factor (E1)  
   - Exact gauged (Q,J) state (T7)  
   - 2×2 VK matrix on (Q,J)  
   - n=2, collisions (E5)

2. **η–gauge coexistence**  
   - Map lifetime(η,g) for rings and balls (E2)  
   - Decide whether η is only a transient/spin-orbit channel or a viable dressed-matter coupling

3. **QFI as science, not slogan**  
   - hrange and operator suite (E3–E4)  
   - Δφ taxonomy with `gen_sfa_pair`  
   - Surrogate: phase-randomized / time-shuffled fields  
   - Publishable statement: what is *proven* classically vs what is *analogous* to Mazza-style QFI

### 4.3 Medium-term track (structure)

- Multi-frequency `eta_qflow` → stationary flavored multiplet + GDR on clean states  
- Fusion + f_Q tracking  
- Statistics probes (exchange of two identical rings/balls; braid phase if any)  
- Spectrum: discrete modes with selection rules beyond one GDR line

### 4.4 Hard gates before “unification” language

Do **not** escalate classical→quantum claims until several hold:

1. **Witness robustness** — multipartite-style bound stable under normalization and operator family, with documented classical meaning.  
2. **Universal ħ** — same ℏ_eff scale from independent processes (soliton E/ω, radiation, scattering).  
3. **Statistics** — exchange phases or exclusion-like signatures, not only classical J.  
4. **One dimensionless match** to a real constant without free retuning (still open since CONCEPT §8).  
5. **Controlled continuum limit** — results stable in N, L, sponge; not box artifacts (v70 lesson).

### 4.5 Explicit non-goals (for now)

- Rewriting gravity into this framework (CONCEPT §8 still open/null).  
- Treating process ontology as a substitute for the field equations.  
- Chasing (n,m) twist taxonomy before magneton + η+gauge.  
- Expanding QFI narrative without E3–E4.

### 4.6 Suggested document updates when committing

- `CONCEPT.md`: keep η-stationary particle + QFI status; add one sentence on **spin carriers** (polarization free J vs topological ring J = nQ) only after a short cohesive § draft — not a lab dump.  
- `FUTURE.md`: v73 open list already good; promote E1–E3 to top.  
- `DISCOVERIES.md`: chronological OK; ensure fixJ implementation note is not contradictory.  
- Version dirs: `v72/FINDINGS.md`, `v73/PROCESS.md` remain the detailed records.

---

## 5. One-page checklist (print this)

**Before any η>0 production run**

- [ ] Seed from `eta_qflow` (or proven stationary SFA), not Θ=0  
- [ ] Profile f = per-component  
- [ ] SFA semantics 0/1/2  
- [ ] Report drift vs η=0 same-grid floor  
- [ ] gauss_max floor if gauged  

**Before claiming QFI / entanglement**

- [ ] Single null + pair cross on **same** pipeline  
- [ ] hrange scan reported  
- [ ] ≥2 operators  
- [ ] Phase / separation control  
- [ ] Language: diagnostic vs Hilbert entanglement  

**Before claiming new spin physics**

- [ ] L_z/Q and mode decomp (not winding alone)  
- [ ] η=0 vs g>0 vs η>0 comparison  
- [ ] Yrast vs isomer energy if comparing to balls  
- [ ] Kernel drop of constrained state (exact vs transient)  

**Next three builds**

1. Magneton analysis on gring  
2. Gauged fixJ / exact gauged ring  
3. QFI skeptical suite  

---

## 6. Bottom line for the project

The **engine** (kernel, fixed-Q/J relaxers, gauges, diagnostics) is producing real, often beautiful classical field results: absolute stability with torsion dressing, topological spin carriers, free polarization spin, symmetry-limited fission, and clean correlation diagnostics.

The **metaphysics** (fabric process as QM, QFI as entanglement, ℏ_eff as quantum ħ) should stay one step **behind** the engine, written as analog hypotheses with falsifiers — which is already the project’s best habit when it is at its best.

**Direction to keep:** forced structure from stability + symmetry in a classical continuum.  
**Direction to avoid:** declaring unification when the measurements are still soliton phenomenology done carefully.
