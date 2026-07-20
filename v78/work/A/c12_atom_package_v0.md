# C₁₂ Atom Package v0 — End State, Pass Bar, Path, Readiness

**Agent:** A (atom)  
**Campaign:** v78 Phase 5  
**Stamp:** **CP-A ADOPT** (package ready for U co-agree)  
**As of:** 2026-07-19  
**Sources of truth:**  
- Goal: [`v75/C12_ATOM_GOAL.md`](../../v75/C12_ATOM_GOAL.md)  
- Freeze: [`v75/STATE.md`](../../v75/STATE.md)  
- Lab: [`v75/FINDINGS.md`](../../v75/FINDINGS.md) F11–F19, F14 atom ladder  
- Ladder: [`v75/ATOM_LADDER.md`](../../v75/ATOM_LADDER.md)  
- Nuclear map: [`v74/CARBON_MAP.md`](../../v74/CARBON_MAP.md)  

This document is the **atom readiness package**: what counts as done, what P/N+Z cloud means, pass bars, phased path P1–P3, current readiness, and what is **BLOCKED** (kernel / auth / GPU). It is not a claim that a time-stable C₁₂ atom has been simulated.

---

## 1. Ideal end state (product)

**Demonstrate a time-stable multi-particle C₁₂ analog** on the full gauged multi-fabric:

| Property | Requirement |
|----------|-------------|
| **Nuclear core** | Carbon-scale: **Z≈6** proton-analogs + optional **N** neutron-analogs; target **A=Z+N≈12** or documented Z-carbon droplet + path to A=12 |
| **P/N distinction** | Firm proton- vs neutron-analogs (see §2) — not two labels on one ball |
| **Light sector** | Multi-L opposite-charge cloud; **cloud count tracks Z**, not A |
| **Structure** | Bound atomic composite: core + cloud remain **distinct sectors** (private bags) |
| **Stability** | Long-T: **no bag merge**, **no unintended fission**, **no cloud/core dispersion** |
| **Fields** | Matter C/Q/L + shared \(A\); Gauss floor ~1e-13 |
| **Visual** | Headless volview package: nuclear charge + L sector + \|E\| / charge over time |

**Primary product name:** fabric **C₁₂ atom** = Z-carbon-class nuclear core + L cloud with \(n_L=Z\).  
**Stretch:** same Z, larger N (isotope) → metastability + decay calc/sim.

---

## 2. P/N + Z cloud (frozen bookkeeping)

From v75 **B2 unlock** + F17–F19 freeze (`STATE.md` §3):

| Symbol | Meaning | Seed content | EM |
|--------|---------|--------------|-----|
| **p-analog** | Proton | \(\Phi_C\) bag + \(\Phi_Q\) co-located | \(\rho_{\mathrm{em}}\approx q_Q\rho_Q\neq 0\) |
| **n-analog** | Neutron | \(\Phi_C\) bag only (\(\Phi_Q=0\) at site) | \(\rho_{\mathrm{em}}\approx 0\) |
| **e-analog / L** | Light opposite | \(\Phi_L\), \(q_L=-1\) | opposite nuclear Q |
| **Z** | Proton count | # p seeds (or \(Q_{\mathrm{em}}/Q_{\mathrm{em,1p}}\)) | fixes nuclear EM + L count |
| **N** | Neutron count | # n seeds | mass / binding; **no** net EM |
| **A** | Z+N | carbon-scale ≈12 ideal | history; droplet may forget multi-center |

**Fabric charges (locked):**

\[
q_C=0,\quad q_Q=+1,\quad q_L=-1,\quad g_{\mathrm{gauge}}=0.05
\]

**Isotope control:** vary \(N\) at fixed \(Z\) and fixed \(n_L=Z\).  
Must hold: \(Q_{\mathrm{flux}}\) (nuclear EM proxy) identical; massL / Ql evolution identical for L; bag \(Q_\phi\) and \(E\) grow with N.

**What does NOT define P/N:** flavored \(\Delta\omega\) alone (still same-sign EM); B1 lock alone (every nuclear lump sources EM).

**Seed tool:** `bin/gen_pn_core`  
```text
gen_pn_core N L profNuc omega out_C out_Q out_L \
  nZ  xz yz zz ... \
  nN  xn yn zn ... \
  nL  xl yl zl ... \
  [profL omegaL]
```
C gets all nuclear balls (Z+N); Q gets **only** Z proton centers; L shell optional with hierarchy `profL omegaL` (e.g. ω=1.46).

---

## 3. Pass bar (scorecard)

### 3.1 Nuclear park (no L / nuclear twin)

\[
c_{Q,\mathrm{park}}=\frac{|Q_{\mathrm{mid}}-Q_{\mathrm{end}}|}{|Q_{\mathrm{mid}}|},\quad
c_{Q\mathrm{em},\mathrm{park}}=\frac{\big||Qf_{\mathrm{mid}}|-|Qf_{\mathrm{end}}|\big|}{\max(|Qf_{\mathrm{mid}}|,1)}
\]

| Gate | Threshold |
|------|-----------|
| PASS_nuc | \(c_{Q,\mathrm{park}}\le 0.15\), \(c_{Q\mathrm{em},\mathrm{park}}\le 0.20\), Gauss floor, \(\|Q_{\mathrm{end}}\|>0.5\|Q_{\mathrm{mid}}\|\) |

### 3.2 Atom (nuclear + L)

| Gate | Threshold |
|------|-----------|
| Nuclear | PASS_nuc (or documented soft park with morphology hold) |
| L hold | \(c_L\le 0.15\), massL_end > 0.5 massL_0 |
| Isolation | L mass stays nonzero; **no bag merge** into C |
| Gauss | `gauss_max` ~ 1e-13 floor |
| Isotope | fixed Z: \(Q_{\mathrm{flux}}\) nuclear identical under +N; L track identical |

**Hard fails:** bag merge (L into C), unintended core fission, cloud dispersion (massL collapse), Gauss drift.

**Soft fails:** energy drift without morphology loss (document; retune later).

**Do not** treat late **net** \(Q_{\mathrm{flux}}\) drop alone as L death when massL/Ql hold (net = nuclear EM − L).

**Ideal long-T bar:** T ≳ several nuclear / orbital timescales (multi-rev if orbiting, or long rest with no secular L loss); clear two-sector morphology in renders; Z and N reportable.

**Tool:** `v75/analysis/score_pn_park.py`

### 3.3 Atom ladder (historical rungs — F14)

| Rung | Setup | Status |
|------|--------|--------|
| A1 soft orbit H1 | C@1.42 + L@1.46, D=20 | **PASS** (inspiral, no merge) |
| A2 dual L | C + 2×L rest | **PASS** (no L absorb) |
| A3 Z6 atom (B1) | 6C octa + 6L shell R≈15.6 | **FAIL** persistent (L −35%, nuc evaporates) |
| F16 B4 + park | Z6+L6 R=22 park-aware | **PASS_park** (pre-P/N B1) |
| F19 B2 Z6N6+L6 | P/N + isotope | isotope EM **PASS**; strict PASS_atom **False** |

---

## 4. Phased path P1–P3

Aligned with `C12_ATOM_GOAL.md`. Work in order; do not skip to chemistry geometry.

### Phase 1 — Bound hydrogenoid (single-center atom)

| ID | Task | Success | Status |
|----|------|---------|--------|
| P1.1 | Retune circular vt (~0.05–0.06) from F16 bracket | Near-flat D(t) vs sub/super | **OPEN** (F16 bracket exists) |
| P1.2 | Multi-rev orbit T ≳ 2000–4000 on B4 θ | No secular L loss; Gauss floor | **OPEN** |
| P1.3 | Multi-cluster / shell-radius diagnostic | Shell not only COM D≈0 | **OPEN** (needed) |
| P1.4 | Binding proxy (E vs free C+L) | Documented binding signal | **OPEN** |

**P1 done when:** multi-rev (or clear bound) single-C + L cloud, visualizable.  
**Evidence so far:** A1 soft orbit PASS short-T; B4 single-C + L6 rest PASS (F15).

### Phase 2 — Parked multi-Z nucleus + light cloud (+ P/N)

| ID | Task | Success | Status |
|----|------|---------|--------|
| P2.0 | Firm p vs n definition | Two species; \(Q_{\mathrm{em}}\propto Z\) | **DONE** (F17) |
| P2.1 | Parked templates (c6 / flavored p+n) | Qc park band; report Z,N | **DONE** Z2 (F18); **partial** Z6 (F19) |
| P2.2 | L shell from **Z** not A | PASS_park + massL; L tracks Z | **DONE** Z2; **sector PASS** Z6 (L −12.5%) |
| P2.3 | Soft kinematics around droplet | Ordered D_shell; no L absorb | **not run** |
| P2.4 | Long-T rest + volview package | No merge/fission/disperse | **partial** (T=400 only; droplet merge) |
| P2.5 | Isotope smoke: fixed Z, +ΔN | Core stability Δ; L charge fixed | **DONE** Z2 & Z6 EM |

**P2 done when:** Z-carbon-class core + L time-stable **and** Z/N independently reportable.  
**Today:** Z/N reportable; Z6N6 parks; Z6N0 soft; L hold soft at carbon class; multi-ball → **single droplet** by t∼400.

### Phase 3 — A≈12-class core + light cloud → ideal C₁₂

| ID | Task | Success | Status |
|----|------|---------|--------|
| P3.1 | Q_max / g survey for A=12 | Parked or quasi-stable A=12-class | **NOT STARTED** (v74: A=12 super-critical at g=0.05) |
| P3.2 | Born–Oppenheimer: freeze nuc, relax L | L in nuclear multipole field | **NOT STARTED** |
| P3.3 | Assemble (Z≈6, A≈12) + multi-L | PASS_park; no L death | **NOT STARTED** |
| P3.4 | Long-T + visual C₁₂ package | Ideal product | **NOT STARTED** |
| P3.5 | Stretch: neutron-rich / decay | Rate + channel vs stable control | **NOT STARTED** |

**v74 nuclear fact (g=0.05):** \(Q_{\max}=921\); \(12\times Q_N^{\min}\gtrsim 1080>Q_{\max}\) → free A=12 always super-critical. Z-carbon (A=6 light) parks. A=12 atom core may use Z6N6 droplet (A=12 inventory as bag+charge split) not 12 free nucleons on branch.

---

## 5. Current readiness (honest)

### 5.1 What is ready

| Capability | Evidence |
|------------|----------|
| Multi-fabric isolation (private bags) | F11–F14: no C–L bag merge |
| Coulomb / opposite charge | F11–F13; U^q; seed same-sign ω |
| Single-C + multi-L rest package | F15 B4 **PASS** |
| Soft orbit hierarchy (H1) | F14 A1 **PASS** |
| P/N firm (B2) | F17: n \(Q_{\mathrm{flux}}=0\); p charged |
| Isotope EM at Z=2 and Z=6 | F18/F19: identical \(Q_{\mathrm{flux}}\) under +N |
| Z6N6 nuclear park | F19 \(c_Q=0.046\) **PASS** |
| Sector L survival at Z=6 | F19 massL −12.5% (not sponge death) |
| Seed + score + grid freeze | `gen_pn_core`, `score_pn_park.py`, `PN_GRID.md` |
| Kernel surface for B2 | already in sim (auth used in v75); configs frozen |

### 5.2 What is not ready

| Gap | Detail |
|-----|--------|
| Strict PASS_atom at Z=6 | F19 PASS_atom False (N0 soft park + L −12.5%) |
| Multi-center nuclear crystal | Multi-ball seeds **fuse to droplet** by t∼400 |
| Long-T (T≫400) visual C₁₂ | Not claimed |
| Multi-rev hydrogenoid | P1.2 open |
| Shell-radius diagnostic | P1.3 open (COM D≈0 useless for concentric shell) |
| A=12 free-nucleon park at g=0.05 | Impossible on static branch (charge budget) |
| Full atom campaign re-run under v78 | Needs multi-fab binary + GPU time |

### 5.3 Readiness matrix (atom product)

| Layer | Status | Notes |
|-------|--------|--------|
| **Theory / freeze** | **READY** | STATE.md + this package |
| **Seed recipe** | **READY** | gen_pn_core; F19 geometries |
| **Config recipe** | **READY** | `v75/cfg/pn/z6_*.cfg` |
| **Scorecard** | **READY** | park-aware + sector L |
| **Z6N6 short-T smoke** | **PARTIAL PASS** | park yes; L soft; not ideal bar |
| **Ideal long-T C₁₂ atom** | **NOT READY** | P3 open |
| **Stretch isotope decay** | **NOT READY** | needs P/N + long-T first |

**One-liner readiness:**  
*Package and short-T Z6 isotope atom are science-ready; strict time-stable visual C₁₂ atom is **not** claimed — next is park/L retune then P3 long-T on multi-fab GPU.*

---

## 6. BLOCKED items (kernel / auth / GPU)

| Block | Type | What is needed | Who |
|-------|------|----------------|-----|
| **Multi-fab production re-run** | **GPU** | V100-16GB+ via scp-runner; N=192 multi-fab ~8 GB proven; ~48 min/run × matrix | O / ops; no kernel change |
| **Long-T (T≳2000–4000)** | **GPU** | Same binary; longer wall; storage under `/space/scp/` | ops |
| **N≳256 / larger shell** | **GPU / margin** | V100-32GB or careful staging; only if R_L must grow (`PN_GRID.md`) | ops |
| **New kernel physics** | **AUTH** | e.g. bag-only-on-C (Q no self-bag), ε_CQ portal, Option C triple Cosserat, sfa.h changes | **Human explicit** — do not edit `scp_sim` / `sfa.h` without OK |
| **Shell-radius diagnostic** | **tooling** | Analysis only (allowed); not a kernel block | A/L analysis |
| **g&lt;0.05 for true A=12 free park** | **profiles + GPU** | New radial profiles at lower g; not blocked on auth if using existing kernel knobs | N/C + GPU |
| **B2 already in tree** | — | v75 authorized B2 path; **re-use**, do not re-litigate | — |

**Non-blocks (available now without new auth):**  
local profile shoot, seed generation, score scripts, docs, replaying archived F19 diags/tracks from `/space/scp/v75/pn/z6/`.

---

## 7. Measured snapshot (F19 freeze numbers)

Grid: `N=192`, `L=48`, `damp_width=5`, p octa R=8, n R=5.5, L shell R=22, T=400.

| Metric | Z6N0 | Z6N6 |
|--------|------|------|
| Q_phi end | 1057 | 4921 (tracks A) |
| Q_flux end (nuclear) | **990.3** | **990.3** (identical) |
| c_Q_park | 0.184 soft | **0.046 PASS** |
| Atom massL 0→end | 421→369 (−12.5%) | **same** |
| Atom Q_flux end (net) | 670 | **670** identical |
| Gauss | floor | floor |
| PASS_atom strict | False | False |

Data: `/space/scp/v75/pn/z6/` · scores `v75/results/pn_z6/` · sheet `images/C6_atom_sheet.png`

---

## 8. Recommended next work (when multi-fab GPU available)

1. **Retune Z6N0 park** (spacing / soft θ / ω) → \(c_Q\le 0.15\).  
2. **Tighten L hold** at carbon shell (R, ω_L hierarchy, soft vt).  
3. **Shell-radius diagnostic** (P1.3).  
4. Prefer **Z6N6 + L6** as primary atom seed (better park than N0).  
5. **Long-T rest** T≳2000 + volview package → claim or fail ideal bar.  
6. Only then P3.1 g-survey / A=12 free vs Z6N6 droplet path.  
7. Kernel auth only if private-bag + retune still blocks atoms.

---

## 9. Cross-links (campaign)

| Doc | Role |
|-----|------|
| This file | Atom end state + readiness + blocks |
| [`recipe_checklist_multifab.md`](recipe_checklist_multifab.md) | Operator checklist when multi-fab ready |
| `v78/GOALS.md` §3.5 | Campaign atom target |
| `v78/CAMPAIGN_MAP.md` | CP-A co-agree with U |
| Agent log | `v78/logs/A_atom.log` |

---

## 10. CP-A stamp

**STAMP: CP-A ADOPT**

| Clause | Decision |
|--------|----------|
| End state defined | Yes |
| P/N + Z cloud frozen | Yes (v75 B2) |
| Pass bar machine-checkable | Yes |
| P1–P3 path ordered | Yes |
| Readiness honest | Partial product; package complete |
| Blocks explicit | GPU for runs; auth only for new kernel |

**U co-agree required** per CAMPAIGN_MAP before unlocking U final solely on CP-A.  
**Does not claim** ideal long-T C₁₂ atom PASS — claims **package ready** for that work.
