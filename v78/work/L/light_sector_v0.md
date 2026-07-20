# Light Sector v0 — L fabric prep (Phase L)

**Campaign:** v78  
**Owner:** agent L  
**Status:** **COMPLETE** — CP-L **ADOPT** (2026-07-19)  
**Sources:** `v75/STATE.md`, `v75/LIGHT_L_PLAN.md`, `v75/FIRST_STEPS.md`, `v75/FINDINGS.md` F11–F19, `v75/ATOM_LADDER.md`, `v75/C12_ATOM_GOAL.md`  
**Does not modify:** kernel, SFA format, or v75 sources (read-only synthesize)

---

## 1. L fabric definition (frozen)

| Item | Value |
|------|--------|
| Fabric | \(\Phi_L^a\), \(\Theta_L^a\) — third Cosserat sector (`n_fabrics=3`) |
| EM weight | \(\boxed{q_L=-1}\) (with \(q_C=0\), \(q_Q=+1\)) |
| Role | Light **opposite-charge** sector (electron analog) |
| Coupling | Shared abelian \(A_\mu\); matter transport \(U^{q_L}=\mathrm{e}^{i q_L\theta}\) |
| Bag | **Private** \(s_L=\prod_a|\Phi_L^a|^2\) — never mixes with \(s_C\) |
| Potential | Same \(V_t(s)\) package as nuclear (nuclear-class \(V_t\); not a separate leptonic potential yet) |
| Default profile | \(\omega_L=+1.46\), `f_w146_g005` → \(Q_L\approx 114\), \(E\approx 173\), mass\(_\rho\approx 78\) |
| Soft min (g=0.05) | \(\omega_L\approx 1.485\), `f_w1485_g005` → \(Q_L\approx 90\) (hierarchy ceiling for single ball) |
| Seed rule | **Same-sign Noether \(\omega\)** on C and L; fabric \(q\) supplies EM sign flip (opp \(\omega\) + opp \(q\) double-flips → null force) |

**EM source:**

\[
\rho_{\mathrm{em}} = q_C\rho_C + q_Q\rho_Q + q_L\rho_L
= 0\cdot\rho_C + (+1)\rho_Q + (-1)\rho_L.
\]

**Electron analog bookkeeping (F17 freeze):** \(\Phi_L\) with \(q_L\rho_L\) opposite nuclear \(Q\)-sector charge. Cloud size for atoms is **L count \(= Z\)** (proton count), **not** \(A=Z+N\).

---

## 2. Isolation (private bag)

**Claim:** L and nuclear sectors can contact in real space without bag merge.

| Gate | Evidence | Verdict |
|------|----------|---------|
| E-lite contact | F10 head-on **merge** (same fabric) | motivates multi-fabric |
| B1 G3 equal-mass head-on | F11: D→1.02 then recede to 60; mass C/L both ~78 intact; **pass-through** | **PASS** |
| Hierarchy head-on H1 | F12: D_min≈0.15; massC 221→225, massL 78→79 | **PASS** |
| Soft-branch head-on H2 | F13: two centers; L mass −22% after contact | **PARTIAL** (no full merge; shed) |
| Dual-L rest | F14 A2: massL 157→158; Ql stable | **PASS** |

**Freeze:** Isolation is **proven** for \(\omega_L=1.46\) under B1/B2 private bags. Soft-branch \(\omega=1.485\) survives as a separate center but is **contact-fragile**. Do not rely on softer single-ball ω alone for atom contact.

---

## 3. Hydrogenoid path (single heavy + light)

Path from force → orbit → packaging (pre multi-Z carbon core):

### 3.1 Same-fabric precursor (E-lite)

| Step | Result | ID |
|------|--------|-----|
| Opposite force D=20 | Attract ΔD class −0.3; 2 clusters T=100 | F8–F9 **PASS** |
| Orbit scan | sub inspiral / circ arc / super expand; head-on merge | F10 **PARTIAL** |
| Multi-rev immortality | period ~3300; T=400 only ~¼ period | **OPEN** |

E-lite = force/orbit calibration; **not** hydrogenoid at contact (merge).

### 3.2 Multi-fabric force + hierarchy

| Run | Setup | ΔD (D₀=20, T=100) | Verdict |
|-----|--------|-------------------|---------|
| F11 G2 | equal 1.46+1.46 | −0.336 | **PASS** (matches E-lite) |
| F12 H1 | C@1.42 + L@1.46 | −0.656 | **PASS** (hierarchy) |
| F13 H2 | C@1.42 + L@1.485 | −0.616 | **PASS** force |

Mass ratio M_C/M_L ≈ 2.6 (H1) … 3.3 (H2). First rung toward \(m_e\ll M_N\); **not** electron-scale (Q_min≈90 wall at g=0.05).

### 3.3 Soft orbit + packaging

| Rung | Setup | Result |
|------|--------|--------|
| A1 soft orbit | C@1.42 + L@1.46, D=20, vt hierarchy | **PASS** — inspiral, massL held, no merge (F14) |
| B4 packaging | n_C=1 ω=1.42 + n_L=6 R_shell=18 | **PASS** scorecard (F15) — massL/Ql flat T=400 |
| F16 pair kinematics | B4 θ, vt sub/vc/super @ R=18 | **PASS** ordered ΔD; L stable; true \(v_c\) **below** naive 0.071 (bracket ~0.05–0.07) |

**Hydrogenoid readiness:** force, isolation, soft orbit arc, and multi-L rest shell around **single** heavy center are **ready**. Multi-rev bound orbit (P1.1–P1.2) and shell-radius diagnostic (P1.3) remain **open** for Phase A / P1 stretch — not blockers for CP-L prep stamp.

---

## 4. Z-matched light cloud

**Rule:** number of L lumps \(n_L = Z\) (proton count); isotope \(N\) must **not** change L package.

| Scale | Setup | massL hold | L tracks Z not A | Verdict |
|-------|--------|------------|------------------|---------|
| Z=1 + L6 | B4 rest R=18 (F15) | 0% loss | single-C packaging | **PASS** |
| Z=2 + L=2 | B2, `gen_pn_core`, R~18 (F18) | −0.6% | identical for N=0 and N=2 | **PASS** |
| Z=6 + L=6 | B2, R=22, N=192 L=48 (F19) | −12.5% | identical for N=0 and N=6 | **PASS sector** / strict atom soft |

**Diagnostics freeze (STATE §5):**

- Atom: nuclear park + \(c_L\le 0.15\), massL_end > 0.5 massL_0.  
- **Do not** treat late net \(Q_\mathrm{flux}\) drop alone as L death when massL/Ql hold (net = nuclear EM − L in flux cube).

**Seed:** `gen_pn_core` … `nL` centers + optional `profL omegaL` (default light ω=1.46).  
**Grid safety (Z6):** L shell R=22 needs L=48 damp_width=5 (\(R_\mathrm{damp}=43\)); **unsafe** in F18 L=32 box (sponge eats L) — `PN_GRID.md`.

---

## 5. What PASS exists (L sector)

| ID | Claim | Status |
|----|-------|--------|
| **L-P1** | \(q_L=-1\) + \(U^q\) → opposite EM force vs nuclear Q | **PASS** F11 (null without \(U^q\)) |
| **L-P2** | Private bag: head-on no C–L merge | **PASS** F11–F12 |
| **L-P3** | Hierarchy light L (1.46) force + no-merge | **PASS** F12 |
| **L-P4** | Soft-branch L (1.485) force | **PASS** F13 |
| **L-P5** | Soft-branch L contact mass hold | **FAIL soft** F13 (−22%) |
| **L-P6** | Soft hydrogenoid orbit (pair) | **PASS** A1 / F16 kinematics |
| **L-P7** | Multi-L shell rest around single C | **PASS** A2, B4 F15 |
| **L-P8** | Z-matched cloud Z=2, L hold + isotope | **PASS** F18 |
| **L-P9** | Z-matched cloud Z=6, sector hold + isotope L fixed | **PARTIAL** F19 (−12.5% massL) |
| **L-P10** | Multi-rev circular hydrogenoid | **OPEN** |
| **L-P11** | Multi-L shell radius diagnostic | **OPEN** (COM D≈0 useless) |
| **L-P12** | True \(m_e\)-scale single ball | **BLOCKED** at g=0.05 (Q_min≈90) without package change |

**Phase L complete bar (GOALS §4):** *L-sector prep + positronium/hydrogenoid path* → **met**.

---

## 6. Open blockers for 6L cloud (carbon atom handoff)

These block **Phase A** full C₁₂ atom claim, not CP-L ADOPT:

| # | Blocker | Severity | Owner next |
|---|---------|----------|------------|
| B1 | **Z6+L6 massL −12.5%** over T=400 (F19); strict PASS_atom False | High for ideal atom | A: retune R_shell / soft θ / park quality |
| B2 | **Z6N0 nuclear park soft** (c_Q=0.184) couples to atom score | High | C/A: spacing/g; prefer Z6N6 core |
| B3 | Multi-ball nuclear seed → **single droplet** by t∼400 (not multi-center nucleons) | Medium (morphology) | C: document droplet as core; A uses parked droplet multipole |
| B4 | **Shell-radius / multi-cluster L tracker** missing | Medium (structure claim) | A: P1.3 diagnostic |
| B5 | Multi-rev / true circular vt retune (~0.05–0.06) not closed | Medium | A: P1.1–P1.2 |
| B6 | Single-ball light hierarchy ceiling Q≈90 | Low for prep; High for real m_e | Future package (multi-L already used); optional m_L / V_L later |
| B7 | Dual nuclear-class \(V_t\) on L (not leptonic potential) | Theory debt | Document; not required for Coulomb cloud |
| B8 | Long-T visual C₁₂ (T≫400, volview package) | Goal | A Phase 5 |

**Not a blocker for 6L:** isolation physics, q_L sign, Z-matching rule, B2 isotope control of L (F18/F19 identical L for +N).

---

## 7. Recipes (handoff to A)

### 7.1 Physics block (atom-ready)

```
n_fabrics=3
mf_lock_CQ=0          # B2 when P/N core
mf_stage=2
q_C=0  q_Q=1  q_L=-1
complex_phi=1 complex_gauge=1 g_gauge=0.05
m=1.5 m_theta=1.6 eta=0
# same-sign ω on C and L
```

### 7.2 Light profiles

| Use | Profile | ω |
|-----|---------|-----|
| Default L cloud | `v74/profiles/f_w146_g005.txt` | +1.46 |
| Soft hierarchy scout | `v74/profiles/f_w1485_g005.txt` | +1.485 (contact-fragile) |
| Heavy single-C hydrogenoid | `v74/profiles/f_w142_g005.txt` | +1.42 |

### 7.3 Seeds

| Tool | Role |
|------|------|
| `gen_mf_pair_boost` | C–L pair force/orbit (dual profiles) |
| `gen_mf_shell_orbit` | B4-class shell + optional vt |
| `gen_pn_core` | Z protons + N neutrons + **nL=Z** lights |

### 7.4 Proven packages

| Package | Config / data | Score |
|---------|---------------|-------|
| B4 single-C + L6 R=18 | F15 B4_full | PASS |
| Z2 + L2 | `v75/cfg/pn/p2_a_*.cfg` · F18 | PASS |
| Z6 + L6 R=22 | `v75/cfg/pn/z6_a_*.cfg` · F19 | sector partial |

### 7.5 Score tools

- Tracks: `mf_pair_track` (pair SoT)  
- Park: `v75/analysis/score_pn_park.py`  
- Atom: massL/Ql tracks + nuclear park; ignore net Q_flux death alone  

---

## 8. CP-L stamp

| Checkpoint | Stamp | Rationale |
|------------|-------|-----------|
| **CP-L** Light sector prep | **ADOPT** | L fabric \(q_L=-1\), isolation, hydrogenoid force/orbit path, and Z-matched cloud recipes are measured (F11–F19). Residual work is atom packaging quality (A), not missing L definition. |

**Co-agree:** U should co-stamp CP-L on board.  
**Unlocks:** Phase A cloud recipes (GOALS Phase 5).

---

## 9. Checkin / handoff

### Checkin (L → O/U)

- Phase L **COMPLETE**.  
- Deliverable: this file `work/L/light_sector_v0.md`.  
- Log: `logs/L_light_sector.log` entry L-001.  
- Stamp: **CP-L ADOPT**.

### Handoff → A (atom)

1. Use **Z-matched L**: \(n_L=Z\), ω_L=1.46, shell R≥18 (Z2) / R=22 (Z6), grid from `PN_GRID.md`.  
2. Prefer **Z6N6** parked core over Z6N0 for atom first claim.  
3. Score with park-aware + massL/Ql; do not fail solely on net Q_flux.  
4. Close P1 multi-rev + shell diagnostic when claiming bound structure.  
5. 6L carbon cloud: **sector-viable** at F19; **tighten massL hold** before ideal C₁₂ atom PASS.

### Handoff → U

- Relation R6 (Coulomb) + R8 (isolation) for L: evidence closed.  
- R10 for multi-L at Z=6: partial — board as soft for atom, not reject L prep.

### Explicit non-work (this phase)

- No new GPU runs in this stamp.  
- No kernel / sfa.h edits.  
- No claim of finished C₁₂ atom (that is A).
'''