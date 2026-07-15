# Priority B — Lighter L on B1 (run design)

**Status:** design ready (2026-07-13)  
**Prerequisite:** B1 G2/G3 PASS (F11). Kernel unchanged except seeds/configs.  
**Priority:** hierarchy first; B2 unlock deferred.

---

## 1. Goal

F11 used **equal** balls (ω_C=ω_L=1.46, Q≈114, E≈173). That proved private-bag
isolation and multi-fabric Coulomb, but is **not** an atom mass hierarchy.

**Goal of this campaign:** attach a **lighter L** (higher ω, smaller Q, smaller E)
to a **heavier C/Q** under B1 lock, and show:

1. Coulomb still attracts (force sign).  
2. Head-on / soft contact still **no merge**.  
3. Mass ratio M_C/M_L ≳ 2.5 with L remaining a stable lump.

This is the first step toward m_e ≪ M_N. True electron-scale may need further
ω→ soft-branch min or multi-L later — **g=0.05 min single-ball Q ≈ 90**
(ω≈1.485), so hierarchy is **bounded** on this package.

---

## 2. Branch numbers (g=0.05, from `v69/theory/gscan.tsv`)

| ω | Q | E_total | r_half | Role |
|---|---|---------|--------|------|
| 1.420 | **311** | **456** | 3.80 | **Heavy C** (existing `f_w142_g005`) |
| 1.460 | 114 | 173 | 2.63 | **Light L baseline** (`f_w146_g005`) — F11 |
| 1.485 | **90** | **137** | 2.32 | Soft-branch min Q (profile **to shoot** if H2) |
| 1.406 | 921 | 1316 | — | Q_max / nuclear heavy (not this sprint) |

**H1 mass ratio** (C@1.42 / L@1.46): M_C/M_L ≈ **456/173 ≈ 2.63**, Q ratio ≈ **2.73**.  
**H2** (C@1.42 / L@1.485): ≈ **456/137 ≈ 3.3** — needs new gauged profile.

---

## 3. Physics package (frozen)

```
n_fabrics=3  mf_lock_CQ=1  mf_stage=1
q_C=0  q_Q=+1  q_L=-1
complex_phi=1  complex_gauge=1  g_gauge=0.05
m=1.5  m_theta=1.6  eta=0
# SAME-SIGN ω on C and L (fabric q supplies EM opposite)
```

Private bags: C self-bag, L self-bag, no shared V — G3 isolation must hold.

---

## 4. Campaign gates

### B-H1 — Hierarchy force (rest)  ★ first GPU run

| Item | Value |
|------|--------|
| Seed | C: `f_w142_g005` ω=+1.42 @ +D/2; L: `f_w146_g005` ω=+1.46 @ −D/2 |
| D | **20** (compare F11); optional D=24 if cores large |
| Grid | N=192, L=48 (same as F11) or L=56 if heavy r_half=3.8 needs pad |
| T | 100, snap_dt=5, diag_dt=0.5 |
| v | all 0 |
| IDs | `mf_h1_force_D20` |

**Naive a_rel** (order-of-magnitude):

\[
a_{\mathrm{rel}} \approx \frac{g^2 Q_C Q_L}{4\pi D^2}\left(\frac{1}{M_C}+\frac{1}{M_L}\right)
\sim 1.4\times 10^{-4}
\]

(F11 equal ≈ 7×10⁻⁵; expect **~2× larger** |ΔD| if calibration holds.)

| Pass | Fail |
|------|------|
| ΔD < 0 (attract) | repel / null (ΔD ≲ 0.01 like broken q) |
| 2 centers; massC, massL stable to ~5% | L dissolves or merges |
| gauss_max ~1e-13 | drift |
| \|Q_C\|, \|Q_L\| Noether stable | charge dump |

Analysis: `mf_pair_track` → D(t); compare ΔD to F11 −0.336.

### B-H1b — Hierarchy head-on (no-merge check)

| Item | Value |
|------|--------|
| Same seeds | + vr: C −0.05…−0.1, L +same (softer than F11 vr=0.1 optional) |
| T | 400, snap_dt=10 |
| ID | `mf_h1_headon_D20` |

| Pass | Fail |
|------|------|
| D_min contact then recede **or** soft capture **without** single merged bag | massC+massL fuse to one fabric blob; L swallowed into C bag |
| Both masses remain O(seed) | L evaporates |

**Note:** unequal mass → COM drifts; track **both** centroids and D. Soft capture of L into Coulomb orbit around C would be a **bonus** (hydrogenoid), not required for PASS.

### B-H2 — Soft-branch L (optional second round)

- Shoot `f_w1485_g005` (ω≈1.485, Q≈90) via `gauged_shooter` continuation.  
- Repeat H1 force + soft head-on with L@1.485, C@1.42.  
- Only if H1 PASS and hierarchy still too mild for atom story.

### B-H3 — Deferred

- Multi-L (Z electrons)  
- Nuclear heavy multi-ball C + light L cloud  
- B2 unlock  

---

## 5. Seed tool

`gen_mf_pair_boost` now supports **dual profiles**:

```bash
# H1 force (rest)
bin/gen_mf_pair_boost 192 48 20 \
  v74/profiles/f_w142_g005.txt 1.42 \
  v74/profiles/f_w146_g005.txt 1.46 \
  0 0 0  0 0 0 \
  /space/scp/v75/seeds/mf_h1_force_D20_C.sfa \
  /space/scp/v75/seeds/mf_h1_force_D20_L.sfa

# H1 head-on (vr=0.08 each → v_rel=0.16, slightly softer than F11 0.2)
bin/gen_mf_pair_boost 192 48 20 \
  v74/profiles/f_w142_g005.txt 1.42 \
  v74/profiles/f_w146_g005.txt 1.46 \
  -0.08 0 0  0.08 0 0 \
  /space/scp/v75/seeds/mf_h1_headon_D20_C.sfa \
  /space/scp/v75/seeds/mf_h1_headon_D20_L.sfa
```

Single-profile mode (F11 style) still works with 14 args after binary name.

---

## 6. Configs

See `v75/cfg/mf_h1_force_D20.cfg`, `v75/cfg/mf_h1_headon_D20.cfg`.

Same physics block as F11; only seeds/output names change. Consider
`L=56` if heavy ball feels boundary (damp_width=4, r_half_C=3.8, D=20 → wall
clearance ~28 — **L=48 OK** for first try).

---

## 7. Execution order (GPU)

```
1. Build gen_mf_pair_boost; make H1 seeds
2. Upload seeds + ensure q-fix binary on instance
3. mf_h1_force_D20   T=100  (~17 min @ 128 ms/step)
4. Track D(t); decide PASS
5. mf_h1_headon_D20  T=400  (~50 min) if force PASS
6. Track + remote D(t) plot; archive diags via rclone→scpsfa
7. FINDINGS F12
```

Disk: delete or finish archiving F11 SFAs before headon if free space tight
(`v75mf` had ~35–40G free after F11).

---

## 8. Risks

| Risk | Mitigation |
|------|------------|
| Heavy+light profile mismatch → radiation | Watch E drift; Gauss; mass loss |
| Coulomb too weak / strong at D=20 | Optional D=16 force; scale vs naive |
| volview useless for MF | Rely on `mf_pair_track` + D(t) PNG |
| Hierarchy still ≪ atom | H2 then multi-L; document Q_min≈90 wall |
| COM motion confuses merge call | Require **both** massC and massL intact |

---

## 9. Success → next

| Outcome | Next |
|---------|------|
| H1 force+headon PASS | Soft orbit / multi-L scout; or H2 lighter |
| Force FAIL | Debug mass ratio / D; not B2 yet |
| Head-on L swallowed | Isolation issue with unequal bags — probe L bag off |
| Hierarchy “enough” for carbon story | Z-carbon template + N light L |

---

## 10. Open questions (do not block H1)

1. Should C use ω=1.42 single ball or a **multi-ball nuclear template** (v74)?  
   → H1 uses single heavy ball; multi-C later.  
2. Is M ratio 2.6 enough to call “light”?  
   → Yes as **first rung**; label as hierarchy step not m_e.  
3. Dual bag potentials both nuclear-class V_t — L not “leptonic” potential.  
   → Still true; separate m_L / V_L is a later package change.
'''