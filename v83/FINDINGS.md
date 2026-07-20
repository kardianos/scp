# v83 FINDINGS — kernel force fix + E1/E2 campaign

**Date:** 2026-07-20  
**Plan:** `v83/PLAN.md` (Fable-informed)  
**Stop:** **E1b FAIL** (no multi-rev hybrid orbit); **E2 NULL** (no static multi-center)

---

## Kernel fix

**Bug:** Lock Lorentz force was `F = −(g q) E`. With Gauss `div E = g ρ` and outward E for +ρ, that **inverts** continuum Coulomb: hybrid Cosserat Q-ball (+Q) **repelled** opposite locks and **attracted** same-sign locks (E1a FAIL).

**Fix** (`sfa/sim/scp_locks.h`):
- `F_em = +(g q) E` (continuum: opposite attract, same-sign repel)
- Pinned locks now record full EM force in `f[]` (was soft/bag only)

**Dual regression after fix:**

| Test | Result |
|------|--------|
| Hybrid q=−1 at D=8 | ⟨F_x⟩=−7.2e−4 (inward) **PASS** |
| Hybrid q=+1 at D=8 | ⟨F_x⟩=+8.9e−4 (outward) **PASS** |
| Lock–lock opposite (medium_only) | mild **recede** (self-force/CIC dominates; not a clean continuum pair) |

Prior K-L “opposite attract” under `F=−(gq)E` was consistent with inverted force + self-force, not hybrid-external truth.

---

## E1a — hybrid force (re-run) → **PASS**

N=48, L=12, Q-ball Q≈315, g=0.05, free heavy opposite lock, Gauss ~7e−14.

| D | ⟨F_x⟩ (q=−1) | continuum −g²Q/(4πD²) | ratio |
|---|--------------|------------------------|-------|
| 5 | −2.12e−3 | −2.51e−3 | 0.85 |
| 6 | −1.47e−3 | −1.74e−3 | 0.85 |
| 8 | −7.22e−4 | −9.79e−4 | 0.74 |
| 10 | −2.86e−4 | −6.27e−4 | 0.46 |

**Gate:** attract sign + Coulomb-class magnitude + Gauss floor → **PASS**.  
H1 strong form still open (needs durable orbit), but **coupling is no longer broken**.

---

## E1b — free orbit hybrid → **FAIL**

Seed: r=8, m=2, v=√(|F|r/m)≈0.0537, soft core (2.5,0.15), T=1200 (~1.3 orb periods).

| Metric | Value |
|--------|-------|
| revs | **0.28** |
| r min / max | 3.03 / 12.5 (boundary L=12) |
| mean r / rms | 7.2 / 2.0 |
| Gauss | floor held |
| Q_ball | 315→252 (shedding / interaction) |
| E drift | −20% |

**Not** multi-rev band. Trajectory: mild inspiral → soft-core encounter → fling to sponge.  
Fable absorb branch: free hybrid orbit **not** a Stage-4 pilot yet — needs better second scale / box / mass ratio / radiation handling.

---

## E2 — gauged Cosserat phase map (minimal) → **NULL** (no park)

N=64, L=16, g=0.05, two ω=1.42 balls D=10, T=80. Cluster track:

| config | clusters t=0 → t=80 | morphology |
|--------|---------------------|------------|
| **co-phase** | 2 → **1** | merge (liquid-drop) |
| **all-anti** | 2 → **2** | separate (sep~10→18) |
| **mix 1-anti** | 2 → **1** | contact → merge |

No static multi-center standoff. Consistent with v71 interlock moral (frustration not load-bearing at equilibrium). **Not** a full D×flavor sweep; closes “easy park at D=10, g=0.05, T=80” as NULL.

---

## Scoreboard

| Step | Result | Implication |
|------|--------|-------------|
| Kernel F sign | **FIXED** | Hybrid Coulomb orientation correct |
| E1a force | **PASS** | Matter↔lock coupling works |
| E1b orbit | **FAIL** | No Stage-4 free-orbit pilot yet |
| E2 multi-center | **NULL** | No static 2-center park in smoke |

**Stopped** per plan (orbit FAIL + multi-center NULL). Did **not** open E3 anti-lock unit this session.

---

## Next (authorized follow-ons)

1. E1b hardening: larger box, heavier lock / lighter effective coupling, capacity surface scale, N=96–128, longer T.  
2. E2 denser scan only if multi-center remains priority (user 2B).  
3. E3 high-lock:1-anti with N_anti=0 control (Fable) after E1b path choice.

---

## Files

| Path | Role |
|------|------|
| `sfa/sim/scp_locks.h` | force sign + pinned f[] |
| `v83/e1/results/e1a_fixed_scan.tsv` | force table |
| `v83/e1/results/track_e1b.tsv` | orbit track |
| `v83/e2/results/clusters_*.tsv` | cluster counts |

---

## Parallel options push (E1b B1–B5, E2 C1–C3) — 2026-07-20

See `v83/e1b_e2_EVAL.md`.

- **E1b:** 0/31 band PASS. Best revs 0.58 (inspiral); flattest r ~0.34 rev (B1 r=10 vf=1). Soft/bag/mass ineffective.
- **E2:** No true multi-center park. Co merges (small D) or stalls (large D); anti separates; flavored interlock no static molecule.
- **Stop:** both tracks FAIL/NULL under option matrix (C4 skipped).
