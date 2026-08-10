# CALLING — v93 theory/sustain campaign (re-read every hour)

**Authority:** User 2026-08-08 — continue experiments until **11:00 AM
2026-08-09 local (CDT)**; modify experiments, develop theory, refine method,
judgment calls, **update kernel** as needed. Over-invest. Prefer **Go
`fabrun`**. Files primarily under `v93/theory/`; kernel in `v93/fab/` (+ C
mirror if physics changes).

**Deadline:** 2026-08-09 **11:00 America/Chicago**.  
**Stop new multi-hour arms after ~09:00 AM 08-09** so they can finish by 11:00.  
Last ~2 h = scoring, final theory writeup, no new launches.

**Stable paths (survive compaction):**

| path | role |
|---|---|
| `v93/theory/CALLING.md` | **this file — re-read first every hour** |
| `v93/theory/HOURLOG.md` | append-only hourly status + decisions |
| `v93/theory/MECHANISM.md` | isolation → mechanism |
| `v93/theory/SUSTAIN_CAMPAIGN.md` | T=3000 matrix design |
| `v93/theory/R1_ENDURANCE.md` | T=10000 endurance on HOLD winners |
| `v93/theory/NUMERICS.md` | scoreboard |
| `v93/theory/runs/ACTIVE_SUSTAIN` | symlink → primary T=3000 dir |
| `v93/theory/runs/ACTIVE_R1` | symlink → R1 endurance dir (when live) |
| `v93/theory/bin/probe` | Go harness |

---

## 1. Mission

Determine whether a condensate **sustains**, **why** (which switches), and
whether that mechanism can become **internal** dynamics. Secondary: shear
vs genesis (mostly closed: shear not the birth switch; workfn is the pure-field
wall). Bohmian-method = progressive amplitude guidance + link currents +
non-tag identity — not pilot particles.

**By 11:00 08-09 deliver:**

1. Full T=3000 matrix scored (`SUSTAIN_REPORT.md` final).  
2. R1 T=10000 on compact HOLD winners (validate endurance).  
3. Optional: monodispersity / rad fixed-point / progressive deeper probes.  
4. Updated MECHANISM + NUMERICS + HOURLOG trail.  
5. Clear next-internal-dynamics recommendation.

---

## 2. Standing facts (do not re-discover)

1. Pure field + `wf_on=1` → **cond≡0** (door presence gate).
2. `wf_on=0` → massive birth; **C1≈C2 / G01≈G02** (shear ≠ genesis switch).
3. Progressive **Jx** / unitary **fd** real; gated τ → STRUCT_MELT at long T.
4. Tags lie; use Em_tot / nHot / Em_max / Rg — not Em_tag alone.
5. Global R₂ inverted for traveling waves.
6. **T=3000 compact HOLD:** S00, S04 (deep), S10≡S00, S14 (kx, extended), S16 (seed_mw).  
   Requires **unitary linear + q_detune trap**. Additive DRAIN; qd=0 DIFFUSE.
7. Radiance → SLOW_DRAIN on quiet clumps (S05≡S18); deep resists somewhat.
8. Field **noise fuel** → box-fill plateau Em~2200–2370, nHot~1300+ — **not**
   the same object as quiet HOLD. Do not call box-fill a particle.
9. Genesis long: open-door pure field → **diffuse** Em~640 plateau (not compact).
10. Em-only automated HOLD is insufficient — always check structure.

---

## 3. Campaign queue (priority order)

### A. Finish primary (if not `.done`)
`theory/runs/sustain_T3000_20260808_104548` — 29 arms T=3000.

### B. R1 endurance (after A, or when HOLD list solid)
`probe -suite r1 -kernel fabrun -T 10000 -j 4`  
Arms: S00, S04, S14, S16 clones at T=10000 + optional S00+rad, S00 seed-robust.  
See `R1_ENDURANCE.md`.

### C. Refinement (overnight capacity)
- Monodispersity scan: k_rad × qd on quiet unitary (shorter T=3000 if needed).  
- S00 last-frame hoard spectrum vs Face C hierarchy.  
- Ring R01 tiny HOLD — optional longer vacuum ring.  
- Document what is **internal** (trap+U) vs **external fuel** (noise).

### D. Kernel (only if measured need)
Already: Go diag `Em_tot=`. Streaming logs in probe for new launches.  
Battery ALL GREEN at amp_nat=0 after any physics edit.

---

## 4. Hourly loop

1. **Read this file** fully.  
2. **Append** to `HOURLOG.md` (time, fabrun count, arms done, decisions).  
3. Check campaigns / score finished arms.  
4. **Act:** finish A → launch B → C as time allows.  
5. After **09:00 CDT 08-09:** no new multi-hour arms; finalize docs only.  
6. At **≥11:00 CDT 08-09:** stop; final judgment only if not already written.

---

## 5. Judgment priorities

1. Sustain yes/no + which switch > more genesis re-runs.  
2. Compact HOLD ≠ box-fill ≠ diffuse genesis plateau — keep classes separate.  
3. Honest negatives.  
4. Go fabrun preferred.  
5. Internalize: unitary+qd already law-compatible; amp_nat adoption needs R0 safety later.

---

*Updated 2026-08-08 ~17:00 CDT — deadline extended to 2026-08-09 11:00 CDT.*
