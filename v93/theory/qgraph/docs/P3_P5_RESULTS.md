# P3 / P5 results (qgraph predictions vs fabrun)

**Run dir:** `theory/runs/qgraph_P5_20260808_175524`  
**Horizon:** T=500 (short; trend only) · L=16 · unitary amp_nat=2 · seed=20260802  
**Scored:** 2026-08-08 ~19:00 CDT

---

## P3 — does `p_rad` set the tax severity as Stefan law claims?

**Setup:** quiet blob + `k_rad=0.05`, compare `p_rad=2` vs `p_rad=4` (adopted V3a).

| arm | Em₀ | Em₅₀₀ | Δ% | dEm/dt (last-half slope) |
|---|---:|---:|---:|---:|
| **P3_prad2** | 148.8 | **99.2** | **−33%** | −0.077 |
| **P3_prad4** | 148.8 | **130.7** | **−12%** | −0.031 |

**Reading:** Lower power **taxes harder** at this amp (p=2 drains ~3× faster than p=4). Consistent with p high → only tall peaks radiate; p=2 hits the multi-hot clump more broadly. Adopted p_rad=4 is the *milder* quiet tax (matches S05 SLOW_DRAIN class, not catastrophic).

**vs UNIFIED_LAW P1/P3 claim** (p_rad = d_s+1): matter-support d_s~1.4 would prefer p_eff~2.4 — and indeed p=2 is harsher on HOLD matter. Full-foam d_s~2.9 supports p=4 as *field continuum* exponent, not matter-graph exponent. **P3 does not falsify dual I/O; it sharpens that p_rad is not set by compact-matter d_s.**

---

## P5 — door coherent vs random fuel preference?

**Claim (SIM_PREDICTIONS):** coherent intake / random ≥ 1.35 if door is selective.  
**Incumbent:** ratio = 1.00 (noise just fills).

**Setup:** blob + `noise_amp=1.5`, `amp_door=0` vs `amp_door=1`, quiet unitary otherwise.

| arm | Em₀ | Em₅₀₀ | nHot trend | class |
|---|---:|---:|---|---|
| **P5_door0** | 148.8 | **2327** | box-fill | capacity plateau |
| **P5_door1** | 148.8 | **2321** | box-fill | capacity plateau |

**Ratio Em(door1)/Em(door0) ≈ 1.00** at T=500 — **no coherent preference measured**. Both saturate the same box-fill Em~2320.  
**Verdict at this apparatus:** **FAIL** for door-selectivity claim (matches incumbent). Door knob does not convert noise fuel into a distinct particle class on short T. Longer T or phase-matched (non-noise) drive still open for Face B.

---

## Carry into sustain campaign

1. Radiance tax severity: p_rad=4 is mild SLOW_DRAIN; monodispersity needs **k_rad ladder** more than p_rad flip.  
2. Noise + any door → box-fill; do not call it metabolic particle.  
3. R1 `*_rad` arms (T=10000) still the long-horizon Em* park test.
