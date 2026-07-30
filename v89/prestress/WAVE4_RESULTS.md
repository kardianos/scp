# PRESTRESS Wave 4 results — P19-lite / P20 m-family / tube harden (2026-07-30)

**Kernel:** integer default `cellfab` (`ledger_mode=3`).  
**Foam:** **single** foam `foam_s20260727.tsv` (full multi-seed P19 deferred).  
**Apparatus:** `kappa_plast=1.0`, `tau_harden=50` on all arms.  
**Scores:** `runs/SCORES.tsv` + `runs/WAVE4_SCORED.tsv`.

---

## 1. What Wave 4 was for

| arm | pre-reg |
|-----|---------|
| P19-lite cube/ring8/c8 harden T=4000 | harden pins mass; tighter identity vs frozen (cluster prep, not full C3) |
| P20 m-family under harden | ≥2 coexisting m-classes with distinct late masses |
| tube harden T=4000 | exception re-test under pin |

**Not in this pass:** multi-foam formfind ensemble (true P19), P21 kicks.

---

## 2. Scored table

| run | role | nv | x50 | Law A | t_death | vs frozen/κ | hard_n end | verdict |
|-----|------|---:|----:|------:|--------:|------------:|-----------:|---------|
| tube_harden | exception re-test | 12 | 0.431 | 2179 | **1437** | W2 tube_k1 1613 (−11%) | 8931 | MID |
| p20_ring12_m6 | m-family | 12 | 0.403 | 2027 | **1401** | W3 1402 (≈same) | 8463 | ON_BAND |
| p19_c8_harden | P19-lite | 12 | 0.405 | 2037 | **1102** | W2 c8_k1 1269 (−13%) | 9277 | MID |
| p19_cube_harden | P19-lite | 8 | 0.389 | 1952 | **1017** | = W2 harden | 8245 | MID |
| p20_ring10_m4 | m-family | 10 | 0.579 | 2980 | 908 | W3 846 (+7%) | 7542 | EARLY |
| p20_ring9_m3 | m-family | 9 | 0.580 | 2985 | 571 | W3 761 (−25%) | 7878 | EARLY |
| **p20_ring12_m5** | m-family | 12 | 0.474 | 2408 | **26** | W3 **1940 → 26** | 8002 | **EARLY kill** |
| p19_ring8_harden | P19-lite | 8 | 0.545 | 2794 | **14** | W1 1631; same as κ kill | 8623 | EARLY |
| p20_ring8_m3 | m-family | 8 | 0.545 | 2794 | **14** | same | 7590 | EARLY |

Hardening **engages** (hard_n thousands by run end) on every arm, including those that die at t=14 (post-death foam still hardens vacuum-adjacent links).

---

## 3. Decisions vs pre-registration

### P19-lite (exact mass / pin)

| Claim | Result |
|-------|--------|
| Harden extends life / parks mass on this foam | **No** — cube 1017, c8 1102, ring8 14; all die; no C1 |
| Tighter mass identity than frozen | **Not measurable** — no object lives as a stable mass attractor; deaths span 14–1102 with no late plateau |
| Hardening is inert if object dies in anneal window | Ring8 / m5 die before pin can matter for the *lump* |

**Verdict:** On **single foam 20260727**, P19-lite does **not** produce exact mass. Multi-seed formfind ensemble remains required for a real M-R1 test; this pass only shows harden+κ is not a self-contained particle recipe for the pre-registered seeds.

### P20 m-family under harden

Frozen ladder lives (W3): m5 1940, m6 1402, m4 846, m9 761, m8 1631 (W1).  
Harden lives: m6 **1401**, m4 908, m9 571, m5 **26**, m8 **14**.

| Claim | Result |
|-------|--------|
| ≥2 coexisting stable m-classes | **Fail** — nothing stable |
| Distinct late masses beyond scatter | **Cannot claim** — most classes die early; survivors are still dying, not coexisting |
| Harden preserves m-ladder ranking | **No** — m5 (best frozen) is **destroyed**; m6 survives ≈ frozen life |

**Verdict:** m-family coexistence under harden is **not observed**. The harden+κ apparatus **reorders** the ladder (selectively kills m5/m8). Any future P20 needs either different κ/τ or a foam where frozen m-classes already live past ~5k.

### Tube exception re-test

t_death **1437** under harden < W2 κ-only **1613** and ≪ 4600. **Not** the exception class. Harden does not promote tube into comp12-class life.

---

## 4. Cross-cutting

1. **No C1, no M-R1, no P20 coexistence** on this single-foam pass.  
2. Harden is **implemented and active** but **not life-saving** on these seeds; sometimes life-**destroying** (ring12_m5, ring8).  
3. Closest ON_BAND objects in W4: **ring12_m6** (1401, 0.69×LawA).  
4. **Do not** fold κ_plast / tau_harden into `laws_V2g`.  
5. Honest next for mass spectrum: multi-foam formfind + only carry nets that already beat ~0.8×LawA frozen (free1, ring12_m5) — and **re-pin κ** before assuming harden.

---

## 5. Artifacts

| path | content |
|------|---------|
| `runs/WAVE4_SCORED.tsv` | full metrics |
| `runs/w4_*.log` | T=3000–4000 diags + `# PLAST` hard_n |
