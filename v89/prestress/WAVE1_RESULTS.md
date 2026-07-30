# PRESTRESS Wave 1 results — frozen foam showdown (2026-07-29)

**Kernel:** integer default `cellfab` (`ledger_mode=3`).  
**Foam:** `foam_s20260727.tsv`, seed 20260727, L=24.  
**Laws:** `battery/laws_V2g.cfg`. Apparatus: edge_sink=1.6, lump_diag=1.  
**Fleet:** `prestress/fleet.go` (max 2 concurrent × 4 threads).  
**Scores:** `runs/SCORES.tsv` (score_net) + `runs/WAVE1_SCORED.tsv` (x50 / Law A).

---

## 1. What Wave 1 was for

Frozen foam (κ_plast=0). Pre-registered questions (E regression + RESUME):

| Question | Pre-reg |
|----------|---------|
| Does exact-ish retune save a cube? | Dies on load line t≈1875 (band 1250–2810) despite good gates; skin real only if ≥4700 or plateau |
| Is tube the structural exception? | ≥4600 if wound-mutual protection is real |
| Do chords without winding help? | 2000–3400 if closure redundancy is enough |
| Back-gate ladder (ring8_m3)? | Direction: wind down, load up → longer life |
| c8 = comp12 twin | Validation / spectra |

**Exact mass (M-R1) was not expected to pass** on frozen foam.

---

## 2. Scored table (measured x50)

Definitions (this campaign):

- **x50** = Emfree(t≈50) / (nv · cap), cap=2.5  
- **Law A** = 274 · (x50/0.0617)^1.066  
- **t_death** = first t with largest lump mass &lt; 0.3·m0 (score_net)  
- **c_eff** = cap · (x50 − x_skirt) / t_death  (whole-life per-voice form)  
- **c0 ref** = 4.25e−4  

| run | nv | x50 | Law A | t_death | t_d / LawA | c_eff | c_eff/c0 | gmin | verdict |
|-----|----|----:|------:|--------:|-----------:|------:|---------:|-----:|---------|
| **c2_cube150** | 8 | 0.385 | 1928 | **908** | **0.47** | 8.9e−4 | 2.1 | 0.401 | **EARLY** vs load line; not skin |
| **c2_cube150_ctrl** | 8 | 0.381 | 1907 | **698** | 0.37 | 1.1e−3 | 2.7 | 0.469 | control dies sooner |
| c8_ring12 | 12 | 0.413 | 2080 | **449** | 0.22 | 2.0e−3 | 4.6 | 1.000 | gates perfect; short life |
| c8_spectra (T=200) | 12 | 0.413 | 2080 | &gt;200 | — | — | — | 1.000 | spectra only; not a death run |
| **ring8_m3** | 8 | 0.547 | 2802 | **1631** | 0.58 | 7.4e−4 | 1.7 | 0.000 | **longest** W1; still dies |
| c5_tube6 | 12 | 0.427 | 2153 | **1042** | 0.48 | 8.8e−4 | 2.1 | 0.055 | **not** ≥4600 exception |
| c4_chords | 20 | 0.522 | 2671 | **630** | 0.24 | 1.8e−3 | 4.3 | 0.219 | no chord rescue |

Integer ledger: `max_sum_err=0` on scored full runs (conservation exact in int).

---

## 3. Decisions vs pre-registration

### Cube showdown (sharpest null)

| Claim | Result |
|-------|--------|
| On load line (1250–2810) despite gates | **No** — death **908**, ratio LawA **0.47** (below ×/1.5 band) |
| Skin if ≥4700 or plateau | **No** |
| Gates buy nothing vs ctrl | **Partial:** retune **1.30×** ctrl (908/698) — real but small; not a rate-level fix |

**Verdict:** On this foam + integer kernel + top-lump death definition, the retuned cube is an **early death**, not a parked consonant skin. Phase retuning still beats a mis-phased twin modestly. **P15 plasticity remains the only path past the frozen-foam gate ceiling** (formfind mean ~0.87, not ≥0.95).

### Tube (structural exception bet)

| Claim | Result |
|-------|--------|
| t_death ≥ 4600 | **Fail** (1042) |
| c_eff ≤ ~0.4·c0 | **Fail** (c_eff ~ 2.1·c0) |

**Verdict:** Wound tube is **not** the comp12-class exception on seed 20260727 under this apparatus. Exception remains **unconfirmed** outside historical comp12.

### Chords (closure vs chirality)

t_death 630 ≪ 2000–3400 band → **chirality / winding still required** for the long-lived class; redundant cycles alone did not protect.

### Back-gate ladder (ring8_m3)

Longest life in the wave (**1631**) at high load (x50≈0.55) despite poor mean gate (0.27) and gmin=0. Consistent with **load/back-gate direction** (E+C), not with “open gates = life.” Still not a particle.

### c8_ring12

Gates 1.0000 but t_death **449** — same qualitative lesson as regression: **seeded gate quality does not set lifetime** in this corpus. Parasite gpar_max=0.48 remains a watch item. Short T=200 spectra run is available under `w1_c8_spectra` for chiral-pump analysis (not scored here as death).

---

## 4. Cross-cutting conclusions

1. **No C1 plateau** in Wave 1. Exact mass (M-R1) correctly blocked.  
2. **Vacuum bleed / load-line story** still dominates directionally (all die; rough_share ~0.885 on rings), but **absolute Law A calibration is off** for this scorer/kernel: every full run dies **early** (ratios 0.22–0.58). Possible causes (do not pick one without a follow-up):
   - death definition (top-lump 0.3·m0 vs corpus rolling M_sum),
   - integer ledger + poly trig vs historical FP64 corpus,
   - x50 from Emfree vs historical x50 definition,
   - edge_sink / L=24 apparatus mismatch to corpus boxes.  
3. **Retune &gt; ctrl** is real (~1.3×) and should not be overstated as “gates don’t matter at all.”  
4. **Plasticity (Wave 2)** is the next mechanism test, not more frozen topology hunting.

---

## 5. Artifacts

| path | content |
|------|---------|
| `runs/SCORES.tsv` | score_net rows |
| `runs/WAVE1_SCORED.tsv` | x50, Law A, c_eff, verdicts |
| `runs/w1_*.log` | full diags + `# LEDGER` |
| `runs/w1_*.{slab,3d}.png` | last-frame renders |
| `fleet.go` | orchestrator used to finish the wave |

---

## 6. Next

- **Wave 2:** PLAST-1 cube ±κ* (gates ≥0.9 by t≈200, leak vs control; pin hunt).  
- Optional: re-score death with corpus M_sum definition to reconcile Law A ratios.  
- Optional: T≥5000 only if a Wave-2 object is still alive near band.
