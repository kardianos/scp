# v81 Path 3 FINDINGS — Token / Update-Budget CA

**Date:** 2026-07-19  
**Agent:** P3 sandbox spike  
**Binary:** `python3 src/run_all.py` (CPU, NumPy)  
**Wall time (full T0–T2):** ~20 s on workstation CPU  

---

## Scorecard

| Exp | Result | One-line |
|-----|--------|----------|
| **T0** Token conservation + hop cap | **PASS** | Drift \(=0\); max transfer \(=c=3\), never exceeded |
| **T1** Vortex–antivortex signed force | **PASS** | Clear **sign-dependent** interaction; **not** EM-like Coulomb |
| **T2** Long-lived bound / coherent pattern | **PASS** | \(\|C\|\) ratio \(0.70\), KE ratio \(0.43\), track frac \(1.0\) over \(T=4000\) |

**Aggregate:** all_passed = **true** (MVP gates).  
**Council two-week kill:** not fired on this session’s metrics — but see honesty on force law and atom readiness below.

---

## T0 — Conservation + hop cap

| Metric | Value |
|--------|-------|
| Grid | \(64^2\), \(\rho=20\), \(c=3\), \(T=500\) |
| Tokens initial / final | 135168 / 135168 |
| Max abs drift | **0** |
| Observed max bond transfer | **3** (\(=\;c\)) |
| Cap hit under overload | yes (east band `f_E += 10c`) |

**Interpretation:** The hop cap is a true engine law: overloaded directed occupations cannot push more than \(c\) tokens across a bond per tick; excess stays. Token total is an exact ledger (integer CA).

---

## T1 — Signed interaction of circulation defects

| Metric | Opposite \(\pm\) | Same \(++\) |
|--------|-----------------|-------------|
| \(d\) initial → final | 17.2 → 40.3 (\(\Delta d=+23.1\)) | 12.1 → 1.0 (\(\Delta d=-11.1\)) |
| Mean approach rate | \(-0.020\) (recedes) | \(+0.009\) (approaches) |
| Midpoint travel | 0.88 | 3.28 |
| \(\|\Delta\theta\|\) (pair axis) | 0.057 | 0.175 |

**Pass path used:** *viscous-gas signed force* — **same-sign merges / spirals in; opposite recedes.**  
Also: approach rates differ by circulation sign.

### Honesty — not Coulomb / not classical Euler dipole

| Expectation | Observed |
|-------------|----------|
| EM-like: opposite attract, same repel | **No** — opposite **recedes**, same **merges** |
| 2D Euler: opposite → translating dipole; same → co-rotate | **Weak** — little dipole mid-travel; same shows partial co-rotation while merging |
| Clear signed force | **Yes** — relative circulation sign controls separation history |

So T1 passes OP language **“opposite attract *(or clear signed force)*”** on the parenthetical branch. It does **not** yet deliver a positronium-class opposite attractor.

Likely mechanism: BGK + hop-cap lattice gas is a **viscous / under-resolved hydro** medium. Same-sign cores viscously merge; opposite pair spreads under residual pressure/noise rather than locking into a clean inviscid dipole. No pairwise force was inserted.

---

## T2 — Longevity (anti-evaporation)

| Metric | Value |
|--------|-------|
| Grid / \(T\) | \(96^2\), \(T=4000\), measure every 50 |
| \(\|C\|\) mean | 0.0111 → 0.0077 (**ratio 0.700**) |
| KE | 2448 → 1047 (**ratio 0.428**) |
| Centroid track fraction | **1.00** |
| Separation mean ± std | 44.3 ± 7.7 (pair slowly separates; stays trackable) |
| Token drift | **0** |

**Interpretation:** Coherent circulation patterns **do not evaporate** on this window the way multiplet field humps did in v80 soft orbits. Vorticity and kinetic energy relax toward a plateau rather than nulling. Tokens remain exact.

Caveat: “bound gyration” in the strong sense (closed multi-rev opposite orbit) is **not** claimed — the opposite pair slowly separates while remaining identifiable. Durability of *defects* is the win; *binding* is not.

---

## Kill-gate map (OP §1 + P3 two-week)

| Gate | Status |
|------|--------|
| State reduces to multiplet \(\phi\) | **Clear** — tokens + directed occupations only |
| Particle = field hump by definition | **Clear** — vortex = circulation defect |
| Inserted pairwise Coulomb | **Clear** — none |
| Hand-placed multiplet L | **Clear** |
| \(c\) as GRIN gravity | **Clear** — hard hop cap only |
| First metric = carbon spectroscopy | **Clear** |
| No clean opposite interaction | **Soft** — signed force yes; opposite *attract* **no** |
| Patterns evaporate under own radiation | **Clear on \(T=4000\)** — retained \(\|C\|\), KE plateau |

**Two-week kill (OP):** would fire only if *both* (no signed interaction) **and/or** (evaporation). This MVP avoids that kill. It does **not** graduate P3 toward kernel or Stage-3 atom delivery.

---

## Comparison intent (for orchestrator `COMPARE.md`)

| Axis | P3 token CA (this) |
|------|---------------------|
| Ledger honesty | Strong (exact integer tokens; hop cap audited) |
| Durability of “charge” carriers | Strong for circulation defects on short–medium \(T\) |
| F-like opposite attract | **Fail / not demonstrated** |
| Implementability | High (≤300 LOC Python, CPU seconds) |
| Kernel readiness | **None by design** — never into `scp_sim` this OP |
| Atom path | Research probe only (council rank 3) |

---

## What worked

1. Literal update-budget monism (shape 6): \(c\) binds; tokens exact.  
2. Charge-as-circulation is measurable and trackable with windowed centroids.  
3. Sign of circulation changes interaction (merge vs recede).  
4. Patterns outlive multiplet-style soft evaporation on the tested window.

## What did not work (yet)

1. EM-like opposite attraction / positronium bound state.  
2. Clean inviscid Euler dipole translation at usable SNR.  
3. Quantitative \(F(D)\) force law or Coulomb scaling.  
4. Any bridge to Cosserat nuclear stack or SFA export (out of scope).

## Next experiments (only if P3 is continued)

1. Lower viscosity / larger \(N\) / weaker hop-cap relative to \(\rho\) to hunt Euler dipole branch.  
2. Measure \(a_{\mathrm{rel}}(D)\) for both signs; publish sign-resolved force curves.  
3. Sequestration: park tokens in `rest` at cores; test inertia vs free budget.  
4. Hard kill still applies at two weeks if opposite attract never appears **and** the program needs that for Stage 3 — P1 remains the atom candidate.

---

## Reproducibility

```bash
cd v81/path3_token
python3 src/run_all.py
# expects results/summary.json all_passed=true
```

Key result files: `results/T0/t0_result.json`, `results/T1/t1_result.json`, `results/T2/t2_result.json`.
