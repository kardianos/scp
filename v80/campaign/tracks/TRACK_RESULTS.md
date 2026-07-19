# Tracks campaign results — force (2) + orbit (3)

**Campaign:** v80 tracks (keep-SFA)  
**Instance:** Vast 45305492 · V100-32GB · **destroyed** after analysis  
**Rulebook:** `../SCORECARD_AND.md`  
**Updated:** 2026-07-19 (queue complete; gates scored)

---

## Layer 1 integrity

| Gate | Status | Evidence |
|------|--------|----------|
| L0 Gauss | **PASS** | `gauss_max` ~9e-14 all physics runs |
| L1a Diags | **PASS** | `results/*_diag.tsv` + runlogs for all 10 jobs |
| L1b Tracks | **PASS** | `mf_pair_track` / `sfa_qball_track` TSVs for force, elite, both orbits |

---

## Gate F — Force (opposite multi-fab C↔L) — **PASS**

**Method:** rest pair, centers ±D/2, T=150, snap 2.5.  
**Metric:** \(D(t)\) from `mf_pair_track`; quadratic fit middle 70%: \(D = D_0 + v_0 t + \tfrac12 a_{\mathrm{rel}} t^2\).

| \(D_{\mathrm{nom}}\) | \(D(0)\) | \(D(T)\) | \(\Delta D\) | \(a_{\mathrm{rel}}\) | sign |
|----------------------|----------|----------|--------------|----------------------|------|
| 8 | 8.000 | 1.764 | −6.24 | −6.21×10⁻⁴ | attract (capture) |
| 12 | 12.000 | 9.694 | −2.31 | −2.13×10⁻⁴ | attract |
| 16 | 16.000 | 14.896 | −1.10 | −9.90×10⁻⁵ | attract |
| 20 | 20.000 | 19.461 | −0.54 | −4.67×10⁻⁵ | attract |
| 24 | 24.000 | 23.802 | −0.20 | −1.59×10⁻⁵ | attract |

**Scaling:** \(|a_{\mathrm{rel}}| \sim D^{-3.2}\) (full set); \(\sim D^{-3.6}\) for \(D\in\{12..24\}\) — steeper than Coulomb \(D^{-2}\), still monotone attractive.

**Ledger:** \(Q_\phi\) held at 114.72; Gauss floor. D=8 is merge-scale (\(E_{\mathrm{em}}\) 0.55→0.086); exclude from soft-force law quotes.

**Gate F: PASS** — attract at all five \(D\); \(|a_{\mathrm{rel}}|\) monotone ↓ with \(D\).

---

## Gate R — Same-sign control (elite) — **SOFT**

| \(D_{\mathrm{nom}}\) | \(D(0)\) | \(D(T)\) | \(a_{\mathrm{rel}}\) | note |
|----------------------|----------|----------|----------------------|------|
| 12 | 11.94 | 7.61 | −1.23×10⁻³ | attract + merge (2→1 clusters) |
| 16 | 15.97 | 15.77 | −8.67×10⁻⁵ | weak attract |
| 20 | 19.99 | 20.28 | **+1.00×10⁻⁴** | **repel** |

At \(D=20\): same-sign repels while multi-fab opposite attracts (desired polarity).  
At \(D=12,16\): same-sign still attracts (bag / overlap).  
Control is single-fabric elite, not multi-fab same-\(q_L\) — architecture mismatch noted.

**Gate R: SOFT**

---

## Gate O — Low-\(v_t\) orbit — **SOFT**

| Job | \(v_t\) | \(D(0)\!\to\!D(T)\) | revs | \(E_{\mathrm{em}}\) | note |
|-----|---------|----------------------|------|-------------------|------|
| `mf_orbit_R16_vt0p03` | 0.03 | 16.00 → **13.46** | **0.13** (47.5°) | 0.706 → 0.666 (94%) | slow **inspiral** |
| `mf_orbit_R16_vt0p05` | 0.05 | 16.00 → **20.92** | **0.16** (58.8°) | 0.706 → 0.657 (93%) | **flyby** expansion |

Naive orbital period \(T_{\mathrm{orb}} = 2\pi R / v_t\):

| \(v_t\) | \(T_{\mathrm{orb}}\) | \(T_{\mathrm{run}}/T_{\mathrm{orb}}\) |
|---------|----------------------|----------------------------------------|
| 0.03 | ~3350 | 0.12 |
| 0.05 | ~2010 | 0.20 |

**Interpretation:**

- Soft \(v_t\) with \(T=400\) covers only ~12–20% of one naive orbit — multi-rev was **kinematically unreachable**.
- vt=0.03: continuous attract + inspiral arc; masses/Q held; **not** v79 \(E_{\mathrm{em}}\) death.
- vt=0.05: tangential excess → separation grows; mild Ql drop (−2.5%); still not multi-rev.
- Neither meets PASS (≥2 revs **or** clear capture with bounded radius). Both rule out “instant soft-kill.”

**Gate O: SOFT** — re-test with \(T \gtrsim 2\,T_{\mathrm{orb}}\) or \(v_t \gtrsim 0.5\) for multi-rev; or longer soft-vt capture run.

---

## Notes / caveats

1. **Elite init=output symlink** wrote trajectories into seed paths — recoverable; fix cfgs next time.
2. **Force SFAs** tracked remotely; local archive has D8 + D12 (~19 G under `/space/scp/v80/tracks/`). D16–24 optional re-download if needed.
3. **Orbit SFAs** (~15 G each) left on remote until teardown; tracks retained locally.
4. Instance **45305492 destroyed** after scoring.

---

## Final verdict (gates 2+3)

```text
L0 / L1a / L1b: PASS
F Force:        PASS
R Same-sign:    SOFT
O Orbit:        SOFT
H Long-T:       PASS (prior S4)
P Pair dwell:   SOFT (prior)
N Multi-L:      N/A — do not re-park Z6+L6

VERDICT = CONTINUE_GATE_ONLY
  product two-body force is real (F PASS)
  orbit needs longer T or larger vt for multi-rev claim
  do not ADVANCE to staged multi-L
```
