# Parallel approaches A/B/C — results

**Date:** 2026-07-20  
**Binary:** `v82/parallel/bin/three_approaches`  
**Logs:** `results/run_all.log` (first pass), `results/A_fixed_run.log` (A with circular fix)

## Scorecard

| Approach | Thesis | Result | One-line |
|----------|--------|--------|----------|
| **A** capacity well | Repel small \(D\), EM attract large \(D\) + correct circular | **Well PASS / analytic orbit PASS / live orbit PASS** | \(r_0\!\approx\!8.6\); continuum circular family parks; live PIC multi-rev band at \(r_c\!\approx\!9.8\) |
| **B** magnetic \(B_z\) | Gyration sets scale | **FAIL** | No sep band; high revs only with wild sep (not park) |
| **C** composite | One object, internal scale | **PASS\*** | Rigid \(\pm\) dipole: durable, \(E_{\mathrm{em}}\) holds, spin persists |

## A — capacity well (fixed circular)

### Force shape

Continuum design: \(F_{\mathrm{EM}}=1/(2\pi D)\), \(F_{\mathrm{core}}=-k e^{-D^2/(2s^2)}\), \(k=0.2\), \(s=4\).

- Force zero \(r_0\approx 8.644\) (repel → attract). **Well PASS.**
- Lattice `measure_Fem` at **integer** \(D\): ratio to continuum \(\sim0.92\)–\(1.07\) (good). Off-grid dense sampling is **poisoned by CIC self-force** — do not use for circular design.

### Circular condition (two bugs fixed)

1. **Not** at \(F_{\mathrm{net}}=0\): need \(F_{\mathrm{along}}(r_c)=\mu v_{\mathrm{rel}}^2/r_c>0\) (attract side).
2. Equal-mass COM: seed each lock with \(v_{\mathrm{each}}=v_{\mathrm{rel}}/2\), not \(v_{\mathrm{rel}}\).

Old seed at \(r_0\) used bare \(F_{\mathrm{EM}}\) as if reduced → \(v_{\mathrm{each}}\sim0.20\). Correct first shell: \(r_c\sim9.04\), \(v_{\mathrm{each}}\sim0.034\) (ratio \(\sim0.17\)).

### Orbit results (`Afix_orbit.tsv`)

| Mode | Best | Result |
|------|------|--------|
| **analytic** central \(F_{\mathrm{net}}\) | all 6 seeds | **PASS** — sep flat to machine noise, revs \(2.4\)–\(4.5\), rms \(\approx0\) |
| **live** PIC re-Poisson + core | \(r_c=9.82\), \(v=0.063\) | **PASS** — sep\([6.33,13.19]\), revs \(=2.00\), rms \(=1.59\) (band) |
| **CONTROL** wrong seed analytic | \(r_0\), \(v=0.20\) | expands sepmax \(\sim68\) as expected |

Live is noisier than analytic (self-force, SOR lag, 2D log + periodic): several seeds fail band on sepmax/rms, but **at least one multi-rev band** promotes live.

### Interpretation

- Capacity well + correct circular family is a **real** bound shell under the force law (analytic).
- Live monist PIC can hold multi-rev with the same seed — not only a continuum fantasy.
- Remaining live fragility is integrator / self-force (GAP A3), not wrong equilibrium formula.

## B — magnetic

Uniform \(B_z\in\{0,0.05,0.15,0.3\}\) × \(v_t\): no band with \(\mathrm{sepmax}<2.5 D_0\) and \(\mathrm{revs}\ge0.5\) simultaneously.  
**Magnetic alone is not the second scale.** Demoted.

## C — composite

Rigid internal \(D_{\mathrm{int}}=6\), free COM + orientation:

- Cannot pass through by construction  
- \(E_{\mathrm{em}}\) final \(\sim0.34\) (live medium)  
- \(\omega\) remains finite  

**PASS\*** for “durable opposite-charge content without two free radial DOF.”  
Not yet positronium *dynamics* (relative \(r\) free); structural Stage-3 candidate.  
Open: COM self-force budget on long runs (see GAP C).

## Recommendation

| Priority | Action |
|----------|--------|
| **1** | **Keep A** as primary free-pair scale: port capacity into monist kernel only after live band is robust (more steps, better J/Gauss, 3D). |
| **2** | **Keep C** as durable bookkeeping baseline if free-pair live stays fragile. |
| **3** | **Drop B** as primary scale; optional \(B_z\) later. |
| Harmonics | Secondary once \(r_c\) held under full PIC. |

## Reproduce

```bash
make -C v82/parallel
./v82/parallel/bin/three_approaches A v82/parallel/results   # fixed A
./v82/parallel/bin/three_approaches all v82/parallel/results # A+B+C
```
