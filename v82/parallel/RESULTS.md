# Parallel approaches A/B/C — results

**Date:** 2026-07-20  
**Binary:** `v82/parallel/bin/three_approaches`  
**Log:** `results/run_all.log`

## Scorecard

| Approach | Thesis | Result | One-line |
|----------|--------|--------|----------|
| **A** capacity well | Repel small \(D\), EM attract large \(D\) | **Well PASS / orbit FAIL** | Clean \(r_*\!\sim\!7{-}9\) for \(k=0.1{-}0.2\); free seed expands (sepmax~49) |
| **B** magnetic \(B_z\) | Gyration sets scale | **FAIL** | No sep band; high revs only with wild sep (not park) |
| **C** composite | One object, internal scale | **PASS\*** | Rigid \(\pm\) dipole: durable, \(E_{\mathrm{em}}\) holds, spin persists |

## A — capacity pair-overlap well

Pinned \(F_{\mathrm{net}}(D)=F_{\mathrm{EM}}-k e^{-D^2/2s^2}\):

- \(k=0.1,s=4\): **well** \(r_*\!\sim\!7\) (repel→attract)  
- \(k=0.2,s=4\): **well** \(r_*\!\sim\!9\)  
- Orbit at \(r_*\) with \(v_t\approx\sqrt{F_{\mathrm{EM}} r/\mu}\): **sepmin=9, sepmax~49, revs=0.22** — leaves the well (radiation / 2D soft dynamics / electrostatic re-solve lag)

**Keep A as primary scale mechanism** (force shape is right). Orbit needs better integrator (full PIC step + capacity) or 3D + softer radiation.

## B — magnetic

Uniform \(B_z\in\{0,0.05,0.15,0.3\}\) × \(v_t\): no band with \(\mathrm{sepmax}<2.5 D_0\) and \(\mathrm{revs}\ge0.5\) simultaneously.  
**Magnetic alone is not the second scale** (matches prior theory note). Keep as optional support later.

## C — composite

Rigid internal \(D_{\mathrm{int}}=6\), free COM + orientation:

- Cannot pass through by construction  
- \(E_{\mathrm{em}}\) final \(\sim0.34\) (live medium)  
- \(\omega\) remains finite  

**PASS\*** for “durable opposite-charge content without two free radial DOF.”  
Not yet positronium *dynamics* (relative \(r\) free); it is a **structural** Stage-3 candidate.

## Recommendation

| Priority | Action |
|----------|--------|
| **1** | **Merge A+C:** composite *or* free pair with **A’s well** + full charge-conserving PIC (not cheap Poisson re-solve each step) |
| **2** | Keep **C** as baseline “atom of bookkeeping” if free-pair multi-rev stays hard |
| **3** | Drop **B** as primary scale; optional \(B_z\) after well exists |
| Harmonics | Secondary: once \(r_*\) held, scan \(v_t\) for low-radiation survivors |

## Reproduce

```bash
make -C v82/parallel run
```
