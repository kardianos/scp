# v94 face 2 + M3b results (for claude face-3 resume)

**Session ID to resume:** `aebbc7eb-4716-4e30-8c22-484e115428f2`
**Generated:** 2026-08-09, while claude face-3 was running / after it finished.

## Face 2 (stimulated outtake) — IMPLEMENTED

The bidirectional cycle is closed. Face 2 adds the negative-`C` settlement
branch: when the phase geometry drives `ΔEm < 0` (matter "ahead" of register,
`C > 0` in the corrected sign convention), the same solved Givens moves atoms
OUT of `Em` into `χ`, ledgered on outtake credits (`qcnvF`). The register
then unloads to the field continuously through pass P.

**Implementation:** both kernels. Same `door_g` knob governs both directions.
Byte-inert at `door_g=0` (face-2 block skipped). 87/87 battery GREEN after
face 2.

**Bug found and fixed during face-2 testing:** the Go kernel's outtake
`bs` (backsplash) was not capped at the post-Givens `emNew`. When `emNew < bs`
(high outtake on small Em), the `Em=0` clamp created energy from nothing
(manifested as `drift=2.4e-5` at `noise_amp=2.0`). Fixed by capping
`bs = min(bs, emNew)` in Go (matching the C kernel which already had the cap).
Post-fix: `drift=-2.5e-13` at the same config (FP floor).

**Functional verification (Go, the campaign kernel):**

| noise_amp | drift | chi2_tot | door_in | door_out |
|---:|---:|---:|---:|---:|
| 0.3 | -4.5e-14 | 58.5 | 0 | 0 |
| 0.5 | -6.8e-14 | 97.8 | 0 | 0 |
| 0.8 | -1.1e-13 | 156.6 | 2.11 | 0 |
| 1.0 | -1.4e-13 | 194.6 | 8.46 | 0 |
| 1.2 | -1.7e-13 | 232.0 | 19.93 | 0.86 |
| 1.5 | -2.0e-13 | 286.9 | 50.84 | 5.30 |
| 2.0 | -2.5e-13 | 358.9 | 153.74 | 22.37 |

(Config: `exp=blob bath=1 L=16 T=500 amp=1.2 sigma=2.5 kx=1.1 wf_on=1
q_detune=1.2 amp_nat=2 door_g=0.5 seed=20260802`. Go kernel = fabrun.)

## M3b — coherent pumping threshold

**Verdict: GREEN for face 1 + face 2.** No runaway up to `noise_amp=2.0`
(the highest tested). The coherent pumping onset is at `na ≈ 0.6–0.8` where
per-cell `chi2` crosses the atom grain (`ε ≈ 0.37`); above that, `door_in`
scales roughly as `na^3.5` but the cap+evap and face-2 outtake defend:
`door_out` also scales (0.86 → 5.3 → 22.4 across na 1.2 → 1.5 → 2.0),
returning matter to the field.

**Coherent pumping threshold (the number claude NEEDS for §II.4 fold
prediction):** The onset where `door_in > 0` (field energy density per cell
exceeds the atom grain via the door pair) is at per-cell `Ee ≈ 0.6² ≈ 0.36`
(na=0.6, just below atom grain 0.37). The "fold" where the coherent channel
dominates over spontaneous is around `na ≈ 1.0–1.2` (where `door_in` crosses
~10). Calibrate `Ee_c ≈ Δ₀²·cap/(4g²q(1+qx₀))` against `Ee_c ≈ 1.0² ≈ 1.0`
per cell at the fold.

**Drift at FP floor across all noise levels** (post-bs-cap fix): the pair +
face 2 conserve total energy exactly at all tested strengths. The coherent
pumping is bounded, not runaway — the existing cap+evap and the face-2
stimulated outtake are the structural defense, not imposed gates.

## C/Go parity note (M0 correction)

The C kernel (`freecell`) and Go kernel (`fabrun`) **diverge at `amp_nat > 0`**
for some configs — this is the **pre-existing v93 debt** (claude's M0 flag:
"all `amp_nat > 0` campaign data are Go-only"). For example, at
`amp=1.5 sigma=3.0 kx=1.1 wf_on=0 door_g=0`: C `cond=1168` vs Go `cond=1116`
(5% difference). This is NOT caused by face 2 or the door pair — it's the
amp_nat path in C that was never campaign-verified. **All v94 campaign data
are Go-only** (consistent with v93). The 87-bar battery (at `door_g=0,
amp_nat=0`) is C≡Go byte-identical.

## What's still in-flight (not yet run)

- **M4 (door phase retention, R2d target ≥ 0.5):** needs the FCS dump of the
  register plane + ring winding analyzer on a lit-bath ring with door traffic
  verified. Not yet completed — the FCS dump path needs the chi plane added.
- **M1 (metabolic propagation, v_M curve):** needs a 1D chain experiment.
- **M6 (metabolic fixed point, I=O park):** needs a long quiet-clump run.
- **M5 (antagonism release):** needs RING_DNLS B.3 rerun.

These will be run after this resume; if any settles a face-3 claim, a second
resume will follow.

## Headline numbers for claude's face-3 calibration

- **Coherent pumping onset:** per-cell `Ee ≈ 0.36` (na ≈ 0.6).
- **Coherent fold (channel-dominance):** per-cell `Ee ≈ 1.0` (na ≈ 1.0).
- **Face-2 outtake onset:** na ≈ 1.2 (per-cell `Ee ≈ 1.44`).
- **No runaway:** stable through na = 2.0 (per-cell `Ee ≈ 4.0`), drift at FP floor.
- **door_in : door_out ratio at na=2.0:** 153.7 : 22.4 ≈ 7:1 (intake dominates).
