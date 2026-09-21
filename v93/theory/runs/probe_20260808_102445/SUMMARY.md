# theory/probe SUMMARY

seed=20260802  n=22  dir=`/home/d/code/scp/v93/theory/runs/probe_20260808_102445`

## Shear / field organization (genesis + progressive field meters)

| arm | family | cond | qcnvF | Jx | \|J\|_f | R2f | Ee | note |
|---|---|---:|---:|---:|---:|---:|---:|---|
| A10_long_kx | shear | 0 | 28.39 | -0.4173 | 0.7030 | 0.005 | 1248.8 | longer wind — more door opportunity |
| A11_vac_kx | shear | 0 | 0 | +0.0000 | 0.0000 | 0.000 | 0.0 | vacuum pulse wind (no bath) |
| A12_slit_beam | shear | 0 | 0.5474 | -0.2683 | 0.7965 | 0.016 | 62.5 | directed slit field beam |
| A13_spin_kx | shear | 0 | 0 | -0.2852 | 0.4924 | 0.082 | 548.3 | spin + wind combined progressive field |
| A1_kx06 | shear | 0 | 0 | -0.1570 | 0.2578 | 0.117 | 548.3 | mild phase wind |
| A2_kx12 | shear | 0 | 0.1989 | -0.2775 | 0.4538 | 0.002 | 548.3 | phase wind kx=1.2 |
| A3_kx20 | shear | 0 | 7.067 | -0.1608 | 0.4297 | 0.010 | 548.3 | steep phase wind |
| A4_spin2 | shear | 0 | 0 | +0.0026 | 0.2191 | 0.002 | 548.3 | azimuthal field winding spin_m=2 |
| A5_imp_kx | shear | 0 | 71.5 | -0.5381 | 0.9273 | 0.075 | 548.3 | radial imp_k + phase wind |
| A6_hiamp_kx | shear | 0 | 229.7 | -1.0370 | 1.7252 | 0.001 | 1644.8 | high amp + wind — door threshold hunt |
| A7_hiamp_k0 | ctrl | 0 | 67.93 | -0.0003 | 0.3425 | 0.725 | 1644.8 | high amp standing — energy control for A6 |
| A8_noise_plus_kx | shear | 0 | 114.3 | -0.1216 | 1.1329 | 0.006 | 1905.4 | noise + organized pulse wind (lit + shear) |
| A9_noise_only | ctrl | 0 | 21.94 | +0.0048 | 0.7962 | 0.017 | 1372.4 | noise alone vs A8 |
| G0_dark | ctrl | 0 | 0 | +0.0000 | 0.0000 | 0.000 | 0.0 | dark bath — door floor |
| S0_pulse_k0 | ctrl | 0 | 0 | +0.0004 | 0.0802 | 0.763 | 548.3 | standing field pulse equal-Ee baseline (no shear) |
| S0b_noise | ctrl | 0 | 369.6 | -0.0257 | 1.5120 | 0.012 | 2744.9 | isotropic field noise (structureless energy) |

## Progressive dense / Bohmian-method transport

| arm | fd | fl | Em_tag | Em_tot | Em_max | PR | F_tag | Jd | note |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| B1_e3b_uni | +284 | +0 | 64.9 | 93 | 2.79 | 145.8 | 0.698 | 0.3027 | unitary linear dense hop — guidance transport |
| B2_e3b_add | +0 | -109 | 75.9 | 77.7 | 0.762 | 336.2 | 0.829 | 0.2262 | additive want control |
| B3_e3b_gated | +5.53 | +0 | 83 | 93 | 0.739 | 412.9 | 0.893 | 0.2506 | unitary gated τ |
| B4_uni_melt | +0 | +0 | 67.6 | 93 | 0.863 | 390.0 | 0.727 | 0.1614 | unitary + q_detune=0 melt |
| B5_uni_deep | +0 | +0 | 160 | 186 | 2.39 | 306.0 | 0.859 | 0.5524 | unitary + deep q_detune trap |
| B6_ring_mw | +0 | +0 | 0 | 4.88 | 0.813 | 6.0 | 1.000 | 1.4162 | vacuum ring matter winding |

## Genesis scoreboard (cond>0?)

- **No arm produced cond>0.** Shear genesis remains negative under this knob surface.
- Progressive field structure (Jx ≠ 0) can still be positive without birth.

## How to read

- **Genesis**: cond, Em birth — primary shear-genesis bar.
- **Progressive field**: signed Jx high + R2f low ⇒ traveling organization; R2f high + Jx~0 ⇒ standing sterile.
- **Progressive dense**: fd (unitary hop current) is the guidance meter; ignore blob_drift.
