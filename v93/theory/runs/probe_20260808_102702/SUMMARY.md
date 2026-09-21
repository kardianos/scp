# theory/probe SUMMARY

seed=20260802  n=5  dir=`/home/d/code/scp/v93/theory/runs/probe_20260808_102702`

## Shear / field organization (genesis + progressive field meters)

| arm | family | cond | qcnvF | Jx | \|J\|_f | R2f | Ee | note |
|---|---|---:|---:|---:|---:|---:|---:|---|
| C1_nowf_kx | shear | 488.1 | 231.6 | -0.1645 | 0.3104 | 0.034 | 338.5 | pulse wind with wf_on=0 (door not presence-gated) |
| C2_nowf_k0 | ctrl | 472.4 | 235.2 | -0.0015 | 0.1395 | 0.355 | 351.6 | standing pulse wf_on=0 energy control for C1 |
| C3_blob_field | shear | 498 | 311.4 | -0.0157 | 1.2935 | 0.022 | 2293.3 | dense blob + ambient noise (matter presence for workfn door) |
| C4_blob_quiet | ctrl | 0 | 13.83 | +0.0000 | 0.0000 | 0.000 | 0.0 | dense blob no field noise control for C3 |
| C5_uni_noise | shear | 220.1 | 65.89 | +0.0018 | 1.0646 | 0.016 | 1856.5 | unitary dense + noise field (organized matter, lit medium) |

## Progressive dense / Bohmian-method transport

| arm | fd | fl | Em_tag | Em_tot | Em_max | PR | F_tag | Jd | note |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|

## Genesis scoreboard (cond>0?)

- **HIT** C1_nowf_kx cond=488.078
- **HIT** C2_nowf_k0 cond=472.402
- **HIT** C3_blob_field cond=497.963
- **HIT** C5_uni_noise cond=220.078

## How to read

- **Genesis**: cond, Em birth — primary shear-genesis bar.
- **Progressive field**: signed Jx high + R2f low ⇒ traveling organization; R2f high + Jx~0 ⇒ standing sterile.
- **Progressive dense**: fd (unitary hop current) is the guidance meter; ignore blob_drift.
