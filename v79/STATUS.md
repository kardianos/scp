# v79 Status — GPU long-T multi-fabric atom

**Started:** 2026-07-18  
**Version:** v79 (execution after v78 package freeze)  
**Updated:** 2026-07-19 — both long-T runs COMPLETE

## Instance

| Item | Value |
|------|--------|
| Vast contract | **45259142** |
| Label | `v79atom` |
| GPU | **Tesla V100-SXM2-32GB** |
| Disk | **100 GB** |
| SSH | `ssh -p 19142 root@ssh7.vast.ai` |
| $/hr | ~0.21 |
| Remote work | `/root/v79/work` |
| Binary | `/root/v79/scp_sim_cuda` |

## Runs

| ID | T | Status |
|----|---|--------|
| **z6_a_n6_T800** | 800 | **COMPLETE** — wall 6160 s, Gauss OK, Q_phi −15%, Q_flux −68%, E_em collapse |
| **z6_n6_nuc_T800** | 800 | **COMPLETE** — wall 6128 s, Gauss OK, Q_phi −15% (same as atom), Q_flux +17%, E_em stable ~6 |

Full tables: `RESULTS.md`.

## Downloads

| File | Local |
|------|--------|
| Atom diag/runlog | `v79/results/` ✓ |
| Core diag/runlog | `v79/results/` ✓ |
| `z6_a_n6_T800.sfa` | `/space/scp/v79/` ✓ (4.2 GB) |
| `z6_n6_nuc_T800.sfa` | `/space/scp/v79/` ✓ (2.1 GB) |

## Headline finding

**Q_phi / Q_core / r_core END values match atom vs L=0 control to printed precision.**  
L sector does **not** save the charge ledger; it **nulls E_em and Q_flux**. Not a cold atom. See RESULTS §C.

## Teardown

After core SFA rsync finishes and sizes verified:

```bash
# optional: vastai destroy instance 45259142
```

Do not teardown while rsync is open.
