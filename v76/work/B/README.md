# Approach B work folder

Owner: Approach B agent only.  
Log: `v76/logs/B_forward_numeric.log` (append-only).

## Round 4 (current) — inertia triad / J5

| Path | Role |
|------|------|
| `sandbox_r4_inertia.py` | U, boost model, m_ray, pair force, renorm |
| `outputs/round4_result.json` | full package + J5 gates |
| `outputs/round4_inertia_triad.tsv` | triad for D |
| `outputs/round4_summary.txt` | human summary |

**J5:** raw three-way FAIL (form factor); ray=ledger PASS; boost=field PASS; renorm PASS.

```bash
python3 /home/d/code/scp/v76/work/B/sandbox_r4_inertia.py --skip-grid
python3 /home/d/code/scp/v76/work/B/sandbox_r4_inertia.py --N 24 --iters 350
```

## Round 3 — 3D free-capacity F1

| Path | Role |
|------|------|
| `medium_design_r3_3d.md` | F1 3D design |
| `sandbox_r3_3d_free.py` | N³ SOR + Born rays + exports |
| `offline_round3_3d.py` | analytic 3D Green + mini SOR |
| `run_round3.sh` | convenience runner |
| `outputs/round3_result.json` | package + gates |
| `outputs/round3_rays.tsv` | D-ingest (sector_tag, phi_origin) |
| `outputs/round3_path_cost.tsv` | radial ψ ~ 1/r |
| `outputs/round3_summary.txt` | human summary |

**Verdict:** `monist_3d_1r_pass=True` — exterior ψ∝1/r, rays without second gravity solver.

```bash
python3 /home/d/code/scp/v76/work/B/offline_round3_3d.py
python3 /home/d/code/scp/v76/work/B/sandbox_r3_3d_free.py --N 32 --iters 450
```

## Round 2

2D free Laplace (log) + monist_kernel_failed for hand 1/r. See `outputs/round2_*`.

## Round 1

Local GRIN + postulated kernel. See `outputs/results*.json`.
