# Approach D work folder

Owner: Approach D agent only.  
Log: `v76/logs/D_reverse_numeric.log` (append-only).

## Round 1 artifacts

| Path | Role |
|------|------|
| `scoring_skeleton.md` | Inversion/scoring design (D1–D3) |
| `invert_demo.py` | Full reverse demo: kernel / local-optics / dualist / softE |
| `gen_results_data.py` | Offline pure-stdlib table writer |
| `results/round1_*` | R1 tables |

## Round 2 artifacts

| Path | Role |
|------|------|
| `congruence_score_r2.py` | Ingest B maps; non-iso dualist; Occam; NC checklist |
| `offline_compute_r2.py` | Multi-start offline writer |
| `congruence_verdict_r2.md` | Goal-(2) / C-checklist verdict |
| `results/round2_*` | R2 scores, winners, NC, summary |

## Round 3 artifacts (3D free-response)

| Path | Role |
|------|------|
| `congruence_score_r3.py` | 3D free Green vs dualist 3D Poisson; multipole gate |
| `congruence_verdict_r3.md` | goal2_minimal vs goal2_PC3D verdict |
| `results/round3_*` | R3 scores |

## Round 4 artifacts (final)

| Path | Role |
|------|------|
| `congruence_score_r4.py` | R3 rays + inertia triad (if B ships) |
| `congruence_verdict_r4.md` | Final D goal-(2) verdict |
| `results/round4_*` | Final scores + J5 |

## Re-run

```bash
python3 /home/d/code/scp/v76/work/D/congruence_score_r3.py
python3 /home/d/code/scp/v76/work/D/congruence_score_r4.py
```

No `scp_sim` / kernel changes. Monist-eligible toys only.
