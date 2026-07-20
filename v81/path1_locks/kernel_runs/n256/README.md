# N=256 lock orbit campaign (lock-Higgs = anti-lock)

**Not Cosserat multiplet Higgs.** Anti-lock is a lock species (`type=1`, `q=0`, no EM ρ/J).

| Run | Locks | Bag | Goal |
|-----|-------|-----|------|
| **base** | ± charge only | off | Large-box soft+\(v_t\) baseline |
| **anti** | ± + pinned anti-lock @0 | `r=12 k=0.8` | Bag-assisted park |
| **antib** | ± + anti-lock | stronger bag | Stress bag + larger radius |

Grid: `N=256 L=50 T=80`, `g_gauge=1`, vacuum Φ, periodic.

```bash
# launched sequential by master.sh / nohup
OMP_NUM_THREADS=16 bin/scp_sim base.cfg
```

Config: `lock2=0,m,x,y,z,...,pinned=1,type=1`  
Keys: `lock_bag_r`, `lock_bag_k`.
