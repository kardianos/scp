# R1 endurance — T=10000 on structural HOLD winners

**Deadline context:** continue through **2026-08-09 11:00 CDT**.  
**Launch:** after primary `sustain_T3000_*.done` (or immediately if done).  
**Harness:** `theory/bin/probe -suite r1 -kernel fabrun -T 10000 -j 4`

## Arms (8)

| id | clone | question |
|---|---|---|
| R1_S00 | S00 uni+qd1.2 quiet | compact HOLD to 10k? |
| R1_S04 | deep qd=3.6 | deep peaks endure? |
| R1_S14 | kx progressive | HOLD_ext Rg bounded? |
| R1_S16 | seed_mw=2 | winding-compatible? |
| R1_S00_rad | S00 + k_rad=0.05 | Em* park / slow drain? |
| R1_S00_s111 | S00 seed=111 | seed-robust? |
| R1_S04_rad | deep + rad | deep resists tax at 10k? |
| R1_S00_amp12 | heavier seed | mass dependence? |

Shared: L=16 dt=0.02 bath=1 amp_nat=2 amp_logate=1 hop_order=1 wf_on=1 (except noted),  
snap_every=10000 (~200 t.u.), diag_every=2000.

## Pass bars

- Em last-half |slope| ≲ 0.02 and Em_end ≳ 100  
- nHot ≳ 20, Em_max ≳ 1.0  
- Rg_hot not runaway (S14: prefer ≲ 8)  
- R1_S00_rad: if Em parks near fixed point, note Em* / n_hot product  

## Wall time

~3.3× T3000 ≈ **~3 h/arm** → j=4 ≈ **6 h** for 8 arms. Fits overnight before 09:00 stop-new-arms on 08-09.

## Launch

```bash
cd /home/d/code/scp/v93
go build -o theory/bin/probe ./theory/probe
OUT=theory/runs/R1_T10000_$(date +%Y%m%d_%H%M%S)
ln -sfn "$(pwd)/$OUT" theory/runs/ACTIVE_R1
setsid env GOMAXPROCS=4 ./theory/bin/probe -suite r1 -kernel fabrun -T 10000 -j 4 \
  -seed 20260802 -snap-every 10000 -out "$OUT" \
  > "${OUT}.launch.log" 2>&1 &
```

## Explicit non-arms

No noise box-fills, no pure genesis G* as R1 primary.
