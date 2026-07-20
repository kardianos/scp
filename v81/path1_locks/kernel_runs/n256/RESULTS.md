# N=256 campaign results — large box + anti-lock (lock-Higgs)

**Date:** 2026-07-19  
**Grid:** `N=256 L=50 T=100`, `g=1`, vacuum Φ, `locks_medium_only=1`, periodic  
**Wall:** ~53 min/run (CPU 16 threads, medium-only Maxwell)

## Design

| Species | `type` | EM ρ/J | Role |
|---------|--------|--------|------|
| Charge lock | 0 | yes | ±q=\pm1\) light carriers |
| **Anti-lock** | 1 | **no** | Lock-Higgs bag seed; attracts charge locks for \(r<r_{\mathrm{bag}}\) |

**Not Cosserat multiplet Higgs** (`higgs_v` unused). Bag force:  
\(F = -k_{\mathrm{bag}}(1/r-1/r_{\mathrm{bag}})_+\,\hat r\) toward each type-1 lock.

## Runs

| Run | Setup | Soft | Bag |
|-----|--------|------|-----|
| **base** | \(\pm\) only, \(D=16\), \(v_y=\pm0.12\) | r=2.5 k=0.25 | off |
| **anti** | \(\pm\) + pinned anti-lock @0 | same | r=12 k=0.8 |
| **antib** | \(\pm\) \(D=20\) \(v=0.15\) + anti | r=3 k=0.35 | r=16 k=1.2 |

## Results

| Run | `gauss_max` | alive | \(E_\star\) | sep \(0\to T\) | min–max | joining revs |
|-----|-------------|-------|------------|---------------|---------|--------------|
| **base** | \(2.7\times10^{-3}\) | 2 | 32 | 16.0 → **26.7** | 16–27 | 0.15 |
| **anti** | \(1.9\times10^{-3}\) | 3 | 72 | 16.0 → **23.4** | 16–23 | 0.13 |
| **antib** | \(1.9\times10^{-3}\) | 3 | 90 | 20.0 → **30.4** | 20–30 | 0.14 |

### Interpretation

1. **Large box works:** sep stays well inside \(L=50\); no wrap/image shred of the N=48 era.  
2. **Anti-lock slows expansion:** base \(\Delta\mathrm{sep}=+10.7\); anti \(+7.4\) at same initial \(D=16\). Bag is *real* (not decorative).  
3. **Still not multi-rev park:** ~0.15 joining-line revolutions; pair slowly unwinds outward (soft Coulomb + radiation + soft core).  
4. **Stronger bag/antib** does not buy multi-rev at these \(v_t\); larger radius still expands.  
5. **Gauss** \(O(10^{-3})\) under free motion with `locks_medium_only` (floor at \(t=0\); residual from medium-only path — better than old \(O(1)\), not Cosserat \(10^{-13}\)).  
6. **Locks durable** by construction; \(E_\star\) exact.

## Verdict

| Question | Answer |
|----------|--------|
| Does N=256 help? | **Yes** — room to measure orbit physics cleanly |
| Does anti-lock (lock-Higgs) help? | **Mildly** — reduces expansion; not enough for bound multi-rev |
| Cosserat Higgs? | **Not used** (correct) |
| Multi-rev positronium? | **Not yet** |

## Next (suggested)

1. Medium-**sourced** bag: anti-locks deposit a grid bag density; charge locks feel \(-\kappa\nabla B\) (true co-field, not pairwise \(1/r\) form).  
2. \(v_t\) scan at fixed bag for circular orbit (angular momentum match).  
3. Or sequestration capacity (\(n_{\mathrm{seq}}\)) as rest-energy soft well.

Reproduce:

```bash
make -C sfa sim/scp_sim && cp sfa/sim/scp_sim bin/
cd v81/path1_locks/kernel_runs/n256
OMP_NUM_THREADS=16 bin/scp_sim anti.cfg   # ~50 min
```
