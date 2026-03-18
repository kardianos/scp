# T3: Late-Time Radiation Rate

## Question
How fast does the braid radiate in open space? Is it oscillon-like
(power-law, quasi-stable) or genuinely unstable (rapid decay)?

## Method
1. Large domain: L=60, N=192 (dx=0.63)
2. Damping ONLY at r > 54 (0.90L) — very far from core
3. Initialize bimodal braid at center
4. Run to T=2000
5. Track at intervals of 25:
   - E_core(t) (energy within R=8)
   - E_shell(t) (energy in R=8..30, the "radiation zone")
   - |P|(t), fc(t), winding(t)
   - dE_core/dt computed by finite difference
6. After T=500 (initialization shedding complete), fit:
   dE/dt ~ -C/t^α (power-law = oscillon)
   dE/dt ~ -C·exp(-γt) (exponential = resonance decay)
   dE/dt ~ const (steady radiation = blackbody-like)

## Key Observable
- α > 2: effectively stable (radiation decreases faster than 1/t²)
- α ≈ 1: marginal (logarithmic energy loss)
- α < 1 or constant: the braid is genuinely radiating and will die

## IMPORTANT: Dynamics Mass
The BIMODAL params have mass=1.50 (phys[14]=1.50). The V28 validated result
used this mass in the EOM (mass²=2.25). Use mass2_override=-1 (default from
params) to replicate V28. Do NOT use mass²=0 — that was a V27 test, not V28.

## Grid & Runtime
- N=192, L=60: 7.1M points, dx=0.63, dt=0.126
- T=2000: ~15900 steps
- Single-threaded: ~30 min. With 16-core OMP: ~3 min.
- Single run.

## Build
```
cd v29/T3_radiation && gcc -O3 -fopenmp -o t3 src/t3.c -lm
./t3
```
