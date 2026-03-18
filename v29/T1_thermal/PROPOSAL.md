# T1: Thermal Bath Equilibrium

## Question
Does the bimodal braid reach dynamic equilibrium in a radiation bath,
absorbing as much as it radiates?

## Hypothesis
In a universe filled with these fields, radiation is not loss — it's the
interaction mechanism. The braid should reach steady state where P_in = P_out.
The "absorbing boundary" in V28 was a crude model of this thermal bath.

## Method
1. Initialize bimodal braid (t=0.85 params) at N=96, L=20
2. Add low-amplitude thermal noise to ALL grid points:
   δφ_a(x) = A_noise × randn() at each point
   δv_a(x) = A_noise × ω_avg × randn()
3. Use PERIODIC BC in ALL directions (no absorbing layer)
4. Scan A_noise ∈ {0.001, 0.005, 0.01, 0.05, 0.1}
5. Run to T=1000
6. Measure at intervals of 50:
   - E_core(t), E_total(t), fc(t), |P|(t), winding(t)
   - FFT of φ₀ at 4 far-field probe points → spectral energy ε(ω)
   - dE_core/dt (radiation rate)

## Key Observables
- If dE_core/dt → 0 at late times: dynamic equilibrium achieved
- If ε(ω) ∝ ω³/(exp(ω/T)-1): Planck spectrum → well-defined temperature
- If winding survives in thermal noise: robust topological protection

## Success Criteria
- Braid survives T=1000 in noise bath with fc > 0.3
- Energy stabilizes (|dE/dt| < 0.01 × |E|/T at T > 500)

## Grid & Runtime
- N=96, L=20, periodic all: dx=0.42, dt=0.084
- T=1000: ~12000 steps, ~3 min per eval with OMP
- 5 noise levels + 1 control (no noise): 6 runs, ~18 min total

## Build
```
cd v29/T1_thermal && gcc -O3 -fopenmp -o t1 src/t1.c -lm
./t1
```
