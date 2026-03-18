# T5: Self-Consistent Metric Coupling

## Question
Does the bimodal braid produce spin-2 radiation when coupled to the
metric it creates?

## Motivation
The strain/torsion are currently "proxies." This test promotes them to
dynamical fields by making the fields propagate on the metric they create.
The bimodal braid is non-breathing and aspherical — the best candidate
for genuine spin-2 gravitational wave emission.

## Method
V25's self-consistent metric framework applied to V28 braid:

    g_ij = δ_ij + α_g × h_ij
    h_ij = ∂_i φ_j + ∂_j φ_i    (strain tensor, symmetrized)

Modified EOM:
    ∂²φ_a/∂t² = g^{ij} ∂_i ∂_j φ_a - ∂V/∂φ_a

where g^{ij} ≈ δ^{ij} - α_g h^{ij} (linearized inverse).

The force computation needs the full gradient tensor at each point.
This is a MODIFIED SOLVER — does not use braid_core.h force computation.

## IMPORTANT: Dynamics Mass
The BIMODAL params have mass=1.50 (m²=2.25 in EOM). V28's validated result
used this mass. The modified EOM should be:
    ∂²φ_a/∂t² = g^{ij} ∂_i ∂_j φ_a - m²φ_a - ∂V/∂φ_a
with m²=2.25 from BIMODAL[14]². Do NOT use mass²=0.

## Scan
- α_g ∈ {0, 0.0001, 0.001, 0.005, 0.01, 0.05}
- α_g = 0 is the control (standard V28)
- At each α_g: measure l=2 content of h_ij at far field
- Also track stability: too large α_g → tachyonic instability

## Key Observable
- l=2 fraction of h_ij at R=10: ratio of quadrupolar to monopolar
  metric perturbation in the far field
- If l=2/l=0 > 0.1 at any α_g: first genuine spin-2 signal

## Grid & Runtime
- N=96, L=20, T=300
- Modified force is ~2x slower (needs gradient tensor)
- ~6 min per eval, 6 values = ~36 min

## Build
```
cd v29/T5_metric && gcc -O3 -fopenmp -o t5 src/t5.c -lm
./t5
```
