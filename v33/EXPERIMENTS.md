# V33 Characterization Campaign

## The Model (no modifications)
    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a
    V(P) = (μ/2)P²/(1+κP²), P=φ₀φ₁φ₂
    m²=2.25, μ=-41.3, κ=50, A_bg=0.1
    Single alloc, periodic BC, symplectic Verlet

## Confirmed: attraction D=20→18.2 (ΔD=-1.8), energy conserved (1.0000×)

## Experiments

### C1: Force Law (D=5 to D=80)
N=128, L=D+20, T=200
D ∈ {5, 8, 10, 12, 15, 18, 20, 25, 30, 40, 50, 60, 80}
Fit: F ∝ 1/D^n

### C2: Single Braid Steady State + Radial Profile
N=256, L=20, T=1000
Measure: ρ(r), energy flux through shells, intake/outtake
Save profiles every T=50

### C3: Two Braids Long Run
N=256, L=40, T=2000
Track D(t): orbit, merge, scatter?
Save field snapshots every T=200

### C4: Mass Dependence
N=128, L=30, T=200, D=20
m ∈ {0.0, 0.5, 1.0, 1.5, 2.0, 2.5}
Does force range change with mass?

### C5: Five Braids
N=256, L=80, T=1000
Pentagon at D=30
Clustering? Hierarchy? Bound pairs?
