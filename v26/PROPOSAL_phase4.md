# V26 Phase 4+: Torsion Constraint and Teleparallel Gravity

## Reference

Algebraic topology of knots/braids: https://aeb.win.tue.nl/at/algtop-5.html
Teleparallel gravity: frame field formulation where torsion = gravity.

## Thesis

The braided soliton (V26 mode 2: twisted tube, triple product only) creates
a STRONGLY CURVED local spacetime from its own strain field. This is not
Newtonian gravity (weak field, scalar) — it is naturally General Relativistic
(strong field, tensor) because:

1. The strain ε_{ij} is O(1) inside the braid (not a perturbation)
2. The torsion ω_{ij} is nonzero (the braid TWISTS the frame)
3. The metric g_{ij} = δ_{ij} + 2ε_{ij} is self-consistent by construction
4. Causality (ds = c·dt) is automatic (field propagation speed)

The torsion term κ_T serves double duty:
- Constrains the braid structure (EM-like, prevents over-twisting → tightens fc)
- Contributes to the gravitational connection (Einstein-Cartan formulation)

## The Lagrangian

    L = Σ_a [½(∂_t φ_a)² - ½(∂_i φ_a)(∂_i φ_a)]
        - (μ/2)P²/(1+κP²)
        - ½κ_T (∂_i φ_j - ∂_j φ_i)(∂_i φ_j - ∂_j φ_i)/4
        - ½κ_S (∂_i φ_j + ∂_j φ_i)(∂_i φ_j + ∂_j φ_i)/4

Decomposition of the gradient:
    ∂_i φ_j = ε_{ij} + ω_{ij}
    ε_{ij} = ½(∂_i φ_j + ∂_j φ_i)  → strain → gravity (spin-2)
    ω_{ij} = ½(∂_i φ_j - ∂_j φ_i)  → torsion → EM (spin-1) + braid constraint

The forces on φ_a:
    acc_a = ∇²φ_a                           [from kinetic gradient]
            + κ_S ∂_j(∂_j φ_a + ∂_a φ_j)/2  [symmetric strain force]
            + κ_T ∂_j(∂_j φ_a - ∂_a φ_j)/2  [antisymmetric torsion force]
            - ∂V/∂φ_a                        [triple product]

Simplifying: the symmetric part gives ½κ_S(∇²φ_a + ∂_a(div φ))
The antisymmetric part gives ½κ_T(∇²φ_a - ∂_a(div φ))

Combined: acc_a = (1 + ½κ_S + ½κ_T)∇²φ_a + ½(κ_S - κ_T)∂_a(div φ) - ∂V/∂φ_a

When κ_S = κ_T = η: this reduces to (1+η)∇²φ_a → just a speed rescaling.
When κ_S ≠ κ_T: compression and shear have DIFFERENT stiffnesses.

## What to Test (ALL 3D, N=128)

### Phase 4a: Torsion Constraint on the Braid

1. Take the V26 mode 2 configuration (twisted tube, triple-product only)
2. Add the torsion term κ_T with κ_T = {0.0, 0.1, 0.5, 1.0, 2.0}
3. Measure: does κ_T TIGHTEN the braid (increase fc)?
4. Measure the torsion field ω_{ij} at equilibrium
5. Compute the vorticity Ω_k = ε_{ijk}ω_{ij} — is it localized?

### Phase 4b: Separate Strain and Torsion

6. Scan κ_S and κ_T INDEPENDENTLY:
   (κ_S, κ_T) = {(0,0), (0.5,0), (0,0.5), (0.5,0.5), (1.0,0.5), (0.5,1.0)}
7. At each: measure fc, breathing (DFT), l=2 content
8. Does separating κ_S ≠ κ_T improve the l=2 fraction?
9. Does the torsion term prevent breathing while strain allows it?

### Phase 4c: Self-Consistent Teleparallel Metric

10. Compute g_{ij} = δ_{ij} + 2ε_{ij} at each grid point
11. Use g_{ij} in the Laplacian: ∂²φ_a/∂x_i² → g^{ij}∂_i∂_j φ_a
12. This is FULL nonlinear self-consistent metric, not weak-field
13. Does the braid survive in its own strongly curved metric?
14. Is the metric STABLE (converges) or unstable (runs away)?

### Phase 4d: Torsion as EM — Quantized Flux

15. Compute the total torsion flux through a surface bisecting the braid:
    Φ_T = ∫∫ Ω_k dS_k
16. Is it quantized (integer × 2π from the braid crossing number)?
17. How does it depend on the braid parameters (twist rate, tube radius)?

### Phase 4e: Two-Braid Interaction

18. Two braided solitons at separation D=30
19. They are NON-BREATHING → no monopole radiation contamination
20. The strain field mediates gravity (symmetric ε)
21. The torsion field mediates EM (antisymmetric ω)
22. Measure: is the force attractive? At what angular pattern (l=2)?
23. Is there ALSO a torsion-mediated force (EM-like, 1/r²)?

## Success Criteria

**Minimum**: κ_T tightens the braid (fc > 0.5, up from 0.37)
**Moderate**: torsion flux is quantized, l=2 > 10%
**Full**: self-consistent teleparallel metric stable, two-braid attraction
with quadrupolar angular pattern
**Breakthrough**: three forces (gravity from ε, EM from ω, strong from P)
all from one Lagrangian with correct spin (2, 1, 0)

## Parameters

μ=-20, κ=20 (NO mass term m=0)
κ_T scan: {0.0, 0.1, 0.5, 1.0, 2.0}
κ_S scan: {0.0, 0.5, 1.0}
Braid: twisted tube, L_twist = L, R_tube=3, A0=0.8
Grid: N=128, L=20, t=500

Compile: `gcc -O3 -fopenmp -Wall -o v26p4 src/v26_phase4.c -lm`
