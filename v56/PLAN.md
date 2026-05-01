# V56 — Multivector field in ℝ(3,1,1) on the Voronoi foam

## Premise

V54 concluded that the 6-field Cosserat equations don't support stable
particles across 750+ parameter configurations. Two structural reasons
the V55 foam kernel doesn't fix:

1. **No spinor structure.** The Cosserat θ field is a pseudovector
   under SO(3); a 2π rotation maps it to itself. Fermions need a
   double cover that picks up a sign — and our current fields don't.
2. **Trivial vacuum manifold.** A_bg cos(k·z) is a single field
   configuration; there's no topological charge protecting any
   particle-like state.

V56 adopts Lengyel's projective spacetime geometric algebra ℝ(3,1,1)
as the field's algebraic home. This brings (a) Lorentz covariance
built into the algebra, (b) spinor structure via the relativistic
quaternion subalgebra (8 real components, double cover of the proper
Lorentz group), and (c) a Higgs-like vacuum manifold whose homotopy
admits topological solitons.

See `LENGYEL_TAKEAWAYS.md` for what we get from the GA framework, and
`DERIVATION.md` for the concrete equations of motion + discrete
operators.

## Architecture

V56 starts from a copy of V55's foam kernel and changes the **field
representation and dynamics**, not the simulation infrastructure:

| Layer | Status |
|-------|--------|
| Voronoi mesh (gen_foam_mesh.cpp + voro++) | ✓ unchanged |
| Foam kernel scaffold (foam_sim.c) | ✓ unchanged structure, modified dynamics |
| SFA cell-native I/O (FMSH/FCEL/FCEP) | ✓ unchanged |
| volview cell-native renderer | needs minor update for n_columns=8 |
| Field count per cell | 6 → 8 |
| Mass / potential terms | Cosserat V(P) → Higgs (\|M\|²−v²)² |
| Curl coupling | dropped (no φ-θ exchange in stage A) |
| Initial conditions | new seed types for soliton candidates |

## Roadmap

### Stage A: Klein-Gordon multivector with Higgs potential

Free + quartic Lagrangian:
```
ℒ = ½ ∂_μ M^a ∂^μ M^a g_{ab} − ¼ λ (|M|²_bulk − v²)²
```
Equation of motion (per component):
```
□ M[i] = − λ (|M|²_bulk − v²) g_{ii} M[i]
```

Tasks:
1. Modify `foam_sim.c` to allocate 8 fields per cell, use the metric
   tensor `g = diag(−,−,−,−,+,+,+,+)`, replace `compute_forces` to use
   the Higgs term.
2. Update `gen_foam_mesh.cpp` — no change needed (mesh is generic).
3. Update `cellsfa_to_voxel.c` — bump `n_columns` to 8.
4. Update `volview/main.go` — handle 8-column FCEL/FCEP with sensible
   RGB mapping (probably show |M|² and the rotor amplitude).
5. Bench with `bench/run_bench.fish` — confirm dynamics still fit in
   the same ms/step envelope.

Experiments:
- A1: Vacuum verification (M = (v,0,...,0), should be stationary)
- A2: Tachyonic instability check (M=0 → fall to vacuum)
- A3: Q-ball candidate (localised perturbation in rotor part)
- A4: Skyrme hedgehog seed (winding number 1)

### Stage B: Dirac multivector (if Stage A produces solitons)

Replace the second-order Klein-Gordon with first-order Dirac:
```
γ^μ ∂_μ ψ + m ψ = ∂U/∂ψ̄
```
where ψ is the 8-component relativistic quaternion (a spinor in PGA
form). This requires:
- A new time integrator (split-operator or Crank-Nicolson)
- Encoded gamma matrices as algebra elements
- Reformulated potential `U(ψ̄ψ)` instead of `U(|M|²)`

This is where actual fermions could appear — Pauli exclusion, half-
integer spin, electron-like solitons.

### Stage C: Topological refinement

If A or B produces soliton candidates, add the Skyrme term:
```
ℒ_skyrme = (1/32 e²) Tr([L_μ, L_ν]²),  L_μ = q⁻¹ ∂_μ q
```
on the rotor part `q = M[0..3]`. This guarantees winding-number-stable
solitons.

## Files

```
v56/
  PLAN.md                — this file
  LENGYEL_TAKEAWAYS.md   — what the GA framework provides (and doesn't)
  DERIVATION.md          — equations of motion + discrete operators
  foam/                  — kernel + format + tools (copied from v55)
    foam_sim.c           — to be modified for 8-component dynamics
    gen_foam_mesh.cpp    — unchanged
    cellsfa_to_voxel.c   — bump n_columns
    voro_src/            — unchanged
    bench/               — same harness, fresh results.tsv
```

## Open theoretical questions

1. Which 8 components? Lengyel's relativistic quaternion picks
   `(eᵢⱼ₀, 𝟙, eᵢ, e₃₂₁)` — temporal-translation rotors plus spatial
   pseudoscalar. Alternatives: full 16-component even subalgebra
   (gives pure Lorentz boosts on equal footing); 10-component
   bivector (most physical — momentum + EM tensor structure).
2. What's the right metric on M? Lengyel's bulk norm signs depend on
   grade; we should verify the (+,+,+,+,−,−,−,−) split for our
   chosen basis.
3. Does the relativistic quaternion's double cover give real spin-1/2
   under our discretisation? Need to test with a known Dirac seed.
4. What replaces the η × curl coupling that gave EM-like behaviour
   in v55? It might emerge naturally as the bivector self-coupling
   term in the algebra; need to derive.
