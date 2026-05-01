# v56 derivation — multivector field in ℝ(3,1,1) on the foam mesh

## Goal

Concrete equations of motion + discrete operators, mappable directly
onto v55's `foam_sim.c`. The result is a working stage-A kernel: Klein–
Gordon dynamics for an 8-component multivector field with a Lorentz-
invariant self-coupling.

## Step 1: Field choice — the even subalgebra of ℝ(3,1,1)

Lengyel's "relativistic quaternion" picks 8 components (rotors that
combine spatial rotations + temporal translations). The full even
subalgebra has 16 components and includes both rotors and Lorentz
boosts as group elements. We have a choice:

| Field                | Components | Captures                     | Cost |
|----------------------|------------|------------------------------|------|
| Vector (grade 1)     | 5          | spacetime points             | low  |
| Bivector (grade 2)   | 10         | screws + EM tensors          | med  |
| Relativistic quat    | 8          | rotation + temporal trans    | med  |
| Even subalgebra      | 16         | all Lorentz transformations  | high |
| Full multivector     | 32         | everything                   | high |

The **bivector field** (10 components) is the most physically natural
choice — it carries the structure of an electromagnetic field tensor
plus angular momentum, and Lengyel's paper on the relativistic inertia
tensor builds momentum-position into a 10-component bivector. Drop to
**relativistic quaternion (8 comp)** for stage A as the closest analog
to our current 6 scalar fields.

We pick **8 components** for stage A.

## Step 2: Component basis

The relativistic quaternion in ℝ(3,1,1) has basis (per Lengyel):

```
e₄₁₀, e₄₂₀, e₄₃₀     (3 components — temporal translations)
𝟙 = e₁₂₃₀₄          (1 component — antiscalar, identity element)
e₁, e₂, e₃           (3 components — spacelike vectors, contracted)
e₃₂₁                 (1 component — spatial pseudoscalar)
```

We label them `M[0..7]`:

| Index | Basis    | Interpretation           |
|-------|----------|---------------------------|
| 0     | e₄₁₀     | temporal translation, x   |
| 1     | e₄₂₀     | temporal translation, y   |
| 2     | e₄₃₀     | temporal translation, z   |
| 3     | 𝟙        | scalar identity           |
| 4     | e₁       | spatial, x                |
| 5     | e₂       | spatial, y                |
| 6     | e₃       | spatial, z                |
| 7     | e₃₂₁     | spatial pseudoscalar      |

A useful split: `M = (q, s)` where `q = M[0..3]` is the "rotor part"
(close to ordinary quaternion) and `s = M[4..7]` is the "boost +
pseudoscalar part". Lengyel's bulk norm:

```
|M|²_bulk = s_x² + s_y² + s_z² + s_w²  −  m_x² − m_y² − m_z² − m_w²
         = M[4]² + M[5]² + M[6]² + M[7]²  −  M[0]² − M[1]² − M[2]² − M[3]²
```

This has Minkowski signature `(+ + + + − − − −)` — four spacelike,
four timelike. Subluminal motion has `|M|²_bulk > 0`.

## Step 3: Linear Lagrangian

The free field Lagrangian (Lorentz invariant by construction):

```
ℒ_free = ½ ∂_μ M^a ∂^μ M^a g_ab − ½ m² |M|²_bulk
```

where `g_ab` is the diagonal metric `(+,+,+,+,−,−,−,−)` matching the
bulk-norm signs. Equation of motion (Klein–Gordon for each component
with the appropriate sign):

```
(□ + m²) M[i] = 0   for i ∈ {4, 5, 6, 7}    (spacelike, +mass²)
(□ − m²) M[i] = 0   for i ∈ {0, 1, 2, 3}    (timelike, tachyonic free; will be stabilised by interaction)
```

The negative-mass-squared on the time-like components is exactly the
condition that gives the Higgs mechanism a nontrivial vacuum: the field
spontaneously settles to a configuration where `|M|²_bulk = const ≠ 0`.

## Step 4: Self-interaction

The simplest Lorentz-invariant nonlinear potential:

```
U(M) = ¼ λ (|M|²_bulk − v²)²
     = ¼ λ (M[4]² + M[5]² + M[6]² + M[7]² − M[0]² − M[1]² − M[2]² − M[3]² − v²)²
```

This is the **double-well Higgs potential**. It has minima on the
hyperboloid `|M|²_bulk = v²` — a 7-dimensional vacuum manifold. The
field naturally settles onto this hyperboloid, and the residual modes
are massless (Goldstone) plus a single massive (Higgs) mode.

Full Lagrangian:

```
ℒ = ½ ∂_μ M^a ∂^μ M^a g_ab − ¼ λ (|M|²_bulk − v²)²
```

Equation of motion:

```
□ M[i] = − λ (|M|²_bulk − v²) × g_{ii} × M[i]    (no sum over i)
```

For `i ∈ {4,5,6,7}`: `g_{ii} = +1` so RHS pulls toward |M|² = v².
For `i ∈ {0,1,2,3}`: `g_{ii} = −1` so the same potential pushes away
from origin — exactly what we want for spontaneous symmetry breaking.

## Step 5: Topological charge

The vacuum manifold is the hyperboloid `H⁷₃ = SO(4,4)/SO(4,3)`. Its
homotopy group `π₃(H⁷₃) = π₃(S⁷) = 0` (trivial) — so this isn't going
to give us topological solitons by itself.

If we restrict to the unit sphere `|M|² = 1` with positive-definite
metric, we get `S⁷` again — `π₃(S⁷) = 0`. Still trivial.

If we go down to `|M|² = 1` with one signature flip (e.g. drop one
component, work on `S⁶`), we have `π₃(S⁶) = 0`. Trivial.

`π₃(S³) = Z` (Skyrme winding) requires the field to be 4-component.
The natural sub-field is the *quaternionic part* `q = M[0..3]` with
its own induced norm. Imposing `|q|² = 1` gives a unit quaternion field,
which IS the Skyrme baryon number.

**Strategy**: use the full 8-component dynamics, but include a Skyrme-
like topological term that couples to the rotor part `q` only. The
field is rich enough for spinor structure; the topological charge is
inherited from the quaternionic Skyrme model that's already
well-understood.

## Step 6: Discrete operators on the foam mesh

The Voronoi foam discretisation from v55 carries over verbatim. Per
cell `c`:

**Laplacian** (unchanged from v55, two-point flux):

```
∇²M[i](c) = (1/V_c) Σ_{f∈faces(c)} A_f × (M[i](nb_f) − M[i](c)) / d_f
```

**Bulk norm** (per cell):

```
|M|²_bulk(c) = Σ_{i∈{4,5,6,7}} M[i](c)² − Σ_{i∈{0,1,2,3}} M[i](c)²
```

**Velocity-Verlet update** for each component i:

```
v[i](c) ← v[i](c) + ½ dt × a[i](c)
M[i](c) ← M[i](c) + dt × v[i](c)
recompute bulk_norm[c] and a[i](c) for all i
v[i](c) ← v[i](c) + ½ dt × a[i](c)
```

where `a[i](c) = ∇²M[i](c) − λ (|M|²_bulk(c) − v²) × g_{ii} × M[i](c)`.

**No curl coupling** in stage A — the multivector self-couples through
the bulk norm, not through curl exchange. (Curl-style exchange would
appear naturally in stage B with the Dirac equation, where it becomes
the γ^i ∂_i ψ structure.)

## Step 7: Initial condition (seed)

For the Higgs-like field with vacuum on `|M|²_bulk = v²`, a natural
seed is a **defect** — a region where the field is forced off the
vacuum manifold:

- **Q-ball seed**: localised oscillation in the rotor part
  `q(x) = (1, 0, 0, 0) × A · exp(−r²/2R²) · cos(ωt)`
  with the rest of M at the vacuum.
- **Vortex seed**: `q(x) = (cos(θ), sin(θ), 0, 0)` where θ is the
  azimuthal angle in xy-plane. This has unit Hopf-like winding.
- **Skyrme hedgehog**: `q(x) = (cos(f(r)), n̂(x) sin(f(r)))` with
  `f(0) = π`, `f(∞) = 0`. Winding number 1.

Stage A starts with the Q-ball seed (simplest, no topology) to verify
the dynamics work. Stage B adds the Skyrme hedgehog once the kernel
runs.

## Step 8: Observables

Same as v55 + new ones from the algebra:

- **|M|²_bulk(c)** per cell — should oscillate around v² for vacuum,
  deviate inside the soliton
- **Energy density** `½ Σ_i (∂_t M[i])² + ½ Σ_i |∇M[i]|² + ¼ λ (|M|²_bulk − v²)²`
- **Skyrme winding** `B = (1/24π²) ∫ ε_ijk Tr(L_i L_j L_k)` where
  `L_i = q⁻¹ ∂_i q` (only meaningful when |q|² > 0)
- **Bulk-norm distribution histogram** — to see whether the field
  stays on the vacuum manifold globally

## Step 9: What's the same vs different from v55

| Aspect | v55 (Cosserat) | v56 stage A (mv field) |
|--------|---------------|------------------------|
| Components per cell | 6 (φ₃ + θ₃) | 8 (relativistic quat) |
| Field representation | scalar SO(3) vector | spinor-capable Lorentz multivector |
| Laplacian operator | unchanged FV | unchanged FV (same code) |
| Mass term | +m² φ + 0 × θ | g_{ii} × m² × M[i] (mixed signs) |
| Self-coupling | V(P) = μ P²/(1+κP²), P = φ₀φ₁φ₂ | λ (|M|²_bulk − v²)² |
| Curl coupling | η ∇×θ in φ EOM | none |
| Vacuum | A_bg cos(k·z) carrier | uniform |M|²=v² hyperboloid |
| Particle candidate | breathing oscillon (unstable) | Q-ball, Skyrme hedgehog |
| Time integrator | velocity Verlet | velocity Verlet (same) |
| Kernel changes | — | bump n_fields 6→8, change potential, drop curl |

The kernel structure is essentially preserved. The biggest change is
the **potential** — replacing V(P)=(μ/2)P²/(1+κP²) with the Higgs
quartic (|M|² − v²)². The physics interpretation is completely
different (broken symmetry, Goldstone modes, topological solitons).

## Step 10: First experiment design

Once the v56 kernel runs:

1. **Vacuum verification**: start at |M|² = v² uniformly. Should be
   stationary (only zero-amplitude oscillations).
2. **Tachyonic instability check**: start at M=0. Should grow
   exponentially in the rotor part until |M|² ≈ v².
3. **Q-ball test**: localised perturbation of the rotor part. Should
   either persist as a Q-ball or radiate away.
4. **Skyrme hedgehog**: hedgehog seed with winding number 1. Should
   converge to a Skyrme-like profile if the topological term is
   present.
5. **Two-soliton interaction**: same as v55 stage but with the new
   field type. Compare attraction/repulsion to fermion expectations.

## Why this might succeed where v54 didn't

v54's null result was 750+ configurations of the 6-field Cosserat
producing no stable particles. Two specific reasons:

1. **No spinor structure** — the 6-field Cosserat has SO(3) rotation
   symmetry but no spinor representations. v56 has spinor structure
   built into the algebra (relativistic quaternions are Pin(3,1)
   double covers).
2. **No vacuum manifold with nontrivial topology** — v55's "vacuum"
   is a single field configuration (A_bg carrier wave). v56's vacuum
   is a 7D hyperboloid; field configurations can wind around it,
   carrying topological charge that physics can't dissipate.

Stage A only delivers (2) (Higgs-like vacuum, Q-balls if any).
Stage B (Dirac) is needed for (1).
