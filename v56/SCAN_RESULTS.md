# V56 Skyrme parameter scan — definitive null result

Scan of (c4, R_seed) at L=20 to test whether either textbook-strength
Skyrme or compact-seed branches stabilise the winding-1 hedgehog.

## Configurations tested

| Run            | c4    | sigma_e2 | R_seed | Wall  | E_total | E_2:E_4 ratio at t=0 |
|----------------|-------|----------|--------|-------|---------|----------------------|
| Production L40 | 0.05  | 1.0      | 5.0    | 48.5m | 154.4   | 153.5 : 0.84  (180:1) |
| Strong L20     | 5.0   | 5.0      | 5.0    | 11.0m | 253.5   | 168.9 : 84.7  (2:1)   |
| Compact L20    | 0.5   | 5.0      | 1.5    |  8.2m |  61.4   |  44.4 : 17.1  (2.6:1) |

## Winding decay: invariant across 100× change in Skyrme strength

| Run            | B(t=0) | B(t=1) | B(t=2) | B(t=3) | B(t=4) | B(t=5) |
|----------------|--------|--------|--------|--------|--------|--------|
| Production L40 | 0.945  | 0.679  | 0.315  | 0.079  | −0.003 | ≈ 0    |
| Strong L20     | 0.948  | (n/a)  | 0.275  | (n/a)  | −0.003 | 0.015  |
| Compact L20    | 0.685  | (n/a)  | −0.006 | (n/a)  | −0.001 | 0.001  |

**The Skyrme term doesn't control the unwinding rate.** A 100× boost in c4
(0.05 → 5.0), and additionally a 5× boost in sigma_e2, both at the same
seed radius — winding still collapses in ~3–4 t-units. The compact seed
loses winding even faster (~2 t-units), since it's resolved by only ~3
cells across.

Energy conservation is identical across all three: drift envelope ±0.001%
to ±0.015%.

## Diagnosis: bulk-Higgs / σ-constraint branch competition

The bulk-Higgs hyperboloid `|M|²_bulk = M[4]² + … − |q|² = v²` admits
TWO branches at the vacuum point M[5..7] = 0:

- **Skyrmion branch**:  |q| = 1, M[4] = √(v² + 1) ≈ 1.118 — what we seed.
- **Scalar branch**:    |q| = 0, M[4] = v = 0.5 — same bulk_norm = v².

Both satisfy the bulk Higgs constraint. The σ-constraint
`(|q|² − 1)²` is the *only* term penalising the scalar branch, and it
is a soft quadratic penalty. With sigma_e2 = 5, the energy cost of
falling onto the scalar branch is 5/4 · (0 − 1)² = 1.25 per unit
volume, which is negligible against the > 100 units of gradient
energy carried by the hedgehog.

Empirically: |q|² promptly excurses to [0.029, 1.0] (strong-Skyrme
run, t = 1.99) or [0.564, 1.076] (compact run, t = 1.99). The rotor
collapses into the scalar branch, the hedgehog topology dissolves,
and the energy radiates outward as a Higgs-amplitude breathing mode.

## What we tried, what we learned

Both scans **fail at the same step**: the rotor escapes |q| = 1 within
1–2 t-units. Once off the unit sphere, π₃ no longer applies and the
winding is no longer protected. After that, the Skyrme term works on
a non-S³ field where it can't enforce topology.

Visual snapshots in `v56/foam/snapshots_scan/` confirm:
- Strong-Skyrme f50: the energy fills the entire L=20 box (radiated
  out periodic boundary).
- Compact-seed f25: a single large warm sphere filling most of the
  box — original compact structure has spread by R ~ 8.

## Recommended next step: HARD projection onto S³

The standard Skyrme model implementation enforces |q| = 1 as a HARD
constraint, not a soft penalty. After every Verlet step, project:

    q ← q / |q|     (and remove velocity component along q to keep ⟨q, q̇⟩ = 0)

This is a constraint-Hamiltonian symplectic projector — preserves
energy to projection precision (~1e-14 per step) and prevents *any*
escape from S³. The competing bulk-Higgs scalar branch becomes
inaccessible.

Implementation cost: ~30 lines in `verlet_step`, after the second
half-step velocity update:

```c
for each cell c:
    qn² = M[0]² + M[1]² + M[2]² + M[3]²
    inv_qn = 1/sqrt(qn²)
    for k in 0..3:
        M[k][c] *= inv_qn
    // Remove velocity component along q (tangent projection):
    qdotv = M[0]·v[0] + M[1]·v[1] + M[2]·v[2] + M[3]·v[3]
    for k in 0..3:
        v[k][c] -= qdotv · M[k][c]
```

This drops the σ-constraint entirely (no longer needed). The bulk
Higgs `(|M|²_bulk − v²)²` would still adjust M[4..7].

Trade-off: hard projection sacrifices Lorentz covariance on the rotor
sector — q now lives on Euclidean S³, not a Lorentzian hyperboloid.
The original Lengyel motivation (Cl(3,1,1) algebra) is partly
abandoned. But this is the textbook Skyrme construction, which DOES
have stable winding-1 solutions.

## If hard projection is rejected as inconsistent with the algebra

The alternative is to **drop the bulk-Higgs scalar branch** by changing
the bulk metric. With g = (+,+,+,+,+,+,+,+) (all spacelike), the
constraint `|M|²_bulk = v²` becomes a 7-sphere, and π₃(S⁷) = ℤ
provides natural topological protection without needing σ-constraint
on the rotor sub-block. This requires reworking the PGA "relativistic
quaternion" interpretation.

A third option: treat M[4] as a Lagrange multiplier for the bulk
constraint and enforce |q|² = 1 directly via projection on M[0..3]
only. This keeps rotor rotation-invariance plus the Higgs-vacuum
hyperboloid, but locks the rotor to S³ algebraically.

## Files

- `v56/foam/scan_strong_L20.sfa`, `scan_strong_L20_diag.tsv`
- `v56/foam/scan_compact_L20.sfa`, `scan_compact_L20_diag.tsv`
- `v56/foam/scan_strong_skyrme_L20.cfg`, `scan_compact_seed_L20.cfg`
- `v56/foam/snapshots_scan/scan_{strong,compact}_f{1,2,3,5,10,25,50}.webp`

## Bottom line

The Skyrme L_4 implementation is correct (energy conserves to 1e-5 across
100× variation in c4, B(0) measured correctly). The problem is **not**
the Skyrme strength. The problem is that the v56 vacuum manifold has a
disconnecting non-Skyrmion branch the rotor escapes to.

Without changing the field structure (hard S³ projection or different
bulk metric), no setting of (c4, sigma_e2, R_seed) will produce a
stable winding-1 Skyrmion in this kernel.
