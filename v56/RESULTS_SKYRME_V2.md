# V56 Stage C — Real Skyrme L_4 on L=40 mesh

First production run with the genuine Skyrme term `(1/32 e²) Tr([L_μ, L_ν]²)`
implemented as a face-divergence of the conjugate momentum

    J^{a,b} = 2 c4 [Tr(g) ∂_b q^a − Σ_j g_bj ∂_j q^a]

with c4 = 0.05 (≈ 1/(2 e²) for e ≈ 3.16, near textbook). Energy density
`e_4 = (c4/2) [(Tr g)² − ‖g‖²_F]` where `g_ij = Σ_a ∂_i q^a ∂_j q^a` over
rotor M[0..3]. Seed: |q|=1 hedgehog `q = (cos f, n̂ sin f)`,
`f(r) = π exp(−r/R)`, R_seed = 5.0. M[4] = √(v² + 1) so bulk-Higgs
vacuum holds at every cell at t=0. T = 50, dt_factor = 0.020.

Wall: 48.5 min on 8 CPU threads, 3.26M cells.

## Energy conservation: PERFECT

| t | E_total | drift |
|---|---------|-------|
| 0     | 1.5436e+02 | +0.0000% |
| 0.99  | 1.5436e+02 | −0.0003% |
| 5     | 1.5436e+02 | +0.0003% |
| 25    | 1.5436e+02 | −0.0003% |
| 50    | 1.5436e+02 | −0.0004% |

Drift envelope across 50 time units: **±0.0008%** — confirms the
discrete divergence of J on faces is the exact adjoint of the
cell-centered gradient operator (summation by parts). This is the
tightest energy conservation we've achieved on the foam kernel.

## Winding number: COLLAPSES IN ~4 t-units

| t | B | |q|² range | E_skyrme | E_grad |
|---|---|-----------|----------|--------|
| 0    | **+0.945** | [1.000, 1.000] | 0.84   | 153.5 |
| 0.99 | +0.679 | [0.573, 1.000] | 0.57   | 136.3 |
| 1.99 | +0.315 | [0.124, 1.000] | 0.25   | 103.0 |
| 2.98 | +0.079 | [0.148, 1.000] | 0.08   | 67.8  |
| 3.98 | −0.003 | [0.976, 1.101] | 0.01   | 37.9  |
| 5–50 | ≈ 0    | [0.6, 1.5]      | < 0.01 | 60–70 |

The B = 0.945 at t=0 (vs expected 1.0) reflects the finite-volume
discretisation of the determinant integrand at the box boundary
(f(L=40) = π·e⁻⁸ ≈ 0.001 → tiny S³ deficit at the boundary).

The hedgehog **unwinds completely within 4 t-units**. This is robust
across both the diagnostic (B → 0) and the visual (rotor structure
gone by frame 5).

## Why it unwinds: Derrick virial mismatch

At t=0:
- E_2 (gradient) = 153.5
- E_4 (Skyrme)   = 0.84

The Skyrme L_4 contribution is **180× too small** to balance E_2.
Skyrmion stability requires the Derrick virial E_2 = E_4 at the
soliton scale. The current parameters have:

- c4 = 0.05 (likely needs ~50× more, or ~2.5 to balance)
- R_seed = 5.0 (gradients are shallow → small E_4); halving R doubles
  the gradient and quadruples E_4 (since L_4 ∝ (∂q)⁴)

The σ-constraint (sigma_e2 = 1.0) is also too weak — |q|² excurses to
[0.124, 1.978] during the tachyonic collapse window (t = 1–5).
The constraint penalty `(|q|² − 1)²` can't hold the rotor to S³ when
the Skyrme term is small.

## Visual: a coherent late-time breather

Despite winding loss, the energy stays **localised** at the seed site
throughout T = 50. By t = 48 the system has settled into a clean
concentric-ring pattern (snapshots in `foam/snapshots_skyrme_v2/`):

- **t = 0**: warm hedgehog core (rotor M[0..2] active).
- **t = 3**: still localised, structure simplified.
- **t = 8**: incipient ring structure forming.
- **t = 13–23**: green/yellow concentric rings around a bright core.
- **t = 48**: clean radial standing-wave pattern, M[3]+M[4] dominant
  (green/cyan), with warm M[0..2] halo at the outer rings.

This is **qualitatively different** from the OLD broken-Skyrme run
(no real L_4, |q|=A=0.2 partial S³, `(|q|²−1)²` penalty as the only
"Skyrme" term). That run dispersed into a chaotic sprawl filling the
entire box by t = 35. The new run holds the energy in a coherent
localised breathing mode — the σ-constraint + L_4 together keep the
system from thermalising.

We do not have a stable winding-1 Skyrmion yet. We do have:
- A tightly conserved Hamiltonian dynamics with the new infrastructure.
- A localised non-topological "Higgs breather" oscillating around
  bulk vacuum.

## Diagnostic side notes

- bulk_norm stays in [0.20, 0.30] for t > 5 (vs v² = 0.25) — bulk
  Higgs vacuum is well-respected outside the soliton site.
- E_kin / E_grad oscillate around 77 / 63 with no decay — system is
  in a long-lived bound oscillation, not relaxing toward vacuum.
- M_max for M[0..3] stays bounded (no rotor blow-up).

## Next-step parameter recommendations

For a stable winding-1 Skyrmion the Derrick virial must be satisfied
at *some* preferred radius R*. The standard estimate (textbook
Skyrme) gives R* such that E_2(R*) = E_4(R*); fixing E_2 ≈ 153 at
R = 5 implies E_4(R*) needs c4 ≈ 9 (factor ~180 boost), OR R needs
to drop to ~1.5 with c4 ≈ 0.5.

Two targeted follow-up experiments:

1. **Strong Skyrme**: c4 = 5, sigma_e2 = 5, R_seed = 5.0. Tests
   whether a textbook-strength L_4 stabilises the long-wavelength
   seed.
2. **Compact seed**: c4 = 0.5, sigma_e2 = 5, R_seed = 1.5. Tests the
   alternative branch (small Skyrmion at lattice resolution edge).

Both should be run at L = 20 first (smaller mesh, ~10 min wall) to
scan c4 and R_seed before committing to L = 40.

## Files

- `v56/foam/v56_skyrme_L40.sfa` (2.3 GB, 50 cell-native frames)
- `v56/foam/v56_skyrme_L40_diag.tsv` (50 rows)
- `v56/foam/snapshots_skyrme_v2/` (8 PGA-spectrum frames at t=0..48)
- `v56/foam/v56_skyrme_L40_OLD.sfa` (3.3 GB, original broken-Skyrme
  baseline kept for reference)
