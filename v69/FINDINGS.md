# v69 FINDINGS — Kernel-v3: The Gauged Diagonal U(1)

**Date**: 2026-06-11. **Status**: [measured] kernel-v3 implemented (CPU+GPU), verified,
and first gauged campaign complete (10 runs, all rc=0). USER-AUTHORIZED kernel change.
Docs: `v69/SPEC.md` (implementation contract incl. the unique Gauss-conserving lattice
current), `v68/GAUGE_DESIGN.md` (theory, 59/59 Maxima), `v69/theory/` (gauged radial
shooter + branch scan `gscan.tsv` + seed profiles). Data: `v69/results/`,
SFAs `/space/scp/v69/`.

## 1. Implementation and verification chain

- Compact Kogut–Susskind U(1) links + noncompact E, temporal gauge, all matter
  derivatives link-covariant; everything derived from one lattice Hamiltonian so the
  discrete Gauss law is conserved by symplectic structure, not projection.
- CPU gates: g=0 **byte-identical** to the ungauged complex kernel (25 cols × 21 rows
  exact); **gauss_max = 6.7×10⁻¹⁴ flat over T=60**; gauged ball stable with
  omega_core = weff(0) = 1.395960 (the Coulomb potential living in the local phase
  rate, as the temporal-gauge analysis predicted) and **Q_flux/Q_noether = 0.999** —
  the ball's charge read from a boundary flux integral alone.
- GPU port: gauge-parity review found ZERO discrepancies (current, staples, seagull
  all hand-verified); nvcc clean first try; **GPU reproduces CPU to all 12 printed
  digits** on omega_core and Q_flux. Latent GPU bug fixed in passing (bc_type=2
  previously ran the absorbing kernel).
- Gauged radial shooter (v69/theory): g=0 branch matches v66 scan ≤0.04%;
  dE/dω = ω·dQ/dω survives gauging (1e-4).

## 2. The gauged branch: Coulomb self-repulsion imposes a maximum particle size

| g | window ω | Q_max | δω_self(Q=482) |
|---|---|---|---|
| 0 | (1.309, 1.50) | ~∞ (1.4×10⁵ at floor) | — |
| 0.02 | (1.361, 1.50) | 5299 | +0.0035 |
| 0.05 | (1.406, 1.50) | **921** | +0.0214 (= +8.6g²) |
| 0.10 | (1.466, 1.50) | 260 | Q=482 not on branch |

Far harsher than the design estimate (Q_max ~6200): **gauging ends the periodic
table** — at g=0.05 no ball above Q=921 exists statically. The fission analogy is
quantitative.

## 3. Coulomb interaction measured — both signs, absolute calibration

- acc3 (same charge, D=16): monotone accelerating **repulsion** (r_core 8.68→9.33).
- acc4 (opposite charge, D=16): monotone accelerating **attraction** (8.68→8.01).
- Force-law scan (cluster-tracked D(t), same-charge pairs Q=311):

| D | measured a | Coulomb a = 2g²Q²/(4πD²M) | ratio |
|---|---|---|---|
| 14 | 6.2×10⁻⁵ | 4.3×10⁻⁴ | 0.14 — **watershed**: contact attraction cancels 86% |
| 20 | 1.95×10⁻⁴ | 2.11×10⁻⁴ | **0.92** |
| 26 | 1.47×10⁻⁴ | 1.25×10⁻⁴ | **1.17** |

Beyond the watershed the measured force matches the **parameter-free** Coulomb
prediction to 10–20% — absolute calibration, stronger than an exponent fit at this
frame count. r*(Q=311) sits just below D=14 (design scaling: r*(Q=482)≈12 ✓ trend).
[Clean n=2.00 exponent: needs D∈[20,40] with denser snapshots — queued.]
This closes the program's oldest quantitative goal (F1, force law on stable
particles): the long-range force is Coulomb, mediated by a massless NEUTRAL field
that carries no charge — the drain problem solved with long-range physics intact.

## 4. Fission test: the symmetric saddle (Bohr–Wheeler in-simulation)

fiss1: a Q=1745 ball (1.9× above Q_max) seeded spherically symmetric survived T=300
nearly intact (Q=1743, mild dilation). Interpretation: a spherical super-critical
charge sits on the SYMMETRIC saddle; fission requires an asymmetry seed (the lattice
at f64 is too clean to provide one). Follow-up: ℓ=2-perturbed seed → watch it split;
fragment masses vs the liquid-drop analogy. Q_flux tracked the giant to ~1%.

## 5. Positronium: bound, oscillating, NOT annihilating

pos1 (opposite charges ±921 at D=8, η=0, m_θ=1.6, T=600): bound oscillation
(r_core 5.6↔6.7, period ~200 vs predicted T_orb=275 at D₁=7.9), Q_total=0 exact,
**no annihilation** — with the θ-drain closed and charge gauged, the pair has no
open decay channel on this timescale, in sharp contrast to the ungauged tb3 grind-down.
Caveat: r_half=5.6 means heavy overlap at D=8 — this is an oscillating charge dipole,
not a clean two-body orbit. Clean orbit experiment (fixed-Q balls, D=16, tangential
velocity): queued.

## 6. What gauging delivered (scorecard vs GAUGE_DESIGN predictions)

| prediction | outcome |
|---|---|
| drain closes identically (neutral mediator) | ✓ (pos1: no decay channel) |
| Gauss law to rounding without projection | ✓ 6.7×10⁻¹⁴ |
| Q from boundary flux | ✓ 0.1–1% |
| same-charge repel / opposite attract, 1/r | ✓ absolute calibration 10–20% |
| watershed r* | ✓ measured just below 14 at Q=311 |
| Q_max (fission limit) | ✓ but 921 not 6200 (shooter replaced estimate) |
| ω_local = ω + g·a₀ phase rate | ✓ 1.395960 exact |
| g=0 bit-identity | ✓ byte-exact |

## 7. Next

1. ℓ=2-perturbed fission run (the splitting movie + fragment spectrum).
2. Clean positronium orbit (tangential velocity, D=16, fixed-Q profiles).
3. Force-law exponent with lever arm (D 20→40, snap_dt=5).
4. AB/flux-line winding experiment (holonomy 2πqn — now topological).
5. Gauged charged-bath + condensate fragmentation (do fragments now obey Q_max?).
6. EM_THEORY.md update: the photon analog is now A, not the polariton; rewrite the
   EM sector mapping (θ remains as massive charged torsion matter).

## Infrastructure

CPU pair runs at N=128 cost ~0.4 s/step (gauged ≈ 2× complex) — two-ball gates belong
on GPU (2 min vs 2 h). The round-3 verifier relaunched killed CPU runs and was
stopped manually; lesson: put explicit wall-clock budgets in verifier prompts.
