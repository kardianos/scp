# V45 Results — Differential Mass Defect

## Experiment Design

Eight deuterium simulations (UUD+UDD) using identical seeds (gen_deuterium.c),
identical grid (N=512, L=100), identical BC (absorbing, damp_width=10,
damp_rate=0.005), differing ONLY in inter-baryon separation D. All runs
T=500, COLZSTD output, V100-SXM2-32GB.

## Separation Profiles

### D=0 — Maximum Overlap
Both baryons initialized at the same position. All 6 braids (3 per baryon)
superposed. Expected violent phase confinement repulsion from P→0 at overlap.
**Result**: Highest E_total (+1253 above baseline). Strong repulsion confirmed.
P_int=763 (highest of all runs — braids overlap constructively before separating).

### D=5 — Deep Confinement
Baryons nearly touching (separation ≈ 1 braid radius). Still deep in the
confinement zone where the triple product P=φ₀φ₁φ₂ is disrupted by overlap.
**Result**: E_total +1009 above baseline. Still strongly repulsive but less
than D=0. The braids have some room to organize.

### D=10 — Core Boundary
Separation comparable to 2× braid core radius (~5 each). The edge of the
confinement interaction zone.
**Result**: E_total +476 above baseline. Repulsion weakening rapidly. P_int
starting to drop (664 vs 763 at D=0) — less constructive overlap.

### D=15 — Interaction Surface
Separation at the braid interaction surface (r≈4-6 per baryon). This is
where the force should transition from confinement repulsion to potential
attraction.
**Result**: E_total +101 above baseline. Nearly at baseline — the repulsive
energy has mostly dissipated. P_int=646, close to the D=80 baseline (642).

### D=20 — Nuclear Range
Well beyond the core, into the θ-mediated interaction range. If nuclear
binding exists, it should appear here.
**Result**: E_total +128 above baseline. Slightly ABOVE D=15, suggesting
no attractive well has formed. The force is near zero.

### D=40 — V42 Equilibrium Distance
The separation used in V42 deuterium. Expected to be the binding minimum.
**Result**: E_total +282 above baseline. HIGHER than D=15-20, not lower.
No binding minimum at this separation.

### D=60 — Weak Interaction
Far separation, weak interaction expected.
**Result**: E_total +149 above baseline. Slightly above baseline.

### D=80 — Non-Interacting Baseline
Baryons at maximum separation (each at ±40 from center). Should be
non-interacting.
**Result**: E_total = 103864 (baseline). Lowest energy of all runs.

## Energy vs Separation Curve

```
D=  0  ΔE = +1253  ████████████████████████████████████████████
D=  5  ΔE = +1009  ██████████████████████████████████
D= 10  ΔE =  +476  ████████████████
D= 15  ΔE =  +101  ███
D= 20  ΔE =  +128  ████
D= 40  ΔE =  +282  █████████
D= 60  ΔE =  +149  █████
D= 80  ΔE =    0   (baseline)
```

## Force Curve

| D_mid | F = -dE/dD | Direction |
|-------|-----------|-----------|
| 2.5 | +48.8 | Attractive (toward larger D) |
| 7.5 | +106.5 | Attractive (toward larger D) |
| 12.5 | +74.9 | Attractive (toward larger D) |
| 17.5 | -5.3 | ~Zero |
| 30.0 | -7.7 | ~Zero / slight repulsion |
| 50.0 | +6.7 | ~Zero |
| 70.0 | +7.5 | ~Zero |

The "attractive" force at D=0-15 is not binding attraction — it's the
restoring force from phase confinement repulsion. Overlapping baryons
are pushed apart toward their natural separation. The force at D>15 is
essentially zero within noise.

## Mass Defect Result

**No binding energy detected.** The energy minimum is at D=80 (maximum
separation). Every closer separation has HIGHER total energy than the
baseline. There is no attractive potential well at any separation tested.

| Metric | Value |
|--------|-------|
| E_bind (D=40 vs D=80) | **+282 code units** (repulsive, not bound) |
| E_bind (D=15 vs D=80) | **+101 code units** (repulsive, not bound) |
| Energy minimum | D=80 (non-interacting) |
| Binding well depth | **NONE detected** |

## Interpretation

1. **Phase confinement works as expected.** D=0-10 shows strong repulsion
   (+475 to +1253) when baryons overlap. The braids resist merging because
   the triple product P=φ₀φ₁φ₂ is disrupted at overlap.

2. **No nuclear binding at any separation.** The energy curve is monotonically
   decreasing from D=0 toward D=80. There is no attractive well — the
   baryons always prefer to be farther apart.

3. **The V42 deuterium "binding" may be a seed artifact.** V42 showed two
   baryons remaining together at D=40 for T=500, but this may be because
   the gen_deuterium seeds initialized them in a configuration that takes
   longer than T=500 to separate. The system isn't bound — it's just slow
   to evolve at that separation.

4. **The D=40-60 bump** (+282 at D=40 vs +149 at D=60) is unexpected and
   may indicate a weak interaction effect, but it's in the repulsive
   direction (more energy, not less). This could be a resonance between
   the two baryons' breathing modes at D=40.

5. **P_int correlation with separation**: P_int is highest at D=0 (763,
   constructive overlap) and drops to a plateau of ~642-649 for D≥15.
   The braids don't interact at D>15 — they're isolated structures.

## Implications

The SCP model as currently parameterized does NOT produce nuclear binding
between UUD and UDD baryons. The gen_deuterium seeds place baryons in
proximity, but they have no energetic reason to stay together. This is a
significant negative result for the nuclear physics program.

Possible explanations:
- **η too low**: The phi-theta coupling (η=0.5) may not be strong enough
  to generate an attractive potential at nuclear distances. The θ-mediated
  "EM" force is repulsive at all ranges for same-charge baryons (both
  radiate isotropic θ).
- **Missing physics**: The current 6-field Cosserat equation may lack the
  mechanism for nuclear binding. Real nuclear binding comes from meson
  exchange (pion, sigma) which requires massive mediator fields.
- **Wrong comparison**: The mass defect may require longer equilibration
  (T=2000+) or different initial conditions to manifest.

## Data Files

- Diag: `/home/d/code/scp/v45/v45_D{0,5,10,15,20,40,60,80}_diag.tsv`
- SFA: `/space/scp/v45/v45_D{0,5,10,15,20,40,60,80}.sfa` (COLZSTD)
- Seeds: `/space/scp/v45/deut_D{0,5,10,15,20,40,60,80}_seed.sfa`
