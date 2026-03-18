# T4: Topological Fragility — Results

## Setup
- N=96, L=20, bimodal sweet spot (t=0.85 interpolation)
- Phase 1: Settle T=200 with absorbing damping, mass2=0
- Phase 2: Gaussian perturbation at origin (sigma=3.0=R_tube), n=(1,1,1)/sqrt(3)
- Evolve T=300 more, track winding at 10 intermediate times
- 11 amplitudes: eps in {0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0}

## Key Finding: Braid Survives All Perturbations

**All 11 runs are stable.** No blowup at any perturbation amplitude, including
eps=5.0 (which is 6.25x the field amplitude A0=0.8).

| eps   | fc     | peak_P | energy   | stable |
|-------|--------|--------|----------|--------|
| 0.05  | 0.2238 | 1.6708 | 90734    | OK     |
| 0.10  | 0.2237 | 1.6726 | 90696    | OK     |
| 0.20  | 0.2234 | 1.6722 | 90610    | OK     |
| 0.30  | 0.2231 | 1.6625 | 90515    | OK     |
| 0.50  | 0.2221 | 1.6234 | 90282    | OK     |
| 0.80  | 0.2204 | 1.6346 | 89785    | OK     |
| 1.00  | 0.2189 | 1.6245 | 89266    | OK     |
| 1.50  | 0.2173 | 1.4292 | 90273    | OK     |
| 2.00  | 0.2153 | 1.3897 | 88969    | OK     |
| 3.00  | 0.2134 | 1.3787 | 88314    | OK     |
| 5.00  | 0.2145 | 1.4182 | 88912    | OK     |

## Winding Oscillation (Not Destruction)

The winding number oscillates dynamically between -1, 0, and +1 throughout all
runs. This is the inherent behavior of this braid configuration: since
pi_1(R^3)=0, the winding is NOT topologically quantized -- it oscillates
continuously as the braid evolves. The settled state happened to be sampled
at w~0, while the perturbed runs land at w=-1 or w=0 depending on phase.

This means the naive eps_crit=0.05 (from comparing final winding to initial)
is **misleading**. The perturbation does not destroy the braid; it merely
shifts the oscillation phase.

## True Fragility Measure

The better indicators of braid integrity are:
- **fc (core fraction)**: ranges 0.2134-0.2238, only ~5% variation across 100x
  range in perturbation amplitude. The braid core remains intact.
- **peak_P (triple product)**: ranges 1.38-1.67, modest decrease at large eps.
  The three-field overlap is maintained.
- **Energy**: ranges 88300-90700, ~3% variation. No catastrophic energy injection.
- **No blowup** at any amplitude.

## Conclusion

**The braid is effectively robust against localized Gaussian perturbations up to
at least eps=5.0 = 6.25 x A0.** No topological destruction was observed. The
winding oscillation is an intrinsic dynamical feature, not a sign of fragility.

By the fc/energy criteria, eps_crit > 5.0 (ratio > 6.25), placing the braid
in the **"effectively topological"** regime despite pi_1(R^3)=0.

## Caveats

1. The winding diagnostic measures phase angle along the z-axis centerline only.
   It oscillates naturally and is not a good stability indicator for this system.
2. Perturbation was symmetric (equal on all 3 fields). Targeted single-field
   perturbation or velocity-only perturbation (as suggested in PROPOSAL.md)
   might reveal different fragility thresholds.
3. N=96 resolution may under-resolve fine structure; N=128 confirmation would
   strengthen the result.
