# v88 — choosing the discrete fabric

> **SUBORDINATE TO `v88/PRINCIPLE.md`.** That document is foundational:
> energy is never destroyed, space is a *mode* of energy, matter is
> converted space, and curvature follows from conservation. Anything here
> that assumes a background — a permanent index set with evolving contents
> — is superseded by it, including every instrument built so far.

**Date:** 2026-07-27 · Instruments: `v86/n_battery/fabric_aniso.c`,
`v88/lattice_iso.c`, `v88/fabric_pn.c`, `v88/fabric_iso.c`, `v88/fabric_proc.c`.

Tags: **[D]** derived · **[M]** measured · **[C]** conjecture

---

## 0. Why the fabric has to change

The census's objects are fluid Q-balls: continuous families, no mass
quantisation, zero bound internal modes for ω ≥ 1.36, and — measured in
`v86/NEXT_PROGRAM_RESULTS.md` §10 — a stability criterion that selects an
*interval*, not a discrete set. Homogeneity is intrinsic to a continuum PDE, so
continuous families are what it must produce.

The project has always run at **structure width / fabric spacing ≈ 10–20**,
where the lattice is invisible: D8b measured the Peierls–Nabarro pinning at
R/E ≈ 1.2×10⁻⁸ there. So the fabric has been *effectively continuous*, and
"refine dx until converged" — the criterion gating most of v86 — pushes away
from discreteness rather than toward it. Discrete-fabric physics lives at
width/spacing ≈ 1–3.

---

## 1. A cubic fabric cannot be made discrete without becoming crystalline  **[M]**

`fabric_aniso.c`, analytic bump slid through one lattice period along [100],
[110], [111]. (An earlier version used the tabulated radial profile; its
interpolation nodes move with the bump and produced a spurious dx² law that
swamped the real signal. With a closed-form bump the barrier is properly
**exponential** in width/dx — 9.1×10⁻¹⁰ at ratio 7.3 to 1.7×10⁻³ at ratio 1.5.)

| width/dx | PN[100]/E | PN[111]/E | anisotropy |
|---|---|---|---|
| 7.33 | 9.1×10⁻¹⁰ | 2.7×10⁻⁹ | 66.5% |
| 3.30 | 6.9×10⁻⁶ | 3.3×10⁻⁵ | 78.9% |
| 2.20 | 5.1×10⁻⁴ | 9.8×10⁻⁴ | 48.5% |
| 1.50 | 1.7×10⁻³ | 1.9×10⁻³ | 40.6% |

**The anisotropy never goes away — it is ~½ to ⅔ of the barrier at every
ratio.** Discreteness and crystallinity arrive in fixed proportion. There is no
regime where a simple-cubic fabric is discrete but isotropic.

---

## 2. Crystalline tuning fixes dispersion, not pinning  **[D,M]**

Anisotropy of the *linear* operator lives in the 4th moment
`T4 = Σⱼ wⱼ dⱼ^⊗4`; a neighbour set is 4th-order isotropic iff it is a
spherical 4-design. The octahedron (SC's 6) and cube (BCC's 8) are only
3-designs. But a cubic-symmetric T4 has exactly **one** independent anisotropy,
so **one** shell-weight ratio cancels it exactly while keeping periodicity
(`lattice_iso.c`):

| fabric | nbrs | A4 (3 = isotropic) | spread @ \|k\|a=2 |
|---|---|---|---|
| SC 6 | 6 | ∞ | 22.69% |
| BCC 8 | 8 | 1.000 | 17.73% |
| FCC 12 | 12 | 2.000 | 7.38% |
| **BCC 8+6, w₂/w₁ = 2/3** | **14** | **3.0000** | **1.49%** |
| FCC 12+6, w₂/w₁ = 1/4 | 18 | 3.0000 | 1.55% |
| icosahedral 12 (not periodic) | 12 | 3.0000 | 0.68% |

The tuned weights are exact rationals — the cancellation is linear in the
weight, so they are analytic, not fitted. **BCC 8+6 is the Kelvin foam:** the
truncated-octahedron cell has exactly those 14 faces. "Foam lattice" and
"crystalline grid" are the same fabric.

**But it does not transfer.** `fabric_pn.c`, same fabrics, nonlinear pinning:

| width/a | SC anisotropy | Kelvin anisotropy |
|---|---|---|
| 2.00 / 1.59 | 72.9% | 13.3% |
| 1.80 / 1.43 | 34.3% | 25.5% |
| 1.60 / 1.27 | 38.9% | **56.9%** |
| 1.45 / 1.15 | 53.9% | 49.2% |
| 1.30 / 1.03 | 62.0% | 18.8% |

15× on dispersion, ~1.6× on pinning, and at one width Kelvin is *worse*. The
reason is structural: pinning is set by coherent Bragg reflection off the
**reciprocal lattice**, and stencil weights cannot move the reciprocal lattice —
only the geometry can. (A simple |G_min| estimate predicts BCC *below* SC and
the data shows the opposite, so the mechanism is probably the reciprocal-star
geometry but the simple version does not reproduce the magnitudes. **Not
established.**)

---

## 3. The procedural irregular fabric  **[M]**

`fabric_proc.c`. **Regular computation, irregular geometry:**

```
    x(i,j,k) = a·(i,j,k) + A·δ(i,j,k mod M)          ← the replaceable core
```

Connectivity is a fixed 18-point *logical* stencil, so memory layout, neighbour
lookup and parallelism are as cheap as a cubic grid; only the physical vectors
differ per site. Everything downstream (neighbour vectors, weights, volumes) is
derived from `x()` at build time, so making δ dynamical later — driven by local
field energy or regional strain — touches nothing else.

**Weights are solved, not guessed:** per site, `Σⱼ wⱼ dⱼ⊗dⱼ = 2I` exactly
(6 constraints, minimum-norm correction to a 1/|d|² prior). Second-order
consistency *and* isotropy therefore hold at every cell however irregular it is,
so residual anisotropy is genuinely 4th-order fabric structure rather than an
artifact of a lazy weight choice.

### 3.1 Seven moduli, resolved

Displacement amplitude 0.28a. 64 origins × 12 directions.

| fabric | vol sd | vol mx/mn | disp @\|k\|=2 | corrugation | systematic | incoherent |
|---|---|---|---|---|---|---|
| M=2 | **4.75%** | 1.22 | 18.53% | 1.20×10⁻³ | **31.4%** | 24.0% |
| M=3 | 14.44% | 1.83 | 10.55% | 3.26×10⁻³ | 36.8% | 32.5% |
| M=5 | 13.69% | 2.27 | 10.37% | 6.60×10⁻³ | 42.7% | 38.0% |
| M=7 | 16.58% | 3.31 | 10.00% | **9.10×10⁻³** | 48.4% | 49.4% |
| M=11 | 16.58% | 3.31 | 9.85% | 8.32×10⁻³ | 46.3% | 47.9% |
| M=13 | 16.58% | 3.31 | 9.63% | 7.84×10⁻³ | 37.0% | 49.5% |
| hash | 16.58% | 3.31 | **9.63%** | 7.43×10⁻³ | 44.6% | 50.4% |

**Cell volumes are regular**, as required: sd 4.75% → 16.6%, max/min 1.22 →
3.31. Same volume, different edges. The spread is set by the displacement
amplitude and is a free parameter.

**Coherent vs incoherent is the distinction that matters.** A crystal makes
every particle everywhere feel the *same* preferred direction — a broken
rotational symmetry, and fatal. A disordered fabric makes pinning random from
place to place — noise, statistically isotropic, harmless. Averaging the
corrugation over origins *first* and taking the directional spread of those
averages isolates the first; the origin-to-origin scatter is the second. An
earlier pass with 8 origins reported 78–119% "systematic" anisotropy that was
**pure sampling noise** — it collapsed to 37–48% at 64 origins. With incoherent
scatter ~50% and 64 origins the noise floor on a 12-direction range is ≈21%, so
the measured 31–48% sits about 2× above it: there is a real systematic
component, of order 30–40%, but it is not sharply resolved.

### 3.2 The prediction failed

Stated before the runs: *pinning anisotropy falls with M, because a
superlattice of period M splits each Bragg peak into ~M³ satellites.* It does
not. Systematic anisotropy **rises** 31% → 48% from M=2 to M=7, then plateaus.
What does fall monotonically with M is the **dispersion** anisotropy
(18.5% → 9.6%), and what rises is the **pinning strength** (1.2×10⁻³ → 9.1×10⁻³,
7.6×) and the incoherent scatter (24% → 50%).

So irregularity buys *localisation strength* and *dispersion isotropy*, and
costs *randomness*. It does not buy pinning isotropy.

---

## 4. Where this leaves the design

| | dispersion aniso | pinning aniso | periodicity / momentum / memory |
|---|---|---|---|
| SC 6 | 22.7% | ~52% | all preserved |
| **Kelvin BCC 8+6 tuned** | **1.5%** | **~33%** | all preserved |
| procedural M=13 / hash | 9.6% | ~37–45% | preserved (logical lattice) |
| amorphous point set | ≲5% (noise-limited) | incoherent | **all lost** |

**On the measurements taken, the tuned Kelvin crystal is ahead on both
channels** — 6× better dispersion isotropy than the procedural fabric and
comparable-or-better pinning anisotropy — while keeping exact periodicity,
discrete translation symmetry and regular memory.

**Caveat that blocks a firm ranking:** §2's Kelvin pinning number and §3's
procedural pinning number come from **different statistics** (peak-to-peak
barrier over one lattice period from a single origin, vs origin-averaged RMS
corrugation over 12 directions). They are not strictly comparable. The
`fabric_proc` framework at displacement amplitude 0 reduces exactly to simple
cubic and so provides an internal control in the *same* statistic; that run is
what settles the ranking, and until it is in, "Kelvin wins" is **provisional**.

### Open, and genuinely undecided

* Is ~30–40% systematic pinning anisotropy disqualifying, or tolerable at this
  stage? Under the discrete-fabric programme pinning is the *desired* effect
  (it localises); its anisotropy means orientation-dependent mobility.
* Displacement amplitude is unexplored — all seven runs used 0.28a. It trades
  volume regularity against pinning strength and is the obvious next knob.
* A quasicrystalline (icosahedral) fabric is the only construction that is both
  deterministic and 5th-order isotropic (0.68% dispersion spread, §2). It is not
  periodic, but it is not random either — untested here and arguably the
  strongest remaining candidate.
