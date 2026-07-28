# v86 — running down the NEXT_PROGRAM list
**Date:** 2026-07-26 · Results for the rungs `NEXT_PROGRAM.md` left ready to run.
Instruments: `n_battery/hc3_volume.py`, `n_battery/hc1_gauged.py`,
`n_battery/d8b_symbolic.py`, `n_battery/hc6_profiles.py`, `n_battery/hc1_bdg.py`.

Tags: **[D]** derived · **[M]** measured · **[P]** postulated · **[C]** conjecture

---

## 0. Headline — the list's own gate flipped

`NEXT_PROGRAM.md` opened HC-6 with a live possibility:

> It is therefore a live possibility that **the flavoured branch has no
> GSS-unstable region at all**, in which case HC-6 is unrunnable — and that is
> itself a result worth stating.

**It is not the case. HC-6 has targets, and they were hiding in plain sight.**
HC-3's `n(D) = 1` everywhere was an artifact of scanning a single detuning ray
at a single base frequency (ω̄ = 1.42), far below the branch's VK turning
point. The volume scan finds `n(D) = 0` on a whole region near the thin-wall
end, and **two independent code paths** put the turning point at the same
place.

| rung | status | one-line result |
|---|---|---|
| **HC-3-volume** | **done** | n(D) = 0 on 20 window-valid partitions near ω ≈ 1.48 (11 interior, 9 near-edge) → **HC-6 is runnable**; the cause is the ordinary VK turn, not a flavour effect |
| **HC-1-gauged** | **done, scoped** | gauged operators built; n(H_ω)**^(ℓ=0)** = 1 on the branch; Coulomb term is PSD so gauging cannot add a negative direction; **gauged VK turn measured at ω ≈ 1.485** |
| **D8b** | **done (derivation)** | exact link-centred discrete flux derived, verified to 2×10⁻¹⁶; pinning residual 1.2×10⁻⁸ ⇒ the 4.4% defect is **not** pinning. Link flux **not yet** re-integrated on N7 data |
| **HC-1 completion debt** | **done** | catalog extended 1.47 → 1.495, all 7 ω accepted, n(H_ω) = 1, zero bound modes |
| **HC-6** | **running** | target/control pair seeded and integrating |
| HC-4-lite, N8, N10-spectral | queued | CPU-bound behind HC-6 |

---

## 1. HC-3-volume — the gate, and it opens  **[M]**

### 1.1 What the old scan actually covered

Fixed total charge means moving in the traceless plane of (ω₀,ω₁,ω₂). Every
traceless detuning is

```
    delta_a(theta, rho) = rho * cos(theta - 120° * a),    a = 0,1,2
```

which sums to zero identically. The S₃ permutation action is rotation by 120°
plus a reflection, so the **fundamental domain is θ ∈ [0°, 60°]**. The old HC-3
ray, δ ∝ (−1, +½, +½), is θ = 180° ≡ θ = 60° — **one edge of that domain.** The
opposite edge θ = 0° (one flavour detuned *up*, two *down*) had never been
sampled, and it is physically a different object.

### 1.2 The scan

424 converged points over ω̄ ∈ {1.39, 1.42, 1.45, 1.48} × θ ∈ [0°,60°] × ρ ≤ 0.10.
(ω̄ = 1.36 failed at the symmetric seed and is not reported.)

**Window filter, applied before any claim.** A converged Newton solution is not
automatically physical: the continuation can walk a component past the window
edge ω = 1.5, where Q jumps by an order of magnitude and one eigenvalue blows
up to ~10⁶. 23 points were dropped on `ω_a ∈ (1.30870, 1.5)` — **including the
only n(D) = 2 point**, which was a continuation artifact and is *not* claimed
as a target.

| | window-valid | dropped |
|---|---|---|
| n(D) = 0 | **20** | 1 |
| n(D) = 1 | 381 | 21 |
| n(D) = 2 | 0 | 1 |

By base frequency:

| ω̄ | points | n(D) seen |
|---|---|---|
| 1.390 | 134 | {1} |
| 1.420 | 133 | {1} |
| 1.450 | 94 | {0, 1} |
| 1.480 | 40 | {0, 1} |

So the old ray's region (ω̄ = 1.42) really is uniformly n(D) = 1 — HC-3 was not
wrong, it was *local*. The signature change lives at higher ω̄. **[M]**

### 1.3 Why: the VK turning point, confirmed twice

The symmetric branch's total charge has a **minimum**, i.e. dQ/dω changes sign:

| ω | Q (from `n4_hc3_flavored`) | Q (from `hc1_bdg`, independent) |
|---|---|---|
| 1.470 | 91.93 | 91.9 |
| 1.475 | 88.53 | 88.5 |
| 1.480 | 86.73 | 86.7 |
| 1.483 | — | **86.7** |
| 1.485 | 87.27 | 87.3 |
| 1.490 | 92.34 | 92.3 |
| 1.495 | 111.35 | 111.3 |

dQ/dω runs −118 at ω = 1.480 to +531 at ω = 1.485, so the turning point is at
**ω ≈ 1.481–1.483**. Two solvers with different discretisations and different
code paths agree to 0.1%. This is the classical Van der Waals / Vakhitov–
Kolokolov turning point of the Q-ball branch, and it is the reason n(D) drops
from 1 to 0: past it, dQ/dω > 0 in every eigen-direction. **[M]**

### 1.4 The targets

20 window-valid partitions, all at ω̄ = 1.48 except one at ω̄ = 1.45. Smallest
|eigenvalue| among them is 4.3 against a matrix scale of ~1600, so these are
genuine sign changes, not near-zero flips on finite-difference noise. Overall
`max|D−Dᵀ|/max|D| = 2.0×10⁻³`.

**The GSS reading.** The criterion is n(H_ω) = n(D). The right anchor here is
the **ungauged** HC-1 (`hc1_bdg`, which sweeps ℓ, now extended to ω = 1.495 —
§4), not the gauged run of §2, which is ℓ = 0 only and a different sector.
Ungauged n(H_ω) = 1; these partitions carry n(D) = 0 ≠ 1; indices mismatch;
decay candidates. The inference direction is correct (checked against
`GROUNDING.md` §1 by the reviewer, Finding 2.4: not inverted).

**Three things this is NOT, all raised by the grok-4.5 review (Finding 2.2,
MAJOR — "turning real; interpretation overstated"):**

1. **It is not a partition-space discovery.** At ω̄ = 1.48 the *symmetric*
   point still has n(D) = 1; the 19 detuned points cross zero because the
   symmetric eigenvalue is already only O(10)–O(50) from the turn. The detuned
   partitions **inherit** the crossing, they do not create it.
2. **It is not novel flavour physics.** Past the turn, n(D) = 0 with n(H) = 1
   is the ordinary monochromatic VK-unstable branch — textbook Q-ball
   behaviour. Framing all 20 as a "partition-volume index boundary" oversells
   flavoured novelty.
3. **Not all 20 are clean.** **9 of the 20** have some ω_a > 1.495, and several
   show soft pathology (large FD eigenvalues, Q jumping as ρ increases). The
   hard window cut at 1.5 is *necessary but not sufficient*; a residual/edge
   quality diagnostic is needed before any of those are seeded.

What survives as new: the crossing exists in this project's flavoured sector,
it was invisible to the single-ray scan, and the **interior** (non-near-edge)
n(D) = 0 points give HC-6 something legitimate to seed.

---

## 2. HC-1-gauged — the census's biggest gap, closed at ℓ = 0  **[D,M]**

Full derivation in the file header. The three structural results:

### 2.1 The Goldstone survives exactly — and validates the operator

The background equation *is* `(−∇² + m² − wt(r)² + P₀) f = 0`, so `L₀ f = 0`
identically even with `wt(r) = ω − χ(r)` radial. Measured residual
3.18×10⁻⁵ relative, against a lowest |eigenvalue| of L_x^sym of 1.15×10⁻¹ —
**3623× above the discretisation floor**, so the index is resolved. This is a
sharp check that the gauged operator was constructed correctly, not merely the
ungauged limit re-labelled. **[D,M]**

### 2.2 The flavour argument now holds AT g = 0.05, not by inheritance

`L_x^flav = L₀ − A` with `A = 2P₀` and `P₀ = μf⁴/(1+κf⁶)² < 0` for μ < 0, so
`L_x^flav ≥ L₀ ≥ 0`. **The argument never used g.** Measured minimum of
L_x^flav over the scan: +3.48×10⁻². The census statement "flavour channels
contribute no negative directions" upgrades from *heuristic at g = 0.05* to
*established at g = 0.05*. **[D,M]**

### 2.3 The Coulomb back-reaction is positive semi-definite

Eliminating the constrained A₀ perturbation gives

```
    L_x^sym  ->  L_x^sym + 12 g² (wt f) K⁻¹ (wt f),   K = -lap + 3 g² f²
```

K is positive definite (screening inside the ball), so the added operator is
PSD: **gauging can only raise eigenvalues, never add a negative direction.**

| ω | L_x^sym min | with Coulomb | shift |
|---|---|---|---|
| 1.41 | −0.11517 | −0.04758 | **+0.06759** |
| 1.42 | −0.15530 | −0.09926 | +0.05604 |
| 1.45 | −0.28287 | −0.24273 | +0.04014 |
| 1.48 | −0.46840 | −0.43703 | +0.03136 |
| 1.49 | −0.56617 | −0.53740 | +0.02877 |

Not a small effect: at ω = 1.41 the Coulomb term cancels **59%** of the
negative eigenvalue, and the fraction grows as the branch edge is approached.

### 2.4 The two window edges are NOT the same failure mode  **(corrected)**

> **This section originally claimed the narrowed window is "the background
> ceasing to exist, NOT a GSS index change". That was half wrong, and the
> grok-4.5 review (Finding 1.4, MAJOR) caught it. The corrected version, with
> the missing measurement now taken, is below.**

The gauged window has **two** edges and they fail differently.

**Low-ω edge — existence.** The continuation finds no gauged solution at all
for ω ≤ 1.405, reproducing the v69 window bottom 1.406 from an independent
path. Q_max = 921 lives here. On every surviving point n(H_ω)^(ℓ=0) = 1. This
edge genuinely is an existence limit, not an index change. **[M]**

**High-ω edge — a real VK turn, i.e. a genuine GSS index change.** The gauged
branch does *not* stop; it turns. Measured directly (`hc1_gauged_vk.log`):

| ω | Q (gauged, g=0.05) | dQ/dω | n(D) |
|---|---|---|---|
| 1.4800 | 90.9473 | −477.2 | 1 |
| 1.4820 | 90.1497 | −302.7 | 1 |
| 1.4840 | **89.7363** | −90.5 | 1 |
| 1.4860 | 89.7878 | +173.1 | **0** |
| 1.4900 | 91.8590 | +997.3 | 0 |
| 1.4950 | 101.9092 | +2497.1 | 0 |

**Gauged Q_min = 89.7363 at ω ≈ 1.484–1.486.** The independent grok re-solve
put it at Q ≈ 89.70, ω ≈ 1.485 — agreement to 0.04%. (Ungauged turn, for
contrast: Q ≈ 86.59 at ω ≈ 1.482.)

So past ω ≈ 1.485 the gauged branch has n(H_ω)^(ℓ=0) = 1 while n(D) = 0 — the
indices do **not** match, and that *is* a GSS mismatch, in the gauged case as
well as the ungauged one. The corrected statement:

> The **low-ω** termination and Q_max = 921 are existence limits with
> n(H_ω)^(ℓ=0) = 1. The **high-ω** end is a Vakhitov–Kolokolov turning point
> where the GSS change is carried by n(∂Q/∂ω), which the original hc1_gauged
> run did not compute. Both edges are real; only one is an existence
> statement. **[M]**

**Scope, and it bounds every n(H) claim above:** this is the **ℓ = 0 sector
only** (`base = M2 - wt^2 + P0`). The index is therefore
**n(H_ω)^(ℓ=0) = 1**, not the full n(H_ω). ℓ = 1 should be the translational
Goldstone; ℓ ≥ 2 can host multipole/fission negatives, and precisely at the
large-Q end (ω = 1.41, Q ≈ 529 rising to Q_max = 921) is where one would
expect them. Closing this is a bounded extension — add the ℓ(ℓ+1)/r²
centrifugal term and re-scan — not a new build, but it is **not done**.

**Also corrected (Finding 1.7):** the file's docstring calls this a "BdG
spectrum". It is not — it diagonalises the static operators L₀, L_x^flav,
L_x^sym + Coulomb, which is the right object for the GSS **Morse index** but
not the dynamical BdG frequency spectrum.

---

## 3. D8b — the discrete momentum balance  **[D,M]**

7/7 symbolic identities verified, plus two numeric arms.

### 3.1 The gradient sector has an exact LINK-centred flux

With `G^d_n = f_{n+d} − f_n`:

```
   d = i :  (lap_i f)(Dc_i f) = Δ_i⁻[ (G^i_n)² / (2 dx³) ]
   d ≠ i :  (lap_d f)(Dc_i f) = Δ_d⁻[ G^d_n B^i_n / (2 dx³) ]
                              − Δ_i⁻[ G^d_{n−d} G^d_{n−d+i} / (2 dx³) ]
```

Verified site-by-site symbolically (sympy, 3D 4³ periodic) and at scale on a
random N=24 field — **max site-by-site error 2.09×10⁻¹⁶, total sum 9.13×10⁻¹⁹**.
A random field exercises every lattice wavenumber, so this is a far stronger
test than a smooth soliton. The flux lives on **links**; sampling the continuum
T^{ji} at cell centres is the error. **[D,M]**

### 3.2 The mass sector contributes exactly zero

`Σ_n m² f_n (f_{n+i} − f_{n−i}) = 0` by relabelling, for any lattice, any dx.
No flux term is needed for it at all. **[D]**

### 3.3 The only obstruction is lattice pinning, and it is tiny

The nonlinear potential does **not** telescope:

```
   R_i = −1/(2dx) Σ_n [ U'(f_n) f_{n+i} − U'(f_{n+i}) f_n ]
```

which vanishes iff U' is linear. This is the Peierls–Nabarro pinning force.
Measured by sliding a real ω=1.430 Q-ball through one lattice cell:

| dx | max\|R_x\| | E_grad | **R/E** |
|---|---|---|---|
| 0.60 | 7.87×10⁻⁴ | 18.40 | 4.27×10⁻⁵ |
| 0.40 | 3.53×10⁻⁶ | 18.51 | 1.91×10⁻⁷ |
| 0.30 | 2.23×10⁻⁷ | 18.54 | **1.20×10⁻⁸** |
| 0.24 | 1.66×10⁻⁷ | 18.56 | 8.94×10⁻⁹ |
| 0.20 | 3.77×10⁻⁸ | 18.57 | 2.03×10⁻⁹ |

Steeply falling, as an exponentially-small-in-width/dx pinning force must.

**What is established:** across the whole scanned dx range the pinning term is
6+ orders below the 4.4% flux defect, so **pinning is ruled out as the source
of that defect**. The ordering is robust; "at production dx" is the loosest
part of the wording, since N7's defect was measured on a specific mesh. **[M]**

**What is NOT established (grok Finding 3.5, MAJOR).** `sfa_momentum.c` still
samples the continuum T^{ji} at cell centres with central differences — it has
**not** been rewritten to the link flux, and the link flux has **not** been
re-integrated on the N7 run to show the residual collapse. So:

> "The defect is the surface object and §3.1 gives the corrected one" is a
> **derivation claim, not a closed re-measurement.** Pinning is excluded; that
> the link flux restores the balance on N7 data is **untested**. Closing it
> means porting §3.1 into `sfa_momentum.c` and re-running N7's plane scan.

### 3.4 The structural statement that should go in the theory document

There is **no exactly conserved lattice momentum** here, and the reason is not
sloppiness. Charge conservation is exact by construction because its U(1) is a
*continuous* internal symmetry that survives discretisation intact. Translation
is broken by the lattice to a *discrete* subgroup, and discrete symmetries carry
no Noether current. So the hierarchy D8a exposed —

| quantity | conservation | why |
|---|---|---|
| Q | 1.3×10⁻⁹ | continuous internal symmetry, exact by construction |
| Gauss residual | 10⁻¹³ | constraint, exact by construction |
| E, P | ~10⁻⁴ | integrator quality only, no structural guarantee |

— is **structural, not a bug**, and it bounds what any re-derivation can buy.
**[D]**

---

## 4. HC-1 completion debt — closed at the thin-wall end  **[M]**

The catalog stopped at ω = 1.47. Extended to 1.47, 1.475, 1.48, 1.483, 1.485,
1.49, 1.495: **all 7 accepted** (Goldstone check passed at every one, none
rejected), n(H_ω) = 1 throughout, and **zero bound internal modes** — consistent
with HC-1's finding that none exist for ω ≥ 1.36. Log: `n_battery/hc1_thinwall.log`.

This matters more than bookkeeping now: the newly interesting region (§1.3) is
exactly the range the catalog did not cover, and n(H_ω) = 1 there is what makes
the n(D) = 0 points a genuine index *mismatch*.

---

## 5. HC-6 — now runnable, and running

Target/control pair chosen at matched charge and energy, which is the whole
point of a converse test:

| seed | ω = (ω₀,ω₁,ω₂) | Q_tot | E | eig(dQ_a/dω_b) | n(D) | role |
|---|---|---|---|---|---|---|
| `target` | (1.4950, 1.4950, 1.4500) | 116.44 | 177.40 | [+61, +3605, +5464] | **0** | GSS mismatch |
| `control` | (1.4650, 1.4650, 1.4200) | 117.95 | 178.70 | [−545, +347, +773] | 1 | GSS match |

1.3% apart in Q, 0.7% in E. Both integrating at N=128, L=22, ungauged, absorbing
BC, T=300, with per-component charges `Q_p0/Q_p1/Q_p2` in the diag stream.

**Pre-registered reading, stated before the runs finish:**

* target redistributes charge between components while control does not →
  **GSS converse confirmed**
* both persist → n(D) = 0 is not sufficient for decay on the timescale probed;
  report the bound, do **not** claim confirmation
* both decay → the seeding or the box is the instability, not the partition;
  the test is void and must be redesigned

### 5.1 The pair is CONFOUNDED — a third arm was required  (grok Finding B1, MAJOR)

Reviewed **before** reading the result, which is the only time it counts. The
target and control differ in three ways that all bias the target to look worse:

| | target | control |
|---|---|---|
| mean ω_eff | **1.480** | 1.450 |
| vs the VK turn (≈1.482) | **on/past it** | below it |
| tail extent r(f<0.01) | **15.0** | 9.5 |
| margin to the sponge (L−damp = 19) | **~4** | ~9.5 |

So "target decays, control does not" would establish nothing: past the turn even
a *monochromatic* ball is classically VK-unstable, and the target's tail sits
much closer to the absorber.

**The discriminator, now running:** a **monochromatic** seed past the same turn,
matched to the target in charge and extent — ω = (1.4955, 1.4955, 1.4955),
Q = 115.43 vs 116.44 (0.9%), E = 175.64 vs 177.40 (1.0%), extent 14.64 vs 15.0
(2.5%) — carrying **no flavour structure at all**. Then:

* target ≈ mono, both ≠ control → the effect is **VK + box, not flavoured GSS**
* target ≠ both → flavour structure is doing something real
* nothing moves → report a bound only

### 5.2 Partial reading at t ≈ 101 (target/control only — NOT a verdict)

| run | dQ_tot/Q | max\|Δx_a\| |
|---|---|---|
| target (n(D)=0) | **−1.328%** | 1.59×10⁻³ |
| control (n(D)=1) | −0.005% | 9.25×10⁻⁶ |

The target loses 260× more charge and its partition fractions move 170× more —
but the fraction drift is only 0.16% while the charge loss is 1.33%, i.e. the
loss is **mostly uniform**, which is the radiation signature grok predicted and
*not* redistribution to a sector minimum. **No verdict until the mono arm
lands**; `n_battery/hc6_verdict.c` applies the three-way test automatically.

---

## 6. What this changes in NEXT_PROGRAM.md

| item | before | after |
|---|---|---|
| HC-6 | "may be unrunnable — live possibility of no unstable region" | **runnable**, 20 targets, running |
| HC-1-gauged | "the weakest brick in the Part-0 wall" | built; n(H_ω)=1 at g=0.05; Coulomb PSD; l=0 scope flagged |
| D8b | "the target" | derived, verified, and the 4.4% attributed |
| HC-1 debt | ω range incomplete | closed to 1.495 |
| Q_max interpretation | implicitly a stability statement | **an existence statement** |

Still open from the list: **HC-4-lite**, **N8**, **N10 spectral**, **D8c/D9**
(Lean), **EX-2 resolution debt**, and the l ≥ 1 extension of HC-1-gauged that
§2.4 opens.

---

## 8. D8c — closing D8b on real data (grok Finding 3.5)

`n_battery/d8c_balance.py`. Three tests; the first is the deliverable, the third
is the honest limit.

### 8.1 The link flux IS the exact discrete stress — proven on real field data

Pure algebra, no time derivative, no EOM, no integrator. On D8a's actual frames
(N=64, L=16, periodic, ungauged), for every slab boundary X:

```
   sum_slab [ (lap f)(Dc_x f) + pi(Dc_x pi) - m^2 f(Dc_x f) ]  ==  F(X) - F(wrap)
```

**Max relative error 3.3×10⁻¹⁰**, and 0–10⁻¹⁵ wherever cancellation is not
dominant — the residual is the f32 storage floor, not the formula. **[M]**

Two bookkeeping errors had to be fixed to get there, both of which produced
plausible-looking wrong answers first:

1. **Periodic boxes have two boundaries.** The region x < X is bounded by the
   plane *and* by the wrap at −L. Counting only the plane is an error the same
   size as the answer.
2. **The dx power.** The identity gives A_n = F_n − F_{n−1} as a *density*, so
   the per-area flux is F·dx. Missing that produced a constant relative error of
   exactly 1/dx at every X — which is what identified it.

With the factor right, the flux is exactly the continuum stress with two-point
link averaging: ½(∂ₓφ)² − ½(∂_yφ)² − ½(∂_zφ)² + ½φ̇² − ½m²φ², **and no potential
term**. That absence is structural: Vt(s) does not telescope, so on a lattice
its momentum transfer is a *bulk* residual, not a surface flux.

### 8.2 On N7 data the flux helps, but does not close the 4.4%

Plane-to-plane spread of the matter-channel flux over |x| ≤ 2, N7 frame 0:

| object | spread |
|---|---|
| naive cell-centred (what `sfa_momentum.c` integrates) | **16.469%** |
| exact link flux | **11.328%** (1.45× flatter) |
| link flux + bulk PN | 13.549% (1.22× flatter) |

**So D8b's attribution was too strong and is corrected here.** Flux placement
accounts for roughly a third of the plane-to-plane variation, not all of it, and
adding the bulk PN term does not recover the rest.

The likely reason is that N7's premise does not hold on this configuration:
"a divergence-free flux must be flat in the vacuum gap" needs an actual vacuum
gap, and at D = 10 with tails reaching r ≈ 9 the region |x| ≤ 2 is **not**
vacuum. Part of the variation is therefore physical, not a stencil error.

### 8.3 What could not be tested, and why

The *dynamical* closure — dP/dt against the balance — is **untestable on D8a**.
Its snapshots are 1.0 apart while the clock period is 4.4, so the field rotates
1.43 rad between frames and any two-frame average of a quadratic quantity is
meaningless. Closing that needs snap_dt ≲ 0.1 (40+ samples/period), which is a
new run, not a re-analysis.

**Status of Finding 3.5: partially closed.** The flux is now proven exact on
real data (8.1) and its effect on N7's number is measured rather than asserted
(8.2) — but the full attribution of the 4.4% is still open, and the residual is
now diagnosed as a non-vacuum gap rather than a stencil defect.


---

## 9. HC-6 — the observable was void, and the corrected one is positive

### 9.1 The original HC-6 observable cannot work in this kernel  **[D,M]**

HC-6 was designed to watch charge **redistribute between components** — "decay
to the sector minimum". That is impossible here, and not for a subtle reason.

The potential depends on the fields only through `s = Π_a |Φ_a|²`, which is
invariant under **independent** phase rotations `Φ_a → e^{iθ_a}Φ_a`. The kernel's
force term is `−2·Vp(s)·Φ_a·Π_{b≠a}|Φ_b|²` — proportional to Φ_a with a real
coefficient built from moduli alone. So the model carries **U(1)³, not one
U(1)**, and every Q_a is its own conserved Noether charge. (η ≠ 0 would break
it through the curl coupling; these runs are η = 0, where it is exact.)

Measured, in a box with the absorber switched off:

```
   Q_p0 = 50.846808956130  at t = 0
   Q_p0 = 50.846808956130  at t = 201      unchanged to 13 significant figures
   worst relative drift of any Q_a, any arm: 6.97e-13
```

**So partition fractions x_a = Q_a/Q_tot are pinned by a conservation law.**
Every verdict built on them — including the null I first reported, and the
"flavour is doing something real" my own tool first printed — was measuring a
quantity that cannot move. The monochromatic arm's frozen fractions were
likewise pinned, by permutation symmetry, at ~1e-17.

**This is a defect in the rung as specified, not in the run.** HC-6 in
`PROPOSAL.md` / `GROUNDING.md` asks for an observable the model forbids.

### 9.2 The channel that can respond — and it does  **[M]**

A GSS index mismatch does not move charge between U(1)³ sectors; it destabilises
the object **at fixed (Q₀,Q₁,Q₂)**, i.e. in the *shape* channel. With the
absorber off (so absorption cannot masquerade as decay) and E conserved:

| arm | Δ s_max | Δ Q_core | Δ r_core | Δ P_int | ΔE/E |
|---|---|---|---|---|---|
| **target, n(D)=0** | **+43.86%** | **+2.87%** | **−8.89%** | **+30.39%** | 9.6×10⁻⁶ |
| **control, n(D)=1** | +0.17% | +0.002% | −0.36% | +0.54% | −1.7×10⁻⁵ |
| mono, n(D)=0 (sponged) | +134.7% | −3.17% | −21.4% | +45.6% | −5.9×10⁻² |

**Four independent channels move coherently for the n(D)=0 object and are flat
for the n(D)=1 control at matched charge.** Peak density up, core charge up,
core radius down, binding up — the object is *contracting and densifying at
fixed charge*, with energy conserved to 10⁻⁵. That is the classical
Vakhitov–Kolokolov picture: a soliton past the turn migrating toward the stable
branch carrying the same Q. The second n(D)=0 arm (monochromatic, no flavour
structure at all) moves the same way; its Q_core falls instead of rising only
because its sponge is eating the tail.

**The 714× rate ratio I first quoted is retracted** (grok Finding 4, MAJOR).
A log-linear rate is the wrong statistic: the control's fit has R² = 0.001, so
its "rate" is undefined and the ratio was signal/noise, not signal/signal. The
runs also carry heavy breathing (~200 sign changes in ds over the target run),
so a single exponential is not the right model at all. Corrected statistics:

| arm | s_end/s₀ | secular amplitude | breathing (p-p) | R² |
|---|---|---|---|---|
| target n(D)=0 | 1.437 | **+0.4039** | 0.2785 | 0.898 |
| control n(D)=1 | 1.001 | +0.0006 | 0.0372 | 0.001 |

The honest comparison is **amplitude against the control's own noise floor**:
the target's secular change (+0.40) is 11× the control's entire breathing
excursion (0.037), and exceeds its own breathing (0.28). That is a real trend,
but on one run it is not an overwhelming margin.

### 9.3 Review verdict and what is not yet established

Independent grok-4.5 pass (`n_battery/REVIEW_session3_grok45.md`):

| piece | verdict |
|---|---|
| old observable void under U(1)³ | **SOUND** |
| census HC-6 "redistribution" language malformed | **SOUND** (design bug, not a run bug) |
| GSS criterion still right, applies at fixed Q_a | **SOUND** |
| shape channel is the right place to look | **SOUND** |
| **"result flips to positive / GSS confirmed"** | **NOT YET ESTABLISHED — overclaimed** |


* **dx convergence.** dx = 0.3465 gives ~6.6 cells per core half-radius. An
  N = 160 (dx = 0.2767) target run is in flight. If the collapse is
  discretisation rather than physics it should change materially; if physical,
  the four-channel signature should persist.
* **The mono arm ARGUES AGAINST a flavoured reading, not for it** (grok
  Finding 5, MAJOR). It deforms *more* than the flavoured target (3.4 vs 1.44)
  with no flavour structure at all, and is n(D)=0 only because it sits past the
  monochromatic VK turn. So the strongest effect in the table is consistent
  with "past the turn ⇒ shape runs away", which mono and target share and the
  control (below the turn) does not. **Nothing here yet attributes anything to
  flavoured index structure.**
* **The controlled pair is now running.** Two *monochromatic* balls, no sponge,
  matched charge, differing only in which side of the turn they occupy:
  ω = 1.4955 (past, Q = 115.43, n(D)=0) vs ω = 1.4520 (below, Q = 114.52,
  n(D)=1). Same object type, same box, same resolution, no flavour, no
  absorber — if the first collapses and the second does not, VK is demonstrated
  cleanly; only then can the flavoured pair be asked whether flavour adds
  anything on top.
* **Seed-representation confound** (grok Finding 6): the target's lattice E(0)
  sits 0.37% below its continuum value against the control's 0.09%, i.e. the
  target is the less faithfully represented object on this grid.
* The instability **saturates** rather than running away, so what is measured is
  a finite rearrangement, not a blow-up.


---

## 10. HC-6 resolved — VK demonstrated cleanly, flavour excluded

The controlled experiment grok Finding 5 demanded. Two **monochromatic** balls
(no flavour structure at all), no absorber, matched charge, same box, same
resolution — differing **only** in which side of the VK turn they sit:

| arm | ω | Q | n(D) | Δ s_max | Δ r_core | Δ P_int | secular amp | R² | ΔE/E | max ΔQ_a/Q |
|---|---|---|---|---|---|---|---|---|---|---|
| **past turn** | 1.4955 | 115.43 | **0** | **+197%** | **−10.17%** | **+78.6%** | **+1.730** | 0.937 | −3.8×10⁻⁴ | 2.6×10⁻¹³ |
| **below turn** | 1.4520 | 114.52 | 1 | +0.19% | −0.25% | +0.42% | +0.0013 | 0.036 | −1.9×10⁻⁵ | 2.6×10⁻¹³ |

**Every confound is closed here.** Same object type, matched Q to 0.8%, identical
grid and box, charge conserved to 10⁻¹³, energy to 10⁻⁴, no sponge to absorb
anything, no flavour structure to argue about. The only variable is the sign of
dQ/dω. Past the turn the object collapses — peak density triples, core radius
falls 10%, binding rises 79%. Below the turn nothing happens at all.

**The below-turn arm is also the resolution control.** It sits on the same
dx = 0.3465 grid and moves 0.19%, so the discretisation is demonstrably capable
of holding a Q-ball static; a 197% effect on the same grid is not a dx artifact.
An N = 160 confirmation of the past-turn arm is running.

### 10.1 Flavour is excluded, not confirmed

Comparing the two n(D) = 0 objects at the same position relative to the turn:

| n(D)=0 object | secular amplitude |
|---|---|
| monochromatic, past turn | **+1.730** |
| flavoured target | +0.404 |

**The flavoured object deforms ~4× LESS than the plain monochromatic one.** So
flavoured index structure adds nothing to the instability — if anything it damps
it. HC-3-volume's 20 "flavoured HC-6 targets" are ordinary VK points that happen
to carry a flavour label, exactly as the review said (Finding 2.2, Finding 5).

### 10.2 Final status of HC-6

| claim | status |
|---|---|
| original observable (charge redistribution) | **VOID** — forbidden by U(1)³, verified to 13 digits |
| GSS/VK instability exists in the fabric at fixed Q | **DEMONSTRATED** [M] — controlled monochromatic pair |
| the instability is *flavoured* | **EXCLUDED** [M] — flavoured object deforms 4× less |
| "714× rate ratio" | **RETRACTED** — wrong statistic (control R² = 0.001) |

What HC-6 delivers, restated so it is not oversold: the fabric's Q-balls obey the
Vakhitov–Kolokolov criterion, demonstrated with all confounds removed. That is a
real and load-bearing check on the soliton sector — but it is textbook Q-ball
physics recovered, not new flavour physics, and the rung's original flavoured
framing does not survive.
