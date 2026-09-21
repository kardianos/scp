# The quantized harmonic graph — one law for emission, absorption, shells and charge

**Status:** symbolic + Go only. No kernel change, no battery run, no adoption.
Written 2026-08-08 under `v93/theory/qgraph/` against the brief
`qgraph/BRIEF_claude.md`.

**Companions:** `MAXIMA_NOTES.md` (what was proved and how to re-run it),
`SIM_PREDICTIONS.md` (falsifiable bars), `VERDICT.md` (queue / do-not-queue).
Every number quoted below is reproduced by `go run ./cmd/qgdemo` and gated by
`go test ./...` (35 assertions, green).

**Ontology held fixed:** energy only changes mode; no background; no imported
U(1) field and no electron species; free-cell is the substrate, this is its
symbolic shadow; nothing before v89 is consulted.

---

## 0. The claim in one paragraph

A live contact set is a finite graph. Its Laplacian has harmonics; the
harmonics fall into **degeneracy classes forced by the graph's own symmetry**,
and those classes are shells with **integer capacities**. Quanta occupy them
under a cap. **Emission and absorption are one operation** — a quantum moves
between two harmonics and the difference is paid to or from the field — so the
quantum count `N` is conserved and only the field is created and destroyed.
**Charge is the deficit** `Q = Z − N` between the graph's bound capacity and
the quanta present; it is conserved exactly because `N` is conserved and `Z` is
spectral. **Circulation is the same defect seen dynamically**: a full shell's
counter-rotating members cancel, an unclosed shell's do not. **Stefan's fourth
power is the field sector's spectral dimension plus one**, so `p_rad = 4` is
the assertion `d_s = 3` and nothing else. **Monodispersity is the shell gap
over the temperature**, never the absence of radiation.

---

## 1. Objects

Let `G = (V, E)` be the contact graph *as it stands now*. It is a chart, not a
stage: nothing below assumes `V` survives to the next instant, and every
quantity used is a spectral or topological invariant of `G` at that instant.

| object | definition | code |
|---|---|---|
| harmonics | eigenpairs of `L = D − W`, `omega_a = sqrt(lambda_a)` | `Spectrum` |
| action atom | `eps_a = A0 * omega_a / (2*pi)` (the programme's door quantum) | `Mode.Eps` |
| occupation | `n_a` in `0..Nu`, `Nu = 1` (cap refusal / PAULI-0) | `Occ` |
| shell | a degeneracy class of `L`; capacity `c_k = dim(eigenspace) * Nu` | `Shells` |
| well | the shells below the binding pitch (the `Es` depression depth) | `Well` |
| bound capacity | `Z = sum_{k bound} c_k` | `Well.Z` |
| cycle rank | `b1 = |E| − |V| + comps`: where winding can live | `CycleRank` |

`Nu = 1` is not fermion statistics imported; it is the cap that `laws_V3a`
already carries, written at mode grain.

---

## 2. The law table

```
L1  HARMONICS   the modes of a contact set are the eigenvectors of its own
                Laplacian; nothing else is available without extra structure.

L2  ATOMS       occupation is integer; a mode carries at most Nu quanta;
                a quantum of mode a carries action eps_a = A0 omega_a / 2 pi.

L3  ONE DOOR    the only matter operation is a PROMOTION a -> b.  It pays
                (eps_a - eps_b) to the field if positive, and takes it from
                the field if negative.  Emission and absorption are this one
                operation read in its two directions.  N = sum n_a is
                therefore conserved by every matter event.

L4  ONE RATE    a single coupling A_ab drives both directions:
                  O(a->b) = A n_a (Nu - n_b) (1 + N_w)
                  I(b->a) = A n_b (Nu - n_a)  N_w
                with w = eps_a - eps_b.  No second constant, no second law.

L5  SELECTION   the door is a vertex profile d_v; its matrix element is
                  M_ab = sum_v conj(phi_a) d_v phi_b,
                i.e. the transform of d at the winding difference.  The door
                profile IS the selection rule.
```

Charge, shells, temperature and Stefan's law are **consequences** of L1–L5,
derived in §4–§8. No further postulate is added.

---

## 3. Emission and absorption are one law

**Theorem 3.1 (exact duality).** With the single coupling `A` of L4, Fermi
occupancy on matter and Planck occupancy on the field *at one temperature*,

```
O(a->b) − I(b->a) ≡ 0.
```

*Proof:* `detailed_balance.mac` §CLAIM 1 returns identically `0`; the
underlying identity is

```
[n_a (1−n_b)] / [n_b (1−n_a)]  =  N_w / (1 + N_w)  =  exp(−w/T),
```

both sides checked to `0` symbolically. The Boltzmann factor is **forced** by
the shape of L4, not assumed. `TestDetailedBalanceIsExact` reproduces this to
`1e−12` across a grid of `(T, mu, eps)`; `TestDetailedBalanceBreaksOffEquilibrium`
confirms the identity is not vacuous — a hotter field drives net absorption. ∎

**Corollary 3.2 (a cold closed object is dark).** In an aufbau ground state
every downward promotion is blocked by cap refusal at the destination. The
object cannot emit at all until something excites it. (`TestEmissionAndAbsorptionAreOneOperation`.)
Darkness is derived, not switched off.

### 3.3 The coarse form, and how it maps onto `laws_V3a`

Coarse-graining L4 over a saturated stack gives, in the adopted notation,

```
O = k_rad * cap * x^p                                 (adopted radiance)
I = k_rad * cap * ( xb^p  +  g * x^(p−1) * xb )       (its dual — SAME k, p)
```

The first intake term is the spontaneous dual of radiance. The second is the
**stimulated cross term**, nonzero only when the incoming field is phase-matched
to what the object already carries; `g` is the coherent fraction.

**Theorem 3.4 (the box-fill theorem).** With `g = 0` the balance root is
`x* = xb`, exactly, for every `k`, `cap`, `p`.

*Proof:* `Balance = k cap (xb^p − x^p)`, whose only positive root is `x = xb`
(`detailed_balance.mac` §CLAIM 2; `TestRandomBathFixedPointIsTheBathItself`). ∎

This is the derivation of a measured null. A random-phase bath can never hold
matter above ambient — the best it can do is fill the box to the bath's own
load. `MECHANISM.md` records exactly that for S07–S15 and S08/S12/S19:
*"open I fills volume"*, *"I ≈ O at filled capacity ≠ particle"*. Under this
law that outcome was not a tuning failure; it was the theorem.

**Theorem 3.5 (coherence lifts the fixed point).** With `g > 0` the balance
condition in the load ratio `u = x/xb` is `u^4 − g u^3 − 1 = 0`, giving

| `g` | 0 | 0.5 | 1 | 2 | 3 |
|---|---|---|---|---|---|
| `u* = x*/xb` | 1.000 | 1.153 | 1.380 | 2.107 | 3.036 |
| `Phi*` | 0.018 | 0.033 | 0.067 | 0.364 | 1.569 |

(`k=0.05, cap=2.5, xb=0.62, p=4`; `TestCoherentBathHoldsMatterAboveAmbient`.)
In every case `d(I−O)/dx < 0`: the fixed point is restoring. And `Phi* > 0`
always — **the balanced object radiates**. `I = O` is a live cycle, never a
sealed bag. That is `METABOLIC_MATTER.md` §2.3 made exact.

**What this says to the programme:** the missing intake is not *more* energy.
It is *phase-matched* energy. `g`, not `noise_amp`, is the knob that separates
a particle from a filled box.

---

## 4. Stefan's fourth power is a dimension, not an axiom

**Theorem 4.1.** For a field sector with density of states `g(w) ~ w^(d_s−1)`
and Planck occupancy,

```
u(T) = Gamma(d_s+1) zeta(d_s+1) * T^(d_s+1),     j = (c/4) u  ~  T^(d_s+1).
```

Hence `p_rad = d_s + 1`.

*Proof:* substituting `w = T y` in `g(w) w nbar(w) dw` gives
`T^(d_s+1) y^d_s/(e^y − 1) dy` — the `T` power is exact and geometry-free —
and the Bose sum evaluates the constant (`stefan.mac` §1–§3). For `d_s = 3` the
constant is `pi^4/15 = 6.493939`, recovered numerically to `1e−9`
(`TestStefanConstantIsPi4Over15`). ∎

**Therefore `p_rad = 4` is the claim `d_s = 3` and nothing else** — and on a
graph `d_s` is a *measurement*. Two estimators are implemented and calibrated
against exact tori:

| torus | `d_s` (random walk) | `d_s` (eigenvalue count) | implied `p_rad` |
|---|---|---|---|
| 1D, L=60 | 0.994 | 1.211 | 2 |
| 2D, L=30 | 2.001 | 2.165 | 3 |
| 3D, L=14 | 3.020 | 3.310 | 4 |

`ReturnDimension` (lazy-walk return probability, `p(t) ~ t^(−d_s/2)`) is
accurate to ~0.02 and needs no diagonalisation — it is the estimator that can
run on a real free-cell contact graph of any size.
`SpectralDimensionFromEigenvalues` carries a stable `+0.25` finite-size bias
and is documented as such; it must not be quoted to better than `±0.3`.

**This is the sharpest test in the note.** A free-cell foam is not obviously a
3D continuum: near jamming, contact networks are famously anomalous. If the
measured `d_s` of the freecell contact graph is not 3, then `p_rad = 4` is
mis-set and the law's own outtake exponent is wrong. See `SIM_PREDICTIONS.md`
P1.

### 4.2 The `x ↔ T` dictionary, and where it fails

The law writes `O ~ x^4`, and §4.1 gives `O ~ T^4`. The two agree **only if
`x ~ T`**, which requires `Em` linear in `T` — a stack with a *fixed, fully
live* mode count (Dulong–Petit). For a small object with a cap that is the
normal regime, so the adopted law is self-consistent at high load.

At low load, modes freeze out. If the object's own stack has spectral dimension
`dm`, then `Em ~ T^(dm+1)`, so `T ~ x^(1/(dm+1))` and

```
O ~ T^4 ~ x^(4/(dm+1)):     dm=1 -> p_eff = 2,  dm=2 -> 1.33,  dm=3 -> 1.
```

**The adopted `p_rad = 4` is the saturated limit only.** At small `x` the true
exponent is much softer, which means the current table *under-radiates small
hoards by a large factor*. Since radiance-as-selection works by eroding a
hierarchy from above, this directly changes the Face C selection story.
`LoadIsLinearInT` is the in-code regime check; `TestLoadIsLinearInTOnlyWhenSaturated`
exhibits the freeze-out signature. See `SIM_PREDICTIONS.md` P3.

---

## 5. Shells are graph symmetry

A shell is a **degeneracy class of the Laplacian**. Its capacity is the
dimension of that eigenspace, an integer nobody chose.

```
C6   ring        b1= 1  caps=[1 2 2 1]      magic=[1 3 5 6]
C10  ring        b1= 1  caps=[1 2 2 2 2 1]  magic=[1 3 5 7 9 10]
Q3   cube        b1= 5  caps=[1 3 3 1]      magic=[1 4 7 8]
Q4   tesseract   b1=17  caps=[1 4 6 4 1]    magic=[1 5 11 15 16]
icosahedron      b1=19  caps=[1 3 5 3]      magic=[1 4 9 12]
P8   path (tree) b1= 0  caps=[1 1 1 1 1 1 1 1]
```

Two of these deserve to be read slowly.

**The icosahedron gives capacities 1, 3, 5.** Those are `2l+1` for
`l = 0, 1, 2`: an s/p/d ladder. Its magic numbers are therefore **1, 4, 9** —
exact squares, the hydrogenic `n²` filling — produced by the automorphism
group of a 12-vertex contact graph, with no orbital, no central potential and
no angular momentum operator imported anywhere.
(`TestIcosahedronShellLadder`, `TestRingMagicNumbers`.)

**The hypercube gives binomial capacities `C(d,k)` at eigenvalue `2k`,**
verified for `d = 2..5` (`TestHypercubeShellCapacitiesAreBinomial`).

**A tree gives capacity 1 everywhere.** Its spectrum is non-degenerate, so
there is no multiplet to half-fill — a fact §6 turns into a theorem.

---

## 6. Charge

### 6.1 Definition

```
Z = bound capacity of the well  (spectral invariant of the contact graph)
N = quanta present              (conserved by L3)
Q = Z − N
```

Worked example — the icosahedral well bound at shell 2, `Z = 9`:

| `N` | 6 | 7 | 8 | 9 | 10 | 11 |
|---|---|---|---|---|---|---|
| `Q` | +3 | +2 | +1 | **0** | −1 | −2 |

`Q = 0` exactly at closure; both signs occur (a quantum held above closure is
the negative branch). (`TestChargeIsTheShellDeficit`.)

**This is not renamed `Em`.** `Z` counts eigenspace dimensions — a symmetry
fact about the contact graph, independent of how much energy is present — while
`N` counts atoms. Two objects with the same `Em` and different graphs have
different `Q`; two objects with the same graph and different `Em` at fixed `N`
have the *same* `Q`.

### 6.2 Conservation

**Theorem 6.1.** `Q` is invariant under every promotion, emission and
absorption; it changes only when a quantum crosses between objects, and then
oppositely on the two, so `sum Q` is invariant. It changes otherwise only by
graph surgery (`Z` moving when contacts are made or broken).

*Proof:* L3 makes every matter event a promotion, which leaves `N` fixed; `Z`
depends on `G` alone. Hence `Q = Z − N` is fixed. Transfer moves one quantum
from one object to another: `N` rises by one and falls by one, so `Q` falls by
one and rises by one. (`TestPromotionConservesCharge`, `TestTransferConservesTotalCharge`.) ∎

The physical reading: **light is neutral by construction** in this law, because
radiation is not a quantum leaving the count — it is the *difference* paid
during a promotion. Charge transfer requires contact, i.e. chemistry.

### 6.3 Charge is circulation, and circulation needs cycles

The occupancy defect has a dynamical face. On `C_N` the harmonics
`chi_m(j) = exp(2 pi i m j/N)/sqrt(N)` carry link current `circ_m = sin(2 pi m/N)`.

**Theorem 6.2 (shell neutrality).** `circ(t) + circ(2 pi − t) ≡ 0`, so a **full
`±m` doublet carries exactly zero net circulation**; and `circ = 0` for both
singlets `m = 0` and `m = N/2`. Hence net circulation is carried **exactly by
the unclosed shell**.

*Proof:* `shells.mac` §S2/S2b, both identities returning `0` symbolically;
confirmed numerically at `1e−16` for `N = 5..12`
(`TestFullShellIsNeutralHalfShellIsCharged`, `TestZoneBoundaryModeIsNeutral`). ∎

Measured circulations, `C10`: full shells `~1e−16`; half-filled doublets
`+0.588, +0.951, +0.951, +0.588`. The integer label is the winding `±m`; the
magnitude is the current.

Note the zone-boundary subtlety, which the theory gets right rather than
papering over: the `m = N/2` mode is self-conjugate (real, standing), its
winding label is ambiguous mod `N`, and its physical current is exactly zero.
It is therefore a **neutral singlet like `m = 0`** — the sign of a charge needs
a genuine counter-rotating partner to be defined against.

**Theorem 6.3 (charge needs cycles).** A forest has `b1 = 0`, a non-degenerate
Laplacian spectrum, capacity 1 in every shell, and therefore **no half-filled
multiplet at any occupancy** — no residual circulation is possible.

(`TestTreeSpectrumIsNonDegenerate`, `TestTreeHostsNoCirculatingCharge`.)

This is why the ring is the programme's charge-capable object and a chain is
not, and it says something concrete about free-cell: **a hoard whose contact
graph is tree-like cannot be charged, however much `Em` it holds.**

### 6.4 What charge does not yet do

`Q` is quantized, signed, additive and conserved. It does **not** yet have a
force law. The sketch — an unclosed shell can lower its energy by completing
against a complementary defect, so opposite `Q` should bind and like `Q`
should not — is plausible and *not derived here*. Until it is, this is a
conserved quantum number, not yet electromagnetism. Named again in §10.

---

## 7. One exact size

**The coarse fixed point is not a species.** §3.4 gives `x* = xb`: a continuum,
one root for every bath. Discreteness enters through the shell ladder only.

**Theorem 7.1.** The neutral occupancies are exactly the magic numbers
`Z_j = sum_{k<=j} c_k`. A closed configuration is stationary at mid-gap, its
one-quantum stability window is exactly the shell gap `Delta_j`, and its
residence time is `~exp(Delta_j / 2T) / 2`.

*Proof:* `monodispersity.mac` §M2–M3. ∎

| `Delta/T` | 1 | 4 | 8 | 16 |
|---|---|---|---|---|
| `tau` | 0.8 | 3.7 | 27.3 | 1490.5 |

And the mass distribution `P(N) ~ exp(−(E(N) − mu N)/T)` on a `C10` well
(magic `1,3,5,7,9,10`, `mu = 1.7`):

| `T` | peak `N` | `P(peak)` | reading |
|---|---|---|---|
| 0.05 | **7** | 0.800 | **species** — one magic size dominates |
| 0.30 | 7 | 0.301 | broad |
| 1.50 | 7 | 0.153 | hierarchy, no line |

(`TestColdDistributionCollapsesOntoAMagicNumber`, `monodispersity.mac` §M1/M4.)

**Consequence.** *Monodispersity is `Delta/T`, and nothing else.* A soft outtake
ceiling can erode a hierarchy from above; it can never manufacture a line. Face
C's broad long-tailed hoard spectrum is the **hot branch of this same
distribution** — the prediction is that cooling it must reveal magic peaks at
the contact graph's own partial capacity sums (`SIM_PREDICTIONS.md` P4).

---

## 8. Temperature

`T` is the parameter of the max-entropy occupancy that reproduces the **bound
stack's own** energy at fixed quantum count. It is recovered by inversion and
round-trips to `1e−3` over `T = 0.05..6`
(`TestStackTemperatureRoundTrips`, `TestStackEnergyIsMonotoneInTemperature`).

This is `T_HARMONICS_AND_WIND.md` §1 made exact, and it settles the correction
recorded there: **`T` is an attribute of the harmonic stack, not the flux.**
`I` and `O` are what the open world sees *once* the stack has a `T`.

**A field with energy but no bound stack has no temperature at all** — there is
no stack to fit, so no `T` exists to be assigned. That is the formal reason a
"hot bath" arm box-fills rather than warming a particle, and it is the same
statement as Theorem 3.4 from the other side.

---

## 9. What this law already explains in the measured record

| measured | this law says |
|---|---|
| noise arms **box-fill**, `I≈O` at capacity but no particle (`MECHANISM.md` S07–S19) | Theorem 3.4: `g=0` ⇒ `x* = xb` exactly. A random bath *cannot* do better. |
| **door traffic scrambles winding** (Face B) | L5: a point-contact door has flat support over every `Δm` — mean transfer 0, sd 3.49 on `C12`. Zero mean is *not* neutrality. |
| **organized field is different from hot field** (SHEAR-G) | L5: a wound door with `k` allows exactly `Δm = k`, mean `+2.000`, sd `0.000`. It is a clean circulation pump. |
| **dark bath ⇒ no door fire** (Face B, G03 `empty`) | Corollary 3.2: a closed cold stack cannot emit; with no field it cannot absorb. Darkness is the ground state. |
| **hierarchical hoards, no single size** (Face C) | §7: that is the `Delta/T <~ 1` branch of one distribution. |
| **radiance selects by erosion, not by minting** (Face C, `Em*≈1.55`) | §7: a ceiling truncates; only a gap makes a line. |
| **the C6 ring holds `W=+2`; the medium dephases it** (RING_DNLS) | §6.3: the ring is the smallest charge-capable graph; contact with a scrambling door is exactly the operator that destroys `m`. |
| **π-flatness of stable matter** (g4) | consistent: a closed shell has zero net circulation, so a neutral object sources no residual current at range. Not derived here, but not contradicted. |

Four independently measured nulls fall out of the same two theorems. That is
the main evidence that the unification is not cosmetic.

---

## 10. What this law does NOT do — the honest list

1. **No force law for `Q`.** Charge is conserved and quantized; attraction
   between opposite defects is sketched (§6.4), not derived. Without it this
   is a quantum number, not electromagnetism.
2. **The graph is given, not selected.** Monodispersity is conditional on
   *which* contact graph a hoard adopts. The theory says magic numbers are
   partial capacity sums; it does not yet say which graph forms. This is the
   single largest gap, and it is exactly the programme's standing "geometry is
   not quantized" problem (`METABOLIC_MATTER.md` §3.2).
3. **`g` is coarse.** The stimulated term exists microscopically (the `1+N_w`
   of L4); the mapping from a field's phase-matching to a number `g` is
   postulated, not computed.
4. **`Nu = 1` is imposed.** It is the cap the adopted table already carries,
   written at mode grain — not derived from anything deeper.
5. **The binding threshold is imported** from the `Es` well depth. `Z` is only
   as well-defined as that threshold.
6. **No spin, no statistics.** The `±` circulation gives a two-valued label;
   nothing here derives exchange behaviour from it.
7. **No dynamics.** This is a law of states, rates and invariants. It is not
   a stepper, and it makes no claim about transport, mobility, or the
   binding/mobility antagonism that `RING_DNLS.md` records.
8. **`d_s = 3` is assumed for the field sector, not measured.** Until the
   freecell contact graph is measured (P1), `p_rad = 4` remains an adopted
   number that this note has merely *explained*, not *confirmed*.

---

## 11. Constraint compliance

| constraint | how it is met |
|---|---|
| one law, energy changes mode | L3: every event is a promotion plus a field payment; `N` conserved, energy booked |
| no background | the graph is re-derived each instant; every invariant used (`c_k`, `Z`, `b1`, `d_s`) is a property of the *current* contact set, and §6.2 explicitly allows `Z` to change under surgery |
| no imported U(1) / electrons | charge is `Z − N`, built from eigenspace dimensions; the `s/p/d` ladder of §5 comes from an automorphism group, not from an orbital |
| free-cell is the substrate | the graph is the contact graph; `ReturnDimension` is designed to run on it directly |
| no pre-v89 material | none consulted; no Q-ball, no continuum PDE, no fixed lattice |

---

*End. The theory is coherent and it retrodicts four measured nulls, but it is
gated by one unmeasured number (`d_s` of the freecell foam) and one unclosed
mechanism (a force law for `Q`). Both are stated as bars in `SIM_PREDICTIONS.md`,
and the queue recommendation is in `VERDICT.md`.*
