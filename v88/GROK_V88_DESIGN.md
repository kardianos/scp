# v88 — Multi-dimensional harmonic fabric (design)

**Date:** 2026-07-27 · Seat: grok-4.5 · Working dir: `v88/`
**Supersedes as mechanism candidate:** `fabric_trap.c` (DNLS / discrete breathers)
**Builds on:** `ONTOLOGY.md`, `FABRIC_DESIGN.md`, measured facts therein
**Code:** `v88/fabric_harmonic.c` (instrument), this file (theory + design)

Tags: **[D]** derived · **[M]** measured · **[P]** postulated · **[G]** guess · **[C]** conjecture

---

## 0. What this document is for

The Q-ball programme (v66–v87) is rejected for measured reasons: no bound
internal spectrum, continuous Q, stability of an *interval*. The first v88
trap model (`fabric_trap.c`) recreated localisation as discrete breathers —
one amplitude, one frequency shift, a continuous dial with a threshold. That
mimics trapping without supplying species, mass ratios, or integer labels.

This document formulates the ontology’s actual mechanism:
**multi-dimensional harmonics + complete-cycle transfer**, as a precise state,
a matching condition where integers enter, a lattice model that can be coded,
candidate species, and annealing schedules that would *find* them rather than
seed them.

---

## 1. Theory — the multi-dimensional harmonic cell

### 1.1 What a cell is  **[P, structured]**

A cell is not a point holding a scalar. It has:

1. **External geometry** — the fabric displacement that *is* the tessellation,
   ```
   x(i) = a·i + ξ(i),     θ(i) = (div ξ)_i
   ```
   with θ < 0 = densified / shrunk. Measured: cyclic energy tightens
   (r_half 29.3 → 2.3 as ω 1.32 → 1.4955; reorganisation past the VK turn). **[M]**

2. **Internal multi-dimensional cyclic process** — the cell has an *internal*
   configuration space of dimension `m ≥ 2` (not the fabric dimension `d`).
   The free motion on that space has a **discrete harmonic spectrum** because
   the internal space is compact (or bounded with zero boundary data): mode
   indices are integers by topology of the internal manifold, not by tuning.

**Postulate (internal space).** Each cell carries an internal compact manifold
`M` of dimension `m`. The free cyclic process is a field (or point particle)
on `M`. Expanding in the Laplace–Beltrami eigenbasis of `M` gives modes
labelled by integer multi-indices `α ∈ ℤ^m` (or a discrete spectrum
`{λ_α}` with integer degeneracies). **[P]**

The implementable reduction used below keeps a finite truncation
`α = 1..m_modes` of that spectrum (the lowest harmonics). That truncation is
a *computational* cut, not a new species list.

### 1.2 State of a cell  **[P → concrete]**

For cell `i`:

| symbol | meaning |
|---|---|
| `ξ_i ∈ ℝ^d` | fabric displacement (geometry DOF); may freeze to scalar `θ_i` for first instruments |
| `θ_i = (div ξ)_i` | local density; mass-bearing |
| `z_{iα} ∈ ℂ`, `α = 1..m` | complex amplitude of internal harmonic `α` |
| `I_{iα} = \|z_{iα}\|²` | action / cyclic energy in mode `α` |
| `φ_{iα} = arg z_{iα}` | phase of mode `α` |

No imported scalar “Higgs field.” The interior configuration of a particle is
a *pattern of occupied modes and locked phases* across a cluster of cells —
emergent, intra-particle, kernel-free. **[P, ontology]**

### 1.3 Isolated-cell free motion  **[P, with measured contact]**

Frequencies depend on cell size (tightening):

```
ω_α(θ) = ω̄_α · (1 − σ θ)          # θ < 0 ⇒ densified ⇒ faster cycle  [P]
```

with `ω̄_α` the free spectrum at θ = 0. The first honest choice is a pure
harmonic tower from an internal circle/interval:

```
ω̄_α = α · ω₀                         # α = 1,2,...,m   [P, simplest]
```

A more geometric choice (same structure, less special): `ω̄_α` = eigenvalues of
a fixed internal Laplacian on a reference cell of unit size; scaling with θ
implements “size maintains the cycle / cycle maintains the size.”

**Cycle ↔ size force (tightening law).** Cyclic energy drives densification
until a hard core, matching the measured branch and the ontology’s rule:

```
θ_eq(I) = −θ_max · tanh( γ Σ_α w_α I_α )     # more cyclic energy → tighter  [P]
```

with weights `w_α = α²` (higher harmonics cost more geometric strain — **[G]**;
any positive monotone weights preserve the sign of the effect). Dynamics for θ
is either overdamped relaxation toward `θ_eq` plus elastic neighbour terms, or
second-order with mass (breathing). Past a threshold of total cyclic energy,
`θ_eq` saturates at `−θ_max`; further energy cannot densify and must reorganise
into *which modes are occupied* — release into configuration. **[P]**; the
VK-turn reorganisation on the old branch is the measured analogue of that
saturation. **[M, analogy only]**

### 1.4 What a “complete cycle” is between two cells  **[D from averaging]**

Consider two neighbouring cells `i, j` with mode sets and frequencies
`{ω_{iα}}`, `{ω_{jβ}}`. Any bilinear or multi-linear coupling of the form

```
H_int = ε · z_{iα}^{*p} z_{jβ}^{q} + c.c.     (p,q positive integers)
```

produces a *secular* (net, long-time) energy exchange if and only if the
**resonance condition**

```
p · ω_{iα}(θ_i)  =  q · ω_{jβ}(θ_j)           (1)
```

holds to within a detuning smaller than the coupling bandwidth. Off resonance
the relative phase `p φ_{iα} − q φ_{jβ}` winds at a nonzero beat frequency;
the time average of power vanishes. That average is taken over the common
period of the two motions — i.e. over a **complete joint cycle**.

**Definition (complete-cycle transfer).** Net energy transfer on channel
`(α,β; p:q)` between cells `i` and `j` is the resonant (secular) part of the
Hamiltonian flow generated by `H_int`. Partial-cycle instantaneous power may
be nonzero; *net* transfer over a joint period is zero unless (1) holds.
**[D]** from standard averaging / rotating-wave theory.

This is the precise content of the ontology’s governing rule, without a
non-Hamiltonian “hard gate.” A hard stroboscopic gate (transfer only at
phase-match instants) is a possible discrete approximation; it is not required
for integers to enter, and it risks artificial dissipation. Prefer resonant
Hamiltonians. **[P, design choice]**

### 1.5 Closure / matching condition — where integers enter  **[D]**

Three integer structures, independent of amplitude thresholds:

**(A) Mode indices.** Internal spectrum labels `α ∈ ℕ` (or multi-indices in
`ℤ^m`) come from eigenfunctions on a compact internal space. **[D from (P)]**

**(B) Resonance vectors.** A channel is labelled by integers `(p,q)` (or a
vector for multi-mode). Secular transfer ⇔ (1). The *order* of the channel is
`r = p + q`; low-order channels dominate if couplings fall with order
(generic smooth interaction). **[D]**

**(C) Internal phase lock (the emergent “Higgs”).** Inside a cell or a
cluster, modes can lock when there exists a nonzero integer vector `n_α` with

```
Σ_α n_α ω_α(θ) = 0                         (2)
```

and a stable fixed relative phase for `Σ n_α φ_α`. That locked multi-mode
state is a **configuration**, not a field. Different integer vectors `n` are
different species candidates. **[D structure; stability is dynamical]**

**(D) Cell count.** A trapped region occupies an integer number of cells
`N_c`. That is a discrete label no continuum soliton has. **[D from discrete
fabric]**

**Matching condition for a particle.** A cluster `C` of cells is a particle
candidate when:

1. Its interior admits a locked configuration (2) with total cyclic energy
   above the densification knee (so θ is near hard core inside).
2. Every low-order resonance channel from interior modes to *exterior*
   fabric modes is detuned outside the coupling bandwidth — no complete cycle
   connects interior to exterior.
3. The exterior carries a sharp elastic / frequency distortion (θ and ω
   jump) that maintains that detuning.

Then energy reflects internally. Trapping is a **spectral mismatch of a
multi-frequency object against the exterior band**, not an amplitude going
outside a single linear band. **[D from 1.4 + ontology]**

### 1.6 Why this is not a discrete breather  **[D]**

| | DNLS breather (`fabric_trap`) | multi-harmonic cell (this) |
|---|---|---|
| DOF per cell | 1 complex mode | `m ≥ 2` modes + geometry θ |
| localisation | amplitude shifts one ω outside linear band | multi-mode lock + exterior detuning of *channels* |
| “discreteness” | threshold on a continuous dial | integer resonance vectors + cell count |
| species | one family parametrised by amplitude | distinct occupied-mode sets / lock vectors |
| mass ratio | continuous | ratios of locked energies / cell counts |
| 3–7 dim claim | none (works all d) **[M]** | can couple to elastic sector where d>2 is required |

A breather *mimics* localisation. It has nothing for integers to label except
maybe an action that is still continuous.

### 1.7 Dimension restriction 3–7  **[SUPERSEDED — see v88/ONTOLOGY.md C6]**

> **This section answers the wrong question.** It reads "3–7" as a RANGE of
> workable fabric dimensions and scans d. The user's structure is **3 emergent
> dimensions from 7 degrees of freedom**, tied to the cross product existing
> only in 3 and 7 dimensions. Scans over d measure nothing. The material below
> is retained only as a record of the misreading.

**[original, superseded: D lower; open upper]**

Two different “dimensions” must not be confused:

| symbol | meaning |
|---|---|
| `d` | fabric / embedding dimension (lattice lives in `ℤ^d`) |
| `m` | internal dimensionality of one cell |

**Lower bound on fabric `d`:** elastic self-energy of a dilational inclusion.
Strain ~ r^{1−d}, self-energy integral ∫ r^{1−d} dr converges at infinity only
for `d > 2`. So localised density lumps with finite elastic cost need
**`d ≥ 3`**. **[D]** (stated already in `fabric_mass` / ONTOLOGY). Static
`fabric_mass` still *coarsens* in d≥3 — the bound is necessary, not
sufficient. **[M]**

**Internal `m`:** `m = 1` collapses to single-mode / breather territory.
`m ≥ 2` is required for nontrivial integer resonance vectors. **[D]**
Nothing in the present construction *forces* an upper bound on `m`.

**Upper bound near 7:** **not derived.** Possible [G] reasons that are
testable, not theorems:

* High fabric `d` raises coordination number and the number of exterior
  channels; trapping may become harder (more leakage). Measure fraction of
  energy retained in seeds vs `d`.
* Internal `m` large → dense resonance web → almost everything couples to
  everything (Arnold web); locked islands shrink. Measure lock lifetime vs `m`.

**Honest status:** the 3–7 *window as a hard fact* does not follow from the
multi-harmonic + complete-cycle rules alone. What follows is: **`d ≥ 3`**
(elastic) and **`m ≥ 2`** (nontrivial harmonics). The upper edges are
empirical questions for the instrument below. The DNLS all-d result **[M]**
is expected and irrelevant once density coupling and multi-mode content are
present.

### 1.8 Tensions / unworkable points  **[honest]**

1. **Hard vs soft complete-cycle rule.** A literally hard rule (“transfer is
   exactly zero except at phase-match events”) is non-generic and usually
   non-Hamiltonian. The resonant/averaged form is the mathematically clean
   reading and is what this design uses. If the user requires a hard gate, say
   so — it changes the integrator and conservation laws.

2. **Imported mode count.** Truncating to `m` modes looks like importing `m`
   fields. Defence: `m` is a cutoff of a geometric spectrum, not a species
   list; species are *integer lock vectors and occupancy patterns*, not mode
   slots. Still a risk: if results depend sensitively on `m`, the cutoff is
   physical and must be derived from cell geometry.

3. **Commensurate free spectrum.** Choosing `ω̄_α = α ω₀` makes every
   harmonic pair resonant at uniform θ. That is convenient and may be *too*
   resonant: the vacuum itself is a resonance web. Prefer a weakly detuned
   tower `ω̄_α = α ω₀ (1 + ε_α)` with small irrational ε_α, so resonance is
   special and requires densification or configuration to tune into. **[P]**

4. **Variable size vs fixed stencil.** Large θ changes physical neighbour
   distances; a fixed logical stencil with solved weights (procedural fabric)
   is still valid, but topology of contact does not change. True foam with
   neighbour birth/death is a later stage.

5. **Mass = cell count** is clean only if particles have sharp boundaries and
   little exterior weight. If most energy lives in a long-range elastic skirt,
   mass is not purely combinatorial. Measure both `N_c` and total energy.

---

## 2. Simulation design

### 2.1 Architecture  **[P]**

- **Geometry:** start on a cubic logical lattice for mechanism tests; production
  geometry is procedural irregular (`fabric_proc`) or Kelvin BCC 8+6
  (`FABRIC_DESIGN.md`), with fixed logical stencil and solved weights
  `Σ_j w_j d_j ⊗ d_j = 2I`. Spacing `a` is a physical constant — do not
  “converge dx → 0.” Target particle width / `a` ≈ 1–3 (discrete regime).
- **No kernel change.** Standalone C under `v88/`. No `scp_sim`, no `sfa.h`
  changes, no v71 Higgs sector.
- **Regular computation:** neighbour offsets fixed; only weights/vectors may
  vary per site when irregular geometry is on.

### 2.2 Degrees of freedom per cell

```
θ_i ∈ ℝ                         # dilation (geometry; mass-bearing)
z_{iα} ∈ ℂ,  α = 1..m           # internal harmonics (m = 3 default)
```

Optional later: full `ξ_i ∈ ℝ^d` with `θ = div ξ` and the geometric identity
`Σ_i θ_i = 0` on periodic domains.

Default parameters (all dimensionless with `a = 1`, `ω₀ = 1`):

| symbol | default | role |
|---|---|---|
| `m` | 3 | number of retained harmonics |
| `ω̄_α` | `α(1 + 0.03·α)` | weakly detuned tower |
| `σ` | 0.4 | frequency–density coupling |
| `θ_max` | 0.8 | hard-core densification |
| `γ` | 1.0 | cyclic energy → tightness |
| `w_α` | `α²` | mode weight in θ_eq |
| `ε_1` | 0.15 | 1:1 neighbour coupling |
| `ε_2` | 0.05 | 2:1 / 1:2 neighbour coupling |
| `g_int` | 0.08 | intra-cell 3-wave (1+2→3 style) |
| `μ` | 0.5 | elastic stiffness for θ |
| `η_θ` | 0.3 | θ damping (overdamped option) |

### 2.3 Hamiltonian and update rule  **[P, concrete]**

**Energy:**

```
H = Σ_{i,α} ω_α(θ_i) I_{iα}
  + Σ_i ½ K (θ_i − θ_eq(I_i))² + (Λ/4)(θ_i² − θ_max²)²_{optional hard wall}
  + (μ/2) Σ_{〈ij〉} (θ_i − θ_j)²
  + Σ_{〈ij〉} Σ_α  ε₁ (z_{iα}^* z_{jα} + c.c.)
  + Σ_{〈ij〉} Σ_{α,β: 2α=β} ε₂ (z_{iα}^{*2} z_{jβ} + z_{jα}^{*2} z_{iβ} + c.c.)
  + Σ_i g_int (z_{i1} z_{i2} z_{i3}^* + c.c.)     # internal lock channel
```

(First instrument may drop the optional hard-wall double-well if `θ_eq` plus
elastic already saturates.)

**Complex equations** (canonical: `i ż = ∂H/∂z*`):

```
i ż_{iα} = ω_α(θ_i) z_{iα}
         + ε₁ Σ_{j∼i} z_{jα}
         + (2:1 terms from ε₂)
         + (intra-cell 3-wave terms from g_int)
```

**Geometry** (overdamped, first instrument):

```
θ̇_i = −η_θ [ K (θ_i − θ_eq(I_i)) − μ Σ_{j∼i} (θ_j − θ_i) ]
```

Integrator: RK2 or RK4 on `z`, forward Euler or RK2 on `θ`, `dt ≈ 0.01–0.02`.
Conserve `Σ I` only approximately under θ dynamics (geometry is not purely
canonical in the overdamped reduction); track total `H` and total action as
diagnostics. A fully Hamiltonian second-order θ is the upgrade if conservation
must be machine-tight.

### 2.4 Conserved / tracked quantities

| quantity | status |
|---|---|
| total cyclic action `Σ_{iα} I_{iα}` | ≈ conserved if ε-couplings are pure exchange and θ slow |
| total energy `H` | conserved if θ is Hamiltonian; drifts if overdamped |
| `Σ_i θ_i` | conserved if θ dynamics is pure exchange Laplacian + odd forces; with `θ_eq(I)` source terms it is **not** — densification is paid by local physics, not zero-sum, unless ξ-based div is used |
| Gauss / charge | N/A — no U(1) gauge here |

**Preferred geometry upgrade:** evolve `ξ` with energy depending on `θ = div ξ`,
so `Σ θ = 0` is geometric. That restores the fabric_mass zero-sum rule inside
this theory. First instrument may use scalar θ for speed; flag the missing
constraint in results.

### 2.5 What would falsify this model

1. **Continuous family only.** After annealing, localised objects fill a
   continuum of actions / radii with no clustering in the fingerprint space
   (defined in §3). → integers are not selecting; collapse to breather physics.
2. **No recurring species.** Independent anneals never produce the same
   fingerprint bins (within tolerance). → no particle types.
3. **m-sensitivity.** Species set changes qualitatively with every change of
   mode cutoff `m`. → cutoff is not a truncation of a geometric spectrum;
   design is wrong or needs derived `m`.
4. **Identical to DNLS.** Setting `m=1`, `g_int=0`, `ε₂=0` recovers
   `fabric_trap` behaviour *and* `m≥2` runs show the same amplitude-only
   clustering with no mode-content structure. → multi-harmonics are inert.
5. **Coarsening.** Number of objects falls without bound in time at fixed
   energy budget (the outstanding `fabric_trap` check). → no preferred size.
6. **No lock.** Intra-cell relative phases never stabilise; (2) never holds
   for macroscopic times. → “Higgs = configuration” does not occur.

### 2.6 Diagnostics (required outputs)

Per snapshot / end state:

* total `H`, total action, max `I`, max `|θ|`
* connected clusters of cells with `Σ_α I_{iα} > I_thr` (or `|θ| > θ_thr`)
* per cluster fingerprint (see §3)
* histogram of fingerprints across clusters and across seeds
* lock quality: min over cluster of `|Σ n_α ω_α|` for candidate `n`
* for coarsening runs: `n_lumps(T)` at T = 200,400,800,1600

### 2.7 Code

Instrument: **`v88/fabric_harmonic.c`**

```
gcc -O3 -march=native -fopenmp -o fabric_harmonic fabric_harmonic.c -lm
./fabric_harmonic [schedule] [d] [seed]
```

Schedules: `quench`, `thermal`, `cyclic`, `driven`, `staged`, `coalesce`,
`smoke` (short self-test). See §4.

---

## 3. Emergent particles — candidates, labels, mass ratios

### 3.1 Fingerprint of a cluster  **[P]**

A species label is the discrete part of:

```
F = ( N_c,  s_α = round(Ī_α / Ī_tot, k bits),  n̂,  p_ext )
```

| piece | meaning | type |
|---|---|---|
| `N_c` | number of cells in cluster | integer |
| `s_α` | fractional action in mode α (coarsely binned) | rational / bin index |
| `n̂` | integer lock vector of the dominant internal resonance (normalised, gcd 1) | integer vector |
| `p_ext` | lowest exterior resonance order that is near-detuned (leak channel) | integer |

Amplitudes enter only through *which bins are occupied*, not as continuous
coordinates of the species. Continuous leftovers (total energy within a bin,
orientation) are moduli of a given species, not new species — if the design
works.

### 3.2 Candidate species (minimal list)  **[P / G]**

Assuming `m = 3` and channels 1:1, 2:1, internal 1+2↔3:

| id | occupancy pattern (schematic) | lock | expected size scale |
|---|---|---|---|
| **S1** | mostly α=1 | trivial single-mode; *suspect breather* | small; continuous family risk |
| **S2** | α=1+2 locked 1:2 | n = (2,−1,0) | discrete if lock required for trap |
| **S3** | α=1+2+3 locked via 1+2↔3 | n involving (1,1,−1) | “first genuine multi-mode” |
| **S2′** | α=2 dominant with 2:1 to exterior blocked | higher harmonic core | heavier / tighter than S1 at same action |
| **C_k** | bound cluster of k of the above | composite lock across cells | annealing / coalesce products |

**S1 is the danger state:** if S1 dominates and fills a continuum, the model
has failed into breathers. Success requires S2/S3 (or composites) to recur
with *narrow* `N_c` and `s_α` histograms.

### 3.3 What sets mass ratios  **[D structure; numbers are G]**

Within a locked configuration, actions are constrained by the lock and by
θ_eq saturation. Schematic:

```
E ≈ Σ_α ω_α(θ*) I_α ,   θ* ≈ −θ_max inside,
I constrained by lock manifold and total action budget of the anneal.
```

Then mass ratios of two recurring species A,B are predicted as

```
M_A / M_B ≈ E_A / E_B   or   N_c(A) / N_c(B)
```

whichever diagnostic is stable across seeds. Prefer `N_c` if the ontology’s
“mass = cell count” holds; report both.

### 3.4 One falsifiable numerical prediction  **[P → testable]**

**Prediction P1.** Under quench and thermal schedules at fixed total action
density `ρ_I = Σ I / N` in a window where multiple lumps form, the joint
histogram of `(N_c, bin(s_1,s_2,s_3))` develops **at least two isolated peaks**
that appear in ≥ 3 independent random seeds, with each peak’s `N_c` having
coefficient of variation `< 0.15` inside the peak.

If only one broad ridge in amplitude appears (CV ≳ 0.25, no multi-mode
structure), P1 fails.

**Prediction P2 (dimension).** With density coupling on, stable multi-lump
states with CV(`N_c`) < 0.15 exist for `d = 3` and do **not** exist for
`d = 1` (elastic / coordination). `d = 2` is marginal. This is the first
place a 3-ish lower bound can show up — absent in pure DNLS. **[P]**

**Prediction P3 (mass ratio).** If S2 and S3 both recur, their mean cell counts
satisfy `N_c(S3)/N_c(S2) ∈ [1.2, 2.5]` under the default parameters above —
order-unity, not 10× — because both sit near the same hard-core θ and differ
by mode content, not by a deep hierarchy. A measured ratio ≫ 10 or ≪ 1
without parameter retuning falsifies the “same knee, different configuration”
picture. **[G bounds; order-unity is the claim]**

---

## 4. Annealing / organisation

Success criterion (agreed): **a species that recurs across independent anneals
is the first real evidence of a particle type.**

### 4.1 Common protocol

* Lattice: `d ∈ {1,2,3,4}`, `L` so `N ≈ 10^4–4·10^4` (match fabric_trap sizes).
* Seeds: ≥ 5 independent RNG seeds per schedule.
* Record fingerprints every `ΔT`, final histogram, `n_lumps(T)`.
* Control: `m=1` run of the same schedule must be reported side-by-side.

### 4.2 Schedules

| schedule | procedure | expected products | discrimination |
|---|---|---|---|
| **smoke** | short T, one seed | code sanity: finite H, some localisation | crash / NaN / flat zero |
| **quench** | random z, θ=0; evolve T_cold with no noise | metastable multi-mode traps; broad then freezing | frozen `N_c` histogram; compare to thermal |
| **thermal** | Langevin noise on z (and weak on θ) with T_noise(t) = T0·e^{−t/τ} | near-ground configurations; lightest species | should favour lowest-order locks (S2/S3) over exotic |
| **cyclic** | K periods of heat (raise noise) / cool (anneal) | defect shake-out; survivors = hard species | species that vanish after cycles are metastable only |
| **driven** | weak coherent pump on α=1 at frequency ω_p, or broadband action injection at rate Γ | larger objects if growth needs feed; threshold Γ* | plot `N_c^max(Γ)`; step = growth threshold |
| **staged density** | raise mean action density ρ_I in stages; anneal to steady at each | heavier / higher-harmonic species only above thresholds | ladder steps in fingerprint vs ρ_I |
| **seeded coalesce** | place 2–4 preformed S2/S3 seeds at separations 2..6; anneal | bound composites C_k or elastic molecules | binding if final `N_c` = sum and fingerprint stable; scatter if not |
| **detune control** | same as quench but `ε₂ = g_int = 0` | *should* collapse toward S1 continuum | if still “species”, labels are fake |

### 4.3 How to tell species apart in data

1. Cluster cells with `I_i > I_thr`.
2. Compute `N_c`, `s_α`, best lock vector `n̂` (small integer search `|n_α|≤3`).
3. Hash into bins; build histogram over clusters and over seeds.
4. A **species** = a bin with count ≥ n_min in ≥ 3 seeds, stable under cyclic
   schedule, absent or deformed under detune control.

### 4.4 Priority order (compute budget)

1. **smoke** + **detune control** — instrument validity.
2. **quench** vs **thermal** at `d=3`, `m=3` — first look at recurrence.
3. **coarsening check** T→1600 — preferred size real?
4. **cyclic** — hard vs soft species.
5. **staged** + **driven** — ladder / growth.
6. **coalesce** — composites.
7. Geometry upgrade: procedural / Kelvin stencil; full `ξ`.

---

## 5. Relation to measured fabric geometry

Mechanism tests may use cubic neighbours. Production claims about particles in
an isotropic world need the fabric ranking in `FABRIC_DESIGN.md`:

* Simple cubic excluded for isotropy at discrete widths. **[M]**
* Kelvin BCC 8+6 zeroes T4; pinning still anisotropic. **[M]**
* Procedural irregular: regular compute, solved weights; best current
  engineering path to dynamical δ later. **[M]**

None of that replaces multi-harmonics; it is the *stage* they act on. Pinning
is desirable (localisation); crystallinity of mobility is the cost.

---

## 6. What is derived vs postulated vs guessed (summary)

| claim | tag |
|---|---|
| Continuum PDE ⇒ continuous families | **[D]** + **[M]** (Q-ball census) |
| Discrete fabric needed for discrete mass | **[P]** (ontology); cell count integer is **[D]** once fabric is discrete |
| Complete-cycle ⇔ resonant secular transfer | **[D]** (averaging) |
| Integers from mode indices + (p:q) + lock vectors | **[D]** given compact internal space **[P]** |
| Higgs = locked multi-mode configuration, not a kernel field | **[P]** (ontology) |
| Tightening law sign | **[M]** (branch) + **[P]** (implementation) |
| `d ≥ 3` for finite elastic self-energy | **[D]** |
| Upper bound ~7 | **not derived** |
| Default parameters / P3 ratio window | **[G]** |
| Recurring fingerprint bins = particle types | **[P]** success criterion (agreed) |

---

## 7. Files

| file | role |
|---|---|
| `v88/GROK_V88_DESIGN.md` | this document |
| `v88/fabric_harmonic.c` | multi-harmonic resonant fabric instrument + anneal schedules |
| `v88/fabric_trap.c` | **cautionary** DNLS breather — do not extend as the particle mechanism |
| `v88/ONTOLOGY.md` | user ontology (authority) |
| `v88/FABRIC_DESIGN.md` | geometry measurements |

---

## 8. Immediate next measurements (for whoever runs the instrument)

1. `smoke` then `quench` at d=3, m=3, 5 seeds — look for multi-peak fingerprints.
2. Same with detune control (`ε₂=g_int=0`).
3. Coarsening curve `n_lumps(T)`.
4. Only if (1) shows recurrence: cyclic + staged + coalesce.

If (1) fails cleanly, revise the intra-cell lock (g_int, spectrum detuning)
before touching geometry. Do not add a Higgs field.

---

## 9. Instrument status (2026-07-27, this seat)

**Built:** `v88/fabric_harmonic.c` → binary `v88/fabric_harmonic`

```
gcc -O3 -march=native -fopenmp -o fabric_harmonic fabric_harmonic.c -lm
./fabric_harmonic [smoke|quench|thermal|cyclic|driven|staged|coalesce|detune] [d] [seed]
```

**Smoke (d=3, seed=1, T≈6):** action conserved to ~3×10⁻⁴ relative (301.39 → 301.31).
Free+geometry energy falls as overdamped θ densifies (expected with η_θ > 0).
Sparse random init → hundreds of small peaks; not yet organised.

**Quench (d=3, seed=42, T=400), single seed — not a claim of species:**

| t | n_lumps | mean N_c | sd/mean | act | notes |
|---|---|---|---|---|---|
| 300 | 10 | 3.6 | 0.90 | 519 | mixed mode content |
| 400 | 3 | 4.7 | 0.44 | 534 | fewer objects; peak I large |

Observations (provisional, **not** success):

* Laplacian 1:1 conserves action on short times; longer quench still shows
  action/energy growth and large peak I — treat as an **integrator / 2:1
  channel bug or stiffness** until a Manley–Rowe audit is clean. Do not read
  late-time quench as physics until E is stable.
* `n_lumps` falling 10 → 3 by T=400 means **coarsening is not excluded**
  (same outstanding check as for `fabric_trap`).
* Mode-3-dominated cores appear often; a few objects keep multi-mode weight.
  Lock residuals ~0.07 are mostly the built-in free-spectrum detuning, not
  proof of dynamical phase lock — need a phase-variance diagnostic.
* **No multi-seed recurrence test has been run.** Prediction P1 is open.

**Code debt (before trusting anneals):**

1. Full interaction energy in the diagnostic; track Manley–Rowe invariants
   per channel.
2. Phase-lock order parameter (variance of `Σ n_α φ_α` over a cluster window).
3. Smaller `dt` or RK4 when peak I ≳ 1; optional action-projection.
4. `ξ`-based geometry with `Σ θ = 0` (upgrade from scalar θ).
5. Procedural / Kelvin stencil swap (geometry is currently cubic).

**What this seat claims:** a precise multi-harmonic theory, an implementable
lattice model distinct from DNLS breathers, falsifiers, species fingerprints,
annealing schedules, and a runnable first instrument. **Not** a measured
particle spectrum.
