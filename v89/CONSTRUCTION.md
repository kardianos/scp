# Construction of the no-background universe

**Proponent construction.** Subordinate to `v88/PRINCIPLE.md`. Takes the
principle as given and builds outward: regular method, conversion direction
and rate, mode structure, quantitative curvature, and particles.

Tags used throughout:

| tag | meaning |
|-----|---------|
| **[D]** | derived from the principle (and, where noted, from measured facts read under it) |
| **[P]** | postulated to make the structure work — free to revise under audit |
| **[G]** | guess: plausible, not load-bearing |
| **[M]** | measured in project instruments (background instruments; numbers kept, ontology discarded) |

Anything in `DESIGN_CLEAN.md` that assumes a permanent index set, fixed
connectivity, or immortal cell labels is superseded here. Recoverable content
from that design is re-derived without those assumptions, or explicitly
re-postulated in a form that does not reintroduce a stage.

---

## 0. The absolute constraint, made operational

**Energy is never destroyed. It only changes mode. Space is one of the modes.**
**[D, PRINCIPLE §§1–2]**

Therefore there is nothing for anything to move through, sit in, or be located
against. **[D, PRINCIPLE §4.1]**

### Operational test for background

A construction reintroduces a stage if it contains any of the following:

1. A set of labels that persist while their values change (`ψ[i]`, immortal
   cell IDs, fixed arrays of sites).
2. A permanent reference geometry from which physical geometry is a
   deviation (`x(i) = a·i + ξ(i)` — the term `a·i` is the stage).
3. Fixed connectivity that is not itself a dynamical product of the current
   energy structure.
4. A quantity whose conservation *requires* a fixed ambient space (ordinary
   momentum).

**What is allowed to persist:** only scalars and algebraic types that do not
locate anything — total energy, the vocabulary of modes, the integer lattice
of harmonic labels as an algebraic object (not a place).

**Audit of this document is in §8.** Every permanent-looking object is checked
there.

---

## 1. The regular method

This is the content that replaces equations of motion. Nothing has a
trajectory. What exists is a rule for which **successions of structure** are
permitted.

### 1.1 State: the energy complex

**[P]** A **world-state** is a single finite **energy complex** \(C\).

An energy complex is not values on an index set. It is a combinatorial object:

```
C = ( V, E, m, e, h, χ, φ, e_mid )
```

| symbol | type | meaning |
|--------|------|---------|
| \(V\) | finite set | **parcels** of energy currently present |
| \(E\) | finite set of pairs from \(V\) | **channels** — open conversion opportunities |
| \(m(v)\) | mode label | current mode of parcel \(v\) (see §3) |
| \(e(v) > 0\) | real | energy in that parcel |
| \(h(v) \in \mathbb{Z}^{k}\) | integer word | internal harmonic content |
| \(\chi(\varepsilon) \in \{+1,-1\}\) | sign | handedness of channel \(\varepsilon\) |
| \(\phi(\varepsilon) \in \{0,1,\ldots,N_{\phi}-1\}\) | discrete phase | progress through the channel cycle |
| \(e_{\mathrm{mid}}(\varepsilon) \ge 0\) | real | energy mid-cycle on the channel |

**No element of \(V\) or \(E\) is named across transitions.** The sets are the
structure. A transition replaces one complex by another; it does not update
fields on surviving labels. Asking “which parcel is this parcel later?” is a
meaningless question — the same way “where is it now?” is meaningless for a
particle under PRINCIPLE §4.4.

**Why parcels at all?** Energy is never an amorphous total alone: conversion
is local and cycle-gated **[M, link mediation]**, so energy must be present as
structurally distinct portions that can couple through channels. Parcels are
those portions. They are not sites of a tessellation; they *are* currently
existing quanta of mode-bearing energy. When their energy fully converts into
another mode’s structure, they cease.

**Why channels are not a fixed lattice:** A channel exists only when two
parcels are in a conversion relationship *now*. Channels are created and
destroyed by rewrites (below). Fixed neighbour tables are forbidden.

**Internal harmonics \(h(v)\):** **[P, from DESIGN_CLEAN §1, stripped of
background]** Energy in a parcel has a compact internal configuration; free
motion on a compact space has a discrete spectrum, so mode indices are
integers by topology of that internal space. The integer lattice
\(\mathbb{Z}^{k}\) is algebraic structure of the *kind of energy*, not a
spatial background. Dimension \(k\) is fixed by the internal space; the
working choice aligned with DESIGN/ONTOLOGY is \(k\) large enough to carry
chiral plane data (see §3.3).

### 1.2 Total energy and the only law

```
E(C) := ∑_{v ∈ V} e(v)  +  ∑_{ε ∈ E} e_mid(ε)
```

**[D, PRINCIPLE §1]** Every permitted transition \(C \to C'\) satisfies

```
E(C') = E(C).
```

There is no other fundamental law. Everything else is the *form* of allowed
conversion — which successions exist — not an extra dynamical equation.

### 1.3 Atomic rewrites (the generators of succession)

A **permitted transition** is any finite composition of the following atomic
rewrites, each required to obey \(E(C')=E(C)\) and the cycle gate.

#### R1 — Deposit into a channel  **[P; structural form of measured mediation]**

Preconditions:

* Channel \(\varepsilon=\{v,w\}\) has \(\phi(\varepsilon)=0\) and \(e_{\mathrm{mid}}=0\),
  or is already mid-cycle and can accept more only if resonance still holds.
* **Resonance:** the harmonic word \(h(v)\) matches the channel’s drive
  condition (integer relation; precise form in §2.2).
* Available energy \(e(v)\) is at least the deposit quantum \(\delta > 0\).

Effect: \(e(v) \leftarrow e(v)-\delta\), \(e_{\mathrm{mid}}(\varepsilon)\leftarrow e_{\mathrm{mid}}+\delta\),
and \(\phi\) advances off 0 if it was 0. Parcel \(v\) is **removed from \(V\)**
if \(e(v)\) hits 0 (it ceases; no ghost label remains).

#### R2 — Advance a channel cycle  **[P; “transfer takes time”]**

Preconditions: \(e_{\mathrm{mid}}(\varepsilon) > 0\) and \(\phi < N_{\phi}-1\).

Effect: \(\phi \leftarrow \phi+1\). No energy moves off the channel. Partial
cycles deliver nothing — this is the gate, not a suppression after the fact.
**[M: detuned link suppresses delivery 227× vs bilinear hop’s ~1.65]**

#### R3 — Deliver from a completed cycle  **[P]**

Preconditions: \(\phi = N_{\phi}-1\) and \(e_{\mathrm{mid}} > 0\).

Effect: all (or a quantum of) \(e_{\mathrm{mid}}\) moves to the far parcel \(w\).
If \(w\) does not exist as a suitable receiver, a **new** parcel is created in
the mode determined by §2–§3 (birth of structure). Then \(\phi \leftarrow 0\),
\(e_{\mathrm{mid}}\leftarrow 0\). If the channel has no remaining role, it is
**deleted**.

Delivery is the only way energy leaves a channel. Incomplete \(\phi\) never
delivers. **[D from PRINCIPLE §3 + measured cycle gate]**

#### R4 — Split / merge of parcels  **[P]**

Energy in one mode may re-partition among parcels of that mode when channels
complete in a way that demands structural rearrangement (e.g. one receiver
becomes two coherent parcels, or two parcels lock into one). Split/merge
conserves \(E\), conserves total harmonic content as an additive integer word
where applicable, and does not reference any external coordinate.

#### R5 — Mode conversion of a parcel  **[P; details in §2–§3]**

A parcel’s mode \(m(v)\) may change when a structural condition on its
harmonic word and its incident completed cycles is met. Energy amount is
unchanged by the relabeling; only the mode changes. This is conversion in the
sense of PRINCIPLE §1: the same energy, different mode.

#### R6 — Channel birth / death  **[P]**

A channel between two existing parcels may appear when their harmonic words
admit a resonant chiral coupling (transverse-plane condition, §2.2). A channel
dies when it is idle (\(\phi=0\), \(e_{\mathrm{mid}}=0\)) and resonance fails, or
when either endpoint ceases.

### 1.4 What “dynamics” means

There is no vector field on a manifold of configurations with immortal
coordinates. There is a **binary relation** on the class of energy complexes:

```
C  ⟹  C'    iff    C' is reachable from C by a finite sequence of R1–R6
                   each obeying E-conservation and the cycle gate.
```

A **history** is a chain \(C_0 \Rightarrow C_1 \Rightarrow C_2 \Rightarrow \cdots\).

**Motion, persistence, and change** are all patterns in such chains — never
properties of a labeled object sliding on a stage. **[D, PRINCIPLE §4.4]**

### 1.5 Branching and regularity

The relation \(\Rightarrow\) need not be functional: several rewrites may be
fireable at once. **[P — selection among concurrent rewrites]** When several
atomic rewrites have satisfied preconditions, all that **commute** (disjoint
support in the complex) may be regarded as jointly fireable; when they
conflict, the rewrite that most reduces the open-cycle deficit \(\Delta(C)\)
(defined in §2.1) is selected. If deficit reduction ties, both futures are
physically allowed (branching). No hidden variable chooses; the theory’s
predictions are about **which structures recur**, not about a single
trajectory.

This is enough to define a **regular method**: the content of physics is the
characterization of \(\Rightarrow\) and of the invariants of its recurrent
structures. Differential equations on index sets are not an approximation
target; they are a different ontology.

### 1.6 What conservation constrains

Conservation alone does **not** pick the next state. It only forbids
transitions that create or destroy energy. The cycle gate forbids transitions
that move energy without complete chiral cycles. Together:

| filter | forbids |
|--------|---------|
| \(E(C')=E(C)\) | creation/destruction of energy; “teleport” of energy as accounting fiction |
| cycle gate (R1–R3) | instantaneous hop; partial-cycle leakage; non-resonant deposit |
| mode rules (§2–§3) | arbitrary mode flips; species without integer closure |

What remains is a sparse rewrite system. That sparseness is the theory’s
dynamical content.

---

## 2. Direction and rate of conversion

Conservation says energy changes mode. It does not say which way, or when.
This section constructs both.

### 2.1 Open-cycle deficit  **[P]**

Define, for a complex \(C\),

```
Δ(C) = ∑_{ε ∈ E}  ω_φ · 1_{e_mid(ε)>0 and φ incomplete}
     + ∑_{v ∈ V}  ω_h · mismatch(h(v); incident channels)
     + ∑_{v: m(v)=S} ω_s · saturation_excess(v)
```

Schematic weights \(\omega_\bullet > 0\). Precise coefficients are not needed for
structure: \(\Delta\) is a **nonnegative structural count** of unfinished
conversion work.

* **mismatch** vanishes when harmonic words stand in the integer resonance
  relation required by their channels; otherwise positive.
* **saturation_excess** is positive when a space-mode parcel holds more cyclic
  content than its tightening capacity allows (see §2.3).

**Gate (when):** A deposit or delivery rewrite is permitted only when it is a
step that can lead to a decrease of \(\Delta\) upon cycle completion — i.e.
only resonant, completeable cycles fire. **Conversion proceeds where chiral
cycles can close and stalls where they cannot.** **[P, developing PRINCIPLE
§7’s candidate; reinforced by measured 227× gate]**

**Direction among fireable rewrites:** select rewrites that minimize
\(\Delta(C')\) (steepest closure). **[P]**

### 2.2 Resonance and handedness  **[P + M]**

A channel actuates in two transverse degrees of freedom relative to the
conversion it mediates. Writing the two transverse amplitudes as complex
(or as a pair of harmonic occupations),

```
u± = (u₁ ∓ i u₂)/√2
χ  = 2 Im(u₁* u₂)          signed handedness
```

**[M: χ well-defined only where a transverse plane exists; strong at emergent
d=3, suppressed in d=1,2]**

**Resonance condition (integer form) [P]:**

Channel \(\varepsilon\) between \(v\) and \(w\) is resonant when there exist
integers \(n_1, n_2\) such that the transverse content of \(h(v)\) and \(h(w)\)
satisfies

```
n_1 · Ω_1(v) + n_2 · Ω_2(v)  =  Ω_channel
```

and the relative phase of the two transverse occupations matches \(\chi(\varepsilon)\).
Non-integer mismatch ⇒ no deposit (R1 fails). This is the structural gate, not
a rate constant.

Same-handed vs opposite-handed channels couple through different integer
relations — measured as a pair-energy split in background instruments **[M]**;
here retained as: **handedness is part of the resonance word**, so different
\(\chi\) open different rewrite channels.

### 2.3 Direction of mode change  **[D + P]**

Three directed rules, in order of logical force:

**(A) Saturation forces space → dense pattern. [P, from DESIGN_CLEAN §4
without fixed volume budget]**

A space-mode parcel can hold cyclic content only up to a structural capacity
\(e_{\mathrm{cap}}(h)\). Beyond that, further energy cannot remain as free
space-mode cyclic content: conservation requires it change mode. The available
mode that does not require an open exterior cycle is **dense pattern** — energy
locked into a closed integer configuration of harmonics (§3, §5).

So: if \(e(v) > e_{\mathrm{cap}}(h(v))\) and \(m(v)=S\), then R5 to \(m(v)=D\)
(or split + lock) is not optional. \(\Delta\)’s saturation term enforces this
under steepest closure.

**(B) Completed exterior cycles allow dense → space. [P]**

If a channel completes delivery into a dense-pattern parcel in a way that
**breaks** the integer lock (mismatch becomes nonzero), the lock fails and R5
to space mode is fireable. Without such a completed exterior cycle, unlock
does not occur — spectral isolation of particles (§5).

**(C) Open completable chains produce field pattern, not mass. [D from
topology, §3]**

If energy can advance by R1–R3 along a chain without ever satisfying a closed
self-reproduction condition, steepest closure keeps it in open succession:
that is field pattern. It is not a third independent rule so much as what
saturation-free, lock-free conversion *is*.

### 2.4 Rate, without a time stage  **[D + P + M]**

There is no external clock. **Time is the partial order of rewrite events**
(the chain index along a history, refined by causal dependence when rewrites
share support).

**Rate** is therefore always **relative**:

```
c ≡ maximum number of successive R3 deliveries along an open chain
    per unit of succession depth, in natural units of one parcel advance.
```

**[D, PRINCIPLE §4.5: c is the conversion rate, not a travel speed]**

Measured in background instruments: \(c \sim 1.05\, g\, a\), constant over a
4× range in coupling **[M]**. Reading under this ontology:

| background symbol | meaning here |
|-------------------|--------------|
| \(a\) | characteristic energy (or inverse structural density) of one space parcel — the natural unit of “one step of structure” |
| \(g\) | resonance acceptance / cycle-coupling strength entering the number of succession steps needed per delivery |

So \(c \sim g \cdot a\) is the same relation: cycle coupling strength times
structural step size. The numerical factor 1.05 is instrument-dependent; the
**constancy over \(g\)** is the structural claim that \(c\) is set by cycle
completion, not by a posted speed limit. **[M read as D]**

Detuning suppresses delivery by destroying resonance (R1 fails), measured
227× **[M]** — here: almost no R3 events, so almost no advance, so almost no
effective \(c\) along that channel.

### 2.5 Summary of direction and rate

| question | answer |
|----------|--------|
| When does conversion happen? | When a chiral channel is resonant and a full cycle can complete (deficit can fall). |
| Which way does mode change? | Saturation → dense; unlock by exterior cycle → space; open chains → field. |
| How fast? | Relative rewrite depth; max open-chain advance rate is \(c\). |
| What stalls conversion? | Non-resonance, incomplete cycles, spectral mismatch at a boundary. |

---

## 3. Mode structure

### 3.1 The established three, plus one forced fourth

PRINCIPLE §2 lists space, dense pattern (mass), field pattern (EMF). Transfer
that takes time, with energy never destroyed, **forces a fourth**:

**[D]** While a channel holds \(e_{\mathrm{mid}} > 0\) with incomplete \(\phi\),
that energy is not in space, not in dense pattern, and not in a delivered
field pattern at a parcel. It is **transfer mode** \(T\): energy mid-conversion.

Without \(T\), either transfer is instantaneous (contradicts mediation **[M]**)
or energy vanishes mid-cycle (contradicts PRINCIPLE §1). So:

| mode | symbol | what it is |
|------|--------|------------|
| space | \(S\) | energy not locked into closed or open pattern; free internal harmonics |
| dense pattern | \(D\) | energy in a **closed** conversion topology (self-sustaining) |
| field pattern | \(F\) | energy in an **open** conversion topology (advancing chain) |
| transfer | \(T\) | energy mid-cycle on a channel (\(e_{\mathrm{mid}}\), incomplete \(\phi\)) |

\(T\) is not an imported field. It is the only place mid-cycle energy can be
once “energy is never destroyed” and “partial cycles deliver nothing” are both
kept.

### 3.2 What forces the mass / EMF distinction  **[D]**

Both mass and EMF are fabric patterning. The distinction is not substance and
not an extra label. It is the **topology of the conversion graph** that the
energy participates in:

```
closed, self-isomorphic succession  →  dense pattern (mass)
open, support-advancing succession →  field pattern (EMF)
neither (free harmonics only)      →  space
mid-cycle on a channel             →  transfer
```

**Forced, not chosen:** given that conversion is cycle-gated rewrite of
structure, every parcel’s energy either (i) participates in a closed recurrent
pattern, (ii) participates in an open advancing pattern, (iii) sits free, or
(iv) is mid-channel. Those are the logical possibilities for how energy
relates to the rewrite relation \(\Rightarrow\). The names mass/EMF/space/T
are those four.

### 3.3 What distinguishes the *patterning* of mass from EMF  **[P]**

Beyond topology, the harmonic content differs:

* **Dense pattern:** harmonic word \(h\) satisfies a **closed integer system**
  — mutual resonance among a finite set of parcels such that the only
  complete cycles are internal. Handedness may be locked as a net \(\chi\)
  (charge-like invariant; see §5.4). Exterior channels are off-resonance by
  construction of the lock (spectral mismatch).

* **Field pattern:** harmonic word supports **chiral open resonance** along a
  chain: each completed cycle shifts the support (birth of a parcel “ahead,”
  death or mode-return “behind”). Net handedness \(\chi\) can propagate as a
  signed signal without locking into a closed multi-parcel relation.

* **Space:** free spectrum; no saturated lock; channels may open if neighbors
  resonate, but nothing forces a particular topology.

### 3.4 Is the list complete?  **[P]**

Claim: **yes, at the level of conversion topology.** Any further “species” of
energy is a *substructure* of \(D\) or \(F\) (different integer locks or chain
types), not a new mode. Composites are several \(D\)-locks sharing channels.
Vacuum structure is inhomogeneity of \(S\), not a fifth mode.

If a future effect cannot be expressed as rewrite behavior of
\(\{S,D,F,T\}\), the construction has failed — the remedy is not to add a
mode by hand, but to revise the rewrite rules. **[constraint: no imported
fields]**

### 3.5 Retired: fixed volume budget  **[D, PRINCIPLE §5]**

\(\sum_i \theta_i = 0\) treated total space as pinned. Here total **energy** is
pinned; \(E_S\) is dynamical. Space converts to \(D\) and \(F\). Rarefaction is
not required to pay for densification globally; what “pays” is the mode
ledger. Locally, removing space energy still forces the remaining space
complex to rewire — that is curvature (§4), not a global \(\sum\theta=0\).

---

## 4. Curvature, quantitatively

### 4.1 Qualitative derivation restated  **[D, PRINCIPLE §4.3]**

Matter is converted space. Converting space to dense pattern at a locus
removes space-mode energy there. Energy is never destroyed, so the removal is
not a gap in a container — there is no container. The **remaining space-mode
structure** must account for the conversion: its incidence changes. That
change is curvature. Gravity is not a force and not a separate field; it is
conservation expressed in the structure of remaining space.

### 4.2 Intrinsic space complex

Let \(C_S\) be the subcomplex of parcels with \(m(v)=S\) and channels that
presently couple two space parcels (or space to transfer returning to space).
\(C_S\) is an ordinary finite graph (plus harmonic decorations). It has **no
embedding**. Graph distance \(d_S(u,v)\) is the channel-count of the shortest
path in \(C_S\).

When energy \(m\) converts \(S\to D\), some space parcels cease or shrink in
energy until removal, and channels rewire. \(C_S\) changes combinatorial type.

### 4.3 Combinatorial Gauss law  **[P, aimed at recognisable continuum limit]**

**Postulate (defect = converted content):**

For any partition of \(V\) into interior \(I\) and exterior \(X\), let
\(B\) be the set of channels with one end in \(I\) and one in \(X\). Define
the **space flux** through \(B\) by

```
Φ(B) := ∑_{ε ∈ B} σ(ε)
```

where \(\sigma(\varepsilon)\) is the space-structural weight of the channel
**[P: default \(\sigma=1\) per idle space–space channel; reduced when channel
is mid-transfer or off-resonance]**.

Let \(M(I)\) be total dense-pattern energy in \(I\), and let \(E_S(I)\) be
space-mode energy in \(I\). Let \(\varepsilon_0\) be the characteristic energy
of one space parcel at unit structural role **[P: sets units]**.

**Conservation identity for a region that has reached local rewrite rest
(no incomplete cycles touching \(B\)):**

```
M(I)/ε_0  +  E_S(I)/ε_0  +  E_T(I)/ε_0  =  N_slots(I)
```

\(N_{\mathrm{slots}}(I)\) is the number of structural roles interior energy
currently occupies — a pure count from the complex, not from an ambient
lattice. When conversion increases \(M(I)\) at the expense of \(E_S(I)\),
either \(N_{\mathrm{slots}}\) drops (parcels cease) or the boundary flux
reorganizes.

**Gauss form [P]:** the **combinatorial curvature** enclosed by \(B\) is

```
K(I) := d_* · |I_S| − ∑_{v ∈ I_S} deg_{C_S}(v)
```

with \(I_S = I \cap V_S\) and \(d_*\) the characteristic degree of a free space
foam **[P: fixed by the internal chiral structure that realises three
emergent dimensions; not a background lattice constant — it is a property of
how space-mode channels form in the absence of dense pattern]**.

Then the load-bearing relation is

```
K(I) = κ · M(I) / ε_0
```

for regions in rewrite rest with \(E_T\) negligible on \(B\). One constant
\(\kappa\) absorbs the packing of chiral space channels. **[P]**

**Reading:** enclosed dense energy equals total angular defect of the
surrounding space complex. This is discrete Gauss–Bonnet bookkeeping applied
to conversion, not a field equation on a manifold.

### 4.4 What is recognisable  **[D from §4.3 as effective description]**

If one forms continuum averages over many parcels (an *effective* description
only — the continuum is not fundamental and must not be reified as a stage):

1. \(K(I) \propto M(I)\) → integral curvature proportional to enclosed mass
   (Einstein–Hilbert structure with matter source, schematic).
2. Weak, slow variation → Poisson equation \(\nabla^2\Phi \propto \rho_M\) for a
   potential \(\Phi\) built from local cycle-rate shift (next item).
3. Inverse-square free-fall of conversion patterns as the continuum shadow of
   steepest-closure bias on a defect-satisfying space complex (see §4.6).

None of these continuum forms is postulated as fundamental. They are
**recoveries**: what the combinatorial relation looks like when parcel counts
are large and smooth on the scale of interest.

### 4.5 Clock rates and gravitational time dilation  **[P + D]**

Cycle completion rate at a locus depends on local space-mode structure
(DESIGN’s “cyclic energy shifts rates,” rephrased without cell size on a
lattice): where \(E_S\) is depleted and \(K\) is concentrated, the resonant
step count per delivery changes.

**[P]** Define the local **conversion tick rate** \(\nu(v)\) as the number of
R2 advances per unit succession depth at channels incident on \(v\). Then
\(\nu\) is lower (or shifted) where space is converted away. Clocks built from
cyclic processes near dense pattern run at different rates from clocks in
free space. That is gravitational time dilation as pure conservation
geometry — no metric field imported.

### 4.6 Attraction without a force  **[P]**

A particle (§5) advances by open conversion at its boundary: space ahead →
pattern, pattern behind → space. Under steepest closure (§2.1), rewrites that
most reduce \(\Delta\) are selected.

Near a mass, the space complex already carries defect \(K \propto M\). The
particle’s exterior harmonic mismatch is not independent of that defect: the
same integer relations that isolate the particle also **half-match** the
distorted space near another dense region. Converting “toward” the mass
closes more deficit per rewrite than converting “away.”

Therefore the self-sustaining pattern’s succession is biased toward the mass.
Nothing pulls. The regular method, applied in a curved space complex, *is*
free fall.

**Tried alternative that failed:** “less space near mass means harder
conversion ahead, so particles avoid masses.” That would give repulsion and
contradict the geometric reading of PRINCIPLE §4.3. The repair is that what
matters is not raw \(E_S\) abundance but **deficit reduction given existing
defect and spectral match**. Scarcity of space is already encoded as \(K\);
steepest closure on \(K\)-laden complexes yields attraction, not repulsion.
**[P — this is the load-bearing postulate linking Gauss law to free fall]**

### 4.7 Relation between converted energy and distortion (summary)

```
δK = κ · δM / ε_0
```

for a region at rewrite rest: an increment of converted dense energy produces
a proportional increment of combinatorial curvature of remaining space.
Flux form: change in boundary space-flux accounts for the same \(\delta M\)
when interior space energy is held fixed by the ledger.

**No ambient coordinates enter.** \(K\) is degree defect in \(C_S\); \(M\) is
mode energy; \(\varepsilon_0,\kappa\) are structural constants of space-mode
channel formation.

---

## 5. What a particle is

### 5.1 Not a trapped lump

The background picture: a region of a permanent lattice whose amplitudes are
high, held by a potential or a spectral gap *on that lattice*. That picture
is void here — there is no lattice to be trapped on.

### 5.2 Definition  **[P, implementing PRINCIPLE §4.4 + DESIGN §5 without stage]**

A **particle** is a **self-sustaining pattern of conversion**: a finite
substructure \(P\) of an energy complex such that

1. **Closure:** all complete cycles that carry the dense energy of \(P\) are
   internal to \(P\) (or return to \(P\) after a fixed finite rewrite loop).
2. **Self-reproduction:** there exists a permitted succession
   \(C \Rightarrow \cdots \Rightarrow C'\) under which the substructure
   isomorphic to \(P\) recurs (same integer invariants, same mode topology),
   even though no parcel of \(P\) need survive as a labeled entity.
3. **Spectral boundary:** every channel from \(P\) to the exterior fails
   resonance or fails to complete; energy does not leak. The boundary is a
   **failure of cycle completion**, not a surface drawn in a space.
4. **Lock:** the harmonic words of \(P\) satisfy a closed integer system that
   is stable under the internal rewrite loop (saturation completed into
   configuration, §2.3A).

Existence is a threshold: the integer system either closes or it does not.
There is no continuous family of particles. **[D from integer closure; same
moral as DESIGN §5 without immortal cells]**

### 5.3 What makes it persist

Persistence is **recurrence under \(\Rightarrow\)**, not survival of objects.

Internally, R1–R3 fire in a loop that returns the same integer lock. Exterior
R1 fails (off-resonance). Steepest closure is satisfied by the internal loop:
\(\Delta\) inside stays at a minimum compatible with the lock; opening to the
exterior would raise \(\Delta\) or is precondition-blocked.

If an exterior channel ever completes in an unlocking way (§2.3B), the
particle **ceases** as that species — decay, absorption, or conversion to
field/space. Persistence is conditional on isolation of the lock, not eternal.

### 5.4 What makes it discrete  **[D]**

The lock is an integer relation among harmonic words and channel handedness.
The set of solutions in \(\mathbb{Z}^{k} \times \{\pm 1\}^{E_P}\) is discrete.
Species are labeled by:

| invariant | role |
|-----------|------|
| integer harmonic word class \([h]\) | internal configuration / “flavor” |
| parcel count \(N_c\) of the lock | structural size |
| net handedness \(\chi_P = \sum \chi\) | signed structural charge |
| total dense energy \(E_D(P)\) | mass energy |

**Charge is a property of structure** and does not reference position — so it
is conserved across destruction and creation of parcels, matching the
measured charge–momentum asymmetry **[M, PRINCIPLE §6 retrodiction]**.
Momentum is not fundamental: it required a background to define.

### 5.5 What sets mass  **[P]**

**Mass = dense-pattern energy \(E_D(P)\) locked in the self-sustaining
structure.**

When the lock fixes energy per parcel (saturation capacity filled to a
discrete level set by \(h\)), equivalently

```
mass ∝ N_c
```

with a species-dependent factor from the harmonic word. Mass ratios are then
ratios of integer configurations, not values of a continuous dial.
**[aligns with DESIGN §5 “mass is a cell count,” but cells are lock-parcels
of the pattern, not cells of a tessellation]**

### 5.6 Motion of a particle  **[D, PRINCIPLE §4.4]**

A “moving particle” is a history in which the self-reproducing substructure
\(P\) recurs with **advancing support**:

* ahead: space parcels convert (via \(T\)) into the dense lock pattern;
* behind: dense lock energy unlocks into space;
* the integer invariants of \(P\) are preserved along the history.

Nothing traverses anything. Identity across time is reconstruction from
invariants, not a fact about a continuing object.

The maximum advance of support per succession depth is \(c\) (§2.4).

### 5.7 Multi-particle structure and binding  **[P / G]**

Two locks \(P_1, P_2\) in one complex may share exterior channels. If a joint
integer system closes with lower \(\Delta\) than the separated locks, a
**composite** is a new self-sustaining pattern (binding channel). If not,
they are separate particles; steepest closure may still bias their successions
toward each other via the curvature of \(C_S\) (§4.6) — that is gravitational
proximity without a binding lock — or via chiral field-pattern exchange
(open chains of mode \(F\) between them).

**[G]** Attractive vs repulsive non-gravitational interaction is the difference
between joint locks that reduce \(\Delta\) and field-mediated rewrite chains
that raise the cost of approach for given \(\chi\) combinations. Measured
same-vs-opposite handed pair split **[M]** is a seed of this, but was mixed
with self-energy in background instruments; the clean prediction here is:
**only separation-dependent deficit reduction counts as interaction.**

### 5.8 Comparison to DESIGN_CLEAN §5

| DESIGN_CLEAN | this construction |
|--------------|-------------------|
| region on a cell complex whose cycle rate mismatches exterior | self-reproducing substructure under rewrites; no immortal region |
| interior configuration of cell harmonics | closed integer lock of harmonic words |
| sharp exterior density jump | spectral failure of exterior channels (no completed cycle) |
| mass = cell count | mass = \(E_D\), often \(\propto N_c\) of lock parcels |
| trapped by spectral mismatch | same mechanism, stated as rewrite preconditions failing at the boundary |

The physics kept: spectral isolation, integer species, mass from discrete
structure. The ontology dropped: permanent cells, fixed tessellation, density
as a field on a lattice.

---

## 6. Emergent dimension and chirality (supporting structure)

Not in the user’s priority list, but required for channels to have transverse
planes.

**[P, from DESIGN/ONTOLOGY C6]** Handedness needs a cross product. Cross
products of the relevant kind exist in 3 and 7 dimensions. The fabric’s free
algebraic degrees of freedom are 7; the chiral conversion structure realises
**3 emergent dimensions** as the stable open-chain directions of field-pattern
succession.

**Under this construction:** “dimension” is not an ambient \(\mathbb{R}^n\).
It is the count of independent **advance directions** for open conversion
chains — how many inequivalent ways support can propagate under chiral R1–R3.
That count is forced to 3 by the algebra of transverse actuation, not by
embedding in a pre-space.

No scan over \(d\) is part of the theory. Reduced-dimension models remove the
transverse plane and therefore remove the gate that makes the regular method
work; they are not approximations of this construction.

---

## 7. Worked consequences and small derivations

### 7.1 Why charge can be exact while momentum cannot  **[D + M]**

Charge invariants (\(\chi_P\), integer harmonic classes) are properties of the
lock’s combinatorial type. They are unchanged by any rewrite that preserves
the lock, and when locks form or decay they redistribute by integer
bookkeeping. No reference to a stage.

Momentum would be generator of translations through a background. There is no
background, so there is nothing for a Noether momentum to generate. Apparent
momentum conservation in instruments is residual approximate symmetry of a
fixed lattice — integrator-quality, ~1e-4 **[M]** — not a principle.

### 7.2 Why c is universal for free field pattern  **[D]**

Open chains in free space (\(K=0\), uniform resonance) share the same cycle
step count per delivery. Hence one \(c\). Near curvature, tick rates shift
(§4.5): effective propagation of field pattern changes — refraction /
gravitational delay — without changing the local cycle gate.

### 7.3 Schematic mass defect of a composite  **[P]**

If \(P_1, P_2\) bind into \(P_{12}\) with

```
Δ(P_{12}) < Δ(P_1 ⊔ P_2)
```

and energy conserved, part of \(E_D(P_1)+E_D(P_2)\) may convert to field or
space on formation (radiation of \(F\) or release to \(S\)), leaving

```
E_D(P_{12}) < E_D(P_1) + E_D(P_2).
```

Binding energy is mode conversion required to reach the lower-deficit joint
lock — not a potential evaluated at coordinates.

### 7.4 What computation is possible without a stage

Allowed:

* Enumeration of small integer locks: finite search over harmonic words and
  channel handedness for self-reproducing closed systems (see code listed in
  §9). **Done:** 772 closed / 502 isolated species at \(N_c\le 3\), mode bound 2.
* Rewrite graphs: given a finite complex, list fireable R1–R6, compute
  \(\Delta\), step steepest closure.
* Combinatorial \(K(I)\) vs \(M(I)\) on abstract graphs with designated dense
  subsets — check linear Gauss relation in toy foams. **Done:** \(r^2\approx 0.90\)
  for \(K-K_0 \propto M\) on a triangular disk foam.

Forbidden (reintroduces stage):

* Fixed \(N^3\) arrays with evolving `psi[i]`.
* `x(i) = a*i + ξ(i)`.
* PN-barrier scans of a bump through a crystal.
* Any “particle position” time series as fundamental output.

Effective continuum equations may be *written* as asymptotic descriptions of
parcel-count averages; they must not be *simulated as the ontology*.

---

## 8. Background audit

| object | persists? | merely re-valued? | verdict |
|--------|-----------|-------------------|---------|
| parcel set \(V\) | no — created/destroyed | no | OK |
| channel set \(E\) | no — birth/death with resonance | no | OK |
| mode labels \(S,D,F,T\) | as vocabulary only | N/A | OK — types, not places |
| total energy \(E\) | yes, as scalar | no location | OK — conservation |
| harmonic vocabulary \(\mathbb{Z}^k\) | as algebra | occupations change | OK if not spatial |
| characteristic degree \(d_*\), \(\varepsilon_0\), \(\kappa\) | constants of structure | no | OK — not a lattice |
| rewrite relation \(\Rightarrow\) | law | no | OK |
| history index (succession depth) | order of events | not a container | OK — time as rewrite order |
| fixed neighbour table | would persist | would re-value fields | **forbidden** — not used |
| immortal cell IDs | would persist | yes | **forbidden** — not used |
| ambient \(\mathbb{R}^3\) | would persist | yes | **forbidden** — not used |
| momentum | N/A | N/A | **not fundamental** |

**How background was avoided:** state is a finite complex replaced as a whole
by rewrites; identity of parts across rewrites is undefined; geometry is
incidence of currently existing space-mode energy; curvature is degree defect
of that incidence; motion is advancing support of a recurrent pattern;
species are integer lock classes.

**Residual risk:** any software implementation that uses stable integer IDs for
parcels for bookkeeping can *look* like a background. The theory treats those
IDs as ephemeral implementation handles, recycled on death, never as physical
coordinates. The species enumerator in §9 uses only isomorphism classes of
locks — no spatial index.

**Honest limit:** continuum recoveries (§4.4) use language of fields on
manifolds. That language is scaffolding for recognition (“this looks like
Poisson”), not part of the ontology. If it creeps back in as fundamental, the
construction has failed.

---

## 9. Code created

| file | role |
|------|------|
| `v88/construct_species.c` | Enumerate small closed integer locks (harmonic words + handedness) that are self-consistent under the resonance/closure rules of §§2–5. No lattice, no coordinates, no immortal cells — pure combinatorial species count and mass labels \(N_c\). |
| `v88/construct_gauss.c` | Toy check of combinatorial Gauss law: build abstract space graphs, designate dense energy, measure degree defect \(K\) vs enclosed \(M\). No embedding. |

Build:

```
gcc -O2 -o v88/construct_species v88/construct_species.c -lm
gcc -O2 -o v88/construct_gauss   v88/construct_gauss.c -lm
```

### 9.1 Species enumeration result (ran 2026-07-27)

```
./construct_species 3 2
# total_closed_species 772
# isolated_species     502
```

With harmonic components in \(\{-2,\ldots,2\}\) and \(N_c \le 3\), the schematic resonance rule of §2.2 produces a **finite discrete set** of closed locks, of which a majority are spectrally isolated from a fixed exterior probe set. Species sort by \((N_c,\, e_{\mathrm{proxy}},\, \chi_{\mathrm{net}})\).

This does **not** claim the correct particle spectrum of nature. It claims that under a pure integer-lock rewrite ontology — no lattice — discreteness and isolation are generic, not tuned. Changing the resonance polynomial changes the census; finiteness at fixed mode bound remains.

### 9.2 Combinatorial Gauss result (ran 2026-07-27)

```
./construct_gauss
# fit: (K-K0) = kappa*M + b  with kappa≈-0.874  b≈-2.82  r²≈0.90
```

On a triangular-disk abstract graph (\(d_*=6\), 91 vertices), converting \(M = 1\ldots 12\) interior parcels to dense (remove + reseal boundary) yields **approximately linear** total degree-defect change versus \(M\). Sign of \(\kappa\) is convention/boundary-dependent on a finite patch; the structural claim is linearity \(K \propto M\), foam-dependent slope. Coordinates were used only once to *generate* a neighbour list; the measured quantities are purely combinatorial on that list.

---

## 10. What was tried and what still sits open

### Tried while constructing

1. **Parcel-free pure global mode ledger** — only totals \(E_S, E_D, E_F\).
   Failed: cannot express local cycle gates or spectral boundaries; conversion
   became ambient accounting again.
2. **Fixed foam with dynamical occupation only** — closest to DESIGN_CLEAN.
   Failed the operational background test (immortal cells + fixed
   connectivity).
3. **Repulsion from space scarcity** (§4.6 alternative) — inconsistent with
   gravity-as-conservation; replaced by steepest closure on defect.
4. **Three modes only** — inconsistent with mid-cycle energy under non-instant
   transfer; added \(T\) as derived, not imported species of particle.

### Still open (not failures of principle; unfinished construction)

* Exact integer resonance polynomial (the concrete form of \(h(v)\) matching)
  is postulated schematically; a unique algebraic choice (octonionic 7→3 or
  other) is not derived here.
* The constant \(\kappa/\varepsilon_0\) is not computed from a first-principles
  channel algebra — only the form \(K \propto M\) is constructed.
* Confluence of the rewrite system (does steepest closure always terminate
  locally?) is not proved; branching is allowed but not classified.
* Binding channels (§5.7) need explicit joint-lock examples from the
  enumerator before mass ratios can be claimed.
* Quantum Ring Theory was not examined (pointer only, per brief).

### Load-bearing postulates (audit list)

These are the joints. If the construction fails, it fails here:

1. State = finite energy complex with ephemeral parcels/channels.
2. Dynamics = rewrite relation with \(E\)-conservation + cycle gate.
3. Steepest closure on open-cycle deficit \(\Delta\).
4. Saturation ⇒ \(S\to D\); exterior complete cycle ⇒ possible unlock.
5. Modes exhausted by conversion topology \(\{S,D,F,T\}\).
6. Combinatorial Gauss law \(K(I)=\kappa M(I)/\varepsilon_0\).
7. Free fall = steepest-closure bias on defect-laden \(C_S\).
8. Particle = self-reproducing integer lock with spectral boundary.

Everything tagged **[D]** stands if the principle stands. Everything tagged
**[P]** is available for revision without abandoning the principle — but
revising it requires a replacement that still fills the gap, not an empty
objection.

---

## 11. Closing statement

The universe constructed here has one law: **energy is never destroyed; it only
changes mode through complete chiral cycles.** The state is the current finite
pattern of that energy. The dynamics is which patterns may succeed which.
Space is what energy looks like when it is not locked and not mid-transfer;
mass is locked closed conversion; light is open conversion; gravity is the
bookkeeping of remaining space when some has been converted; particles are
recurrent locks; \(c\) is how fast open conversion can advance; charge is
structure; momentum is a ghost of a stage that is not there.

No background is used. No field is imported to carry an effect. Where the
construction is unfinished, the unfinished piece is named and the postulate
that currently holds the joint is labeled for audit.
