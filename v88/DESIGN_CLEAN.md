# The cycle-trapped fabric — design

> **SUBORDINATE TO `v88/PRINCIPLE.md`.** That document is foundational:
> energy is never destroyed, space is a *mode* of energy, matter is
> converted space, and curvature follows from conservation. Anything here
> that assumes a background — a permanent index set with evolving contents
> — is superseded by it, including every instrument built so far.


A statement of the model as it stands. No history, no superseded alternatives,
no record of how it was arrived at.

---

## 1. Substance

**Space is a discrete cell complex.** The cells are the substance, not a
numerical sampling of something continuous. The cell spacing `a` is a physical
constant of the theory, not a resolution parameter. There is no continuum limit
in which the theory is supposed to be recovered; the discreteness is the physics.

**Cells have variable size.** Cell size is a dynamical quantity, not fixed by the
tessellation. Writing the mapping from logical index to physical position as

```
    x(i) = a·i + ξ(i)              ξ is the fabric's own degree of freedom
    θ(i) = (∇·ξ)(i)                local density;  θ < 0 = shrunk, densified
```

the displacement that *defines* the tessellation is the same object that carries
mass. Nothing lives "on" the fabric; the fabric's geometry is the only field.

**Cells are internally multi-dimensional.** A cell is not a point holding a
number. It carries an internal configuration space, compact, whose free motion
has a discrete harmonic spectrum. Mode indices are integers by the topology of
that internal space, not by tuning or truncation.

---

## 2. Links and the transfer of energy

**Cells do not exchange energy directly.** Every link between adjacent cells
carries its own harmonic. Energy is deposited into the link, carried through the
link's cycle, and only then delivered to the far cell.

**Transfer therefore takes time.** It is not an instantaneous term evaluated at a
point in time; it is a process with a duration set by how fast the link's cycle
can complete.

**Energy crosses only in complete cycles.** A partial cycle delivers nothing. A
link whose harmonic cannot be driven resonantly by the cell modes never accepts
the energy in the first place — the gate is structural, a property of what the
link can be excited into, not a suppression applied after the fact.

**A link actuates in two planes, along a chiral harmonic.** A link along
direction `d̂` has a plane transverse to it, and the actuation occupies both
transverse directions with a relative phase. Writing the two transverse
amplitudes as `u₁, u₂`:

```
    u± = (u₁ ∓ i u₂)/√2                left / right circular
    χ  = 2 Im(u₁* u₂)                  signed handedness
```

Handedness is a property carried by the geometry of the link, not an added
label. Same-handed and opposite-handed configurations couple through different
channels.

**c is the rate at which cycles couple.** The speed limit on the propagation of
energy is not postulated. It is the fastest that a link harmonic can complete a
cycle and hand energy across, so

```
    c ~ a × (rate of cycle coupling)
```

is a derived constant of the fabric.

---

## 3. Dimension

Handedness requires a cross product. The cross product exists only in **3 and 7
dimensions**. The fabric carries **7 degrees of freedom**, from which **3
dimensions emerge** — the three that the chiral link structure realises.

The dimension count is therefore structural, not a parameter. It is not to be
scanned over or tuned.

---

## 4. Cyclic energy and cell size

**Cyclic energy tightens a cell.** This is opposite to a fluid, in which internal
energy expands. A cell holding more cyclic energy runs at a shifted rate and
occupies less volume.

**The tightening has a limit.** Beyond some cyclic energy the cell cannot tighten
further. The energy must then be released into *configuration* — into the
arrangement of the cell's internal harmonics and of its neighbours — rather than
into further compression.

**Densification is globally scarce.** On a closed fabric the total dilation is
pinned: `Σᵢ θᵢ = 0` identically, a geometric fact rather than an imposed
constraint. Not every cell can be shrunk. Any densified region must be paid for
by rarefaction elsewhere, so the ground state cannot be uniform.

---

## 5. What a particle is

A particle is a **region whose cycle structure cannot complete a cycle with its
surroundings**. Energy inside is reflected internally at the boundary. The object
is held by that spectral mismatch, not by a balance of forces.

It has three parts, all necessary:

1. **An interior configuration** — a locked arrangement of the cells' internal
   harmonics, holding the dense core against its own compression. This
   arrangement is what plays the stabilising role; it is a *configuration* of
   degrees of freedom the fabric already has. It is emergent and
   intra-particle, and it is never a field.
2. **A sharp exterior distortion** — the density and rate jump at the boundary
   that maintains the mismatch.
3. **Internal reflection** — the consequence: no complete cycle connects
   interior to exterior, so energy cannot leave.

**Mass is a cell count.** A particle occupies an integer number of cells.

**Existence is a threshold, not a dial.** A configuration either closes its
cycles — internally, and against the exterior — or it does not exist. There is
no continuous family of particles.

**Species are configurations.** Different particles are different internal
harmonic arrangements satisfying closure, distinguished by integer labels: mode
indices, resonance vectors, handedness, and cell count. Mass ratios follow from
those configurations.

---

## 6. Constraints on any implementation

1. **Emergent from the space-time fabric alone.** Particles, their species,
   their masses and their interactions must arise from the fabric's own degrees
   of freedom and its own dynamics. Nothing may be placed by hand to stand in
   for a physical effect, and no structure may be assumed that the fabric does
   not itself produce.

2. **No imported field or species.** No additional field may be introduced to
   carry an effect, and no particle species may be put in by hand. In particular
   the stabilising interior arrangement of §5.1 is a *configuration* and must
   never be represented as a field with its own potential. If an effect requires
   a new field to appear, that is a failure of the design, not a licence to add
   one.

3. **No reduced-dimensionality substitutes.** The structure is three-dimensional
   and chiral. A one- or two-dimensional reduction removes the transverse plane
   and therefore removes handedness, which is a load-bearing part of the model.
   Reduced-dimension tests do not probe this design.

4. **No linearised proxies.** The mechanism is nonlinear — tightening,
   saturation, release into configuration, and trapping are all nonlinear. A
   linear model can exhibit propagation and gating but cannot bind, and must not
   be taken as evidence about binding.

---

## 7. What the model must produce to be right

* Localised objects that are **multi-cell**, with a definite interior and a
  sharp exterior.
* A **discrete** spectrum of such objects — distinct species with integer
  labels, recurring across independent preparations rather than varying
  continuously with initial conditions.
* **Binding**: an interaction between objects that vanishes with separation and
  is attractive in some channel, so composites can form.
* **Mass ratios** fixed by configuration, not by a continuous parameter.
* A derived **c**, and a transfer that is genuinely gated by cycle closure.

Failure of any of these is a failure of the design, not of a parameter choice.
