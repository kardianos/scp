# VERDICT — should the simulation be queued?

**Answer: YES, and the first three items cost nothing.**

Queue P1, P7 and P4 (`SIM_PREDICTIONS.md`) immediately — they are offline
analyses of frames that already exist. Queue P5 and P6 next (config only).
P2 is the only item that needs a kernel change and therefore user
authorization; it is also one of the two decisive discriminators, so it should
be *asked for*, not slipped in.

**Nothing is proposed for adoption.** No law-table change, no bar retired, no
knob defaulted differently.

---

## Against the brief's success criteria

| # | criterion | met? | evidence |
|---|---|---|---|
| 1 | one law table with `I` and `O` dual in the same language | **yes** | `UNIFIED_LAW.md` §2 (L1–L5) and §3. One coupling `A` drives both directions; `fullratsimp(O − I) ≡ 0` symbolically, `1e−12` numerically over a `(T, mu, eps)` grid. |
| 2 | shells as closed harmonics with discrete capacity | **yes** | §5. Shells are Laplacian degeneracy classes; capacity = eigenspace dimension. Verified: hypercube `C(d,k)` for `d=2..5`; icosahedron `1,3,5,3`; ring `1,2,2,...`. |
| 3 | charge as incomplete-shell defect + conservation statement | **yes** | §6. `Q = Z − N`; Theorem 6.1 (conserved under every promotion, emission, absorption; changes only by transfer or graph surgery); Theorem 6.2 (a full shell's circulation cancels identically); Theorem 6.3 (a forest cannot be charged). |
| 4 | ≥2 sharp predictions distinguishing this from "radiance only + noise fuel" | **yes** | P2 (`x*` tracks the bath with slope 1 — the incumbent predicts slope 0) and P5 (`x*(coherent)/x*(random) ≥ 1.35` — the incumbent predicts exactly 1.00). P1 and P6 are two more. |
| 5 | `go test ./...` green | **yes** | 35 assertions pass; `gofmt` clean; `go vet` clean. |

## Against the brief's NO conditions

| NO condition | status |
|---|---|
| "dual absorption collapses to handwaved `T_bath^4`" | **avoided.** Absorption is not postulated at all. It is the *same operation* as emission (L3) with the *same coupling* (L4); the `T_bath` dependence is derived, and the Boltzmann factor is forced rather than inserted (`detailed_balance.mac`). |
| "charge is just renamed `Em`" | **avoided.** `Z` counts eigenspace dimensions — a symmetry fact about the contact graph, independent of energy content. Two objects with equal `Em` and different graphs carry different `Q`; two with the same graph and different `Em` at fixed `N` carry the same `Q`. And Theorem 6.3 gives `Q` a property `Em` cannot have: it is impossible on a tree. |
| "shells require imported orbitals" | **avoided.** The icosahedral contact graph yields capacities `1, 3, 5` and magic numbers `1, 4, 9` — the `2l+1` ladder and the `n²` filling — from its automorphism group alone. No central potential, no angular momentum operator, no orbital. |

---

## What actually earned the YES

Three things, in order of weight.

**1. It retrodicts four measured nulls from two theorems.** This was not
designed in; it fell out.

- Box-fill (`MECHANISM.md` S07–S19) is Theorem 3.4: with a random-phase bath
  the fixed point is `x* = xb` *exactly*, so a hot bath can only fill the box to
  its own load. That was never a tuning failure.
- Door scrambling (Face B) is L5: a point-contact door is precisely the
  operator with flat support over every winding transfer — mean 0, sd 3.49.
- Dark-bath silence (Face B, G03) is Corollary 3.2: an aufbau ground state has
  every downward promotion blocked by cap refusal.
- Face C's broad hoard spectrum is the `Delta/T <~ 1` branch of one
  distribution (§7).

**2. It converts an adopted constant into a measurement.** `p_rad = 4` stops
being a number and becomes the claim `d_s = 3` about the foam's own contact
graph — testable offline, today, with an estimator calibrated to `±0.02`
against exact tori.

**3. It names the missing knob.** The programme's standing gap is "phase-safe
self-restoring intake" (`METABOLIC_MATTER.md` §5.2). This law says the variable
is the **coherent fraction `g`**, not `noise_amp`, and gives the curve:
`u* = 1, 1.380, 2.107, 3.036` at `g = 0, 1, 2, 3`. `g = 0` is provably a filled
box; `g > 0` provably holds matter above ambient, with `Phi* > 0` throughout —
a live cycle, not a sealed bag.

---

## What the YES does not cover

These are real and are stated in `UNIFIED_LAW.md` §10. Repeated here so the
verdict is not read as broader than it is.

1. **There is no force law for `Q`.** Charge is conserved, quantized and
   signed, but nothing here derives attraction between opposite defects. Until
   that closes, this is a quantum number, not electromagnetism — and the
   programme's charge lane is *not* finished by this note.
2. **The graph is given, not selected.** Monodispersity is conditional on which
   contact graph a hoard adopts; the theory says the magic numbers are partial
   capacity sums, not which geometry forms. This is the same wall
   `METABOLIC_MATTER.md` §3.2 names ("geometry is not quantized"), relocated
   rather than removed. **It is the largest remaining gap.**
3. **`d_s = 3` is assumed, not measured.** Until P1 runs, `p_rad = 4` has been
   *explained*, not *confirmed*.
4. **No dynamics.** This is a law of states, rates and invariants. It says
   nothing about transport, mobility, or the binding/mobility antagonism in
   `RING_DNLS.md`.
5. **`Nu = 1` and the binding threshold are imported** from the adopted table
   and the `Es` well, not derived.

A fair summary: the note unifies **bookkeeping** — what states exist, what
transitions are allowed, what is conserved, what sets a scale. It does not
unify **dynamics**, and it does not close the charge lane.

---

## Recommended next action

1. **Run P1 offline now.** It is the cheapest item in the list and the only one
   that can invalidate an adopted constant.
2. **Run P7 and P4 offline** against existing Face C / RING frames, with the
   shuffled-graph control pre-registered.
3. **Then decide on P5/P6** (config only, no kernel touch).
4. **Ask the user about P2** — it adds a term to the door and is a kernel
   change. It is the sharpest test in the set, and it should be requested
   explicitly rather than folded into an analysis pass.

If P1 returns `d_s ≈ 3` and P5 returns a ratio near 1.00, the interesting half
of this law is dead and the incumbent picture stands. That is the outcome the
predictions are shaped to produce quickly if it is the true one.
