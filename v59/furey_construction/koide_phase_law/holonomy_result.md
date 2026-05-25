# Result: the inter-generational holonomy test — the geometric mechanism is ruled out

*Built the inter-generational connection and ran `arg(ξ³) = Q`.  Verdict: **negative** —
the geometric/holonomy mechanism (M2/M3, the former lead) cannot derive `φ = Q/3`.*
See `intergenerational` computation in this note and `witt_map.py`; machine-checked pieces
in `lean/PhaseExclusions.lean`, `lean/ColorSU3.lean`.

## The connection that was built

The genuinely-motivated inter-generational connection is the **Brannen 3-site ring**: three
generations on a ring with complex hopping `ξ = t e^{iφ}` (the off-diagonal of
`M = a(I + ξS + ξ̄S²)`).  Its **Wilson loop** around the closed generation cycle `0→1→2→0`
is `W = ξ³`, with `arg(W) = 3φ` — a **Peierls / Aharonov-Bohm flux** through the generation
ring.  So `arg(ξ³) = 3φ` holds *by construction*: the holonomy **is** the phase.

The real test is therefore not "is `arg(ξ³)=3φ`" (tautological) but **"is that flux *forced*
to `Q` by a structural curvature?"**

## The test: no structural curvature forces the flux to Q

A flux equals `∮A = ∫∫F`; to force `3φ = Q` needs a curvature `F` integrating to `Q`.  Every
structural candidate fails:

| candidate curvature source | holonomy | = Q = 2/3? |
|---|---|---|
| v58 connection `ω ∝ ∇log(M M̃)` (exact 1-form) | `0` | no (gradient ⇒ zero loop holonomy) |
| discrete / center `Z₃` flux | `{0, 2π/3, 4π/3}` | no |
| octonion associator `[eₐ,e_b,e_c]` | integer coeff `{−2,0,2}` ⇒ `arg ∈ {0,π}` | no |
| color `Z₃` on the lepton (machine-checked) | `0` (fixes the singlet) | no |

**The clincher (machine-checked, `PhaseExclusions.koide_not_pi_rational`):** a
structural/geometric holonomy — from a discrete group, a root system, or the octonion
associator — is a **rational multiple of π**.  But `Q = 2/3 rad` is **not** π-rational
(`2/3 = π·(2/3π)`, `2/3π` irrational, via `irrational_pi`).  So **no structural holonomy can
equal `Q`.**  The flux `3φ` is a *free Aharonov-Bohm phase*, unforced by the algebra
(consistent with `BrannenPhase.Q_phase_independent`).

## Verdict — the leading mechanism is dead

The geometric-phase mechanism (M2/M3 in `physical_mechanisms.md`, previously rated "most
promising") is **ruled out**:
1. its connection (the v58 `ω`) is exact ⇒ zero holonomy;
2. any structural holonomy is π-rational, but the required value `Q` is not π-rational;
3. the color `Z₃` fixes the lepton singlet (so even the intra-ideal version is trivial).

So `φ = Q/3`, though empirically Koide-tight (`CONCEPT.md`), is **not a geometric/holonomy
phenomenon.**  Building the connection — as requested — closed the route rather than opening it.

## What survives, and where the search goes now

With the geometric route excluded, the surviving options for `φ = Q/3` are:
- **non-geometric algebraic/spectral origin** — a relation among the kernel's *real* spectral
  data (e.g. a second symmetric-function invariant of the masses that happens to fix `cos 3φ`
  to `cos(2/3)`), with no angle/holonomy interpretation.  This is the only structural avenue
  left; it must produce the transcendental `cos(2/3)` (`PhaseAmbiguity.invariant_is_cos_two_thirds`),
  not a "nice" number, so it faces the same hurdle.
- **genuine numerical coincidence** at `m_τ` precision — `φ = Q/3` is a real `~10⁻⁵` regularity
  (as tight as Koide) but with no mechanism, like Koide itself for decades.

**Honest bottom line:** the holonomy test was run and answered — `arg(ξ³) = Q` is **not
forced** (the flux is free; structural holonomies are π-rational, `Q` is not).  The most
attractive physical mechanism is excluded, leaving `φ = Q/3` as a precise but mechanism-less
empirical law whose only structural target (the transcendental `cos 2/3`) matches no candidate.
The negative result is itself the deliverable: it removes the geometric hypothesis from the
design space for good.
