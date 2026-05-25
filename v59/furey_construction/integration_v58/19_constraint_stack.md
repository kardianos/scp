# Stacking the constraints on the OFE: the right idea — and it reduces the lepton sector to ONE number

*2026-05-24.  Per the steer "what if we added our other constraints to the OFE."  The framing is
right — the goal is the **OFE plus all constraints as one well-posed problem**.  But the decisive
move is to **classify** the constraints: most of ours are **Koide in disguise** (circular).  Once
the circular ones are removed, stacking the genuinely-independent constraints **collapses the
lepton sector to a single free dimensionless number `t²`** — with the phase `φ` *following* from the
phase law.  The sole remaining job is for the OFE dynamics to fix `t²=½`.*

## The classification (the thing that was being missed)

The 2 dimensionless unknowns are `t²` and `φ`.  Sorting every constraint we have:

| constraint | effect on `(t²,φ)` | class |
|---|---|---|
| `Q=(1+2t²)/3` | defines `Q` | **proven** (not a constraint on `t²`) |
| **`t²=½` (Koide)** | `t²=½` | **TARGET** |
| **grade-balance** `|diag|=|off|` | `t²=½` | **= target (circular)** |
| **`|ξ|²` = L-grade complex norm** | `t²=½` | **= target (circular)** |
| phase law `φ=Q/3` | `φ=(1+2t²)/9` | **real relation** (empirical 10⁻⁵) |
| vacuum `MM̃=v²` | fixes scale `a` | independent (no `t²,φ`) |
| mass ∈ L (`u²=−1`) | which grade | proven (no value) |
| S₃ circulant | `M=a(I+ξS+ξ̄S²)` | proven form (leaves `t²,φ` free) |
| chiral: electron light | `φ < φ_chiral` | bound only |
| **OFE F-channel dynamics** | extremize `‖⟨Ω²⟩₄‖²+…` | **INDEPENDENT (dynamical)** |

**The trap:** grade-balance, the L-grade complex norm, and `t²=½` are the **same statement** — the
target.  Stacking *those* "constraints" proves nothing (it assumes the answer).  This is what makes
"adding constraints" feel like it should work yet not close: three of our strongest-looking
constraints are circular.

## What stacking the INDEPENDENT constraints actually buys

Remove the circular ones and impose the rest.  Then:

- **the phase is no longer independent:** `φ = Q/3 = (1+2t²)/9` — at `t²=½` this gives `φ=2/9` ✓.
  So **fixing `t²` fixes `φ` for free** (via the empirical phase law, one of our structural
  relations).
- vacuum fixes the scale `a`; `mass∈L` + S₃ fix the *form*; chiral only *bounds* `φ`.

> **Net: stacking collapses the lepton sector to a SINGLE free dimensionless number, `t²`.**  (Down
> from `t²,φ`.)  Everything else is either circular, fixes the scale, fixes the form, or follows
> from `t²` via the phase law.

## The bottleneck — and why it's sharp now

**No constraint we have is both (a) independent of Koide and (b) fixes `t²=½` — except the OFE
dynamics.**  And the OFE's coupling-selecting term, the F-channel `‖⟨Ω²⟩₄‖²` alone, scales as `t⁴`
(monotonic, minimised at `t²=0`) — so it fixes `t²` only when **balanced against** the norm/kinetic
terms (`MM̃=v²`, `|DΩ|²`, the `λ,μ`).  The interior extremum exists, but **its location depends on
the OFE coefficients**, which are not yet pinned.

So the whole lepton-mass problem reduces to one sharp, well-posed question:

> **Does the constrained OFE — the grade-`{0,2,4,6}` energy with the F-channel, on the `mass∈L`,
> `MM̃=v²`, S₃-circulant vacuum — have its extremum at `t²=½`?**  If yes, the phase law delivers
> `φ=2/9`, and the entire lepton sector is fixed (one dimensionful scale `a`, nothing else).

## Honest status & next step

- **Gained (this analysis):** the constraint stack **reduces 2 unknowns to 1** (`t²`; `φ` follows),
  and **isolates the circular constraints** (grade-balance, L-norm — don't count them).  The problem
  is now *"fix one number with the OFE."*
- **Needed:** (i) pin the OFE coefficients (`λ,μ`, the F-channel weight) — from v58's safe band, or
  a normalization principle; (ii) compute the constrained extremum in `t²`; (iii) the generation↔
  internal map (`14_…md`) to read the extremum as the Brannen `t²`.  Then it's a clean yes/no on
  `t²=½`.

**Your instinct was right and sharpened the problem:** adding the constraints *does* help — it
removes the phase as an independent unknown and leaves a single number — but only after discarding
the circular ones, and the last number still needs the OFE dynamics (with pinned coefficients) to
fix.  That is the cleanest the lepton-mass problem has been posed.
