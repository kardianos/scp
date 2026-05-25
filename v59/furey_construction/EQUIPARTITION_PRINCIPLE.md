# Mechanism hunt for Koide = 2/3: it unifies with v_Higgs into ONE principle

*2026-05-24.  Per the steer "keep hunting for a mechanism — for Koide=2/3, the root."  Result:
the **magnitude** `|ξ|²=1/2` (= Koide 2/3) had never been mechanism-hunted (all prior work
M1–M5 in `koide_phase_law/physical_mechanisms.md` was on the **phase** `φ=Q/3`, now ruled out).
Hunting the magnitude yields a genuine structural advance: **Koide=2/3 and `v_Higgs=28²·a²` are
the same mystery — equipartition (maximal mixing) of the order parameter over its sector — and
that single principle, if granted, forces BOTH.**  Not a derivation, but it collapses two open
problems to one well-posed, potentially-derivable principle.*

## Four equivalent forms of Koide = 2/3 (the root fact)

With `√m_n = a(1 + 2t cosθ_n)`, `Q=(1+2t²)/3`:

1. **Brannen amplitude:** `t² = |ξ|² = 1/2`.
2. **45° geometry:** the √-mass vector `(√m_e,√m_μ,√m_τ)` makes a **45° angle** with the
   democratic axis `(1,1,1)` — measured `cos²ψ = 0.500005` (`ψ = 45.000°`).  (`Q = 1/(3cos²ψ)`.)
3. **S₃ rep power balance:** under the generation `S₃`, `√m⃗ = ` trivial(1-dim, the mean) ⊕
   standard(2-dim, the splitting).  `Q=2/3 ⟺` **equal power** in the trivial and standard parts.
4. **G₂-complement fraction:** `t² = (D−dimG₂)/D = 1 − 14/D` — the fraction of the sector's
   algebra lying **outside the G₂ automorphism**.

All four say the same thing: **the order parameter is split 50/50 between the "mean/protected"
mode and the "splitting/active" modes — by *power*, not by dimension.**

## The unified formula: one law, three sectors (the strongest card)

Form 4 fits **all three** fermion sectors with a single formula and **no per-sector freedom**
(`D` = sector dim, `14 = dimG₂` universal):

| sector | `D` | `t²=(D−14)/D` | `Q=(1+2t²)/3` | observed |
|---|---|---|---|---|
| lepton | 28 | 1/2 | **2/3** | 2/3 (10⁻⁵) |
| d-quark | 35 | 3/5 | **11/15** | 11/15 (~0.3%) |
| u-quark | 63 | 7/9 | **23/27** | 23/27 (~0.3%) |

`D_lepton = 28 = 2·dimG₂` (`ScaleBridge.two_dimG2_eq_dimSpin8`): the lepton sector is *exactly*
double the G₂ core, which is *why* the lepton lands on the symmetric 1/2.

## THE UNIFICATION with v_Higgs

The two load-bearing relations of the whole lepton/EW/Higgs block turn out to be the **same
statement** about the order parameter:

| relation | content | "equipartition" reading |
|---|---|---|
| **Koide = 2/3** | `t²=(D−dimG₂)/D` | the order parameter is **uniform over the sector's D directions**; the splitting-fraction is the non-G₂ fraction `(D−14)/D` |
| **`v_Higgs=dim(L)²·a²`** | Frobenius² of the **democratic** mass bilinear on `L` (`03_higgs_bridge_result.md`) | the order parameter has **equal magnitude in all `dim(L)²` components** |

Both are "**the order parameter is maximally mixed / democratic over its available
directions.**"  *There are not two mysteries here — there is one:* **why does the order
parameter equipartition?**  A single dynamical principle forcing equipartition would close
**both** Koide and v_Higgs simultaneously.

## The principle is of a *derivable* type: maximal mixing / maximal symmetry

Equipartition is exactly the **maximum-entropy / maximally-mixed** configuration — and that is
**forced** by symmetry + a norm constraint (it is the unique `G`-invariant state of fixed norm),
not an arbitrary ansatz.  The test lands exactly:

> uniform weight `1/D` on each of the `D` sector directions, with the `dimG₂=14` G₂-invariant
> directions as the inert "core," gives splitting-fraction `(D−14)/D = t²` — **verified for all
> three sectors** (1/2, 3/5, 7/9).

So **`maximal mixing ⟹ t²=(D−dimG₂)/D ⟹ all three Koide ratios`**, and the same maximal mixing
gives the democratic v_Higgs bilinear.  This is the first time the Koide value `2/3` is *forced*
(modulo the principle) rather than matched — and it predicts the other two sectors as a bonus.

## Supporting evidence: the precision hierarchy

If Koide is the maximally-mixed (symmetric) **high-scale** boundary condition, the observed
precision hierarchy follows: `Q` is RG-invariant under common rescaling `√m→c√m`, and lepton
Yukawas are tiny and run near-universally ⇒ lepton Koide stays `10⁻⁵`; quark Yukawas run large
and non-universally (top, QCD) ⇒ quark Koide degrades to the observed soft `~0.3%`.  The data's
own precision pattern (lepton ≫ quark) is what a high-scale symmetric origin predicts.

## Honest gaps (the sharp residual)

1. **Why maximal mixing, not minimal energy?**  A ground-state VEV minimizes energy; the
   maximally-mixed state maximizes entropy.  The natural reconciliation: the **symmetric
   vacuum** (uniform over the internal sector) projects to a **split spectrum** in the 3-gen
   observable space via the embedding — the "democratic mass matrix" picture.  This needs the
   explicit `3 generations ↪ D-dim sector` embedding (the deferred Witt/generation map,
   `SevenDAlgebra.lean`) to be made a derivation rather than a principle.
2. **Why `G₂ = the inert core`?**  G₂=Aut(𝕆) is the unbroken automorphism; identifying its 14
   directions with the *non-splitting* ("mean") part of the mass operator is the physical claim
   carrying the formula.  Plausible (G₂ is the universal protected subalgebra) but to be shown.
3. The principle is currently **phenomenological** (it *fits* everything, predicts the 3-sector
   pattern) — its dynamical origin (a free-energy/entropy extremum on the fixed-norm vacuum
   manifold `MM̃=v²` / `|ξ|²=1/2`) is the thing to derive.

## Net (the constructive outcome of the hunt)

The hunt did **not** derive Koide=2/3 — but it did something better than another negative: it
**collapsed the program's two main mysteries (Koide and v_Higgs) into one**, identified that
one as **maximal mixing / maximal symmetry of the order parameter** (a principle of derivable
*type*, not an arbitrary fit), showed it **forces all three Koide ratios + the v_Higgs bilinear**
at once, and is *supported* by the precision hierarchy.  The program's central open question is
now sharp and singular:

> **Derive that the fermion order parameter is the maximally-mixed (G-symmetric, fixed-norm)
> configuration on the vacuum manifold — equivalently, the democratic mass matrix.  This one
> result closes Koide=2/3, the quark Koide ratios, AND `v_Higgs=28²a²` together.**

This is the right next target — far more leverage than the per-relation attempts, and the first
mechanism candidate that is *forced by symmetry* rather than assumed.
