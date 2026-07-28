# v87 Rung B0-G — the geometric hypothesis
**Supplement to `BELL_BRIEF.md`. Read that first; everything in it still binds.
This file adds ONE specific hypothesis for the GEOM seat to pursue.**

---

## 1. The hypothesis, stated

> Locality produces a **hard stop**. The geometric bound carried by **c²** — the
> quadratic form of the metric, the light cone — is **the same relationship**
> that produces the asymptotic stop at **2√2** in CHSH. What *enables* exceeding
> the classical bound of 2 in the first place is **particle/energy-as-fabric**
> (O2), which permits **recursive non-locality by geometric construction**.

So the hypothesis assigns the two halves of the problem to two different
structures:

| half | assigned to |
|---|---|
| **Reachability** — how does the model get *above* 2 at all? | field monism (O2): "recursive non-locality by geometric construction" |
| **The bound** — why does it stop at 2√2 and not run to 4? | the metric quadratic form; the same geometry that gives c² |

**Keep these two halves separate at every step.** The most likely failure mode of
this rung is to produce an elegant geometric account of the *bound* and quietly
assume the *reachability*. A geometric derivation of 2√2 that does not also say
which Bell assumption is relaxed has solved the second half only, and must say
so in those words.

---

## 2. Why the "bound" half is promising — the geometry is real and elementary

There is a genuinely geometric derivation of Tsirelson, and it is worth taking
as the starting point because it involves nothing but a real inner-product
space.

For unit vectors **a**, **a′** in ℝⁿ, the parallelogram law gives
‖a+a′‖² + ‖a−a′‖² = 2(‖a‖² + ‖a′‖²) = 4, and Cauchy–Schwarz then gives

    ‖a + a′‖ + ‖a − a′‖ ≤ √(2 · 4) = 2√2,

with equality iff **a ⊥ a′**. If the correlation function takes the inner-product
form E(a,b) = −a·b, the CHSH combination is
S = −(a+a′)·b − (a−a′)·b′, and |S| ≤ ‖a+a′‖ + ‖a−a′‖ ≤ 2√2 immediately.

**So 2√2 is a fact about the 2-norm on a real inner-product space** — a
parallelogram identity plus Cauchy–Schwarz. It is the same class of object as a
metric quadratic form. That is the concrete sense in which the hypothesis could
be right, and it is the thing to make precise.

The operator route says the same thing differently: S² = 4·1 + [A,A′]⊗[B,B′]
with ‖[A,A′]‖ ≤ 2, hence ‖S‖ ≤ √8 = 2√2.

**The real question this rung must answer is therefore not "why 2√2 given inner
products" — that is elementary — but: WHY ARE THE CORRELATIONS INNER PRODUCTS OF
UNIT VECTORS?** If fabric geometry forces E(a,b) to be a bilinear form of unit
vectors in a space with a definite (or Lorentzian) quadratic form, the bound
follows for free. Deriving that form from fabric structure *is* the deliverable.

---

## 3. Two hard constraints. Ignoring either invalidates the result.

**Constraint 1 — the light cone alone does NOT give 2√2.** Popescu–Rohrlich
boxes are perfectly no-signalling and reach S = 4. Relativistic causality is
therefore *insufficient*; "locality produces a hard stop" cannot be
"no-signalling produces a hard stop," because that stop is at 4, not 2√2. The
hypothesis survives only in the stronger form: it is not causality *per se* but
the **specific quadratic form / inner-product structure** — the same algebraic
object that appears in s² = c²t² − x² — that supplies the 2-norm bound. Make
that distinction explicitly and early, or the whole construction is void.

**Constraint 2 — field monism alone does NOT escape Bell.** The tempting
argument is: "the two particles are one extended field structure, never two
objects, so factorisability was never justified." This fails as stated, and it
is important to understand exactly why. Given local dynamics, one may always
take λ to be **the entire field configuration on a Cauchy slice** in the common
past. The responses A(a, λ), B(b, λ) are then local functions of it, (R) and (L)
both hold, and CHSH ≤ 2 follows regardless of how holistic the field is. Ordinary
classical field theory is exactly this case.

Therefore **"recursive non-locality by geometric construction" must do something
sharper than assert holism.** The specific thing it must break is the
*Cauchy-slice specifiability of λ*: the fabric state must not be freely
specifiable on a spacelike surface. Candidate mechanisms, all of which are
legitimate and none of which is free:

- a **global self-consistency constraint** on histories (an all-at-once /
  boundary-value structure rather than an initial-value one) — this is Door C
  territory and is the most natural reading of "recursive";
- a **topological or holonomic constraint** linking regions, so that slice data
  is over-determined;
- **measurement dependence** arising because the chooser's fabric state is not
  independent of the source's (Door A), which the geometry might *derive* rather
  than postulate.

State plainly which one the construction uses. "Recursive" is currently a word,
not a mechanism; the rung's job is to turn it into one or to show it cannot be.

---

## 4. Anchors worth working through

Verify each before relying on it; mark **UNVERIFIED** if you cannot confirm it.

- **Cirel'son (Tsirelson)**, Lett. Math. Phys. **4**, 93 (1980) — the bound.
- **Landau**, Phys. Lett. A **120**, 54 (1987) — an elementary derivation.
- **Popescu & Rohrlich**, Found. Phys. **24**, 379 (1994) — no-signalling reaches 4.
- **Cabello, Severini, Winter**, PRL **112**, 040401 (2014) — graph-theoretic
  approach; the quantum bound is the **Lovász ϑ function** of the exclusivity
  graph. *This is the sharpest existing example of "geometry produces a hard
  stop": ϑ is defined by an orthonormal representation in a real vector space.*
  If the fabric hypothesis is right, it should connect here.
- **Pawłowski et al.**, Nature **461**, 1101 (2009) — information causality.
- **Navascués & Wunderlich**, Proc. R. Soc. A **466**, 881 (2010) — macroscopic
  locality: the quantum set is what remains when the coarse-grained limit must
  look classical. Structurally close to a "geometric hard stop from a limit."
- **Fritz et al.**, Nat. Commun. **4**, 2263 (2013) — local orthogonality.
- **Navascués, Guryanova, Hoban, Acín**, Nat. Commun. **6**, 6288 (2015) —
  "almost quantum" correlations. **Read this as a caution:** there is a set
  strictly larger than quantum satisfying all the known physical principles, so
  no principle yet derived picks out quantum uniquely. A construction landing at
  "almost quantum" rather than exactly 2√2 is a *known* outcome, not a failure —
  but it must be identified as such rather than rounded to success.
- **Barrett**, PRA **75**, 032304 (2007) — generalized probabilistic theories;
  the state-space geometry that fixes the correlation set.
- **Howard**, Stud. Hist. Phil. Sci. **16**, 171 (1985) — Einstein's separability
  (*Trennungsprinzip*) as distinct from locality. Relevant to §3 Constraint 2,
  and to stating precisely what field monism does and does not buy.

---

## 5. What would count as success, in descending order

1. **Full**: fabric geometry forces E(a,b) into a unit-vector bilinear form,
   the metric quadratic form yields max|S| = 2√2 as a theorem, AND the
   reachability mechanism is identified and derived (not postulated), with the
   cost quantified.
2. **Strong**: the bound derived geometrically from the c²-form, with the
   reachability mechanism identified but postulated. State the gap.
3. **Partial**: a demonstration that the quadratic form gives *a* bound, but at
   the wrong value (e.g. "almost quantum", or 4, or below 2). Report the value.
4. **Negative, and still valuable**: a proof that the c² quadratic form cannot
   by itself produce 2√2, with the obstruction identified. Given Constraint 1,
   this is a live possibility and would be a real result.

Do the algebra in Maxima or SymPy and commit the worksheets to
`v87/work/geom/`. The brief says the math will not be easy to construct; attempt
it regardless, and if it defeats you, record precisely where and why — a sharp
statement of the obstruction is worth more than a vague construction.

---

## 6. Relationship to the other two seats

FABLE and GROK are working the same problem from the doors-and-costs angle
(`BELL_BRIEF.md` §3), without this geometric hypothesis. You are the **GEOM**
seat. Read their files; they may already have settled which assumption must be
relaxed, which directly constrains your §3 Constraint 2. Where you disagree with
them, quote and argue in your own file. Where their reachability answer plugs
into your bound, say so — the ideal outcome of running three seats is that the
two halves join.
