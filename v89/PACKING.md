# PACKING — mass as a compaction of field (user-directed line, 2026-08-02)

**User direction (2026-08-02), restated as five theses:**

1. Mass cannot survive without also taking into account charge and EMF.
2. The question to ask is: how much field can we densely pack into a single
   harmonic structure.
3. Mass is a compaction of field, built up by compatible harmonics.
4. Quantum teleportation is a momentary alignment through this dense field
   of harmonics — normally impossible.
5. A positive quark in a proton signals that the packing uses not just
   stable (tail-aligned) harmonics, but charged harmonics as well.

Subordinate to `PRINCIPLE.md`, `CHARGE.md`, `MASS.md`. This document walks
each thesis, maps it to standing machinery, and advances it where an
advance is provable. Measurement: `charge_pack.c` (log
`charge/runs/charge_pack.log`).

> **SCOPE NOTICE — 2026-08-02.** A second wave of this document (dynamics,
> 2D/3D long-range force, annihilation, a Gauss meter, a hedgehog) was
> written and has been **withdrawn in full**, along with its artifacts
> (`charge_dyn*.c`, `charge_xy.c`, `charge_gauss.c`, `charge_hedge.c`,
> `charge_breather.c`, `charge_sym3.mac`, `charge_sym4.mac`) and its
> proof table (P7–P14). The reasons are recorded in
> `KIMI_REVIEW_2026-08-02.md` §3; in short, several of those runs
> falsified their own headline in their own logs (a Gauss meter validated
> on a lattice containing no charge; a hedgehog whose index read 0 on
> every shell after relaxation; a "collision" whose log says no collision
> occurred; an evanescent-tail measurement whose amplitude grew with
> distance). They were also stage-built — fixed 2D and 3D lattices, i.e.
> permanent index sets — which `PRINCIPLE.md` constraint 2 forbids, so
> even the parts that were arithmetically fine could not have entered the
> record as fabric results. What remains below is the statics, which is
> a reduction of one closed cycle and is where this document's own
> machinery lives.

---

## 1. Thesis 1 — "mass cannot survive without charge and EMF" is already
## measured; it is the missing intake law

The MASS program's microscopic result (M-B1/M-B2) is that the dense-only
balance I(x)−O(x) has **no zero crossing on any axis** — no object closes
its books at constant mass, everything leaks, best case ~45k t.u. on V2s.
Read the user's sentence against this: the dense-mode books *cannot* balance
alone, because the only channels in those books are dense-mode channels. The
missing intake side is the **field channel** — conversion door F→D firing in
atoms (the standing law 1 machinery), and the charge texture's dressing
books (E10: charge self-energy is a load-bearing part of the nucleon mass
budget in reality — BMW decomposition: +2.5 MeV flavor − 1.0 MeV EM).

So thesis 1 is not a new postulate; it is the correct diagnosis of the
no-particle theorem: **the zero crossing will be found when the balance
books include the field intake and the charge dressing.** Concrete form
(queued as N8): re-measure the balance curve I(x)−O(x) for the M-B5/M-B7
objects with the field-conversion channel's books opened and the winding
content counted — the prediction of this document is that the crossing
appears on an axis that mixes dense load with field throughput, not on any
dense-only axis.

**Status: a diagnosis and a pre-registered prediction, not a result.**

## 2. Thesis 2 — the packing question has two answers, one measured today

"How much field can one harmonic structure hold" splits into two capacities,
because field packs in two distinguishable ways:

**Charge capacity (measured, `charge_pack.c`).** Pack k units of integer
winding into one cycle of N=100 voices (standing gate law, p=8, A=1):

| k | N/k | E/(k·E₁) | max step | resolvable cores | regime |
|---|---|---|---|---|---|
| 1–18 | ≥5.6 | 1.000–1.008 | 1.40–1.43 | k | charged: k distinct kinks |
| 20 | 5.0 | 1.035 | 1.75 | 20 | cores begin to touch |
| 22 | 4.55 | 1.035 | 1.82 | 16 | cores merging |
| **25** | **4.0** | **1.066** | 1.59 | **0** | **the relaxed lattice IS the uniform state** |
| 28–45 | 3.6–2.2 | 1.10–1.39 | 2.0–3.1 | 0–8 | past capacity |
| 50 | 2.0 | 1.462 | π | 0 | maximal strain (every site at the potential top) |

The charge-packing capacity of one cycle is **k\* = 25 at N=100, i.e. ~4
sites per charge**: below it the winding packs as an even kink lattice
(self-organized — see the merged-jump control below); at k\* the relaxed
lattice and the uniform-strain state coincide to 0.1% and the identity of
the individual charges is gone; beyond it there are no charges, only strain.

**Do not compress this to a single "ξ".** `charge_fk.c`'s independent
floor, from the energy excess on small cycles, is N/k ≥ 6, not 4. Both are
real measurements of different questions (when does the identity vanish
entirely vs when does the energy start to depart from additive), and the
capacity is a crossover rather than a sharp edge. CHARGE §8.3 lists all
four numbers that have been called the core width and why they differ.

**Field capacity (cellfab-side, stated).** The FK chain has no saturation —
uniform strain costs (A/2)δ² per site without bound. The real field-packing
limit lives in the standing law table: `cap` (a voice saturates at x→1) and
the vacuum skirt (load within 2Γ of the room pitch dissolves). How much
*field* one structure holds is therefore N·cap·x\* with x\* set by the
skirt — measurable in cellfab, not in FK. Noted honestly as the half of
thesis 2 that today is a law-table statement, not a number.

**The compatibility constraint on both:** not every packing is admitted —
closure (frequency closure Πp/q = 1 over the comb, Tenney ≤ 6) is the
admission ticket. The closure-allowed set is enumerated
(`charge/species_enum.py`, 17922 dihedral classes for N ≤ 8): thesis 2's
real answer is a *discrete set* of packings, not a continuum.

## 3. Thesis 3 — "mass is a compaction of field, built up by compatible
## harmonics": the cost accounting that makes it precise

The packing scan gives the economic content of "compatible". Compare two
**stationary states of the same ring** at fixed total winding W=8:

* the relaxed **kink lattice** (8 separated cores): E = 59.4621;
* the **merged jump** (all the winding on one bond): E = 67.7205.

The merged defect costs **1.14×** more. Kink repulsion (`charge_fk2.c` N1)
is what spreads the winding into the even lattice: a closed lock
concentrates the packed field into a few voices and leaves the rest of the
cycle at the vacuum, and eight such locks prefer to share the ring evenly.
As density rises toward k\* the distinction dissolves, which is exactly
where the charged/strain boundary sits.

> **Retraction.** The first version of this section quoted a **10.8×**
> ratio for "incompatible packing", from E_strain(k=1) = 80.6 against
> E₁ = 7.43. That comparison was against the *unrelaxed uniform ramp*,
> which at k=1 is not a competing packing at all — it relaxes to the kink.
> The ratio has no physical content and is withdrawn. The uniform state is
> a legitimate reference only at large k, where it *is* the relaxed state,
> and that is the only role it now plays (§2's crossover).

E=Mc² in this picture: mass = the packed field total, priced at the
conversion rate C. The packing scan is why the *amount* is what it is —
the structure's energy is its cheapest compatible packing. Whether that
equals the inertial mass is **not** settled here; the inertial-mass
theorem (m_tilt·C² = E_book on the M-B5 object ladder) remains unmeasured
and belongs to the MASS program on the annealed substrate.

## 4. Thesis 4 — teleportation as momentary alignment: framed, and made
## testable

Translate: through a densely packed medium, transport is *normally
impossible* — the gate forbids it (partial cycles convert nothing; the
loaded core is a mirror — g3 opacity is measured). The only channel through
is a **momentary full-cycle alignment**: a conversion process that closes
across the packed structure as one event. This is precisely the ontology e5
already runs on ("a pair is one conversion process" — the CHSH responder);
thesis 4 says that is not an approximation but the general mechanism:
correlation without transport, because the packed field's harmonics align
for one cycle and then unalign.

What this buys concretely — a pre-registered experiment (**N10, winding
transfer**): on a packed chain carrying one kink (charge) at site A,
inject a kink–antikink pair at distance R; the injected antikink
annihilates the resident kink; the injected kink remains — **the charge
state has moved A→B without any texture traversing R**. Bars: (i) winding
books exact at every step (charge conserved, ledger floor); (ii) no kink
detected in the intervening region at any time (the alignment, not a
traveler); (iii) transfer time ≈ the pair's local cycle time, independent
of R below the field crossing time.

**Status: pre-registered, not built.** It must be built in cellfab on a
V2s chain, not as a standalone FK dynamics — a fixed 1-D index set with
its own time integrator is a stage, and the withdrawn second wave is what
that costs.

## 5. Thesis 5 — the positive quark: charge is packing overflow

The advance, stated as a theorem candidate on the corrected monodromy
(`charge_sym.mac` PART 2):

**Closure of a packed structure splits into a payable part and an
unpayable part.** Around any cycle: M(θ)−θ = comma + 2πk/Q. The comma
(the α-sum) is *payable*: annealed by plasticity, dumped on seams, paid in
atoms (ε = A₀ω/2π, credit two — all standing machinery). The branch 2πk/Q
is *unpayable*: no continuous deformation and no atom payment changes k.
A packing whose interval content closes with Q=1 can be neutral; a packing
whose compatible harmonic content forces Q>1 **cannot close its books
without a residue — and the residue IS charge.**

So: charge is not an ingredient added to the packing; it is the packing's
un-closable remainder. "A positive quark signals the packing uses charged
harmonics" becomes: **the proton's harmonic content has no neutral
closure.** Its charge is forced by what it is made of, not attached
afterward.

And the two charge mechanisms separate cleanly:

* **Structural charge** (interval holonomy): carried by the closure
  arithmetic itself — no size floor, no soliton needed; confined to
  interval-rich structures (fractions cannot traverse unison — §7.1).
  *The quark mechanism.*
* **Textural charge** (kinks): needs the core to fit — measured floor
  N/k ≳ 6 at the standing parameters (CHARGE §8.3), integer only, mobile
  (PN barrier 2.2% of E_core), dressed by the far field.
  *The lepton mechanism.*

The electron is the textural integer (delocalized kink + dressing — §7.3);
the quarks are structural thirds. The proton is a packing that must use
both: closed locks for the bulk of its mass and un-closable interval
content for its charge.

**One honest gap.** The uud pattern (+2/3, +2/3, −1/3 = +1) is *stated* as
the minimal Q=3 packing with net k/Q = 1; nothing here derives it, and the
species table deliberately no longer emits slot assignments (CHARGE §8.6).
The energetic half of the structural/textural split also awaits the
fractional core cost on a ℤ_Q interval chain — the number withdrawn in
CHARGE §8.4.

## 6. What this adds to the program (the provable increments)

| # | item | status |
|---|---|---|
| P1 | Charge-packing capacity: k\* ≈ N/4 by the identity-loss criterion, N/k ≳ 6 by the energy criterion; the two differ and both are reported | **measured** (`charge_pack.c`, `charge_fk.c`) |
| P2 | Merged jump costs 1.14× the relaxed kink lattice at equal winding (both stationary) | **measured** (`charge_pack.c`) |
| P3 | Two charge mechanisms (structural holonomy vs textural kink) with different size laws | derived from the corrected monodromy; the fractional energetics are owed |
| P4 | Charge = unpayable closure residue (packing overflow) | theorem candidate; symbolic machinery verified |
| P5 | Zero crossing lives on the field-mixed balance axis | pre-registered prediction (N8) |
| P6 | Winding transfer = fabric-level teleportation test | pre-registered (N10), not built |
| P7 | Field-packing cap = N·cap·x\* (skirt) | law-table statement; cellfab measurement queued |

**Next experiments (queued), in program-value order:**

* **N8** — the balance curve with field-channel books (thesis 1's test;
  highest value; pure analysis on existing M-B5 machinery plus one new
  book).
* **N9** — the fractional core cost on a ℤ_Q interval chain. This is the
  number CHARGE §8.4 withdrew and the one that would let the
  structural/textural split be more than a classification.
* **N10** — winding transfer (§4), in cellfab.
* **N11** — the packing scan reproduced in cellfab on V2s: does k\* and
  the 1.14× survive the full gate physics? They should — the gate law IS
  the FK potential — and if they do not, the FK reduction is not a
  reduction and everything in this document is stage-built too.

## 7. Artifacts

| file | role |
|---|---|
| `PACKING.md` | this document |
| `charge_pack.c` | the packing-capacity measurement (log `charge/runs/charge_pack.log`) |
| `charge_sym.mac` | the corrected monodromy (the §5 theorem's machinery) |
| `charge_fk.c` / `charge_fk2.c` | core energy, size floor, kink repulsion (the §3 lattice self-organization rests on N1) |
| `charge/species_enum.py` | the closure-allowed packing set |
