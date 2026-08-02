# Review of the GLM charge commit, and audit of the KIMI charge wave
# — 2026-08-02

Two passes are recorded here, in order:

* **§1–§3** — the KIMI review of the GLM charge commit (f3678a1). Its
  three defect findings are independently confirmed and are repaired in
  `charge_sym.mac`.
* **§4–§5** — an audit of the KIMI wave itself, which found a further set
  of defects, several of them fatal. That audit is the reason most of the
  KIMI artifacts no longer exist. **Read §4 before trusting anything
  attributed to that wave elsewhere in the tree.**

Subordinate to `PRINCIPLE.md`. Standing table `laws_V2g.cfg`; ratchet rule
governs. Kernel untouched throughout; all work in standalone files.

---

## 1. What reproduces from the GLM commit

* `charge_fk.c` reproduced its own printed tables exactly. (Whether those
  tables measured what they claimed is a separate question — see §2.4.)
* `charge_sym.mac` ran clean; I(1) and I(8) via `quad_qags` are right, and
  the max-torque locus 1.232578 at ψ=0.505360 is right.

The in-kernel Q11–Q13 / Q12 arc (circulator bistability, pump → pinned
strain, the dynamical basin separation, the mobile ℤ₃ wall exiting at
exactly −2π/3) is well-documented and internally consistent, and the
verification pass re-deriving §7.10's numbers from the run logs is exactly
the right discipline. **That arc is not touched by anything below** — it
ran in-kernel on the gated battery, not on a reduced chain.

## 2. Defects found in the GLM artifacts — all confirmed

### 2.1 `charge_sym.mac` PART 1 presented a failed symbolic check as a success

The file claimed I(1)=4 was verified symbolically. The `integrate()`
antiderivative is atan2-branch-broken and **printed 0**, not 4. The
verification that actually held was the numerical quadrature.

Repair, and it is exact: √V(ψ,1) = √((1−cosψ)/2) = |sin(ψ/2)|, so
I(1) = ∫₀^{2π}|sin(ψ/2)|dψ = 4 identically. Also added: I(2) exactly
= 2√2 + 2·asinh(1) = 4.591174298785276.

### 2.2 The monodromy composition was wrong, and it bit the physical cycle

`monodromy()` composed edge maps θ ↦ (p·θ + α + 2πn)/q as

    c := (c + α_l + 2πn_l)/q_l          — missing the factor p_l on c.

The correct composition, from CHARGE §7.1's own definition, multiplies the
accumulated offset by p_l at each subsequent edge:

    c := (p_l·c + α_l + 2πn_l)/q_l.

The printed triangle {3:1,1:3,1:1} is insensitive (p₂=p₃=1), so the bug was
invisible in every case the artifact showed. For the **fifth-triangle
{3:2, 2:3, 1:1}** — the actual Q11/Q12/Q13 cycle — the two differ:

    old code : c = α₁/6 + α₂/3 + α₃ + πn₁/3 + 2πn₂/3 + 2πn₃     (wrong)
    correct  : c = (α₁+α₂)/3 + α₃ + 2π(n₁+n₂)/3 + 2πn₃

The continuous comma was wrong (α₁/6 vs α₁/3) **and the n₁ branch
coefficient was wrong** (πn₁/3 vs 2πn₁/3 — a spurious ℤ₆-looking term).
The headline conclusion (branch = 2πk/Q, charge = k/Q) survives — the
correct branch is 2π(n₁+n₂+3n₃)/3 = 2πk/3, k ∈ ℤ₃ — but only after the
fix. A "general N-cycle monodromy" exercised only on a degenerate case is
a hazard. Repaired and re-verified on {3:1,1:3,1:1}, {3:2,2:3,1:1} (ℤ₃)
and a 4-edge mixed cycle {2:1,1:2,3:1,1:3} (ℤ₆).

### 2.3 The printed branch formula was mis-stated

The old file printed "discrete branch 2*pi*n1*n2*n3 / 3" — the branch is a
**sum** over edges, not the product n₁·n₂·n₃. Print-string error only (the
computed `c` showed the sum), but a reader lifting the formula got
nonsense. Fixed.

### 2.4 Two further GLM defects, found in the audit and not in the first pass

These were missed by the first review and are the more consequential ones:

* **`charge_fk.c` relaxed to the PN saddle, not the minimum.** A single
  uniform-ramp seed on an even-length chain is symmetric about a bond, and
  gradient descent preserves that symmetry, so it converged to the
  bond-centred stationary state: E = 7.598210 at p=8, A=1 (max|grad| ~
  1e-14 — a genuine stationary point, just not the ground state). The
  minimum is the site-centred kink at **7.432761**, 2.2% lower. Every
  ratio in the §8.1 tables was high for this reason; p=8 read 0.994 of the
  Bogomolny bound where the true figure is 0.972. Fixed (multi-seed,
  lowest kept); CHARGE §8.1 corrected.
* **CHARGE §8.4's fractional-charge scan measured an artifact.** See §4.1.

## 3. Assessment against the goal (stable proton & neutron)

**MASS.** Best object: the M-B7 wound composite on V2s. The blocker is
microscopic and named: the balance curve I(x)−O(x) has **no zero crossing
on any axis** (M-B1/M-B2) — the intake-side law is missing, so no object
can close its books at constant mass. Nucleon mass is additionally blocked
on charge (E10).

**EMF.** EM5 is built and gated (leapfrog E-links/B-faces on the link
graph; cone ω = 0.672k + 0.033). The Maxwell operator exists; what is
missing is the **source term** ∇·E = ρ, and there is no ρ because there is
no seeded static charge texture on the link graph.

**CHARGE.** The O3′ ontology has real measured content: quantization (the
branch arithmetic, now with the right composition law), confinement (the
unison-transport argument plus the in-kernel pump experiments: fractional
deposits pin, integer textures move), C-symmetry (exact, 1:10⁶), the
Josephson I–V (all three bars), the integer core energy, the kink force
law, and the PN barrier. What it lacks for the nucleon: a Gauss meter
(E3), a two-charge force at range (E4), the fractional core cost, and a
*balanced* host cycle to put the charge on.

**The conjunction that matters:** the three blockers are one blocker seen
from three sides — *no object yet maintains a steady conversion cycle at
constant books*. GLM's conclusion (livefab is the keystone) is correct and
I endorse it.

---

## 4. Audit of the KIMI wave — what was withdrawn and why

The KIMI wave produced, over three sub-waves, a set of standalone models:
FK statics (kept), FK dynamics, a 2D phase lattice, a 3D O(3) lattice, a
plaquette index meter, and a breather study. **All but the statics are
withdrawn and deleted.** Two independent grounds, either sufficient.

### 4.1 Several runs falsified their own headline in their own logs

| claim as written | what the run log actually said |
|---|---|
| "the quark pattern, measured": ⅓ charge costs 0.60 of the linear share | On a unison ring 2π/3 is not a period of V, so there is no soliton. The relaxed state puts **1.3948 rad on the ring's artificial closing bond** and leaves 398 of 400 sites at the vacuum. The energies scale as w² (strain). |
| "Gauss/index meter PROVEN, exact to 7e−17, surface-independent" | Every ε reported `+1:0 -1:0`, total winding 0.0000 — **no vortex survived relaxation in any run.** The exactness statement was about 2304 zeros. |
| "π₂ index meter PROVEN; hedgehog E ∝ L, lattice factor 0.61" | `deg = 0.0000` on every shell after relaxation, for every R — the hedgehogs unwound. The E ∝ L number came from a log **no committed source could produce**, and its normalization (box side vs radius) was never stated. |
| "k+k̄ collision at v=0.2 → bound breather + ~30% radiation" | The cited artifact's log ends `# D1 note: no collision (Emax=14.8713)` with Ekin ≈ 0.005. The companion artifact ("dyn6") did not exist. |
| "the meson slot is occupied" (breather stable ≥2000 t.u.) | The spectrometer returned ω = 0.0628 = exactly 2π/(FFT window) — the k=1 bin of an unmeaned, still-decaying signal, i.e. the residual trend, not a mode. E_total was still falling at t=5000. |
| "below-gap drive suppressed ≥10³ (evanescent)" | Below-gap RMS **grew 7× with distance**, where κ=1.57 evanescence predicts decay by ~10⁻⁸². The number quoted was an amplitude ratio at fixed drive, not a decay. |
| "gapless ⇒ log-Coulomb, slope 4.5–6.3 vs 6.28" | That range is [one measured value, the prediction]. The (++) branch measured −0.93 against −6.28 — a 6.8× miss — and was not mentioned. |
| PN barrier 1.27% of E_core | The sweep pinned a site at π ± 0.5 rad, but π *is* the minimum and the saddle is 0.79 rad away: it never left the bottom of the well. True barrier **2.23%**, a factor 1.75 higher. |
| species table: 175 classes, "Q=3: 5 classes" | The enumerator used sorted edge multisets, so it saw one cyclic ordering of each — but Q depends on ordering, as its own docstring argued, and the example that docstring cited was among the cases it dropped. Ordered enumeration gives **17922** classes. |

Two items from that wave did survive audit and are folded into
`charge_sym.mac`: the phonon band ω ∈ [2, 2√2] with max group velocity
√2−1 (re-derived independently), and the static screening ratio
r = 3−2√2 (the root of r + 1/r − 2 = V″/A, and the cleanest single
measurement of the wave). Their host artifacts are gone, so they are kept
as the symbolic statements they always were.

### 4.2 Most of it was stage-built, which the principle forbids

`charge_xy.c` was a fixed 2D lattice of phases; `charge_gauss.c` the same;
`charge_hedge.c` a fixed 3D lattice with an S² field; the `charge_dyn*`
series a fixed 1-D index set with its own time integrator. These are
permanent index sets whose values evolve — the exact construction
`PRINCIPLE.md` constraint 2 forbids and that the README warns "comes back
through inherited habits of construction."

The wave's own ontology audit noticed this and proposed a mitigation
(mark stage-built results MECHANISM-ONLY). That mitigation was not
applied: the status tables still read PROVEN. And the mitigation would
not have been enough anyway, because a stage-built result cannot enter
the record as a fabric result at all — which is the whole point of the
constraint.

The 1-D FK reduction that remains is defensible on narrower grounds: it is
a reduction of **one closed cycle** of fixed voice count, with the standing
gate law as its potential, which is what CHARGE §7.2 says the slip sector
is. Even there the reduction leaks — the fractional-twist artifact of §4.1
is precisely the ring's artificial seam being mistaken for physics — so
every FK result now carries the boundary diagnostic that would expose it.

### 4.3 What the FK statics keeps

Repaired and re-run; logs in `charge/runs/`:

* **Integer core energy.** E_core = 7.432761 at p=8, A=1, additive in |k|
  to all printed digits, 0.972 of the continuum bound (the lattice
  correction, in the right direction and magnitude for a sub-lattice core).
* **Size condition.** The energy excess collapses onto N/k: ≥6 costs
  <0.3%, 4 costs 6.6%, 3 costs 15.4%. The k=1,2,3 rows at equal N/k agree
  to every digit.
* **Two-kink force law.** sign(E_int) = sign(q₁q₂), monotone in R,
  antisymmetric where both signs resolve, range ~1 site — consistent with
  the tail rate κ = √(p/2A) = 2/site. Integer windings only; the
  fractional version belongs on a ℤ_Q interval chain.
* **PN barrier.** 0.165449 = 2.23% of E_core at p=8, A=1, from the two
  symmetric stationary states of the translation path. Which of the two
  is the minimum flips with (p, A), and the winner cross-checks against
  `charge_fk.c`'s independently relaxed E_core at every shared (p, A) —
  a consistency test the two artifacts previously failed.
* **Mobility size floor.** The barrier holds 2.226e-2 to four digits from
  N=400 down to N=8, then collapses below the FP floor at N=6: the two
  symmetric states merge and there is no distinct kink position left to
  pin. Floor N ≥ 8 in mobility form, against N/k ≥ 6 from the energy
  excess — different questions, reported separately.
* **Packing capacity.** k\* ≈ N/4 by identity-loss; merged jump costs
  1.14× the relaxed lattice at equal winding.
* **Closure/holonomy census.** 17922 classes for N ≤ 8; every charge
  quantum is 1/(2^a·3^b).

## 5. Method note — the failure mode, stated once

Every defect in §2.4 and §4.1 was invisible in prose and visible in the
run log. A seed that converged to the wrong branch. A sweep that never
reached the saddle. An enumerator that skipped 99% of its domain. A meter
validated on a lattice containing no charge. A collision that did not
happen, in a file whose own last line said so.

The common cause is not carelessness in the physics; it is that results
were written up at PROVEN/MEASURED **before the log was read back against
the claim**. The standing correction:

1. An artifact prints its convergence witness (max|grad|, residual drift,
   iteration count) and its negative controls, and the doc quotes them.
2. A claim cites the log line that carries it. A number without one is not
   a result.
3. Before a status word is assigned, the log is read for the thing that
   would falsify the claim — not for the thing that would confirm it.
4. Where the state lives is declared before the run, not audited after. A
   permanent index set is not a fabric result at any status level.

## 6. Standing queue

1. **N8** — balance curve with field-channel books in cellfab (thesis 1's
   direct test; the zero-crossing hunt). Highest program value.
2. **N9** — the fractional core cost on a ℤ_Q interval chain: the number
   CHARGE §8.4 withdrew, and the only route to a fractional charge
   energetics that is not a unison-ring artifact.
3. **N10** — winding transfer (PACKING §4), in cellfab.
4. EM5 source hook + a Gauss meter on a pinned Q12 wall — in-kernel, on
   the link graph, not on a stand-in lattice. Needs user authorization.
5. livefab cellfab integration (the keystone). Needs user authorization.

## 7. Artifacts

| file | role | status |
|---|---|---|
| `KIMI_REVIEW_2026-08-02.md` | this review + the audit | current |
| `charge_sym.mac` | exact I(1), I(2); corrected monodromy; exact torque locus; V″(0), ξ, tail rate | verified |
| `charge_fk.c` | core energy, size floor, integer additivity, fractional negative control | verified |
| `charge_fk2.c` | two-kink force law; PN barrier | verified |
| `charge_pack.c` | packing capacity; merged vs lattice | verified |
| `charge/species_enum.py` | closure/holonomy census | verified |

Deleted by this audit: `charge_dyn.c`, `charge_dyn3.c`, `charge_dyn4.c`,
`charge_dyn5.c`, `charge_breather.c`, `charge_xy.c`, `charge_gauss.c`,
`charge_hedge.c`, `charge_sym2.mac` (folded into `charge_sym.mac`),
`charge_sym3.mac`, `charge_sym4.mac`, and their logs. They remain in git
history; they are not part of the record.
