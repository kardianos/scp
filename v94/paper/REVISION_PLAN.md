# Revision plan for `scp_theory.md` (response to grok's pedagogical review)

**Session ID:** `aebbc7eb-4716-4e30-8c22-484e115428f2` (resume 4, task B).
**Input:** `REVIEW_grok.md` (full pedagogical review, HS-reader lens).
**This file:** (1) author responses to the four "muddy" law-mapping grades,
(2) section-by-section revision plan, (3) status of the top-3 fixes —
items marked **DONE (this pass)** were executed in this revision; items
marked **PENDING** are planned but not yet applied.

---

## 1. Responses to grok's findings

### Maxwell — grade "muddy": AGREE

The Maxwell box assumed the reader has seen curl-form Maxwell; the stated
audience has not. The structural claim ("two quantities, each sourcing
only the other's *change*") is the real content and survives without any
`∇×`. **Revision:** prose-first restatement of the familiar law ("the
electric field's change is caused by the magnetic field, and vice versa");
curl equations retained but explicitly marked as optional notation for
readers who know them; the `-i` in the exchange rotation grounded in the
new quarter-turn figure (a `-i` = "the neighbour's arrow arrives turned a
quarter-turn," the same 90° relationship the quadrature pair has). The
Schrödinger aside stays but is fenced as an explicitly skippable remark.
**Status: prose-first + skippable fencing DONE (this pass); the
"E drives B / B drives E ↔ fa1/fa2" side table PENDING.**

### E=mc² — grade "muddy (concept clearer than equation)": AGREE

"Ledger identity" lands; the *formula* correspondence (where did c² go?)
does not. **Revision (PENDING, next pass):** (a) a worked ledger example
with law-table numbers — one atom leaves matter, `E_m` drops by
`ε = A₀w₂ₑ/2π ≈ 0.38`, field gains exactly that, drift meter unchanged;
(b) the explicit sentence "in units where the substrate speed is 1,
`E = m`; the `c²` of the SI formula is a unit conversion, nothing more";
(c) the open inertia/momentum caveat moved under its own subhead, "What
E=mc² does *not* yet give us," so the honest gap reads as a labelled gap
rather than an undercut.

### T⁴ — grade "muddy (adopted shape, not derived)": AGREE

The reader should never discover the exponent is adopted *after* seeing
the formula presented in derivation position. **Revision:** the "adopted
form, battery-gated; whether the substrate *must* use 4 is open" sentence
moved to *before* the equation; an explicit caution that the load `x` is
temperature-*like*, not a thermometer reading. The mode-counting sketch
stays, clearly labelled a sketch. **Status: adopted-first reorder DONE
(this pass); the `x` vs `T` caution DONE (this pass, same edit).**

### Poisson/Newton — grade "muddy near field if jargon stays": AGREE

"Green's-function well" and "monopole" are grad-school shorthand doing no
work for this audience; the physics ("a permanent drain makes a permanent
dip; things get pushed down the dip") is fully expressible in plain
language. **Revision:** Green's-function phrasing replaced with the
drain/dip/pressure formulation; "monopole tail" glossed as "a steady
`1/r` tail" at first use; π-flatness renamed in-line to "the flat-medium
result" with the internal label given once in parentheses. The far-field
honesty stays verbatim in spirit — grok correctly identifies it as a
strength. **Status: Green's replacement DONE (this pass); the explicit
map table (deficit ↔ potential dip, conversion ↔ source density,
push ↔ falling) PENDING.**

### Double slit — grade "clear": no change beyond micro-glosses

"Unitary" glossed in six words; Born explained before named. **Status:
DONE (this pass).**

---

## 2. Section-by-section plan

| § | Keep | Rewrite / add | Demote |
|---|---|---|---|
| Abstract | thesis, honesty framing | add "all symbols: Appendix A" pointer early in §2 instead | "three-node mode fiber" → "a three-station conversion path" |
| §1 | one law, four constraints, claim tags (grok: "teaching gift") | add reading-path note + correspondence-table-as-roadmap pointer (PENDING) | — |
| §2 | cells/roster, three sectors, register, ledger | **DONE:** new §2.3 "The arrow and the quarter-turn" + Fig. arrow; **DONE:** flight-vs-register cartoon (Fig. transit); **DONE:** `2π` motivation for ε; Appendix-A pointer; PENDING: everyday anchors for the three modes (light / stuff / room), action toy number, HS-math half-page preamble | — |
| §3 | precession+exchange words; component pair + shape commentary (grok: "best equation pedagogy") | **DONE:** pulse-hop time-step cartoon (Fig. pulsehop); **DONE:** prose-first Maxwell; quadrature defined in body; `-i` tied to quarter-turn | **DONE:** Schrödinger aside fenced skippable; PENDING: 5-step 1-D `c_sub` mini-derivation or demotion to appendix |
| §4 | register-necessity §4.1 (grok: "conceptual heart — protect") ; beat/click | PENDING: door pair introduced as "the §3 primitive aimed at (ψ_e, χ)" *before* the complex ODE form; detuning-vs-load toy table | PENDING: E=mc² worked example + "not yet" subhead |
| §5 | radiance fixed-point story; fold narrative | **DONE:** T⁴ adopted-first + `x`≠thermometer caution; **DONE:** "nonperturbative" glossed; "fold" introduced as "a point of no small reversible dressing" | **DONE:** ω± formula out of main text → Appendix C; main text keeps the avoided-crossing *picture* in plain words ("mixed wave", the word polariton given once) |
| §6 | pressure-push physics; far-field honesty | **DONE:** Green's-function → drain/dip language; monopole glossed | PENDING: explicit map table |
| §7 | everything (model section) | **DONE:** unitary glossed; Born explain-then-name | — |
| §8 | the claim + three irreducibles | **DONE:** plain-English theorem paragraph *first*; pass-F/U/P labels replaced by "field hop / matter hop / door"; "stroboscopic" replaced | PENDING: forbidden-alternatives cartoon; formal proof sketch → appendix (kept in place this pass, shortened) |
| §9 | table | PENDING: also cite it in §1 as roadmap | — |
| §10 | list | PENDING: light gloss pass (synchronous detection, accretion) | — |
| App. | notation, method | PENDING: 1-page glossary; Appendix C added **(DONE)** for the mixed-branch formula | — |

## 3. Top-3 fixes — status

1. **Arrow / i = 90° / energy = radius² figure + prose; complex ODEs
   optional** — **DONE (this pass):** new Fig. (fig_arrow.pdf) + new
   §2.3; `-i` referenced back to it in §3.1; Schrödinger fenced. The
   full "every `i d/dt` equation optional" restructure of §4.2 is
   PENDING (needs the geometry-first door-pair presentation).
2. **Demote polariton algebra + specialist nouns; process cartoons** —
   **DONE (this pass):** ω± → Appendix C; noun pass (Green's function,
   pass-F/U/P, stroboscopic, nonperturbative, unitary, Born, monopole,
   quadrature, polariton); two process cartoons added (pulse hop; flight
   vs register).
3. **§8 plain-English-first** — **DONE (this pass)** in reduced form
   (plain statement precedes the formal theorem; labels cleaned).
   The forbidden-alternatives cartoon and full proof-to-appendix move
   are PENDING.

Remaining PENDING items are ordered for a next pass: (a) §4 geometry-first
door pair + E=mc² worked example; (b) §3.3 `c_sub` mini-derivation or
demotion; (c) HS-math preamble + roadmap + reading path; (d) glossary;
(e) forbidden-alternatives figure; (f) §6 map table; (g) birth/death
cartoon (grok's optional item 7).

---

## 4. Revision log — resume 5

**Format fix (user request): DONE.** Front matter (title, abstract, TOC)
is now one-column; the body switches to two columns via a raw
`\twocolumn` at the start of the markdown body. The pandoc invocation
drops `-V classoption=twocolumn` (README updated). Verified: p.1 is a
full-width title/abstract/TOC page; the body opens two-column at §1.
15 pp.

**Grok items executed this resume:**

- **E=mc² worked ledger example + units sentence + split subheads**
  (from §1 "muddy" response) — DONE: one atom fired outward moves
  exactly `ε ≈ 0.38` between ledger columns; "`c²` is a unit
  conversion, not a physical ingredient"; the open inertia item now
  sits under its own "What it does *not* yet give us" subhead. This
  also covers the "action atom toy number" item from grok §7.B.
- **Reading path + table-as-roadmap** (grok §7.I) — DONE: closing
  paragraph of §1 ("Two ways to read"), pointing at Table 1 as a map.

**Still PENDING after resume 5** (unchanged priority): (a) §4
geometry-first door-pair presentation (door pair as the §3 primitive
aimed at `(ψ_e, χ)` before the complex ODEs) + detuning-vs-load toy
table; (b) §3.3 `c_sub` 5-step mini-derivation or demotion to appendix;
(c) HS-math half-page preamble; (d) one-page glossary; (e)
forbidden-alternatives cartoon for §8; (f) §6 explicit map table (the
prose map "dip ↔ Φ, conversion rate ↔ ρ, push ↔ attraction" landed in
resume 4 — the tabular form is optional polish); (g) birth/death-of-a-
cell cartoon; (h) Maxwell "E drives B / B drives E ↔ f_a1/f_a2" side
table; (i) everyday anchors for the three modes in §2.2; (j) §10 gloss
pass (synchronous detection, accretion).
