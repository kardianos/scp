# S2 — the amplitude-completion frontier (queue #8): entry measurement + design memo

Executed 2026-08-04. Queue #8's deliverable in two parts: (1) the
**entry criterion measured on the free substrate** — kappa_reac=1 (the
raw derived choir term, `v89/s2/choir_pull.c`) across the battery
surface and every acceptance face; (2) the **design memo**: what the
completion must be, given how the rate-level candidate fails. The law
table is untouched — kr=1 is a declared per-run override, printed in
every header. Logs: `runs/s2/`.

## 1. The entry measurement: kappa_reac=1 on the free substrate

v89 measured kr=1 = 16/17, failing only e3b. The v90 apparatus
(composites, FCQ, XSEC, DS) sees more:

| face | kr=0 (gated) | kr=1 (measured) | verdict |
|---|---|---|---|
| bath (vacuum) | drift 0, births 4611 | **identical** (Sm=0 ⇒ term exactly zero) | untouched |
| pulse (field) | v/C 0.5740 | **0.5740** | untouched |
| ds tier-0 (optics) | exposure 71.963056 | **71.963056** | untouched |
| pair (unison bond) | d_final 1.538698, off +0.0329 | 1.538969, off +0.0332; up=dn to 6 decimals | holds (slightly tighter) |
| ring6 (even parity) | dev 0.0350, gg ≥ 0.97 | dev 0.0329, gg ≥ 0.980 | holds |
| **ring5 (odd parity)** | dev 0.1319, **≥3 of 5 gates dead** | dev **0.0551, only 1 gate dead** (others 0.95–0.98) | **signature weakened: the π defect LOCALIZES onto one edge instead of distributing** — 2 gated bars would go red |
| blob (retention) | 0.8052 | 0.786 (bar ≥ 0.78) | holds, degraded 2.4% |
| **e3b (driven translation)** | speed 0.002580, cos 0.9508 | speed **0.001617, cos 0.3010** (kr=0.6: 0.000532 / 0.653) | **DEGRADED — slower and incoherent; the v89 failure reproduced on free cells, with magnitudes** |
| **FCQ UUD (the fifth as binding)** | closed, gates alive | **ggm = 0.0000; D voice bloats to Em 2.06 while the U's drain to 0.50** | **the fifth dies; mass piles onto D** |
| **composite rings nv=6/12 (boundary fifth)** | bond locks gg 0.989 (CQ1) | fifth edge **gg = 0.0000**, ψ ≈ +1.1/+2.7; D-ring x̄ → **0.83**, Em 12.5/25.1 | **the boundary fifth dies; the D droplet becomes a mass sink** |
| xsec headroom absorber | net_tag 7.274 | 8.101 (+11%, still rough=evap=0) | qualitatively same |

**The mechanism, read off the table:** the reactive term is a
*unison instrument*. Its −κ·g·sin(ψ) pull drives phase errors toward
lock, which even ring5's frustrated cycle partially exploits (the π
defect concentrates instead of spreading). But on a FIFTH (3:2 chart)
the same pull suppresses the forward want whenever ψ oscillates — the
gate collapses, the fifth stops transporting, D cannot shed, load
detunes it further, and the channel dies with the mass piled on the
wrong side. The choir correction retires the fifth everywhere the
fifth is load-bearing: FCQ binding, composite boundaries, and the
tilt-driven transport window (e3b's coherence is fifth-adjacent: the
drive rides phase gradients the correction flattens).

**Entry verdict: kr=1 is NOT the completion and is now measured to be
further from it than v89 knew** — it fails e3b (as in v89) *and* kills
flavored structure (new, worse). The entry criterion stands as the
acceptance test precisely because it currently fails: any proposed
completion must pass the full battery with the derived term at its
raw strength, or derive its retirement.

## 2. The design memo: what the amplitude completion must be

Everything measured in v89–v90 converges on one statement: **the dense
sector transports magnitudes; reality transports amplitudes.** The
field sector (two-plane signed amplitude, unitary hops) produces
interference, Huygens wavelets, clicks that rebuild fringes, HOM
exchange signs — every wave phenomenon we have tested. The dense
sector (want-driven magnitude transfers, phase used only to GATE)
produces bonds, composites, absorption — and none of the transport
phenomena. The completion is to give the dense sector what the field
sector already has, without losing what the dense sector does:

1. **Coherent dense channel.** Dense transfers must carry a complex
   amplitude (magnitude AND phase) whose phase composes along paths —
   "translation IS the current" (v89). Rate-level corrections on wants
   (kappa_reac, kappa_freq) are the wrong-level object: both were
   derived, both measured non-load-bearing or harmful. The e3b test:
   a tilted blob's forward window must become seed-robust and
   coherent (cos → 1), not a ±1e-3 residual.
2. **The fifth must survive coherence.** The 3:2 chart must transport
   amplitude with a p:q phase map (θ → (q/p)θ), not be flattened
   toward unison. The kr=1 sweep is the negative template: any
   completion that kills `rint1 gg` on the nv=6 boundary or bloats
   the D-ring past x̄ 0.4 fails the composite face.
3. **Exchange registry.** Two-quantum amplitudes need a sign under
   exchange (the v89 HOM tier measured g_b 0.42 < 0.5 < g_f 0.58 at
   mode-match-limited depth). Exclusion (PAULI-0's successor) is
   testable only with this: identical excitations refusing one state
   in OCCUPANCY, distinct from cap saturation (which is measured
   pitch-blind).
4. **The flux-machine interior.** Five independent composite
   measurements (CQ7, CQ8 A/B, CQ8-long, H-stiffness refutation,
   bath-retune) force it: a stationary composite pitch requires
   steady internal circulation. The amplitude channel must support
   persistent internal currents at zero net transport — the discrete
   analog of a bound state's probability current.
5. **What must NOT change.** The atoms machinery (grains ε = A₀ω/2π,
   two-atom credit — the ħ-analog, linear to 1e-8), the conversion
   door (cap, headroom, evap — PAULI-0 and XSEC live there), the
   space sector (pressure pushes, footprint g1, skirt), and the
   contact rule. The battery's 93 bars are the invariant surface; the
   completion lands as a new mode of the existing channels, gated by
   the same ratchet.

**Acceptance surface (quantified, from this cycle's measurements):**
e3b speed ≥ 2.6e-3 with cos ≥ 0.95 seed-robust; nv=6 boundary fifth
gg ≥ 0.9 at T=100 with D-ring x̄ ≤ 0.35; UUD ggm ≥ 0.9; composite EMF:
a flavored composite must RADIATE (rough > 0 sourced at the boundary
fifth) instead of staying dark; exclusion precursor: identical vs
fifth-consonant occupancy statistics at a Y-junction diverge. Full
battery green throughout, kr at its derived value.

**Status: the completion itself is a law-level construction — v91-scale
work (or explicit sign-off to modify the v90 law).** This memo + the
entry table is the queue-#8 deliverable that makes it actionable.

## Addendum (2026-08-05, post-CANTUS — the §2 item 3 upgrade, measured)

The CANTUS campaign (`v91/CANTUS.md`) measured the churn bath to be a
**persistent-closure medium**: 66% of cells hold a time-averaged
two-sided gate > 0.5 at zero lock force (the null-meter arm).
Consequence, measured twice (cell-borne v1, link-borne v1.1): any
coherence law that self-grows on closure statistics ignites a
bath-wide Kuramoto transition and rewrites the thermodynamics (glow
suppressed 10–50×). §2 item 3 (the exchange registry) is therefore
UPGRADED: from "the door to exclusion" to **the gate on which
coherence self-assembly itself hangs** — the bond-vs-churn
discriminator must be identity, not phase alignment. Two further
measured constraints on any completion design: charts must be CHORD
charts (unison cannot connect at the balance — the derived geometry
wall, CANTUS §3.2), and ring-locking needs winding-sector bookkeeping
(the differential lock conserves loop winding). Active lane:
`v91/REGISTRY.md`.
