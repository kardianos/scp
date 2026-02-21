1. The B=3 proton claim contradicts the established Skyrme identification (Section 13)                                
																													
The document says protons/neutrons are "three linked or intertwined vibrational knots forming a composite with baryon 
number B=3." But in the Skyrme model — which the document elsewhere claims to reproduce exactly — the proton IS the 
B=1 soliton. Baryon number = topological charge, so proton = B=1, neutron = B=1 (different isospin orientation), and
the B=3 soliton is a ³He/³H nucleus. The "three quarks" inside a proton are not three separate solitons in the Skyrme
model; they are internal degrees of freedom of a single B=1 configuration.

This is a fork-in-the-road decision: either HFKT follows the standard Skyrme identification (proton = B=1, which is
well-tested and matches the numerical results), or it proposes something genuinely different (proton = linked triple,
which has no numerical support yet and would change the meaning of all existing B=1-4 results). The document currently
conflates both without acknowledging the tension.

2. The Q identification is inconsistent throughout

Q is called "electric charge" in Sections 9, 14, and the topology discussion, but "baryon number" in Section 13 and
Appendix B. In the actual Skyrme model, Q = B (baryon number). Electric charge is derived from isospin: Q_em = I₃ +
B/2. These are not the same thing. The document needs to commit to one identification and derive the other, or
explicitly state the ambiguity as open.

3. Gravity is presented more confidently than the v2 results support (Section 18)

The section says five properties "emerge automatically" from depletion zones. But v2's extensive numerical
investigations found:

- In the sigma model (λ → ∞), ρ ≡ ρ₀ everywhere — there is no density depletion at all. The soliton is a pure twist
with no radial density variation.
- The BLV effective metric gives P/m = 2 algebraically, meaning no gravitational lensing from L₂ + L₄ alone.
- Only the L₆ sextic term breaks P/m = 2, producing nuclear-scale effects (10³⁷× stronger than Newton), not long-range
1/r.
- Paths 3-6 for gravity generation all returned null results. g_top is a free parameter at all levels tested.

The gravity mechanism described in the document is conceptually appealing but is actually the least numerically
validated part of the theory. The section should honestly flag that the depletion picture has not been confirmed and
that the effective metric investigation produced null results for long-range 1/r behavior.

4. Hopfion enthusiasm should be tempered by the v2 instability result (Section 13)

The document advocates upgrading toward "explicit Faddeev-Skyrme Hopfion linking" as the path for composite particles.
But v2's hopfion Langevin simulations showed that Faddeev-Skyrme hopfions dissolve at ALL temperatures (T=0 to 10⁻³).
The torus shrinks through the lattice via the same saddle-point mechanism that destabilizes Skyrmions. This is a
definitive negative result that should appear alongside the proposal, or the section should explain why the
instability is a lattice artifact rather than a continuum feature.

5. The L₂ blindness to the weight sector should be flagged at point of introduction (Section 6)

The Lagrangian is presented in Section 6, but the critical structural fact — that L₂ evaluates to ½|∂_μq|² and is
completely blind to the weight sector — doesn't appear until Section 17 and Appendix A. A reader encountering the
Lagrangian will assume all 8 components participate. A one-line note at Section 6 would prevent this misunderstanding,
since it is the single most important structural property of the dynamics.

6. Section 23 conflates spatial and degenerate bivectors

The null-rotor section says "the 6 bivector components (3 degenerate e₀ᵢ + 3 spatial eᵢⱼ) map exactly to the 6
components of the electromagnetic field tensor." This is wrong. The EM field is identified specifically with the
spatial bivector F (3 dof, massless). The degenerate bivectors J (3 dof, massive at μ) are the flux/weak sector. Only
the 3 spatial components map to EM. The 6-component statement works for the full spacetime Faraday tensor F_μν in
Cl(1,3), but not in the PGA decomposition where the two bivector types have different physics.

7. "Open waves have no persistence" overstates (Section 2)

"They carry energy but have no rest mass and no persistence — they disperse or are absorbed." Axiom 4 says no
dissipation. In a non-dispersive vacuum, a photon propagates indefinitely without losing coherence. It can be absorbed
or scattered, but it doesn't spontaneously disperse. The distinction from knots is topological protection (knots
can't be continuously unwound; waves can be absorbed), not persistence per se.

8. Section 10 ("What a Knot Looks Like") is poetically vivid but physically imprecise

"Field lines twist, compress, and form closed resonant loops. The internal wave energy recirculates coherently." The
actual B=1 hedgehog is a static twist: q = ρ₀(cos f(r) + sin f(r) r̂·σ). In its rest frame, nothing is oscillating —
it's a time-independent solution. The "vibration" language comes from the process ontology (the knot is a persistent
pattern that must be continuously re-realized as it moves), but in the rest frame it is a purely static twist. The
section should distinguish "static internal twist" from "oscillating/vibrating," or clarify that the vibration is the
process of maintaining the twist rather than literal time-dependent oscillation.

9. Parameter count (Section 7)

Says five parameters. But v2 Phase 8-9 added a sixth parameter g (bulk-degenerate coupling strength), bringing the
actual theory to 6 parameters (ρ₀, λ, e, μ, g, c). If v3 intends to drop or subsume g, that should be stated. If not,
the count should be updated.

10. Missing: lattice topology loss caveat

The "Resolved" section lists "Soliton scattering: topology preserved through deep interpenetration" but doesn't
mention that gradient flow loses topology at ALL tested resolutions (N=128-192), and even the time-evolution code
loses topology after ~2.4 time units. The sigma-model Skyrmion is a lattice saddle point. This is relevant to the
foundational claim that topological protection is absolute — in practice it depends on resolution, and v2 found no
resolution sufficient for indefinite stability.

11. Minor: the EOM in Section 8 omits V_D

The equation of motion shown doesn't include the degenerate mass contribution. The full EOM should have +μ²p terms for
the weight sector (or note that they reduce to μ²p = 0, forcing p → 0).

---
Structural suggestion: Consider adding a short "Caveats and Known Tensions" subsection (perhaps in Part VIII after the
Open list) that collects the internal tensions — gravity null results, lattice instability, Q identification
ambiguity, hopfion dissolution — in one honest place. v2's spec had "Warning: Unresolved Contradictions" at the top
level, and that intellectual honesty was one of its strengths. v3 currently distributes the optimism broadly but
concentrates the caveats narrowly.
