# LIT — physics/math similarity search: harmonic cycles, settling, annealment
(2026-07-31, research agent report, received verbatim; commissioned per user
directive after the Q11h scan. Mapped against §7.9–7.10 results.)

## Master correspondence table

| # | Measured result (this system) | Closest established system(s) | Canonical references | Match quality |
|---|---|---|---|---|
| 1 | Bistable triangle: quiet locked state vs persistent circulator; slips antisymmetric around loop (Kirchhoff); winding localized at one vertex | Fluxoid sectors in superconducting rings (metastable winding states); trapped vortex in discrete 3-junction ring; q-twisted states vs sync state on oscillator rings; discrete sine-Gordon kink pinned at one site | Langer & Ambegaokar 1967; McCumber & Halperin 1970; Wiley, Strogatz & Girvan 2006; Delabays, Coletta & Jacquod 2016; van der Zant et al. 1994; Braun & Kivshar 1998/2004 | Strong, near-quantitative |
| 2 | Circulation rate ~interval-independent; on 6-cycle per-edge rate 0.43× triangle, loop total conserved (0.044 vs 0.038) | Fluxon revolving in an annular/discrete Josephson ring: zero-field steps, per-junction slip rate = v/L ∝ 1/N at fixed circulation speed | Davidson et al. / Ustinov annular-JJ literature; Watanabe, Strogatz, van der Zant & Orlando 1995 | Strong, quantitative (1/N law: predicted 0.50, measured 0.43) |
| 3 | Josephson relation ν = Δω/2π verified; knee; partial-lock 0.44×Δ/2π | RSJ/Adler equation: V = R·sqrt(I²−Ic²); SNIC bifurcation; type-I intermittency | Adler 1946; Stewart 1968 / McCumber 1968; Tinkham; Kye & Kim 2000 | Exact (0.44 ⇒ Δc ≈ 0.90Δ) |
| 4 | 3:2 tongue chatters (867 flips/600 t.u.) under flight load; 2:1 never flips | Arnold-tongue widths shrinking with rational complexity; Kramers escape from narrow vs wide tongues; van der Pol & van der Mark's "irregular noise" between subharmonic steps | Arnold 1961; Jensen, Bak & Bohr 1983/84; van der Pol & van der Mark 1927; Kye & Kim 2000 | Strong, mechanism-level |
| 5 | Closure defect res = −ω_high·(perimeter) mod 2π; ℤ_p branch classes; circulation independent of static defect, controlled by seeded winding | Frustration f in SQUID loops / JJ arrays; fractional flux quanta Φ₀/2, Φ₀/p (tricrystal π-rings, 0–π junctions, charge-4e/6e kagome rings); fluxoid number n set by history, f only tilts sectors | Tsuei & Kirtley 1994/2000; Mints et al. 2002; Goldobin et al. 2004–05; PRX 14, 021025 (2024); LAMH theory | Strong conceptually; mixed-p:q loop arithmetic itself is novel |
| 6 | 40× leak spike at exact interval coincidence with ambient; tongue-free desert with 10× further suppression | Discrete-breather / oscillon decay when a harmonic enters the phonon band; exponentially small radiation off-resonance; Fiske resonances | Segur & Kruskal 1987; Flach & Willis 1998; Flach & Gorbach 2008; Ustinov et al. 1998 | Strong, mechanism-level |
| 7 | 9–18× lower bleed on annealed (σ_d 3%) vs disordered (σ_d 9%) substrate | Roughness-limited scattering loss ∝ σ²; phonon boundary specularity p = exp(−4σ²k²) | Payne & Lacey 1994; Ziman 1960; Maznev PRB 2015 | Quantitative: (9/3)² = 9 vs measured 9–18 |

## Quantitative transfers adopted as pre-registrations (v89-native tests)

1. **RSJ knee law (from #3):** ν/ν_free = sqrt(1 − (Δc/Δ)²); measured 0.44 ⇒
   Δc ≈ 0.90·Δ — the flight-load detune sits ~10% above the 3:2 tongue edge.
   TESTS: (a) reduce flight load 10–15% → slip rate collapses to 0;
   (b) ν² vs Δ² is a straight line with intercept −Δc²;
   (c) inter-slip intervals ∝ (Δ−Δc)^(−1/2) near the knee (type-I
   intermittency).
2. **Zero-field-step 1/N division (from #2):** ν_edge·N flat to ~15% across
   N = 3,4,5,6,8,12 cycles; the hexagon's 14% deficit (0.43 vs 0.50) is a
   fluxon speed measurement (flight-load drag on the longer path). The 0.43
   also favors a near-linear (sawtooth-ish) CPR over sin — consistent with
   the measured Q9 CPR.
3. **Noble-address quiet floor (Greene/Shenker, from #5/#6):** the
   locking-resistant extremum of the band is the noble ratio ω/Ω = 1/φ =
   0.6180 ⇒ x_noble = (φ·... ) = 0.515 on the V2s pair. PARAMETER-FREE:
   the desert minimum should sit AT x=0.515, not merely "between tongues."
   (M-B4 grid measured min at 0.536 with 0.486/0.636 neighbors — probe
   running at 0.515.)
4. **π-ring transplant (Tsuei–Kirtley, from #5):** a cycle built with
   closure defect res = π exactly (perimeter-tunable) should SPONTANEOUSLY
   bifurcate into degenerate ± circulators with no seeding.
5. **LAMH sector lifetime (from #1):** circulators are immortal below a
   sharp transport-noise threshold and die by SINGLE discrete unwinding
   events (never gradual decay); dwell ∝ exp(barrier/flight-variance).
6. **Kink-width law (FK, from #1):** winding texture width ∝
   sqrt(tongue width) — octave windings (γ/2) should spread over more
   edges than fifth windings (γ/6); measurable on larger cycles.
7. **Kibble–Zurek annealment (settling protocols):** ramp the pin release /
   coupling at rate 1/τ_Q → trapped-circulator probability P ∝ τ_Q^(−σ)
   (annular-JJ measured σ ≈ 0.5); slow anneal selects neutral, fast quench
   traps charge. Random-texture seeding: P(q) ∝ exp(−k·q²), σ_q ∝ √N,
   winding decided in t ∝ log N.
8. **Roughness σ² law (from #7):** leak vs σ_d at intermediate anneals →
   log–log slope 2.0; deviations fix the disorder correlation length
   (Payne–Lacey form factor); frequency dependence p = exp(−4σ²k²) gives a
   two-knob (σ, ω) surface with zero free parameters once calibrated.
9. **Golden-rule tongue drain (from #6):** the 40× vacuum-fifth spike has a
   Lorentzian profile of half-width γ/6 in pitch (Q ≈ 116); leak at
   coincidence scales linearly with contact number (DOS factor).

## Novelty verdicts (no clean literature analog found)
1. Per-link LIVE rational arbitration (Tenney-weighted Lorentzian re-vote
   each step) with chatter as an observable.
2. Mixed-interval cycles: heterogeneous p:q locks, per-edge ℤ_p branch
   classes, retardation-frustration res = −ω_h·perimeter — candidate
   theorem territory (covering spaces over cycle graphs with edge-wise deck
   groups).
3. Interval-independence of circulation rate stated as an invariance.
4. Transport-load as the intrinsic noise source closing the loop through
   the droop law ω = Ω/(1+q·x) ("flight-load self-droop").
Established & transferable as-is: RSJ knee/partial lock; 1/N division;
seeded-winding (field-cooling) selection; σ² roughness; resonant drains.

(Full citation list in the session log; primary sources named inline above.)
