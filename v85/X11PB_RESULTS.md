# X11-PB results — the two-closure orbit ("fabric hydrogen-lite"), completed T=2000
**Date:** 2026-07-26 · gpu7 (V100-32GB), 19.7 h, clean completion; 22 frames
@ 100 t.u. (frames t≥1800 unreadable — async-writer tail corruption; diag @
0.25 t.u. covers the full run; trajectory intact through t=1700).
Data: `/space/scp/v85/topo1/gpu/x11pb*`.

## Design
Heavy: v71 flavored ball (Q=+284, partition 76.7/103.7/103.7, clocks
1.38/1.42/1.42) at origin. Light: branch-bottom soliton (Q=−88.5, ω=−1.475,
E=135.5, r_half=2.35) at (0,12,0), tangential v=0.064 (predicted circular-
orbit speed; predicted period ≈1200 t.u.). Protections: separation ≫ tails;
clock detuning beat ≈2.9; flavor-partition mismatch. N=384, L=55, g=0.05.

## The measured arc (negative-charge centroid trajectory, cpos)
| phase | t | r | Q_B | what happened |
|---|---|---|---|---|
| Orbit | 0–340 | 12.1 → 16.3 | −87.5 → −71 | **First two-closure orbit in the program**: elliptical, planar, azimuthal rate matching L-conservation; period from early sweep ≈1270 vs 1200 predicted (6%) |
| Runaway evaporation | 300–500 | 16.3 → 11.2 | −71 → −42 | Tidal + self agitation drove accelerating erosion of the evaporation-metastable satellite (rate doubling ~150 t.u.) — the sustained-agitation runaway SURV could not produce in isolation |
| Deorbit | 400–600 | 15.4 → 6.0 | −42 → −37 | Mass-loss retro-rocket: shed charge carried orbital energy; satellite plunged into the heavy ball's skirt |
| Contact burn | 600–900 | 3.4–6.5 | −37 → −29 | **Annihilation throttled**: −2 to −3 charge/100 t.u. at full overlap (pos1-class slowness; detuning + flavor mismatch); heavy ball lost only ~6 units total |
| Re-ascent + burn-orbit | 900–2000 | 5.6–9.2 (cycling) | −29 → −15 | Burn weakens attraction → remnant climbs; ascent slows burn → shallow cycling around the skirt. Ends as a **bound, orbiting, slowly-burning remnant** (burn ≈ −1/100 and easing) |

Ledger: Gauss at 1.9e-13 floor throughout; θ silent; E drift −21% = the
evaporation/annihilation radiation leaving through the sponge. Heavy ball
intact at +263-ish equivalent (lost ≈6 of 284 across the entire contact era).

## Verdicts
1. **Two-closure Coulomb orbit: ACHIEVED** — the Stage-3 Step-2 target that
   failed throughout v75–v83. Kepler-consistent period to 6%, planar,
   L-conserving.
2. **Carrier fragility is the binding constraint.** The only available light
   carrier (branch-bottom ball) is evaporation-metastable; sustained orbital
   agitation triggers runaway shedding. A fabric atom made from same-branch
   carriers self-destructs on ~10² t.u. — sharpening HIER: the theory needs
   a light *stable* carrier, which one branch cannot provide.
3. **Annihilation protections work as engineered**: contact burn throttled
   ~50× below naive overlap rates (compare pos1); flavor + detuning did
   their job. The heavy partner survives contact with a conjugate object.
4. **The end-state is a new object**: a throttled binary — heavy flavored
   ball + orbiting sub-threshold remnant cycling at the skirt in a
   burn-ascend-slow feedback loop.
5. **Companion to SURV** (see SURV_RESULTS.md): together they establish
   that sustained agitation, not charge, drives erosion — the satellite's
   runaway under tides vs isolated balls' quiet simmer.

## Atom-arc synthesis (X10 pathfinder → X10b → X10c → X10d → X11-PB)
The classical fabric atom now has a complete phenomenology: response-mode
shells at predicted radii (retained, breathing, monopole-protected); cascade
under radiation reaction; ground-state access contact-limited; smeared
clouds cannot orbit (carriers required); carrier orbits achieved but
carrier lifetime bounded by agitated evaporation. The architecture's missing
piece is unchanged and now doubly confirmed from dynamics: **a light,
stable, opposite-charge carrier** — the HIER lock, which one branch cannot
supply and which motivates either the second-sector design (TOE §D.2) or
acceptance of cloud-atoms (retention without orbits) as the fabric's honest
atom.
