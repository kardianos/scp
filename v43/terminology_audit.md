# Terminology Audit: Braid vs Proton

**Date**: 2026-03-26
**Context**: V43 established that braids (single z-aligned helical structures) and
protons (phase-confined 3-braid composites) behave qualitatively differently:

- Braids drift toward HIGH density (EM-dominated, anisotropic)
- Protons drift toward LOW density (gravity-dominated, isotropic)
- Braids are quark analogs; protons are the physical particles

This audit flags every instance where "braid" is used in a way that could mislead
the reader about which object's physics is being described.

---

## Summary of Most Important Fixes

1. **CONCEPT.md Section 3 ("Gravity")** uses "braid" throughout when describing the
   gravitational mechanism. Most of this section (footprint asymmetry, F=ma, force
   law) was written from V33 z-aligned braid experiments. V43 showed the braid
   drifts the WRONG WAY. The gravitational mechanism applies to the *proton*, not
   the braid. The section's V43 caveat (lines 139-144) is good but insufficient --
   the surrounding text still reads as if "braid" = "gravitationally responsive object."

2. **CONCEPT.md Section 10 ("Open Questions")** repeatedly says "braid" where
   "particle" or "proton" would be more accurate, since the questions are about
   physics that the proton exhibits, not the bare braid.

3. **FUTURE.md F2** has already been partially updated with V43 findings but still
   uses "braid" in the mechanism summary (lines 266-272) where the described
   behavior (drift toward depletion = gravity) is WRONG for braids and RIGHT for
   protons.

4. **EM_THEORY.md** is mostly CORRECT -- braids genuinely are the EM objects, and
   most discussion is about single-braid EM properties. However, the Coulomb/charge
   sections (OQ1, OQ2) implicitly assume the charged object exhibiting the Coulomb
   force is a braid, when it should be a composite baryon.

5. **DISCOVERIES.md** is generally well-labeled (V43 entries are accurate). The
   V33 and V34 entries correctly describe braid experiments but do not note the
   V43 caveat that braid gravity measurements are EM-contaminated.

---

## File-by-File Audit

### CONCEPT.md

#### Section 2: "Particles (Braids)" (lines 46-89)

- **Line 46**: "## 2. Particles (Braids) [CONFIRMED]"
  **Classification: B (AMBIGUOUS)**
  The section title equates "particles" with "braids." Post-V43, the physical
  particle is the proton (3-braid composite), not the braid itself. The braid is
  a sub-component (quark analog). Suggest: "## 2. Braids and Particles [CONFIRMED]"
  with a note that braids are the sub-components; baryons are the particles.

- **Lines 49-50**: "The braid is NOT separate from the field. It IS the field..."
  **Classification: C (CORRECT)**
  This describes the braid structure itself, which is accurate.

- **Lines 68-79**: Binding anatomy table and description
  **Classification: C (CORRECT)**
  Describes the single braid's internal structure.

- **Lines 82-88**: Stability discussion
  **Classification: C (CORRECT)**
  Describes single-braid stability experiments (V29).

#### Section 3: "Gravity" (lines 92-377)

- **Lines 96-98**: "The braid depletes the surrounding field... When another braid's
  depletion zone creates a rho gradient, the first braid's spatial perturbation
  profile becomes ASYMMETRIC."
  **Classification: A (WRONG)**
  The gravitational mechanism (drift toward low rho) applies to the proton, not the
  braid. The braid drifts toward HIGH rho.
  **Suggested replacement**: "The particle depletes the surrounding field by binding
  field energy into its structure. When another particle's depletion zone creates a
  rho gradient, the test particle's spatial perturbation profile becomes ASYMMETRIC."
  (Or keep "braid" but add: "This mechanism operates correctly for the proton
  composite; the bare z-aligned braid's response is dominated by EM effects -- see
  the V43 caveat below.")

- **Lines 100-107**: "The braid's perturbation decays as a Yukawa tail... The braid
  therefore REACHES FURTHER into the depleted side than the dense side."
  **Classification: B (AMBIGUOUS)**
  The Yukawa-tail mechanism is generic (applies to both braids and protons), but
  the conclusion that this produces gravitational drift is only correct for protons.

- **Lines 110-116**: Footprint asymmetry measured (V33)
  **Classification: D (NEEDS CONTEXT)**
  The measurement is of a z-aligned braid. The footprint asymmetry is real but the
  net drift direction for a braid is opposite to gravity. Should note: "This
  footprint asymmetry was measured on a z-aligned braid. The bare braid's net drift
  is reversed by EM effects (V43); the proton composite shows the correct
  gravitational response."
  DT: This is in CONCEPT, so we need to describe the concepts as we understand them, NOT individual experiments.

- **Lines 118-121**: "Braid A depletes rho... Braid B's perturbation extends
  further... B's equilibrium shifts toward A -> drift -> GRAVITY"
  **Classification: A (WRONG)**
  This narrative implies braid-braid interaction produces gravity. V43 showed braids
  drift the wrong way. The mechanism produces gravity only at the proton level.
  **Suggested replacement**: Replace "Braid" with "Particle" or "Proton" in this
  diagram, or add explicit caveat that this works for protons but not bare braids.

- **Lines 128-130**: "The braid-field coupling is STRONGER in high-rho backgrounds"
  **Classification: C (CORRECT)**
  This describes a measured property of the braid (V33 drag test), accurately.

- **Lines 132-144**: Force Law section with V43 caveat
  **Classification: D (NEEDS CONTEXT)**
  The caveat at lines 139-144 correctly notes C=186 is contaminated. But the
  preceding text (lines 134-137) presents C=186 as THE gravitational constant
  without qualification until the caveat. Consider moving the caveat earlier or
  restructuring so the proton measurement is primary.
  DT: Correct

- **Lines 146-162**: Force table (D vs DeltaD)
  **Classification: D (NEEDS CONTEXT)**
  This table is from V33 braid-braid experiments. Post-V43, these measurements mix
  gravitational and EM contributions. Should note this is braid-braid (not
  proton-proton) and the values include EM effects.

- **Lines 164-176**: "Gradient Response [CONFIRMED]"
  **Classification: A (WRONG)**
  "A braid placed in an imposed rho gradient... drifts toward LOW rho (V33 gradient
  test)." V43 showed braids drift toward HIGH rho. The V33 result may have been a
  different experimental configuration or a misinterpretation.
  **Suggested replacement**: This paragraph should either be removed (since the
  proton gradient test supersedes it) or clearly labeled as the V33 braid result
  with the V43 correction that the braid's net drift is actually toward HIGH rho
  in the controlled V43 test.

  NOTE: There is a possible contradiction between V33 (braid drifts toward low rho)
  and V43 (braid drifts toward high rho). If the V33 gradient test used different
  conditions, this should be explained. If V43 supersedes V33, the V33 claim should
  be corrected.

- **Lines 178-221**: "Proton Gradient Test [CONFIRMED -- V43]"
  **Classification: C (CORRECT)**
  This section properly distinguishes proton from braid behavior.

- **Lines 264-346**: "F = ma from the Field"
  **Classification: A (WRONG) in multiple places**
  This entire section uses "braid" where "particle" or "proton" is meant:
  - Line 277: "gravitational mass of the SOURCE braid"
  - Line 281: "gravitational mass of the TEST braid"
  - Lines 305-306: "The test braid responds to this gradient"
  - Lines 327-328: "m = the braid's bound energy"
  - Lines 331-340: "Both masses originate from... how much field energy is bound
    into the braid's helical structure"
  - Lines 344-346: "testing requires braids of different sizes"
  The F=ma structure describes gravitational physics, which only works for protons.
  **Suggested replacement**: Replace "braid" with "particle" throughout this section,
  or add a header note: "In this section, 'braid' refers to any solitonic structure
  (braid or composite baryon). The mechanism has been confirmed for the proton
  composite (V43); bare braids are EM-contaminated."

- **Lines 353-364**: "Braid-Field Coupling Scales with rho"
  **Classification: C (CORRECT)**
  Describes a measured braid property (V33 drag test).

- **Lines 366-376**: "Comparison to Newton"
  **Classification: D (NEEDS CONTEXT)**
  The exponent n~1.8 was measured from braid-braid interactions, which are
  EM-contaminated. Should note that the proton-proton force law exponent has not
  yet been measured.

#### Section 4: "Electromagnetism" (lines 380-527)

- **Lines 414-415**: "Sourced by curl(phi) -> helical braids radiate theta"
  **Classification: C (CORRECT)**
  Braids ARE the EM sources.

- **Lines 470-471**: "Winding number IS electric charge for the theta-mediated
  interaction."
  **Classification: D (NEEDS CONTEXT)**
  True for single braids, but V43 clarifies that the physically relevant "charge"
  is the NET chirality of the composite baryon (UUD = +1, UDD = -1), not the
  individual braid's winding.

- **Lines 496-497**: "The net force on two braids is gravity + EM"
  **Classification: A (WRONG)**
  Two braids do not experience clean "gravity + EM." The braid's gravity is
  contaminated by EM. The correct statement is: the net force on two *protons* is
  gravity + EM.
  **Suggested replacement**: "The net force on two particles is gravity + EM, with
  the sign of EM depending on relative winding. For bare braids, the EM component
  can dominate and reverse the net force direction (V43)."

- **Lines 505-507**: "What is the analog of Coulomb's law for static braids?"
  **Classification: B (AMBIGUOUS)**
  The Coulomb law question is really about composite baryons (protons), not bare
  braids. Braids are magnetic dipoles, not Coulomb charges.

#### Section 7: "Dark Matter" (lines 761-773)

- **Line 766**: "V30's M7 analysis showed the braid simultaneously..."
  **Classification: C (CORRECT)**
  Historical reference to specific braid experiments.

- **Lines 771-772**: "The 'dark matter' in this framework would be the field
  depletion around each braid"
  **Classification: B (AMBIGUOUS)**
  Should say "around each particle" or "around each baryon" since dark matter
  halos are macroscopic and the relevant depleting object is the proton, not the
  sub-component braid.

#### Section 8: "Cosmology" (lines 776-788)

- **Lines 780-787**: Uses "braid" throughout
  **Classification: C (CORRECT)**
  Discusses braid formation specifically.

#### Section 9: "Lorentz Invariance" (lines 791-879)

- **Lines 810-817**: "Lorentz Contraction of Braids"
  **Classification: B (AMBIGUOUS)**
  The test is about boosting a soliton. Could be a braid or a proton. The title
  should probably say "of Solitons" or "of Particles."

- **Lines 820-855**: "The Depletion Is Power-Law, Not Yukawa"
  **Classification: D (NEEDS CONTEXT)**
  "The braid's COLLECTIVE depletion profile" -- this was measured on a z-aligned
  braid. Should note that the proton's depletion profile may differ (more isotropic).

- **Line 834**: "The braid's depletion extends as a power law"
  **Classification: B (AMBIGUOUS)**
  Same as above; "particle's depletion" would be more general.

#### Section 10: "Open Questions" (lines 883-968)

- **Line 889**: "static electric force between stationary braids of opposite winding"
  **Classification: A (WRONG)**
  The Coulomb force question is about baryons/protons, not bare braids. Braids are
  magnetic dipoles; the "electric charge" emerges at the composite level.
  **Suggested replacement**: "static electric force between stationary composite
  baryons of opposite net chirality"

- **Lines 892-894**: "Does the theta shell around the braid support bound
  standing-wave modes? ... like an electron shell."
  **Classification: B (AMBIGUOUS)**
  The electron would orbit a proton, not a bare braid.

- **Line 895**: "Multi-braid formation: conditions for stellar nucleosynthesis analog."
  **Classification: B (AMBIGUOUS)**
  Nucleosynthesis forms nuclei from baryons, not bare braids.

#### Section 11: "Resolved Questions" (lines 902-968)

- **Lines 904-905**: "E(D) is monotonically repulsive" (about two braids)
  **Classification: C (CORRECT)**
  Correctly describes the braid-braid experiment.

- **Lines 928-930**: "Same-winding braids attract 27% more, opposite-winding
  attract 57% less"
  **Classification: C (CORRECT)**
  Correctly describes V34 braid experiments.

---

### DISCOVERIES.md

#### V43 Section (lines 1-63)

- **Classification: C (CORRECT)** throughout.
  This section properly distinguishes protons from braids and documents the
  dichotomy. No changes needed.

#### V42 Section (lines 66-95)

- **Classification: C (CORRECT)** throughout.
  Uses "baryons" correctly. No braid/proton confusion.

#### V41 Section (lines 97-117)

- **Line 100**: "Braids attracted at medium range by depletion"
  **Classification: C (CORRECT)**
  Describes braids within the composite (sub-component attraction).

#### V34 Section (lines 213-298)

- **Lines 216-218**: "Same-winding braids attract 27% more... Winding number IS
  electric charge."
  **Classification: D (NEEDS CONTEXT)**
  True for braids, but should note that the physically relevant "charge" is the
  composite baryon's net chirality (V43). The braid's winding is the quark-level
  charge analog, not the particle-level charge.

- **Lines 232-235**: Two-braid force comparison table
  **Classification: C (CORRECT)**
  Accurately describes braid-braid V34 experiments.

- **Lines 251-252**: "Not yet confirmed: Whether theta creates a 1/r^2 field,
  whether it mediates charge-dependent forces, whether W=+1 braids produce
  opposite theta patterns."
  **Classification: C (CORRECT)**
  Historical (later confirmed). Correctly references braids.

#### V34 Depletion Section (lines 254-265)

- **Line 257**: "the braid's depletion profile decays as delta-rho ~ 1/r^1.2"
  **Classification: D (NEEDS CONTEXT)**
  Measured on z-aligned braid. Proton depletion profile may differ.

#### V33 Section (lines 267-283)

- **Line 269**: "The gravitational force is F = -C * nabla-rho"
  **Classification: D (NEEDS CONTEXT)**
  Measured using braids. V43 showed this is EM-contaminated for braids. Should
  note that C=186 includes EM effects and C_proton has not been precisely measured.

- **Line 281**: "The braid's footprint is asymmetric in a density gradient"
  **Classification: C (CORRECT)**
  Correctly describes the braid measurement.

#### V32-V33 and earlier sections

- Generally **C (CORRECT)** -- these are historical records of braid experiments
  and correctly use "braid" to describe what was actually tested.

#### V27-M3c (line 848-857)

- **Lines 854-856**: "This is gravity-like (universal attraction regardless of
  'charge'), not EM-like"
  **Classification: D (NEEDS CONTEXT)**
  Two braids attracting is not cleanly "gravity-like" per V43 -- braids have
  EM-contaminated interactions. The statement is from pre-V43 understanding.

---

### EM_THEORY.md

#### Section 1: "The Field Mapping" (lines 1-83)

- **Line 7**: "Braids act as Amperian current loops"
  **Classification: C (CORRECT)**
  Braids genuinely are the current loops.

- **Lines 40-43**: "A braid's helical twist creates a nonzero nabla-cross-phi"
  **Classification: C (CORRECT)**

- **Lines 58-63**: "The braid is fundamentally a magnetic object... Effective
  electric charge emerges at the composite-baryon level"
  **Classification: C (CORRECT)**
  This already makes the correct distinction.

- **Lines 65-66**: "Same-winding braids attract, opposite-winding repel."
  **Classification: C (CORRECT)**
  Describes braid-braid EM interaction.

#### Section 2: "Light as a Coupled Wave" (lines 86-173)

- **Line 107**: "not a particle -- particles are braids"
  **Classification: A (WRONG)**
  Particles are protons (3-braid composites), not individual braids. Braids are
  quark analogs.
  **Suggested replacement**: "not a particle -- particles are braids and their
  composites (baryons)"

- **Lines 121-122**: "Inside or very close to braids (where phi is large)"
  **Classification: C (CORRECT)**
  Describes the braid's local field structure.

- **Lines 126**: "Photons propagate at exactly c through empty space. The phi
  coupling only affects light when it passes through or near matter (braids)."
  **Classification: B (AMBIGUOUS)**
  "Matter (braids)" should be "matter (braids and their composites)" since the
  reader might think braids are standalone particles.

#### Section 3: "Spatial Force Structure" (lines 176-218)

- **Lines 180-193**: Force decomposition around "a single braid"
  **Classification: C (CORRECT)**
  Describes single-braid force structure from V34 hires data.

- **Lines 202-217**: "Connection to V42 Deuterium"
  **Classification: C (CORRECT)**
  Correctly uses "baryon" and "braid" distinctly.

#### Section 4: "The Meta-Static Coulomb Potential" (lines 221-328)

- **Lines 248-253**: "Each braid emits oscillating theta waves... When theta-waves
  from braid A hit braid B's current loop..."
  **Classification: D (NEEDS CONTEXT)**
  The radiation pressure mechanism describes braid-braid EM interaction. The
  Coulomb force between physical particles (protons/baryons) would involve the
  composite radiation pattern (OQ1), not individual braid dipole radiation.

- **Lines 300-302**: "Braid geometry (R_tube, ellipticity, amplitude)"
  **Classification: C (CORRECT)**
  Describing parameters of the braid structure.

#### Section 6: "Ohm's Law and Conductivity" (lines 413-443)

- **Lines 417-418**: "A 'wire' in this theory would be a lattice of braids (same
  winding) in a line or tube."
  **Classification: A (WRONG)**
  A wire would be a lattice of protons/baryons with orbital electrons, not bare
  braids. Braids are confined inside baryons.
  **Suggested replacement**: "A 'wire' in this theory would be a lattice of baryons
  with mobile electron analogs (theta perturbations in orbital modes)."

- **Lines 427-428**: "~10-100 braids in a line geometry"
  **Classification: A (WRONG)**
  Same issue. Should be "~10-100 baryons."
  **Suggested replacement**: "~10-100 baryons in a line geometry"

- **Lines 430**: "braid drift velocity vs applied theta gradient"
  **Classification: A (WRONG)**
  Should be "charge carrier drift velocity." The current carriers would be electron
  analogs (theta perturbations), not the braid/baryon lattice.
  **Suggested replacement**: "charge carrier drift velocity vs applied theta gradient"

- **Lines 434-435**: "the proportionality constant is the conductivity sigma. The
  resistance R would depend on the braid spacing, amplitude, and coupling eta."
  **Classification: B (AMBIGUOUS)**
  "Braid spacing" should be "baryon spacing" or "lattice spacing."

- **Lines 439-440**: "the simulation scale (100+ braids) requires significantly
  larger grids"
  **Classification: A (WRONG)**
  Should be "100+ baryons."
  **Suggested replacement**: "the simulation scale (100+ baryons)"

#### Open Questions (lines 444-589)

- **Lines 447-449**: "Each braid is a magnetic dipole... A UUD baryon has three
  ORTHOGONAL current loops."
  **Classification: C (CORRECT)**
  Properly distinguishes braid (dipole) from baryon (composite).

- **Lines 467-469**: "the scattering cross-section of a target braid and the
  time-averaged Poynting flux from the source braid"
  **Classification: D (NEEDS CONTEXT)**
  The ponderomotive force calculation is for braid-braid EM interaction. The
  physical Coulomb force is between composite baryons. Should clarify whether
  this calculation is for the sub-component or the composite.

- **Lines 535-536**: "the baseline is 28% (from braid_hires.sfa shell_analysis)"
  **Classification: C (CORRECT)**
  References specific braid data.

#### Next Steps (lines 569-589)

- **Line 577**: "single-braid runs"
  **Classification: C (CORRECT)**
  Describes actual single-braid simulations.

- **Line 586**: "Wire/Ohm's law simulation (100+ braids, N=1024+)"
  **Classification: A (WRONG)**
  Should be "100+ baryons."
  **Suggested replacement**: "Wire/Ohm's law simulation (100+ baryons, N=1024+)"

---

### FUTURE.md

#### F17: Nuclear Binding Energy (lines 14-51)

- **Line 38**: "Each baryon: phase-confined 3-braid"
  **Classification: C (CORRECT)**

#### F1: Mass Parameter and Long-Range Gravity (lines 183-210)

- **Line 193**: "The braid's collective perturbation of the background propagates
  as a PHONON"
  **Classification: B (AMBIGUOUS)**
  This is probably true for both braids and protons, but the gravitational
  significance applies to protons.

- **Lines 195-199**: "The m^2 parameter confines the braid... m^2 >= 1.25 needed
  for braid binding."
  **Classification: C (CORRECT)**
  Describes braid confinement specifically.

#### F2: Gravity Mechanism (lines 212-279)

- **Line 214**: "The V33 braid measurement established F proportional to nabla-rho"
  **Classification: D (NEEDS CONTEXT)**
  Already has V43 context in the same paragraph, but "established" is too strong
  for a measurement now known to be EM-contaminated.

- **Lines 230-243**: Footprint asymmetry, drag test
  **Classification: C (CORRECT)**
  Correctly describes braid experiments as braid experiments.

- **Lines 266-272**: Mechanism summary bullets
  **Classification: A (WRONG)**
  "The braid's oscillation cycle creates an asymmetric spatial perturbation...
  The force direction is toward depletion = toward other braids = GRAVITY"
  This is wrong for braids (they drift toward HIGH rho). The mechanism described
  here is correct for PROTONS.
  **Suggested replacement**: Replace "braid" with "particle" in these bullets, and
  note that the mechanism has been confirmed for the proton composite (V43). Add:
  "For bare z-aligned braids, the EM anisotropy reverses the net drift direction."

- **Lines 276**: "Connection between C and braid parameters (m, mu, kappa, braid
  size)"
  **Classification: B (AMBIGUOUS)**
  Should be "particle parameters" since the gravitational C is for the proton.

#### F3: Lorentz Contraction (lines 284-299)

- **Line 287**: "Boost a braid at v=0.1c, 0.3c, 0.5c"
  **Classification: B (AMBIGUOUS)**
  Could test either a braid or a proton. Should specify.

#### F4: Isotropic Background (lines 301-311)

- **Line 304**: "Do braids survive?"
  **Classification: C (CORRECT)**
  Asking about braids specifically (the sub-component).

#### F5: Equivalence Principle (lines 313-327)

- **Lines 321-325**: "Does it create a depletion zone? Does a braid respond to it
  gravitationally? If NO: only 'matter' (braids) gravitates"
  **Classification: A (WRONG)**
  The gravitationally responsive object is the proton, not the braid.
  **Suggested replacement**: "Does a proton respond to it gravitationally? If NO:
  only 'matter' (baryons) gravitates"

#### F6: Electromagnetism (lines 332-362)

- **Line 336**: "Massless theta sourced by the braid's helical twist"
  **Classification: C (CORRECT)**

- **Lines 340-341**: "Same-winding braids attract 27% more... Opposite-winding
  braids attract 57% less."
  **Classification: C (CORRECT)**
  Describes V34 braid experiments.

- **Line 357**: "Is there a static (Coulomb) regime for non-moving braids?"
  **Classification: B (AMBIGUOUS)**
  The Coulomb regime would be between composite baryons, not bare braids.

#### F7: Multi-Braid Formation (lines 364-379)

- **Line 366**: "Heavier 'atoms' are braids sharing helical structure"
  **Classification: A (WRONG)**
  Heavier atoms are multi-baryon nuclei, not multi-braid structures.
  **Suggested replacement**: "Heavier nuclei are multi-baryon composites, each
  baryon being a phase-confined 3-braid structure."

- **Lines 375-376**: "multi-braid bound states... braid condensation"
  **Classification: B (AMBIGUOUS)**
  "Multi-braid" could mean multi-component composites or multi-baryon nuclei.
  The context (nucleosynthesis) suggests multi-baryon.

#### F8: Binding Energy Accounting (lines 381-390)

- **Lines 383-388**: "When a braid forms, does the total field energy DECREASE?...
  Compare E_total of (braid + depleted background)"
  **Classification: C (CORRECT)**
  Asking about individual braid formation energy.

#### F9: Analytical Effective Potential (lines 429-448)

- **Lines 437-441**: "braid-braid interaction potential... second braid as a test
  perturbation in the first braid's depletion field"
  **Classification: D (NEEDS CONTEXT)**
  The braid-braid interaction includes EM contamination. Should note whether this
  analytical derivation targets the proton-proton or braid-braid potential.

#### F10: Spin and Helical Handedness (lines 452-461)

- **Line 455**: "The braid has a helical twist"
  **Classification: C (CORRECT)**

#### F11: Gravitational Waves (lines 463-469)

- **Lines 464-465**: "Does an accelerated braid produce radiation with spin-2..."
  **Classification: B (AMBIGUOUS)**
  Gravitational radiation would come from an accelerated proton/baryon, not a bare
  braid. The braid's radiation is predominantly EM (magnetic dipole).

- **Line 468**: "Accelerate a braid (give it a kick)"
  **Classification: B (AMBIGUOUS)**
  Should specify proton or braid; GW signature depends on which.

#### F12: Dark Matter Profiles (lines 471-475)

- **Line 473**: "Does the depletion profile around a braid match observed dark
  matter halo profiles"
  **Classification: B (AMBIGUOUS)**
  Dark matter halos surround galaxies worth of baryons, not individual braids.
  Should say "around a particle" or "around an accumulation of baryons."

#### F13: Quantization (lines 477-482)

- **Line 481**: "Are the braids stable quantum mechanically?"
  **Classification: C (CORRECT)**
  Asking about braid stability specifically.

#### X7: Density-Dependent kappa (lines 539-556)

- **Lines 546-549**: "Phi/c^2 ~ -0.02 per braid... ~25 braids concentrated...
  N ~ 10^19 braids"
  **Classification: B (AMBIGUOUS)**
  The gravitational potential per "braid" should be per "baryon" or "proton" since
  these are the gravitationally relevant objects. The Planck mass calculation should
  reference protons.

#### Resolved Composite Questions (lines 565-574)

- **Line 566**: "UDD 3-braid at R=4 survived T=200"
  **Classification: C (CORRECT)**

- **Line 569**: "Can braids be confined without merging?"
  **Classification: C (CORRECT)**
  Asks about braids within a composite.

---

## Summary Statistics

| Classification | Count | Description |
|---------------|-------|-------------|
| A (WRONG) | 14 | Says "braid" but physics is specifically proton/baryon |
| B (AMBIGUOUS) | 18 | Could mean either; would benefit from clarification |
| C (CORRECT) | 35 | Correctly refers to the z-aligned helical sub-component |
| D (NEEDS CONTEXT) | 12 | Correct but should note V43 braid/proton caveat |

The most impactful changes are the A-classified items in CONCEPT.md Section 3
(gravity mechanism described using "braid" when it only works for protons),
EM_THEORY.md Section 6 (Ohm's law using "braids" instead of "baryons"),
and FUTURE.md F2 mechanism summary and F7 nucleosynthesis.
