# V38 Research: The Braid-as-Electron Hypothesis

## Executive Summary

This document investigates the speculation (SPECULATION.md) that the helical
braid soliton in the SCP field theory should be identified as an ELECTRON
analog, not a proton. Five research questions are addressed in depth. The
overall verdict is: **the hypothesis is interesting and structurally motivated,
but currently faces severe quantitative difficulties on three fronts** (theta
field decay law, scale hierarchy, and bound state formation).

---

## 1. Is the Braid's Theta Field Truly EM-Like?

### What V34 Actually Measured

The V34 Cosserat theta characterization (`v34/torsion_coupling/theta_characterize/`)
performed three experiments on the 6-field system:

    d^2 phi_a/dt^2 = lap(phi_a) - m^2 phi_a - V'(P) + eta * curl(theta)_a
    d^2 theta_a/dt^2 = lap(theta_a)                  + eta * curl(phi_a)

with eta=0.5, m_theta=0, N=80, L=25.

**Experiment 1 (radial decay)**: The azimuthal theta component theta_phi^2
decays as r^(-0.50), NOT as r^(-2) (Biot-Savart). The raw data
(`theta_characterize/data/theta_radial.tsv`) shows theta_phi^2 oscillating
between 0.9e-3 and 3.8e-3 across r=1-20, with no clean power law. The
azimuthal fraction theta_phi^2/theta_total^2 is only 0.45-0.58, far from
the >95% expected for a true magnetic field.

**Experiment 2 (charge conjugation)**: The time-averaged DC component of
theta_phi reverses sign with winding (correlation -0.68), but this DC signal
is only 0.2% of the oscillation amplitude. The dominant theta response is
a standing/propagating wave with period ~4 time units, whose phase is
independent of winding sign.

**Experiment 3 (two-braid force)**: This is the STRONGEST result.
Same-winding braids attract 23-27% more than gravity-only baseline;
opposite-winding attract 36-57% less. This charge-dependent force is
confirmed and dt-converged.

### Assessment: Does It Match Biot-Savart?

**No.** The theta field is NOT a static 1/r magnetic field. Specifically:

| Property | Classical EM (B around wire) | Observed theta field |
|----------|---------------------------|---------------------|
| Radial decay | B ~ 1/r | theta_phi ~ r^(-0.25) (amplitude) |
| Angular structure | >99% azimuthal | 45-58% azimuthal |
| Time dependence | Static (DC) | Oscillating (period ~4, DC = 0.2%) |
| Force mechanism | Magnetostatic | Wave-mediated (radiation coupling) |

The theta field is better described as a **wave-exchange mediator**, closer
to QFT photon exchange than to classical Biot-Savart. The 0.2% DC bias
carries the charge information, but the 99.8% oscillation carries the
energy and mediates the force.

### What This Means for the Electron Hypothesis

If the braid is an electron, its "electromagnetic field" is fundamentally
quantum-field-theoretic from the start: photon exchange, not static fields.
This is actually MORE correct than classical EM for a fundamental particle.
However, the r^(-0.5) envelope (not r^(-1) or r^(-2)) means the force law
is NOT Coulomb at the tested separations. At longer range, the massless
theta waves propagate at c, so the asymptotic behavior may improve. This
has not been tested (would require N=256+, L=100+, T=500+).

**Verdict**: The charge-dependent force is real and has the correct sign
structure for EM. The field structure is wave-mediated, not static. The
radial decay does not match Coulomb/Biot-Savart at the tested scales.
Further investigation at longer range is needed.

---

## 2. What Determines the Braid's "Charge"?

### Winding Number as Topological Charge

The braid's charge is its **phase winding number** W:

    W = (1/2pi) sum_{k=0}^{N-1} Delta(theta_k)
    theta_k = atan2(phi_1(z_k), phi_0(z_k))

V27-M5 proved this is EXACTLY conserved (to machine precision) for all
time. The mechanism is topological: periodic BC + smooth dynamics make
W an integer that cannot change under continuous evolution. The Jacobian
density is a total divergence whose integral is invariant.

### Is It Integer or Continuous?

**Strictly integer.** W = -1 for all tested braids at the bimodal sweet
spot. V29-T1b shows W = -1.000 at every timestep for T=500+. The winding
is quantized by the periodic BC in z: the phase must return to itself
after one traversal, so the total winding must be an integer multiple of 2pi.

This is EXACTLY analogous to the quantization of electric charge in
electrodynamics. The integer nature is topological (homotopy class
pi_1(S^1) = Z), not dynamical.

### Charge Conservation in Interactions

V27-M3c tested two-braid interactions:
- Same twist (W=+1, W=+1): both survive, Delta_sep = -7.2 (attract)
- Opposite twist (W=+1, W=-1): both survive, Delta_sep = -5.0 (attract)

Total winding is conserved in both cases. The braids do NOT annihilate
or change winding during the tested interactions (T=200). This is
consistent with charge conservation.

However, there is a key limitation: the periodic BC in z ENFORCES
winding conservation for each individual braid. In a fully 3D simulation
where braids are finite-length tubes (not periodic), winding could in
principle unwind at the endpoints. This has not been tested.

### Implications for the Electron Hypothesis

The integer topological charge (+/-1) maps naturally to electron charge
(+/-e). The conservation law is exact and topological. The sign determines
the theta-mediated force direction (same-sign attract via theta, opposite
repel), which is the correct EM sign convention for current-current
interactions (parallel currents attract).

**Verdict**: Winding number is a robust integer topological charge,
conserved under all tested dynamics. It is the natural candidate for
electric charge in this theory.

---

## 3. Scale Hierarchy: What Would the Proton Mass Be?

### Parameter Fitting (from V24 and MEMORY.md)

The V24 parameter fitting established:
- 1 code energy unit = 9.098 MeV
- 1 code length unit = 0.5624 fm
- 1 code time unit = 1.875e-24 s
- Physical constraint: M_p * r_p / (hbar*c) = 4.00

These were derived by fitting the Skyrme model to the proton (M_p = 938 MeV,
r_p = 0.84 fm). But the Skyrme model identifies the Skyrmion soliton
(the 3D topological object) with the proton. The braid is a DIFFERENT
kind of object.

### Braid Energy in Code Units

From V33-C2 (N=128, L=20, single braid):
- E_total(t=0) = 5357 code units
- E_total(t=500) = 5334 code units (after radiation loss)

From V27-M4 (m=0 braid, N=128):
- E_total = 7915 code units (massless regime, higher)

The V33 braid (m^2=2.25, standard parameters) has E ~ 5350 code units.

### If the Braid Is an Electron

Using the V24 conversion: E_braid = 5350 * 9.098 MeV = 48,674 MeV.

This is 95,000 times the electron mass (0.511 MeV). The braid is FAR
too heavy to be an electron with the V24 conversion factors.

**But the V24 conversion was fit to the proton.** If we instead fit to
the electron:
- m_e = 0.511 MeV, E_braid = 5350 code units
- => 1 code E = 0.511/5350 = 9.55e-5 MeV = 95.5 eV
- => 1 code L = hbar*c / (m_e * c^2 * r_factor) ... this requires a
  different length calibration

The classical electron radius is r_e = 2.82 fm. The braid has spatial
extent ~10 code units. So:
- 1 code L = 2.82 fm / 10 = 0.282 fm (if braid radius ~ r_e)
- OR 1 code L = 0.282 fm / some factor

But this doesn't work cleanly because:
1. The braid is a 1D tube, not a point particle
2. The classical electron radius is not a "size" in QM
3. The Compton wavelength (386 fm) is much larger than the braid

### What Would the Proton Mass Be?

If the braid is the electron (0.511 MeV), the proton must be a composite
at mass 938 MeV = 1836 * m_e.

In the SCP theory, the proton would need to be:
- A 3D composite object (three braids bound together, or a Skyrmion-like
  structure in the phi fields)
- 1836 times heavier than the braid
- Spatially SMALLER than the braid (proton r_p = 0.84 fm << electron
  Compton wavelength 386 fm)

The size relationship is the key problem. The braid tube has radius ~2
code units and length ~10 code units. If 1 code unit ~ 0.282 fm, the
braid is 2.8 fm long and 0.56 fm wide. A proton-sized object (0.84 fm)
would need to fit INSIDE the braid tube or at a comparable scale.

### Resolution Requirements

To see proton-scale structure (if it exists) at the electron-braid scale:
- Current: N=128, L=15, dx=0.24 code units
- Braid core: ~2 code units wide = 8 grid cells
- Proton (if at ~0.84 fm, using electron-scale calibration):
  ~3 code units = 12 grid cells
- This is actually RESOLVABLE at current N=128

The problem is not resolution but ENERGY: the proton needs 1836x the
braid energy concentrated in a similar or smaller volume. The triple-product
coupling V(P) = (mu/2)P^2/(1+kappa*P^2) saturates at large P due to the
kappa regularization. Getting 1836x the energy would require either:
- Much larger field amplitudes (limited by kappa saturation)
- Much higher frequency oscillations (limited by CFL condition)
- A fundamentally different binding mechanism

### Assessment

The scale hierarchy problem is SEVERE. The braid energy (5350 code units)
maps to either:
- 48,674 MeV with V24 proton-based calibration (too heavy for electron)
- 0.511 MeV with electron-based calibration (then proton requires 1836x
  energy concentration not demonstrated)

Neither calibration produces a self-consistent two-particle spectrum.

**Verdict**: The theory currently has no mechanism to produce the observed
electron-proton mass hierarchy (1:1836) from braid dynamics alone. This
is the biggest obstacle to the hypothesis.

---

## 4. Multi-Braid Bound States

### V23 Results: 1D Oscillon Pairs

V23 found inter-oscillon equilibrium at D_eq ~ 15.7 in 1D, with:
- Repulsive core (D < 15.7)
- Attractive well (D > 15.7), depth ~ 4e-4 code energy units
- Spring constant K ~ 1e-3
- Phonon speed c_s = 0.40 +/- 0.02

But the 1D lattice was UNSTABLE (well too shallow by 250x), and 3D
pairs MERGED (no repulsive barrier in 3D due to geometric dilution of
radiation pressure).

### V33 Results: 3D Braid Pairs

V33-C1 measured the force between two braids in 3D:

| D | Force direction | Regime |
|---|----------------|--------|
| 5 | Repulsive | Hard core overlap |
| 8 | Repulsive | Core overlap |
| 12 | Attractive | Onset |
| 15 | STRONG attraction | Peak force |
| 18 | Strong attraction | |
| 20 | Weak attraction | |
| 40 | Detected | Long-range tail |
| 80 | Detected | Very long range |

Power law: F ~ 1/D^1.8 (D=15-30).

V33-C3 (two braids at D=20, T=500) showed SCATTER: approach to D=14,
bounce, escape. No stable bound orbit.

### Can Two Braids Form a "Meson"?

For a bound orbit, the attractive force must exceed the kinetic energy
of the braids. The measured force F ~ -2.4/D^1.8 gives:

    V(D) = -integral F dD ~ -2.4/(0.8 * D^0.8) ~ -3.0/D^0.8

For a circular orbit at radius D:
    F_cent = m*v^2/D = m*omega^2*D
    F_grav = -2.4/D^1.8

Equating: m*omega^2*D = 2.4/D^1.8, so omega^2 = 2.4/(m*D^2.8)

The orbital period T_orb = 2*pi/omega ~ D^1.4 * sqrt(m/2.4).

With the braid mass M ~ 5350 code units:
    At D=15: omega^2 = 2.4/(5350 * 15^2.8) ~ 2.4/(5350*2946) ~ 1.5e-7
    omega ~ 3.9e-4, T_orb ~ 16,000 time units

This is extremely slow. The orbital velocity at D=15:
    v = omega*D = 3.9e-4 * 15 = 0.006 code units/time = 0.006c

At v << c, the orbit would take T ~ 16,000 time units to complete.
V33-C3 ran for T=500 and saw scatter, consistent with the orbital
period being much longer than the simulation time.

### The Theta Force Changes Everything

With the 6-field Cosserat system, same-winding braids attract 27% MORE
than gravity alone. For opposite-winding braids, attraction is 57% LESS.

A "meson" analog would be two OPPOSITE-winding braids:
- Gravitational attraction: F_grav ~ -2.4/D^1.8
- Theta repulsion (opposite winding): partially cancels, net ~ -1.0/D^1.8
- Much weaker binding, even slower orbit

A "baryon" analog would be three SAME-winding braids:
- Total gravitational attraction: enhanced by theta (3 pairs, all same)
- More bound than two-braid system
- But requires 3-body simulation not yet performed

### What About the V23 Equilibrium?

The V23 equilibrium at D_eq=15.7 was for 1D oscillons, not 3D braids.
For braids:
- There IS a repulsive core at D < 12 (V33-C1 data)
- There IS an attractive regime at D > 12
- But V33-C3 shows SCATTER, not equilibrium

The key difference: 1D oscillons have radiation pressure creating a
repulsive barrier at D < D_eq. In 3D, radiation pressure dilutes as
1/r^2 and cannot balance the attractive force. The 3D braids have a
hard-core repulsion from field overlap (D < 12), but this is too steep
to create a smooth equilibrium.

**Verdict**: Multi-braid bound states are plausible in principle (the
force law supports orbits at D ~ 15, T_orb ~ 16,000), but have NOT been
demonstrated. The required simulation times (T > 10,000) far exceed
current runs (T ~ 500). The 6-field theta force provides charge-dependent
corrections but does not qualitatively change the picture.

---

## 5. Literature: Models with Coexisting 1D and 3D Solitons

### The Skyrme-Faddeev-Niemi Connection

The most relevant literature comes from the Skyrme-Faddeev-Niemi (SFN)
model, which is the CLOSEST existing framework to what SCP proposes:

**The Skyrme model** (SU(2) chiral field, pi_3(SU(2)) = Z):
- 3D topological solitons = baryons (protons, neutrons)
- Baryon number B is the topological charge
- Well-established in nuclear physics since 1961

**The Faddeev-Niemi model** (CP^1 or S^2 target, pi_3(S^2) = Z):
- Knot/link solitons = GLUEBALLS (proposed identification)
- Hopf charge Q classifies the knots (Q=1 ring, Q=7 first true knot)
- The solitons are 1D loops/knots, not 3D balls
- Key paper: Faddeev & Niemi, Nature 387 (1997) 58

Cho, Maison, and others showed that the Skyrme theory contains BOTH types:
- Skyrmions (3D): baryon number from pi_3(SU(2))
- Baby Skyrmions (2D vortex sheets)
- Faddeev-Niemi knots (1D loops): from pi_3(S^2) via Hopf map
- All exist as solitons of the SAME SU(2) field theory

Reference: Cho, "New topological structures of Skyrme theory: baryon
number and monopole number," Eur. Phys. J. C 77 (2017) 655.

**THIS IS EXACTLY THE STRUCTURE SCP PROPOSES**: both 1D extended objects
(braids) and 3D compact objects (oscillons/Skyrmions) in the same field
theory. However, in the established literature:
- The 3D objects = baryons (nucleons)
- The 1D objects = glueballs or magnetic vortex rings
- There is NO standard identification of 1D solitons with electrons

### Witten Superconducting Strings

Witten (1985) showed that cosmic strings can carry fermion zero modes
in their cores, making the string a 1D superconductor. The fermions
propagate along the string as massless excitations. Key features:
- The string is a topological defect (1D)
- Fermionic currents flow along it
- The current is quantized and conserved
- Stable loops (vortons) can form

This is structurally similar to the SCP braid: a 1D tube carrying a
conserved current (winding number). However, in Witten's model the
fermions are ADDITIONAL fields trapped on the string, not the string
itself.

### Bosonization / Soliton-Fermion Duality

The Coleman-Mandelstam correspondence (1975) shows that in 1+1D, the
fermion of the massive Thirring model IS the soliton of the sine-Gordon
model. This is the closest theoretical precedent for "a soliton IS an
electron":
- The soliton is a bosonic, extended object
- After quantization, it behaves as a fermion
- The soliton's topological charge = fermion number

In 2+1D and 3+1D, analogous bosonization dualities exist:
- Particle-vortex duality (2+1D): vortex <-> fermion via Chern-Simons
- String-particle duality: conjectured but not rigorously established

### Kelvin Waves on Vortex Lines

In superfluid helium, quantized vortex lines support Kelvin waves --
helical disturbances that propagate along the vortex core. These are
structurally identical to the SCP braid: a helical wave traveling along
a 1D tube. Key papers:
- Kelvin wave cascade on quantum vortex, Kozik & Svistunov, PRL 92 (2004)
- Direct excitation of Kelvin waves, Nature Physics (2024)

The analogy: the SCP background field plays the role of the superfluid,
and the braid plays the role of a Kelvin wave packet on a vortex line.
However, in superfluids, Kelvin waves are low-energy excitations, not
particles. There is no identification of Kelvin waves with electrons.

### Assessment of Literature Support

| Model | 1D Object | 3D Object | 1D = Electron? |
|-------|-----------|-----------|---------------|
| Skyrme + Faddeev-Niemi | Knot soliton (glueball) | Skyrmion (baryon) | NO |
| Witten string | Cosmic string (topological) | None | NO (fermion ON string) |
| Coleman bosonization | Sine-Gordon soliton | N/A (1+1D only) | YES (in 1+1D) |
| Superfluid vortex | Kelvin wave | Roton/vortex ring | NO |
| **SCP** | **Braid** | **Oscillon (if found)** | **Proposed** |

**No established model identifies a 1D extended soliton with an electron
in 3+1 dimensions.** The closest precedent is the Coleman-Mandelstam
bosonization, which works rigorously only in 1+1D. The SCP proposal
would be genuinely novel if it could be made quantitative.

The Skyrme-Faddeev-Niemi connection is the most promising: it already
has both 1D (knot) and 3D (Skyrmion) solitons in the same theory.
However, the 1D objects are identified with glueballs, not electrons.
Reinterpreting them as electrons would require:
1. A mechanism to give them spin-1/2 (not present classically)
2. A mass scale 10^3 times lighter than the Skyrmion
3. A coupling to a U(1) gauge field (EM)

**Verdict**: No existing literature supports the braid-as-electron
identification in 3+1D. The hypothesis is novel but faces the same
obstacles that all soliton-as-fermion proposals face: spin statistics,
mass hierarchy, and gauge coupling.

---

## 6. Overall Assessment

### Strengths of the Hypothesis

1. **Charge quantization**: Winding number W = +/-1 is an exact topological
   invariant, naturally quantized. This is the strongest point.

2. **Charge-dependent force**: The theta-mediated force has the correct
   sign structure (parallel currents attract, antiparallel repel), with
   quantitative measurements showing 27% enhancement / 57% reduction.

3. **Wave-mediated interaction**: The theta force is closer to QFT photon
   exchange than to classical EM, which is actually the correct physics
   for fundamental particles.

4. **Right-hand rule**: The volumetric visualization (V34) confirms the
   theta circulation follows the right-hand rule around the braid axis.
   This emerges from the curl coupling, not from imposed symmetry.

5. **Self-sustaining dynamics**: The braid's "self-reconstruction through
   motion" (V26-V27) is an attractive mechanism for a stable fundamental
   particle that is not a static solution.

### Weaknesses of the Hypothesis

1. **No Coulomb law**: theta_phi ~ r^(-0.25), not r^(-1). The force law
   has not been shown to be 1/r^2 at any tested separation. This is a
   critical failure for EM identification.

2. **No mass hierarchy mechanism**: The braid energy (5350 code units)
   cannot be simultaneously identified with the electron (0.511 MeV) and
   be 1836x lighter than a proton-like composite. No mechanism produces
   the required energy concentration ratio.

3. **No spin-1/2**: The braid is a classical field configuration with no
   intrinsic angular momentum quantum number. Getting fermion statistics
   requires either:
   - Quantization of the soliton (as in Coleman bosonization)
   - Additional structure (Finkelstein-Rubinstein mechanism for Skyrmions)
   - This has not been addressed

4. **No stable bound state demonstrated**: The required orbit time
   (T ~ 16,000) far exceeds simulation capability. Two-braid systems
   scatter rather than orbit.

5. **1D vs 3D topology**: The braid is a 1D tube, but the electron has
   no spatial extent in the Standard Model. The braid's anisotropy (it
   has a preferred axis) is unlike the isotropic electron.

6. **Periodic BC artifact**: The topological protection of winding number
   depends on periodic BC in z. A finite-length braid in open BC could
   unwind. This is a serious concern for the robustness of "charge."

### What Would Strengthen the Hypothesis

1. **Long-range theta test**: Simulate N=512, L=200, T=1000 to measure
   theta decay at r=50-100. If it transitions to 1/r at long range, the
   EM identification gains significant support.

2. **Bound orbit demonstration**: Run two same-winding braids at D=15
   for T=20,000+ to see if a stable orbit forms. This requires GPU
   acceleration (V36 CUDA port).

3. **Spin from quantization**: Compute the semiclassical quantization
   of the braid's rotational modes. If the helical twist gives half-integer
   angular momentum, this would be a major breakthrough.

4. **Compact 3D soliton at a different scale**: Find a stable 3D
   localized configuration in the same theory at a different energy.
   V37's evolutionary search found compact binding (E_pot = -30) but
   it dissolved by T=200. A truly stable compact structure would
   demonstrate scale separation.

5. **Fine structure constant**: Compute alpha_EM = (theta force at D) /
   (gravitational force at D). If this ratio is ~1/137 at some scale,
   the identification becomes compelling.

---

## 7. Connection to the Faddeev-Niemi Framework

The most productive theoretical direction may be to formally connect
the SCP field theory to the Skyrme-Faddeev-Niemi framework. The
parallels are striking:

| SCP | Skyrme-Faddeev-Niemi |
|-----|---------------------|
| 3 scalar fields phi_a | SU(2) chiral field U |
| Triple product V(P) = V(phi_0*phi_1*phi_2) | Skyrme term tr([L_mu, L_nu]^2) |
| Braid (1D helical tube) | Faddeev-Niemi knot (1D closed loop) |
| Oscillon (3D breathing) | Skyrmion (3D static) |
| Winding number (pi_1) | Hopf charge (pi_3(S^2)) |
| Baryon number (pi_3?) | Baryon number (pi_3(SU(2))) |

The SFN theory already has BOTH 1D (knots) and 3D (Skyrmions) solitons
in the same framework, with a rich mathematical structure. The SCP triple
product V(phi_0*phi_1*phi_2) may be relatable to the Skyrme term through
some field redefinition.

Key difference: the SFN knots are CLOSED loops, while the SCP braid is
an OPEN tube (made periodic by BC). Closing the braid into a loop would
create a toroidal soliton, which is a Faddeev-Niemi hopfion.

The suggestion is: pursue the mathematical connection between the SCP
Lagrangian and the SFN Lagrangian. If they are related, the full
machinery of the SFN model (topological classification, quantization,
mass spectrum) becomes available.

---

## References

### Literature

- Skyrme, T.H.R., "A non-linear field theory," Proc. R. Soc. Lond. A 260 (1961) 127
- Faddeev, L.D. & Niemi, A.J., "Knots and particles," Nature 387 (1997) 58
- Cho, Y.M., "New topological structures of Skyrme theory," Eur. Phys. J. C 77 (2017) 655
- Witten, E., "Superconducting strings," Nucl. Phys. B 249 (1985) 557
- Coleman, S., "Quantum sine-Gordon equation as the massive Thirring model," Phys. Rev. D 11 (1975) 2088
- Hietarinta, J. & Salo, P., "Faddeev-Hopf knots: dynamics of linked un-knots," Phys. Lett. B 451 (1999) 60
- Battye, R.A. & Sutcliffe, P.M., "Knots as stable soliton solutions," Phys. Rev. Lett. 81 (1998) 4798
- Kondo, K., "Glueball mass from quantized knot solitons and gauge-invariant gluon mass," Phys. Rev. D 74 (2006) 125003
- Gleiser, M. & Krackow, M., "Resonant configurations in scalar field theories," Phys. Rev. D 100 (2019) 116005

### SCP Project Files

- `/home/d/code/scp/v38/SPECULATION.md` -- the hypothesis
- `/home/d/code/scp/CONCEPT.md` -- master concept document
- `/home/d/code/scp/DISCOVERIES.md` -- complete V21-V37 record
- `/home/d/code/scp/v34/torsion_coupling/theta_characterize/RESULTS.md` -- theta field measurements
- `/home/d/code/scp/v34/torsion_coupling/theta_characterize/data/theta_radial.tsv` -- raw radial data
- `/home/d/code/scp/v34/torsion_coupling/RESULTS_cosserat.md` -- Cosserat system results
- `/home/d/code/scp/v34/torsion_coupling/OBSERVATION_theta_righthand.md` -- right-hand rule
- `/home/d/code/scp/v33/V33_RESULTS.md` -- force law and braid energetics
- `/home/d/code/scp/v23/phonon/RESULTS.md` -- inter-oscillon equilibrium
- `/home/d/code/scp/v37/RESULTS_evolutionary_search.md` -- compact structure search
- `/home/d/code/scp/v35/SUMMARY.md` -- hydrogen spectrum and multi-scale architecture
