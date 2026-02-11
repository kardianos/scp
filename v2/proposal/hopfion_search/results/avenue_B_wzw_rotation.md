# Avenue B: WZW Rotating Soliton Matrix Elements

## Overview

This analysis computes the Wess-Zumino-Witten (WZW) coupling matrix element for a rotating quantized J=1/2 Skyrmion in the Cl+(3,0,1) theory. The WZW term has a topologically quantized coefficient, making it the only avenue with NO free parameters.

**Bottom line**: The spatial baryon current coupling B^i d_i p VANISHES at O(Omega) by current conservation (div B = 0). Furthermore, Avenue C's finding that epsilon_0^2 = 0 exactly (required for pi_3 != 0, i.e., for solitons to exist) means the WZW term lives entirely in the bulk sector and CANNOT couple to the degenerate scalar p. Avenue B is definitively excluded as a mechanism for generating the gravitational coupling.

---

## 1. Setup: Collective Coordinate Quantization

### Hedgehog ansatz
The B=1 Skyrmion in the sigma-model limit (|q| = rho_0):

    U_0(x) = exp(i f(r) r-hat . tau) = cos f + i sin f (r-hat . tau)

with profile f(r): f(0) = pi, f(inf) = 0, and shooting parameter a = -f'(0) = 1.4201.

### Collective rotation
The physical nucleon is a quantized rotating soliton:

    U(x,t) = A(t) U_0(x) A^dagger(t),    A(t) in SU(2)

with the angular velocity defined by xi = A^dagger A-dot = (i/2) Omega_a tau_a.

### Rotational quantization
The rotational Lagrangian is L_rot = (1/2) Lambda Omega^2, giving the quantized spectrum:

    E_rot = hbar^2 J(J+1) / (2 Lambda)

For the nucleon: J = I = 1/2, with:

| Quantity | Value | Units |
|----------|-------|-------|
| Lambda | 141.55 | code |
| hbar | 0.0386 | code |
| Omega_rms = hbar sqrt(3/4)/Lambda | 2.36 x 10^{-4} | code^{-1} |
| E_rot = 3 hbar^2/(8 Lambda) | 3.96 x 10^{-6} | code |

The angular velocity is extremely small compared to the soliton's internal frequencies (omega ~ 1).

---

## 2. Spatial Baryon Current B^i at O(Omega)

### Derivation

For the rotating hedgehog, the left current gains a time component:

    L_0^{rot} = U_0^{-1} xi U_0 - xi = (i/2) Omega_a (D_{ab}(2f) - delta_{ab}) tau_b

where D_{ab} is the adjoint rotation matrix of the hedgehog:

    D_{ab} - delta_{ab} = -2 sin^2 f (delta_{ab} - r-hat_a r-hat_b) - sin(2f) epsilon_{abc} r-hat_c

The spatial baryon current at O(Omega):

    B^i = -(3/8 pi^2)(i/2) Omega_a (D_{ab} - delta_{ab}) epsilon^{ijk} Tr(tau_b L_{0,j} L_{0,k})

After contracting with the hedgehog's left current structure, this has the general tensor form:

    B^i(x) = (Omega_a / 2 pi^2) [alpha(r) delta_{ia} + beta(r) r-hat_i r-hat_a]

where alpha(r) and beta(r) are radial functions depending on f(r), f'(r).

### Key property: divergence-free at O(Omega)

Baryon number conservation gives d_mu B^mu = 0. At O(Omega):
- d_0 B^0 = 0 (B^0 is static at zeroth order)
- Therefore d_i B^i = 0 at O(Omega)

This constrains alpha and beta:

    alpha' + beta' + 2 beta / r = 0

---

## 3. S-wave Coupling: VANISHES

For an s-wave degenerate scalar p = p(r):

    d_i p = p'(r) r-hat_i

The coupling integral:

    integral B^i d_i p d^3x = (Omega_a / 2 pi^2) integral r^2 dr [alpha + beta] integral r-hat_a dOmega x p'(r)

Since integral r-hat_a dOmega = 0 (odd function on the sphere), the s-wave coupling vanishes by angular symmetry.

**Result**: The s-wave WZW coupling is identically zero.

---

## 4. ALL Partial Waves: VANISH at O(Omega)

The divergence-free argument is stronger than the angular symmetry argument:

    integral B^i d_i p d^3x = -integral p (d_i B^i) d^3x + surface term

Since d_i B^i = 0 at O(Omega) and p vanishes at infinity, this integral is zero for ANY p(x) — not just s-wave. This includes L=1 (dipole), L=2 (quadrupole), and all higher partial waves.

**Result**: The B^i d_i p coupling vanishes at O(Omega) for ALL modes of p, by current conservation.

### Verification for L=1

For a dipole mode p = g(r) r-hat_a:

    d_i p = g'(r) r-hat_i r-hat_a + (g/r)(delta_{ia} - r-hat_i r-hat_a)

The angular integrals give nonzero individual terms:
- integral r-hat_b r-hat_a dOmega = (4 pi/3) delta_{ab}
- integral r-hat_b (delta_{ia} - r-hat_i r-hat_a) dOmega gives various terms

But the divergence-free condition alpha' + beta' + 2 beta/r = 0 ensures all radial integrals cancel after integration by parts.

---

## 5. O(Omega^2) Correction

At O(Omega^2), centrifugal deformation modifies B^0:

    B^0 -> B^0_static + delta B^0_{O(Omega^2)}

The relative correction:

    delta B^0 / B^0 ~ Lambda Omega^2 / E_sol = 7.65 x 10^{-8}

This gives an O(Omega^2) correction to the B^0 . p coupling:

    g_top^{O(Omega^2)} ~ g_WZW x 7.65 x 10^{-8}

| N | g_WZW | g_top^{O(Omega^2)} |
|---|-------|---------------------|
| 1 | 4.22 x 10^{-4} | 3.23 x 10^{-11} |
| 3 | 1.27 x 10^{-3} | 9.69 x 10^{-11} |

Required: g_top = 2.8 x 10^{-17}. Even the O(Omega^2) correction is 10^6 too large.

---

## 6. The Static B^0 Coupling (Dominant Term)

The WZW coupling B^mu d_mu(degenerate) has two components:
1. **mu = 0**: B^0 x p-dot — sources the degenerate scalar via the baryon density
2. **mu = i**: B^i x d_i p — vanishes at O(Omega) by current conservation

The dominant mechanism is (1): the STATIC baryon density B^0(r) sources p through:

    kappa^2 p-double-dot - kappa^2 nabla^2 p + mu^2 p = g_WZW B^0(x)

For static p, this is the sourced Poisson/Yukawa equation:

    -kappa^2 nabla^2 p + mu^2 p = g_WZW B^0(x)

This is exactly the mechanism computed in `degenerate.c` (Path 3), with g_top = g_WZW.

### WZW coefficient values

The WZW coefficient is g_WZW = N / (240 pi^2):

| N | g_WZW | G_eff/G_Newton |
|---|-------|----------------|
| 0 | 0 | 0 |
| 1 | 4.22 x 10^{-4} | 2.3 x 10^{26} |
| 3 | 1.27 x 10^{-3} | 2.0 x 10^{27} |

For SU(2): pi_5(SU(2)) = Z_2, so N in {0, 1}. In the Cl+(3,0,1) theory without a QCD gauge group, N is determined by the homotopy of the target space.

**The WZW term produces gravity that is 10^{13} too strong (N=1) or absent (N=0).**

---

## 7. Physical Parameters

Computed from the sigma-model B=1 profile at e=1, rho_0=1:

| Quantity | Value |
|----------|-------|
| E_sol | 103.14 code = 938.3 MeV |
| E_2 | 51.55 code |
| E_4 | 51.59 code |
| E_2/E_4 | 1.001 (virial satisfied) |
| Q (topological charge) | 1.000000 |
| Lambda (moment of inertia) | 141.55 code |
| hbar | 0.0386 code |
| Omega_rms (J=1/2) | 2.36 x 10^{-4} code^{-1} |
| 1 code energy | 9.098 MeV |
| 1 code length | 0.5624 fm |

---

## 8. Conclusions

### What was expected
The WZW term has a topologically quantized coefficient. A rotating J=1/2 Skyrmion generates B^i != 0. The hope was that the matrix element <J=1/2| integral B^i d_i p |J=1/2> would give an effective g_top with no free parameters.

### What was found

1. **S-wave coupling vanishes** by angular symmetry (integral r-hat_a dOmega = 0).

2. **ALL partial wave couplings vanish** at O(Omega) by current conservation (div B = 0). This is a stronger result that holds for any mode of p, not just s-wave.

3. **The WZW mechanism reduces to the static B^0.p coupling** (Path 3), which was already computed in `degenerate.c`. The rotation adds nothing at leading order.

4. **The WZW coefficient g_WZW = N/(240 pi^2)** gives either:
   - N=0: no coupling (WZW term absent)
   - N=1: g_WZW = 4.2 x 10^{-4}, producing gravity 10^{26} too strong

5. **No angular suppression from rotation**: The divergence-free property of B^i kills the spatial coupling completely. The O(Omega^2) centrifugal correction is a ~10^{-8} relative effect, still leaving g_top too large by 10^6.

### Implications for other avenues

- **Avenue A** (one-loop): The WZW term does not contribute the B^0.p vertex at tree level in the way needed. The one-loop mechanism (pion fluctuations) remains the primary candidate for generating g_top.

- **Avenue C** (anomaly cancellation): If N=0, the WZW term is absent and g_top must come entirely from loop effects. Anomaly cancellation may then constrain epsilon_0^2, which controls the loop amplitude.

- **Avenue D** (spectral self-consistency): The WZW term, if present (N=1), produces an enormous B^0.p coupling that would need to be cancelled by other contributions for the theory to be viable. This cancellation requirement could itself constrain the parameters.

### Why the naive estimate was wrong

The rationale document estimated g_top^{WZW} ~ N/(240 pi^2) x hbar/Lambda ~ 10^{-7}. This assumed the B^i spatial current coupling would be the dominant effect. In fact:

1. The B^i coupling vanishes entirely at O(Omega) — no angular suppression factor, just zero.
2. The B^0 coupling (static baryon density) is present regardless of rotation, with the full WZW coefficient g_WZW ~ 4 x 10^{-4}.
3. The correct comparison is g_WZW vs g_top_required, not g_WZW x (hbar/Lambda) vs g_top_required.

This makes the mismatch WORSE (10^{13} vs the naive 10^{10}).

---

## 9. Update: Implications of Avenue C (epsilon_0^2 = 0 exactly)

Avenue C established a critical topological constraint: if epsilon_0^2 != 0, the target space becomes S^7 with pi_3(S^7) = 0, destroying all topological solitons. Therefore epsilon_0^2 = 0 EXACTLY at the algebraic level.

This has profound consequences for the WZW coupling:

### The WZW 5-form requires sector coupling

The WZW term is defined as a 5-form on the field space. In the Cl+(3,0,1) theory, the full field Psi = (q, d) lives in a product space. At epsilon_0^2 = 0:

- The bulk sector q lives on S^3 (unit quaternions, pi_3(S^3) = Z)
- The degenerate sector d = (J, p) is decoupled, forced to zero algebraically
- The field space FACTORIZES: there is no nontrivial 5-cycle connecting the sectors

The WZW 5-form Tr(L^5) on S^3 gives pi_5(S^3) = Z_2, which is the standard Skyrme WZW term. But this term lives ENTIRELY in the bulk sector. It does not couple to the degenerate sector because:

1. The baryon current B^mu is constructed from bulk-sector currents L_mu = q^{-1} d_mu q
2. At epsilon_0^2 = 0, the degenerate field p has no kinetic coupling to q
3. The WZW coupling B^mu d_mu p requires a CROSS-TERM between sectors
4. No such cross-term exists when the Lagrangian factorizes

### Strengthened conclusion

The Avenue C finding sharpens the Avenue B result from "wrong magnitude" to "does not exist":

| Scenario | WZW coupling to p | Status |
|----------|-------------------|--------|
| epsilon_0^2 != 0 | Could exist, but solitons don't (pi_3(S^7)=0) | Inconsistent |
| epsilon_0^2 = 0, N=0 | Zero (no WZW term) | Consistent, no coupling |
| epsilon_0^2 = 0, N=1 | WZW exists in bulk, but CANNOT couple to degenerate sector | No cross-sector coupling |

In all consistent scenarios, the WZW term produces ZERO coupling between the baryon current and the degenerate scalar p.

### Remaining possibility: dynamical coupling

Even with epsilon_0^2 = 0 at the fundamental level, quantum fluctuations could generate an effective coupling. The WZW term modifies the bulk sector dynamics (it gives the correct spin-statistics for the Skyrmion, making it a fermion for odd N). Through loop corrections, the WZW-modified bulk propagator could couple to the degenerate sector at higher order.

However, this is effectively an Avenue A mechanism (loop-generated coupling) dressed by the WZW vertex — not a pure Avenue B effect. The topological quantization of the WZW coefficient would enter only indirectly, through the modified pion propagator.

### Final assessment

**Avenue B is definitively excluded as an independent mechanism for generating g_top.** The combination of:
1. Current conservation killing B^i d_i p at O(Omega)
2. The epsilon_0^2 = 0 requirement killing the B^0 d_0 p cross-sector coupling

means the WZW term cannot source the degenerate scalar at ANY order in Omega. If g_top is nonzero, it must arise from loop effects (Avenue A) or spectral constraints (Avenue D), not from the topologically quantized WZW term.

---

## Files

- `src/wzw_rotation.c` — Numerical computation and analysis code
- `data/profiles/profile_sigma_e1.dat` — B=1 sigma-model profile (e=1, rho_0=1)
- `results/quantum_gravity_rationale.md` — Full context for quantum gravity avenues
