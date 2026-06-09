# v60 Lagrangian — Target Specification (Living Document)

**Date**: 2026-05-25 (initial)
**Status**: Gated. This document will be expanded only after at least one viable candidate from the G9 or G1 tracks survives initial filters. Until then it records the *requirements* that the eventual dynamical Lagrangian must satisfy.

**Gate condition**: A live induced-metric recast (G9) *or* a live two-piece Y + breaking chain (G1) that passes the polarization / rank-tension numerical audits and preserves the v59 precision bridges.

---

## 1. Non-negotiable requirements (derived from v59 CLOSEOUT + trusted skeleton)

The Lagrangian \(\mathcal{L}\) (or first-order action) on the appropriate bundle (after G9 recast) must:

1. **Reproduce the algebraic structures as solutions**:
   - Static (or slowly varying) solutions whose grade projections yield the Brannen kernels on the sector-selected ambients (L for leptons, F for d-quarks, L⊕F for u-quarks).
   - The eigenvalues must be exactly the Brannen form with \(Q = 2/3\) (leptons) and the corresponding structural \(Q_N\) for quarks.
   - The vacuum constraint surface must enforce \(t^2 = 1/2\) (or the sector-specific equilibrium) as a consequence of the potential / EL equations, not an external imposition.

2. **Derive (or exactly match) the structural numbers**:
   - Gauge prefactors 5 and 2 from Killing-form embeddings or equivalent geometric quantities.
   - \(\sin^2\theta_W = 2/9\) from the Brannen phase (or an equivalent dynamical mechanism).
   - The 784 scale for the electroweak vev as the dimension of the operator algebra resolving non-associativity (or its post-G1 two-piece analog).
   - Gravity charge as the second moment of the mass kernel, EP-exact.
   - Radial law \(\square \Omega = f_g \rho_{\rm grav}\) (or its curved-space generalization) in the appropriate limit.

3. **Gravity sector (post-G9)**:
   - The carrier (2-form or projected multivector) plus any soldering/simplicity constraint must yield a symmetric tensor with exactly two transverse-traceless physical degrees of freedom in the linearized spectrum.
   - The long-range force must couple to \(\rho_{\rm grav}\) with the observed (or \(\alpha^{21}\)-class) strength.
   - Equivalence principle must hold (universal charge/mass = 1 across sectors).

4. **Spectrum and rank tension (post-G1)**:
   - The mass bilinear that appears in the Yukawa / mass terms must be the "active" rank-3 piece of a larger (two-piece or chain-broken) object whose total Frobenius norm supplies the 784 bridge scale.
   - The heavy / Goldstone sector must not source long-range gravity (or must be confined / massive).

5. **Gauge content and breaking**:
   - The silent SU(2)_L must be gaugable, with the three Goldstone modes from the Brannen constraint surface eaten (or identified) to give the correct \(W^\pm, Z^0\) masses.
   - The full Spin(7) → Pati-Salam → SM chain must be compatible with the dynamics.

6. **Falsifiability and precision**:
   - The theory must pass all v59 precision tests (lepton masses to \(3\times10^{-14}\) given one scale, EW masses at 0.02–0.04%, bridge at 0.068%, etc.).
   - It must be compatible with LIGO (or provide a principled reason why tensor modes are absent or modified at observable distances).
   - All free parameters beyond the single dimensionful scale \(a_\ell\) (and possibly \(\alpha\)) must be eliminated or shown to be structural.

---

## 2. Current (pre-gate) content

- The tree-level effective form from `../v59/LAGRANGIAN.md` is copied here as the initial target sketch (see `v59_target_lagrangian.md` — to be created on first expansion).
- Every structural identification in that document is now annotated "must be derived from EL" rather than "is postulated."

---

## 3. Evolution rules

- This file is updated whenever a G9 or G1 candidate becomes live enough to constrain the form of \(\mathcal{L}\).
- Once the gate opens, the document becomes the checklist against which the first candidate Lagrangian is measured.
- At v60 closeout, the final version of this spec (plus the actual derived \(\mathcal{L}\)) will be the authoritative reference for the dynamical theory.

---

**This is a requirements document, not a derivation. The derivation lives in the post-gate work under this directory.**