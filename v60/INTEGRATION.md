# v60 Integration — Delta from v59

**Date**: 2026-05-25 (kickoff)
**Parents**: [`../v59/INTEGRATION.md`](../v59/INTEGRATION.md), [`../v59/CLOSEOUT.md`](../v59/CLOSEOUT.md), [`README.md`](README.md)

This is the *delta* document. The full integration picture (v58 multivector force + v59 algebraic skeleton + Furey construction + scale bridges) remains in `v59/INTEGRATION.md`. This file records only what changes or must be re-integrated as a consequence of the G9 + G1 attack and the move to a true dynamical Lagrangian.

---

## 1. What stays the same (the trusted v59 core)

- The single-source algebra \(\mathrm{Cl}(7)_{\rm even} \cong \mathbb{C}\otimes\mathbb{H}\otimes\mathbb{O}\), its grade decomposition, the L⊕F bisection, and the forcing of lepton = L.
- The Brannen kernel on the 3-generation space, the structural origin of \(Q = 2/3\) (both constraint-surface and Lie-algebraic), and the sedenion \(S_3\) origin of the "/3".
- The gauge decomposition \(\mathrm{Spin}(7) = G_2 \times \mathrm{SU}(2)_L \times \mathrm{SU}(2)_R \times \mathrm{U}(1)_{B-L}\), the dual-Coxeter 5, and \(\sin^2\theta_W = 2/9\).
- The 784-dimensional operator algebra on \(\mathrm{End}(L)\) as the resolution of non-associativity (generic for absolutely irreducible adjoint representations).
- Gravity charge = second moment of the mass kernel (EP-exact, Frobenius² structure shared with the bridge at 6/784 for leptons).
- The radial law for a massless kernel with nonzero monopole.
- The exponent 21 in the magnitude (data-forced to \(\dim\mathrm{Spin}(7)\)) and the instanton / top-form reading.
- All Lean theorems that were axiom-clean at v59 closeout (8278-job build, minimal axioms).

These form the *target specification* for the dynamical theory: any proposed \(\mathcal{L}\) or induced-metric recast must reproduce them (exactly or to the documented precision) as derived consequences.

---

## 2. What must be re-examined or extended

### 2.1 Gravity sector (the largest delta)

**v59 picture** (now provisional):
- Long-range force carried by a scalar connection \(\Omega \in \Lambda^2\) sourced by \(\rho_{\rm grav} = {\rm Tr}(M^\dagger M)\).
- Magnitude \(G_e = (21/16)\alpha^{21}\) (or close variant).
- Radial law \(\square \Omega = f_g \rho_{\rm grav}\).

**v60 requirement**:
- The carrier must be recast (or shown to induce) a spacetime tensor mode with exactly two transverse-traceless degrees of freedom (LIGO-compatible) or a principled reason must be given why the observable gravity in this theory is scalar or short-range.
- The simplicity / soldering constraint (if used) must be compatible with the algebraic structure (grade projectors, complex structures \(J \in L\), color \(\mathfrak{su}(3)\), etc.).
- The second-moment charge and the \(\alpha^{21}\) (or revised) magnitude must emerge from the same geometric object that supplies the tensor modes.

**Integration points**:
- The existing polarization decomposition (`g9_polarization_test.py`, `G8G9_Gravity.lean`) becomes the baseline test harness for every induced-metric candidate.
- The OBE radial test (`obe_radial_test.py`) must be re-run on the emergent curved-space equation.
- Magnitude fitting (`gravity_magnitude_test.py`) may need revision once the measure / Jacobian of the emergent metric is known.

### 2.2 Mass bilinear and bridge (rank tension)

**v59 picture**:
- Democratic full-rank \(Y \in {\rm End}(L)\) with \(\|Y\|_F^2 = 784 a_\ell^2\) sets the electroweak scale.
- Physical spectrum is rank-3 with \(\|Y\|_F^2 = 9Qa^2\).

**v60 requirement**:
- A two-piece (or chain-broken) structure that lets the bridge keep its 784 count while the active sector is rank-3 and sources gravity.
- The 0.068% numerical match must either be preserved exactly or derived from the geometry of the split.
- The Goldstone / heavy sector (the ~25 directions) must be shown to be either inert for gravity or to contribute only short-range / gauge forces.

**Integration points**:
- `formalize_bridge.py` and `EwScaleBridge.lean` are the reference implementations of the tension; they will be extended rather than replaced.
- `gravity_charge_test.py` must be updated to accept a candidate split and verify that only the rank-3 piece sources the long-range charge.
- The Higgs vacuum (`XiVacuum.lean`) must be re-interpreted as the vacuum of the *active piece*.

### 2.3 Selection rule and L⊕F dynamics

The Z₂×Z₂ pattern (lepton = L, d-quark = F, u-quark = L⊕F) was an empirical observation in v59. In the dynamical theory it must emerge from the form of the Lagrangian (or from the breaking chain in G1). The additive identity \(D_u = D_e + D_d = 63\) remains a strong structural hint.

### 2.4 The preliminary Lagrangian (v59/LAGRANGIAN.md)

This document is re-interpreted as the **target specification**, not the derived theory. Every structural identification (gauge couplings from Killing indices, \(\sin^2\theta_W = 2/9\) from Brannen phase, Yukawa structure from the Brannen kernel, etc.) must be shown to follow from the EL equations of the new \(\mathcal{L}\) once the G9 carrier and G1 mass bilinear are in place.

---

## 3. New or modified interfaces

- **Geometry ↔ Algebra**: The induced-metric recast must map the internal grades and complex structures of \(\mathrm{Cl}(7)_{\rm even}\) onto spacetime tensors and forms in a way that preserves the already-proved theorems (e.g., lepton = L forced by \(B^2 = -1\), color pinning by \(J_c\)).
- **Geometry ↔ Dynamics**: The OBE equation \(\square \Omega = f_g \rho_{\rm grav}\) becomes (at minimum) the weak-field limit of the EL equation for the 2-form + soldering system.
- **G1 split ↔ Gauge content**: The Goldstone sector from the breaking chain must map onto (or be consistent with) the silent SU(2)_L and the Pati-Salam factors already derived.
- **Full \(\mathcal{L}\) ↔ Algebraic skeleton**: The EL equations must admit the Brannen kernels on the grade-selected sectors as solutions, with the correct eigenvalues, the correct vacuum constraint surface \(t^2 = 1/2\), and the correct phase \(\varphi = 2/9\).

---

## 4. Documentation and version control

- The master cohesive picture remains in `v59/UNIFIED_THEORY.md` until G9 + G1 are resolved; `v60/UNIFIED_THEORY.md` (if created) will be a delta or a full replacement only at closeout.
- All new geometry, representation, and Lagrangian work is documented in the `v60/` subdirectories with the standard `0X_findings.md` + Lean modules pattern.
- At the end of v60 (or at a major checkpoint), a consolidated `v60/CLOSEOUT.md` will re-audit every claim with updated status tags, and the root `CONCEPT.md` / `DISCOVERIES.md` / `FUTURE.md` will be updated from the authoritative version documents.

---

**Bottom line**: v59 gave us a mathematically beautiful and numerically tight *skeleton*. v60 must give it *dynamics and a face that can look at the sky (LIGO) and the spectrum (three light generations) without contradiction*. The integration task is to ensure that every new geometric or breaking construct is compatible with — and ultimately derives — the trusted algebraic results.