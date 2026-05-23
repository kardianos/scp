# v58 + v59 + Furey Integration — One Picture

**Date**: 2026-05-22
**Parents**: [`SESSION_2026-05-22.md`](SESSION_2026-05-22.md),
[`synthesis/FINDINGS_scale_bridge.md`](synthesis/FINDINGS_scale_bridge.md),
[`synthesis/MOTIVATIONS_AND_CONSEQUENCES.md`](synthesis/MOTIVATIONS_AND_CONSEQUENCES.md),
[`furey_construction/ALL_FINDINGS.md`](furey_construction/ALL_FINDINGS.md)

**Purpose**: Now that Step 1 (fermion content) is complete, integrate ALL
prior work — v58 multivector force law, v59 Brannen-Koide kernel, scale
bridge, Furey color algebra, prior Cosserat simulations — into a single
coherent picture. This document is the **map** of how everything fits.

---

## 1. The algebraic foundation

```
                  ┌──── ℂ⊗ℍ⊗𝕆 (64 cplx-dim, full Furey color algebra)
                  │                  ≅ Cl(7)_even
                  │     ┌──────────┼──────────┐
                  │     │          │          │
                  │   Λ⁰=1      Λ²=21       Λ⁴=35       Λ⁶=7
                  │  (idem.)  (Spin(7))  (d-quark)    (S⁷)
                  │
                  │   ── Sector-specific Yukawa contractions ──
                  │
                  │   Lepton (D_L=28):  Λ²⊕Λ⁶   "L content"
                  │   d-quark (D_d=35): Λ⁴       "F content"
                  │   u-quark (D_u=63): Λ²⊕Λ⁴⊕Λ⁶ = D_L+D_d
                  │
                  │   ── 28+35=63: u-quark = L⊕F ──
                  │
                  │   Plus: Z₃ ⊂ S₃ triality of Spin(8) → 3 generations
                  │         (28 = dim Spin(8))
                  │
                  ▼
              ψ_X(x), Φ(x) — physical fields on spacetime
```

**Key identities (mostly axiom-free in Lean):**
- $\dim \mathrm{Cl}(7)_\text{even} = 1+21+35+7 = 64$ — the Furey color algebra
- $D_L = 28 = \dim \mathrm{Spin}(8) = 2\cdot\dim G_2$ — lepton ambient = generation triality
- $D_L + D_d = D_u$ — additive identity (28+35=63)
- $(1-Q_N)\cdot D_N = 28/3$ — universal Koide deviation across sectors

---

## 2. The dynamical fields

| Field | Type | Dimension | Role |
|---|---|---|---|
| $\psi_e^{(i)}(x)$ | Dirac spinor in Cl(3,1) | 4 × 3 gen | Charged leptons + neutrinos |
| $\psi_u^{(i)}(x)$ | Dirac spinor × color | 4 × 3 × 3 | Up-type quarks |
| $\psi_d^{(i)}(x)$ | Dirac spinor × color | 4 × 3 × 3 | Down-type quarks |
| $\Phi(x)$ | Cl(7)_even multivector | 64 cplx | Higgs + sector-specific scalars |
| $\xi(x) = \Pi_L[\Phi]$ | ℍ slice | 4 real | **Lepton Higgs** (current v59 study) |
| $A_\mu^a$ | SU(2)_L gauge | 3 | Weak bosons (silent direction) |
| $B_\mu$ | U(1)_Y gauge | 1 | Hypercharge (origin: TBD) |
| $G_\mu^a$ | SU(3)_c gauge | 8 | Gluons (Furey Cl(6) Witt decomp) |
| $g_{\mu\nu}$ | metric | (4×4 sym) | Spacetime |

**Fermion content** (per generation, from `07_full_generation.py`):
- L-doublet $(\nu_L, e_L)$: $Y = -1$
- L-doublet $(u_L, d_L) \times 3$ colors: $Y = +1/3$
- R-singlets: $\nu_R(Y=0),\ e_R(Y=-2),\ u_R^c \times 3(Y=4/3),\ d_R^c \times 3(Y=-2/3)$
- All anomalies cancel (verified in `07_full_generation.py`).

---

## 3. The Lagrangian (sketch)

The v59 Lagrangian, in schematic form, is now writeable as:

$$
\mathcal{L}_{\rm v59} = \mathcal{L}_\Phi + \mathcal{L}_\text{gauge} + \mathcal{L}_\psi + \mathcal{L}_\text{Yukawa} + \mathcal{L}_\text{grav}
$$

with:

```
L_Φ        = (1/2)⟨(D_μΦ)(D^μΦ̃)⟩₀  − V(|Φ|²)         (Higgs sector kinetic + potential)
L_gauge    = -(1/4)[F^a_μν F^{aμν}]_SU(2) − (1/4)B²μν − (1/4)[G^a G_a]                  (all 3 gauge factors)
L_ψ        = i ψ̄_X γ^μ D_μ ψ_X                       (Dirac kinetic + gauge coupling)
L_Yukawa   = Σ_X y_X ψ̄_X^{(i)} M_X(ξ_X)_{ij} ψ_X^{(j)}      where M_X = a_X(I + ξ_X S + ξ̄_X S²)
L_grav     = R/(16π G_N) + ρ_M · g_μν              with ρ_M = (1/2)(|Φ|² - v_ref²)
```

Brannen kernel $M_X = a_X(I + \xi_X S + \bar\xi_X S^2)$ replaces the SM 3×3 Yukawa
matrix with TWO parameters per sector ($a_X$, $\xi_X$), constrained by $|\xi_X|^2 = 1 - \dim G_2/D_X$.

---

## 4. What this looks like at each "layer" of the SCP project

The SCP project's history is a layered understanding. Here's how v59 sits in it:

| Layer | What it provided | What v59 does with it |
|---|---|---|
| v6–v28 | Skyrme solitons, baryon binding, hopfion | Numerical context; not directly used in v59 |
| v34 | 6-field Cosserat reference (3 phi + 3 theta) | Provides the 6 bivector degrees of freedom of Cl(3,1) (bivector grade ≅ 6-dim) |
| v50–v55 | GPU SFA infrastructure, large-scale lattice | Available for **Stage 2 simulation** of v59 Lagrangian |
| v58 | Multivector force law, ρ_M source, chiral current | Provides the **dynamical bosonic field structure**; v59 ξ ∈ ℍ is the lepton-sector slice of v58's M ∈ Cl(3,1) |
| v59 (this) | Brannen-Koide kernel, Koide Q = 2/3, scale bridge, EW sector closure | **Mass and gauge sector specification** |
| Furey | Color algebra ℂ⊗𝕆 gives SM fermion content | **Fermion content + charges** (Variants B + new Variant H in `07_full_generation.py`) |

**The 6-field Cosserat connection (potentially deep)**:
The v34/v50–v55 simulations used 6 fields (3 phi + 3 theta) with coupling
$\eta \cdot \mathrm{curl}(\theta)$ in the phi-EOM and vice versa. These 6 fields
naturally correspond to the **6-dim bivector grade of Cl(3,1)**:
- 3 phi = 3 spatial bivectors $\{e_{23}, e_{31}, e_{12}\}$ (the quaternion units)
- 3 theta = 3 boost bivectors $\{e_{01}, e_{02}, e_{03}\}$

The Brannen $\xi \in \mathbb{H}$ uses scalar + 3 spatial bivectors = 4 components.
The Cosserat 6-field uses ALL 6 bivectors. So the v59 Brannen sector is a
4-dim slice of the 6-dim Cosserat sector, which is a 6-dim slice of the 16-dim
full Cl(3,1) multivector field of v58. **This is the integration**:

$$
\text{v59}\ \xi \in \mathbb{H}\ \subset\ \text{v34/v55 Cosserat}\ (\phi\oplus\theta\in\Lambda^2(\mathrm{Cl}(3,1)))\ \subset\ \text{v58 }M\in\mathrm{Cl}(3,1)
$$

---

## 5. The full Cl(3,1) multivector content (v58 view)

Cl(3,1) has dim 16, decomposing as:

| Grade | Dim | Components | Physical role |
|---|---|---|---|
| Scalar | 1 | Identity | Higgs radial mode / gravity source ρ_M |
| Vector | 4 | $\gamma^\mu$ | 4-momentum / current direction |
| Bivector | 6 | $\gamma^{\mu\nu}$ | EM tensor $F_{\mu\nu}$ + 6-field Cosserat |
| Pseudovector | 4 | $\gamma^5\gamma^\mu$ | Chiral / axial current |
| Pseudoscalar | 1 | $\gamma^5$ | CP-odd Higgs (or topological term) |

So when v58 writes $M \in \mathrm{Cl}(3,1)$, it encodes:
- Higgs (scalar grade)
- 4-current (vector grade)
- EM/Cosserat (bivector grade)
- Chiral current (pseudovector grade)
- Topological (pseudoscalar)

v59 extracts the *quaternion slice* (scalar + 3 spatial bivectors) and uses it for the
Brannen-Koide mass kernel. The OTHER grades of v58 are gauge fields,
matter currents, and gravity — currently underutilized in v59.

**Open consequence**: the Lagrangian should treat all 16 Cl(3,1) grades
plus the Furey ℂ⊗𝕆 algebra structure simultaneously. The current
Brannen-only formulation captures *only the Higgs sector*.

---

## 6. The simulation ladder — now grounded

With fermion content (Step 1) done, the simulation ladder I outlined
before is now better-defined:

| Stage | Field content | Purpose | Status |
|---|---|---|---|
| **0** | None (analytical) | Verify v59 conjectures against PDG | **COMPLETE** (this session) |
| **1** | Lagrangian only (Maxima symbolic) | Derive v59 conjectures from a candidate L | TODO (next) |
| **2** | $\xi(x) \in \mathbb{H}$ on lattice (Higgs alone) | Verify Brannen vacuum is ground state | Feasible now — `scp_sim.cu` ready |
| **3** | $\xi(x)$ + $A_\mu$ (gauged Higgs) | Verify $m_W, m_Z$ from EW breaking | Stage 2 + SU(2)_L gauging |
| **4** | $+\psi_X$ (Dirac fermions with Brannen-Yukawa) | Verify lepton Koide from spectrum | Needs Step 1 (now done) |
| **5** | Full Furey $\Phi \in \mathrm{Cl}(7)_\text{even}$ + all sectors | Verify cross-sector predictions | Requires Step 3 (CKM check) |
| **6** | Cosmological (curved background) | Test ρ_M as gravity source | Long-term |

---

## 7. The "central spine" of the v59 framework

After this session, the v59 framework can be summarized as:

**Inputs** (truly empirical):
- $a_\ell \approx 17.7\ \sqrt{\rm MeV}$ — Brannen lepton scale (from $\sqrt{m_e}+\sqrt{m_\mu}+\sqrt{m_\tau}$)

**Structural inputs** (geometric, not empirical):
- $\mathbb{C}\otimes\mathbb{H}\otimes\mathbb{O}$ as the central algebra
- Furey color/charge structure from Cl(6) Witt decomposition
- Three generations from $\mathbb{Z}_3 \subset S_3$ triality of Spin(8)
- Brannen kernel $M_X = a_X(I + \xi_X S + \bar\xi_X S^2)$ with $|\xi_X|^2 = 1 - 14/D_X$
- Higgs sector = radial mode of $\xi$ on constraint surface

**Outputs** (with empirical-match precision in parentheses):
- Lepton Koide $Q = 2/3$ (6×10⁻⁶)
- Brannen phase $\varphi_e = 2/9 = \sin^2\theta_W$ (7×10⁻⁶ / 0.4 %)
- $\cos^2\theta_W = 7/9 = t^2_{u\text{-quark}}$ (0.11 %)
- $\alpha(0) \approx 1/137$ from $-\ln\alpha+2\alpha = \pi^2/2$ (0.004 %)
- $\alpha(M_Z) = 25/(324\pi^2)$ from EW-tree consistency (0.03 %)
- $g_W^2 = 5\sqrt{\alpha}$ from Killing-form embedding (0.28 % at α(0))
- $v_\text{Higgs} = 28^2 \cdot a_\ell^2$ (0.07 %)
- $m_W, m_Z, m_H$ at <0.15 %
- $G_e = (21/16)\cdot\alpha^{21}$ for gravity (0.25 %)
- Quark Koide $Q_d = 11/15, Q_u = 23/27$ (~0.3 %)
- $m_\text{top} \approx 172.6$ GeV with $a_u^2 = 72\cdot a_\ell^2$ (0.02 %)
- $y_\text{top} \approx 1$ (the "unit Yukawa" identity, 0.09 %)
- Full SM fermion content with anomaly cancellation
- $m_H \cdot m_Z \cdot \sqrt{n_\text{gen}} = m_W \cdot v$ (0.07 %)
- $(1-Q_N) \cdot D_N = 28/3$ universal (exact algebraic)

**Predictions for new physics** (not yet tested):
- d-quark "Higgs analog" near 800 GeV (from $D_d^2 \cdot a_d^2$)
- u-quark "Higgs analog" near 90 TeV (from $D_u^2 \cdot a_u^2$)
- $\delta\alpha/\alpha \sim 10^{-10}$ in Earth's gravitational potential
- Possibly: cosmological constant from $\langle\rho_M\rangle$ of quark sectors

---

## 8. Open Issues, Ranked

These are what remain. Each is one or more sessions of focused work:

| # | Issue | Frontier | Effort |
|---|---|---|---|
| **A** | Write the candidate Lagrangian (Step 6 of fermion-path) | Cross-cutting | 2–3 sessions |
| **B** | Selection rule (Z₂×Z₂: why each N → which Cl-grade) | Frontier 2 | 3–5 sessions |
| **C** | U(1)_Y geometric origin | Frontier 3 | 1–2 sessions |
| **D** | CKM/PMNS from Brannen eigenvectors | Open | 1 session (computable now) |
| **E** | Lagrangian mechanism for $g_W^2=5\sqrt\alpha$ | Frontier 3 | 2–4 sessions |
| **F** | Lagrangian mechanism for $v=28^2 a_\ell^2$ | Cross-cutting | 1–2 sessions |
| **G** | v58 ρ_M as universal (not sector-dep) gravity | v58⊕v59 | 2–3 sessions |
| **H** | Quark Brannen phases φ_u, φ_d (currently empirical) | Open | Unclear |
| **I** | Cosmological constant from $\langle\rho_M\rangle$ | Speculative | Long |
| **J** | Lattice simulation of Lagrangian (Stage 2 onward) | Numerical | 2-4 sessions, GPU |

---

## 9. The decision tree for next sessions

Given the user's preference for "primary path", the order is:

```
You are here ──→ Step 1 fermion content (DONE)
                     │
                     ▼
              Step 2: Decision — what's the dynamical Higgs field?
                  Option A: ξ ∈ ℍ only (current v59)
                  Option B: Multi-vacuum potential V(|ξ|²)
                  Option C: Effective potential (matter-dependent)
                  Option D: ξ ∈ ℍ is lepton; quarks composite
                  Option E: Φ ∈ Cl(7)_even with sector projections (RECOMMENDED)
                     │
                     ▼
              Step 3: CKM/PMNS from Brannen eigenvectors (testable now)
                     │
                     ▼
              Step 4: Selection rule mechanism (Z₂×Z₂)
                  - try Spin(8) → Spin(7) → G_2 breaking
                  - try Furey Fock-space derivation
                  - try octonion automorphism analysis
                     │
                     ▼
              Step 5: U(1)_Y origin
                     │
                     ▼
              Step 6: Write the Lagrangian
                     │
                     ▼
              Stage 1: Symbolic (Maxima) — derive v59 conjectures from L
                     │
                     ▼
              Stage 2: Lattice simulation of Higgs alone
                     │
                     ▼
              Stage 3-5: Add gauge, fermions, full v59 SM
```

---

## 10. The single insight worth committing to

If one structural reading deserves emphasis, it's this:

> **The v59 Brannen ξ ∈ ℍ is the lepton-sector projection of a richer
> Higgs multivector Φ ∈ Cl(7)_even = ℂ⊗𝕆. Each fermion sector sees its
> own Cl-grade projection, with the L⊕F bisection determining which
> sectors couple to which grades. The same ℂ⊗𝕆 algebra carries:**
>
> - **Fermion charges** (Cl(6) Witt decomposition, Variant B)
> - **Yukawa coupling structure** (sector-specific Cl-grade projections)
> - **Mass scales** (Brannen-Koide kernel on 3 generations)
> - **Three generations** (Z₃ ⊂ S₃ triality of Spin(8) ⊂ ℂ⊗𝕆)
> - **EW boson mixing** (silent SU(2)/U(1) acting on ξ)
> - **Gravity coupling** (constraint deviation, v58 multivector ρ_M)
>
> **This is the unification claim: ONE algebra (ℂ⊗𝕆) sources EVERYTHING.**

The Lagrangian's task is to encode this unification dynamically — to make
the ALGEBRAIC structure into a FIELD-THEORETIC one. The Lagrangian
sketch in Section 3 is the first step; the remaining details (selection
rule, U(1)_Y origin, etc.) need to be filled in before the Lagrangian
is fully specified.

---

## Cross-references

- **Master session**: [`SESSION_2026-05-22.md`](SESSION_2026-05-22.md)
- **Scale bridge findings**: [`synthesis/FINDINGS_scale_bridge.md`](synthesis/FINDINGS_scale_bridge.md)
- **Motivations & consequences**: [`synthesis/MOTIVATIONS_AND_CONSEQUENCES.md`](synthesis/MOTIVATIONS_AND_CONSEQUENCES.md)
- **Furey ALL_FINDINGS**: [`furey_construction/ALL_FINDINGS.md`](furey_construction/ALL_FINDINGS.md)
- **Code (this session)**:
  - `furey_construction/07_full_generation.py` — Step 1
  - `furey_construction/08_brannen_yukawa.py` — Yukawa structure
  - `synthesis/scale_bridge.py` — v_Higgs = 28²·a_l²
  - `synthesis/alpha_consistency.py` — α(M_Z) = 25/(324π²)
  - `synthesis/v59_predictions_all.py` — consolidated table
  - `synthesis/higgs_quartic.py` — Higgs quartic identity
- **Lean** (Mathlib-checked):
  - `furey_construction/lean/ScaleBridge.lean` — new theorems (4 axiom-free)
  - `furey_construction/lean/Predictions.lean` — prior session
  - `furey_construction/lean/AxiomCheck.lean` — audit
