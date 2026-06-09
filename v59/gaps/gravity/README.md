# Gravity gaps — G8 (magnitude mechanism) and G9 (scalar vs tensor)

**Cluster**: gravity · **Gaps**: G8, G9 · **Date**: 2026-05-25
**Folder scope**: all work here is in `/home/d/code/scp/v59/gaps/gravity/`.

Status tags follow `UNIFIED_THEORY.md`:
**[thm]** machine-checkable identity (certain mathematics);
**[emp≈X]** empirical match of a structural number to data at precision X;
**[conj]** ansatz with a structural integer, not derived from a Lagrangian.

---

## 0. Where gravity sits in the v59 OBE

In the unified field equation (`synthesis/NEW_OBE_FORMULATION.md` §6), a test
excitation feels a bivector connection `Ω(x) ∈ Λ²(𝒜)` sourced by three currents.
The gravity piece, after the live results established this session, is:

- **Charge** [thm-of-structure / emp]: gravity is sourced by the **second moment of
  the mass kernel**, `ρ_grav = Tr(M†M) = Σ_k m_k = 9Q a²` (Frobenius² of the 3×3
  Brannen kernel over the generation space). This is **equivalence-principle exact**
  (charge/mass = 1 for all sectors) and *the same* Frobenius² structure as the EW
  bridge `v = dim(L)² a²`. The earlier candidate (Koide-constraint deviation
  `ρ_N ~ a²(|ξ|²−½)`) is **dead**: zero for leptons, scale-free, cannot be mass.
  (`synthesis/gravity_charge_test.py`, `FINDINGS_obe_bridge.md` §4.)
- **Radial law** [thm]: for a **massless** kernel, `□Ω_grav = f_g ρ_grav` gives a
  `1/r` potential (slope −1.000) and a `1/r²` force (slope −2.000), because
  `ρ_grav` carries a nonzero monopole (= total mass). (`synthesis/obe_radial_test.py`.)
- **Magnitude** [conj]: `G_N = f_g²/4π`, target `α_G(e) = (m_e/M_Pl)² = 1.752×10⁻⁴⁵`.

The two remaining gravity gaps are **G8** (what fixes the magnitude `f_g`) and **G9**
(the propagating mode is a scalar, but LIGO needs a spin-2 tensor).

---

## 1. G8 — the gravity magnitude mechanism

### The conjecture
$$
\boxed{\;G_e=\frac{21}{16}\,\alpha^{21},\qquad \frac{m_e}{M_{\rm Pl}}=\sqrt{\tfrac{21}{16}}\;\alpha^{21/2}\;}
\qquad\textbf{[conj, emp}\approx0.25\%\textbf{]}
$$
with `α = α(0)` the IR value (gravity is long-range; `α(M_Z)` misses by 4.2×).

### What is robust vs weak (verified in `g8_exponent_test.py`)

| element | status | evidence |
|---|---|---|
| **exponent = 21** | **strong signal** | the exponent forcing prefactor 21/16 is `n = 21.0005` = dim Spin(7) to **0.0024%**. A falsifier scan over algebra dims (14, 16, 21, 28) shows **only** n=21 lands the prefactor in the O(1) window; n=14/16/28 are off by `10⁻¹⁵`–`10¹⁵`. So the *data demands* exponent exactly 21. |
| **prefactor 21/16** | **weak** | the prefactor that n=21 forces is `1.3092`; nearest low-denominator rationals: **17/13 = 1.3077 (0.12%)** beats 21/16 = 1.3125 (0.25%); 4/3 = 1.333 (1.8%). So "21/16 = dim Spin(7)/dim Cl(3,1)" is **not uniquely selected**. |
| **21 itself** | **over-determined** [thm] | `21 = dim Spin(7) = 28−7 = 14+7 = 35−14` — four v59 readings, so the integer alone is not a unique fingerprint. |
| **mechanism** | **open** | no Lagrangian derives the magnitude. |

### The additive-vs-multiplicative hypothesis
Gauge couplings are **additive** over generator subsets: `g_W² = 5√α`, with
`5 = 21 − 16 = h∨(Spin(7))` the dual Coxeter number [thm]. Gravity is conjectured
**multiplicative** over **all 21** Spin(7) generators:
`g_grav² = (√α)^(2·21) = α²¹`. The picture: gravity is the *unique* coupling that
needs **all** generators acting *coherently* (a product), so it is `α^(21/2)`-
suppressed, whereas gauge forces use additive subsets and stay `O(√α)`. This would
**explain why gravity alone is so weak**, and the V6 historical `~10⁴⁰` overshoot
is *exactly* the missing `α²¹ ≈ 10⁻⁴⁵` suppression. **[conj — no Lagrangian
derives the product.]**

---

## 2. G9 — scalar vs tensor (the decisive LIGO test)

### The gap
The propagating gravity mode in the OBE is a **scalar**: `□Ω_grav = f_g ρ_grav`
with `ρ_grav = Tr(M†M)` a **Lorentz scalar** (a generation-space trace). A massless
scalar carries **one** polarization, helicity `h = 0`. LIGO detects the
**transverse-traceless `h = ±2` strain** of a spin-2 graviton. A scalar (or vector,
`h = ±1`) force cannot produce the observed `h = ±2` waveform.

### Why it is fatal
Helicity is a **spacetime little-group** label. The hope that "gravity lives in the
28-dim `L = so(8)/Λ²` bivector" supplies spin-2 fails because the `Λ²` index is
**internal** (so(8)), and internal indices are **inert** under spatial rotations —
they add no spacetime helicity. The polarization decomposition (`g9_polarization_test.py`,
helicity along z) confirms:

| carrier | spacetime helicities | spin | LIGO? |
|---|---|---|---|
| scalar `S(x)` (the v59 ρ_grav mode) | `{0}` | 0 | no |
| 4-vector `A_μ` | `{±1}` (physical) | 1 | no |
| antisymmetric 2-form `F_μν` (EM-like bivector connection) | `{±1}` | 1 | no |
| **symmetric traceless transverse `h_μν`** | **`{±2}`** | **2** | **yes** |

Only a **symmetric spacetime rank-2 tensor** carries `h = ±2`. A *spacetime*
bivector connection `F_μν` is **antisymmetric** ⇒ spin-1, *also* fails. So even
promoting `Ω` from internal to a spacetime 2-form does not help.

The one place the bivector structure can reach spin-2 is the **quadratic bilinear**
`T_μν = Ω_μ^a Ω_ν^a` (the gauge-field stress tensor, symmetric in μν, TT part = `±2`).
But that is a **source** (the energy-momentum that sources gravity), `O(g²)=O(α)`-
suppressed and short-range — **not** a fundamental massless graviton.

**Verdict:** as currently formulated, **v59 gravity is scalar and fails LIGO.** To
pass, it needs a genuine symmetric-tensor massless field `h_μν` (a graviton in the
*spacetime* sector), which the present internal-bivector OBE does **not** contain.
G9 is the more decisive of the two gaps: magnitude is moot if the polarization is
wrong.

---

## 3. Files in this folder

| file | purpose |
|---|---|
| `README.md` | this background (states both gaps precisely, tagged) |
| `ALTERNATIVES.md` | solution space for G8 (product/det/top-form/instanton) and G9 (routes to a spin-2 mode), each with test + falsifier |
| `g8_exponent_test.py` | numerically verifies α²¹, pins exponent=21, tests the product/det/top-form counting, falsifier scan, instanton/transmutation alternatives |
| `g9_polarization_test.py` | helicity decomposition of scalar/vector/2-form/sym-tensor; shows internal so(8) index is inert; only sym-TT tensor carries h=±2 |
| `g9_spin2_route_test.py` | tests the two spin-2 routes: (A) DOF count — a massless 4D 2-form is dual to a scalar (1 DOF, h=0), provably not a graviton; (B) `Sym²(so(8) adjoint)` *does* contain a spin-2 rep, but it stays internal (no soldering) and composite (quadratic, short-range) |
| `G8G9_Gravity.lean` | Lean formalization (`import Mathlib`, self-contained): 21=dim Spin(7), 5=21−16, 21/16 prefactor, product-exponent identity `(α^½)^{2·21}=α^21`, det realization, helicity multiset classification. One flagged `sorry` (empirical magnitude match). **Written, not built this run.** |
| `FINDINGS.md` | verdicts on G8 and G9, most-promising avenue |
