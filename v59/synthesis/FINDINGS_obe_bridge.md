# OBE bridge variants — test results (radial law, Option 7, Option 6)

**Date**: 2026-05-25
**Scripts**: `obe_radial_test.py`, `obe_geometric_6_7.py` (both self-contained, numpy).
**Context**: making the `v_Higgs = 28²·a_ℓ²` bridge workable in the v59 OBE
(`NEW_OBE_FORMULATION.md`), keeping octonion non-associativity and the L⊕F /
Koide constraints in mind. Pursues the **Yukawa-free** readings (Options 6, 7)
and tests whether the gravity term can carry a **1/r²** (long-range) law.

---

## 1. Can the OBE gravity term give 1/r²?  (`obe_radial_test.py`)

OBE gravity piece: `Ω_grav(x) = f_g ∫ K(x,x') ∇'ρ_N(x') d³x'`, `ρ_N = |ξ|² − ½`.

**Integration-by-parts identity (confirmed numerically):** for a *massless* kernel
`K` (static Green's fn of `−∇²`), the gradient comes out of the convolution,
`Ω = f_g ∇(K * ρ_N) = f_g ∇Φ`, with `−∇²Φ = ρ_N` (residual `|−∇²Φ − ρ_N| < 6×10⁻⁵`).
So the "∇ρ_N" source is *not* a problem — the connection is just the **force of a
potential sourced by ρ_N itself**.

**The radial law is governed entirely by the monopole of ρ_N and the mass of K:**

| case | `∫ρ_N d³x` | massless K result | verdict |
|---|---|---|---|
| Profile M (net deviation, Gaussian) | +15.75 | Φ slope **−1.000**, force slope **−2.000**, `\|Φ(10)\|=0.125` | **1/r² ✓** |
| Profile N (breathing, `= ∇²blob`) | ~0 (7×10⁻¹¹) | `\|Φ(10)\|=10⁻¹³`, ratio 8×10⁻¹³ vs M | short range |
| massive K (m=0.5, 1.0) on Profile M | +15.75 | `\|Φ(10)\| = 9.6×10⁻⁴, 9.4×10⁻⁶` | Yukawa screened |
| v59 `|ξ|²=½ − D·e^{−r²/2σ²}` (dip D=0.2, 0.5) | −3.15, −7.87 | force slope **−2.000** | **1/r² ✓** |

**Verdict — 1/r² works, with two conditions:** (i) `K` massless, (ii) `∫ρ_N ≠ 0`.
A localized `|ξ|²` deviation generically *has* a nonzero monopole (a net
Koide-constraint excess/deficit), so 1/r² is **reachable** — this is the same
mechanism that gave V6 its one long-range result (`□δρ = −½|ω|²`), now with the
structural v59 source `ρ_N`.

**Two caveats that are physical, not radial (unchanged by this test):**
1. A fully relaxed lepton sits *on* the constraint (`|ξ|²=½ ⇒ ρ_N≡0 ⇒` no source).
   So `ρ_N` alone cannot be the mass/gravity charge — the known `ρ_M ≠ mass`
   problem (`FINDINGS_synthesis.md`). The bridge needs the charge that a *massive*
   particle actually carries, which is not simply the constraint deviation.
2. The propagating mode is a **scalar** (the density `ρ_N`), so this is scalar
   long-range force, not tensor — the LIGO `h=±2` problem persists.

So: **the radial obstruction is gone; the source-identification and tensor
obstructions remain.** Changing the OBE from `∇ρ_N` (ambiguous) to an explicit
massless `□Ω = ρ_N` makes the 1/r² intent manifest.

---

## 2. Option 7 — Kaluza-Klein separation  (`obe_geometric_6_7.py`)

A massless field on `ℝ³ × M_int` (compact internal manifold, the 28-dim L-coset):

| internal dim `d` | long-range slope (`r≫R`) | short-range slope (`r≪R`) |
|---|---|---|
| 1 | **−1.000** (4D 1/r zero mode) | −1.875 → −2 (full (3+1)D) |
| 2 | **−1.000** | −2.646 → −3 (full (3+2)D) |

**Zero-mode coupling scales as `1/Vol_int`** (measured `4πrG → 1/Vol`): `1.000,
0.500, 0.250` for `Vol = 1, 2, 4` — exact.

**Conclusion:** the radial law (1/r → 1/r² force) is the **4D massless zero mode**
and is *independent* of the internal manifold; the internal volume only sets the
**coupling** (the "spread over the 28 directions" dilution, `1/Vol` per leg). So
**the 1/r² law and the 28-count bridge live in separate factors of `K(x,x') =
K_spatial · K_internal`** and do not conflict. The short-range steepening to the
full `(3+d)`-D law is the massive KK tower (Yukawa corrections) — exactly the
short-range/long-range split we want.

---

## 3. Option 6 — rectangular generation × ambient overlap (Yukawa-free)

Leptons = 3 modes (rows) overlapping `D` internal directions (columns); masses =
**singular values** of a `3×D` overlap matrix `W` (no `ψ̄Φψ` coupling):

- `Q = 0.66666` (=2/3), `φ = 0.2223` (=2/9 to `|Δ|<10⁻⁴`) — Brannen shape intact.
- **Bridge in overlap form:** `√v = (D/N_gen)·Σ√m`. Holds at **0.034%** for the
  lepton ambient **D=28 only**:

  | sector | `D` | `√v/a` | law wants | verdict |
  |---|---|---|---|---|
  | lepton | 28 | 28.01 | 28 | **match** |
  | d-quark | 35 | 19.47 | 35 | fails (1.8×) |
  | u-quark | 63 | 3.29 | 63 | fails (19×) |

- **SVD realization:** a democratic `3×28` `W` (3 orthonormal Z₃-phase column
  patterns) reproduces the singular values `= √m_{e,μ,τ}` to machine precision;
  `‖W‖_* = Σ√m`, and `√v = (28/3)‖W‖_* = 496.04` vs `√v_obs = 496.21` (0.034%).

**Honest limit:** the SVD recovers whatever singular values are inserted, so
Option 6 **relocates** the Brannen shape `(t, φ)` into the generation rows — it
does not derive them. What it *does* make structural and lepton-specific is the
ambient leg count `D` and the bridge `√v = (D/3)Σ√m`, valid only for `D=28`.

---

## 4. Gravity charge — which scalar map? (`gravity_charge_test.py`)

The radial test left one obstruction: the OBE sources gravity from the
constraint *deviation* `ρ_N = |Π_H[Φ]|² − (1−14/D)`, which is **zero for leptons**.
Head-to-head across the three sectors (equivalence principle ⇒ charge/mass must be
universal):

| candidate | scalar map | lepton charge | charge/mass universal? | verdict |
|---|---|---|---|---|
| A: constraint deviation `a²(\|ξ\|²−½)` | Clifford grade-0 | **0** | no (2.4× spread; u/d=2.8 or 98 ≠ 40.6) | **DEAD** |
| B: second moment `Tr(M†M)=Σm=9Qa²` | generation trace (Frobenius²) | `1883` | **yes (=1)** | **LIVE** |

**Diagnosis of why A is dead:** `ρ_N` is the dimensionless *shape* `|ξ|²`, so it
**drops the amplitude `a`** — but mass `= 9Qa²` lives in `a²`. A scale-free shape
cannot be a mass.

**The live charge is the Frobenius² of the kernel** = the total mass, which is the
*same* Frobenius² structure as the bridge (`v = ‖Y‖²_F = dim(L)²a²`), over the
3-dim generation space instead of the 28-dim ambient. This gives a clean bonus
relation tying gravity to the Higgs VEV:

> `Σm_lepton = (9Q / dim(L)²)·v = (6/784)·v` — empirical 0.07%.

**Net on gravity:** the radial obstruction (1/r²) and the charge obstruction (EP)
are now *both* resolved by the consistent Frobenius²/scalar-map structure. **Still
open** (untouched here): the coupling magnitude (Newton `G`; V6 was ~10⁴⁰× strong)
and scalar-vs-tensor (LIGO).

## 5. Options 5 and 2 for the 28² factor (`obe_options_2_5.py`)

**Option 5 (operator-resolution algebra) — LIVE, and the cleanest derivation.**
Non-associativity forces `|Φ|²` onto the associative left-action algebra; for `L`,
the product *is* the bracket, so the operator algebra is `𝒜_L = ⟨ad_X : X∈so(8)⟩`.
Computed the dimension of the algebra generated:

| rep of L=so(8) | generated algebra dim | reading |
|---|---|---|
| defining (8×8) | **64** (=8²) | wrong reading |
| **adjoint (28×28)** | **784** (=28²) | `End(L)`, the bridge |

`ad(so(8))` is absolutely irreducible ⇒ by Burnside/Jacobson density it generates
the *full* `End(L) = M₂₈(ℝ)`, dim **784**. So **28² is the dimension of the
associative algebra that resolves the non-associative L-action** — a derivation,
not a posit. The adjoint reading is the physical one (the connection `Ω∈L` acts on
`Φ∈L` by the bracket); the defining rep (→64) is the wrong one.

**Option 2 (composite condensate) — LIVE as a counting mechanism.**
`‖Y‖²_F = (# nonzero entries)·a²`: the so(8)-singlet (`Y∝δ`, 28 entries) gives
`28a²` (the excluded linear reading); a symmetry-breaking democratic condensate
(784 entries) gives `784a²` (the bridge); a rank-1 single-field outer product gives
the wrong power (`784b⁴`). So 784 needs **two independent 28-legs** (`⟨ψ̄_L^a ψ_R^b⟩`)
or a broken vacuum. Live for the count, but the dynamical scale `a_ℓ` is unproven
(gap equation) and it implies a collider-testable compositeness scale; the broken
vacuum carries the 25 extra Goldstone directions (the "25 unaccounted" issue).

## 6. Gravity coupling magnitude (`gravity_magnitude_test.py`)

In the OBE, `G_N = f_g²/4π`, so the magnitude question is "what fixes `f_g`?"
The target is `α_G(e) = G_N m_e²/ℏc = (m_e/M_Pl)² = 1.752×10⁻⁴⁵`. v59's conjecture:
`G_e = (21/16)·α²¹` with `α = α(0)` (IR).

- **Numerically live; the exponent is the real signal.** `(21/16)α²¹ = 1.756×10⁻⁴⁵`,
  ratio 1.0025 (0.25%). The exponent that fits prefactor 21/16 is **n = 21.0005** —
  exactly dim Spin(7) to 0.002%, *not* merely "near 21." So `f_g ~ α^(21/2)`, and the
  V6 ~10⁴⁰ overshoot is precisely the missing `α²¹ ~ 10⁻⁴⁵`. Hierarchy form:
  `m_e/M_Pl = √(21/16)·α^(21/2)`.
- **IR-consistent:** uses `α(0)`; `α(M_Z)` misses by 4.2× — correct for long-range.
- **Universal G_N:** μ, τ, p, W, Z, top all match at the same 1.0025 (`α_G ∝ m²`);
  `α²¹` is the electron-electron value, others follow by mass².
- **Combination mechanism (hypothesis):** gauge `g_W²=5√α` is *additive* over
  `5 = 21−16` generators; gravity `g_grav² = (√α)^(2·21) = α²¹` is *multiplicative*
  over all 21 Spin(7) generators — which would explain why gravity is uniquely weak
  (it needs all generators coherently). No Lagrangian derives the product.
- **Numerology caveat (honest):** scanning integer exponents, only `n=21` gives an
  O(1) prefactor — that part is real. BUT the nearest simple rational to the fitted
  prefactor is **17/13 (0.1%), which beats 21/16 (0.25%)** — so "21/16 =
  dimSpin7/dimCl31" is *not* uniquely selected, and 21 is over-determined
  (`=28−7=14+7=35−14`). Exponent 21 = strong signal; prefactor + mechanism = open.

**Verdict:** the magnitude is a *value conjecture* (like α itself), not a derivation
— but the exponent being exactly dim Spin(7) makes it the most law-like of the
gravity pieces, and it quantitatively cures the V6 overshoot.

## 7. Net

**Which mechanisms are live (this round of tests):**

- **Option 5 (algebra dimension) — LIVE, strongest.** `784 = dim End(L)` is
  *generated* by `ad(so(8))` (numerically confirmed; defining rep gives 64). The
  bridge factor is derived as the dimension of the associative algebra that
  resolves the non-associative L-action — no posit.
- **Gravity charge — RESOLVED to the second moment.** The constraint-deviation
  source is DEAD (zero for leptons, fails EP); the live charge is `Tr(M†M)=Σm`,
  EP-exact, and the *same* Frobenius² as the bridge. Bonus: `Σm_ℓ=(6/784)v` (0.07%).
- **1/r² — reachable** (massless `K` + monopole charge); both the radial and
  charge obstructions are now closed by the Frobenius²/scalar-map structure.
- **Options 6, 7, 2 — live, weaker.** 7 separates the 1/r² zero mode from the
  28-count (internal volume); 6 gives masses as `3×28` overlap singular values
  (lepton-specific `√v=(D/3)Σ√m`); 2 supplies 784 via two fermion legs. All three
  reproduce the numbers but leave a residual (6/7: relocate the Brannen shape /
  need explicit `M_L`; 2: dynamical `a_ℓ` + compositeness scale).

**Recommended OBE edits** (folded into `NEW_OBE_FORMULATION.md §3a, §6`): (i) make
the gravity sector an explicit massless wave `□Ω_grav = f_g·ρ_grav` with charge
`ρ_grav = Tr(M†M)` (the second moment / mass density), **not** the constraint
deviation; (ii) ground the EW factor 28² in `dim End(L)` (Option 5) populated by
the mass bilinear, with `a_ℓ` the per-mode quantum `√v/dim(L)`.

- **Gravity magnitude — numerically live, mechanism open.** `G_e=(21/16)α²¹` cures
  the V6 overshoot exactly; the exponent is `21.000 = dim Spin(7)` (strong), the
  prefactor and the additive-vs-multiplicative generator mechanism are not
  established (17/13 beats 21/16; 21 over-determined).

**Still open after this round:** the gravity **mechanism** (a Lagrangian for the
`α²¹` / the 21-generator product), the **scalar→tensor** upgrade (LIGO `h=±2`), and
the dimensionful R1 identification `EW scale = ‖Y‖²_F`.
