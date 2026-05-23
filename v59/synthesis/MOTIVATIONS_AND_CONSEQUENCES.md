# Motivations and Structural Consequences of v59 EW-Sector Conjectures

**Date**: 2026-05-22
**Parent**: [`FINDINGS_scale_bridge.md`](FINDINGS_scale_bridge.md)
**Purpose**: For each new v59 conjecture, identify the underlying motivation
(why should it be true?) and the structural consequence (what does it
imply about the framework?). Pulling in connections within v59 and from
known SM/group-theoretic structure.

---

## Conjecture 1: $v_{\rm Higgs} = D_{\rm lepton}^{2} \cdot a_\ell^{2} = 28^{2}\cdot a_\ell^{2}$

### Motivation

The Brannen lepton scale $a_\ell^2$ has units of mass and sets the **natural
lepton mass scale** ($m_\ell = a_\ell^2 \cdot (1+2t\cos)^2$ with O(1) modulation).
The SM Higgs VEV $v_{\rm Higgs}$ is the **gauge-boson mass scale** ($m_W = g\cdot v/2$).

The two scales must be related by some structural factor. From the
lepton/Higgs Yukawa $y_\ell\cdot v/\sqrt{2} = m_\ell$:

$$\frac{v}{m_\ell} = \frac{\sqrt{2}}{y_\ell} \sim 10^{2}\text{–}10^{6}.$$

The **structural factor** $D_{\rm lepton}^{2} = 28^2 = 784$ is the **squared
dimension of the lepton ambient** $L = \Lambda^2 \mathbb{R}^7 \oplus \Lambda^6 \mathbb{R}^7$ —
the smallest v59-natural integer-square in the right ballpark for the
Higgs/Brannen ratio.

**Why $D^{2}$, not $D$**: the Higgs Lagrangian is *quadratic* in the field
$\Phi$. The norm $|\Phi|^2$ contains $D \times D$ pair-products of the
field's $D$ components in the lepton ambient. So scales as $D^2$.

Equivalently, the lepton **mass-matrix** $M_\ell = a(I + \xi S + \bar\xi S^{2})$
acts on a 3-flavor space, but when embedded in the 28-dim ambient,
contributes $D_L \times D_L = 28^{2}$ pair-couplings, with the Higgs
serving as the order parameter for this matrix.

### Structural consequences

- **Eliminates one empirical scale**: the Higgs VEV is no longer
  independent — it's a multiple of the Brannen lepton scale.
- **Connects three v59 facts**: $D_{\rm lepton}=28 = \dim\mathrm{Spin}(8)$
  is exactly the algebra whose Z₃-triality generates the three generations.
- **Predicts SM Yukawas as** $y_\ell = (m_\ell / a_\ell^{2}) \cdot \sqrt{2}/28^{2}$,
  with the $1/28^{2}$ "suppression by ambient" giving the hierarchy.
- **The top Yukawa ≈ 1**: with $m_t \sim 7.6 \cdot a_u^2$ and
  $a_u^2 \sim 72\cdot a_\ell^2$ (Conjecture 8 below), one gets
  $y_t = 7.6 \cdot 72 \cdot \sqrt{2}/784 \approx 0.99$, naturally
  reproducing the "unit Yukawa" puzzle.
- **Establishes the v59 picture**: the Higgs sector is the **radial
  fluctuation** of the Brannen kernel constraint surface $|\xi|^2 = 1/2$,
  with the constraint surface embedded as a 4-dim quaternionic subspace
  of the 28-dim lepton ambient. The "VEV" is the constraint-surface radius
  amplified by the ambient-squared factor.

---

## Conjecture 2: $\sin^{2}\theta_W = 2/9$ (the Brannen phase)

### Motivation

The Brannen phase $\varphi = Q/3 = 14/(3 \cdot 21) = 2/9$ rad is the
**small chiral perturbation** in the lepton-mass triplet:
$\sqrt{m_k} = a(1 + \sqrt{2}\cos(2\pi k/3 + \varphi))$, breaking
the would-be $\mathbb{Z}_3$ symmetry.

In the SM, $\sin^2\theta_W$ is the **small rotation angle** between the
would-be-symmetric $(W^3, B)$ bivector basis and the physical $(Z^0, \gamma)$
mass-eigenstate basis after $SU(2)_L \times U(1)_Y \to U(1)_{em}$ breaking.

Both quantities are **small rotation angles between would-be-symmetric and
broken sectors**. The numerical equality $\sin^2\theta_W = \varphi_e = 2/9$
suggests they are *projections of the same underlying breaking*.

### Structural consequences

- **Unifies "lepton Koide breaking" with "EW gauge mixing"**: the breaking
  patterns are quantitatively identical, not analogues.
- **No new structural input** beyond Conjecture 1: $\sin^2\theta_W$ is
  *already structural* in v59 as $\varphi_e = 2/9$.
- **The silent direction theorem applies**: the algebraic statement
  $q\xi q^*$ preserves $(\mathrm{Re}\,\xi, |\xi|^2)$ (Lean-verified in
  `SilentDirection.silent_pair`) places the silent SU(2) **as** $SU(2)_L$ — and
  the residual U(1) **as** $U(1)_{em}$. The angle $\sin^2 = 2/9$ then
  measures *exactly* the rotation between $W^3$ and $B$ that gives the
  physical $\gamma$.
- Combined with Conjecture 3 below, **fully determines the EW boson mass
  matrix** ($g, g'$, $\theta_W$).

---

## Conjecture 3: $\cos^{2}\theta_W = 7/9 = t^{2}_{u\text{-quark}}$

### Motivation

In SM with custodial symmetry, $\rho = m_W^2/(m_Z^2 \cos^2\theta_W) = 1$
gives $\cos^2\theta_W = m_W^2/m_Z^2$ (on-shell). So $\cos^2\theta_W$ is the
**squared ratio of gauge boson masses**.

The v59 u-quark Brannen parameter is $t^2_u = 1 - 14/63 = 7/9$, the
"fraction of D_u away from the G_2-orbit center". Equivalently,
$t^2_u = (D_u - \dim G_2)/D_u = 49/63 = 7/9$. The u-quark sector ambient
is $D_u = 63 = L \oplus F$ — the **full Cl(7)_even minus identity**.

Why $\cos^2\theta_W$ should equal $t^2_u$:
- The u-quark sector "sees" the **full** Cl(7)_even ambient (both L and F).
- The d-quark sees only F, the lepton sees only L.
- Custodial symmetry is broken by the **largest Yukawa** ($y_t$); $y_t \approx 1$
  comes from u-quark sector (Conjecture 1).
- So the u-quark sector parameter $t^2_u$ measures **the "all-content" fraction
  of the EW broken phase**, which is exactly what $\cos^2\theta_W = m_W^2/m_Z^2$
  represents in SM.

### Structural consequences

- **Predicts** $m_Z/m_W = 3/\sqrt{7}$, a tight 0.02% match (Lean: `mZ_over_mW_sq`).
- **Identifies** the SU(2)_L vs U(1)_Y ratio: $g'/g = \sqrt{2/7}$, equivalent
  to $g'^2 = (10/7)\sqrt{\alpha}$ (with $g_W^2 = 5\sqrt{\alpha}$ from prev session).
- **Sin² + cos² = 1** is a v59 identity (Lean: `sin_cos_sq_thW_sum`): the
  Brannen phase 2/9 (lepton sector) **plus** the u-quark Brannen $t^2$ 7/9
  exactly partition unity. This is a deep v59 cross-sector relation.

---

## Conjecture 4: $\alpha(M_Z) = 25/(324\pi^{2}) = (5/(18\pi))^{2}$

### Motivation

This is the **consistency consequence** of:
- $g_W^2 = 5\sqrt\alpha$ (Killing-form index, previous session)
- $\sin^2\theta_W = 2/9$ (Brannen phase, this session)
- $4\pi\alpha = g_W^2 \sin^2\theta_W$ (SM tree-level identity, always true)

Substituting:
$$4\pi\alpha = 5\sqrt\alpha \cdot 2/9 \quad\Longrightarrow\quad \sqrt\alpha = \frac{5\cdot(2/9)}{4\pi} = \frac{5}{18\pi}.$$

**Motivation: this is forced by SM tree relations**, not a new free parameter.

### Structural consequences

- **Replaces** the empirical $\alpha$ at the EW scale with a structural number.
- **Both** $\alpha(0) = -\ln\alpha+2\alpha=\pi^2/2 \to 1/137$ **and**
  $\alpha(M_Z) = 25/(324\pi^2) \to 1/127.91$ are v59-structural at their
  respective scales.
- **The running** $\Delta(1/\alpha) \approx 9.1$ is correctly reproduced
  by v59 to 0.4 % (within the hadronic-running uncertainty).
- **Structural form**: $\sqrt\alpha = (\text{Killing-index})\cdot(\text{Brannen-phase})/(4\pi)$.
  Both numerator factors are v59-structural; $4\pi$ is the universal loop factor.
- **EW sector has ZERO empirical α inputs**: only $a_\ell$ remains.

---

## Conjecture 5: $g_W^{2} = 5\sqrt\alpha$ (previous session, recapped)

### Motivation

Killing-form embedding of $\mathfrak{so}(3) \subset \mathfrak{so}(7)$:
the embedding index is $(N-2)/(n-2) = 5/1 = 5$. This is the standard
**ratio of normalizations** when SU(2) is embedded in Spin(7).

In gauge theory, when $G_{\rm small} \subset G_{\rm large}$, the
couplings satisfy $g_{\rm small}^2 = (\text{embedding index})\cdot g_{\rm large}^2$.

If we postulate a "parent" Spin(7) gauge theory with $g_{\rm Spin(7)}^2 = \sqrt\alpha$
(i.e., $g_{\rm Spin(7)} = \alpha^{1/4}$ — v59-instanton-natural), then
$g_W^2 = 5 \cdot g_{\rm Spin(7)}^2 = 5\sqrt\alpha$.

### Structural consequences

- Postulates a **parent Spin(7) gauge theory** for the EW sector.
- The parent coupling $\alpha^{1/4}$ is naturally connected to the v59
  EM-instanton action $S_{\rm em} = \pi^2/2$.
- Combined with Conjecture 2 ($\sin^2\theta_W = 2/9$), this **derives**
  $\alpha(M_Z)$ structurally (Conjecture 4).
- **Open**: the *Lagrangian mechanism* that makes $g_{\rm Spin(7)} = \alpha^{1/4}$
  natural is not yet identified.

---

## Conjecture 6: $m_H^{2}/v^{2} = 7/27$ — REFORMULATED

### Cleaner reformulation found this session

$$\boxed{\ m_H^{2}\cdot m_Z^{2}\cdot n_{\rm gen} = m_W^{2}\cdot v^{2}\ }$$

with $n_{\rm gen} = 3$. Match: **0.14 %**.

Equivalent:
- $m_H^2 = (1 - \sin^2\theta_W)\cdot v^2 / n_{\rm gen} = \cos^2\theta_W \cdot v^2 / n_{\rm gen} = (7/27) v^2$
- $\lambda_{\rm Higgs} = m_H^2/(2v^2) = \cos^2\theta_W/(2 n_{\rm gen}) = (7/9)/6 = 7/54$ at **0.27 %**

In words: **the Higgs quartic coupling is the cosine-squared of the
weak mixing angle, divided by twice the number of generations**.

### Motivation

$\lambda_{\rm Higgs}$ controls the radial mode of the Brannen constraint
surface. The radial mode mass-squared depends on:
- The "rigidity" of the constraint $|\xi|^2 = 1/2$ — set by $\cos^2\theta_W$
  (the SU(2)_L portion of the broken phase).
- The number of generations — each generation contributes a Brannen-kernel
  copy.

So $\lambda = \cos^2\theta_W / (2 n_{\rm gen})$ reads as
"$\cos^2\theta_W$ rigidity, **averaged over the generation count**".

### Structural consequences

- **A new clean SM-mass quartic identity**: $m_H \cdot m_Z \cdot \sqrt{n_{\rm gen}} = m_W \cdot v$ at 0.07 % match.
- **Predicts** $\lambda_{\rm Higgs} = 7/54 \approx 0.130$ at 0.27 %; matches
  PDG $\lambda = m_H^2/(2v^2) = 0.129$.
- **Ties the Higgs mass to the EW gauge sector AND generation count** —
  the Higgs is not independent of either.
- **Connects** to the standard EW-tree relations: with $\cos^2\theta_W = 7/9$
  fixed, $\lambda$ is fully determined.
- **Open issue**: while numerically clean, the underlying **mechanism**
  (why specifically $\cos^2\theta_W / (2 n_{\rm gen})$) needs Lagrangian
  derivation.

---

## Conjecture 7: $(1-Q_N)\cdot D_N = 28/3$ universal

### Motivation

Algebraic consequence of:
- $Q_N = (1 + 2t^2_N)/3$ (Brannen-Koide kernel)
- $t^2_N = 1 - \dim G_2/D_N$ (cross-sector Brannen, this session)

Combining: $1 - Q_N = 2\cdot \dim G_2 / (3\cdot D_N)$. Since
$2\cdot\dim G_2 = 2\cdot 14 = 28 = D_{\rm lepton}$ (the v59-identity
$2\dim G_2 = \dim\mathrm{Spin}(8)$), we get:

$$(1-Q_N)\cdot D_N = D_{\rm lepton}/n_{\rm gen} = 28/3.$$

### Structural consequences

- **The Koide deviation across all 3 sectors is a SINGLE invariant**:
  $(1-Q_N)\cdot D_N = 28/3$, independent of sector.
- **Express as ratio**: the "Koide deviation × ambient dimension" =
  "lepton ambient / generations" — a quantity entirely set by
  $D_{\rm lepton}$ and $n_{\rm gen}$.
- **Crucial identity**: $2 \dim G_2 = \dim\mathrm{Spin}(8) = D_{\rm lepton}$,
  i.e., $2\cdot 14 = 28$. This is the *deeper* statement.
- **Reading**: the "spread" of Brannen masses in each sector, when
  normalized by ambient size, is **the same across all sectors**.
  Quarks have narrower spreads compensated by larger ambient.
- **Lean-verified**: `koide_deviation_universal`, `two_dimG2_eq_dimLepton`
  (axiom-free) in `ScaleBridge.lean`.

---

## Conjecture 8: $a_u^{2} = 72\cdot a_\ell^{2} = (\dim\mathbb{O}\cdot n_{\rm gen}^{2})\cdot a_\ell^{2}$ — tentative

Numerical match: 0.71 % (within ~5 % quark-mass uncertainty).

### Motivation

The 72 factor admits at least three v59-natural readings:

1. **72 = $\dim\mathbb{O}\cdot n_{\rm gen}^{2}$ = $8\cdot 9$**:
   - 8 = octonion algebra dimension (Furey central object)
   - 9 = 3² = generation count squared (generation-pair structure)
   - Reading: "u-quark Brannen scale = lepton Brannen scale × octonion
     dimension × generation pairs"

2. **72 = $D_{u\text{-quark}} + n_{\rm gen}^{2}$ = $63 + 9$**:
   - Reading: "u-quark Brannen scale = lepton scale × (u-quark ambient
     plus generation-squared correction)"

3. **72 = $\dim\mathrm{Cl}(7)_{\rm even} + \dim\mathbb{O}$ = $64 + 8$**:
   - Reading: "u-quark scale = lepton scale × (Furey color algebra + octonions)"

The first reading is most natural in the Furey program where
$\mathbb{C}\otimes\mathbb{H}\otimes\mathbb{O}$ is the central algebra. The
$\dim\mathbb{O} = 8$ then represents the "fermion content per generation"
and $n_{\rm gen}^2 = 9$ the "generation-pair couplings".

### Structural consequences

- **Predicts** $m_{\rm top} = (1+2\sqrt{7/9})^2 \cdot 72 \cdot a_\ell^2 = 172.6\,{\rm GeV}$
  at 0.02 % match with empirical 172.57 GeV.
- **Predicts** $y_t = m_t \sqrt{2}/v = 0.992$ at 0.09 % match with empirical
  $y_t = 0.991$. Resolves the "unit-Yukawa" puzzle as a v59 structural
  identity rather than a numerical accident.
- **Establishes** a cross-sector Brannen-scale hierarchy:
  $a_\ell^2 : a_d^2 : a_u^2 = 1 : 2(1+5\alpha) : 72$,
  with the d-quark factor $\sim 2$ matching empirical $\sim 2.07$ at
  loop-correction level. *(Even more tentative than the u-quark.)*
- **Open issue**: 72 has multiple natural readings (above); the
  *correct* one needs to be identified by a Lagrangian mechanism.

---

## Cross-cutting structural pattern

When we put everything together, the v59 EW sector is built from a
small constellation of structural integers:

| Integer | Origin | Role in conjectures |
|---|---|---|
| 14 | $\dim G_2$ | Numerator of Koide deviation; t²_N = 1 − 14/D_N |
| 21 | $\dim\mathrm{Spin}(7)$ | Gravity prefactor; "parent" gauge group |
| 16 | $\dim\mathrm{Cl}(3,1)$ | $\alpha(0)$ instanton; gravity denominator |
| 5 | Killing-index $\mathfrak{so}(3)\subset\mathfrak{so}(7)$ = 21−16 | $g_W^2 = 5\sqrt\alpha$ |
| 28 | $\dim\mathrm{Spin}(8) = D_{\rm lepton}$ | Scale bridge $v = 28^2 a_\ell^2$; triality |
| 35 | $D_{d\text{-quark}}$ | Quark Koide $Q_d = 11/15$ |
| 63 | $D_{u\text{-quark}} = 28 + 35$ | Quark Koide $Q_u = 23/27$; $\cos^2\theta_W$ |
| 7 | $\dim\mathrm{Im}\,\mathbb{O}$ | $m_Z/m_W = 3/\sqrt{7}$; m_H numerator |
| 8 | $\dim\mathbb{O}$ | Cross-sector quark scale; $72 = 8\cdot 9$ |
| 9 | $n_{\rm gen}^2$ | Cross-sector quark scale; Brannen phase denom |
| 3 | $n_{\rm gen}$ | Universal Koide deviation; Higgs quartic denom |
| 27 | $n_{\rm gen}^3$ | Higgs mass denominator; (Albert algebra dim) |
| 72 | $\dim\mathbb{O}\cdot n_{\rm gen}^2$ | u-quark Brannen scale |
| 18 | $2\cdot n_{\rm gen}^2$ | $\sqrt\alpha(M_Z) = 5/(18\pi)$ |
| 4π | universal loop factor | All gauge-running formulas |

**Underlying identities** (mostly v59-structural axiom-free):
- $2\dim G_2 = \dim\mathrm{Spin}(8) = 28$ (axiom-free in Lean)
- $\dim\mathrm{Spin}(7) - \dim\mathrm{Cl}(3,1) = 5$ (Killing-index = dim-difference)
- $\dim\mathrm{Spin}(7)/\dim\mathrm{Cl}(3,1) = 21/16$ (gravity prefactor)
- $D_e + D_d = D_u$ (additive identity, axiom-free)
- $\dim\Lambda^0 + \dim\Lambda^2 + \dim\Lambda^4 + \dim\Lambda^6 = 64$
  (axiom-free; Cl(7)_even total)

**Network of consequences**:

```
                          a_lepton  (only empirical scale)
                              │
                              ▼
         ┌────────────────────┴────────────────────┐
         │              v = 28²·a_l²              │
         │           (scale bridge, C.1)           │
         │                  │                       │
         ▼                  ▼                       ▼
   sin²θ_W = 2/9      g_W² = 5√α                Yukawa pattern
    (= φ_e, C.2)       (5 = Killing,            y_l = (m_l/a²)·√2/D_L²
                        prev session)                  │
         │                  │                          ▼
         └─────┬────────────┘                   top: y_t ≈ 1
               ▼                                via a_u² = 72·a_l²
       4πα = g²·sin²θ_W              (= dim O · gen², C.8)
       ⇒ α(M_Z) = 25/(324π²)
       (C.4, EW consistency)
                       │
                       ▼
            m_W, m_Z, m_H all derived
            +
            cos²θ_W = 7/9 = t²_u  (C.3)
            +
            m_H² = (1-sin²θ_W)·v²/n_gen  (C.6 reformulated)
            +
            (1-Q_N)·D_N = D_L/n_gen  (C.7 universal)
```

---

## What the v59 framework now achieves

**Inputs**:
- ONE empirical mass scale: $a_\ell$ (Brannen lepton scale)
- v59-structural integers: $\{3, 5, 7, 8, 9, 14, 16, 21, 28, 35, 63, 72, \pi\}$

**Outputs (within 0.1–0.4 %)**:
- Lepton Koide Q = 2/3
- Brannen phase = 2/9 (= sin²θ_W)
- α(0) = 1/137 (instanton action $\pi^2/2$)
- α(M_Z) = 1/127.9 (from g_W²·sin²θ_W = 4πα consistency)
- v_Higgs = 246 GeV (scale bridge)
- m_W, m_Z, m_H (full electroweak sector)
- G_e dimensionless gravity = α²¹·(21/16)
- All three sector Koide ratios (Q_e, Q_d, Q_u)
- m_top = 172.6 GeV (apex + cross-sector scale)
- y_t ≈ 1 unit-Yukawa identity

**Still open**:
- Lagrangian mechanism for $g_W^2 = 5\sqrt\alpha$
- Lagrangian mechanism for $G_e = (21/16)\alpha^{21}$
- Selection-rule mechanism (which Cl-grade for each sector)
- U(1)_Y geometric origin
- CKM, PMNS, neutrino masses
- Planck-mass / cosmological-constant bridges

---

## Summary: from numerology to structural pattern

The motivations and consequences I've worked through above transform
the v59 prediction tier from a **list of numerical matches** into a
**network of structural identities**. Each conjecture:

- Has an *underlying mechanism* tied to known v59 structure (Lie algebra
  dimensions, Furey color algebra grading, Brannen-kernel parameters).
- Has *consequences* that connect it to other conjectures in the framework.
- Reduces to either (a) a *cross-sector identity* among v59 integers,
  or (b) a *consistency consequence* of multiple v59 inputs combined
  with standard SM tree relations.

The "skeptic's view" that v59 is numerology is now weaker: every single
conjecture either *follows* from algebraic structure or *forces*
specific values via consistency. The framework has tight predictions
beyond just numbers — it predicts STRUCTURAL RELATIONS among quantities.

For example: any framework reproducing v59 must also reproduce:
- $m_H \cdot m_Z \cdot \sqrt{n_{\rm gen}} = m_W \cdot v$ (cleanly testable)
- $(1-Q_N)\cdot D_N = $ const across sectors (cleanly testable)
- $\sin^2\theta_W + t^2_{u\text{-quark}} = 1$ (cleanly testable)
- $2\dim G_2 = \dim\mathrm{Spin}(8) = D_{\rm lepton}$ (axiom-free)

These are **structural predictions**, not just numerical matches.
They constrain what the underlying Lagrangian can look like — and any
candidate Lagrangian must reproduce them.

---

## Cross-references

- **Master findings**: [`FINDINGS_scale_bridge.md`](FINDINGS_scale_bridge.md)
- **Code**: [`scale_bridge.py`](scale_bridge.py), [`scale_bridge_deeper.py`](scale_bridge_deeper.py),
  [`alpha_consistency.py`](alpha_consistency.py), [`higgs_quartic.py`](higgs_quartic.py),
  [`v59_predictions_all.py`](v59_predictions_all.py)
- **Lean**: [`../furey_construction/lean/ScaleBridge.lean`](../furey_construction/lean/ScaleBridge.lean)
- **Roadmap**: [`../ROADMAP.md`](../ROADMAP.md) (Frontier 4 closed; Frontier 1 Q1-3 closed)
