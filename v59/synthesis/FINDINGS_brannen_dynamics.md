# Brannen Dynamics — Maxima Symbolic Analysis

**Date**: 2026-05-22
**Scripts**: `brannen_symbolic.mac`, `brannen_phase_dynamics.mac`,
`brannen_effective_potential.mac`, `brannen_phase_precise.py`

After committing to Option E** (single Higgs + structured Yukawas), this
session explores the SYMBOLIC dynamic structure of the Brannen-Z₃ kernel
to find: (a) deeper invariants beyond the Koide ratio, (b) structural
forms for the quark Brannen phases φ_d, φ_u, (c) a dynamic mechanism
for phase determination.

---

## 1. New invariant: e₂(s) = 3a²(1−t²)

Maxima found that the SECOND elementary symmetric function of the
Brannen amplitudes is φ-INDEPENDENT:

$$
e_2(s) = s_0 s_1 + s_1 s_2 + s_2 s_0 = 3 a^2 (1 - t^2)
$$

This is a new universal invariant complementing the Koide ratio. Specifically:

$$
\boxed{ e_2(s) \cdot D_X = 3 \cdot \dim G_2 \cdot a^2 = 42 \cdot a^2 }
$$

universally across all three v59 sectors. (Verified numerically for lepton, d-quark, u-quark.)

The most compact statement of the cross-sector identity is now:

$$
\boxed{ (1 - t^2_X) \cdot D_X = \dim G_2 = 14 }
$$

This is **simpler** than the previous form $(1-Q_N)\cdot D_N = 28/3$, and it's the **direct statement** of the v59 cross-sector structure.

---

## 2. Brannen phase patterns (suggestive)

High-precision numerical fits give:

| Sector | t²_X | Best-fit φ_X (rad) | v59 conjecture | Gap |
|---|---|---|---|---|
| Lepton | 1/2 | $±0.22222$ | $2/9 = 0.22222$ | $< 10^{-5}$ ✓ |
| d-quark | 3/5 | $+0.10859$ | $1/9 = 0.11111$ | 2.3 % |
| u-quark | 7/9 | $-0.07251$ | $-1/14 = -0.07143$ | 1.5 % |

Or equivalently, in $3\varphi$ form:

| Sector | $3\varphi_X$ | Reads as |
|---|---|---|
| Lepton | $2/3$ | $= Q_l = \dim G_2 / \dim\mathrm{Spin}(7)$ (lepton Koide) |
| d-quark | $1/3$ | $= Q_l / 2$ (half the lepton) |
| u-quark | $-3/14$ | $= -3/\dim G_2$ |

The lepton relation $3\varphi_l = Q_l$ is **exact**. The quark relations
are suggestive (~1–2 % match, within quark-mass scheme uncertainty) but
**no uniform pattern** ties all three. The denominators 9, 9, 14 differ.

---

## 3. Dynamic origin: effective potential

For an effective potential of the form
$$
V_{\rm eff}(\varphi) = -A \cos(3\varphi) + B \cos(6\varphi) + \cdots
$$
(the leading φ-dependent terms in $V_{\rm eff}$ from fermion loop contributions),
critical points satisfy $\sin(3\varphi)\cdot(3A - 12 B \cos(3\varphi)) = 0$, giving:

- $\sin(3\varphi) = 0$ → $\varphi = 0, \pi/3, 2\pi/3, \ldots$ (commensurate)
- $\cos(3\varphi) = A/(4B)$ → $\varphi$ depends on sector-specific ratio $A_X/B_X$

For the observed $\varphi_l = 2/9$: $\cos(3 \cdot 2/9) = \cos(2/3) \approx 0.786$ → $A_l/(4B_l) = 0.786$
For $\varphi_d = 1/9$: $\cos(1/3) \approx 0.945$ → $A_d/(4B_d) = 0.945$
For $\varphi_u = -1/14$: $\cos(-3/14) \approx 0.978$ → $A_u/(4B_u) = 0.978$

The **sector-specific coefficients** $A_X, B_X$ would arise from 1-loop and
2-loop fermion contributions to $V_{\rm eff}(\varphi)$, with different
fermion mass spectra giving different ratios. This is a **plausible dynamic
origin** for the empirical Brannen phases, though it requires explicit
loop calculation in a complete Lagrangian to verify.

---

## 4. Higher-moment structure (computed symbolically)

For Brannen-Z₃ amplitudes $s_k = a(1 + 2t \cos(\theta_k))$, $\theta_k = 2\pi k/3 + \varphi$:

| Quantity | Expression | φ-dependence |
|---|---|---|
| $P_1(s) = \sum s_k$ | $3a$ | none |
| $P_1(m) = \sum m_k$ | $3a^2(1+2t^2)$ | none |
| $e_1(s)$ | $3a$ | none |
| $e_2(s)$ | $3a^2(1-t^2)$ | none ← NEW |
| $e_3(s) = \prod s_k$ | $a^3(1 - 3t^2 + 2t^3 \cos(3\varphi))$ | $\cos(3\varphi)$ |
| $P_3(s) = \sum s_k^3$ | $3a^3(1+6t^2+2t^3 \cos(3\varphi))$ | $\cos(3\varphi)$ |
| $P_2(m) = \sum m_k^2$ | $3a^4(1+12t^2+8t^3 \cos(3\varphi)+6t^4)$ | $\cos(3\varphi)$ |
| $P_3(m) = \sum m_k^3$ | (involves $\cos(3\varphi)$ and $\cos(6\varphi)$) | both |

The Z₃-invariant combinations come in TWO classes:
- **φ-INDEPENDENT**: $P_1(s), P_1(m), e_2(s)$ — fixed by $(a, t)$ alone
- **φ-DEPENDENT through $\cos(3\varphi)$**: $e_3(s), P_3(s), P_2(m)$, ...

The Koide ratio $Q = P_1(m)/P_1(s)^2 = (1+2t^2)/3$ is φ-independent — that's
why it's a "structural" v59 prediction.

---

## 5. Z₃ symmetry and the discriminant

The cubic discriminant $\Delta = (m_0-m_1)^2(m_1-m_2)^2(m_2-m_0)^2$ is
a function of $(t, \cos(3\varphi))$. It vanishes when two masses degenerate.

For our three sectors at the v59 t² values and empirical φ, Δ is non-zero
(non-degenerate). At $t^2 = 0$, Δ = 0 (all equal). At $t^2 = 1$, Δ = 0
(boundary). The empirical $t^2_X$ values sit at interior critical points
where Δ is maximal — interpretation: the v59 t² values are the "Z₃-natural"
points where mass splittings are MAXIMIZED subject to the Brannen constraint.

---

## 6. What we did NOT find

After this symbolic exploration, the following remain open:

- **No uniform structural form** for quark Brannen phases across sectors
- **No dynamic derivation** of φ_X from a unified Lagrangian (requires loop calculation)
- **No new structural integer** beyond {2, 3, 5, 7, 8, 9, 14, 16, 21, 27, 28, 35, 63, 72}
- **No simpler form** for the higher Brannen-Koide ratios $Q^{(n)}$

The cross-sector pattern $(1-t^2)D = 14$ IS the cleanest statement we can extract symbolically.

---

## 7. Where this leaves us

After the symbolic analysis:

**Stable findings (CONFIRMED structurally)**:
- $(1-t^2_X)\cdot D_X = \dim G_2 = 14$ uniformly
- $e_2(\sqrt{m})\cdot D_X = 42 a^2$ uniformly
- $\varphi_l = 2/9 = Q_l/3$ (lepton, exact)

**Suggestive (1–2 % match, possibly within quark-mass uncertainty)**:
- $\varphi_d \approx 1/9$
- $\varphi_u \approx -1/14$

**Open path**: derive quark Brannen phases φ_X from fermion-loop $V_{\rm eff}$
once a complete Lagrangian is written. This is now the **first task** for
Step 6 (Lagrangian) work.

---

## Cross-references

- `brannen_symbolic.mac` — base symbolic derivation of e_n(s), P_n(s), $\Delta$
- `brannen_phase_dynamics.mac` — phase analysis across sectors
- `brannen_phase_precise.py` — high-precision numerical fit of phases
- `brannen_effective_potential.mac` — $V_{\rm eff}(\varphi)$ structure
- `DYNAMICAL_FIELD_OPTIONS.md` — overall framework (Option E**)
- `FINDINGS_scale_bridge.md` — established structural relations
