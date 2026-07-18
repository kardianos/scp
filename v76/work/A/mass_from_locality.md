# Deriving \(m = E_{\mathrm{bound}}/c^2\) without dualist matter

**Approach:** A  
**Date:** 2026-07-18  
**Status:** Round 1 sketch — definitions + consistency argument; not a closed theorem  
**Axioms:** `axioms_v0.md` Ax1–Ax9

---

## 0. Goal

Show that, under the locality-first axioms, the inertial mass of a lock **must** be the bound ledger divided by \(c^2\), where \(c\) is free-field locality — not an independent “amount of stuff” and not energy of matter *on* a stage.

---

## 1. What must be related

Three quantities about a stable lock \(L\) in its rest free-frame \(F_0\):

| Symbol | Meaning under monism |
|--------|----------------------|
| \(E_{\mathrm{bound}}(L)\) | Unlock / rest ledger of bound medium (Ax6) |
| \(m_{\mathrm{inertial}}(L)\) | Resistance of \(L\) to acceleration relative to free medium |
| \(c\) | Free-field locality bound in \(F_0\) (Ax3, Ax7) |

**Claim:** \(m_{\mathrm{inertial}} = E_{\mathrm{bound}}/c^2\).

If this holds, writing \(m\) for inertial mass recovers \(E = mc^2\) as a **relation among continuum ledgers**, not a dualist postulate about matter.

---

## 2. Units from locality alone

Free rods and clocks (Ax7) fix:

- time unit = free oscillator period (or free round-trip);
- length unit = free-signal path in one time unit at local \(c\).

Thus \([c] = \mathrm{L}/\mathrm{T}\) is **purely free-field operational**. Energy is the ledger of reconfiguration capacity (Ax6); in units where free budget density is scaled so that free-signal action balances, \([E] = \mathrm{budget}\). Then

\[
\Bigl[\frac{E}{c^2}\Bigr] = \mathrm{budget}\cdot\frac{\mathrm{T}^2}{\mathrm{L}^2},
\]

which is the dimension of an inertial coefficient when “force” is free-budget flux rate into reconfiguring a lock’s exterior (see §4).

No Minkowski metric is required for this unit bookkeeping — only free-signal calibration.

---

## 3. Unlock energy = bound ledger

By Ax5–Ax6, converting \(L\) entirely to free form costs

\[
E_{\mathrm{unlock}}(L) = E_{\mathrm{bound}}(L;F_0)
\]

(in v0, with \(\mathcal{E}[\rho]=\rho\); more generally the integral of \(\mathcal{E}[\rho_{\mathrm{bound}}]\)).

Budget identity (Ax4): that energy is returned to \(\rho_{\mathrm{free}}\). There is no separate “rest mass energy of matter” stored in another sector.

**Monist identity:** rest energy **is** unlockable bound continuum, not energy *of* a particle *on* fabric.

---

## 4. Inertia as free-medium resistance (relational sketch)

### 4.1 Acceleration of a lock

To change the velocity of lock \(L\) relative to the surrounding free medium by \(\Delta v\) (\(|\Delta v|\ll c\)):

1. Free medium must be rearranged around the lock (bound structure’s causal envelope).
2. Free updates that reconfigure that envelope cannot exceed speed \(c\) (Ax3).
3. The amount of free budget that must be reorganized scales with how much medium is locked (size of bound ledger) and with how much the free-path structure must be skewed (kinematic factor involving \(v/c\)).

### 4.2 Soft argument (dimensional + Einstein-style)

In special relativity, the expansion of energy for a free particle is \(E = mc^2 + \tfrac12 m v^2 + \cdots\). Under monism we **reuse the same expansion for the free-frame ledger of a lock**, but justify the coefficients:

- Static lock: free-frame energy is pure bound ledger \(E_0 = E_{\mathrm{bound}}\).
- Slowly moving lock: free signals that maintain the lock’s stability condition must chase a moving bound region; the free-path cost of keeping \(\mathrm{Stable}(L)\) acquires a kinetic surplus \(\sim E_0 v^2/(2c^2)\).

Identifying the coefficient of \(\tfrac12 v^2\) as \(m_{\mathrm{inertial}}\) forces

\[
m_{\mathrm{inertial}} = \frac{E_0}{c^2} = \frac{E_{\mathrm{bound}}}{c^2}.
\]

### 4.3 What is still missing (honest)

A closed derivation needs:

1. a precise free-path cost functional \(C\) of free medium + locks — **prefer nonlocal free response** (C: local \(n(\rho_f)\) alone fails long-range tests);
2. a definition of “boost” of a lock as an automorphism of free-frame data;
3. expansion of free-medium cost of maintaining \(\mathrm{Stable}(L)\) at velocity \(v\) to \(O(v^2)\), yielding \(\tfrac12 (E_\star/c^2) v^2\).

Those are **dynamics + operational geometry** (Ax8–Ax9 + Approach B). C reverse (EQ checklist) requires unlock, rest, and inertia to share one ledger \(E_\star\). Until then, (M) in `axioms_v0.md` is:

- **definition** of mass as bound ledger / \(c^2\);
- **consistency requirement** that inertial coefficient match that definition;
- **exclusion** of any mass not tied to free-budget lock.

This is still monist-eligible: mass is not dualist substance; it is defined only from continuum ledger and locality.

---

## 5. Photon / free-packet limit (sanity)

A pure free radiation packet has \(\rho_{\mathrm{bound}}=0\), \(E_{\mathrm{free}}>0\). Then \(m_{\mathrm{rest}}=0\), while energy and momentum still relate by free-locality \(E = pc\) for null free packets. That matches: **rest mass is only for bound form**.

No second “massless matter field” is required — free field in flight is not a lock.

---

## 6. Link to warp (why mass without warp fails)

From Ax4: \(\rho_{\mathrm{free}} = \rho_{\mathrm{tot}} - \rho_{\mathrm{bound}}\) (strong form).  
From Ax7: observers enforce local free speed \(=c\).  
From Ax8–Ax9: free paths through depleted free budget acquire nontrivial global geometry.

If \(E_{\mathrm{bound}}>0\) on a compact lock, free budget is reduced there; free rays cannot see a globally flat free medium while keeping local \(c\). Hence:

\[
m>0 \;\Longrightarrow\; \text{warp in global free-signal chart}.
\]

Conversely, a theory with \(m>0\) on a fixed flat stage with globally constant free-path cost **contradicts** Ax4+Ax7+Ax8 — that is dualist fiction.

---

## 7. Killable claims (Round 1)

| ID | Claim | Kill if… |
|----|-------|----------|
| K1 | \(m = E_{\mathrm{bound}}/c^2\) with \(c=\) free locality | Inertia of lock fails to track \(E_{\mathrm{bound}}/c^2\) in any monist medium dynamics |
| K2 | Unlock energy = rest ledger | Stable locks require energy not accounted in free+bound budget |
| K3 | \(m>0\Rightarrow\) free-path warp | Medium with locks and local \(c=\mathrm{const}\) remains globally Euclidean for free rays |
| K4 | Strong budget (B) | Weak lensing / inertia force only integral or topological constraints (then weaken Ax4) |

---

## 8. Minimal formal fragment

See Lean: `work/A/lean/LocalityBudget.lean`

- Budget identity free + bound = total  
- Lemma: bound increase ⇒ free decrease (fixed total)  
- Def: mass = boundEnergy / c²  
- Lemma: mass monotone in bound energy (c fixed, c≠0)

That encodes the **ledger half** of the derivation; inertia expansion remains math notes until free-path cost is defined.
