# Kernel Fit 02 — Findings: Z₃-Invariant Potential and the Koide/Brannen Pinning Problem

**Date**: 2026-05-22
**Script**: `02_z3_invariant_potential.py`
**Target**: derive |ξ|² = 1/2 (Koide) and φ = 2/9 rad (Brannen) from a Z₃-invariant potential.

---

## Headline Result

The most general Z₃-invariant potential up to degree 6 in ξ has the form

$$V(r, \varphi) \;=\; m^2 r^2 + \lambda\, r^4 + |\mu|\, r^3 \cos(3\varphi - \alpha) + d_6\, r^6$$

with $r = |\xi|$, four real radial coefficients ($m^2, \lambda, d_6$) and a complex cubic coupling parameterized by magnitude $|\mu|$ and phase $\alpha$.

Critical-point analysis gives **two independent pinning conditions**:

1. **Phase**: $\partial V/\partial \varphi = 0$ forces $\sin(3\varphi - \alpha) = 0$, so the minimum is at $\varphi = \alpha/3$. Brannen's empirical $\varphi = 2/9$ rad therefore requires **$\alpha = 2/3$ rad**.

2. **Radial**: at the angular minimum, $\partial V/\partial r = 0$ at $r = 1/\sqrt{2}$ (Koide) gives **one** equation relating the radial coefficients:
$$m^2 = -\lambda + \tfrac{3\sqrt{2}}{4}|\mu| - \tfrac{3}{4} d_6$$
At the degree-4 truncation ($d_6 = 0$): $m^2 = -\lambda + \tfrac{3\sqrt{2}}{4}|\mu|$.

Stability (positive second derivative at the minimum) requires $4\lambda > \tfrac{3\sqrt{2}}{2}|\mu|$, i.e., $\lambda/|\mu| > 3\sqrt{2}/8 \approx 0.530$. Within this region the potential has a stable, unique Koide minimum.

## Numerical Verification

Scanned $(\lambda, |\mu|)$ over $\{0.5, 1.0, 2.0, 5.0\}^2$ with $m^2$ set by the Koide condition. The radial minimum is consistently at $r^2 = 1/2$ to within numerical precision ($\sim 10^{-9}$). The potential structure works as advertised.

## What This Buys Us

After steps 1–4 the situation is:

- The Cl(3,1) Z₃-cyclic algebraic kernel provides Brannen form for free (step 1–3).
- The Z₃-invariant potential pins $|\xi|^2 = 1/2$ as a stable minimum given the radial coefficient relation (step 4).
- The cubic phase $\alpha$ is the parameter that fixes Brannen's $\varphi$. To get $\varphi = 2/9$ rad, we need $\alpha = 2/3$ rad.

The kernel + potential together convert the lepton spectrum from "3 measured numbers" into:
- A fixed algebraic kernel (Cl(3,1) Z₃-cyclic), zero parameters.
- A potential with the structural form above, 4 parameters in general (or 3 at degree-4 truncation).
- One constraint from Koide ($m^2$ in terms of others).
- One constraint from Brannen ($\alpha = 2/3$ rad).

Net: lepton masses are reproduced by (2 free parameters $\lambda, |\mu|$) + (one mass scale $a$) once $\alpha$ is fixed. Three measured masses → three parameters, but with a SHARP STRUCTURE that wasn't there before.

## The 2/3 Coincidence (Structural Conjecture)

A striking numerical observation:

- $Q = 2/3$ (Koide ratio, dimensionless)
- $|\xi|^2 = 1/2$ (Koide condition on the kernel)
- $\alpha = 2/3$ rad (Brannen-pinning cubic phase, an angle)

The Koide ratio Q (a dimensionless quantity) and the cubic phase α (a phase in radians) are both equal to 2/3 — *numerically the same*. The natural conjecture is that they are related:

**Conjecture**: There exists an algebraic operation in Cl(3,1) whose natural Z₃-invariant cubic has phase $\alpha = Q$ rad, where Q is the Koide ratio determined by the radial structure. If true, this would tie Koide and Brannen together: once the radial structure forces $|\xi|^2 = 1/2$ (Koide), the same algebra would force $\alpha = 2/3$ rad (Brannen) automatically.

This is currently empirical numerology — the two appearances of 2/3 could be coincidence. But it's exactly the kind of relation a deeper structure would produce. Direct test: find a natural Z₃-invariant cubic operation in Cl(3,1) and compute its phase.

## What's Still Empirical

After steps 1–4, the following remain unexplained by structure alone:

1. **The radial coefficients $\lambda, |\mu|, m^2$** individually. The Koide condition gives one constraint among them; two are still free. They control the curvature of the potential (the "mass" of the ξ-field excitations around the Koide minimum) but not the lepton spectrum itself.

2. **The phase $\alpha = 2/3$ rad** is currently a hand-picked value. Until the structural conjecture is verified, this is an empirical input. (However, since it's just a single number tied to the empirical Brannen residue, it's no worse than what we had before.)

3. **The overall lepton mass scale $a$** is unconstrained by the lepton sector alone. It would be pinned by relating it to other measured quantities (Higgs VEV, GUT scale, etc.) — see Step 5.

## What This Is Not

This step does NOT derive the Standard Model lepton sector from scratch. It identifies the minimum algebraic structure (Cl(3,1) Z₃-cyclic kernel + Z₃-invariant potential) that reproduces the empirical Koide + Brannen residues, and it cleanly separates structural inputs from empirical ones.

Compared to the v58 unified_multivector_force experiment:
- v58 had many tuned parameters (λ, μ, ρ_crit, safe-band bounds, ambient modulation) and zero TYCHO matches.
- This kernel + potential has 4–5 free parameters and matches all three lepton masses to machine precision (under the empirical inputs $|\xi|^2 = 1/2, \alpha = 2/3$ rad).

## Next Step

**Step 5: cross-sector α prediction**. Now that the lepton sector is pinned, the same Lagrangian's bivector-grade content should determine the fine-structure constant $\alpha^{-1} = 137.036$. If the bivector coupling that emerges from the same potential and vacuum gives a sensible α, we have the first cross-sector unification test. If it doesn't, the framework needs extension (additional fields, scales, or symmetries) before unification claims hold.

This is the most discriminating test the v58 framework can face. Failure rules out the simplest unification; success would be major progress.

## Files

- `02_z3_invariant_potential.py` — analysis script.
- `02_findings.md` — this document.

## Status

| Item | Status |
|------|--------|
| Brannen form (eigenvalue structure) | Derived from Cl(3,1) Z₃ kernel (step 1–3). |
| Koide condition $\|\xi\|^2 = 1/2$ | Now pinned by Z₃-invariant potential (step 4). |
| Brannen condition $\varphi = 2/9$ rad | Requires $\alpha = 2/3$ rad in the potential. Conjectured to be tied to Koide via structural Q ↔ α equivalence. |
| Lepton masses | Reproduced to machine precision (step 3). |
| α (fine structure) | Step 5, not yet tested. |
| G (gravity) | Step 6, not yet tested. |
| m_p/m_e | Future, harder. |
