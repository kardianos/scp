# Tycho Table — Tight Numerical Mysteries for Substrate Fitting

**Date**: 2026-05-22
**Sources**: PDG 2024 (Particle Data Group, particle masses and mixings) and CODATA 2022 (fundamental constants). Where direct measurement is unavailable, values are quoted from the latest peer-reviewed determinations.

**Purpose**: This is the authoritative reference for what any v59-candidate substrate must reproduce. Each entry states a precise experimental relation, the inputs that feed it, the precision available, and the form a substrate-side prediction must take to count as a fit.

Entries are sorted into three tiers:
- **Tier 1**: dimensionless, theoretical orphan (not derived from the Standard Model), single-relation test of 5+ significant digits. These are the primary fitting targets.
- **Tier 2**: suggestive but with messier dependencies (multiple components, running corrections, or fewer digits). Secondary tests.
- **Tier 3**: open mysteries that remain too qualitative for clean fitting at the first stage. Listed for reference and future use.

For each Tier-1 entry the section also states the **fit prescription**: what a substrate must produce to count as a hit, and the precision the substrate prediction must reach.

---

## Tier 1: Dimensionless, Theoretical Orphans

### T1.1 — Koide Formula (Charged Leptons)

**Relation**:
$$
Q \;=\; \frac{m_e + m_\mu + m_\tau}{\big(\sqrt{m_e} + \sqrt{m_\mu} + \sqrt{m_\tau}\big)^{2}}
$$

**Experimental value**: $Q = 0.666661 \pm 0.000007$
**Conjectured exact value**: $Q = 2/3$.
**Agreement**: $\sim 5$ significant digits.

**Inputs (charged lepton pole masses)**:
- $m_e = 0.51099895069(16)$ MeV
- $m_\mu = 105.6583755(23)$ MeV
- $m_\tau = 1776.86(12)$ MeV

**Geometric content**: The vector $\big(\sqrt{m_e},\sqrt{m_\mu},\sqrt{m_\tau}\big)$ makes an angle $\arccos\sqrt{2/3} = 35.2644^\circ$ with the diagonal $(1,1,1)$. This is the angle between an edge of a cube and the body diagonal through the same vertex. Equivalently, the square-root masses lie on the "isoceles cone" whose axis is the symmetric diagonal direction.

**Why it qualifies as a Tycho number**: no derivation from the Standard Model. The formula is approximately MS-bar–invariant under renormalization but not exactly. The fact that it holds to 5 digits at low energy and not exactly under running is itself a clue: it points at a tree-level structural origin with small radiative corrections.

**Fit prescription**: A substrate passes T1.1 if it produces three natural invariants $\nu_1, \nu_2, \nu_3$ — eigenvalues of a canonical self-adjoint operator, ropelengths of an allowed three-knot family, mode energies of three localized excitations of the same algebraic type — such that
$$
\frac{\nu_1^{2} + \nu_2^{2} + \nu_3^{2}}{(\nu_1 + \nu_2 + \nu_3)^{2}} \;=\; \frac{2}{3} \;\pm\; 0.01
$$
on first pass, tightened to experimental precision ($\pm 10^{-5}$) on subsequent refinement. The $\nu_i$ play the role of $\sqrt{m_i}$. The substrate must select the three invariants without per-particle tuning.

### T1.2 — Brannen's Parametrization (Phase Offset)

**Relation**: With Brannen's ansatz,
$$
\sqrt{m_k} \;=\; \sqrt{m_0}\,\Big(1 + \sqrt{2}\,\cos\!\big(2\pi k/3 + \varphi\big)\Big), \qquad k = 0, 1, 2,
$$
the Koide formula holds identically for any $\varphi$ and any $m_0$. The experimental charged-lepton masses are reproduced at $\varphi \approx 0.2222$ rad.

**Experimental value**: $\varphi = 0.22222(2)$ rad.
**Suggestive value**: $\varphi = 2/9 = 0.22222\ldots$ rad (rational; agreement to 4 digits, possibly exact).

**Geometric content**: the square-root masses are three points on a circle, spaced 120° apart, offset by a small phase. The structure is 3-fold symmetric, with a small breaking parameter $\varphi$. A substrate with exact $\mathbb{Z}_3$ symmetry plus a small chiral perturbation could naturally produce this.

**Fit prescription**: A substrate that already passes T1.1 passes T1.2 additionally if the three invariants take the form $1 + \sqrt{2}\cos(2\pi k/3 + \varphi)$ with $\varphi$ within $\pm 0.01$ of the experimental value (initial pass), tightening to $\pm 10^{-4}$. Two passes against a single substrate is a much stronger result than one.

### T1.3 — Proton-Electron Mass Ratio

**Relation**: $\mu = m_p / m_e$

**Experimental value**: $\mu = 1836.15267343(11)$
**Agreement targeted**: 11 significant digits.

**Inputs**:
- $m_e = 0.51099895069(16)$ MeV
- $m_p = 938.27208816(29)$ MeV

**Why it qualifies**: $m_p$ arises from QCD confinement ($m_p \sim 3\Lambda_{\rm QCD}$, mostly binding, not bare quark mass), while $m_e = y_e v_{\rm Higgs}/\sqrt{2}$ arises from Yukawa coupling to the Higgs field. The ratio crosses two unrelated Standard Model sectors and is determined by independent parameters. There is no known reason it should take this specific value.

**Fit prescription**: A substrate that supports both a stable extended composite (proton analog) and a stable elementary mode (electron analog) passes T1.3 if the ratio of their natural invariants reproduces $1836.15$ to within 10% (initial pass), tightening toward experimental precision. Reproducing T1.3 with a substrate that also passes T1.1 would be a very strong result (it crosses sectors).

### T1.4 — Fine Structure Constant

**Relation**: $\alpha = e^2 / (4\pi\epsilon_0 \hbar c)$

**Experimental value**:
- $\alpha = 7.2973525693(11) \times 10^{-3}$
- $\alpha^{-1} = 137.035999084(21)$

**Agreement targeted**: 12 significant digits.

**Why it qualifies**: $\alpha$ is the most fundamental dimensionless coupling in QED. Many attempted geometric derivations (Eddington, Wyler, Atiyah, etc.) have failed to produce the experimental value without tuning. It is suspiciously close to combinations of $\pi$ and small integers but not exactly any of them. A substrate that produces a coupling constant from its geometric or algebraic invariants must reproduce it.

**Fit prescription**: A substrate passes T1.4 if a natural dimensionless ratio of its invariants (a homotopy index over a symmetry-group dimension, a volume ratio in a fundamental domain, a topological coupling extracted from an action) equals $1/137.036$ to within 1% on initial pass and to experimental precision on refinement.

T1.4 is the hardest single-number target because $\alpha$ has no obvious algebraic intuition. It is listed in Tier 1 because of its extraordinary precision and theoretical-orphan status; substrates that pass T1.1 should be tested against T1.4 as a stretch goal.

### T1.5 — Muon-Electron Mass Ratio

**Relation**: $m_\mu / m_e$

**Experimental value**: $m_\mu / m_e = 206.7682830(46)$
**Agreement targeted**: 9 significant digits.

**Why it qualifies**: same flavor, same charge, same conserved quantum numbers; difference is mass alone. Why this particular ratio?

**Fit prescription**: A substrate that supports at least two stable modes in the same algebraic sector (differing only by some discrete index — radial quantum number, level on a finite-group representation tower, knot complexity) passes T1.5 if the ratio of their natural invariants squared (treating the invariant as $\sqrt{m}$) reproduces 206.77 to within 5% on initial pass.

T1.5 is the lowest-cost first test: just two natural invariants of the same type, no need to identify a triple or solve for $\varphi$.

---

## Tier 2: Suggestive but Messier

### T2.1 — Tau-Electron Mass Ratio

$m_\tau / m_e = 3477.23 \pm 0.23$. Combined with T1.5 fully determines Koide for the charged leptons.

### T2.2 — Cabibbo Angle

$\sin\theta_C = 0.22500 \pm 0.00067$.

Suggestive Gatto-Sartori-Tonin relation (1968): $\tan\theta_C \approx \sqrt{m_d/m_s}$. Approximately true; not exact.

### T2.3 — Weinberg Angle

$\sin^2\theta_W = 0.23121(4)$ in MS-bar at $M_Z$.

Conjectured GUT-scale value $3/8$, runs down under RG flow. Connects U(1) and SU(2) gauge couplings.

### T2.4 — Top Quark Yukawa

$y_t = 0.9356 \pm 0.0036$ at $\overline{m}_t$ scale.

Suspiciously close to unity. Top quark mass sets the natural Higgs sector scale.

### T2.5 — Quark Mass Ratios

Six numbers (running masses, PDG 2024):
- $m_u = 2.16(7)$ MeV (MS-bar at 2 GeV)
- $m_d = 4.67(7)$ MeV (MS-bar at 2 GeV)
- $m_s = 93.4(8)$ MeV (MS-bar at 2 GeV)
- $m_c = 1.27(2)$ GeV (MS-bar at $m_c$)
- $m_b = 4.18(3)$ GeV (MS-bar at $m_b$)
- $m_t = 172.69(30)$ GeV (pole)

Less clean than leptons due to QCD running and non-perturbative confinement effects. Useful only after lepton sector is in hand.

### T2.6 — Lepton Anomalous Moments

- $a_e = (g_e - 2)/2 = 1.15965218059(13) \times 10^{-3}$
- $a_\mu = 1.16592061(41) \times 10^{-3}$

Leading order QED: $a = \alpha/(2\pi) \approx 1.16141 \times 10^{-3}$. Substrate-side: if $\alpha$ is fit, $a_e$ is constrained.

The muon $a_\mu$ has a $\sim 4\sigma$ tension with Standard Model prediction, possibly due to new physics.

### T2.7 — CKM Matrix Magnitudes

$|V_{us}| = 0.2243(8)$, $|V_{cb}| = 0.0410(14)$, $|V_{ub}| = 0.00382(20)$. PMNS matrix similarly. Pattern of magnitudes (Wolfenstein parametrization) suggests hierarchical structure.

---

## Tier 3: Exploratory / Open Mysteries

### T3.1 — Hadronic Mass Spectrum

- $m_p = 938.272$ MeV, $m_n = 939.565$ MeV, $m_n - m_p = 1.29333$ MeV (mostly EM + isospin breaking).
- $m_\pi^\pm = 139.57$ MeV, $m_\pi^0 = 134.98$ MeV, $m_K = 494$ MeV, $m_\eta = 548$ MeV, $m_\rho = 775$ MeV, etc.
- Nucleon mass arises ~95% from QCD binding, not constituent quarks. Substrates that produce binding-energy mass would be tested here.

### T3.2 — Cosmological Constants

- $\Lambda \sim 10^{-52}\,{\rm m}^{-2}$
- Cosmological constant problem: naive QFT estimate $\sim M_{\rm Planck}^4$ exceeds observation by 120 orders of magnitude.

### T3.3 — Three Generations

The qualitative observation that there are three families of fermions. Why three? Topological? Algebraic (e.g., octonion structure)? Anthropic? Open.

### T3.4 — Hydrogen Spectrum

- Rydberg constant $R_\infty = 10\,973\,731.568160(21)\,{\rm m}^{-1}$.
- Ground state binding $13.6056931229(3)$ eV.
- Predicted by Schrödinger equation given $m_e$, $e$, $\hbar$, $c$. A substrate deriving the hydrogen spectrum from first principles, without inserting $\hbar$ or QM, would be a definitive result.

---

## Priority Order for First Experiments

1. **T1.5** (muon/electron ratio) — only two invariants needed; lowest-cost first test.
2. **T1.1** (Koide formula) — primary target. Three invariants of the same type.
3. **T1.2** (Brannen $\varphi$) — refines T1.1 if it hits. Tests 3-fold symmetry hypothesis.
4. **T1.3** (proton/electron) — cross-sector test. Pursue only after a substrate passes T1.1 or T1.5.
5. **T1.4** ($\alpha$) — single-number stretch goal. Pursue against substrates that already pass T1.1.

Candidate substrate families (initial list — to be expanded in `first_experiments/`):

- **Clifford algebra eigenvalues**: eigenvalues of canonical self-adjoint elements in $\mathrm{Cl}(3,0)$, $\mathrm{Cl}(3,1)$, $\mathrm{Cl}(0,7)$.
- **Knot ropelengths**: minimum-energy lengths of unit-thickness rope tied in prime knots up to crossing number 9 (data from Pieranski, Cantarella, et al.).
- **Knot invariants**: Jones polynomial evaluations at roots of unity, hyperbolic volumes of knot complements, signature, Alexander polynomial roots.
- **Finite-group spectra**: irreducible representation dimensions and Casimir eigenvalues of $A_4$, $S_4$, $A_5$ (icosahedral), and their natural representations.
- **Quasicrystal vertex frequencies**: ratios of vertex-type populations in Penrose tilings and other aperiodic structures.
- **Octonion / Cl(0,7) invariants**: traces and norms of canonical octonionic structures, following the Furey program.
- **Lie group root-system invariants**: ratios of root lengths and weights of exceptional Lie groups $G_2$, $F_4$, $E_6$, $E_7$, $E_8$.

Each substrate family is one short Python session per test. Results — including null results — are logged in `first_experiments/`.

---

## A Note on Numerology

Numerology has a bad name because most of it is post-hoc fitting with many free parameters or no falsifiable predictions. v59's program is numerology in the Kepler sense: identify a tight, multi-digit, dimensionless relation in advance; fit a substrate to it with zero or one free parameter; require that the substrate be defined independently of the target. A substrate that happens to reproduce one Tycho number is interesting; one that reproduces two is a discovery; one that reproduces three is the start of a theory.

The right time to do this is exactly when many substrate proposals exist (the project has at least a dozen) and none of them have produced quantitative particle physics. Locking the targets first stops the search from being open-ended.
