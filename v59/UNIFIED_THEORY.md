# v59 Unified Theory — Complete Description, Equations, and Gaps

**Date**: 2026-05-25 · **Status**: consolidated theory document (replicable description,
not a chronological notebook). Sources: `SUMMARY.md`, `furey_construction/RIGOR_AUDIT.md`,
`PARAMETER_BUDGET.md`, `synthesis/NEW_OBE_FORMULATION.md`, `synthesis/FINDINGS_obe_bridge.md`,
and the Lean modules in `furey_construction/lean/`.

> **Reading rule (status tags used throughout):**
> **[thm]** machine-checked Lean identity (mathematics, certain);
> **[emp≈X]** empirical match of a structural number to data at precision X;
> **[conj]** ansatz with a structural integer, not derived;
> **[free]** genuine input.
> A meta-point that governs the whole theory: *even `Q=2/3` is not "derived physics"* —
> it is an empirical fact matched to a structural ratio. "Derived" applies to the algebraic
> identities; every physical claim is ultimately a match of data to structure. The honest
> questions are **(a) how tight, and (b) mechanism or fitted formula?**

---

## 0. One-paragraph statement

A single 64-dimensional algebra, the Furey color algebra
$\mathcal{A}=Cl(7)_{\text{even}}\cong\mathbb{C}\otimes\mathbb{H}\otimes\mathbb{O}$, is posited as
the source of all Standard-Model structure. Its grade decomposition fixes the fermion sectors
and their mass-space dimensions; an exceptional-group chain $G_2\subset \text{Spin}(7)\subset
\text{Spin}(8)$ fixes the lepton mass ratios (Koide $Q=2/3$) and the weak mixing angle
($\sin^2\theta_W=2/9$); the dimension of the operator algebra of the lepton grade fixes the
electroweak scale ($v=\dim(L)^2 a_\ell^2$); and gravity couples to the second moment of the mass
kernel with a strength $G_e=(21/16)\alpha^{21}$ set by a power equal to $\dim\text{Spin}(7)$. The
genuine inputs are **one dimensionful scale $a_\ell$ and the fine-structure constant $\alpha$**;
the rest of the lepton+EW+Higgs block is dependent. The quark-flavour, CKM, neutrino, and strong
sectors are **not** reduced.

---

## 1. The source algebra and its grading  **[thm]**

**Postulate 1.** The fundamental field is multivector-valued, $\Phi(x)\in\mathcal{A}$,
$$
\mathcal{A}=Cl(7)_{\text{even}}\;\cong\;Cl(6)\;\cong\;\mathbb{C}\otimes\mathbb{H}\otimes\mathbb{O},
\qquad \dim_\mathbb{R}\mathcal{A}=64 .
$$

**Grade decomposition** (of $Cl(7)_{\text{even}}=\Lambda^0\oplus\Lambda^2\oplus\Lambda^4\oplus\Lambda^6$ on $\mathbb{R}^7$):
$$
\underbrace{\dim\Lambda^0}_{1}+\underbrace{\dim\Lambda^2}_{21}+\underbrace{\dim\Lambda^4}_{35}+\underbrace{\dim\Lambda^6}_{7}=64,
\qquad \tbinom{7}{2}=21=\dim\text{Spin}(7),\;\;\tbinom74=35,\;\;\tbinom76=7=\dim\text{Im}\,\mathbb{O}.
$$

**The $L\oplus F$ bisection.** Define
$$
L=\Lambda^2\oplus\Lambda^6\ \ (\dim 28),\qquad F=\Lambda^4\ \ (\dim 35).
$$
*Derivation that this is the complex/real split (the lepton-grade forcing):* a simple even
$k$-blade squares to $(-1)^{k(k+1)/2}\mathbb{1}$. Hence
$$
B\in\Lambda^2,\Lambda^6\Rightarrow B^2=-\mathbb{1}\ (\text{complex structures}),\qquad
B\in\Lambda^4\Rightarrow B^2=+\mathbb{1}\ (\text{real/involutions}).
$$
Equivalently $L=$ skew $=\mathfrak{so}(8)$, $\{\Lambda^0\}\oplus F=$ symmetric. A genuine complex
structure $J$ ($J^2=-\mathbb1$) carrying the Brannen phase therefore must live in $L$.
(`BladeSquareSign.prod_sq`, `LeptonComplexStructure.*`, `LeptonRealityForcing.*` — all **[thm]**.)
Note the additive identity $28+35=63$ and $2\dim G_2=28=\dim\text{Spin}(8)$.

---

## 2. Fermion sector: three generations and the selection rule

**Three generations** = the $\mathbb{Z}_3\subset S_3$ cyclic subgroup of the triality outer
automorphism of $\text{Spin}(8)$; concretely the sedenion automorphism group
$\text{Aut}(\mathbb{S})=G_2\times S_3$ supplies an order-3 automorphism $\psi$ (verified:
$\psi^3=\mathbb1$, $\psi$ an automorphism on all $16\times16$ basis pairs — **[thm-verified]**),
the "$/3$" that appears in the phase. 

**The $\mathbb{Z}_2\times\mathbb{Z}_2$ selection rule** assigns each fermion sector $N$ a pair of
grade-bits $(B_L,B_F)$, i.e. a mass-space ambient $D_N$:
$$
\text{lepton}=(L)\ \Rightarrow D_0=28,\qquad
d\text{-quark}=(F)\ \Rightarrow D_1=35,\qquad
u\text{-quark}=(L\oplus F)\ \Rightarrow D_2=63 .
$$
The $L\oplus F$ bisection is the grade-mod-4 eigenspace split $\mu|_{\Lambda^k}=(-1)^{k/2}$
**[thm]**; *which* $N$ takes *which* bits is **not derived** (see Gap G2).

**lepton $=L$ is forced** (not assumed): the color $\mathbb{Z}_3$ and the color complex structure
$J_c=\gamma_0\gamma_5\in\Lambda^2$ split $\mathbb{R}^8=\mathbb{C}^4=\mathbf1\oplus\mathbf3$ (lepton
singlet $\oplus$ 3 quark colors); an explicit $\mathfrak{su}(3)$ (8 generators, full $A_2$
closure) annihilates the singlet and commutes with $J_c$, pinning the lepton complex structure to
$\pm$canonical (`ColorSU3.*` — **[thm]**).

---

## 3. Mass spectrum — the lepton "Kepler ellipse"

**The Brannen mass kernel.** On the 3-generation flavour space,
$$
\boxed{\,M=a\,(\,\mathbb{1}+\xi S+\bar\xi S^2\,)\,},\qquad S=\text{cyclic shift }(S^3=\mathbb1),
\quad \xi\in\mathbb{H},\ |\xi|^2=t^2 .
$$
$M$ is circulant; its eigenvalues are $\lambda_k=a\,(1+2t\cos\theta_k)$, $\theta_k=\varphi+\tfrac{2\pi k}{3}$, where the phase is $\varphi=\arg\xi$. Identify $\lambda_k=\sqrt{m_k}$:
$$
\sqrt{m_k}=a\,(1+2t\cos(\varphi+\tfrac{2\pi k}{3})),\qquad k=0,1,2 .
$$

**Derivation of Koide** (using $\sum_k\cos\theta_k=0$, $\sum_k\cos^2\theta_k=\tfrac32$):
$$
\sum\sqrt{m_k}=3a,\quad \sum m_k=3a^2(1+2t^2)\ \Rightarrow\
\boxed{\,Q\equiv\frac{\sum m_k}{(\sum\sqrt{m_k})^2}=\frac{1+2t^2}{3}\,}\quad\textbf{[thm: }\texttt{BrannenKernel.Q\_value}\textbf{]}.
$$
The maximal-mixing / constraint value $t^2=\tfrac12$ gives $Q=\tfrac23$ **[thm:
`koide_iff_constraint`]**, matched to the Lie-algebraic ratio
$$
\boxed{\,Q=\frac{\dim G_2}{\dim\text{Spin}(7)}=\frac{14}{21}=\frac23\,}\quad\textbf{[thm for the ratio; emp}\approx10^{-5}\textbf{ for the match]}.
$$

**The phase.**
$$
\boxed{\,\varphi=\frac{Q}{3}=\frac{2}{9}\,}\quad\textbf{[emp}\approx10^{-5}\textbf{]},\qquad
\text{the “}/3\text{” = sedenion }S_3\ \textbf{[thm-verified]};\ \text{the magnitude }Q\ \text{is residual}.
$$
*Falsification/limits:* $\varphi$ is fixed only mod the generation $S_3$ (only $\cos3\varphi$ is
observable); the invariant $\cos(2/3)$ is **not** $\pi$-rational, so $\varphi$ is **not** a
geometric/holonomy angle (`PhaseExclusions`, `PhaseAmbiguity` — **[thm]**). It holds at Koide
precision ($|2/9-\varphi_{\text{meas}}|\approx7\times10^{-6}$, $0.9\sigma$).

**Quark Koide** (extends the same form):
$$
t_N^2=1-\frac{\dim G_2}{D_N}=1-\frac{14}{D_N}\ \Rightarrow\ Q_N=\frac{1+2t_N^2}{3},\qquad
Q_d=\tfrac{11}{15},\ Q_u=\tfrac{23}{27}\quad\textbf{[emp}\approx0.3\%,\ \text{RG/scheme-dependent — weak]}.
$$

**Total mass (used by gravity, §7).** From the same kernel,
$$
\boxed{\,\sum_k m_k=\mathrm{Tr}(M^\dagger M)=3a^2(1+2t^2)=9Q\,a^2=6a^2\ (\text{lepton})\,}.
$$

---

## 4. Gauge sector

**Group from the chain.** $\text{Spin}(7)$ decomposes as
$$
\boxed{\,\text{Spin}(7)=G_2\times SU(2)_L\times SU(2)_R\times U(1)_{B-L}\,},\qquad
\underbrace{14}_{G_2}+\underbrace{3}_{L}+\underbrace{3}_{R}+\underbrace{1}_{B-L}=21\quad\textbf{[thm]} .
$$
$U(1)_Y$ is the Pati–Salam diagonal $Y=2T_3^R+(B{-}L)$.

**Weak mixing angle** (from the Pati–Salam matching with Killing indices $(c_W,c_{B-L})=(5,2)$):
$$
\boxed{\,\sin^2\theta_W=\frac{c_{B-L}}{c_W+2c_{B-L}}=\frac{2}{9}\,}\quad\textbf{[thm given (5,2); }c_{B-L}=2\text{ pinned]},\qquad
\cos^2\theta_W=\frac79=\frac{\dim\text{Im}\,\mathbb{O}}{9}=t_u^2\quad\textbf{[thm]} .
$$

**Couplings.** Each gauge coupling$^2$ is an embedding (Killing) index times $\sqrt\alpha$:
$$
g_W^2=5\sqrt\alpha,\quad g_R^2=5\sqrt\alpha,\quad g_{B-L}^2=2\sqrt\alpha,\quad
\frac1{g'^2}=\frac1{g_R^2}+\frac1{g_{B-L}^2}\Rightarrow g'^2=\tfrac{10}{7}\sqrt\alpha .
$$
Here $5=h^\vee(\text{Spin}(7))=\dim\text{Spin}(7)-\dim Cl(3,1)=21-16$ is the dual Coxeter number
**[thm]**; the $\sqrt\alpha$ *form* is **[conj]** (see Gap G3). $g_W=0.6535$ vs $0.6517$ (0.28%).

**Fine structure.** Two conjectured structural forms (consistent under RG running):
$$
\text{IR: } -\ln\alpha+2\alpha=\frac{\pi^2}{2}=\frac{8\pi^2}{\dim Cl(3,1)}\ \textbf{[conj]},\qquad
\text{EW: } \alpha(M_Z)=\frac{25}{324\pi^2}=\Big(\frac{5}{18\pi}\Big)^2\ \textbf{[conj, emp}\approx0.03\%]} .
$$
*Derivation of $\alpha(M_Z)$ (given the conjectures):* the SM tree identity $4\pi\alpha=g_W^2\sin^2\theta_W$
with $g_W^2=5\sqrt\alpha$ and $\sin^2\theta_W=2/9$ gives $\sqrt\alpha=5\cdot(2/9)/(4\pi)=5/(18\pi)$
**[thm arithmetic]** — load-bearing on $g_W^2=5\sqrt\alpha$.

---

## 5. The electroweak scale bridge  $v_{\text{Higgs}}=\dim(L)^2 a_\ell^2$

**The relation** (the lepton mass scale sets the EW scale):
$$
\boxed{\,v_{\text{Higgs}}=\dim(L)^2\,a_\ell^2=28^2\,a_\ell^2=784\,a_\ell^2\,}\quad\textbf{[emp}\approx0.07\%]},
\qquad\Longleftrightarrow\ \sqrt v=\dim(L)\,a_\ell,\quad \frac{\sum\sqrt{m_\ell}}{\sqrt v}=\frac{N_{\text{gen}}}{\dim(L)}=\frac{3}{28}.
$$
Equivalently the lepton Yukawa sum rule $\sum\sqrt{y_\ell}=(3/28)\,2^{1/4}$ **[emp≈0.03%]**.

**Why the count is $\dim(L)^2$, not $\dim(L)$** (the Frobenius² reading; `03_higgs_bridge_result.md`):
a single operator trace $\mathrm{Tr}(L_\Phi L_{\tilde\Phi})=\dim\cdot|\Phi|^2$ is *linear* and gives
the excluded $\sqrt{28}$ scaling; the bridge needs the **Frobenius² of a rank-2 mass bilinear**
$Y$ on $L$, whose $\dim(L)^2=784$ entries each $\sim a_\ell$ give $\|Y\|_F^2=784\,a_\ell^2$.

**Why the relevant space has dimension $784$ (Option 5 — derived this session, `obe_options_2_5.py`):**
non-associativity forces $|\Phi|^2$ onto the associative left-action algebra; for $L=\mathfrak{so}(8)$
the product *is* the bracket, so the operator algebra is $\mathcal{A}_L=\langle\mathrm{ad}_X:X\in\mathfrak{so}(8)\rangle$.
Since $\mathrm{ad}(\mathfrak{so}(8))$ is the absolutely irreducible 28-dim adjoint, Burnside/Jacobson
density gives
$$
\boxed{\,\mathcal{A}_L=\mathrm{End}(L)=M_{28}(\mathbb{R}),\quad \dim\mathcal{A}_L=28^2=784\,}\quad\textbf{[thm, numerically confirmed: defining rep}\to64,\ \text{adjoint}\to784]} .
$$
This **upgrades** the $784$ from "noticed coincidence" to "dimension of the algebra the L-action
generates." Residual (R1) why the EW scale $=\|Y\|_F^2$, and (R2) why each entry $=a_\ell$, remain
**open**.

**Downstream EW/Higgs masses** ($v\equiv v_{\text{Higgs}}$):
$$
m_W=\tfrac12 g_W v=\tfrac12\sqrt{5\sqrt\alpha}\,v\ (0.04\%),\quad
m_Z=\frac{m_W}{\cos\theta_W}=\frac{3}{\sqrt7}\,m_W\ (0.02\%),\quad
m_H\approx\sqrt{\tfrac{7}{27}}\,v\ (0.07\%) .
$$

---

## 6. The unified field equation (OBE)

**Postulate 2 (the connection).** A test excitation feels the effective bivector connection
$\Omega(x)\in\Lambda^2(\mathcal{A})$:
$$
\boxed{\ \Omega(x)=\int K(x,x')\Big[\,f_g\,\rho_{\text{grav}}(x')\;+\;f_W\,\mathcal{J}_L(x')\;+\;f_C\,\mathcal{J}_F(x')\,\Big]\,d^3x'\ }
$$
with
- $K(x,x')$ the propagation kernel — **massless** ($-\nabla^2K=\delta$) for the long-range
  (gravity, EM) sectors; massive $\Rightarrow$ Yukawa $e^{-mr}/r$ for the weak/strong sectors;
- $\rho_{\text{grav}}=\mathrm{Tr}(M^\dagger M)=\sum_k m_k$ the **second-moment mass charge** (§7);
- $\mathcal{J}_L=\langle\Phi\tilde\Phi\rangle_L$ (28-dim weak/spin current),
  $\mathcal{J}_F=\langle\Phi\tilde\Phi\rangle_F$ (35-dim color current);
- $f_W,f_C$ the gauge couplings of §4 ($g_W^2=5\sqrt\alpha$, …).

**Postulate 3 (the force law).** A test fermion in sector $N$ feels only its allowed grades:
$$
a_N=\big\langle\Pi_N[\Omega]\,\Psi_t\big\rangle_1+\kappa\big\langle\Pi_N[\Omega]\wedge\Psi_t\big\rangle_{\text{higher}},
\qquad
\Pi_0=\Omega_L,\ \Pi_1=\Omega_F,\ \Pi_2=\Omega_L+\Omega_F,
$$
enforcing the $\mathbb{Z}_2\times\mathbb{Z}_2$ selection rule (lepton→L, d-quark→F, u-quark→L⊕F).

The silent gauging $D_\mu\Phi=\partial_\mu\Phi+\tfrac{g}{2}[A_\mu^a\tau^a,\Phi]$ on the $\mathbb{H}$
sub-algebra reproduces $SU(2)_L$, with the 3 Goldstones of the $|\xi|^2=\tfrac12$ vacuum eaten to
give $m_W,m_Z$ at $v=784 a_\ell^2$.

---

## 7. Gravity

**(a) Charge — the second moment, not the constraint deviation** (`gravity_charge_test.py`).
The earlier candidate $\rho_N=|\,|\Pi_\mathbb{H}\Phi|^2-(1-14/D)\,|\sim a^2(|\xi|^2-\tfrac12)$ is
**dead**: it vanishes for leptons ($|\xi|^2=\tfrac12$) and is scale-free (drops $a$), so it cannot
be mass. The **live** charge is
$$
\boxed{\,\rho_{\text{grav}}=\mathrm{Tr}(M^\dagger M)=\sum_k m_k=9Q a^2\,}\quad(\text{Frobenius}^2\text{ of the kernel}),
$$
which is EP-exact (charge/mass $=1$ for all sectors) and the *same* Frobenius² structure as the
EW bridge (§5), over the 3-dim generation space. Bonus relation:
$$
\sum m_{\text{lepton}}=\frac{9Q}{\dim(L)^2}\,v=\frac{6}{784}\,v\quad\textbf{[emp}\approx0.07\%]} .
$$

**(b) Radial law — $1/r^2$ is reachable** (`obe_radial_test.py`). For massless $K$, integration
by parts gives $\Omega_{\text{grav}}=f_g\nabla(K*\rho_{\text{grav}})=f_g\nabla\Phi$, $-\nabla^2\Phi=\rho_{\text{grav}}$.
Because $\rho_{\text{grav}}$ has nonzero monopole ($=$ total mass), $\Phi\sim1/r$ (slope $-1.000$)
and the force $\sim1/r^2$ (slope $-2.000$). Recommended OBE form:
$$
\boxed{\,\Box\,\Omega_{\text{grav}}=f_g\,\rho_{\text{grav}}\,}\quad(\text{explicit massless wave}).
$$

**(c) Magnitude** (`gravity_magnitude_test.py`). With $G_N=f_g^2/4\pi$, the target is
$\alpha_G(e)=G_N m_e^2/\hbar c=(m_e/M_{\text{Pl}})^2=1.752\times10^{-45}$. The conjecture:
$$
\boxed{\,G_e=\frac{21}{16}\,\alpha^{21}\,}\quad\textbf{[conj, emp}\approx0.25\%]},\qquad
\frac{m_e}{M_{\text{Pl}}}=\sqrt{\tfrac{21}{16}}\;\alpha^{21/2},
$$
using $\alpha=\alpha(0)$ (IR; $\alpha(M_Z)$ misses by $4.2\times$). *The exponent fitting the
prefactor is $n=21.0005=\dim\text{Spin}(7)$ to $0.002\%$* — the strongest signal. So
$f_g\sim\alpha^{21/2}$, and the historical V6 $\sim10^{40}$ overshoot is exactly the missing
$\alpha^{21}$ suppression.

**(d) Mechanism hypothesis (combination).** Gauge couplings are *additive* over generator subsets
($g_W^2=5\sqrt\alpha$, $5=21-16$); gravity is conjectured *multiplicative* over **all** 21
$\text{Spin}(7)$ generators, $g_{\text{grav}}^2=(\sqrt\alpha)^{2\cdot21}=\alpha^{21}$ — which would
explain why gravity alone is so weak. **Not derived** (no Lagrangian; and the prefactor is
ambiguous — $17/13$ fits better than $21/16$).

---

## 8. Parameter budget

| block | optimistic (accept conjectures) | rigorous (theorems + tight emp only) |
|---|---|---|
| charged-lepton + EW + Higgs | **1 input: $a_\ell$** | **2 inputs: $a_\ell$ and $\alpha$** |
| + gravity scale | $G_N$ from $(21/16)\alpha^{21}$ (conj) | $+1$ (the $\alpha^{21}$ conjecture) |
| full Standard Model | **~10–12 free** (quark scales ×2, quark phases ×2, CKM ×4, neutrinos ≥7, $\theta_{QCD}$, $\alpha_s$ — all untouched) | same |

Dependent (turned structural this program): the 2 non-scale lepton d.o.f. ($\to Q$),
$v_{\text{Higgs}}$, $\sin^2\theta_W$, $\cos^2\theta_W$, $m_W$, $m_Z$, $m_H$, the gauge integers
$(5,2)$.

> **Gap-attack update (2026-05-25, `gaps/SYNTHESIS.md`):** two of these are demoted. The
> $g_W^2=5\sqrt\alpha$ "law" is **not independent** — it is $\alpha(M_Z)=25/(324\pi^2)$ in
> disguise (so the rigorous input set is exactly $a_\ell+\alpha$; the $\sqrt\alpha$ adds nothing).
> The quark Koide ratios are a **scale-convention artifact** (2–5% at a physical scale, not
> 0.3%) and are dropped from "dependent." See §9.

---

## 9. Where the gaps live

Ordered by how load-bearing they are. **Statuses below incorporate the 2026-05-25 gap attack
(six agents, one per cluster; full detail and Lean/Python backing in `gaps/` + `gaps/SYNTHESIS.md`).**

### A. Load-bearing conjectures (the theory's "1–2 input" headline rests on these)

| # | gap | current status | what it needs |
|---|---|---|---|
| **G1** | $v_{\text{Higgs}}=\|Y\|_F^2$ (R1) and entry $=a_\ell$ (R2) | Option 5 derives the space is $\dim\mathrm{End}(L)=784$ **[thm, shown generic]**; R1/R2 still **[conj]**. **New (rank tension):** the EW count needs a *full-rank* $Y$ (784 entries $\sim a$), but the 3-generation spectrum is *rank-3* $Y$ with $\|Y\|_F^2=\sum m=$ the **gravity charge**, not $784a^2$. The same $Y$ can't be both; the "25 unaccounted directions" is the symptom, and no single-step $so(8)\to H$ gives 3 light directions. | resolve the rank tension (full-rank EW vacuum vs rank-3 spectrum); identify the `XiVacuum` Higgs with the active block of $Y\in\mathrm{End}(L)$ |
| **G2** | the $\sqrt\alpha$ in $g_W^2=5\sqrt\alpha$ | **DEAD as an independent law** — given $4\pi\alpha=g_W^2\sin^2\theta_W$ with $\sin^2\theta_W=2/9$, $g_W^2=5\sqrt\alpha\iff\alpha(M_Z)=25/(324\pi^2)$ (one point). The "5"$=h^\vee$ survives **[thm]**; the $\sqrt\alpha$ carries no new info. | nothing — **G2 collapses into G3** |
| **G3** | $\alpha$ itself (both conjectured forms) | **[conj]**; the only dimensionless input. Tight within a form-family, but large look-elsewhere across families (precision alone isn't evidence). | best avenue = an **RG / β-function fixed point** that standard running passes through at $M_Z$ (would explain away the $5\sqrt\alpha$ coincidence); test: does one scale jointly fit $\sin^2\theta_W$ and $\alpha(M_Z)$ |

### B. Undriven mechanisms (structurally observed, not derived)

| # | gap | status | what it needs |
|---|---|---|---|
| **G4** | the $\mathbb{Z}_2\times\mathbb{Z}_2$ selection rule ($N\to(B_L,B_F)$) | pattern observed; lepton=L forced **[thm]**, but the quark assignments undriven | a dynamical reason each Furey $N$ takes its grade-bits (7+ hypotheses tried, none worked) |
| **G5** | quark Koide $Q_d=11/15,Q_u=23/27$ | **DEMOTED to soft pattern** — the "0.3%" is a *scale-convention artifact* (PDG mixed-scale only); at a single physical scale the gap is 2–5%, the RG spread is 16–27× the gap, and 12 rationals fit each band. | an RG-invariant reformulation, if one exists; otherwise honestly *not a prediction* |
| **G6** | quark phases $\varphi_u,\varphi_d$, quark scales $a_u,a_d$ | **[free]** ($\varphi_q\neq Q_q/3$). The $n\cdot\alpha$ numerology only tests on the weak $\cos3\varphi$; scale ratios (35, 72) are convention artifacts. | a quark analog of the lepton phase mechanism |
| **G7** | the magnitude $Q=2/3$ of the lepton phase | **Magnitude now LIVE** (two non-geometric readings: max-mixing G₂-content $t^2=(D-\dim G_2)/D$, and equipartition $t^2=\tfrac12$). $\cos(2/3)$ as an algebra invariant is **NULL**. **Reframing:** $\cos3\varphi$ *is* the standardized skewness of $\sqrt m$. | the **radian-insert** residual — why $Q$ re-enters as the phase *argument*; find a v59 character/eigenvalue giving $\cos(2/3)$ with $Q$ as a *weight* |

### C. Gravity-specific gaps

| # | gap | status | what it needs |
|---|---|---|---|
| **G8** | gravity **magnitude mechanism** | $G_e=(21/16)\alpha^{21}$ **[conj]**; exponent **data-forced to $21.0005=\dim\text{Spin}(7)$** (0.0024%). $\alpha^{21/2}=\det(\sqrt\alpha\,I_{21})=$ a $\Lambda^{21}$ top-form $=$ instanton $e^{-S}$, $S=21\ln(1/\alpha)$; the prefactor is the 0.27 sub-leading correction. Not derived; $21/16$ beaten by $17/13$. | a Lagrangian for the coherent 21-generator product. **Note: moot until G9.** |
| **G9** | **scalar vs tensor** | **FATAL & decisive — the top blocker.** Rigorous helicity decomposition: a scalar source ($\rho_{\text{grav}}$ a Lorentz scalar) gives $h=0$ only; the internal $so(8)/\Lambda^2$ index is **inert under spacetime rotations** (helicity is a spacetime label), so it supplies no $h=\pm2$. $\mathrm{Sym}^2$ has spin-2 reps but they are stranded (no soldering to spacetime Lorentz). **Fails LIGO.** | an induced/emergent-metric (Plebański-style) recast with $\Omega\in\Lambda^2$ the fundamental 2-form $B$ and $h_{\mu\nu}$ derived; tests = a simplicity/soldering constraint + a 2-TT-DOF count |

### D. Untouched sectors

| # | gap | status |
|---|---|---|
| **G10** | CKM (3 angles + phase) | only $\sin^2\theta_C\approx7\alpha$ conjectured (0.45%); rest **[free]** |
| **G11** | neutrino masses + PMNS (≥7) | **[unaddressed]**; $D_\nu$ fits no v59 ambient |
| **G12** | strong sector ($\theta_{QCD}$, $\alpha_s$) | **[unaddressed]**, with one **clean null** — $g^2=h^\vee\sqrt\alpha$ does *not* generalize ($g_s^2=3\sqrt\alpha\approx0.26$ vs measured 1.48), and one **live lead**: $\theta_{QCD}\approx0$ may be *forced* by the color complex structure $J_c$ **[thm]** (explains a zero, not fits a number) |
| **G13** | dynamical $\xi(x)$ → executable Lagrangian | framework only (3 Goldstones + 1 radial); the OBE (§6) is a *target structure*, not yet an Euler–Lagrange-derived field theory |

### Honest bottom line

The **algebraic skeleton is theorem-grade**: the grade structure, lepton=L forcing, the
exceptional-chain dimensions, the dual Coxeter $5$, $\sin^2\theta_W=2/9$ given $(5,2)$, the
sedenion $S_3$, and $\dim\mathrm{End}(L)=784$ as the algebra the L-action generates. Every
**physical** reduction — Koide, $v_{\text{Higgs}}$, $\alpha$, $G_e=(21/16)\alpha^{21}$ — is an
empirical match (tight: $10^{-5}$ for Koide/phase) **hanging on the conjectures G1, G3, G8**. The
"$\approx1$ input" headline is *conditional on those conjectures*; rigorously it is **one
dimensionful scale $a_\ell$ + $\alpha$ + a stack of unproven structural ansätze**.

**Revised frontier (post gap-attack), in order:**
1. **G9 (tensor mode)** — decisive; v59 gravity is scalar and fails LIGO. Test the induced-metric recast. *Gates the whole gravity sector (G8 is moot until this).*
2. **G1 rank tension** — the EW-bridge $Y$ (full-rank, 784) vs the spectrum/gravity-charge $Y$ (rank-3, $\sum m$) cannot be the same object; this is the real content of R1/R2 and links the EW and gravity sectors.
3. **G3 ($\alpha$)** via an RG fixed point (subsumes G2).
4. **G7** radian-insert; **G12** $\theta_{QCD}$-via-$J_c$.

Quark Koide (G5), CKM (G10), neutrinos (G11) are soft/open — patterns, not predictions.
