# Key Findings — Per-Paper Notes

Quick reference for the most important results from each paper.
Cross-references to INDEX.md by [number].

---

## Composition & Topology

### [01] Zheng et al. — Hopfion rings in FeGe (Nature 2023)
- **Key result**: Hopfion charge depends on LINKING between components, not just individual windings
- **Classification**: Ordered pairs (Q, H) where Q=skyrmion string charge, H=Hopf charge
- **Construction**: Hopfion ring = closed twisted skyrmion string (the string IS the building block)
- **Implication**: A single "particle" is a composite: skyrmion string + hopfion ring threading it
- **Observation method**: TEM + electron holography; emergent magnetic field gives H directly

### [02] Sutcliffe — JNR Skyrmions (2024)
- **Key result**: N+1 weighted points → charge-N Skyrme field (algebraic, not dynamical)
- **Energies within few percent** of numerically computed solutions for B=1–19
- **Symmetries preserved**: tetrahedral (B=3), octahedral (B=5), cubic (B=7), icosahedral (B=11)
- **Advantage over product ansatz**: respects broken symmetries; our scatter.c product ansatz does not
- **Implication**: better initialization for multi-soliton simulations

### [05] Chen, Yang, Li — Hopfion knot shapes (EPJC 2023)
- **Key result**: Hopf charge Q = W₁ × W₂ (toroidal × poloidal winding)
- **Constraint**: If W₁ < W₂, writhing CANNOT fully convert to twisting
- **Shape varies** with Euler angle θ — more writhed at smaller θ
- **No simple composition algebra**: knot shape ≠ hopfion shape at same Q

## Non-Local Effects

### [07] Ornelas et al. — Non-local quantum skyrmions (Nature Photonics 2024)
- **Key result**: Topology exists ONLY in entanglement between separated photons
- **Neither photon individually carries topological information**
- **Skyrmion number more robust to noise than standard entanglement measures**
- **Paradigm**: "topology has traditionally been thought to exist in a single local configuration; now it is nonlocal"
- **Implication for CHPT**: The degenerate sector (J,P) mediating bulk correlations could be the mechanism

### [09] Canfora et al. — Integral identities for solitons (PRD 2024)
- **Key formula**: ∫[f(x)T + (∂_j f)T_{ji}x_i]d^d x = 0 for ANY function f(x)
- **Standard Derrick (f=1)**: E₂ = E₀ + E₄ → our E₂ = E₄ (sigma model, confirmed)
- **Higher identities (f=r²)**: link moment of inertia, D-term, and rms radius
- **BPS solitons**: ALL identities become trivial (more constraints = less structural info)
- **Implication**: Non-BPS solitons carry richer non-local structural information

## BPS / Sextic Term

### [11] Adam et al. — Near-BPS Skyrmions (JHEP 2022)
- **Key result**: L₆ + V alone → E = E_FB exactly (Bogomolny saturation)
- **BPS solutions are COMPACTONS** (field reaches vacuum at finite radius!)
- **Zero classical binding**: multi-soliton composition is additive at BPS level
- **Binding appears at O(ε²)** when L₂ + L₄ added as perturbation
- **Near-BPS expansion**: small ε = (L₂+L₄)/L₆ gives controlled perturbation theory
- **Implication**: L₆ is the natural starting point; L₂+L₄ are corrections

### [12] Gudnason & Halcrow — Quantum binding energies (PLB 2024)
- **Key result**: ~7B vibrational modes per B-Skyrmion; zero-point energy CANCELS much of classical binding
- **B=4**: 6.99 MeV/nucleon (theory) vs 7.08 MeV (experiment) — remarkable!
- **Critical insight**: Making classical binding zero (BPS) is NOT sufficient — spin quantization alone gives 10× too much binding
- **Vibrational spectrum** is what our normal_modes.c computes for B=1 (1 mode); multi-B needs all ~7B modes
- **Implication**: Must compute full vibrational spectrum of composed hopfions

### [13] Carbon-12 in generalized Skyrme (2024)
- **Key result**: L₆ preserves alpha-clustering in B=12 Skyrmions
- **Energy ordering correct**: D₃ₕ (triangular) vs D₄ₕ (chain) with L₆ + pion mass
- **Alpha particles visible** as distinct substructures within the B=12 configuration
- **Implication**: L₆ doesn't destroy composition structure — composition IS preserved

## Vector Mesons / Effective Non-Locality

### [15] Huidobro et al. — Rho-meson Skyrme variant (JHEP 2025)
- **Key result**: K₀ drops from 1080 → 351 MeV (empirical ~240 MeV)
- **Only ONE rho-pi interaction term needed** (motivated by holographic reduction)
- **Mechanism**: ρ-meson mediates effective Yukawa-type non-local interaction between pions
- **Implication**: Vector mesons are the simplest way to add effective non-locality

### [16] Gudnason & Speight — Omega-meson crystals (JHEP 2024)
- **Key result**: Omega couples to baryon current via Wess-Zumino term
- **Effect**: REPULSIVE short-range force that "puffs up" soliton, reducing peak baryon density
- **Implication**: Core density depletion can be enhanced by vector meson couplings

## Hopfion Energy Scaling

### [VK bound — Vakulenko & Kapitanskii (1979, still foundational)]
- **E ≥ C|Q|^{3/4}** for Faddeev-Skyrme (Hopf charge Q)
- **Compare Skyrmions**: E ≥ C|B| (linear in baryon number)
- **Sublinear scaling** → composite hopfions ALWAYS energetically favored
- **Numerical**: H=1 unknot, H=2 linked unknots, H=5 first composite link, H=7 first trefoil knot
- **Implication**: Multi-hopfion states have STRONG binding (stronger than Skyrmions)

### [18] Garcia Salcines — Eidic Field Theory (SSRN 2025)
- **Key result**: Field theory with Hopf sectors + Faddeev-Skyrme rigidity
- **VK bound recovered**: E ≥ C|Q|^{3/4}
- **SM compatible**: Embeds into SU(3)_c × SU(2)_L × U(1)_Y
- **"Eidic section"**: orientational variable (not positional) → CP¹ map → S² target
- **Implication**: Closest existing framework to what CHPT needs for Hopf sectors

## Condensed Matter Realizations

### [21] Knapman et al. — Spacetime hopfions from braiding (CommPhys 2024)
- **Key result**: TWO routes to higher Hopf charge:
  1. **Braiding** skyrmion positions in time (inter-soliton non-local)
  2. **Internal excitations** (periodic breathing/elliptical modes of single skyrmions)
- **Hopf index tunable** by applied AC field
- **Implication**: Hopf charge can emerge from DYNAMICS of lower-order solitons — a hierarchy

### [20] Hopfion crystals (PRL 2025)
- **Construction**: Hopf map ∘ rational map (same technique as our verify3d.c B=2–4 init)
- **Lattice types**: simple cubic, FCC, BCC hopfion crystals
- **High charge**: torus links (Q_H=6), trefoil knots (Q_H=7)
- **Implication**: Our rational map infrastructure extends naturally to hopfion crystals

## Gravity & Solitons (from gravity search)

### Barcelo-Liberati-Visser — Effective metric from soliton background
- **Key formula**: g_eff^{μν} = -∂²L/∂(∂_μφ₀)∂(∂_νφ₀)
- **Pure L₂ (sigma model)**: g_eff = η (FLAT) — no gravity from quadratic kinetic term!
- **L₄ (Skyrme)**: g_eff is nontrivial — pions see curved spacetime near soliton core
- **Implication**: Need nonlinear kinetic terms (L₄, L₆) for effective gravity; L₂ alone gives nothing

### Gudnason & Speight — Derrick's theorem in curved spacetime (PRD 2019)
- **Key result**: Gravity weakens Derrick instability but CANNOT eliminate it
- **Virial with gravity**: I₃ = I₁ + 2I₂, forces I₂ = 0 for pure scalar → no stable config
- **Skyrme term still essential** for stability even with gravity
- **Implication**: "Gravity as field depletion" alone doesn't stabilize — topology still needed

### Bizon & Chmaj — Critical gravitational coupling
- **α²_max ≈ 0.0404** for B=1 gravitating Skyrmion
- **Two branches**: Skyrmion branch (lower mass) and Bartnik-McKinnon branch (higher mass)
- **Beyond α²_max**: no regular solution — soliton destroyed by its own gravity
- **Implication**: Self-gravity has a maximum; too deep a well destroys the soliton

### Dvali & Gussmann — Black hole Skyrmion hair (2016)
- **Key result**: B=1 Skyrmion outside a black hole = FIRST counterexample to no-hair for global charges
- **Baryon number measurable** via Aharonov-Bohm phase shift
- **Topology survives**: If Skyrmion is swallowed, B must resurface as classical hair when BH evaporates
- **Implication**: Topological charge transcends even black hole formation/evaporation

### Blazquez-Salcedo et al. — Gravitating Skyrmions with fermions (2024)
- **NEGATIVE ADM mass solutions** appear with Yukawa-coupled fermions
- **Anti-gravitating configurations**: regular, asymptotically flat, B=1 preserved
- **Spectral flow**: fermion eigenvalues cross at bifurcation points
- **Implication**: Self-gravitating soliton-fermion systems have richer behavior than simple "field depletion"

### Marino, Vocke et al. — Emergent geometry in photon fluid (Nature Sci Rep 2016)
- **EXPERIMENTAL**: Nonlinear density waves create self-induced background flow
- **Key quote**: "the initially flat spacetime geometry is curved by the waves, which are then distorted by the spacetime metric"
- **Acoustic line element**: ds² = (u²-c_s²)dt² + 2u dt dx + dx²
- **Shock = curvature singularity** of emergent metric
- **Strongest experimental evidence** for "gravity from field density variation"

### BPS Skyrme neutron stars (Adam, Sanchez-Guillen, Wereszczynski 2014)
- **BPS energy-momentum = perfect fluid**: T^{μν} = (ε+p)u^μ u^ν + pg^{μν}
- **Maximum NS mass**: few solar masses, R ~ 10 km (parameter-dependent)
- **L₆ crucial**: stiffens EoS at high density (incompressible BPS matter)
- **Without L₆**: EoS too soft, cannot support observed ~2 M_solar neutron stars
- **Implication**: L₆ needed for realistic gravitational structure of compact objects

---

## Critical Equations to Keep on Hand

### Energy scaling
- **Skyrmion**: E ≥ 6√2 π² ρ₀³/e × |B| (Faddeev-Bogomolny, linear in B)
- **Hopfion**: E ≥ C|Q|^{3/4} (Vakulenko-Kapitanskii, sublinear in Q)

### Virial relations
- **Standard (flat)**: E₂ - E₄ + 3E_V = 0
- **With gravity**: I₃ = I₁ + 2I₂ where I₃ = -E₀ - 2E_G
- **Mass formula**: Mc² = 2E₄ - 2E_V

### Effective metric from soliton background
- **g_eff^{μν} = -∂²L/∂(∂_μφ₀)∂(∂_νφ₀)** (Barcelo-Liberati-Visser)
- L₂ alone → flat; needs L₄ or L₆ for curvature

### Critical gravitational coupling
- **α² = 4πGF_π²**: dimensionless self-gravity parameter
- **α²_max ≈ 0.0404**: beyond this, no B=1 soliton exists

### BPS limit
- **L₆ + V only**: E = E_FB exactly, zero classical binding
- **Compacton**: field = vacuum for r > R_compact (finite support)
- **Near-BPS**: binding at O(ε²) from L₂+L₄ perturbation
