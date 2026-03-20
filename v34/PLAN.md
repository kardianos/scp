# V34 Planning: EMF, Mass, Electrons

**Status**: Planning phase. These three problems are interrelated and may
need to be addressed together.

---

## The m=0 Collapse: What Actually Happens (V33-C4 Analysis)

The m=0 failure mode is NOT a localized collapse (black hole) or
dissipation. It is a **thermal explosion** — a Big Bang.

### Timeline (from V33-C4, m=0, two braids at D=20):

    t=0:   E_total = -136     Cold, low-energy initial state
    t=15:  E_total = -2551    Falls into potential well (D: 20→6.2)
    t=25:  E_total = +3650    BANG — potential energy → kinetic
    t=200: E_total = +62,816  Hot soup, E growing monotonically

### Phases:

**Phase 1 (t=0→15): Runaway attraction**
E_kin grows 200×. The attractive V(P) with μ<0, unchecked by any
m² spring, feeds energy into the field. Braids rush together.

**Phase 2 (t≈15→25): Explosion**
Braids collide, scatter violently (D: 6→50). E_total flips from
negative to positive. The potential is an unlimited energy source.

**Phase 3 (t>25): Hot chaotic thermalization**
No braids survive. E_kin ≈ 100,000, E_pot ≈ -88,000. E_total
grows at ~100/time. The field is a turbulent mess. Braid detection
picks up random hot spots (D jumps wildly: 6, 50, 30, 51, 33...).

### Where the instability lives:

    V''(P) = μ(1 - 3κP²) / (1 + κP²)³

    At P ≈ 0 (vacuum):  V'' = μ = -41.3   → TACHYONIC (m_eff² < 0)
    At P > 1/√(3κ):     V'' > 0            → STABLE

The instability is at SMALL P — the empty vacuum. κ saturation
only helps at large P. Increasing κ does NOT fix the vacuum.

### Cosmological interpretation (speculative):

The m=0 instability resembles the Big Bang:
1. Cold initial state → falls into potential → reheats → hot soup
2. V(P) plays the role of the inflaton potential
3. The tachyonic instability IS the false vacuum decay
4. The t≈25 transition IS reheating

If this is not a failure but a FEATURE:
- Early universe: m ≈ 0 → tachyonic instability → Big Bang
- Expansion/cooling → background forms → m_eff increases
- Once m_eff > stability threshold → vacuum stabilizes
- Braids condense from the cooling field → matter formation

The m=0 problem might be the ORIGIN STORY, not a bug.

### Metastability: m doesn't need to be zero

The vacuum doesn't need to be perfectly stable — just stable for
cosmological timescales. This is physically legitimate:
- Standard Model: our vacuum may be metastable (lifetime >> 10^100 yr)
- Proton: may decay, but not for > 10^34 years
- "Stable enough" IS the standard in real physics

V33-C4 shows a smooth stability transition:

    m² = 0.00:  E_drift = +46,000%  (explodes in T~20)
    m² = 0.25:  E_drift = -133%     (unstable, decays in T~200)
    m² = 1.00:  E_drift = -3.4%     (marginal — physical or numerical?)
    m² = 2.25:  E_drift = -0.3%     (stable, drift is numerical only)

The instability growth rate decreases smoothly with m². There exists a
METASTABILITY WINDOW: a range of m² where:
1. The vacuum is "stable enough" (instability timescale >> braid lifetime)
2. The gravitational range 1/m is much longer than at m=1.5
3. The force law is closer to 1/r² over relevant distances

**What we need**: A fine m² scan in the range [0.5, 2.0] that measures:
- Instability growth rate (separate numerical from physical drift)
- Braid survival time
- Force exponent n in F ∝ 1/D^n
- Yukawa range 1/m vs effective range from depletion (Track B)

If m²=1.0 gives n≈1.9 with braids lasting T=10,000+, that may be
sufficient. The force doesn't need to be EXACTLY 1/r² — just close
enough over the observable range, with Yukawa corrections at very
large D. This is TESTABLE and PRACTICAL.

Key: run the m²=1.0 case at different dt values (0.5×, 1×, 2× current).
If the -3.4% drift scales with dt² → numerical. If constant → physical
instability. This determines whether m²=1.0 is truly marginal or just
under-resolved.

---

## The Three Problems

### Problem 1: EMF / Photons (no Maxwell equations yet)

The theory has braids (matter) and gravity, but no electromagnetic field.
Need: photons, Coulomb force, magnetic fields, Maxwell's equations.

### Problem 2: m → 0 (Yukawa vs Newtonian gravity)

With m²=2.25, gravity is Yukawa (exponential decay, range ~0.67).
Newtonian 1/r² requires m→0, but m=0 gives vacuum instability (V33-C4).

### Problem 3: Electrons / Orbital Structure

A hydrogen atom needs: a heavy braid (proton), a light thing (electron)
orbiting it, quantized energy levels, and a single stable ground-state
orbital for the simplest case.

---

## Conceptual Framework

### EMF as a Stretched Braid

**The idea**: A braid is a localized, self-reinforcing helical pattern
(standing wave with topological winding). What if a photon is the SAME
type of self-reinforcing pattern, but EXTENDED rather than localized?

A braid, conceptually stretched to infinity along its propagation axis,
becomes a traveling helical wave — a circularly polarized electromagnetic
wave. The key differences:

| Property | Braid (particle) | Stretched braid (photon) |
|----------|-------------------|--------------------------|
| Localization | Compact (R~3) | Extended (λ = wavelength) |
| Topology | Winding W=±1 | No net winding (or W=0) |
| Propagation | Stationary/slow | Travels at c |
| Depletion | Binds field → δρ<0 | No net binding → δρ≈0 |
| Mass | m_grav > 0 (from binding) | m_grav = 0 (no depletion) |
| Energy | E_braid ∝ bound field | E = hν (from frequency) |

The helical structure is the SAME — three fields with phase offsets
oscillating along an axis. The difference is: the braid wraps this
helix into a compact envelope (Gaussian), while the photon lets it
propagate freely.

**What this gives us**:
- Photons are massless (no depletion → no gravitational mass) ✓
- Photons propagate at c (they ARE perturbations of the c=1 field) ✓
- Photons carry energy (field oscillation amplitude) ✓
- Photon polarization = helical handedness (left/right circular) ✓
- Photon frequency = oscillation frequency of the traveling helix ✓

**What this does NOT give us** (yet):
- Maxwell's equations (need U(1) gauge structure)
- Coulomb 1/r² force between charges
- Photon emission/absorption by braids
- Quantization of photon energy (E = hν)

**The emission question**: When a braid transitions between states,
it must release the energy difference as a propagating wave packet.
This wave packet IS the photon — a finite-length burst of helical
oscillation that travels outward at c. The braid's internal modes
(K≥1 angular modes, see V12 normal mode analysis) set the allowed
frequencies. The photon wavelength equals the spatial period of
the emitted helix.

**Testable**: Initialize a helical wave packet (Gaussian envelope ×
traveling helix) and check: does it propagate at c without dispersing?
Does it carry definite momentum? Does it interact with a braid?

### The m → 0 Problem: Reframing in Terms of δρ

**The standard framing**: The equation has a bare mass parameter m².
Setting m=0 makes gravity long-range but destabilizes the vacuum.

**The alternative framing**: What if "mass" isn't a bare parameter
at all, but emerges from the BACKGROUND FIELD?

Consider the linearized equation around the background φ_a^(0):

    φ_a = φ_a^(0) + δφ_a

Expand the equation of motion to linear order in δφ:

    ∂²δφ_a/∂t² = ∇²δφ_a - m²δφ_a - V''(P_bg) × (...)

The effective mass for perturbations is:

    m_eff² = m² + V''(P_bg) × (background-dependent factor)

With V(P) = (μ/2)P²/(1+κP²) and μ < 0:

    V''(P) = μ(1 - 3κP²)/(1 + κP²)³

At small P (weak background): V'' ≈ μ < 0 (NEGATIVE contribution).
At large P (strong background): V'' → 0 (saturates).

**The current situation (m²=2.25)**:

    m_eff² = 2.25 + μ × (small) ≈ 2.25 - 41 × O(A_bg⁴) ≈ 2.25 - small

The bare m² dominates. Perturbations have mass ~1.5. Yukawa range ~0.67.

**The question: can the background itself provide stability at m²=0?**

If m²=0:

    m_eff² = V''(P_bg) × (...)

With μ < 0, this is NEGATIVE → tachyonic → unstable. This is exactly
what V33-C4 observed.

**But what if μ > 0?** With positive μ:
- V'' > 0 → positive effective mass from background → STABLE at m=0
- But: positive μ means V(P) is a repulsive potential, not attractive
- Braids might not form (they need the attractive V to self-bind)

**Or: what if the background is strong enough that 1-3κP² < 0?**

When κP² > 1/3, V'' flips sign (even with μ < 0). A sufficiently
strong background could provide a POSITIVE m_eff² even with negative μ.
This requires |P_bg| > 1/√(3κ) = 1/√150 ≈ 0.082.

With A_bg=0.1, P_bg ≈ A_bg³ × cos³(...) ~ 0.001. This is much too
small. Would need A_bg ~ 0.4-0.5 for this to work.

**The δρ reframing — is it useful?**

The user asks: "does the equation need to be re-framed in terms of δρ
rather than m?"

Rewrite the equation as perturbations around the background:

    ∂²δφ_a/∂t² = ∇²δφ_a - m_eff²(ρ_bg) × δφ_a - nonlinear(δφ)

Now m_eff depends on the LOCAL background density. In the depleted
zone around a braid, ρ is lower → m_eff is lower → perturbations
propagate further. In the undisturbed background, m_eff is higher →
perturbations are confined.

This gives a POSITION-DEPENDENT effective mass:

    m_eff²(x) = m² + V''(P_bg(x))

This is EXACTLY what we measured in the footprint asymmetry test!
The Yukawa range varies with position because m_eff varies with ρ.

**The key insight**: The "mass" that determines the gravitational
range is NOT the bare parameter m², but the local m_eff² which depends
on the background density. The bare m² sets the MAXIMUM effective mass
(in undisturbed background). The depletion LOWERS the local effective
mass, EXTENDING the gravitational range near the braid.

**Does this resolve the m→0 problem?**

Not directly. You still need m² > 0 (or equivalently, a positive
m_eff² everywhere) to keep the vacuum stable. But it reframes the
question:

- The GRAVITATIONAL range is set by m_eff in the depleted zone,
  which is SMALLER than the bare m
- If the depletion is deep enough, m_eff → 0 locally → the force
  range extends dramatically near the braid
- The VACUUM is stabilized by m_eff in the undisturbed background,
  which equals the bare m

So the system has TWO scales:
1. Vacuum stability: requires m_eff(background) > 0 → m > some threshold
2. Gravitational range: set by m_eff(depleted zone) < m → LONGER than
   the bare Yukawa range

**This might be sufficient**: if the depletion is deep enough, m_eff
near the braid could be much smaller than the bare m, giving
near-1/r² gravity at intermediate ranges while keeping the vacuum
stable at long distances where ρ = ρ_bg.

**The question is quantitative**: How much does the depletion reduce
m_eff? The depletion is δρ/ρ ≈ -15% at the braid surface (V29-T11).
Does this give a measurably longer range?

### Field Quantization and Collapse Prevention

**Speculation level**: High. Worth recording, hard to test directly.

**The idea**: If the field has a minimum quantum ε (a smallest possible
δφ per "cell" of the physical substrate), then:

1. The collapse feedback at m=0 (δρ → stronger V → deeper δρ → ...)
   cannot proceed continuously — it advances in discrete steps
2. The depletion has a FLOOR — you cannot remove more field than exists
3. The Yukawa tail doesn't decay to true zero, it hits the quantum floor
4. m_eff² = m² + V''(P_bg) is bounded below by the quantized background

**Why the simulation grid doesn't already do this**: The grid discretizes
SPACE (limits gradient steepness = UV cutoff on k), but field values are
continuous float64. A physical field quantum discretizes the FIELD VALUES
themselves (limits well depth = IR cutoff on δρ). The m=0 collapse in
V33-C4 may be driven by the well getting arbitrarily deep, not by
gradients getting arbitrarily sharp. A field-value quantum would arrest
this where the spatial grid does not.

**Connection to ds = c·dt**: If c is the propagation speed on a discrete
substrate, the substrate has a fundamental cell size. Each cell carries a
finite field value that changes in quantized steps. This provides:
- Maximum frequency ω_max (Debye-like cutoff)
- Minimum field step ε per cell
- Natural regularization of the m=0 instability

**Connection to κ saturation**: The κ parameter in V(P) = (μ/2)P²/(1+κP²)
already provides a CLASSICAL cap on the interaction energy (V → μ/2κ as
P → ∞). A field quantum would provide an ADDITIONAL cap at a different
scale. If κ saturation alone can prevent m=0 collapse (Track C), then
a physical field quantum would do the same at the substrate scale.

**Practical test**: Track C (enhanced κ at low m) is the classical proxy.
If κ=1000 stabilizes m=0, the mechanism is well-depth saturation — and
a physical field quantum would provide the same effect naturally.

### Electrons as Orbital Modes

**The idea**: An "electron" is not a second braid. It is a STANDING
WAVE MODE of the field around the braid (the "nucleus").

The braid creates a depletion well (δρ < 0 around it). This well
acts as a potential for perturbations of the field. If the well
supports BOUND MODES, those modes are the "electron orbitals."

From the V12 normal mode analysis:
- K=0 (breathing/radial) modes: NO bound states found
- K≥1 (angular/dipole) modes: NOT YET TESTED
- The depletion well has depth ~15% and width ~6 code units

For hydrogen analog (single orbital): need ONE bound K=1 mode in
the depletion well. The mode's frequency would set the "electron
energy." Transitions between modes would emit photons (traveling
helical wave packets).

**Connection to EMF**: The electron (bound orbital mode) couples to
the photon (traveling wave) through the braid's nonlinear V(P) term.
When the orbital mode is excited, it can decay by emitting a
propagating helical perturbation = photon emission.

---

## Experimental Plan

### Track A: Photon Propagation
1. Initialize a traveling helical wave packet (no braid, just a pulse)
2. Does it propagate at c without dispersing?
3. Does it carry definite momentum and energy?
4. Does it interact with a stationary braid? (absorption, scattering)

### Track B: m_eff Mapping
1. Compute m_eff²(x) = m² + V''(P_bg(x)) around a settled braid
2. Map the effective mass profile: m_eff(r) from core to far field
3. Compute the effective Yukawa range 1/m_eff(r) at each radius
4. Does the range extend significantly in the depleted zone?

### Track C: Vacuum Stabilization at Low m

**Critical insight**: The instability at m=0 lives at SMALL P (the vacuum),
not at large P. V''(P=0) = μ = -41.3 regardless of κ. Increasing κ only
caps the energy at large P — it does NOT fix the vacuum instability.

To stabilize the vacuum at m=0, need V''(0) > 0. Options:

C1. **Strong background**: A_bg large enough that P_bg > 1/√(3κ)
    → V'' flips positive. Needs A_bg ≈ 0.4+ (currently 0.1).
    Tests: m²=0, A_bg=0.3, 0.4, 0.5 with κ=50

C2. **Repulsive φ⁴ term**: Add λ(Σφ²)² to the potential.
    This makes m_eff² = V''(P) + 4λΣφ² → positive for large enough λ.
    Tests: m²=0, λ=0.1, 1.0, 10.0

C3. **Positive μ regime**: μ > 0 makes V'' > 0 everywhere. But braids
    need attractive V to self-bind. Test: do braids survive with μ > 0?

C4. **κ sweep** (for completeness, likely won't help):
    m²=0, κ=100,200,500,1000. Predicted: still unstable at all κ
    because the instability is at P≈0, not at large P.

C5. **Field-value quantization**: discretize φ to multiples of ε.
    May prevent tachyonic growth if the minimum step is larger than
    the instability growth rate per timestep. (See Track F)

### Track D: Orbital Modes
1. Linearize the equation around a settled braid's field configuration
2. Solve for K=1 (dipole) angular eigenmodes in the depletion well
3. Are there bound states? What are their frequencies?
4. If bound states exist: these are the "electron orbitals"

### Track E: κ-Sweep Force Law
1. At fixed m²=2.25, vary κ=10,50,100,500
2. Measure force exponent n and range for each κ
3. Does higher κ give longer range (by deepening the depletion)?

### Track G: Metastability Window (HIGHEST PRIORITY)
1. Fine m² scan: 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0
   - Single braid, N=128, T=500 at each m²
   - Measure: braid survival time, E_drift rate, force exponent n
2. dt convergence at m²=1.0: run at dt×0.5, dt×1.0, dt×2.0
   - If drift scales as dt² → numerical artifact → m²=1.0 may be stable
   - If drift is constant → physical instability → find the threshold
3. Two-braid force law at each stable m²
   - Measure n in F ∝ 1/D^n
   - Map: n vs m² → find the sweet spot (n closest to 2.0)
4. Gravitational range: compute 1/m and compare to effective range
   from the depletion (m_eff from Track B)

This is the most practical path to resolving F1 from FUTURE.md.

### Track F: Field-Value Quantization (Speculative)
1. Modify the simulation to enforce a minimum field step: δφ rounds
   to nearest ε (e.g., ε=0.001, 0.01, 0.1)
2. Run m²=0 with field quantization: does it prevent collapse?
3. If stable: measure braid survival, force law, depletion profile
4. Compare to Track C results: is quantization equivalent to high κ,
   or does it produce different physics?

---

## Dependencies

    Track A (photons) → independent, can start immediately
    Track B (m_eff map) → independent, can start immediately
    Track C (κ at low m) → independent, can start immediately
    Track D (orbitals) → needs a settled braid profile from Track B
    Track E (κ-sweep) → independent, can start immediately
    Track F (quantization) → independent, but do C first (cheaper test)
    Track G (metastability) → INDEPENDENT, START FIRST — most practical

    Tracks A + D together → photon emission/absorption test
    Tracks B + C together → understanding whether δρ reframing resolves m→0
    Tracks B + G together → does depletion extend range beyond bare 1/m?
    Tracks C + E together → parameter space for viable long-range gravity

---

## Combined Experiments

### GB: Metastability + Effective Mass (G+B)

At each stable m² found by Track G, compute the full m_eff(r) profile
around the braid. This answers: **does lowering the bare m² amplify
the depletion's range-extending effect?**

At m²=2.25 (current), m_eff ≈ 2.25 - small_correction → range barely
extends. But at m²=1.0, the bare mass is closer to the V'' correction,
so the fractional reduction could be MUCH larger:

    m²=2.25: m_eff² ≈ 2.25 - 0.1 = 2.15  (5% reduction, negligible)
    m²=1.00: m_eff² ≈ 1.00 - 0.1 = 0.90  (10% reduction, modest)
    m²=0.50: m_eff² ≈ 0.50 - 0.1 = 0.40  (20% reduction, significant)

The RATIO of range extension grows as m² approaches the V'' correction.
Near the braid core where δρ is deepest, m_eff could drop to near zero
even with a nonzero bare m².

**Method**:
1. Take settled braid snapshots from G Phase 1 at each stable m²
2. Compute m_eff²(x) = m² + V''(P(x)) at every grid point
3. Average in cylindrical shells → m_eff(r) profile
4. Compare: bare range 1/m vs minimum m_eff range 1/m_eff_min
5. Plot: m_eff(r) for each m² on same axes

**Key prediction**: At lower m², the depletion's "range tunnel" (where
m_eff << m) extends further from the braid. The effective gravitational
range may be 2-5× the bare Yukawa length 1/m at marginal m² values.

This would mean gravity has TWO regimes:
- Near the braid: m_eff ≈ 0 → ~1/r² (no exponential suppression)
- Far field: m_eff → bare m → Yukawa cutoff
The crossover distance depends on the depletion depth.

### GBF: Metastability + Effective Mass + Field Quantization (G+B+F)

The most speculative but potentially most interesting combination.

**The hypothesis**: At marginal m² (from Track G), the m_eff profile
(from GB) shows m_eff → 0 near the braid core. Adding field-value
quantization (Track F) prevents the tachyonic instability that would
otherwise develop where m_eff < 0. The quantization acts as a "floor"
on m_eff, allowing the bare m² to be reduced BELOW the classical
stability threshold.

**Method**:
1. Find the lowest classically stable m² from G Phase 1
2. Compute m_eff profile from GB — identify where m_eff² first goes
   negative (this is the classical instability boundary)
3. Add field quantization: after each Verlet step, round each
   φ_a to the nearest multiple of ε
4. Run at m² values BELOW the classical threshold (into the
   classically unstable regime)
5. Does the quantization stabilize the vacuum?
6. If YES: measure the force law — is it closer to 1/r²?

**Three regimes to test**:

    ε = 0.001 (fine quantum — nearly continuous, should behave classically)
    ε = 0.01  (moderate — comparable to background amplitude fluctuations)
    ε = 0.1   (coarse — comparable to A_bg itself)

**What we learn**:
- If ε=0.01 stabilizes a previously unstable m²: the field quantum
  IS the missing stabilization mechanism, and m² can be lowered further
  than classical analysis suggests
- If NO ε stabilizes it: the instability is not about field-value
  depth but about spatial spreading (mode count), and quantization
  doesn't help
- If ε=0.1 stabilizes but kills the braid: the quantum is too coarse
  to support structure — need a finer substrate

**Why this matters**: If GBF works, it suggests the physical field
quantum (from ds=c·dt) IS what prevents the m=0 Big Bang from
destroying the universe permanently. The universe STARTS with a
Bang (m_eff ≈ 0 everywhere), the expansion dilutes the field,
and the field quantum prevents re-collapse. Braids form when the
background cools enough that m_eff > 0 in the depleted zones.
The bare m² could be ZERO, with all stability from the quantum.

---

## Key Unknowns

1. Can a traveling helical wave packet propagate without dispersing?
   If NO: photons need a different mechanism.
   If YES: EMF is naturally the traveling mode of the same field.

2. How much does the depletion reduce m_eff?
   If negligibly: the δρ reframing doesn't help.
   If significantly: the gravitational range extends naturally near braids.

3. Does the depletion well support K=1 bound modes?
   If NO: electrons need a different mechanism (complex fields? gauge?).
   If YES: atomic structure emerges from the standard equation.

4. Can κ alone stabilize m=0?
   If NO: m² is truly fundamental, Yukawa is unavoidable.
   If YES: long-range gravity is achievable with parameter tuning.

5. How much does lowering m² amplify the depletion range extension? (GB)
   If the ratio m_eff_min/m shrinks significantly at lower m²:
   gravity range grows FASTER than 1/m — depletion and bare mass
   work together synergistically.

6. Can field-value quantization stabilize below the classical threshold? (GBF)
   If YES: the physical field quantum is the missing piece, and
   m² = 0 may be achievable with quantization alone.
   If NO: the instability is spatial (mode proliferation), not
   depth-limited, and a different stabilization is needed.

---

## Resolved Experiments

### Track G: Metastability Window (COMPLETE)
- **Stability threshold**: m² ≥ 0.50 for vacuum, m² ≥ 1.25 for braids
- **Braid survival is the binding constraint**, not vacuum stability
- **dt convergence**: drift at m²=0.50 is physical (constant across dt)
- **Best constant m²**: 1.50 (97% E_pot retention, range 0.82 = 1.83× current)
- Results: `G_metastability/RESULTS.md`

### Track GB: Field-Dependent Mass (COMPLETE — NEGATIVE)
Two approaches tested, both fail to extend range:

**φ⁴ mode (m_eff² = α×Σφ²)**: Mass too strong at core → overwhelms
V(P) binding → braid dissolves. E_pot → 0 at ALL α values.

**Inverse mode (m_eff² = α/(1+β×Σφ²))**: Two regimes:
- m²(bg)=2.25 (β=1-5): Braid survives, no range improvement
- m²(bg)≤1.0 (β=5-10): V(P) overwhelms reduced core mass → κ-saturated
  collapse → blob slowly drains. Not viable.

**Root cause**: Constant m² is a GOLDILOCKS — it simultaneously provides
vacuum stability and braid binding at the same value. Field-dependent
mass via any monotonic function of Σφ² cannot decouple these because
core Σφ² is only ~60× background Σφ², insufficient dynamic range.
Results: `GB_field_mass/RESULTS.md`

### Phonon Test (COMPLETE — MAJOR POSITIVE RESULT)
The braid's depletion profile is POWER LAW (δρ ∝ 1/r^1.2), NOT Yukawa.
Yukawa m=1.5 is excluded at 500,000× discrepancy. The long-range force
is carried by a massless collective mode (phonon) of the background, not
by massive φ_a excitations. m² does NOT limit gravitational range.

Results: `phonon_test/RESULTS.md`

**This changes the picture fundamentally**: The m→0 problem is NOT a
problem. The gravitational range was never m²-limited. Tracks G and GB
were solving the wrong problem — the range is already infinite (power law).
The remaining question is the EXPONENT: n≈1.2 instead of 2.0.

### Track A Phase 1a: Mode Analysis (COMPLETE)
All three linearized modes around the uniform background are massive
(m_eff ≈ 1.5, splitting only Δm² = 0.004 at A_bg=0.1). No Goldstone
mode exists. The power-law depletion is a NONLINEAR collective effect,
not a single linear mode. Mode splitting grows near the braid surface
(amplitude mode becomes nearly massless at A≈0.5), but all modes are
massive in the far field. Results: `A_emf/RESULTS_phase1a.md`

### Track A Phase 2: Torsion Waves as EM (COMPLETE — NEGATIVE)
Four experiments testing torsion (field-space rotation) as EM carrier:
1. Torsion profile differs from depletion (peaks at shell, not core)
2. Torsion waves do NOT propagate — disperse in ~5t, convert to amplitude
3. Braid response: torsion and amplitude are qualitatively the same
4. Opposite winding braids: NO opposite response — winding ≠ charge

Root cause: all far-field modes are massive. No massless carrier exists
for a second long-range force in the current real-field theory.
Results: `A_emf/RESULTS_phase2.md`

### Cosserat 6-Field Extension (COMPLETE — MAJOR POSITIVE RESULT)
Added 3 angle fields θ_a to the 3 position fields φ_a, coupled through curl.
The equation becomes (Eq. 10):
    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - V'(P) + η×curl(θ)_a
    ∂²θ_a/∂t² = ∇²θ_a - m_θ²θ_a       + η×curl(φ)_a

Results:
- η scan: braid survives at η=0.1-1.0, dissolves at η=2.0, explodes at η=5.0
- θ_rms ≈ 0.065 — braid helical curl sources the angle fields
- **m_θ²=0 (MASSLESS) is STABLE** — no tachyonic instability (θ has no V(P))
- Massless θ helps the braid survive 3× better (provides radiation channel)

Code: `torsion_coupling/src/v33_cosserat.c`
Data: `torsion_coupling/data/cosserat_eta_*/`, `cosserat_mt*/`
Results: `torsion_coupling/RESULTS_cosserat.md`

### θ Field Characterization (COMPLETE)

Three experiments + dt convergence test:

1. **Radial decay**: θ_φ² ∝ r^(-0.5), NOT r^(-2). The field is an
   oscillating wave (period ~4t), not a static 1/r circulation.
   DC component is 0.2% of oscillation amplitude.

2. **Winding reversal**: Time-averaged DC θ_φ reverses sign between
   W=+1 and W=-1. Far-field ratio ≈ -1.0. Correlation: -0.68.
   CONFIRMED: winding = charge.

3. **Two-braid force**:
   - 3-field baseline (gravity): ΔD = -5.05
   - Same winding (gravity + θ): ΔD = -6.41 (27% more attraction)
   - Opposite winding (gravity + θ): ΔD = -2.16 (57% less attraction)
   CONFIRMED: charge-dependent force with correct sign.

4. **dt convergence**: θ_rms varies by 2% across dt×{0.25, 0.5, 1.0}.
   CONFIRMED: θ structure is physical, not numerical artifact.

Results: `torsion_coupling/theta_characterize/RESULTS.md`

### Remaining open questions:
1. Static Coulomb regime: is there a θ force between non-moving braids?
2. Fine structure constant analog: what sets the θ/φ force ratio?
3. Electron orbitals: do the θ shells support bound standing-wave modes?
4. Depletion exponent: why n≈1.2? Isotropic background + longer runs needed.
5. Visualization: the SFA archive format + volume viewer are built and working.
   The right-hand-rule θ pattern is visually confirmed.
