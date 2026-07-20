# Gap analysis — approaches A / B / C

**Date:** 2026-07-20  
**Purpose:** Theory completeness + analytic reading of numerics. Where each path is incomplete or wrong.

---

## Approach A — capacity well

### What was claimed

\[
F_{\mathrm{net}}(D)=F_{\mathrm{EM}}(D)-k\,e^{-D^2/(2s^2)},
\quad
F_{\mathrm{EM}}\approx\frac{1}{2\pi D}\ (2\mathrm{D}).
\]

Label “well” = first zero where \(F_{\mathrm{net}}\) crosses **repel → attract** (\(r_*\sim7\)–\(9\)).  
Seed circular orbit with \(v_t=\sqrt{F_{\mathrm{EM}}\,r_*/\mu}\) at that \(r_*\).

### Numeric facts (from `A_Fnet.tsv` / log)

| Check | Result |
|-------|--------|
| \(F_{\mathrm{EM}}\) vs \(1/(2\pi D)\) | Ratio \(\sim0.99\)–\(1.07\) for \(D\in[3,12]\) — medium Coulomb OK |
| Zero of \(F_{\mathrm{net}}\) | Exists for \(k\gtrsim0.1\), \(s=4\) |
| Free orbit at claimed \(r_*\) | \(\mathrm{sepmin}=9\), \(\mathrm{sepmax}\sim49\), \(\mathrm{revs}=0.22\) — leaves “well” |

### GAP A1 — **Wrong circular condition (theory bug)**

Relative 1D radial problem with conserved angular momentum \(L=\mu\,r\,v_\theta\):

\[
\mu\ddot r = F_r(r)+\frac{L^2}{\mu r^3}.
\]

Here \(F_r=-F_{\mathrm{along}}\) if \(F_{\mathrm{along}}>0\) means attract (force that decreases \(r\)).  
**Circular** \(\ddot r=0\), \(\dot r=0\):

\[
F_{\mathrm{along}}(r_c)=\frac{\mu v_\theta^2}{r_c}\;>\;0.
\]

So the pair must sit where **net force is still attractive**, equal to the centripetal requirement — **not** at \(F_{\mathrm{net}}=0\).

At \(F_{\mathrm{net}}=0\) there is **no** force left for \(\mu v^2/r\). The code then seeded

\[
v_t=\sqrt{F_{\mathrm{EM}}\,r/\mu}
\]

using **bare EM**, not \(F_{\mathrm{net}}\). At \(r_*\) where \(F_{\mathrm{net}}=0\), true \(F_{\mathrm{net}}\ll F_{\mathrm{EM}}\), so \(v_t\) is **systematically too large** → super-circular → climbs out of any shallow well → expand. That matches sepmax~49.

**Correct program:**

1. Choose \(L\) (or \(v_\theta\)).  
2. Solve \(F_{\mathrm{along}}(r)=\mu v_\theta^2/r\) for \(r_c(L)\).  
3. Require **stability**: \(\partial_r\!\left(F_{\mathrm{along}}-\mu v_\theta^2/r\right)<0\) in the appropriate sign convention (effective potential minimum).  
4. Seed at that \((r_c,v_\theta)\), not at \(F_{\mathrm{net}}=0\).

### GAP A2 — **Is there a stable circular family at all?**

Define \(V(r)\) by \(\mathrm{d}V/\mathrm{d}r=F_{\mathrm{along}}\) with \(V(\infty)=0\) (work to separate).  
Effective potential \(V_{\mathrm{eff}}=V+L^2/(2\mu r^2)\).

Analytic sketch for \(k=0.2\), \(s=4\), \(\mu=4\):

- \(V(r)\) has a **minimum** near the force zero (\(r\sim8.7\)).  
- Circular points with \(F_{\mathrm{along}}>0\) exist only for \(r>r_{\mathrm{zero}}\) (attractive side).  
- For large \(L\), the centrifugal barrier may destroy the bound pocket; for small \(L\), circular radius sits just outside \(r_{\mathrm{zero}}\) with small \(F_{\mathrm{along}}\).

So a bound circular family **can** exist, but:

- It lives in a **thin shell** outside \(r_{\mathrm{zero}}\), with **small** \(v_\theta\), not the EM-only \(v_t\sim0.2\).  
- Numerics never seeded that family → orbit FAIL does **not** yet kill the well idea.

### GAP A3 — **Dynamics incompleteness (numeric method)**

Orbit step used: deposit ρ → **few SOR iterations** → force → non-rel update. That is **not** Gauss-preserving Maxwell+current dynamics. Radiation, retardation, and self-force are wrong or missing. Escape may mix:

- wrong \(v_t\) (A1), and  
- numerical heating / force lag (A3).

Until A1 is fixed under a **faithful** stepper, “orbit FAIL” is under-determined.

### GAP A4 — **Ontology vs pair proxy**

Pair Gaussian overlap is a **proxy** for depletion \(n_{\mathrm{free}}\), not the co-field. It reintroduces a **pairwise** form (allowed as regularisation of overlap, but not the monist end state). True capacity must remain **gathered from a medium field** so Gauss/ledger stay clean in 3D.

### A — status after gaps

| Piece | Status |
|-------|--------|
| EM medium force | Complete enough (numeric≈continuum) |
| Short-range repel shape | Complete enough for a well |
| Circular / stability theory in docs | **Wrong / incomplete (A1)** |
| Orbit numeric | **Invalid as falsifier until A1+A3 fixed** |
| Co-field capacity | Incomplete (pair proxy) |

**Not chasing the wrong object** — chasing with the **wrong equilibrium formula** and a **weak dynamical probe**.

---

## Approach B — magnetic / angular

### What was claimed

Uniform \(B_z\) supplies a cyclotron scale \(\omega_c=|q|B/m\); with soft core + Coulomb, multi-rev might park.

### Numeric facts

- \(B_z=0\): already “revs” up to ~1.9 with sep in \([3,31]\) — relative **angle** winds while separation **excurses** (not a radial park).  
- Larger \(B_z\): often **fewer** revs, still huge sepmax, sometimes min sep → soft-core floor.  
- No run with \(\mathrm{sepmax}<2.5 D_0\) **and** \(\mathrm{revs}\ge0.5\).

### GAP B1 — **No radial potential from uniform \(B\) alone**

For a single charge, \(|v_\perp|\) is fixed by energy; radius \(\rho=mv_\perp/(|q|B)\).  
For **two** opposite charges with Coulomb:

- Relative motion in a **uniform** \(B\) does **not** add a scalar \(V(r)\) that balances \(1/r\) (or \(1/r^2\)) into an isolated minimum independent of velocity.  
- Magnetic force is \(\propto v\), odd under \(v\to-v\); it does not define a static radial well.  
- Quantum Landau+Coulomb is a different story (ℏ, Landau levels); we are classical.

So B **cannot** by itself provide the missing **length** the way capacity can. Numerics failing is **expected**, not a surprise bug.

### GAP B2 — **Wrong success metric (angle vs radial bound)**

Counting joining-line revs while \(r(t)\) spans half the box confuses **orientation wind** with **orbital park**.  
A valid B-test needs e.g. time-averaged \(r\), variance of \(r\), or Poincaré section — not revs alone.

### GAP B3 — **Self-consistent \(B\) never tested**

Only **external** \(B_z\). Live medium \(B\) from currents (pair’s own motion) is a different regime (possible support **after** a radial well exists). That gap is optional research, not a fix for B1.

### B — status after gaps

| Piece | Status |
|-------|--------|
| Cyclotron scale for one particle | OK |
| Two-body Coulomb + uniform \(B\) classical park | **Theoretically incomplete as primary scale** |
| Numeric | Consistent with theory; metric muddy (B2) |

**Verdict:** B is the **wrong primary**. Useful later as **angular support / selection**, not as the second scale.

---

## Approach C — composite / rigid dipole

### What was claimed

One sequestration unit with fixed internal \(D_{\mathrm{int}}\); no free relative \(r\) → no pass-through; medium couples to multipole; durable + live \(E_{\mathrm{em}}\).

### Numeric facts

- Rigid \(D_{\mathrm{int}}=6\), free COM + \(\theta,\omega\).  
- Final \(E_{\mathrm{em}}\sim0.34\neq0\), \(\omega\sim0.08\) persists → **PASS\***.  
- COM drifted \((48,48)\to(50,41)\) — **net force on a neutral pair** should vanish in infinite free space.

### GAP C1 — **PASS\* is almost by construction**

“Durable” is **enforced**: internal \(r\) is not a DOF. That **defines away** the Stage-3 problem (free opposite charges that neither evaporate nor pass through while bound by the medium).  
It is a valid **structural** object, but it does **not** demonstrate medium-generated binding of two free carriers.

### GAP C2 — **Self-force / image force on COM**

Neutral \(\pm\) pair: total force \(\sum q_i E_i\) should be \(O(\)discretisation\()\). Observed COM drift ⇒

- CIC self-force asymmetry, and/or  
- periodic images, and/or  
- incomplete Poisson convergence under motion.

So “medium couples cleanly to multipole” is **not** cleanly shown; part of the motion is **lattice artifact**.

### GAP C3 — **No torque/energy theory of free dipole**

A free electrostatic dipole in 2D:

- Has **self-energy** depending on \(D_{\mathrm{int}}\) and shape.  
- Rigid constraint must be maintained by **constraint forces** (Lagrange multipliers / stiff springs) that were **not** put in the medium ledger.  
- Spin–field coupling (radiation of rotating dipole \(\sim \ddot p\)) not measured; “harmonics” undefined without a radiated power vs \(\omega\) curve.

### GAP C4 — **Stage-3 goal mismatch**

Standing goal: light **opposite-charge** content that is medium-bound and eventually combinable with a nuclear bag.  
C gives opposite charge **labels on one rigid body**. Next gaps:

- How does a free \(e^+\) + free \(e^-\) **form** the composite (capture)?  
- How does composite **ionize** / share charge with a nucleus?  
Without a formation channel, C is a **bookkeeping atom**, not a process-form particle.

### C — status after gaps

| Piece | Status |
|-------|--------|
| Durability | Trivial if rigid |
| Live \(E_{\mathrm{em}}\) | True but weak criterion |
| Force integrity (COM) | **Gap (C2)** |
| Formation / ionization | **Missing theory** |
| Relation to free multi-rev | **Does not solve A** |

**Verdict:** C is a **useful degenerate limit** (internal scale fixed), not a complete Stage-3 theory.

---

## Cross-cutting gaps

### G1 — Success metrics mixed three different objects

| Metric | What it tests |
|--------|----------------|
| \(F(D)\) zero / well shape | Static force law |
| Multi-rev + sep band | Free two-body dynamics |
| Durable + \(E_{\mathrm{em}}\neq0\) | Object existence / medium live |

A “PASS” on one does not transfer. A well PASS + orbit FAIL is coherent once A1 is seen.

### G2 — 2D Coulomb is logarithmically special

2D potential \(\sim\log r\) (force \(\sim1/r\)):

- Bound states and virial relations differ from 3D \(1/r\).  
- Radiation and “escape to infinity” are softer/harder depending on IR.  
Orbit failure in 2D sandbox may **not** map 1–1 to 3D kernel. **Gap:** no 3D reduced test of the **corrected** circular condition.

### G3 — Harmonics remain secondary — and untested

No measurement of:

- radiated power vs \(v_t\),  
- multipole content of \(E\),  
- survivor fraction over an ensemble of seeds.

So “orbital harmonics select quantum-like states” is still **philosophy**, not a closed loop on data.

### G4 — Integrator / Gauss / force consistency

Cheap Poisson re-solve ≠ Maxwell + Esirkepov. Capacity in kernel 3D still deferred. Any dynamical claim is provisional until one consistent stepper is used for A (and B).

---

## Where the real GAP is (one sentence each)

| Path | The gap |
|------|---------|
| **A** | Equated **force zero** with **circular orbit**; seeded **too much** \(v_t\); dynamics stepper too crude to falsify the well. |
| **B** | Uniform \(B\) **cannot** supply a classical radial well for Coulomb; metrics counted angle wind as binding. |
| **C** | **Removed** the hard problem (free \(r\)); COM drift shows force/self-force not yet clean; no formation theory. |

---

## What is *not* the gap

- Medium Coulomb: numerically healthy.  
- Need for a short-range **repulsive** scale: still right (A shape works when defined correctly).  
- Multi-fab L / multiplet light: still closed.  
- “Must invent full QM”: not required for the next falsifier.

---

## Tight next steps (if continuing)

1. **A repair (highest leverage):**  
   - Document circular condition \(F_{\mathrm{along}}(r_c)=\mu v^2/r_c\).  
   - Scan \(L\) (or \(v_\theta\)) → \(r_c(L)\); seed those; measure \(r(t)\) variance.  
   - Prefer full PIC + core force (or frozen Poisson only if forces frozen too for debug).  

2. **B demote:** one analytic note “uniform B is not \(V(r)\)”; optional later as \(B_z\) support on top of A.  

3. **C tighten:**  
   - Measure COM acceleration vs lattice resolution (self-force budget).  
   - Write formation/ionization criteria or mark C as **structural only**.  

4. **Harmonics:** only after a seed stays in a radial band for \(\gtrsim1\) rev — then plot radiation vs \(v_t\).

---

## Bottom line

The three options look incomplete because they **are**:

- **A** is the most complete *idea* but was **tested with an incorrect circular-orbit theory** and a weak dynamical scheme.  
- **B** is **theoretically incomplete** as a primary scale.  
- **C** is **definitionally durable** and under-specified as process-form physics.

The highest-value gap to close first is **A1+A3** (right \(r_c,v_\theta\) + honest stepper). Until then, free multi-rev “FAIL” does not kill capacity; and C’s “PASS\*” does not replace it.
