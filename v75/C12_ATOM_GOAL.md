# Goal: Multi-fabric path to time-stable C₁₂ atom

**Set:** 2026-07-15  
**Scope:** Phases **1–3** of the multi-fabric roadmap (hydrogenoid → parked multi-Z + L → A≈12-class core + light cloud).  
**Architecture:** B1 multi-fabric (private bags, shared \(A\), U^q); park-aware scoring for multi-ball nuclei. B2 / Option C only if blocked.  
**Baseline:** F11–F16 (isolation, Coulomb, B4 packaging, pair kinematics, Z6+L6 PASS_park).

---

## Ideal end state (primary goal)

**Demonstrate a time-stable multi-particle C₁₂ analog** on the full gauged fabric:

| Property | Requirement |
|----------|-------------|
| Nuclear core | Multi-ball / droplet at **carbon-scale** (target A≈12 or documented Z-carbon + path to A=12) |
| **P/N distinction** | **Firm proton- vs neutron-analog species** in the nuclear sector (see §P/N below) |
| Light sector | Multi-L opposite-charge cloud (private bag); cloud charge tracks **Z** (proton count), not N |
| Structure | **Bound** atomic composite: core + cloud remain distinct sectors |
| Stability | Long-T evolution: **no bag merge**, **no fission**, **no dispersion** of cloud or core |
| Fields | Retain related field content (matter + shared \(A\); Gauss floor) |
| Visual | Headless `volview` (or equivalent) snapshots/movie: nuclear charge density + light sector + \|E\| / charge view over time |

**Pass bar (ideal):** T ≳ several nuclear / orbital timescales (multi-rev if orbiting; or long rest with no secular loss), massL and core Q within park-aware bounds, clear two-sector morphology in renders; **Z** and **N** both reportable and independently adjustable in seeds.

---

## P/N distinction (hard requirement for C₁₂ + isotopes)

Real carbon chemistry: **atomic number Z** = proton count fixes the electron cloud; **neutron number N** changes mass/stability (**isotopes**) without changing Z.  
To build **stable C₁₂** and later an **unstable variant by adding neutrons**, the fabric must support:

| Symbol | Meaning in sim | Must control |
|--------|----------------|--------------|
| **Z** | Count / charge of **proton-analogs** in the core | EM charge of nucleus (and thus L cloud matching) |
| **N** | Count of **neutron-analogs** in the core | Nuclear mass / binding / stability **without** changing net nuclear EM charge |
| **A** | Z + N | Carbon-scale ≈ 12 for the ideal package |

**Firm distinction means:**

1. **Two nuclear species** (p-analog and n-analog) that are **dynamically stable** as separate lump types (or as two branches of the flavored multiplet), not merely two labels on the same ball.  
2. **Net EM charge of the nucleus ∝ Z only** (n-analog contributes **≈ 0** to gauged charge / Gauss source).  
3. **Adding n-analogs** can move the core toward metastability (super-criticality, deformation, evaporation) **while Z (and L cloud charge) stay fixed** — the isotope control knob.  
4. Diagnostics report **Z, N, A** (or proxies: charged vs neutral nuclear cluster charges / flavor partitions).

**Theory hook (existing):** CONCEPT flavor multiplet — frequency partition \((\omega_a,Q_a)\) as p/n-analog distinction.  
**Gap today:** v75 atom runs use generic same-sign heavy balls; **net-neutral neutron-analog + charged proton-analog under shared \(A\)** is not yet a frozen package. Closing this gap is **on the critical path before stretch decay-by-N**.

**Candidate implementation paths (choose by experiment, not yet locked):**

| Path | Idea | Risk |
|------|------|------|
| **Flavored multiplet** | p-branch vs n-branch of \(\omega_a\) partition; n-branch engineered for \(Q_\mathrm{em}\approx 0\) | Same-sign gauged charge on components — neutrality nontrivial |
| **Composite n** | Neutral bound of opposite internal charges or cancelling multipoles | May need B2 / extra structure |
| **Two nuclear charge species under A** | Explicit \(q=+1\) vs \(q=0\) nuclear matter (beyond B1 lock) | Kernel/config extension; may need B2 |

Do **not** fake isotopes by only changing total Q or ball count without a Z-preserving N knob.

---

## Stretch goal: radioactive C₁₂ variant (isotope path preferred)

Primary stretch path: **same Z, larger N** (neutron-rich) metastable package — not only super-Q same-composition seeds.

Produce a package that:

1. Shares **Z** (and L cloud charge) with the stable C₁₂ analog.  
2. Differs by **N** (extra neutron-analogs) and/or controlled deformation / super-criticality.  
3. **Calculates** a decay rate or lifetime estimate (evaporation, fission barrier proxy, radiation power, mean escape time).  
4. **Simulates** the decay channel with time series.  
5. Contrasts with the **stable** (Z,N) package under the same multi-fabric laws and scorecard.

Fallback stretch (if P/N not yet firm): super-Q or deformed seed at fixed composition — documented as **not** a true isotope until P/N exists.

---

## Phase plan (1–3 only for this goal)

### Phase 1 — Bound hydrogenoid (single-center atom)

| ID | Task | Success |
|----|------|---------|
| P1.1 | Retune circular vt (~0.05–0.06) from F16 bracket | Near-flat D(t) vs sub/super |
| P1.2 | Multi-rev orbit T ≳ 2000–4000 on B4 θ | No secular L loss; Gauss floor |
| P1.3 | Multi-cluster / shell-radius diagnostic for multi-L | Shell not only COM D≈0 |
| P1.4 | Binding proxy (E vs free C+L) | Documented binding signal |

**Phase 1 done when:** multi-rev (or clear bound) single-C + L cloud, visualizable.

### Phase 2 — Parked multi-Z nucleus + light cloud (+ P/N foundation)

| ID | Task | Success |
|----|------|---------|
| P2.0 | **Firm p- vs n-analog definition** (seed + diagnostics for Z, N) | Two species; nuclear \(Q_\mathrm{em}\propto Z\); n adds A not Z |
| P2.1 | Standard parked templates (c6 / flavored p+n mixes) | Qc ≈ parked band; report Z, N |
| P2.2 | L shell co-design from **Z** (not A) | PASS_park + massL stable; L tracks Z |
| P2.3 | Soft kinematics around **droplet** (not point C) | Ordered D_shell; no L absorb |
| P2.4 | Long-T rest stability + volview package | No merge / fission / disperse |
| P2.5 | Isotope smoke: fixed Z, +ΔN vs baseline | Core stability changes; L charge unchanged |

**Phase 2 done when:** Z-carbon-class core + L cloud is time-stable and visual, **and** Z/N are independently reportable.

### Phase 3 — A≈12-class core + light cloud → ideal C₁₂

| ID | Task | Success |
|----|------|---------|
| P3.1 | Q_max / g survey for A=12 (lower g or sub-critical recipe) | Parked or quasi-stable A=12-class core |
| P3.2 | Born–Oppenheimer: freeze nucleus, relax L; optional back-react | L cloud in nuclear multipole field |
| P3.3 | Assemble **(Z,N) with Z≈6, A≈12** + multi-L; park-aware score | PASS_park; no L death; Z fixed under N variation tests |
| P3.4 | Long-T stability + **visual C₁₂ package** (primary ideal) | No merge / fission / disperse |
| P3.5 | Stretch: **neutron-rich (same Z, higher N)** and/or deformed → decay calc + sim | Rate estimate + observed channel; control = stable (Z,N) |

**Phase 3 done (ideal) when:** time-stable visual C₁₂ composite with definite (Z,N).  
**Phase 3 stretch when:** radioactive **isotope-class** variant (prefer +N at fixed Z) + calculated and simulated decay.

---

## Scorecard (mandatory)

- **Multi-ball nuclear:** park-aware \(c_Q\) (mid→end), not seed→end alone.  
- **Hard fails:** bag merge (L into C), core fission (unintended), cloud dispersion (massL collapse), Gauss drift above floor.  
- **Soft fails:** large E drift without morphology loss (document, may retune).

---

## Explicit non-goals (this goal)

- Full B2 ε_CQ unlock unless Phase 2–3 blocked by Q←C lock  
- Option C triple Cosserat unless engagement graph is the wall  
- Cosmology / gravity  
- L3 n-bar chemistry populations before L0 C₁₂ stability  

---

## Tracking

| Doc | Role |
|-----|------|
| This file | Goal definition + phase checklist |
| `STATUS.md` | Live progress |
| `FINDINGS.md` | Measured results (F17+) |
| `FUTURE.md` | Pointer to this goal |

**Start state (2026-07-16 checkpoint):** Phase 1 open. **P2.0 DONE (F17).** **Z2 package DONE (F18).** **Z6 isotope EM PASS (F19)**; Z6N6 park PASS; Z6N0 soft; L −12.5%; multi-ball **nuclear cores merge to a single droplet** by t∼400. Phase 3 not started.

**Canonical physics snapshot:** `STATE.md` (equations, charge table, P/N freeze).

**Critical path note:** Bottleneck is **droplet park quality + L hold + shell diagnostics**, then P3 long-T visual C₁₂ — not inventing Z/N.
