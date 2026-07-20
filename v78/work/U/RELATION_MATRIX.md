# v78 Relation Matrix R1–R10 (U freeze / W gap-fill)

**Status:** FROZEN for CP-W  
**Sources:** `GOALS.md` §2; v76 F1-3D; v77 WORLD/Maxwell/RC1; v69–v75 kernel  
**Role:** who couples to whom — monist language ∥ SCP particle engine

---

## Primitive stack

| Layer | Source | Particle role |
|-------|--------|---------------|
| Complex Cosserat Φ_a, Θ_a | kernel-v3 | Matter continuum / bags |
| Gauged diagonal U(1) A,E | v69–v71 | EM / Coulomb / Gauss |
| Conserved Q, Q_a, ω | U(1) Noether | Inventory, flavor |
| Free capacity ψ | v76 monist | Path-cost / mass-form geometry (**parallel**) |
| Free gauge Maxwell | v77 | Sibling free-gauge channel (**parallel**) |
| Multi-fabric C/Q/L | v75 | Bag / nuclear EM / light opposite |

**Primitive claim:** particles are **stable localized field configurations** (Q-balls / multi-fabric composites), not point masses on empty space.

---

## Matrix

| ID | Relation | Fabric expression | Particle consequence | Owner layer | Status |
|----|----------|-------------------|----------------------|-------------|--------|
| **R1** | Locality \(c\) | Free update bound; monist \(c=1/\sqrt{\varepsilon\mu}\); kernel wave speed | Shared causal structure | monist + kernel | LIVE |
| **R2** | Mass ↔ energy | \(m=E_\star/c^2\); soliton \(E=\int T_{00}\) | Binding / mass defect (E/Q vs free N) | monist + kernel | LIVE |
| **R3** | Charge ↔ Gauss | \(\nabla\cdot E \propto \rho_{\mathrm{em}}\); discrete Gauss floor ~1e-13 | p EM charge; n ~0 EM | kernel (v75 q weights) | LIVE |
| **R4** | Path-cost / gravity-like | ψ from ρ_b (v76 F1-3D) | Large-scale geometry; **optional** at nuclear T | monist sandbox | LIVE (parallel) |
| **R5** | Short-range binding | \(V_t(s)\), phase coherence, private bags | Nuclear fusion / liquid-drop; C–L isolation when multi-fab | kernel | LIVE |
| **R6** | Coulomb | Shared A; opposite q | Atom binding; P–P repulsion; ± force | kernel | LIVE (F8/F11) |
| **R7** | Flavor / multiplet | \((ω_a,Q_a)\) | p/n **not** from flavor alone; B2 Q presence is EM knob | kernel | LIVE (F17 kill flavor-n) |
| **R8** | Multi-fabric isolation | Private \(s_f\); no C–L product mix | No core–cloud bag merge | kernel multi-fab | LIVE (F11 G3) |
| **R9** | Budget free/bound | Monist \(\rho_f+\rho_b\) | Mass-formation language for soliton energy | monist | LIVE (parallel) |
| **R10** | Stability | VK branch dQ/dω<0; park c_Q; long-T drift | Proton/neutron/C₁₂ survival gates | kernel + analysis | LIVE |

---

## Coupling graph (schematic)

```text
                    ┌── R4 ψ path-cost (optional nuclear) ──┐
                    │                                      ▼
  Φ bags ──R5──► fusion / park          free continuum ψ / Maxwell
    │                                      ▲        ▲
    │ R7 flavor                            │ R1 c   │ R9 budget
    ▼                                      │        │
  Q, Q_a ──R3──► ρ_em ──R6──► A,E ─────────┴────────┘
    │                         │
    │ R8 private bags         │ R6 Coulomb (atoms)
    ▼                         ▼
  C / Q / L fabrics         light opposite cloud
    │
    └── R10 park / VK / Gauss floor
```

---

## Explicit non-identifications

| Forbidden | Why |
|-----------|-----|
| ψ ≡ Φ_EM | Dual-channel: path-cost ≠ electrostatic (v77) |
| Flavor multiplet alone = neutron | Residual EM (v75 F17) |
| One same-sign fabric = atom | v74: nuclei only (F1) |
| Local free-density GRIN as gravity | v76 dead branch |

---

## Particle map of relations

| Target | Dominant relations |
|--------|-------------------|
| Nucleon | R2, R5, R7, R10 (+ R3 if gauged) |
| Proton analog | R3, R5, R8 (C+Q) |
| Neutron analog | R5, R8 (C only); R3 null EM |
| Z-carbon nucleus | R5 fusion, R2 mass defect, R10 park |
| Light L / e-analog | R6, R8, R10 |
| Atom | R6 bind + R8 no-merge + R3 Z-count + R10 |

---

## Freeze stamp

**U ADOPT CP-W** — matrix complete for campaign scoring.  
Peer W may publish root `PHYSICS_RELATIONS.md` isomorphic to this freeze.
