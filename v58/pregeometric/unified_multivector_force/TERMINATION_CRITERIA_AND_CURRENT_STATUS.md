# Termination Criteria and Current Status

**Experiment**: Unified Multivector Force Law  
**Location**: `v58/pregeometric/unified_multivector_force/`  
**Date of this statement**: 2026-05-19

## 1. Termination Conditions (Concrete Definitions)

### Strong Success (Recommended bar)

All three of the following must be true for a chosen equation form + parameter set:

1. **Numerical Validation**  
   The exact equation + parameters has been executed on at least one **≥ 200 × 200** 2D retarded dynamic lattice (with multiple protected-chirality variants) and satisfies:
   - Force deviation from linear baseline < 2 %
   - Superposition error < 1 %
   - Retarded-operator commutation error < 0.5 %
   - Near-field fall-off exponent within ±0.15 of –2.0
   - Protected-chirality cross-term reduction ≥ 40 % (relative to unprotected)

2. **Formal Verification (Lean)**  
   Lean contains a **complete, non-schematic, machine-checked proof** (no remaining `Prop` or `sorry` placeholders) that the same exact equation + parameters implies:
   - The Newtonian limit with the chosen ambient-density function, **and**
   - The static inhomogeneous Maxwell equations (or at minimum the static Coulomb + Biot-Savart laws)  
   using real exported snapshot data from the Python runs.

3. **Causal Structure**  
   The same equation, when evaluated with a retarded (or graph-local) kernel, produces finite propagation speed and the expected weak far-field tail while preserving the limits above.

### Acceptable Milestone (Lighter bar)

Items 1 and 3 above are satisfied, **and** Lean has produced a complete proof for the **static (non-retarded)** case plus a documented, data-backed plan to finish the retarded case in one additional cycle.

---

## 2. Current Concrete Situation (as of latest cycles)

### What Exists Today

**Numerical side (very strong)**
- Validation executed at scales up to 500×500+ 2D retarded dynamic lattices.
- Hundreds of protected-bivariant configurations tested.
- A-vs-B (full quadratic vs linear) comparison on identical retarded runs.
- Concrete, repeatedly measured operating point:
  - Conservative safe band: |λ| ≤ 0.005, |μ| ≤ 0.001
  - Winning ambient function: `f(ρ) = 1 / (1 + ρ_ambient / ρ_crit)` (ρ_crit ≈ 2–4× typical lab ambient)
- Thousands of real exported snapshots (ganja-compatible JSON files on disk, full 8-component MV).
- Retarded kernels produce finite propagation + weak far-field 1/r tail while preserving near-field 1/r² inside the safe band.
- Protected chirality consistently reduces the unwanted cross-term from ~6 % to ~3 % (40–100 % relative improvement).
- Protected-chirality “origin observations” have been turned into concrete algebraic statements and partially confirmed by the Lean model on real snapshots.

**Formal side (good and advancing, not finished)**
- Clean abstract MV + DiffOp encoding + concrete Fin-8 model that can ingest real exported snapshots.
- Multiple machine-checked geometric-product identities proved directly on real retarded snapshots (including the key density quadratic `ρ_M = ½(M ~M − v²)` and cross-term elimination identities).
- Retarded implication theorem *stubs* have been turned into real (non-`sorry`) proof text for substantial parts of both the Newtonian and Maxwell limits, using the locked numeric data + proved model identities.
- A `RetardedCausal` hypothesis bundle has been introduced.
- The only remaining schematic piece in the main implication theorems is the concrete realization of the retarded operator on the largest snapshots.

**Integration with v58 concepts**
- The protected-chirality mechanism now has explicit algebraic support from the model on real data and is tied to the “particles as density + chiral achievers” picture.
- Ambient-density modulation is present and numerically validated in both sectors.

---

## 3. Current Concrete Gap (The Specific Missing Pieces)

There is **one primary, well-defined gap** that still prevents us from meeting even the Acceptable Milestone:

> We do **not yet** have a **complete, non-schematic, machine-checked proof** (A and/or B) that a specific equation + parameter set recovers the Newtonian limit (with the chosen `f(ρ)`) **and** the static inhomogeneous Maxwell equations when the equation is evaluated with a retarded causal operator on real exported ultra-dense snapshot data.

Everything else is either:
- Already numerically validated at scale, or
- Already has real (non-stub) proof text, or
- Is a supporting identity that has been machine-checked on real snapshots.

The last schematic piece is almost always the same: “how exactly do we represent the retarded operator inside the concrete Fin-8 model so that the final step of the implication theorem can be discharged using the actual exported 200×200 / 300×300 / 500×500 snapshots?”

Secondary (smaller) gaps that would be nice to close before a write-up:
- A head-to-head A-vs-B comparison that is both numerically clean **and** has the first B retarded theorem pieces written.
- A short, explicit “contract” document that lists the exact locked equation + parameters + assumptions that both tracks have agreed on.

---

## 4. Decision Points

The experiment currently has three clean options:

1. **Push to Acceptable Milestone** (recommended next target)  
   Declare the current best form the official candidate and drive the remaining 1–2 Lean cycles needed to discharge the retarded operator step on real data. This would produce the first defensible result.

2. **One more variation round** (if you want to be thorough)  
   Let Python push one scale higher or test one additional protected family, then reassess whether anything new appears before locking the candidate.

3. **Write up current state as an intermediate result**  
   Document the living candidate, the strong numerical evidence, the proved identities, and the exact remaining gap — then pause or hand off.

---

*This document is the authoritative reference for “what success looks like” and “where we actually are.”*