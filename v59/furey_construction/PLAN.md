# Furey Construction Plan — Multi-Variant Attack on the Standard Model

**Date**: 2026-05-22
**Parent**: [`../README.md`](../README.md), [`../SUMMARY.md`](../SUMMARY.md)
**Goal**: Complete the cross-sector unification by building the full $\mathbb{C} \otimes \mathbb{H} \otimes \mathbb{O}$ (or analogous octonionic) construction, identifying the Standard Model gauge group as its natural symmetry breaking, and deriving α and G individually.

This plan lists multiple distinct approaches to be executed in sequence. **No giving up — every variant runs to completion or honest documented failure.**

---

## Context

After 11 steps in v59 we have:
- **Structurally**: Brannen form, Koide Q = 2/3, Brannen φ = 2/9, cross-sector ratio 21, three generations from triality — all derived from $G_2 \subset \text{Spin}(7) \subset \text{Spin}(8)$.
- **Empirically remaining**: overall lepton mass scale $a$, individual values of α and G, quark sector, full SM gauge structure.

The Furey construction $\mathbb{C} \otimes \mathbb{H} \otimes \mathbb{O}$ (64 real dimensions) is conjectured to provide the natural framework where these all emerge from one algebraic structure.

## Variants in Execution Order

### Variant A — Build the C⊗H⊗O algebra explicitly (Python)

**Goal**: Construct the 64-dimensional algebra explicitly; verify the standard multiplication table; demonstrate basic operations (multiplication, conjugation, inversion).

**Approach**: Tensor product of three normed division algebras. Use 4D quaternions and 8D octonions; the tensor product is a 64D non-associative algebra. Compute structure constants.

**Output**: `01_choh_algebra.py` + `01_findings.md`. Tabulated multiplication, verification of associators, and structure-constant statistics.

### Variant B — Identify SM idempotents

**Goal**: Find the projection operators in $\mathbb{C} \otimes \mathbb{H} \otimes \mathbb{O}$ whose stabilizer in the algebra's natural symmetry group is $\text{SU}(3) \times \text{SU}(2) \times \text{U}(1)$.

**Approach**: Following Furey (2014, 2018), look for primitive idempotents whose left-multiplication ideals decompose into SM representations. Use the Witt decomposition of $\mathbb{O} \otimes \mathbb{C} = \text{Cl}(6)$.

**Output**: `02_sm_idempotent.py` + `02_findings.md`.

### Variant C — Derive U(1)_em from algebra normalization

**Goal**: Compute the natural $U(1)_{em}$ coupling from the algebra without external normalization. Check if it gives α ≈ 1/137.

**Approach**: $U(1)_{em} \subset \text{SU}(2) \times \text{U}(1)_Y$ via the Standard Model symmetry breaking. In the algebra, the $U(1)$ orbit through the SM idempotent has a specific volume. The coupling is determined by this volume.

**Output**: `03_alpha_from_u1.py` + `03_findings.md`.

### Variant D — Test α prediction with the π²/2 conjecture

**Goal**: Verify whether $S_{\rm em} = \pi^2/2$ emerges as an instanton action on the natural Furey-construction quotient (likely $G_2 \backslash \text{Spin}(7)$ or similar).

**Approach**: Compute the instanton action $S = (8\pi^2)/g^2$ for $g = $ U(1)_em coupling derived in Variant C. Check if this equals $\pi^2/2$ or what fractional correction is needed.

**Output**: `04_alpha_prediction.py` + `04_findings.md`.

### Variant E — Quark sector via left-ideals

**Goal**: Extend the lepton-sector results to quarks. Compute mass ratios and charges for the up, down, strange, charm, bottom, top quarks.

**Approach**: The Furey construction places quarks in a 3-dim irrep of $\text{SU}(3)_c$ (color). Their mass operator is a Z₃-cyclic structure analogous to leptons but with color triality. Compute the predicted spectrum.

**Output**: `05_quark_sector.py` + `05_findings.md`.

### Variant F — Lean formal verification

**Goal**: Encode the C⊗H⊗O algebra and key identities in Lean 4, ensuring the structural claims are machine-checked.

**Approach**: Build minimal Lean files for:
- The algebra structure (multiplication tables, associators).
- The SM idempotent and its stabilizer.
- The Koide and Brannen identities from $G_2/\text{Spin}(7)$ ratios.

**Output**: `lean/ChOAlgebra.lean`, `lean/Triality.lean`, `lean/SMEmbedding.lean` + `lean/README.md`.

### Variant G — Gravity sector via Spin(7) coset

**Goal**: Identify gravity with a specific sector of the Furey algebra, derive G individually.

**Approach**: The scalar-density grade of the full algebra corresponds to gravity. Compute the natural coupling from the Spin(7)-orbit volume of the density-grade trace.

**Output**: `06_gravity_sector.py` + `06_findings.md`.

### Variant H — Consolidation and final assessment

**Goal**: Bring together all variant outputs into a coherent picture; identify what's been derived and what remains empirical.

**Output**: `ALL_FINDINGS.md` — consolidated final document.

## Termination Criteria

A variant is **complete** when:
- Its script runs to completion (no errors).
- Its findings document is written.
- The result (positive or null) is honestly recorded.

The whole plan is **complete** when all eight variants have status "complete" — regardless of how many produced positive numerical results. The point is to exhaust the obvious paths.

## Falsifiability

Each variant can:
- **Confirm a target**: produces a TYCHO match with ≤1 empirical parameter.
- **Fail with information**: rules out a specific identification, narrowing the next attempt.
- **Fail without information**: indistinguishable from random — flagged and discussed.

A variant that fails without information is the lowest-value outcome, but still tracked.

## Note on Scope

The full Furey program is large (decades of literature). This plan is a **focused first-pass execution** — building the algebra, identifying the SM, computing the natural couplings — not a complete reconstruction of mainstream octonionic-SM theory. Where this falls short, it's documented; where it succeeds, it's the v59 program's contribution beyond what existed before.
