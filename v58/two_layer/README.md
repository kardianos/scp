# Two-Layer Medium + Excitations Route (Placeholder)

**Date**: 2026-05-18  
**Relation to LAYER_CRITIQUE.md**: This folder corresponds to Route 1 ("Explicit Two-Layer 'Medium + Excitations' Model").

## Purpose

An explicit dynamical "medium state" layer (local node density ρ(x), stiffness, anisotropy, causal parameters, or effective metric coefficients) whose evolution is sourced by the presence and motion of excitations. Small-amplitude test waves, rays, or collective modes then propagate according to the instantaneous local effective geometry derived from the medium state.

## Status in v58

Light placeholder only. No implementation, no kernels, no equations of motion yet.

The physics *target* that any such two-layer construction must meet is defined in `../first_principles/EXPECTED_BEHAVIOR.md`. In particular:

- The medium-state → propagation mapping must be geometric and non-dispersive (same effective metric for all frequencies and species, |c_GW − c| ≲ 10^{-15}).
- Back-reaction must be self-consistent and universal.
- The construction remains grid-based and is therefore provisional; its primary value is as a diagnostic that forces the geometric-propagation question to be explicit and measurable.

## Relation to Other Folders

- `first_principles/` supplies the acceptance criteria.
- `pregeometric/` is the conceptually preferred long-term direction (no presupposed grid or metric).
- This route is retained as a pragmatic, lower-risk bridge for testing whether an explicit medium layer can produce both soliton forces *and* genuine geometric effects before committing to a full pre-geometric substrate.

Further work here will be opened only after the first-principles requirements are stable and a concrete minimal medium dynamics proposal has been scored against them. See DIRECTION.md for prioritization.