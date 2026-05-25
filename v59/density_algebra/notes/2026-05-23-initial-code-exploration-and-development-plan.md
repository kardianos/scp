# 2026-05-23: Initial Exploration of v59/density_algebra/python/ and Framework Development Kickoff

**Agent**: Grok Build subagent (scientific Python / numerical simulation specialist for density algebra research)
**Task**: Aggressively mature the Python code into production-quality full-theory analysis framework incorporating v58 living candidate + f(ρ) ambient modulation + quadratic terms + full 8-comp MV/octonion-like content + explicit protection projectors + direct ρ_M wells + rich observables + sweeps + feedback to HYPOTHESES.md.

## Exploration Summary (Thorough Review Performed)

### Files Inspected
- `v59/density_algebra/python/improved_density_protection_scan.py` (latest 8-comp version):
  - Good skeleton: explicit 8-comp masks for PROTECTION_CHOICES (unprotected_full, protected_L=s+vecs, protected_F=s+bivs+tri, protected_LF=full).
  - Uses ga.MV directly for make_8component_lump (dict terms for components).
  - Has rho_M(mv, v_sq) = 0.5 * ( (mv * ~mv)[0] - v^2 ) — correct direct definition matching Maxima.
  - make_mask, COMPONENT_NAMES = ['s','e1',..,'e123'] (Cl(3) full 8D proxy for octonion).
  - **Critical flaws identified**:
    - compute_rho_at_probe is a *toy non-physical hack*: uses delay logic that returns 0 for most t/probe, raw=source['mv'][0] (index error! MV __getitem__ expects blade tuple like (), not int 0; would crash or wrong), protection_factor crude, no real retarded Omega, no f(ρ), no quadratic, no actual dynamics.
    - Consequently all reports show peak/avg_density=0.0 — useless for testing wells.
    - run_improved_scan uses fixed toy 1D, no real integration with v58 retarded logic.
    - No plots, limited observables (only avg/peak dens + toy cross/comm).
- `v59/density_algebra/python/targeted_density_protection_scan.py` (earlier):
  - Better: actually calls `run_retarded_simulation` from v58 `retarded_dynamic_scan.py` (full living-candidate retarded 1D with history, f_winning, quadratic λ corrections, protected_biv flag).
  - Extracts density_proxy from avg_force_c (indirect, not direct ρ_M).
  - Good report structure (JSON with raw_stats).
  - **Limitations**: Still only scalar rho + biv chi (not full 8-comp MV algebraic content for lumps), protection only binary protected_biv (not general 8-comp masks/projectors), density proxy is force magnitude not ρ_M(M·~M), no multi-probe well profiling, no 2D/3D, no superposition tests, reports also had 0s in the snapshot (likely run timing or param where forces small? but structure good).
- Core dependency `v58/pregeometric/unified_multivector_force/python/ga.py`:
  - Solid lightweight Cl(3,0) MV: dict[blade_tuple, coeff], full geom product, grade(), reverse()/~ , norm2(), scalar/vector/bivector factories, inner/outer.
  - 8 components exactly: grades 0+1+2+3 map to octonion 1+e1..e7.
  - Used throughout v58 experiments.
- `v58/.../retarded_dynamic_scan.py` (and 2d_retarded_grid_scan.py, quadratic_f_scan.py):
  - Living candidate implementation: retarded causal history, rho_ambient_retarded, f_winning(ρ) = 1/(1+ρ/ρ_crit) for both g and em channels, compute_retarded_omega with grav (vec from rho * f_g) + biv (from chi), then quadratic iteration (λ Ω² terms + μ norm), force_on_test(Ω * M_t).grade(1), commutation/leakage measures, protected_biv variant (restricts test biv plane).
  - 1D moving source + fixed probes; 2D grid version exists for better spatial.
  - Produces cross_term_pct, comm_err, exponent proxies, safe band for |λ|≲0.005-0.01.
  - **Key for extension**: Can be generalized — replace scalar rho/chi with full MV "algebraic state" per lump, apply component masks as projectors before contributing to omega or to M_t, accumulate local effective M(x,t) = sum retarded contribs, then direct rho_M(local_M) for well measurement at grid of probes.
- Reports (2 JSONs dated 2026-05-23):
  - Structure good for feeding hypotheses (timestamp, description, results list, links to HYPOTHESES/CONSTRAINTS).
  - But data degenerate (dens=0) due to above issues — first task: generate *non-trivial* quantitative data.
- `v59/density_algebra/HYPOTHESES.md` + `CONSTRAINTS_AND_TROUGHS.md`:
  - Excellent framing: particles = stable high-ρ_M density wells protected by algebraic tech (L=28d, F=35d, L+F=63d graded selection from Cl(7)_even); S³ |ξ|=1/√2 as protection budget sweet spot; 21-factors as internal cost of density gradients; forcing via density-well depth vs protection cost vs force-separation budget vs relationship trough.
  - Direct asks for: vary protection/algebraic content in multivector code, measure peak ρ_M + cleanliness; quantitative signatures for each hyp.
  - Vocabulary: Density Well (depth/width/stability), Protection Budget/Cost (fraction of DOF spent), Relationship Trough, Geometric/Causal Stability Margin, Force-Channel Separation Budget.
  - Maxima `full_octonion_perturbation.mac` provides analytic target: ρ_M = ½(M M̃ - v²), protection projectors P_L/P_F on 8 (or 7 imag) comps, 2nd variation/Hessian eigenvalues for stability of wells under different P.
- Other: `maxima/perturbation_analysis_living_candidate.mac`, lean/, NEXT_STEPS.md (calls for quantitative proxy in v58 code for forcing, budget accounting).

**Memory recall**: Prior session notes on full 8D octonion in Maxima + Python driver on ga.MV for algebraic content/protection as explicit masks; scaling Maxima symbolic + Python numerical for perturbation/Hessian of ρ_M; integrate with living candidate f(ρ), quadratic.

### Theory Alignment Requirements (from task + docs)
Must incorporate:
- Full 8-comp MV / octonion-like algebraic content for sources/lumps/tests.
- Protection = explicit projectors/masks on the 8 components (L/F/ full / custom subsets; maps to graded sector selection Hyp1).
- Direct ρ_M measurement at *multiple probe locations* inside/around protected lumps (not force proxy).
- Full/high-fid living candidate: retarded causal dynamics (history, light cones), ambient f_g(ρ), f_em(ρ) modulation (f_winning or variants), quadratic self-interaction λ,μ terms.
- Retarded for causality.
- Sweeps: λ,μ, protection mask/strength, initial algebraic content (which comps in make_lump), grid/resolution, v_crit for f, etc.
- Rich observables (directly testable vs hyp forcing):
  - Density well: peak_ρ_M, avg_ρ_M (bg-subtracted), well width (FWHM or std of profile), depth vs bg, stability (var over t or response to pert).
  - Protection cost: num_active /8 , or effective DOF spent.
  - Cross-term suppression (neutral vs "charged" i.e. different masks on test), commutation error (leakage between grades/channels).
  - Fall-off exponents (fit ρ_M(r) or force(r) ~ r^{-α} near/far).
  - Superposition quality (two lumps: ρ_M(sum) vs sum ρ_M; deviation measures interaction cleanliness or trough).
  - Force cleanliness, geometric stability proxies.
- Continuous loop: implement, run, measure, propose enhancement (e.g. 2D grid, better f forms, numerical Hessian matching Maxima, adaptive masks, bootstrap from analytic), test, output structured (JSON + tables + plots + summary for hyp validation).
- Clean/extensible: docstrings, config dicts, logging, modular (core sim vs analysis vs reporting), use SFA? (but for now since density algebra not full field sim, np arrays + JSON ok; later if needed).
- Output usable for HYPOTHESES validation/falsification (e.g. "protected_L achieves 2.3x deeper well than full at same λ, with 30% better comm_err — supports L as 'cheapest' tech per Hyp1").

**Current state gaps vs mature**:
- No working direct ρ_M wells in reports.
- Protection too coarse (binary or 4 toy choices).
- Dynamics not fully 8-comp aware for algebraic content.
- No spatial profiling of wells (1 probe only).
- No 2D/3D, limited stats, no viz, no Maxima cross-check.
- Not production: minimal error handling, no config, hard-coded, no plots.

## Development Plan (Aggressive, Iterative)
**Phase 1 (today)**: 
- Create notes/ (done).
- Write this exploration note.
- Design & implement core `density_algebra_lab.py` (or enhance/fork improved_ into production lib): 
  - Flexible ProtectionMask class (from list of active indices or name, projector fn on MV: zero inactive comps).
  - Enhanced make_protected_lump(mask, params) returning MV with realistic amplitudes per grade (or user spec).
  - Direct rho_M, also local field builder: retarded sum of projected MV contribs → effective M_field at (x,t), with f(ρ_amb) where ρ_amb from |M| or scalar, quadratic corrections.
  - Reusable simulator: run_protected_retarded_scan( mask, lambda_nl, mu_nl, n_steps, probes=[...], use_2d=False, ... ) returning dict of rich metrics + time series + spatial profiles.
  - Adapt logic from retarded_dynamic_scan + 2d_... but MV-native for algebraic content (source 'mv' instead of rho/chi; contrib = P(mask, mv) * kernel * f(ρ from scalar(M) or rho_M) ).
  - Observables computer: well_depth, well_width (from probe grid ρ_M(r)), peak_loc, stability, cross_suppression (run with test_mask variants), comm_err (extend measure), falloff_fit (scipy curve_fit on log ρ_M vs log r), superposition_dev (add two sources).
  - Config system (dataclass or dict), logging.
- Fix/run the existing improved_ as smoke test, generate first non-zero report.
- Create driver `run_mature_scan.py` for systematic sweeps (λ grid, all PROTECTION_CHOICES + more realistic e.g. even-grades-only, biv-only for protected chirality, vector-only, full-even, random subsets).
- Generate plots: ρ_M profile vs r for different masks; trade-off scatter (peak_ρ vs cross_term, colored by num_active); stability metrics.
- Write first rich JSON report + summary table in report.
- Date next note with results vs hypotheses (e.g. "L mask gives deeper wells at low λ supporting Hyp1 cheapest protection").

**Phase 2 (immediate follow-ups)**:
- Integrate 2D grid for true well imaging (radial profiles, anisotropy?).
- Numerical perturbation/Hessian: for fixed background M0 in mask, small random δM (all 8), finite-diff d²ρ_M / dδ_i dδ_j , eigvals of Hessian → stability spectrum per protection (compare to Maxima taylor).
- More f(ρ) variants (const, soft, sqrt, learned), μ sweeps, initial content amp sweeps.
- Superposition, multi-lump stability, retarded causality tests (lag effects on well depth).
- Bootstrap: use analytic from Maxima (load or hardcode low-order) to validate numerics.
- Extensibility: plugin for new protection types, new observables, export to ganja-like or sfa if scales to field.
- Auto-generate hypothesis feedback snippet: "Quantitative support for Hyp1: protected_L achieves peak_ρ_M = X ±σ > protected_F at λ=0.005 with comm_err 40% lower, consistent with L as low-budget tech for lepton-like."

**Phase 3+**: 
- 3D if needed (but 2D sufficient for wells), time-dependent well formation dynamics, link to ξ dynamic or S³ constraint (by making "radius" of MV norm part of protection).
- Statistical ensembles (many random initial phases), error bars.
- Visualization dashboard (multi-panel figs saved), pandas summary tables.
- Continuous: after each run, auto-analyze "does this falsify/support any hyp? propose next exp".
- Polish: type hints, pytest stubs if time, comprehensive doc, README in python/.

**Success Criteria**:
- Reports contain non-zero, physically plausible ρ_M wells with variation by protection mask (e.g. L better than full for cleanliness).
- Code is clean, importable as `from density_algebra_lab import DensityProtectionScanner`, configurable sweeps in <20 LOC.
- At least 3-5 rich experiments run, data in reports/, plots generated, progress notes written.
- Direct tie-back: numbers cited in future hyp updates or Lean.
- No breakage of existing v58 (pure additive in v59/python/).

**Risks/Mitigations**:
- Retarded + full MV + multiple probes per step may be slow for sweeps → start 1D few probes, vectorize where poss, limit n_steps=20-50, use numba? (check if avail), or precompute.
- Mapping 8-comp Cl(3) exactly to Cl(7) graded → treat as faithful low-dim analog for protection projector tests; note in docs.
- Getting ρ_M positive/deep: tune initial MV amps in make_lump, v^2 offset, or use |M|^2 scalar part directly as proxy if v=0 for wells; ensure f modulation uses local density estimate.

**Next Immediate Action**: Implement and test the core lab module + first production run of sweeps. All code changes via search_replace or write (prefer edit existing for improved_ first? but new module better for maturity). Update this note with results. Use todo? but since focused, direct execution.

This establishes the baseline. Aggressive development begins now. All future work references this exploration for continuity.
