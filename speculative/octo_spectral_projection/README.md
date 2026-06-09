# Speculative — Octo-space frequency/phase projection

This directory contains bounded, documented experiments exploring the reframing that emerged after the honest negatives in v59–v61:

**Core idea**: The octonionic algebra (octo-space / Cl(7)_even ≅ ℂ⊗ℍ⊗𝕆) is the primary information substrate. What we perceive as 3+1 physics, Brannen kernels, Koide Q, the Brannen phase φ=2/9, gauge integers, etc. is the output of a frequency/spectral transform + phase (or information) projection applied to algebraic data. The "before it comes real" layer is where the actual dynamics live; the configuration-space field theory of v59–v61 was asking the algebra to do the wrong kind of work.

Experiments follow the project's 0X_findings.md discipline: every cut gets a script + a findings file that states the exact operator, the numerical result, the precise obstruction, and what the next minimal upgrade would be. Clean negatives are first-class outputs.

Current cuts (as of 2026-05-28):

- `01_*` — Classical DFT on the Z3 generation cycle + real algebra-derived weights. Clean negative: real scalars cannot manufacture the required phase offset.
- `02_*` — Concrete frequency analysis in the 7 imaginary dimensions using the repo's own Fano table (L_ei left-multiplication operators with spectrum {0, ±i}, Fano graph Laplacian gap of 7, G2 Casimir C₂(7)=2). This is the "freq analysis you can do in 7 dims" the user asked for.
- `03_*` — First wiring of the 7D spectral content through the L-grade complex structures J (J² ≈ −I) to produce a phase + three masses. Partial signal (phase in ballpark) with clear obstruction.
- `04_*` — Explicit compression step (projection onto top-k dominant Fano Laplacian eigenvectors) inserted before the J-rotation. Tests whether "compressing" the 7D frequency content improves the output.

See the individual `0X_findings.md` files for the detailed obstruction and upgrade path for each cut.

This line is a direct response to the derivation obstruction in `v60/gravity_recast/09` (Plebański-style parent is an independent posit) and the rank-tension deflation in `v60/gaps/rank_tension/01`. If the algebra is primarily spectral and only its projection is "real," then many of the configuration-space contradictions become category errors rather than technical gaps.

Next cut (if pursued): replace the toy frequency filter with a G2-highest-weight or actual Shulga S^7 harmonic kernel, couple the compressed 7D layer to the Cl(3,1) spacetime factor (v60/gravity_recast/07), and re-run the same test harness against the known Brannen numbers + structural integers.