# BRIEF — Unified quantized harmonic graph law

**To:** Claude (fable/max), headless review + symbolic derivation  
**From:** SCP v93 programme  
**Date:** 2026-08-08  
**Mode:** **Symbolic only** (Maxima + math). Write **Go** for the
quantized harmonic graph formalization. **No freecell/fabrun simulation
campaign** in this brief — simulation is a later gate if the theory is
coherent.

**Working directory / write all outputs under:**  
`/home/d/code/scp/v93/theory/qgraph/`

You have `maxima` on PATH. Use it for algebraic checks (saved scripts in
`qgraph/maxima/`). Implement the discrete theory as a small Go module under
`qgraph/` (types, graphs, occupation, shell operators, charge deficit —
unit tests with `go test`, no link to freecell required).

---

## 1. Programme constraints (non-negotiable)

1. **One law:** energy never destroyed; only changes mode. Space is a mode.  
2. **No background** — no immortal lattice of passive sites as ontology.  
3. **No imported U(1) charge / electrons as species** — charge and shells must
   *emerge* from harmonic structure.  
4. Free-cell is the physical substrate later; **here** you abstract to a
   **finite undirected (or oriented) graph** with complex amplitudes on
   vertices and/or edges — a **quantized harmonic graph**.  
5. Pre-v89 Q-ball / continuum PDEs: do **not** reimport.

---

## 2. What must unify under ONE law

The user postulate: **absorption, emission, electron shells, and charge
(via incomplete harmonic shells)** are the same mechanism viewed at
different moments / closures.

| phenomenon | empirical / law hook in v93 |
|---|---|
| **Emission** | Radiance: \(O = k_{\mathrm{rad}} C x^{p}\), \(p=4\), \(x=E_m/C\). Matter→field, atom-quantized. |
| **Absorption** | Door condensation field→matter (workfn, beats, atoms). **Not** yet dual T_bath⁴. |
| **Shells** | Es well + field modes; cap refusal (PAULI-0 seed); no charge yet. |
| **Charge** | Absent on free-cell. User: incomplete harmonic shells ⇒ charge-like asymmetry. |
| **Sustain** | T=3000: quiet unitary + q_detune HOLDs multi-hot clump; noise box-fills; rad drains. |

**Goal:** write a single discrete law on a harmonic graph such that:

- Emission and absorption are **time-reversed / dual** book operations on
  the same mode ladder (Stefan-like powers derived or motivated, not
  only postulated).  
- Bound shells are **closed progressive harmonics** (cycles / integer
  winding / integer occupation of graph eigenmodes).  
- **Charge** is a **topological or occupancy defect**: incomplete shell
  (missing occupation, broken cycle, half-filled multiplet) with a
  conserved defect quantum — not an imported U(1) field.  
- Temperature of matter ≈ excitation / cycle rate of the bound stack;
  vacuum field is **not** automatically thermometric.

---

## 3. Suggested formal objects (you may refine)

Graph \(G=(V,E)\). Complex amplitude \(\psi_v\) or \(\psi_e\).  
Occupation \(n_\alpha \in \mathbb{Z}_{\ge 0}\) of graph Fourier / Laplacian
eigenmodes \(\phi_\alpha\), \(\omega_\alpha\).

Action atom: \(\varepsilon_\alpha = A_0 \omega_\alpha / 2\pi\) (matches programme).

Candidate dual rates (you must derive or reject):

\[
O_\alpha = \kappa\, n_\alpha^{p} \quad\text{or}\quad
O_\alpha = \kappa\, (n_\alpha \omega_\alpha)^{p}/something
\]

\[
I_\alpha = \kappa\, n_{\mathrm{bath},\alpha}^{p}
\quad\text{or interference form involving } J\sim \mathrm{Im}(\psi^*\nabla\psi)
\]

Shells: partitions of mode set into closed subsets \(S_k\) with
\(\sum_{\alpha\in S_k} n_\alpha = N_k\) magic or capacity \(c_k\).

Charge candidate:

\[
Q = \sum_k \big(N_k^{\mathrm{full}} - N_k^{\mathrm{occ}}\big) \sigma_k
\]

or winding deficit on cycles. Show conservation of \(Q\) under your
update if absorption/emission only move complete shell quanta or
paired defects.

**Stefan \(p=4\):** either derive from 3D mode density + detailed balance on
the graph continuum limit, or show it is an independent geometric axiom
and how absorption must match.

---

## 4. Deliverables (write these files)

| file | content |
|---|---|
| `qgraph/docs/UNIFIED_LAW.md` | Main theory: definitions, law, dual I/O, shells, charge, T, Stefan |
| `qgraph/docs/MAXIMA_NOTES.md` | What was checked in Maxima + commands |
| `qgraph/maxima/*.mac` | Reproducible Maxima scripts |
| `qgraph/*.go` + `go.mod` if needed | Graph, modes, occupation, shell fill, charge defect, rate laws |
| `qgraph/*_test.go` | Unit tests for conservation identities, shell capacity, Q |
| `qgraph/docs/SIM_PREDICTIONS.md` | **Falsifiable** freecell/fab predictions (knobs, bars) if theory is adopted |
| `qgraph/docs/VERDICT.md` | Success / gaps / whether sim should be queued (YES/NO + why) |

---

## 5. Success criteria for VERDICT.md

**YES (queue sim)** if you have:

1. One law table with I and O dual in the same language.  
2. Shells as closed harmonics with discrete capacity.  
3. Charge as incomplete-shell defect with a conservation statement.  
4. At least two sharp, implementable predictions distinguishing this law
   from “radiance only + noise fuel.”  
5. Go tests pass (`go test ./...` under qgraph).

**NO** if dual absorption collapses to handwaved T_bath⁴, or charge is
just renamed Em, or shells require imported orbitals.

---

## 6. Context paths (read as needed)

- `/home/d/code/scp/v93/theory/METABOLIC_MATTER.md`  
- `/home/d/code/scp/v93/theory/T_HARMONICS_AND_WIND.md`  
- `/home/d/code/scp/v93/theory/PROGRESSIVE_CELL.md`  
- `/home/d/code/scp/v93/theory/BOHMIAN_HYBRID.md` (§10 method morph)  
- `/home/d/code/scp/v93/theory/MECHANISM.md`  
- `/home/d/code/scp/v93/laws/laws_V3a.cfg` (k_rad, p_rad=4, workfn)  
- `/space/scp/v93/README.md` (T⁴ status note)

---

## 7. Explicit non-goals for this turn

- Do not run freecell/fabrun batteries.  
- Do not edit freecell.c / fab step physics (Go theory package only).  
- Do not invent multi-fabric C/Q/L from v75.

Work until the deliverables above exist and `go test` is green (or document
why tests are blocked). End with a clear VERDICT.md.
