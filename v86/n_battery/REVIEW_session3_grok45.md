I'll attack the reversal hard: verify U(1)³ from the kernel, then re-derive whether the shape-channel flip is real physics or an artifact.# Adversarial review — HC-6 reversal

No files modified. Kernel re-read, diags re-analysed, rates recomputed.

---

## Part 1 — U(1)³ and the void observable

### Finding 1 — U(1)³ is exact at η=0, ungauged  
**Severity: OK**  
**Claim:** \(V=V_t(s)\), \(s=\prod_a|\Phi_a|^2\) ⇒ independent phase rotations ⇒ each \(Q_a\) conserved.

**Kernel (`compute_forces_complex`, η=0, no gauge):**
- Potential force: \(-2 V_t'(s)\,\phi_a\cdot\mathrm{prod\_rest}[a]\) with \(\mathrm{prod\_rest}[a]=\prod_{b\neq a}|\Phi_b|^2\)
- Force ∥ \(\Phi_a\); \(V_t'\) and \(\mathrm{prod\_rest}\) are phase-invariant
- Mass and kinetic terms are per-component
- η curl terms are **off** (`ETA=0`); no gauge links

Independent phases \(\Phi_a\to e^{i\theta_a}\Phi_a\) are symmetries. Three Noether charges follow.

**Measured (your no-sponge runs, rechecked):**  
worst \(|\Delta Q_a/Q_a|\): target \(6.97\times 10^{-13}\), control \(3.69\times 10^{-13}\). Matches exact conservation to integrator noise.

**What breaks U(1)³ (not in these runs):**
| Term | Effect |
|------|--------|
| η>0 curl(Φ,Θ) | mixes sectors; can spoil pure Φ-phase symmetry depending on Θ content |
| single U(1) gauge, equal charges | still conserves each \(Q_a\) globally (diagonal gauging); relative U(1)² remain |
| multi-fabric / non-product \(V\) | can break |

**Correction:** claim is right for **this** run class. Do not globalise to “the full Cosserat always has U(1)³.”

---

### Finding 2 — HC-6-as-redistribution was malformed; GSS was not  
**Severity: MAJOR** (census design)  
**Claim:** original HC-6 (“decay by charge redistribution / partition fractions”) is void.

**Agreed.** Under exact \(Q_a\) conservation, \(x_a=Q_a/Q_\mathrm{tot}\) is frozen up to \(O(10^{-13})\). Any null or positive built on \(\max|\Delta x_a|\) was measuring a conservation identity, not dynamics.

**But the deeper point:**  
`GROUNDING.md` §1 already states *“three commuting U(1)s at η=0”* and defines multi-charge GSS with \(D_{ab}=\partial Q_a/\partial\omega_b\). The theory side **knew** U(1)³. The rung language did not:

- `PROPOSAL.md`: “decay to sector minima”  
- `NEXT_PROGRAM.md` / HC-5: “decay to the sector minimum”  
- early `hc6_verdict.c`: partition drift as the discriminator  

That is a **malformed operationalisation**, not a failure of multi-charge GSS.

**Corrected statement:**  
GSS still applies. Instability means the stationary profile is **not a constrained minimizer of \(E\) at fixed \((Q_0,Q_1,Q_2)\)** → deformation / radiation / fission **at fixed \(Q_a\)**. It does **not** mean charge flows between U(1) factors.

---

### Finding 3 — U(1)³ makes GSS more applicable, not less  
**Severity: OK**  
**Claim / Q(c):** does U(1)³ change the criterion \(n(H_\omega)=n(D)\)?

**No change to the criterion; clearer scope.**

- Multi-charge GSS (Grillakis–Shatah–Strauss / Comech–Pelinovsky style) is built for several conserved charges and chemical potentials \(\omega_a\).
- \(D_{ab}=\partial Q_a/\partial\omega_b\) is the Hessian of the reduced action on the charge manifold.
- Orbital stability ⇔ \(n(H_\omega)=n(D)\) (with the usual caveats: orbital not asymptotic; continuum radiation; box).

U(1)³ means the charge lattice is three-dimensional and **true** conserved coordinates of the reduced dynamics. That makes \(D\) **more** fundamental, not less.

**What U(1)³ kills** is only the story “unstable partitions relax by reshuffling \(Q_a\).”  
**What it leaves** is “unstable partitions deform at fixed \(Q_a\).”

---

## Part 2 — Shape-channel “positive” result

Recomputed from diags:

| arm | sponge | \(s_\mathrm{end}/s_0\) | rate (60% fit) | \(R^2\) | \(\max\|\Delta Q_a/Q\|\) |
|-----|--------|----------------------|----------------|-------|-------------------------|
| target | no | **1.435** | \(2.32\times 10^{-3}\) | **0.89** | \(7\times 10^{-13}\) |
| control | no | **1.001** | \(2.7\times 10^{-6}\) | **0.001** | \(4\times 10^{-13}\) |
| mono | yes | 2.35 | \(5.47\times 10^{-3}\) | **0.52** | \(6\times 10^{-2}\) (sponge) |
| target | yes | 1.23 | \(2.94\times 10^{-3}\) | 0.88 | \(2\times 10^{-2}\) |

Rate ratio target_ns/control_ns ≈ **855×** (your 714× is the same comparison, different window).  
Control \(R^2\approx 0\): its “rate” is a fit to noise on a flat line.

---

### Finding 4 — \(s_\mathrm{max}\) is a legitimate *channel*, not a clean Lyapunov probe  
**Severity: MAJOR**  
**Claim:** \(s_\mathrm{max}\) growth diagnoses GSS instability; fitted rate is a growth rate.

**What \(s_\mathrm{max}\) can do:** peak of \(s=\prod|\Phi_a|^2\) responds to core compression / rarefaction. Target no-sponge also shows coherent companions:
- \(r_\mathrm{core}\): 6.11 → 5.57 (−9%)  
- \(\phi_\mathrm{max}\): 0.568 → 0.617  
- \(P_\mathrm{int}\): 4.31 → 5.63  
- \(E\) conserved to \(\sim 10^{-4}\)  

So this is not a single-voxel spike with energy blown up.

**What it cannot do cleanly:**

1. **Non-monotonic start (your worry is right).**  
   target_ns: \(s/s_0\) min **0.902 at \(t\approx 29.6\)**, then rise. That is seed relaxation / breathing kick from lattice projection, **not** pure exponential instability from \(t=0\).

2. **Oscillation load.**  
   target_ns: **199 sign changes** in \(\mathrm{d}s\) over the run; residual about linear trend peak-to-peak **0.23** on top of secular ~0.42. Growth is **secular + large breathing**, not a clean mode.

3. **Log-linear “rate” over first 60%.**  
   - Target \(R^2=0.89\): usable as a crude slope, not a Lyapunov exponent.  
   - Control \(R^2=0.001\): rate **undefined**; ratio “714×” is **rate / noise**.  
   - Mono \(R^2=0.52\) with clear 3.4↔2.3 swings: rate **not meaningful**.

**Correction:** report  
- secular amplitude \(\langle s\rangle_\mathrm{late}/\langle s\rangle_\mathrm{early}\) or envelope,  
- \(\Delta r_\mathrm{core}\), \(\Delta P_\mathrm{int}\),  
- and **do not** headline a rate ratio against a flat control.  
My recomputed amplitudes: target late−early **+0.42**, control **+0.002** (~200× in *amplitude*, not 714× in rate). Still a large separation; the 714× number is the wrong statistic.

---

### Finding 5 — Mono arm undercuts a *flavoured*-GSS reading  
**Severity: MAJOR**  
**Claim (implicit):** n(D)=0 vs n(D)=1 separation ⇒ GSS index mismatch (flavoured).

Mono past-turn (sponged) reaches \(s/s_0\) up to **~3.4** and ends at **2.35**, **larger** deformation than the flavoured target. It has **no** flavour structure and \(n(D)=0\) only because it sits past the monochromatic VK turn.

So the strongest effect in the table is consistent with **“past the VK turn ⇒ shape runs away”**, which both mono and target share. Control is below the turn. That is exactly the confound you still have in the shape channel.

**Correction:** without a **no-sponge monochromatic** arm at matched \(Q\) and extent, you cannot attribute the target’s shape motion to *flavoured* index structure rather than ordinary VK.

---

### Finding 6 — Resolution / seed confounds remain  
**Severity: MAJOR**  
**Claim / fear:** dx=0.3465, ~6.6 cells per half-radius ⇒ discretisation drift.

**Evidence for confound:**
| | Target | Control |
|--|--------|---------|
| continuum \(E\) | 177.40 | 178.70 |
| lattice \(E(0)\) | **176.75** (−0.37%) | **178.54** (−0.09%) |
| continuum \(Q\) | 116.44 | 117.95 |
| lattice \(Q(0)\) | **116.03** (−0.35%) | **117.95** (~0) |
| \(r(f>0.01)\) | ~15 | ~9.5 |
| ω high component | 1.495 (κ≈0.12) | 1.465 (κ≈0.32) |

Target is the softer, worse-projected seed. Extra shape motion can be **relaxation of projection error** along a flatter energy landscape (near/ past VK), not necessarily a linear unstable eigenmode.

**N=160:** only \(t\lesssim 2\) in `hc6_target_ns160_diag.tsv` — **no convergence information yet**.

**What would convince me:**

| Outcome at matched \(t\sim 100\)–\(200\) | Reading |
|----------------------------------------|---------|
| \((s/s_0-1)_{160}/(s/s_0-1)_{128}\in[0.85,1.15]\) and same sign of \(\Delta r_\mathrm{core}\) | consistent with continuum instability |
| ratio \(\lesssim 0.5\) (much less growth when refined) | mostly discretisation / seed transient |
| growth **increases** under refinement | under-resolved, not trusted |
| early dip to ~0.9 disappears or shrinks a lot at N=160 | dip was lattice transient |
| dip persists at same depth | physical breathing / mode content |

Also demand **same seed path** (shooter → same generator) so you are not comparing different projection errors.

---

### Finding 7 — Target/control still confounded in the shape channel  
**Severity: MAJOR** (unchanged from B1)  
**Claim:** pair isolates \(n(D)\).

They still differ in \(\bar\omega\), distance to VK turn, extent, and seed quality. Shape-channel rewrite **does not** remove that.  
Best control remains: **same mean ω / same extent / same distance to turn**, differing only in partition — or mono vs flavoured **both** past turn, **both** no sponge.

---

### Finding 8 — Strongest non-GSS alternative and how to kill it  
**Severity: MAJOR**

**Best counter-hypothesis:**  
> Both arms that move sit **at or past the monochromatic VK turning region** of the soft thick-wall branch. The control sits **below** it on a stiffer, better-resolved object. Observed \(s_\mathrm{max}\) growth is **ordinary VK shape instability / breathing of a soft Q-ball**, amplified by seed projection error and marginal dx — not a multi-charge GSS phenomenon that requires \(n(D)=0\) flavoured structure.

**Why it fits:** mono (no flavour) moves as much or more; target has worse seed residual; non-monotonic early \(s\); large oscillatory component; control is simply deeper in the stable VK region.

**How to kill it (pre-register):**
1. **No-sponge mono** at matched \(Q,E\), extent ≈ target: if mono ≈ target in secular \(\Delta s\) and \(\Delta r_\mathrm{core}\), counter-hypothesis wins.  
2. **Flavoured control at same \(\bar\omega\)** but \(n(D)=1\) if any such point exists interior to the volume (hard near the turn).  
3. **N=160 (and ideally N=192)** secular \(\Delta s\) convergence within ~15%.  
4. **Mode content:** project early evolution onto the known l=0 dilational BdG mode; require growth of that mode, not only broadband \(s_\mathrm{max}\).  
5. **Energy landscape:** if \(E\) at fixed \(Q_a\) decreases while shape moves (true gradient descent on the constrained manifold), that supports constrained instability; pure oscillation at fixed \(E\) is breathing.

Until (1)+(3) pass, **do not** claim GSS converse confirmation.

---

### Finding 9 — What *is* established  
**Severity: OK** (limited)

With no sponge and \(Q_a\) conserved to \(10^{-13}\):
- Control shape is **stable** at the \(10^{-3}\) level in \(s\) for \(t\sim 200\).  
- Target shape **secularly deforms** (\(+43\%\) in \(s_\mathrm{max}\), core shrinks, \(P_\mathrm{int}\) up) at **fixed** \(Q_a\).  
- That is **compatible with** constrained instability past the VK / index-mismatch region.  
- It is **not yet** a clean demonstration that the cause is \(n(H)=n(D)\) mismatch as opposed to VK+resolution+seed.

---

## Verdict on the REVERSAL

| Piece | Verdict |
|-------|---------|
| Old observable void under U(1)³ | **SOUND** |
| Census HC-6 “redistribution” language malformed | **SOUND** (MAJOR design bug, not a run bug) |
| GSS criterion still right; applies at fixed \(Q_a\) | **SOUND** |
| Shape channel is the right place to look | **SOUND** |
| “Result flips to positive / 714× rate = GSS confirmed” | **NOT YET ESTABLISHED** (overclaimed) |

### Overall: **SOUND WITH CORRECTIONS** on the *void*, **NOT YET ESTABLISHED** on the *positive flip*

**What you should write down now:**
1. Partition-fraction HC-6 is void by U(1)³ — stand by that.  
2. GSS converse, if tested, is **shape dynamics at fixed \(Q_a\)**.  
3. Present no-sponge data show a **real shape difference** target vs control, but  
4. **Do not** claim confirmation until no-sponge mono + dx convergence close the VK/resolution counter-hypothesis.  
5. Retire “714× growth rate”; use secular amplitude / \(r_\mathrm{core}\) / mode projection, and never divide by a control rate with \(R^2\approx 0\).

**One line for the lab book:**  
*The reversal correctly kills the old observable; it does not yet earn a GSS converse confirmation.*
