# v93 Critical Review — theory, code, results, conclusions

Independent senior review. I read `REVIEW_REQUEST.md`, `v93/README.md`,
`v93/L1_FINDINGS.md`, the named regions of `v93/kernel/freecell.c`, the Go
mirror `v93/fab/step.go`, and `v92/consult/SUBQUARK_synthesis.md`, and I
re-ran the winding analyzer on the QUENCH `.fcs` outputs myself. Derivations
below are mine; where they coincide with the parallel grok review that simply
makes those code-level findings consensual. The places I add something new are
the discrete-time-Schrödinger identification of pass U, the imprint-wall
diagnostic at t=20, the Gauss-Seidel ordering confound, and the
no-C≡Go-at-strength gap.

## 0. Executive verdict

| Claim (as stated by the authors) | Verdict | One line |
|---|---|---|
| Pass U is a norm-preserving Givens hop, cousin of pass F | **Sound** (algebra) | Verified term-by-term; it is in fact a discrete-time Schrödinger hop |
| Conservation dissolves into a theorem (L1-B) | **Confirmed** | drift_rel 0 at amp_nat 0/2.5/5; pairwise norm-preservation is exact |
| Anti-ignition is structural (L1-C, self-gating) | **Confirmed for the *unarmed* bath** | The armed (ρ_coh≈0.77) case was NOT measured — see §5 |
| In-pass precession "exactly mirrors pass F" | **False on two counts** | Opposite chirality vs F; and τ is phase-dependent (nonlinear) where F is linear |
| L1-A: "the thesis is CONFIRMED … coherent dense translation emerges" | **Overclaimed** | Direction bias is real and a first; the bar (cos→1, speed, seed-robust) is not met |
| Narrow resonance is "a dispersion signature" | **Partly** | Clean dispersion *would* give exactly this — but only if τ were linear; τ here is state-dependent, so it is a nonlinear resonance until proven otherwise |
| arg(ψ) door is "the designed §II.7 retention fix" | **No** | It writes the cell clock `th2[i]`, not a slot-borne atom; IV.6 is unimplemented |
| arg(ψ) door is energy-conserving | **Yes** | Only th2 is written; the Em ledger is untouched. But it is a phase-only merge, not a unitary amplitude merge |
| Early-time "R2d~0.1, ~3× retention" | **Weak** | At t=6 nlive≈21 ⇒ 1/√N≈0.22, so 0.1–0.28 is noise-floor; the matter m-peak is essentially never +2 |
| QUENCH wall = "matter churn" | **Refuted by their own data** | h2 halves births (8733→4580), h3 quadruples them (→33581), R2d flat at 0.02–0.03 in all four arms. Churn is not rate-limiting |

Bottom line: the **engineering move is right and the structural wins (B, C) are
real**, but the two "thesis CONFIRMED" claims (L1-A translation, arg(ψ) door
retention) are one to two notches stronger than the evidence supports, and the
"matter churn" diagnosis for QUENCH is contradicted by the arm-to-arm birth
counts. There are also two genuine code defects (empty-cell phase reset,
one-way evap door) and one large undocumented gap (no C≡Go at amp_nat>0).

---

## 1. Theory / code soundness

### 1.1 Pass U is a discrete-time Schrödinger hop (stronger than "Givens")

The hop (`freecell.c:1561-1565`, identical in Go at `fab/step.go:488-491`) is

```
dm1[i] = cc*m1i + ss*m2j;  dm2[i] = cc*m2i - ss*m1j;
dm1[j] = cc*m1j + ss*m2i;  dm2[j] = cc*m2j - ss*m1i;
```

Writing ψ = m1 + i·m2 and cc=cos τ, ss=sin τ, this is **exactly**

```
ψ_i' = cc·ψ_i − i·ss·ψ_j
ψ_j' = cc·ψ_j − i·ss·ψ_i          ... (A)
```

which is `exp(−i τ (c_i† c_j + c_j† c_i))` applied to (ψ_i, ψ_j) — i.e. the
**one-step discrete-time Schrödinger hop for a tight-binding (adjacency)
Hamiltonian on the link**. This is more than "a Givens rotation": it identifies
pass U as a **quantum walk** on the cell roster, so the entire body of
tight-binding phenomenology (ballistic spreading ∝ t, Bragg reflection, band
revivals) applies by name. That is the right way to read L1-A.

Pairwise norm-preservation is then the statement ⟨ψ|ψ⟩ is invariant under
`e^{−iHτ}`, which is a theorem:

```
|ψ_i'|²+|ψ_j'|² = (cc²+ss²)(|ψ_i|²+|ψ_j|²)
                 + 2·cc·ss·[m1i·m2j − m2i·m1j + m2i·m1j − m1i·m2j] = original.
```

Cross-terms cancel pairwise. **L1-B is structural, confirmed.** I agree with
the authors here.

### 1.2 The link current — derivation, and the Re/Im doc bug

From (A), the energy gained by cell i to leading order in τ is

```
Δ_i = |ψ_i'|² − |ψ_i|² = sin²τ·(|ψ_j|²−|ψ_i|²) + sin(2τ)·Im(ψ_i* ψ_j)
                       = 2τ·Im(ψ_i* ψ_j) + O(τ²).            ... (B)
```

So the **energy current** is `J_s = 2 sin τ · Im(ψ_i* ψ_j)` (exact for this
hop), and the code comment at `freecell.c:1522-1523` and `L1_FINDINGS.md:60`
are correct. The **momentum** is the first moment of this current,
`p ∝ Σ_links (r_j − r_i)·J_s` — the README §II.6 framing is right (this is the
center-of-energy theorem; the monopole telescopes to a boundary flux).

**Bug (documentation, not code):** the README is self-contradictory on the
current. `v93/README.md:176-179` writes `2·Re(ψ_i*·w·ψ_j)`; `v93/README.md:211`
writes `J_s = Im(ψ_i* · w · ψ_j)`; `v92/consult/SUBQUARK_synthesis.md:25` writes
Re again. (B) settles it: **Im, not Re.** The Re form would be the *even*
channel that the additive ledger could already absorb — the whole point of the
unitary channel is that the *odd* (Im) channel is no longer discarded. Fix the
§II.4 line.

### 1.3 Precession ordering — right for the dense clock, but not "exactly F"

The in-pass precession (`freecell.c:1544-1550`, Go `fab/step.go:465-471`) is

```
dm1 = cc·a − ss·b ;  dm2 = ss·a + cc·b   ⇒   ψ_m → e^{+i w2e dt} ψ_m.
```

Pass F's precession (`freecell.c:2057-2062`) is

```
fa1 = cc·a1 + ss·a2 ;  fa2 = −ss·a1 + cc·a2   ⇒   ψ_field → e^{−i w1e dt} ψ_field.
```

**Opposite chirality.** This is not a bug — the dense clock has *always*
advanced by `th2 += w2e·dt` (`freecell.c:2143-2144`, the line face 3 skips), so
`+w2e` is the correct dense-clock sign, self-consistent with the historical
advance. The field happens to use the opposite convention. But the claim that
face 3 "mirrors pass F exactly" (`L1_FINDINGS.md:107-109`,
`freecell.c:1538-1542`) is **false on sign**: it mirrors the *structure*
(precess-then-hop) and the *hop form*, not the chirality. If the programme ever
wants a single complex-structure law across the two sectors (a long-term goal),
this sign needs an explicit decision, not a gloss.

The Trotter logic for face 3 is otherwise sound: applying the on-site rotation
*before* the hops and skipping the out-of-pass advance removes the O(dt)
split-error that face-2 showed as centroid wobble. No double-precession (the
`if (amp_nat==0)` guard at `freecell.c:2143` is correct). Byte-inertness at
amp_nat=0 holds (the whole pass-U block and the pass-6 skip are amp_nat-gated).

### 1.4 Where pass U genuinely diverges from F — τ is state-dependent

This is the most important theory/code gap and it bears directly on §2.

- Pass F hop angle (`freecell.c:2081-2082`): `τ_F = field_J · (sA/Aref)(dref/d) / √(s_i s_j) · dt` — **geometry only, no phase.** Pass F is a *linear* unitary channel with a clean dispersion.
- Pass U hop angle (`freecell.c:1477-1479`): `τ_s = amp_nat · base · √(g_ij g_ji)`, with
  `g_ij = [½(1+cos(bq·th_i − bq·w_i d/C − bp·th_j))]^p_gate`, p_gate=8 (`freecell.c:649-659`, `freecell.c:1380-1381`).

So **τ_s depends on the cell phases th_i, th_j through a steep (p_gate=8)
closure gate.** Pass U is therefore a *nonlinear* (state-dependent) unitary
map, not the linear cousin of pass F the README sells. Two consequences:

1. There is **no clean dispersion relation** for pass U in the sense pass F has
   one. What you measure is a self-consistency resonance of the
   phase-dependent hopping. (See §2.)
2. The gate is *even* in the closure phase (gsym = √(g_ij g_ji)). All directed
   transport therefore comes from Im(ψ_i*ψ_j) inside the hop (the odd channel),
   while the *asymmetric* part of the old want `(g_ij − g_ji)` — the directional
   information the additive ledger *did* carry — is discarded. That may be fine
   (it's the part that broke bookkeeping), but it should be a conscious
   decision, not an unnoticed side-effect of choosing the geometric mean.

### 1.5 arg(ψ) door — energy-safe, but neither unitary-merge nor slot-borne

The door (`freecell.c:2178-2193`, Go `fab/step.go:1301-1317`):

```
aph   = atan2(fa2[i], fa1[i]);
Em_old = Em[i] − (d1 + dsp);           // Em[i] is already post-increment (line 2177)
ro = √Em_old;  rc = √d1;
s1 = ro·cos(th2[i]) + rc·cos(aph);
s2 = ro·sin(th2[i]) + rc·sin(aph);
th2[i] = atan2(s2, s1) mod 2π;
```

**What I verified and approve:**

- **`Em_old` is computed at the right point.** The increment `Em[i] += d1+dsp`
  is at `freecell.c:2177`, *before* the `if (P.amp_door>0)` block, so subtracting
  `(d1+dsp)` back at `:2188` recovers the pre-add value correctly. (The brief
  asked specifically about this; it is right.)
- **Energy-conserving.** The block writes only `th2[i]`. Em/Ee/Es are untouched
  by the door, so `Em+Ee+Es` is invariant by construction. L1-B survives.
- `aph` is read after the real scale `fac=√((Ee−d1)/Ee)` (`:2173-2174`); since
  fac>0, `atan2(fa2·fac, fa1·fac) = atan2(fa2,fa1)`. Phase preserved. OK.

**Why it is *not* the designed §II.7 / IV.6 fix:**

1. **It is a cell-clock write, not a slot-borne atom.** §IV.6
   (`README.md:390-394`) and §II.7 (`README.md:236-239`) are explicit: *"Phase
   is carried slot-borne (protected from delivery churn); the fired atom carries
   arg(ψ), not m·th2."* The code writes `th2[i]` — the same cell clock that
   §IV.6 forbids writing at the door. The forbidden `m·th2` pattern is replaced
   by a *better* cell write (a coherent-arg merge instead of a partial pull),
   but it is still a cell write. The slot-borne infrastructure the kernel
   already has (`sre_/sim_` shadow at `freecell.c:1616-1620`, `slph`/`slgid_` at
   `:1609,:1614`) is unused. **IV.6 is not implemented; it is re-labelled.**

2. **It is a phase-only merge, not a unitary amplitude merge.** A true coherent
   sum would give `Em_new = Em_old + d1 + 2√(Em_old·d1)·cos(aph−th2)`, which
   **breaks** the `Em += d1` ledger (that's the |ΣA|²≠Σ|A|² trap of §IV.4 — and
   the authors correctly refuse it). But refusing it means the interference
   cross-term is **kept for the phase direction and discarded from the energy**.
   So the comment at `freecell.c:253-255` — *"carries the field phase + the
   interference cross-term = the current"* — overstates it: the current is not
   deposited anywhere; only a phase estimate borrows its sign. Honest name:
   *"coherent-arg, scalar-magnitude door."*

3. **The space footprint `dsp` has no phase partner.** The door composes `√d1`
   (field-derivative) but the matter actually gains `d1+dsp` (`:2177`). The
   `dsp` (space→matter, the g1 footprint) is magnitude-only. Consistent with
   "space is a mode," but it means arg(ψ_m) is defined on a *slice* of the new
   matter, not the whole parcel.

4. **`amp_door` is a boolean**, not a knob. It appears only in `>0` tests
   (`freecell.c:2178`, `fab/params.go:468`); the value is never a mix weight.
   Naming implies continuous; it is a switch.

### 1.6 Two genuine code defects

**(a) Empty-cell phase reset — a global phase sink.** At `freecell.c:1576-1582`
the writeback is unconditional:

```
Em[i] = a*a + b*b;                     // 0 for an empty cell
ph = atan2(b, a);                      // atan2(0,0) = 0
if (ph < 0) ph += TWO_PI;
th2[i] = ph;                           // → th2 forced to 0
```

Any cell that briefly held then lost Em (common in a churning cloud) gets its
clock slammed to 0 every step under amp_nat>0. This is an **artificial phase
reference / phase-slip machine** that the additive path does not have, and it
matters directly for QUENCH (every reset is a winding-destroying event). Fix:
`if (a*a+b*b > eps) th2[i] = ph;` (leave the stale clock otherwise).

**(b) The reverse door is not phase-faithful.** `field_inject`
(`freecell.c:922-934`):

```
if (Ee>1e-20) { fac = √((Ee+dE)/Ee);  fa1*=fac; fa2*=fac; }   // real scale: field phase kept
else          { fa1 = √dE·cos(th2[i]); fa2 = √dE·sin(th2[i]); } // vacuum: inject at matter phase
```

Matter→field evaporation (`freecell.c:2218, 2246`) composes the new field
energy as a *real scale* of the existing field amplitude — it does **not** bring
the matter phase `th2` across the door unless the field cell was empty. So the
arg(ψ) door is **structurally one-way**: field→matter is coherent-arg,
matter→field is scalar. In a steady-state churning cloud (births≈deaths≈9000)
any winding imprinted at condensation is erased at the next evaporation. This
alone makes steady-state retention in QUENCH impossible by construction,
independent of all other factors. Either make the reverse door symmetric
(compose `√d2·e^{i th2}` into fa1,fa2 with Ee-ledger renorm), or stop expecting
retention while evap is phase-blind.

### 1.7 P1 meter

`freecell.c:1566-1572` records the pairwise energy displacement `ΔE_j·sd[s]·û`.
This is the exact pairwise-conserved energy moment (matches (B) at leading
order; deviates at O(τ²) after a finite Givens, as it must). No bug.

### 1.8 Minor / hygiene

- `v93/README.md:1-7` still says *"docs-only, NOT opened … awaiting sign-off"*
  while `CLAUDE.md` asserts v93 is the active programme with full
  authorization and the kernel is already edited. The status banner is stale.
- `L1_FINDINGS.md` ends at *"arg(ψ) door is the registered next step"*
  (`:296-299`) while `REVIEW_REQUEST.md:71-77` already reports door results.
  The findings file is behind the brief.
- `quench_argdoor.sh:13` parses `$(NF-3)`, which is the `pkd=m±N` *index*, not
  R2d — see the garbled `runs/quench_argdoor.out`. The suite's stdout cannot be
  trusted; numbers must come from `report/analyze_winding.py` on the `.fcs`
  files (I did this independently in §3).

---

## 2. The L1-A narrow resonance — is it a real dispersion relation?

### 2.1 What a *linear* pass U would predict

If τ were geometry-only (like pass F), pass U on a uniform 1-D chain with
on-site `w2e` and link hop τ is the discrete-time Schrödinger map

```
ψ_n(t+dt) = e^{+i w2e dt} [ cos τ · ψ_n − i sin τ (ψ_{n−1}+ψ_{n+1})/2 … ]   (Trotter-split)
```

whose quasi-energy band is `Ω(k) = w2e − 2τ cos(kd)` (sign per the +w2e
chirality), and whose **group velocity is `v_g(k) = dΩ/dk = 2τ d sin(kd)`.**
A wave packet centered at the seed carrier `kx` drifts at `v_g(kx)`. Three
predictions follow immediately:

1. `v_g` changes sign when `kx·d` crosses π/2 (the Brillouin-zone edge) —
   **alternating-direction windows** as the *effective* kd shifts.
2. As τ (≡ amp_nat) grows, the band widens and folds (discrete-time
   quasi-energy lives in (−π,π]); the carrier crosses zone edges repeatedly —
   **multiple sign reversals vs amp_nat.**
3. Packet diffraction: a localized packet spreads ∝ t (ballistic, the quantum-
   walk signature), so its *centroid* still moves at `v_g` while its *peak
   density* falls and its rms grows.

This is a clean, real dispersion story, and it qualitatively matches the
observed table at `L1_FINDINGS.md:150-166` (sign reversals between amp_nat
1.5/1.75/2.0/2.25, robust +x at 2.0, seed-variant magnitude). **So far so
good — except for §1.4.**

### 2.2 Why I do not yet buy "dispersion" as the *cause*

τ here is **not** the linear tight-binding parameter §2.1 assumes. It is
`τ_s = amp_nat·base·√(g_ij g_ji)` with the p_gate=8 closure gate (`:1477`,
`:1380`, `:649`). Four independent effects can each fake the table in §2.1
without true continuum dispersion:

| Confound | Mechanism | How to rule it in/out |
|---|---|---|
| **(C1) Phase-dependent τ (the real one)** | The gate makes the map a *nonlinear* unitary; "resonance" is a self-consistency window of the phase-locked subgraph, not a band edge | **Linearize τ**: set `gsym ← ⟨gsym⟩` (time-frozen) or `p_gate=0`; if the window persists, it is linear dispersion; if it shifts/vanishes, it is the nonlinearity |
| **(C2) Gauss-Seidel sweep bias** | Hops are applied in-place, in roster order from the i-side (`freecell.c:1551-1555`). A directed in-place sweep on a tilted packet can itself produce a net drift independent of any dispersion | **Reverse the `i` loop** (`for i=NC−1; i>=0; i−−`) at fixed seed; if cos flips sign, ordering is a confound. Cleaner long-term: checkerboard (even/odd) split, the standard symmetric quantum-walk decomposition |
| **(C3) τ saturation** | `if (tau>0.5) tau=0.5` (`:1479`) is a hard nonlinearity; at amp_nat≳few on deep-overlap links the map saturates and amp_nat stops being a linear hopping strength | Histogram per-link τ; report the clamp fraction at each amp_nat. If many links sit at 0.5 near the "resonance edges," the reversals are aliasing of the clamp, not band-folding |
| **(C4) Centroid-fit aliasing** | `blob_drift` fits a straight line to the **second half** of the tag-centroid time series (`freecell.c:4524-4541`) on a packet that the authors themselves show *wobbles* (com_dev 0.06→0.10→0.07, `L1_FINDINGS.md:67-68`) and *melts* 3.6% (`:67`). A sloshing, fragmenting packet aliased into a linear fit can flip vx between seeds and between amp_nat without any bulk translation | Fit the **full** trajectory; report the linear-fit R²; plot COM(t), Em_tag(t), rms(t), peak Em. Require rms bounded and Em_tag retained for "translation" |

None of C1–C4 has been controlled for. The seed-robust *sign* at fixed amp_nat
(5/5 +x at 2.0, 5/5 −x at 1.5) is genuine and rules out *pure* noise — but it
is consistent with all four of C1–C4 (a sweep bias, a clamp edge, a nonlinearity
window, and a fit on a wobbling packet all give seed-robust signs). So:

> **Verdict on §2:** "a real dispersion signature" is *plausible* and would be
> the right answer *if τ were linear*. With the phase-dependent τ actually in
> the code it is **not established** — it is one of four live possibilities, and
> the decisive experiments (C1–C4 above) are cheap and have not been run.

### 2.3 The tag-centroid is the wrong meter for a diffracting packet

`tag_stats` (`freecell.c:2537-2571`) sums only over `tag[i]` cells — the
*seed-labelled* set. Under amp_nat>0 the unitary hops move ψ_m across the
tag/bath boundary freely (untagged bath cells are eligible receivers), so
energy *leaks* the tag silently; the "Em_tag melts 3.6%" is the visible part.
The centroid therefore measures "where the *original* cells' share of the
energy went," not "where the packet went." For a clean dispersion measurement
you want the **energy-weighted COM over all cells**, or — better — the
**Fourier peak of `Em·e^{i th2}`**, which reads the carrier k directly and
gives v_g without any centroid fit at all.

### 2.4 How to widen / pin magnitude (ranked)

1. **Run C1 (linearize τ) first.** This is the single most informative
   experiment in the whole programme: it tells you whether the resonance is
   physics of the *linear* channel (then widen by tuning kd via kx and the
   chart factor) or physics of the *nonlinear* gate (then widen by softening
   p_gate or replacing √(g_ij g_ji) with a smoother phase-matching envelope).
2. **kx sweep at amp_nat=2.** If `v_g(kx)` traces a sin branch, it is tight-
   binding; if it is scatter, it is C2/C3/C4.
3. **dt-halving at fixed `amp_nat·dt`.** True rate physics is invariant;
   Trotter/clamp artifacts move. (The τ>0.5 clamp and the
   precess-then-hop split are both dt-sensitive.)
4. **Uniform clock (`q_detune=0`) on e3b.** The load-detuning
   `w2e=w2/(1+q_detune·x)` (`freecell.c:1311-1314`) shears the on-site
   frequency across the blob (`L1_FINDINGS.md:185-189`); if removing the shear
   sharpens cos and lowers seed-variance, inhomogeneous broadening was
   magnitude-limiting. (Caveat: on QUENCH, q_detune=0 *quadrupled* births with
   no retention gain — see §3 — so it is not a free win.)
5. **Softer odd channel.** gsym is even; if you want directed transport to be
   more robust, a real rotation-angle contribution from an odd closure piece
   (κ_reac-shaped, `freecell.c:1395-1410`, currently zero under V3a) is the
   controlled knob — but it must stay a real angle to preserve norm.

---

## 3. The QUENCH-3 retention wall

### 3.1 I re-ran the analyzer myself (the suite script is broken)

`quench_argdoor.sh:13` extracts the wrong column (§1.8), so I ran
`report/analyze_winding.py` on `runs/quench/h2_door_uni.fcs` directly.
Steady-state (t=20..300), the door+unitary arm:

```
R2d ∈ [0.012, 0.050],  mean ≈ 0.027;
matter m-peak roams: m−3, m+2, m+2, m+0, m−3, m−2, m+1, m+4, m+3, m+0, m−3, m−4, …
                     (essentially never locked to +2);
field RA2: 1.0 → 0.79 (t=20) → 0.05–0.31 (t≥60);
nlive grows 539 → 1870 (the cloud is still condensing, not steady).
```

**Two facts the authors' summary does not surface:**

1. **The matter winding is not merely *small*, it is *not at m=+2*.** R2d is
   `|⟨e^{i(th2 + kx·x − 2φ)}⟩|`, but the spectrum's *peak* is at a random m
   each frame. R2 (no demod) is also ~0.02–0.03. So matter has **no coherent
   winding at any mode** — it is phase noise, not a weak +2 imprint. The meter
   `R2d` is mildly misleading because its noise floor (≈1/√nlive, here
   ≈0.023–0.043 across t=20..300) is *exactly* the reported value.
2. **At t=20 the field still carries the winding (RA2=0.788) but matter has
   none (R2d=0.018).** This is the decisive diagnostic: the field has not yet
   decohered, yet the door has already failed to imprint. **The wall is the
   door/imprint, not field-transit-decoherence** (the latter is real and
   dominant by t≥60, but it is not the t=20 limit). The L1_FINDINGS diagnosis
   ("the field winding decoheres in transit … so by the time matter condenses
   the imprint is already speckle," `:240-242`) is wrong for the early window
   the door was supposed to win.

### 3.2 "Matter churn is the wall" is refuted by the four arms

Births/deaths from `runs/quench_argdoor.out`:

| arm | config | births | R2d (t=20..300) |
|---|---|---|---|
| h0 | no door, amp_nat=0 | 8733 | ~0.025 |
| h1 | arg(ψ) door, amp_nat=0 | 9111 | ~0.025 |
| h2 | arg(ψ) door, amp_nat=2 | **4580** (≈½) | ~0.025 |
| h3 | arg(ψ) door, amp_nat=2, q_detune=0 | **33581** (≈4×) | ~0.025 |

If churn (the door-event rate ∝ births) were rate-limiting, R2d should *rise*
from h3→h0→h2 as births fall. **It is flat.** Halving the churn (h2) buys
nothing; quadrupling it (h3) costs nothing. **Churn is a non-limiting
multiplier, not the wall.** This is the single cleanest experimental fact in
the QUENCH suite and it directly contradicts the authors' conclusion 2
(`REVIEW_REQUEST.md:84-87`).

### 3.3 So what *is* the wall?

Ranked by the evidence above:

1. **The carrier is wrong.** Phase is stored on the cell clock `th2`, which the
   door overwrites in place (`:2193`) and pass U resets to 0 when Em empties
   (§1.6a). A cell clock is not a protected carrier; every overwrite/reset is a
   phase-slip site. This is the §IV.6 design lesson, and the current door does
   not implement it.
2. **The reverse door is phase-blind** (§1.6b). Even a perfect imprint is erased
   on the next evaporation, and in a churning cloud that is within a few beats.
3. **No topological support.** §II.8 (`README.md:250-252`) itself states winding
   is protected *only if |ψ|>0 on a cycle*. The QUENCH "matter" is a sparse
   growing cloud (nlive 539→1870, Em_LIVE=0.05 threshold) riddled with
   sub-threshold holes — phase-slip sites by construction. You are trying to
   hold winding in an object that is topologically incapable of holding it.
4. **Field transit decoherence** (RA2 1→0.15) — real and dominant at t≥60, but
   *not* the t=20 limit (§3.1).
5. **Churn** — present, but shown non-limiting by §3.2.

### 3.4 Experiments that distinguish these (and the one that could yield steady retention)

| Experiment | Predicted outcome if wall is… |  |
|---|---|---|
|  | **carrier/imprint** (3.3.1–2) | **topology** (3.3.3) | **field** (3.3.4) |
| Hand-seed `th2 = 2φ` on a dense ring at t=0, **no field, no door**, amp_nat=2 | R2d should *stay* high if transport retains; if it dies, transport itself loses phase | needs a closed ring first | n/a |
| Slot-borne fired atom (write phase into `sre_/sim_`/`slph`, decay with amp_tau; pass U reads it; **don't** slam th2) | should help — this is the actual §II.7 | neutral | neutral |
| Symmetric reverse door (compose `√d2·e^{i th2}` into fa, Ee-ledger renorm) | should help under churn | neutral | neutral |
| Freeze roster / suppress birth-death (or restrict parcel identity to long-lived cells) | small effect (matches h2) | **decisive**: if R2d rises, topology was it | small |
| Keep RA2 pinned with a continuous spin drive | if R2d still 0, imprint/hold fail | — | if R2d follows RA2, field was it |
| **Seed a closed dense ring** (nv=48/6-style parcel), not a slit cloud, *then* drive | — | **decisive** for topology | — |

The honest path to *real* steady retention is the **closed-ring substrate + a
slot-borne phase atom + a symmetric reverse door**, run with a phase-slip-rate
meter (the pre-registered §II.8 acceptance language). The current slit-cloud +
cell-clock-door combination cannot retain by construction; no amount of churn
reduction will change that.

---

## 4. Packet binding — should L1-A be reframed?

**Agree with the authors' free-packet diagnosis** (`L1_FINDINGS.md:74-79`,
`:192-194`). A pure unitary hop with no confining potential is a discrete
Schrödinger evolution: packets spread ∝ t (ballistic, per §2.1). The additive
model's cap+inflight+head was an *accidental* binder; face 2–3 removed all
three from the hop (`freecell.c:1514-1515` bypasses passes 3–5 entirely when
amp_nat>0). So diffraction of a free packet is the *expected* outcome, not a
surprising failure — and "cos→1 on a bound blob" is the wrong first bar.

Two reframings, ranked:

**A. Reframe L1-A as wave-packet group velocity (recommended, short-term).**
The clean meter is the **Fourier peak of `Em·e^{i th2}`** (carrier k and its
drift), or equivalently the dense dipole of the P1 moment
(`freecell.c:1566-1572`). Drop "bound blob cos→1, seed-robust" as the first
bar. Adopt staged bars:
- **L1-A0:** seed-robust sign of v_g in a window (already met at amp_nat=2.0
  for direction; extend to |v_g| with a Fourier meter).
- **L1-A1:** |v_g| matches the dispersion prediction `2τ d sin(kd)` (after C1
  settles whether τ is linear).
- **L1-A2:** bound translating object (needs a binder).

**B. Add a binding companion, but go to the *nonlinear* binder the substrate
already has.** The obvious candidate is **not** "put head back in τ" (face-1
lesson: head≈0 in the core freezes transport) and **not** a hand-added cohesive
bond. It is **self-trapping via the load-detuning that is already there**: the
on-site frequency `w2e = w2/(1+q_detune·Em/cap)` (`:1311-1314`) is exactly the
`ω(|ψ|²)` nonlinearity of the **discrete nonlinear Schrödinger equation
(DNLS)**, which supports self-trapped discrete breathers above a critical
norm. The narrow "resonance at amp_nat≈2" may be the **edge of a self-trapping
existence region**, not a linear band feature — which would simultaneously
explain (i) the narrowness, (ii) the seed-variance (critical-norm phenomena are
initial-condition-sensitive), and (iii) why a free packet diffracts just
outside it. **Map the existence of a localized steady state (peak Em / rms
stable over T) in the (amp_nat, q_detune) plane.** This is the most
physics-rich experiment available and it reuses machinery that already exists.

**One missing control that should gate any "translation IS the current" claim
(`README.md:216-219`):** `bath=0 → exactly zero motion`
(`L1_FINDINGS.md:69-70`). A current dipole with balanced monopole should be
visible **in vacuum or a thin bath** — if translation requires a dense bath to
push against, the mechanism is bath-drag, not the unitary dense current. Run
e3b at bath=0 and bath_frac∈{0.25,0.5,1} and report v_g vs bath density.

---

## 5. Missed aspects / alternative explanations

1. **L1-C measured the *unarmed* bath; the armed case (ρ_coh≈0.77) was not
   run.** The byte-identical bath stats (`L1_FINDINGS.md:200-214`) show
   self-gating on **random** phases. The design's *named* risk
   (`README.md:258-273`, §II.9) is a **coherent** bath with ρ_coh≈0.77
   (Phase-M's measurement), for which √(g_ij g_ji) is *not* small. The armed
   L1-C — a ρ_coh-seeded or cantus-primed bath at amp_nat=2–5, watching
   births/glow — is the actual anti-ignition test, and it is absent.

2. **No C≡Go comparison exists at any amp_nat>0.** Every abx log in
   `runs/abx_*.log` carries `amp_nat=0`, and there is no Go run with
   amp_nat>0 anywhere in `runs/`. The entire face-3 + QUENCH-argdoor body of
   results is **C-kernel-only**. Design III.3 and §IV.12 (`README.md:325-328`,
   `:429-432`) explicitly warned that compounded rotations would amplify the
   1-ulp C/Go divergence and that tolerance-based abx would be needed at
   strength. That tolerance abx has not been built. Until it is, every face-3
   number is uncorroborated by the Go mirror, and the C≡Go identity claim in
   the status table is only validated at the inert point.

3. **Driven-beam births *drop* at amp_nat>0** (`L1_FINDINGS.md:216-218`:
   8800→3100–4400). The unitary channel reshapes condensation metabolism
   without igniting — fine, but "the channel is inert on structure" is no
   longer strictly true; it actively re-routes the field→matter flow. Document
   this regime; it may matter for ASTRO and the XSEC bars.

4. **The L2 fifth shows the same resonance shape as L1-A** (ggm 0→0.132 at
   amp_nat=2, back to 0.015 at 3; `L1_FINDINGS.md:267-272`). Either one
   coherent mechanism drives both (beautiful) or one amp_nat accident
   underlies both (red flag). A **joint amp_nat sweep of e3b-cos and
   tri-ggm on the same law table** is mandatory before either is claimed.

5. **The early-time "3× retention" is noise-floor.** At t=6, nlive≈21 ⇒
   1/√N≈0.22, so R2d≈0.1–0.28 is statistically indistinguishable from zero
   winding. The honest statement is "the arg(ψ) door does not produce a
   statistically significant matter winding at any time," which is also what
   the steady-state data says. The "door helps retention" claim
   (`REVIEW_REQUEST.md:74-76`, conclusion 2) is not supported.

6. **Alternative explanation for the amp_nat 2→2.5 direction flip**
   (complementing §2.2-C4): the second-half linear fit on a packet that
   *fragments* late (Em_tag melting, rms growing) can flip fitted vx when the
   late-time centroid is dominated by a stray fragment. Without COM(t)
   trajectories this is not excluded; the sign-flip could be a fit artifact on
   a diffracting packet, not a band edge.

---

## 6. Concrete next experiments / code changes, ranked

### P0 — correctness & honesty (before any further physics claims)

1. **Fix the empty-cell `th2` writeback** (`freecell.c:1576-1582`): only update
   `th2` when `a²+b² > eps`. Removes the global phase sink. (Also in Go,
   `fab/step.go:502-510`.)
2. **Fix `quench_argdoor.sh:13`** (parse R2d, not `$(NF-3)`); republish the
   QUENCH table with nlive, the full m-spectrum, and R2d's noise floor
   (1/√nlive) shown alongside.
3. **Update `L1_FINDINGS.md`** with the arg-door results and downgrade
   "thesis CONFIRMED" (`:172-173`) → "mechanism engaged; L1-A bar open."
4. **Fix the design doc:** `J = Im` not `Re` (`README.md:176-179`); STATUS
   banner → active; state explicitly that `amp_door` is a cell-clock write,
   not the slot-borne atom §IV.6 specifies.
5. **Build the tolerance abx and run C≡Go at amp_nat=2** on the e3b and QUENCH
   arms (closes the §5.2 gap).

### P1 — decide what L1-A is and localize the resonance

6. **e3b instrumentation pack** per run: COM(t) over *all* cells, Em_tag(t),
   rms(t), peak Em, linear-fit R², mean τ, τ-clamp fraction, ρ_coh, gsym mean.
7. **C1 — linearize τ** (freeze gsym at its time-mean, or p_gate=0); rerun the
   amp_nat sweep. **This is the single most informative experiment in v93.**
8. **C2 — reverse the pass-U sweep order** (or switch to even/odd checkerboard)
   at fixed seed; check whether cos flips.
9. **kx sweep at amp_nat=2** (v_g(kx) — sin branch vs scatter).
10. **dt-halving at fixed amp_nat·dt** (rate physics vs Trotter/clamp).
11. **bath-density sweep** (bath=0, 0.25, 0.5, 1.0): does v_g survive without a
    bath to push against?
12. **Reframe the bars:** L1-A0 = seed-robust v_g sign + |v_g| from a Fourier
    meter; L1-A1 = matches dispersion; L1-A2 = bound object.

### P2 — retention, the real §II.7

13. **Slot-borne fired-atom phase**: on condensation, write phase into
    `sre_/sim_`/`slph` (decaying with `amp_tau`); pass U reads it; stop
    overwriting `th2` per atom. This is what §IV.6 actually says.
14. **Symmetric reverse door** on evaporation (compose `√d2·e^{i th2}` into fa
    with Ee-ledger renorm).
15. **Hand-seeded `th2=2φ` dense ring, no field, no door, amp_nat=2** —
    isolates *hold* from *imprint*. The decisive QUENCH experiment.
16. **Closed-ring substrate** (nv=48/6-style parcel) instead of the slit cloud,
    + phase-slip-rate meter (the §II.8 acceptance language).

### P3 — binding & fifth

17. **(amp_nat, q_detune) existence map** for a localized steady state (peak Em,
    rms over T) — search for DNLS self-trapping rather than only "add a binder."
18. **Joint L1-A / L2 amp_nat fine sweep** on one law table (§5.4).
19. **Armed L1-C** — ρ_coh-seeded / cantus-primed bath at amp_nat=2–5
    (the actual anti-ignition test, §5.1).

### P4 — cleanup

20. Optional odd (κ_reac-shaped) contribution to τ as a *controlled* experiment.
21. DEC `δJ` refactor only after L1-A0 is clean (design III.4 step 5 is rightly
    last).

---

## 7. Direct answers to the six brief questions

1. **Theory/code soundness.** Pass U is sound and is (more sharply) a
   discrete-time Schrödinger hop `ψ_i' = cc·ψ_i − i·ss·ψ_j`; pairwise
   conservation is a theorem (verified algebraically). In-pass precession order
   is right for the *dense* clock, but the "exactly mirrors pass F" claim is
   false: opposite chirality (+w2e vs −w1e), and τ is phase-dependent
   (nonlinear) where pass F's is geometry-only (linear). The arg(ψ) door is
   energy-conserving (only th2 written; `Em_old` is recovered at the correct
   point, after the increment), but it is a phase-only merge, not a unitary
   amplitude merge, and it writes the cell clock — so §IV.6 is unimplemented.
   No sign error in the hop or in `J = 2τ·Im(ψ_i*ψ_j)`; the README's `Re` at
   §II.4 is wrong. Two real defects: empty-cell `th2→0` reset
   (`freecell.c:1576-1582`), and a one-way evap door (`field_inject:922-934`).
   Go matches C exactly (Sincos arg order, `sFree=0=S_FREE`, precession sign).

2. **Narrow resonance.** A *linear* pass U would give exactly the observed
   alternating-direction windows (Bragg folding of the quasi-energy band
   `Ω = w2e − 2τ cos(kd)`). But τ here is state-dependent via the p_gate=8
   gate, so "dispersion" is one of four live hypotheses (with Gauss-Seidel
   sweep bias, τ-clamp saturation, and centroid-fit aliasing on a wobbling
   packet). Discriminate by linearizing τ (C1), reversing the sweep (C2),
   histogramming τ-clamp (C3), and fitting full COM trajectories (C4). Widen
   via C1's outcome + kx tuning + a possible DNLS self-trapping region.

3. **QUENCH wall.** "Matter churn" is **refuted** by the four arms (h2 halves
   births, h3 quadruples them, R2d flat at 0.02–0.03). The wall is the
   **carrier** (cell clock, overwritten and reset) + **one-way evap door** +
   **no topological support** (sparse cloud, holes = phase-slip sites); field
   transit decoherence is real but only dominant at t≥60 (at t=20 the field
   still has RA2≈0.79 while matter has R2d≈0.02). Decisive experiments: a
   hand-seeded matter-only `th2=2φ` ring (hold vs imprint), and a closed-ring
   substrate. Code path to real retention: slot-borne fired atom + symmetric
   reverse door + closed ring.

4. **Packet binding.** Necessary for a *bound translating blob*; not necessary
   for a first success. Reframe L1-A as wave-packet group velocity (Fourier
   peak of `Em·e^{i th2}` / dense P1 dipole), staged L1-A0/A1/A2. The most
   promising binder is the load-detuning already in the kernel (`w2e(Em)` =
   the DNLS nonlinearity): map the self-trapping existence region in
   (amp_nat, q_detune). And run bath=0 — if v_g needs a bath, it is drag, not
   the dense current.

5. **Missed.** Armed L1-C (ρ_coh≈0.77) never run; no C≡Go at any amp_nat>0;
   driven-beam births *drop* (metabolism reshaping, not inert); L2 and L1-A
   resonances co-located (one mechanism or one accident?); early-time R2d is
   noise-floor (1/√N); evap door is one-way; flight sector fully disabled
   under amp_nat; the matter m-peak is essentially never +2.

6. **Next, ranked.** P0 honesty/correctness (empty-cell reset; script parse;
   doc fixes; C≡Go at strength) → P1 (e3b instrumentation; **linearize-τ
   experiment**; sweep-order reversal; kx & dt sweeps; bath-density; reframe
   bars) → P2 (slot-borne atom; symmetric reverse door; hand-seeded matter
   ring; closed-ring substrate) → P3 (DNLS existence map; joint L1-A/L2 sweep;
   armed L1-C).

---

## 8. Bottom line

The core engineering move — **replace the additive Em want with pairwise
unitary hops** — is sound, byte-inert at amp_nat=0, conservation-theoretic
(L1-B), and self-gating on the incoherent bath (L1-C, unarmed). Face 3's
in-pass precession is a real improvement and produces the first seed-robust
**directional** dense bias, which I read as a quantum-walk group velocity.

It does **not** yet establish "coherent dense translation" as an accepted L1-A
object (direction yes; magnitude and speed seed-variant; binding absent;
resonance origin unresolved between four hypotheses). It does **not** implement
the designed §II.7/§IV.6 retention fix (the door writes the cell clock, not a
slot-borne atom; the reverse door is phase-blind). And the "matter churn is the
QUENCH wall" diagnosis is **contradicted by the authors' own four-arm birth
counts** (R2d flat while births vary 4×–9×); the real wall is carrier +
one-way door + no topological support, with field decoherence joining late.

**Disagreement with the authors' conclusions:**
- (1) "thesis CONFIRMED" — **partial**: structural wins on L1-B/L1-C yes;
  L1-A is *mechanism-engaged, bar open*, and the "dispersion" label is
  unproven while τ is phase-dependent.
- (2) "door helps retention; wall is churn" — **disagree**: door help is
  noise-floor; churn is *refuted* as the wall by the four arms.
- (3) "narrow resonance, seed-variance, no binding are the L1-A gaps" —
  **agree, and add**: the resonance origin is confounded (§2.2), the meter is
  wrong for a diffracting packet (§2.3), and bath-dependence is unmeasured
  (§4/§5.6).

Tighten the meters, run the four C1–C4 discriminants, close the design–code
gap on the door, and build the C≡Go tolerance abx before claiming structural
victory beyond conservation and (unarmed) anti-ignition.
