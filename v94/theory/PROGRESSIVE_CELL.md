# Progressive cell — symbolic treatment and numerical probes

**Status:** theory + initial numerics (2026-08-08).  
**Ontology:** unchanged (no background, Em/Ee/Es books, free roster).  
**Method:** Bohmian/Madelung as *computation* — morph the cell from rigid
harmonic capsule to local sample of **line-like progressive** structure.

Companion notes: `BOHMIAN_HYBRID.md` §10, `SHEAR_GENESIS_DESIGN.md`,
`METABOLIC_MATTER.md`. Results of probes: §6–§8 below; scripts in
`v93/report/analyze_progressive.py`, `v93/runs/prog/`.

---

## 1. Symbolic framework

### 1.1 Two descriptions of the same substance

**Capsule description (default mental model, often wrong for identity):**

\[
\text{cell } i:\quad
(E_m^{(i)},\, \theta_2^{(i)},\, \psi_f^{(i)},\, E_s^{(i)},\, \mathbf{x}_i,\, r_i)
\]

Identity of an object ≈ fixed set of indices \(I_{\mathrm{tag}}\) seeded at \(t=0\).

**Progressive description (method target):**

On the contact graph \(G_t=(V_t,E_t)\) at time \(t\),

\[
\begin{aligned}
\rho_i &= E_m^{(i)} && \text{(dense density sample)} \\
\phi_i &= \theta_2^{(i)} && \text{(dense phase sample)} \\
\psi_i &= \sqrt{\rho_i}\, e^{i\phi_i} && \text{(polar reconstruction)} \\
J_{ij} &= 2\,\mathrm{Im}(\psi_i^{*}\psi_j)
       = 2\sqrt{\rho_i\rho_j}\,\sin(\phi_j-\phi_i)
       && \text{(link current / guidance bilinear)} \\
(\nabla\phi)_{ij} &= \mathrm{atan2}\!\big(\sin\Delta\phi,\cos\Delta\phi\big)/d_{ij}
\end{aligned}
\]

Field sector analog: \(\psi_f = f_{a1}+i f_{a2}\), \(E_e=|\psi_f|^2\), same \(J^f_{ij}\).

**Guidance velocity (discrete Madelung / Bohm method, not a new particle):**

\[
\mathbf{v}_i \;\sim\;
\frac{1}{\rho_i+\varepsilon}
\sum_{j\sim i} J_{ij}\,\hat{\mathbf{e}}_{ij}
\quad\text{(dense)}\,,
\qquad
\mathbf{v}^f_i \;\sim\;
\frac{1}{E_e^{(i)}+\varepsilon}
\sum_{j\sim i} J^f_{ij}\,\hat{\mathbf{e}}_{ij}\,.
\]

Pass U / pass F are discrete realizations of transporting \((\rho,\phi)\)
consistently with these currents (Givens hops preserve pairwise norms).

### 1.2 Progressive identity (what “the object” is)

Define a **progressive object** at time \(t\) as a connected set of samples
with coherent progression, not a tag set:

\[
\begin{aligned}
S_\tau(t) &= \{\, i : \rho_i \ge \tau \,\}
  && \text{(density support)} \\
C_k(t) &= k\text{-th connected component of } S_\tau \text{ under } G_t \\
R_2(C) &= \Bigl|\sum_{i\in C} \rho_i e^{i\phi_i}\Bigr|
         \Big/ \sum_{i\in C}\rho_i
  && \text{(phase coherence of support)} \\
\Phi_\gamma &= \sum_{(ij)\in\gamma} (\nabla\phi)_{ij}\, d_{ij}
  && \text{(path progression; winding if }\gamma\text{ closed)}
\end{aligned}
\]

**Claim (method):** physical identity tracks \((C_k, R_2, \{J_{ij}\})\) in
time, not \(I_{\mathrm{tag}}\). Capsule IDs are mesh labels.

**Tag fidelity** (diagnostic of capsule thinking):

\[
F_{\mathrm{tag}}(t)
  = \frac{\sum_{i\in I_{\mathrm{tag}}} \rho_i}
         {\sum_i \rho_i}\,.
\]

If \(F_{\mathrm{tag}}\) falls while \(\sum\rho\) and \(R_2(C)\) hold, the
capsule description has failed and the progressive description is required.

### 1.3 Line-like harmonics

A **harmonic line** is a 1-chain (path or cycle) on \(G_t\) carrying
consistent phase progression:

\[
H[\gamma] = \bigl(\Phi_\gamma,\; J_\gamma,\; \bar\rho_\gamma\bigr),\quad
J_\gamma = \mathrm{mean}_{(ij)\in\gamma}|J_{ij}|\,.
\]

- **Bound core harmonic:** high \(\bar\rho\), high \(R_2\), closed \(\gamma\)
  with stable \(\Phi_\gamma = 2\pi m\).
- **Outer / shell harmonic:** progressive structure primarily in
  \(\psi_f\) on the \(E_s\) depression around a dense support
  (field lines / phase rings), not extra Em bags.
- **Unbinding:** \(R_2\) or closed \(\Phi_\gamma\) lost while energy may
  remain — delocalization of *progression*, not only of mass.

### 1.4 Continuity equation (what conservation *means* progressively)

Capsule bookkeeping: \(\sum_i \rho_i = \mathrm{const}\) (pass U theorem).

Progressive form (discrete continuity):

\[
\frac{d\rho_i}{dt}
  = -\sum_{j\sim i} \mathcal{F}_{ij} + \mathcal{S}_i^{\mathrm{door}}
\]

where \(\mathcal{F}_{ij}\) is the hop flux implied by the Givens angle
(related to \(J_{ij}\) and \(\tau_s\)), and \(\mathcal{S}_i^{\mathrm{door}}\)
is conversion (field↔dense). **Identity of a progressive object** can move
even when \(\mathcal{S}=0\): pure \(\mathcal{F}\) rearranges support across
cells.

### 1.5 Relation to temperature and shear (kept distinct)

| quantity | progressive reading |
|---|---|
| **Shear** | spatial structure of \(\mathbf{v}\) or \(\nabla\phi\) (organized \(J\)) |
| **Temperature** | excitation of bound progressive stacks (not \(\mathbf{v}\) itself) |
| **Genesis** | door source \(\mathcal{S}>0\) where field progression converges / beats |
| **Melt** | loss of \(R_2\) / spreading of \(S_\tau\) under \(J\) without binding |

### 1.6 What to compute (meter list) — refined after probes

| symbol | name | use | avoid as sole bar |
|---|---|---|---|
| \(F_{\mathrm{tag}}\) | tag fidelity | capsule failure diagnostic | as mass |
| \(R_2(C)\) global | phase alignment of support | **stationary / locked** objects only | traveling waves, DNLS hoards |
| \(\langle|J|\rangle\) | current activity | “something is flowing” | alone: noise also high |
| \(\langle J_x\rangle\) signed | **directional** progressive current | **primary field progressive meter** | — |
| \(\mathrm{Jaccard}(S_t,S_{t+\Delta})\) | support continuity | progressive identity | tau-sensitive |
| \(v_{\mathrm{COM}}^{\mathrm{all}}\) vs tag | pattern vs capsule motion | transport reading | tagged-only COM |
| \(f_d\) (p1_meter) | dense hop current (online) | L1-A progressive transport | blob_drift alone |
| \(n_{\mathrm{comp}},\,E_{m,\mathrm{cmax}}\) | fragmentation / main clump | condensation spectrum | single-component assumption |

**Corrected claim:** \(R_2\) is *not* a universal progressive quality. A pure
phase ramp \(\phi=-kx\) has \(R_2\approx 0\) globally while being *maximally*
progressive. Prefer **path progression** \(\Phi_\gamma\) and **signed \(J\)**.

---

## 2. Numerical programme

### 2.1 Tools

| tool | role |
|---|---|
| `report/analyze_progressive.py` | contact-binned \(J\), \(F_{\mathrm{tag}}\), Jaccard, \(R_2\), COMs |
| inline `dirJ` (in this campaign) | signed field \(J_x\) on pulses |
| existing Face C fcs | long dense condensation / melt |
| `runs/prog/*` new short arms | field progressive + e3b unitary vs additive |

### 2.2 Pre-registered expectations vs outcome

| E | expectation | outcome |
|---|---|---|
| E1 | \(F_{\mathrm{tag}}\downarrow\) while Em fixed under condensation | **CONFIRMED** |
| E2 | melt: PR↑; condense: Em_max↑ PR↓ | **CONFIRMED** |
| E3 | \(R_2\) dies on DNLS hoards | **CONFIRMED** (and sharpened: high \(R_2\) on melt is *misleading*) |
| E4 | organized \(k_x\) ⇒ directed \(J^f\) at equal \(E_e\) | **CONFIRMED** (\(J_x\)) |
| E5 | Jaccard holds better when bound | **CONFIRMED** (deep > law > melt) |

---

## 3. Results — dense progressive structure (Face C reanalysis)

Source: `runs/faceC/{c_law05,c_melt05,c_deep2,c_law05_V3a}.fcs`  
Meters: `runs/prog/*.prog` (τ=0.15, stride 2).  
Common: amp_nat=2, amp_logate=1, Em_tot≈92.97 conserved (except V3a tax).

| arm | F_tag 0→300 | PR 0→300 | Em_max | mean Jaccard | J_mean end | n_comp | R2 end |
|---|---|---|---|---|---|---|---|
| **melt** (qd=0) | 0.89→**0.26** | 456→**1246** | 0.56→0.52 | **0.11** | 0.057 | 30 | 0.79* |
| **law** (qd=1.2) | 0.89→**0.61** | 456→**203** | 0.56→**2.17** | **0.46** | 0.21 | 28 | 0.15 |
| **deep** (qd=3.6,amp=2) | 0.89→**0.84** | 456→250 | 2.25→**7.24** | **0.75** | 0.89 | 18 | 0.03 |
| **law+V3a** | 0.89→0.61 | 456→241 | 0.56→**1.57** | 0.46–0.55 | 0.14 | 24 | 0.12 |

\*Melt’s high end-\(R_2\) is on a **tiny** residual component (size_cmax→2,
Em_cmax→0.7) — accidental / seed-phase leftover, **not** a healthy object.

### 3.1 Interpretation

1. **Capsule identity fails wherever dynamics are live.** Even the deep
   trap loses 16% of tag fidelity; melt loses 71% while total Em is flat.
2. **Progressive support continuity (Jaccard) ranks binding correctly:**
   deep (0.75) > law (0.46) > melt (0.11). This is the right “does the
   object persist?” meter under cell turnover.
3. **Link current \(J\) lights up under condensation** (0.21–0.89) and
   stays low under pure melt (0.06) — progressive activity, not just mass.
4. **Global \(R_2\) is anti-useful for DNLS hoards:** it collapses precisely
   when energy self-traps (Face C: binds energy, not phase). Do **not**
   use \(R_2\) as particle quality here.
5. **V3a tax** pulls Em_max to the radiance fixed point (~1.57) without
   restoring \(R_2\) or tags — tax is progressive-mass filter, not a
   phase machine.

### 3.2 Symbolic lesson

\[
\text{identity}_{\mathrm{capsule}}(I_{\mathrm{tag}})
\;\not=\;
\text{identity}_{\mathrm{progressive}}(S_\tau, J, \text{COM}_{\mathrm{all}}).
\]

Conservation \(\sum\rho=\mathrm{const}\) is necessary but silent about
which description is true. Jaccard + COM_all + directed \(J\) speak.

---

## 4. Results — field progressive structure (new runs)

Config: `exp=pulse bath=1 L=16 T=24 amp=2 σ=3 wf_on=1 k_rad=0`, seed fixed.
**Equal \(E_e\):** pulse_kx and pulse_k0 both \(\sum E_e = 548.27\) at all \(t\).
Noise control is **not** energy-matched (\(\sum E_e\approx 2745\)) — qualitative only.

| arm | organization | \(J_x\) (signed, t=0→24) | \(|J|\) | global \(R_2^f\) |
|---|---|---|---|---|
| **pulse_kx** (k=1.2) | phase wind | **−0.57 → −0.27** | 0.91→0.44 | **≈0.00** |
| **pulse_k0** (k=0) | standing Gaussian, φ=0 | **≈0** | 0→0.08 | **1.00→0.76** |
| noise_amp=2 | random | ≈0 (±0.02) | ~1.5 | ≈0.01 |

### 4.1 Interpretation

1. **Directional current is the progressive signature.** At fixed energy,
   \(k_x\) produces large signed \(J_x\); \(k=0\) does not. Absolute \(|J|\)
   alone is weaker (noise can be large).
2. **Global field \(R_2^f\) is inverted for progressive waves.** Traveling
   phase wind has \(R_2^f\approx 0\); sterile standing packet has
   \(R_2^f\approx 1\). Using \(R_2\) as “coherence quality” would rank the
   **wrong** arm first.
3. **Energy density is not structure.** Same \(E_e\), radically different
   \((J_x, R_2^f)\). Shear/genesis designs must match energy *and* report
   directed \(J\).
4. Pulse arms birthed **no Em** on this short dark-ish door window
   (cond≈0 in logs for k0; tiny field credits only). Progressive field
   structure **exists without genesis** — method meters work before
   SHEAR-G birth bars.

### 4.2 Symbolic field progressive order parameter

Prefer

\[
\mathcal{P}_x
  = \frac{\sum_{\langle ij\rangle} J^f_{ij}\,(\hat e_{ij}\!\cdot\!\hat x)\,w_{ij}}
         {\sum w_{ij}},
\quad w_{ij}=\sqrt{E_e^{(i)}E_e^{(j)}}
\]

over global \(R_2^f\). Path progression \(\Phi_\gamma=\int_\gamma\nabla\phi\cdot d\ell\)
along high-\(E_e\) ridges is the continuum limit of the same idea.

---

## 5. Results — dense transport progressive (e3b short)

Config: blob amp=0.5 σ=2.5 kx=1.1 bath=1 T=40 p1_meter=1 hop_order=1.

| | unitary linear (amp_nat=2) | additive (amp_nat=0) |
|---|---|---|
| **fd** (dense hop current) | **+285.5** | 0 |
| **fl** (inflight) | 0 (pass U bypass) | **−145** |
| F_tag 0→40 | 0.89→**0.66** | 0.89→**0.81** |
| PR | 456→**154** | 456→281 |
| Em_max | 0.56→**2.82** | 0.56→0.63 |
| mean Jaccard | **0.45** | **0.69** |
| Em_tot | ~93 (small evap) | 93→**62** (in-flight ledger + dynamics) |
| dcom (all vs tag) end | 0.43 | 0.33 |
| blob_drift cos | 0.19 (misleading) | 0.43 (also weak) |

### 5.1 Interpretation

1. **Unitary progressive transport is real:** large positive `fd` with
   seed-matched theory direction — confirms item 2 in progressive
   language (current is the beable of transport, not the tag COM).
2. **Unitary condenses while transporting** (PR↓, Em_max↑, Jaccard lower):
   progressive identity is a **moving multi-component pattern**, not one
   blob of fixed cells.
3. **Additive keeps tags better** but moves mass into flight (`fl`) and
   loses cell Em_tot — different progressive book (parcel flight vs hop
   current). Comparing Em_tag between laws is **invalid**.
4. **blob_drift remains a bad progressive meter** (low cos even when fd
   is huge) — do not revive it.

---

## 6. What works (adopt)

| approach | evidence | action |
|---|---|---|
| **Morph cell → sample of progressive amplitude** | F_tag fails; Jaccard ranks binding | analysis default; docs language |
| **Link current \(J\) / signed \(J_x\)** | pulse_kx vs k0; e3b fd | primary progressive / shear meter |
| **Support Jaccard + COM_all** | Face C ranking; e3b | identity & motion |
| **Pass F/U as discrete guidance** | conservation + currents | conceptual + teaching |
| **Multi-component census** (n_comp, Em_cmax) | Face C hierarchy | never assume one object |
| **Energy-matched controls for field structure** | pulse Ee identical | SHEAR-G mandatory |

---

## 7. What to avoid (retire or demote)

| approach | why it fails | replacement |
|---|---|---|
| **Tag sets as object mass** | F_tag 0.89→0.26 with Em fixed (melt) | Em_tot / support Em |
| **Global \(R_2\) as quality** | high for sterile standing field; ~0 for progressive wind; dies on good DNLS hoards | signed \(J\), path \(\Phi_\gamma\), local coherence |
| **\(|J|\) alone** | noise has large |J| with \(J_x\approx 0\) | directed / structured \(J\) |
| **blob_drift / tagged COM for L1-A** | cos low while fd huge | `p1_meter` fd / COM_all |
| **Single Em_max story** | n_comp→20–30 under live unitary | spectrum + components |
| **Equating high \(R_2\) melt residue with success** | size_cmax→2 | require Em_cmax & PR jointly |
| **Tracers / extra corpuscles** (for this morph) | not the intended methodology | progressive meters on existing ψ |
| **Literal trajectory ontology** | fights PRINCIPLE; not needed for method | guidance of *amplitude* |

---

## 8. Updated symbolic picture (post-numerics)

```
                    free-cell mesh (computational support)
                              │ samples
                              ▼
              ψ_i = √ρ_i e^{iφ_i} ,  ψ^f_i = f_a1 + i f_a2
                              │
              links: J_ij = 2 Im(ψ_i* ψ_j)     ← progressive harmonics live HERE
                              │
              objects: connected supports S_τ of ρ or E_e
                       with continuity Jaccard(S_t, S_{t+Δ})
                       and motion from COM(ρ) or integrated J
                              │
              doors: S_i source/sink between modes (genesis / radiance)
```

**Cell morph (operational definition):**

> A cell is a **chart element** that stores samples \((\rho,\phi,\psi_f,E_s)\).
> A **harmonic** is a consistent progression of phase/current on a 1-chain
> of contacts. A **material object** is a time-continuous family of
> supports \(\{S(t)\}\) with nontrivial progressive structure — not a
> fixed index set.

Bohmian methodology enters only as: compute as if guiding
\((\rho,\phi)\) with \(v\sim j/\rho\); do not add corpuscles.

---

## 9. Implications for adjacent programmes

| programme | implication |
|---|---|
| **SHEAR-G** | Pass bar must use **signed/structured \(J^f\)** + energy match; global \(R_2^f\) is a **negative** control (high ⇒ standing, not progressive) |
| **Face C / hoards** | Report Jaccard + n_comp + Em spectrum; stop quoting \(R_2\) as quality |
| **L1-A** | Already fixed by fd; progressive doc makes the reason structural |
| **Metabolic T** | Still not shear; progressive stacks can exist at T≈0 (pulse_kx has structure, no Em thermometer) |
| **Shells** | First instrument: progressive \(\psi_f\) structure around dense \(S_\tau\) (field \(J^f\), path phase on the Es well) — offline on existing dumps before new knobs |

---

## 10. Minimal next experiments (if continuing)

1. **Energy-matched phase scramble** of the same pulse amplitude (true S0 for SHEAR-G) — noise_amp is not matched.
2. **Local path-\(R_2\)** on sliding windows along \(+\hat x\) for pulse_kx (should be high locally while global \(R_2\) is low).
3. **Shell probe offline:** for deep hoard frames, map \(E_e\) and \(J^f\) in radial shells about COM_all — no kernel change.
4. Do **not** implement tracer particles for this morph; if anything, store richer **link-native** diagnostics in fcs (optional).

---

## 11. Files

| path | content |
|---|---|
| `theory/PROGRESSIVE_CELL.md` | this document |
| `report/analyze_progressive.py` | progressive meters |
| `runs/prog/*.prog` | Face C + pulse + e3b summaries |
| `runs/prog/pulse_{kx,k0}.fcs` | energy-matched field progressive pair |
| `runs/prog/e3b_{lin,add}.fcs` | transport progressive pair |

---

*End of progressive-cell treatment. Method confirmed numerically: link-
directed currents and support continuity work; capsule tags and global
R₂ systematically mislead.*
