# FREECELL_LOG — running record of the free-cell substance campaign

**Opened 2026-08-02.** Session worklog for the attempt on `FREECELL.md` §5
(the substance half: can a free-cell fabric hold a localised object
together). Subordinate to `PRINCIPLE.md`; standing table `laws_V2g.cfg`;
ratchet rule governs. **`cellfab.c` is not modified** — all work in a new
standalone instrument, `freecell.c`, per the LIVEFAB precedent.

Method rules of `FREECELL.md` §6 and the DESIGN standing rule apply:
*"Any future boundary proposal should state its percolation fraction and
its signal-to-noise against σ_d BEFORE it is coded."* §2 below is that
statement. §3 states every prediction before the first run.

---

## 1. The mechanism, found in the standing law

Reading `cellfab.c` for the port (full extraction in session notes):
**positions are never read by the law.** The law consumes, per link, only
the length `ld` (retardation `d/C`, conductance `dref/d`, flight rate
`C/d`), the direction `û` (collimation), and the live contact area `lA`
(from live radii `cr = cr0·(Es/e_s0)^{1/3}`). A free-cell port therefore
changes *where `d`, `û`, `lA` come from* — live geometry instead of a
frozen scaffold — and touches no law.

Two mechanisms already in the standing kernel become **geometric cohesion**
the moment positions are dynamical:

**(A) The rung bond** — P15 retardation plasticity (`cellfab.c:3068-3103`)
is the S2 odd term on the link's geometric conjugate:

    dd = -kappa · base · Sm · dV/dd,   V = 1 - G(psi_f)·G(psi_b)
    psi_f = q·th_i - q·w_i·d/C - p·th_j,   psi_b = p·th_j - p·w_j·d/C - q·th_i
    psi_f + psi_b = -(q·w_i + p·w_j)·d/C

so V is minimised exactly on the pair-separation ladder
**(q·w_i + p·w_j)·d/C = 2π·m** — the measured CONSONANCE law (E6/E7
tongue). On the frozen foam this walks the link *parameter* `ld` (the
anneal probe: d 1.4999 → 1.5420 to a joint (θ,ω,d) lock, gates → 1.0000,
within 10 t.u. at κ=1). On free cells the same force, applied ±û/2 to the
two endpoint *positions*, is a **two-sided restoring bond at quantized
lengths** d*_m = 2πmC/(q·w_i + p·w_j). Vacuum is exactly inert:
Sm = sqrt(mi_eff·mj_eff) and mj_eff = sqrt(mob_j·mir) = 0 when Em_j = 0.

**(B) Space-deficit re-jamming** — loading consumes space
(`Em += d1+dsp; Es -= dsp` at condensation, s_pull=0.15), radii shrink,
and on a live substrate the bath re-jams inward (LIVEFAB run 2, measured:
far-field σ_d 8.4% → 4.7% live vs → 14.5% frozen). The bath *confines* a
shrunken loaded cluster instead of frustrating around it. This also
predicts a cure for the D3 self-strangulation death (frozen: r shrinks at
fixed d until d > death radius 1.63 and the channel pinches; free: the
bath pushes the pair together as it shrinks, d falls with r).

Two cohesion **regimes** follow from the law's own numbers:

| regime | load | interior gates | held by |
|---|---|---|---|
| **molecule** (pair/ring/truss) | x ≤ ~0.9, on rungs | open (locked) | rung bonds (A) |
| **droplet** (e3a blob) | x ≈ 0.95 | closed (tail mismatch — the e3a seal) | cage + (B) + gate seal; bond force ~0 because G·G' ≈ 0 |

The heavy blob's m=1 rung d* = π(1+1.2·0.95)/2.9 = 2.32 lies *outside*
contact range (1.7) — but its interior gates are closed (ψ ≈ −w·d/C ≈ −3
rad), so the bond force, ∝ G'·G, is negligible in both directions. The
molecule regime is where bonds act; the droplet regime is where the
seal + cage + re-jamming act. Both must be measured, with controls.

## 2. Percolation fraction and signal-to-noise (stated before coding)

* **Percolation.** The bond exists on a link iff BOTH endpoints carry
  dense load (Sm > 0 needs mob_i, mob_j > 0). Inside a seeded object the
  loaded fraction is 1 by construction — the bond graph IS the object's
  own link set, not a matching condition. This is not harmfrust's
  situation: the rung condition is one condition per LINK on the link's
  own free coordinate d (positions supply 3 DOF/cell; a truss with
  E ≤ 3V−6 internal links satisfies all rung lengths generically), not d
  simultaneous commensurability conditions per cell on one shared scalar.
  harmfrust's 5–13% is what the VACUUM shows, and there it is a feature:
  the bath's bond fraction is ~0 (vacuum exactly inert), far below
  p_c ≈ 1/(z−1), so no accidental aggregation.
* **Signal-to-noise vs disorder.** Bond walk rate near a rung (laws_V2g,
  pair at x=0.325, d≈1.55): |dd/dt| = κ·(kd·geo·gpl·res)·Sm·|dV/dd|
  ≈ 2.4·0.107·1·1 · 0.81 · 1.49 ≈ **0.31 length/t.u.** (κ=1), against a
  geometric relaxation velocity ~0.07–0.1/t.u. at typical overlap — the
  bond dominates its own link by ~3–4× and equilibrates ~0.01 off-rung
  (see P-FC1). Crumb-contaminated vacuum (Em ~ 1e-3 both sides):
  Sm ≈ 3.2e-3 → bond velocity ≤ ~1.7e-3/t.u., i.e. ≤ 2% of the repulsion
  scale — the bath is safe from bond-driven densification (the LIVEFAB
  run-1 catastrophe was an unconditional spring; this force carries Sm).

## 3. Pre-registered predictions (before any run)

* **P-FC1 (pair bond).** Two voices at x=0.325 (Em=0.8125), phases seeded
  at lock, d0 off-rung by ±0.05: d converges to d* = π(1+1.2x)/2.9 =
  1.506 within ~5–20 t.u. (κ=1), from BOTH sides, and equilibrates at
  d* + δ_eq with δ_eq ≈ v_rep/K_b ≈ 0.01 (bond stiffness K_b ≈ 7/t.u.
  per unit length, linearised at the rung; repulsion at 8% overlap).
  Unloaded control pair: force exactly 0 (Sm=0); drifts only under
  repulsion. Detuned control (one voice x=0, one x=0.65): no common rung
  in reach → no lock, no walk to 1.506.
* **P-FC2 (truss).** Ring6 (E=6 < 3V−6=12: floppy) holds its BOND LENGTHS
  but deforms freely transverse; octahedron (E=12 = 3V−6: isostatic)
  resists shear up to bond strength and recovers. Both, in open space,
  retain Em to ~1e-3 over 1000 t.u. (no vacuum contact → no bleed; only
  internal roughness, which is ~0 on-comb). In-bath versions bleed at the
  measured universal c₀ ~ 4e-4 Em/t.u./voice (E-regression), i.e.
  load-line lifetimes, and the SKIRT death (x → x_skirt = 0.0617) is the
  expected end state — free cells are not claimed to fix the balance
  curve (KIMI §3: no zero crossing; that is a different, open problem).
* **P-FC3 (blob at density, the §5a test).** At cellfab-standard density
  (φ_nom ≈ 2.0, z_live ≈ 14): e3a-class blob (amp 1.6, x_core ≈ 0.95)
  holds geometrically (centroid speed ≤ 0.002, the e3a bar) and retains
  energy on the load line; κ_bond=0 control ALSO holds there (cage) —
  the discriminator is the **dilute run**: bath expanded until z_live < 4
  (below isostatic): tuned/bonded object still holds; κ=0 object
  disperses. Sub-jamming dispersal must actually run (method rule 6).
* **P-FC4 (footprint & strangulation cure).** The loaded blob maintains
  a graded Es depression (g1 analog; core/far ≤ 0.8 under s_k=0.06,
  s_disp=0.3) and the LIVE bath re-jams into it: blob-region live degree
  stays ≥ 7 (F2 floor) where the frozen scaffold loses contacts
  (strangulation). Two-blob COM attraction is expected BELOW the noise
  floor (g2 precedent: sub-floor even on the frozen foam) — recorded,
  not gated.
* **P-FC5 (conservation & determinism).** Total ledger (ΣEs+Em+Ee+Σlem,
  Kahan) drifts ≤ 5e-14 relative through birth/death events (rule α:
  a channel dies only at e_mid = 0; dying channels flush their cycles).
  Same-seed reruns byte-identical; dt and dt/2 agree to first order.
* **P-FC6 (bath sanity).** Vacuum bath: jam-settles to a static packing
  (max|Δx| → 0), z_live ≥ 7 everywhere, σ_d comparable to S1's band, NO
  densification (mean degree bounded ≤ ~20), and with a leaky object
  present the crumb-driven vacuum contraction stays ≤ 2% of d̄ (the
  anneal-probe surprise, now bounded by the Sm scaling above).

Failure of P-FC1 kills the molecule regime; failure of P-FC3's dilute
discriminator plus P-FC1 kills the whole substance claim and the honest
next step is FREECELL §5(b)/(c) shape work. Either way the §5a test has
then actually been run at density with the configuration printed.

## 4. The instrument — freecell.c

Standalone, `v89/freecell.c`, build `gcc -O2 -fopenmp`, laws read from
`battery/laws_V2g.cfg` VERBATIM (same key=value parser semantics), plus
freecell-only apparatus keys (geometry sector), echoed in the header:

* geometry: `mob_geo` (overdamped position mobility), `k_rep` (contact
  repulsion), `cfac=1.15` (candidate rule), `jam_sweeps` (init settle),
  walls clamped as in cellfab's anneal;
* the bond: `kappa_bond` (=kappa_plast semantics, applied ±û/2 to the
  endpoint positions; ldd buffered Jacobi, applied in pass D);
* channels: persistent identity table keyed by the pair; birth free at
  (e_mid=0, φ=0); death ONLY at Σlem=0 after leaving the candidate set
  (rule α); dying channels advance flight and deliver (flush); a dying
  channel whose receiver is full on completion returns the residue to
  its source (the (β) fallback, counted and reported);
* law passes ported verbatim from cellfab.c order: 0 flload, S space,
  1 radii/pitch, 2 wants+bond (comb/gates/mutual mobility/heads), 3
  outflow, 4 deposit+entrain, 5 flight/delivery/rough(+atoms), D apply
  bond+repulsion to positions, F unitary field hops (symmetric
  normalization, canonical link order), 6 beat conversions (+atoms,
  s_pull), 7 alignment+tumble (serial gaussian stream);
* diagnostics per row: t, rel. drift, φ_nom, z_live, NL alive/dying,
  births/deaths, d̄, σ_d, max|F_net|, max|Δd_bond|, plus per-experiment
  columns; every claim cites a log line.

Known deviations from cellfab (documented, justified): (i) hop order is
canonical link order, not color batches (same physics class, different
integrator ordering); (ii) dying channels advance flight where cellfab
holds pinched channels (rule α requires the flush); (iii) FP64 ledger
(cellfabf-class reference), no integer ledger in v0.

## 5. Status

- [x] Docs absorbed; law extracted code-faithfully (session notes).
- [x] §2 percolation + S/N stated; §3 predictions registered.
- [x] Symbolic derivations (deriv_bond.py; charge/runs/deriv_bond_results.log).
- [x] freecell.c v0 built (periodic box, kick caps, NaN autopsy).
- [x] FC-1 pair bond + 6 controls.
- [x] FC-2 truss suite (rings, octahedron, yield ladder, parity).
- [x] FC-2b bath-embedded ring6 (THE substance result).
- [~] FC-3 blob campaign + FC-2c death watch — 8 runs in flight.
- [ ] Tier-A field pulse; writeup into FREECELL.md; commit.

## 6. MEASURED RESULTS (all logs charge/runs/freecell_*.log)

### 6.1 Derivations (charge/runs/deriv_bond_results.log, sympy + honest 2-cell integration)

* Bond stiffness EXACT: K_b = kappa·(kd·geo·gpl·res)·Sm·p_gate·S²/(4C²);
  at the FC-1 point K_b = 10.998 /t.u. (FD-checked 1e-10). On the
  symmetric branch the kernel's fixed-phase dVdd = (p·S/2C)·h^(2p−1)·
  sin(S·delta/2C) — restoring both sides of the rung, everywhere in the
  basin. ψ_f+ψ_b = −S·d/C is PHASE-INDEPENDENT, so locking needs no
  entrainment (verified: honest integration hits d* to 1e-15 without it).
* Capture: gate-only force peak |δ|=0.17, but the lens area dies at
  contact — hard zero at δ=+0.1319. Bond and repulsion share support.
* Equilibrium: δ_eq(mob·k_rep=1) = 0.0121 first-order / 0.0247 honest
  branch; exists up to mob·k_rep < 3.91 — huge margin.
* Vacuum exactly inert (Sm=0 if either Em=0); crumb pairs cannot bind
  (repulsion wins 5.5–35× at the crumb rung; gate collapses at 10%
  overlap) — no vacuum crystallisation, DERIVED.
* μ-criterion δ/r > (α−1)/2 confirmed symbolically (per-cell derivative;
  the both-cells-grow path gives (α−1) — a factor-2 convention caveat
  for geopress §3.5 readers).
* Rigidity: ring6 floppy=6, cube floppy=6, octahedron/icosahedron
  floppy=0 (isostatic). ring6+3 diagonals: floppy=4 with one self-stress.

### 6.2 FC-1 pair (freecell_pair_*.log) — THE BOND IS REAL

* up/dn (d0 = d*±0.05): both converge to d=1.5387 (identical to 6
  digits), phase offset dl=+0.116 rad ⇒ δ_live ≈ 0.028 vs derived 0.0247
  (breathing ±0.003 about it). Gates 0.99. drift 4e-16. The pair runs a
  closure current (~0.115 Em/voice in flight) and sits at x = seeded x.
* kappa_bond=0 control: pair slides OUT to exact contact d=1.6377 =
  2r(Es) — the STRANGULATION death radius of the frozen kernel; the bond
  is what rescues strangulation on free cells.
* unloaded control: force exactly 0, drift exactly 0.0 (vacuum inert).
* detuned control (x 0.325 vs 0.65): no comb rung in reach — nothing
  moves, no lock, no coupling.
* dt/2: d_final differs 2e-4 (first-order convergent).
* seedlock=0 (random phases): does NOT acquire — repulsion wins the race
  before locking; FREE PAIRS NEED CONFINEMENT (or a seeded lock) TO
  ACQUIRE. Dense crowding catalyses bonding.

### 6.3 FC-2 truss (freecell_ring*, freecell_oct*.log)

* COLLIMATION THEOREM (measured, then explained): bond strength carries
  gpl; each cell has ONE carrying axis (the planes' intersection), so no
  normal assignment serves a 3D deltahedron — the law's natural dense
  objects are CHAINS/RINGS/TUBES. Radial-normal ring6 (gpl=0.56) held
  35 t.u. then escaped via a growing breathing oscillation; octahedron
  (gpl 0.25/0.0625 mixed) distorted, weak edges strangled.
* Coplanar ring6 (n1=n2=ẑ, gpl=1): STABLE ≥300 t.u., edge_dev 0.035,
  gg 0.98, x̄ flat, drift 3e-16. First standing free-cell object.
* PARITY: ring5 DEAD (odd cycle: n·w·d*/C = nπ needs even n; π-frustrated
  locks never form, slides to contact) — ring4 stable, ring6 stable,
  ring8 metastable (breathing escape between t=70 and 270). Species
  rule: even small rings live; odd rings are not species; large simple
  rings are marginal.
* TRIANGLE/π theorem (analysis): at the m=1 rung every 3-cycle sums to
  3π ≢ 0 — unison objects must be bipartite, so no triangles, so no
  generic 3D rigidity in the unison sector. Hexagon chord bracing is
  length-consonant (2d* = m=2 rung exactly) but PHASE-frustrated (chain
  forces Δθ=π across 3 edges, chord wants 0). 3D rigid species require
  MIXED-INTERVAL triangles — the fifth-triangle {3:2,2:3,1:1} of the
  charge sector (prestress theory's odd-cycle closure). Rigidity ⟺ odd
  cycles ⟺ interval content: the ℤ₃/flavor connection, now motivated
  from cohesion.
* Yield ladder (shear_eps at t=80): 0.02, 0.05 recover (bond lengths
  restore; floppy shape memory persists — plastic shape, elastic
  lengths); 0.08, 0.12 break to contact. Elastic limit between 5% and
  8%, derived ψ*/(w·d*) ≈ 0.079. k_rep 1→0.3 scales δ_eq 0.035→0.0093
  (linear in pressure, as derived).

### 6.4 FC-2b bath-embedded ring6 (freecell_ring6_bath.log) — HEADLINE

Ring6 carved into a jammed periodic bath (φ_nom 1.62, bath z≈15):
crushed at insertion (edge_dev 0.105, shape crumpled), bonds RECOVER and
relock by t≈12 (edge_dev 0.038 ≈ vacuum value, gg 0.977), conn 1.000
throughout, z_tag → 14.8 (the object re-embeds; F2 floor holds through
the object), drift EXACTLY 0.0 through ~2200 channel births + deaths
(rule α exercised). Bleed: ret 1.000→0.9781 over 40 t.u. ⇒
**dEm/dt ≈ 4.5e-4 per voice — the frozen-foam universal c₀ = 4.25e-4
within ~6%**. The skirt law is substrate-independent. Death clock
predicts t ≈ cap·(x−x_skirt)/c₀ ≈ 1400 t.u.; free-cell prediction: the
dying ring SHRINKS (d* tracks x down, −23% at the skirt) while staying
gate-locked — running (freecell_ring6_bath_death.log).

### 6.5 FC-3 blob campaign (freecell_blob_*.log) — the §5a answer

* L=24 live (φ_nom 1.67, z_live 16 printed per row): HOLDS. conn 1.000,
  COM speed 1.6e-3 (e3a bar 0.002), retention 0.848@t72; frozen-scaffold
  control 0.833 — LIVE ≥ FROZEN. g1-analog live: es core/far 0.61 ≤ 0.8
  (frozen keeps deeper 0.55 — bath re-jams into the dip). z_tag 13.7–16.
* Half density (φ_nom 0.84): holds (conn 1.000, ret 0.855@t160).
* Dilute (φ_nom 0.43): FRAGMENTS — conn 0.008, z_tag 0.00; energy parked
  on disconnected cells (ret 1.000, no channels, drift 4e-15). Dispersal
  on free cells = connectivity loss, not energy spreading. Droplets need
  the cage; bonded objects do not (ring6 holds in empty space; ring6 in
  the HALF-density bath holds too: dev 0.045, freecell_ring6_bathhalf.log).

### 6.6 The death law (freecell_ring6_bath_death.log, T=2000)

Three phases: (1) LOCKED SHRINK t=50–650: gg ≥ 0.974 (peak 0.998), x̄
0.314→0.198 at c₀-class rate 4.8e-4/voice/t.u., edge length tracks
d*(x̄) within 0.011 — the ring shrinks along the rung ladder in tune.
(2) UNLOCK at x̄ ≈ 0.19 — ABOVE the skirt: cage-tracking failure (the
bond cannot displace the jammed bath to follow the ladder inward), gg
0.98→0.34 in ~100 t.u. Usable life 700 vs naive skirt clock 1400.
(3) HUSK: bleeds to vacuum (x̄ 0.02 by t=1550), conn 1.000 throughout —
death is loss of harmonic identity, not scattering.
ring3 dead by parity like ring5 (freecell_ring3_cop.log, dev 0.145).
SOFT-BATH discriminator (bath_frac 0.5, freecell_ring6_softdeath.log):
never cleanly locks (gg 0.32–0.68, loose-cage breathing) but bleeds 7×
slower (6.8e-5/voice vs 5.0e-4 — leak ∝ vacuum-contact surface, z_bath
4.7 vs 15), projected life ~5000 t.u. Lock-quality/lifetime trade; the
new lever on the no-particle bound is the object's contact surface.

### 6.7 Tier-A + auxiliaries

* Pulse (e2 apparatus): v/C = 0.574 on the live substrate (frozen 0.605,
  bar 0.3), drift exactly 0.0 through ~7500 births + 6400 deaths under
  the packet (freecell_pulse.log).
* Pair in bath (e7-style): walks to its LIVE rung and locks at gates
  1.0000, offset +0.0087 (confinement compresses the bond); bleeds at
  c₀ (x 0.325→0.305 over 100 t.u.) (freecell_pair_bath.log).
* ring8 in bath: confined mush — deformed, partially locked (gg 0.5–0.77),
  conn 1.000; confinement prevents dispersal but does not restore the
  lock (freecell_ring8_bath.log).
* Rule α measured at scale: thousands of births/deaths per run with
  drift at 0/FP-floor; β-returns fire only on dying channels (counted,
  usually 0 — flush by delivery dominates).
* Flux-locked-topology cohesion (the 'immortal channel' idea) is
  RESOLVED BY DESIGN: rule α's flush means in-flight energy delivers or
  returns within one cycle of candidacy loss — channels cannot pin
  geometry; cohesion lives in the bond, not in channel immortality.

### 6.8 e3b free-cell panel — pending fill

amp=0.5 kx=2.08 across SEED_PANEL {20260802, 111, 222222, 314159,
7777777} + frozen-scaffold control (freecell_e3b_*.log). Frozen-foam
standing result: passes 1/5 (speed ≥ 0.003 AND cos ≥ 0.8).

### 6.9 Engineering record

* Periodic box required: walls + motion create a crust (candidate
  explosion). Dart throw must be wrap-aware or seam-close pairs explode
  the jam (NaN cascade — found by autopsy).
* Truss cells must not enter the jam before their positions exist
  (init-order NaN).
* Kick caps (|Δd_bond| ≤ 0.05/step, |Δx| ≤ 0.1/step) = integrator
  guards for deep-overlap encounters (geo ~ 1/d); all measured regimes
  sit far below them.
* Bath: z_live ≈ 15–20 settling, σ_d ≈ 14–17% (dart+quick-jam; S1-class
  annealing would sharpen — improvement path, not blocking).
