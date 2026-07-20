# Recipe Checklist — When Multi-Fabric Is Ready

**Owner:** Agent A  
**Use when:** multi-fab binary available (CPU `scp_sim` or CUDA with `n_fabrics=3` + B2) and GPU/local budget allocated.  
**Package:** [`c12_atom_package_v0.md`](c12_atom_package_v0.md)

Do **not** invent new kernel physics. Prefer frozen B2 package.

---

## 0. Preflight

- [ ] Read `v75/STATE.md` freeze (q, B2, P/N, grid).  
- [ ] Confirm **no** unauthorized `scp_sim` / `sfa.h` edits needed for this run matrix.  
- [ ] Binary: multi-fab B2 path (`n_fabrics=3`, `mf_lock_CQ=0`, `init_sfa_Q`, `init_sfa_L`).  
- [ ] Tools present: `bin/gen_pn_core`, `bin/radial_qball` (or known profiles), `v75/analysis/score_pn_park.py`.  
- [ ] Profiles on disk (typical): light nuc ω≈1.42–1.46 gauged; L soft ω≈1.46.  
  - Nuclear light: e.g. `profiles/f_w146_g005.txt` / campaign standard used in F19.  
  - Confirm column convention: per-component amplitude (v72).  
- [ ] Data root: `/space/scp/v78/atom/` (or `/space/scp/v75/pn/…` for replay).  
- [ ] Runner: `sim_setup` → `sim_build` → `sim_run`; **never** `sleep` inside `sim_exec`.  
- [ ] Status poll: `v75/analysis/remote_status_snippet.sh` (instant cat/ps only).

---

## 1. Config freeze (copy verbatim unless retuning)

```text
n_fabrics=3
mf_lock_CQ=0
mf_stage=2
q_C=0
q_Q=1
q_L=-1
complex_phi=1
complex_gauge=1
g_gauge=0.05
m_theta=1.6
eta=0
bc_type=0
damp_width=5.0
damp_rate=0.01
init=sfa
```

**Z6 grid freeze (F19):** `N=192`, `L=48`, `dt_factor=0.025`, T=400 smoke / T≥2000 long-T.

**Do not** put Z6 L@R=22 in an L=32 box (sponge eats L) — see `v75/PN_GRID.md`.

---

## 2. Seed generation (`gen_pn_core`)

### 2.1 Nuclear-only twins (park baseline)

- [ ] **Z6N0 nuc:** 6 proton centers (octa R≈8); nN=0; L empty or zero file.  
- [ ] **Z6N6 nuc:** 6 p (octa R≈8) + 6 n (octa R≈5.5); L zero.  
- [ ] Write `*_C.sfa`, `*_Q.sfa` (protons only on Q), `*_L.sfa` (zero).  
- [ ] Verify: `sfa_finalize_header` path in tool (gen_pn_core does this); SFA types POSITION/ANGLE/VELOCITY.

Example shape (geometry from F19 — re-derive centers if profile radius changes):

```text
# Z=6 N=6 + optional L later
gen_pn_core 192 48 <profNuc> <omegaNuc> out_C.sfa out_Q.sfa out_L.sfa \
  6  <6 proton x y z> \
  6  <6 neutron x y z> \
  0
```

### 2.2 Atom seeds (L count = Z)

- [ ] **Z6N0 + L6:** same nuclear as N0; 6 L on shell R=22, ω_L=1.46, `profL`.  
- [ ] **Z6N6 + L6:** same nuclear as N6; **identical** L package (isotope control).  
- [ ] Confirm seed rule: **same sign** of Noether ω on C and L (fabric q flips EM sign).  
- [ ] Prefer **Z6N6+L6** as primary long-T candidate (better nuclear park).

Reference configs: `v75/cfg/pn/z6_n0_nuc.cfg`, `z6_n6_nuc.cfg`, `z6_a_n0.cfg`, `z6_a_n6.cfg`.

### 2.3 Hydrogenoid / orbit (P1)

- [ ] Single p (C+Q) + single L; D≈20; tangential vt ~0.05–0.06 from F16 bracket.  
- [ ] Optional: `gen_mf_shell_orbit` / `gen_mf_pair_boost` for B1-style soft packages if B2 single-center not required.  
- [ ] For P/N bookkeeping on H1, still use B2 + gen_pn_core Z=1 N=0 L=1.

---

## 3. Run matrix (minimum)

| Priority | ID | Content | T | Gate |
|----------|-----|---------|---|------|
| 1 | z6_n6_nuc | Z6N6 rest core | 400 | PASS_nuc |
| 2 | z6_a_n6 | Z6N6 + L6 R22 | 400 | sector L + park |
| 3 | z6_n0_nuc | Z6N0 rest | 400 | retune if c_Q>0.15 |
| 4 | z6_a_n0 | Z6N0 + same L6 | 400 | isotope L identical |
| 5 | long_a_n6 | best atom package | ≥2000 | ideal bar |
| 6 | h1_orbit | single C+L multi-rev | ≥2000 | P1.2 |

Campaign pattern: nohup loop of cfgs (see `v75/cfg/pn/run_z6_remote.sh`); poll status; download diags/tracks first, full SFA if morphology needed.

---

## 4. Diagnostics (every run)

- [ ] `diag.tsv`: Q_phi, E, s_max, gauss_max time series.  
- [ ] Track massL / Ql (SFA columns or track tool).  
- [ ] Score:  
  ```bash
  python3 v75/analysis/score_pn_park.py <diag> <track> --Z 6 --N 0|6
  ```  
- [ ] PASS_nuc / PASS_atom JSON archived under `v78/work/A/results/` or `/space/scp/…`.  
- [ ] Gauss floor check: ~1e-13; any secular rise = **bug tripwire**, stop.  
- [ ] Isotope pair: assert nuclear Q_flux end equal (N0 vs N6); L tracks equal for atom pair.

---

## 5. Visual package (required for P2.4 / P3.4 claim)

- [ ] Headless volview: charge view (9), field, gauge (|E|/|A|), optional flavor.  
  ```bash
  bin/volview -snapshot N -view 9 -out f.webp file.sfa
  ```  
- [ ] Sheet: multi-time composite (see `images/C6_atom_sheet.png` / F19 renders).  
- [ ] Note morphology: multi-ball → droplet by t∼400 is **expected** under current binding; document, do not fail solely on multi-center loss.

---

## 6. Pass / fail decision

| Outcome | Action |
|---------|--------|
| PASS_nuc + c_L≤0.15 + Gauss floor + no merge | Mark atom smoke PASS; proceed long-T |
| Soft park only (c_Q ~0.15–0.20) but morphology + L hold | Document soft; retune spacing/θ before ideal claim |
| L mass collapse or bag merge | FAIL hard; check grid sponge, R_L, private bags |
| Gauss drift | Stop; treat as implementation bug |
| Long-T PASS | Claim P2.4 / partial P3; hand U for CONVERGENCE |
| Long-T FAIL after retune | Document; escalate kernel auth only if isolation/Coulomb path exhausted |

---

## 7. Archive

- [ ] Keep diag.tsv, track.tsv, score JSON, runlog locally under v78.  
- [ ] Large `.sfa`: `/space/scp/`; optional `rclone` archive after analysis.  
- [ ] Append result summary to `logs/A_atom.log`.  
- [ ] Update readiness one-liner in package if status changes.

---

## 8. Explicit non-goals on this checklist

- Full B2 ε_CQ unlock unless blocked  
- Option C triple Cosserat  
- Cosmology / gravity  
- L3 n-bar chemistry populations  
- Fake isotopes by changing only total Q without Z-preserving N knob  

---

## 9. Quick command sketch (remote)

```text
# After sim_setup(remote) + sim_build multi-fab CUDA
sim_upload seeds + cfgs
sim_run id=z6_a_n6 config=<z6_a_n6.cfg content>
# poll sim_run_status / remote_status_snippet — no sleep in sim_exec
sim_download diag + track + optional sfa → /space/scp/v78/atom/
# score locally
```

Reference F19 wall: ~48 min per T=400 N=192 multi-fab on V100-16GB.
