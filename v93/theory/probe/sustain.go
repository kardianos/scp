package main

import "strconv"

// buildSustainArms: long-horizon isolation matrix for SUSTAIN (T default 3000).
// Each arm turns one physics switch; shared blob seed so birth is not confounded.
// Genesis long arms (G*) use pulse + wf_on=0 from the C1/C2 line.

func buildSustainArms(seed int, T float64) []Arm {
	ts := strconv.FormatFloat(T, 'f', -1, 64)
	sd := strconv.Itoa(seed)
	// diag every ~20 t.u.; sparse fcs every ~100 t.u. for Em(t)
	// snap_every is in steps: dt=0.02 → 100 t.u. = 5000 steps
	base := func(extra ...string) []string {
		a := []string{
			"L=16", "dt=0.02", "bath=1",
			"T=" + ts, "seed=" + sd,
			"diag_every=1000", // 20 t.u.
			"hop_order=1",
		}
		return append(a, extra...)
	}
	blob := func(extra ...string) []string {
		return base(append([]string{
			"exp=blob", "amp=0.8", "sigma=2.5", "kx=0",
		}, extra...)...)
	}

	return []Arm{
		// ========== dense sustain isolation (pre-seeded Em) ==========
		// Baseline unitary self-trap, door presence-gated, no field, no tax
		{Name: "S00_uni_base", Family: "sust",
			Desc: "unitary linear + law qd; quiet; baseline hold",
			Args: blob("amp_nat=2", "amp_logate=1", "q_detune=1.2",
				"wf_on=1", "k_rad=0", "noise_amp=0", "p1_meter=1"),
			Snap: true},

		// Transport algebra
		{Name: "S01_additive", Family: "sust",
			Desc: "additive want (no unitary)",
			Args: blob("amp_nat=0", "q_detune=1.2", "wf_on=1", "k_rad=0", "noise_amp=0", "p1_meter=1"),
			Snap: true},
		{Name: "S02_uni_gated", Family: "sust",
			Desc: "unitary gated τ (amp_logate=0)",
			Args: blob("amp_nat=2", "amp_logate=0", "q_detune=1.2", "wf_on=1", "k_rad=0", "noise_amp=0", "p1_meter=1"),
			Snap: true},

		// Binding nonlinearity
		{Name: "S03_qd0_melt", Family: "sust",
			Desc: "q_detune=0 — no self-trap",
			Args: blob("amp_nat=2", "amp_logate=1", "q_detune=0", "wf_on=1", "k_rad=0", "noise_amp=0"),
			Snap: true},
		{Name: "S04_qd36_deep", Family: "sust",
			Desc: "deep self-trap qd=3.6",
			Args: blob("amp_nat=2", "amp_logate=1", "q_detune=3.6", "wf_on=1", "k_rad=0", "noise_amp=0"),
			Snap: true},

		// Radiance tax (outtake)
		{Name: "S05_rad", Family: "sust",
			Desc: "V3a radiance tax k_rad=0.05",
			Args: blob("amp_nat=2", "amp_logate=1", "q_detune=1.2",
				"wf_on=1", "k_rad=0.05", "p_rad=4", "noise_amp=0"),
			Snap: true},
		{Name: "S06_deep_rad", Family: "sust",
			Desc: "deep trap + radiance",
			Args: blob("amp_nat=2", "amp_logate=1", "q_detune=3.6",
				"wf_on=1", "k_rad=0.05", "p_rad=4", "noise_amp=0"),
			Snap: true},

		// Field intake (random) — sustain under continuous field energy
		{Name: "S07_noise", Family: "sust",
			Desc: "isotropic field noise intake wf_on=1",
			Args: blob("amp_nat=2", "amp_logate=1", "q_detune=1.2",
				"wf_on=1", "k_rad=0", "noise_amp=1.5", "p1_meter=1"),
			Snap: true},
		{Name: "S08_noise_rad", Family: "sust",
			Desc: "noise intake + radiance tax",
			Args: blob("amp_nat=2", "amp_logate=1", "q_detune=1.2",
				"wf_on=1", "k_rad=0.05", "p_rad=4", "noise_amp=1.5"),
			Snap: true},
		{Name: "S09_noise_deep", Family: "sust",
			Desc: "noise + deep trap",
			Args: blob("amp_nat=2", "amp_logate=1", "q_detune=3.6",
				"wf_on=1", "k_rad=0", "noise_amp=1.5"),
			Snap: true},

		// Open door (wf_on=0) ± field
		{Name: "S10_nowf_quiet", Family: "sust",
			Desc: "wf_on=0 no field — open door, no intake fuel",
			Args: blob("amp_nat=2", "amp_logate=1", "q_detune=1.2",
				"wf_on=0", "k_rad=0", "noise_amp=0"),
			Snap: true},
		{Name: "S11_nowf_noise", Family: "sust",
			Desc: "wf_on=0 + noise — open door with fuel (C5-class)",
			Args: blob("amp_nat=2", "amp_logate=1", "q_detune=1.2",
				"wf_on=0", "k_rad=0", "noise_amp=1.5", "p1_meter=1"),
			Snap: true},
		{Name: "S12_nowf_noise_rad", Family: "sust",
			Desc: "wf_on=0 + noise + radiance — metabolic balance probe",
			Args: blob("amp_nat=2", "amp_logate=1", "q_detune=1.2",
				"wf_on=0", "k_rad=0.05", "p_rad=4", "noise_amp=1.5"),
			Snap: true},
		{Name: "S13_nowf_deep_noise", Family: "sust",
			Desc: "wf_on=0 + deep + noise",
			Args: blob("amp_nat=2", "amp_logate=1", "q_detune=3.6",
				"wf_on=0", "k_rad=0", "noise_amp=1.5"),
			Snap: true},

		// Progressive dense phase structure (Bohmian-method: phase progression on matter)
		{Name: "S14_kx_tilt", Family: "sust",
			Desc: "dense phase wind kx on blob (progressive matter)",
			Args: blob("amp_nat=2", "amp_logate=1", "q_detune=1.2",
				"wf_on=1", "k_rad=0", "noise_amp=0", "kx=1.1", "p1_meter=1"),
			Snap: true},
		{Name: "S15_kx_noise", Family: "sust",
			Desc: "dense kx + field noise",
			Args: blob("amp_nat=2", "amp_logate=1", "q_detune=1.2",
				"wf_on=1", "k_rad=0", "noise_amp=1.5", "kx=1.1", "p1_meter=1"),
			Snap: true},
		{Name: "S16_seed_mw", Family: "sust",
			Desc: "matter winding seed_mw=2 on blob",
			Args: blob("amp_nat=2", "amp_logate=1", "q_detune=1.2",
				"wf_on=1", "k_rad=0", "noise_amp=0", "seed_mw=2"),
			Snap: true},

		// Additive + noise (control for unitary metabolic)
		{Name: "S17_add_noise", Family: "sust",
			Desc: "additive + noise",
			Args: blob("amp_nat=0", "q_detune=1.2", "wf_on=1", "k_rad=0", "noise_amp=1.5", "p1_meter=1"),
			Snap: true},

		// Full V3a-like: rad + workfn + unitary + optional noise
		{Name: "S18_v3a_quiet", Family: "sust",
			Desc: "uni + rad + workfn quiet (adopted law, no fuel)",
			Args: blob("amp_nat=2", "amp_logate=1", "q_detune=1.2",
				"wf_on=1", "k_rad=0.05", "p_rad=4", "noise_amp=0"),
			Snap: true},
		{Name: "S19_v3a_noise", Family: "sust",
			Desc: "uni + rad + workfn + noise (adopted law + fuel)",
			Args: blob("amp_nat=2", "amp_logate=1", "q_detune=1.2",
				"wf_on=1", "k_rad=0.05", "p_rad=4", "noise_amp=1.5", "p1_meter=1"),
			Snap: true},

		// ========== genesis long (C1/C2 line at T) ==========
		{Name: "G01_nowf_kx", Family: "genesis",
			Desc: "pulse wind wf_on=0 long — C1 extended",
			Args: base("exp=pulse", "amp=3", "sigma=3", "kx=1.2", "wf_on=0", "k_rad=0"),
			Snap: true},
		{Name: "G02_nowf_k0", Family: "genesis",
			Desc: "pulse standing wf_on=0 long — C2 extended",
			Args: base("exp=pulse", "amp=3", "sigma=3", "kx=0", "wf_on=0", "k_rad=0"),
			Snap: true},
		{Name: "G03_wf_kx", Family: "genesis",
			Desc: "pulse wind wf_on=1 long — expect no birth",
			Args: base("exp=pulse", "amp=3", "sigma=3", "kx=1.2", "wf_on=1", "k_rad=0"),
			Snap: true},
		{Name: "G04_nowf_kx_rad", Family: "genesis",
			Desc: "pulse wind wf_on=0 + radiance — growth vs tax",
			Args: base("exp=pulse", "amp=3", "sigma=3", "kx=1.2", "wf_on=0", "k_rad=0.05", "p_rad=4"),
			Snap: true},
		{Name: "G05_nowf_k0_rad", Family: "genesis",
			Desc: "standing + rad control for G04",
			Args: base("exp=pulse", "amp=3", "sigma=3", "kx=0", "wf_on=0", "k_rad=0.05", "p_rad=4"),
			Snap: true},

		// ========== ring progressive (topological) ==========
		{Name: "R01_ring_uni", Family: "sust",
			Desc: "vacuum C6 ring seed_mw unitary long",
			Args: []string{"exp=ring", "ring_n=6", "bath=0", "L=16", "dt=0.02",
				"T=" + ts, "seed=" + sd, "diag_every=1000",
				"seed_mw=2", "amp_nat=2", "amp_logate=1", "hop_order=1",
				"k_rad=0", "wf_on=0"},
			Snap: true},
		{Name: "R02_ring_bath", Family: "sust",
			Desc: "ring in bath — contact dephasing long",
			Args: []string{"exp=ring", "ring_n=6", "bath=1", "L=16", "dt=0.02",
				"T=" + ts, "seed=" + sd, "diag_every=1000",
				"seed_mw=2", "amp_nat=2", "amp_logate=1", "hop_order=1",
				"k_rad=0", "wf_on=0"},
			Snap: true},

		// Second seed for key sustain pairs
		{Name: "S00b_uni_base_s111", Family: "sust",
			Desc: "S00 seed robustness",
			Args: func() []string {
				a := blob("amp_nat=2", "amp_logate=1", "q_detune=1.2", "wf_on=1", "k_rad=0", "noise_amp=0")
				// replace seed
				for i := range a {
					if len(a[i]) > 5 && a[i][:5] == "seed=" {
						a[i] = "seed=111"
					}
				}
				return a
			}(),
			Snap: true},
		{Name: "S11b_nowf_noise_s111", Family: "sust",
			Desc: "S11 seed robustness",
			Args: func() []string {
				a := blob("amp_nat=2", "amp_logate=1", "q_detune=1.2", "wf_on=0", "k_rad=0", "noise_amp=1.5", "p1_meter=1")
				for i := range a {
					if len(a[i]) > 5 && a[i][:5] == "seed=" {
						a[i] = "seed=111"
					}
				}
				return a
			}(),
			Snap: true},
	}
}
