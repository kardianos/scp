package main

import "strconv"

// Arm matrix: shear genesis + progressive / Bohmian-method probes.
// Ontology-safe: only freecell knobs. Genesis success = cond>0 with Em
// birth; progressive success = directed J / fd without requiring birth.

type Arm struct {
	Name   string
	Family string // "shear" | "prog" | "ctrl"
	Desc   string
	Args   []string
	Snap   bool
}

func common(seed int, t float64) []string {
	return []string{
		"L=16", "dt=0.02", "bath=1",
		"wf_on=1", "wf_floor=0.01", "wf_far=99",
		"k_rad=0",
		"diag_every=200",
		"seed=" + strconv.Itoa(seed),
		"T=" + strconv.FormatFloat(t, 'f', -1, 64),
	}
}

func buildArms(seed int) []Arm {
	base := common(seed, 24)
	long := common(seed, 60)
	short := common(seed, 16)

	return []Arm{
		// ---- controls ----
		{
			Name: "G0_dark", Family: "ctrl",
			Desc: "dark bath — door floor",
			Args: append(append([]string{}, base...), "exp=bath"),
		},
		{
			Name: "S0_pulse_k0", Family: "ctrl",
			Desc: "standing field pulse equal-Ee baseline (no shear)",
			Args: append(append([]string{}, base...), "exp=pulse", "amp=2", "sigma=3", "kx=0"),
			Snap: true,
		},
		{
			Name: "S0b_noise", Family: "ctrl",
			Desc: "isotropic field noise (structureless energy)",
			Args: append(append([]string{}, base...), "exp=bath", "noise_amp=2"),
			Snap: true,
		},

		// ---- shear / organized field ----
		{
			Name: "A1_kx06", Family: "shear",
			Desc: "mild phase wind",
			Args: append(append([]string{}, base...), "exp=pulse", "amp=2", "sigma=3", "kx=0.6"),
			Snap: true,
		},
		{
			Name: "A2_kx12", Family: "shear",
			Desc: "phase wind kx=1.2",
			Args: append(append([]string{}, base...), "exp=pulse", "amp=2", "sigma=3", "kx=1.2"),
			Snap: true,
		},
		{
			Name: "A3_kx20", Family: "shear",
			Desc: "steep phase wind",
			Args: append(append([]string{}, base...), "exp=pulse", "amp=2", "sigma=3", "kx=2.0"),
			Snap: true,
		},
		{
			Name: "A4_spin2", Family: "shear",
			Desc: "azimuthal field winding spin_m=2",
			Args: append(append([]string{}, base...), "exp=pulse", "amp=2", "sigma=3", "kx=0", "spin_m=2"),
			Snap: true,
		},
		{
			Name: "A5_imp_kx", Family: "shear",
			Desc: "radial imp_k + phase wind",
			Args: append(append([]string{}, base...), "exp=pulse", "amp=2", "sigma=3", "kx=1.2", "imp_k=0.5"),
			Snap: true,
		},
		{
			Name: "A6_hiamp_kx", Family: "shear",
			Desc: "high amp + wind — door threshold hunt",
			Args: append(append([]string{}, base...), "exp=pulse", "amp=6", "sigma=3", "kx=1.2"),
			Snap: true,
		},
		{
			Name: "A7_hiamp_k0", Family: "ctrl",
			Desc: "high amp standing — energy control for A6",
			Args: append(append([]string{}, base...), "exp=pulse", "amp=6", "sigma=3", "kx=0"),
			Snap: true,
		},
		{
			Name: "A8_noise_plus_kx", Family: "shear",
			Desc: "noise + organized pulse wind (lit + shear)",
			Args: append(append([]string{}, base...), "exp=pulse", "amp=2", "sigma=3", "kx=1.2", "noise_amp=1"),
			Snap: true,
		},
		{
			Name: "A9_noise_only", Family: "ctrl",
			Desc: "noise alone vs A8",
			Args: append(append([]string{}, base...), "exp=bath", "noise_amp=1"),
			Snap: true,
		},
		{
			Name: "A10_long_kx", Family: "shear",
			Desc: "longer wind — more door opportunity",
			Args: append(append([]string{}, long...), "exp=pulse", "amp=3", "sigma=3.5", "kx=1.2"),
			Snap: true,
		},
		// A11 bath=0 pulse is INVALID (z_live=0) — dropped after probe_20260808.
		// Door-side genesis hunts (workfn may block pure-field birth):
		{
			Name: "C1_nowf_kx", Family: "shear",
			Desc: "pulse wind with wf_on=0 (door not presence-gated)",
			Args: []string{"exp=pulse", "amp=3", "sigma=3", "kx=1.2",
				"L=16", "dt=0.02", "bath=1", "wf_on=0", "k_rad=0",
				"diag_every=200", "seed=" + strconv.Itoa(seed), "T=40"},
			Snap: true,
		},
		{
			Name: "C2_nowf_k0", Family: "ctrl",
			Desc: "standing pulse wf_on=0 energy control for C1",
			Args: []string{"exp=pulse", "amp=3", "sigma=3", "kx=0",
				"L=16", "dt=0.02", "bath=1", "wf_on=0", "k_rad=0",
				"diag_every=200", "seed=" + strconv.Itoa(seed), "T=40"},
			Snap: true,
		},
		{
			Name: "C3_blob_field", Family: "shear",
			Desc: "dense blob + ambient noise (matter presence for workfn door)",
			Args: append(append([]string{}, base...), "exp=blob", "amp=0.5", "sigma=2.5", "kx=0",
				"noise_amp=2", "wf_on=1"),
			Snap: true,
		},
		{
			Name: "C4_blob_quiet", Family: "ctrl",
			Desc: "dense blob no field noise control for C3",
			Args: append(append([]string{}, base...), "exp=blob", "amp=0.5", "sigma=2.5", "kx=0",
				"noise_amp=0", "wf_on=1"),
			Snap: true,
		},
		{
			Name: "C5_uni_noise", Family: "shear",
			Desc: "unitary dense + noise field (organized matter, lit medium)",
			Args: append(append([]string{}, short...), "exp=blob", "amp=0.5", "sigma=2.5", "kx=1.1",
				"amp_nat=2", "amp_logate=1", "hop_order=1", "noise_amp=1.5", "p1_meter=1"),
			Snap: true,
		},
		{
			Name: "A12_slit_beam", Family: "shear",
			Desc: "directed slit field beam",
			Args: append(append([]string{}, base...), "exp=slit", "amp=2", "sigma=2", "slit_sy=4",
				"kx=1.2", "slit_srcx=4", "slit_mask=0"),
			Snap: true,
		},
		{
			Name: "A13_spin_kx", Family: "shear",
			Desc: "spin + wind combined progressive field",
			Args: append(append([]string{}, base...), "exp=pulse", "amp=2", "sigma=3", "kx=1.2", "spin_m=1"),
			Snap: true,
		},

		// ---- progressive dense / Bohmian-method transport ----
		{
			Name: "B1_e3b_uni", Family: "prog",
			Desc: "unitary linear dense hop — guidance transport",
			Args: append(append([]string{}, short...), "exp=blob", "amp=0.5", "sigma=2.5", "kx=1.1",
				"amp_nat=2", "amp_logate=1", "hop_order=1", "p1_meter=1"),
			Snap: true,
		},
		{
			Name: "B2_e3b_add", Family: "prog",
			Desc: "additive want control",
			Args: append(append([]string{}, short...), "exp=blob", "amp=0.5", "sigma=2.5", "kx=1.1",
				"amp_nat=0", "p1_meter=1"),
			Snap: true,
		},
		{
			Name: "B3_e3b_gated", Family: "prog",
			Desc: "unitary gated τ",
			Args: append(append([]string{}, short...), "exp=blob", "amp=0.5", "sigma=2.5", "kx=1.1",
				"amp_nat=2", "amp_logate=0", "hop_order=1", "p1_meter=1"),
			Snap: true,
		},
		{
			Name: "B4_uni_melt", Family: "prog",
			Desc: "unitary + q_detune=0 melt",
			Args: append(append([]string{}, short...), "exp=blob", "amp=0.5", "sigma=2.5", "kx=0",
				"amp_nat=2", "amp_logate=1", "hop_order=1", "q_detune=0"),
			Snap: true,
		},
		{
			Name: "B5_uni_deep", Family: "prog",
			Desc: "unitary + deep q_detune trap",
			Args: append(append([]string{}, short...), "exp=blob", "amp=1.0", "sigma=2.5", "kx=0",
				"amp_nat=2", "amp_logate=1", "hop_order=1", "q_detune=3.6"),
			Snap: true,
		},
		{
			Name: "B6_ring_mw", Family: "prog",
			Desc: "vacuum ring matter winding",
			Args: []string{"exp=ring", "ring_n=6", "bath=0", "L=16", "dt=0.02", "T=40",
				"seed_mw=2", "amp_nat=2", "amp_logate=1", "hop_order=1",
				"k_rad=0", "wf_on=0", "seed=" + strconv.Itoa(seed), "diag_every=200"},
			Snap: true,
		},
	}
}
