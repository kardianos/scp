package main

import "strconv"

// buildR1Arms: endurance T default 10000 on compact HOLD winners from T=3000 matrix.
func buildR1Arms(seed int, T float64) []Arm {
	ts := strconv.FormatFloat(T, 'f', -1, 64)
	sd := strconv.Itoa(seed)
	base := func(extra ...string) []string {
		a := []string{
			"L=16", "dt=0.02", "bath=1",
			"T=" + ts, "seed=" + sd,
			"diag_every=2000",
			"hop_order=1",
			"exp=blob", "amp=0.8", "sigma=2.5",
			"amp_nat=2", "amp_logate=1",
			"wf_on=1", "k_rad=0", "noise_amp=0",
		}
		return append(a, extra...)
	}
	return []Arm{
		{Name: "R1_S00", Family: "r1",
			Desc: "endurance clone S00 uni+qd1.2 quiet",
			Args: base("kx=0", "q_detune=1.2", "p1_meter=1"), Snap: true},
		{Name: "R1_S04", Family: "r1",
			Desc: "endurance clone S04 deep qd=3.6",
			Args: base("kx=0", "q_detune=3.6"), Snap: true},
		{Name: "R1_S14", Family: "r1",
			Desc: "endurance clone S14 progressive kx",
			Args: base("kx=1.1", "q_detune=1.2", "p1_meter=1"), Snap: true},
		{Name: "R1_S16", Family: "r1",
			Desc: "endurance clone S16 seed_mw=2",
			Args: base("kx=0", "q_detune=1.2", "seed_mw=2"), Snap: true},
		// refinements (overnight capacity)
		{Name: "R1_S00_rad", Family: "r1",
			Desc: "S00 + mild radiance — does Em* park appear at long T?",
			Args: base("kx=0", "q_detune=1.2", "k_rad=0.05", "p_rad=4"), Snap: true},
		{Name: "R1_S00_s111", Family: "r1",
			Desc: "S00 seed robustness",
			Args: func() []string {
				a := base("kx=0", "q_detune=1.2")
				for i := range a {
					if len(a[i]) > 5 && a[i][:5] == "seed=" {
						a[i] = "seed=111"
					}
				}
				return a
			}(), Snap: true},
		{Name: "R1_S04_rad", Family: "r1",
			Desc: "deep + radiance long",
			Args: base("kx=0", "q_detune=3.6", "k_rad=0.05", "p_rad=4"), Snap: true},
		{Name: "R1_S00_amp12", Family: "r1",
			Desc: "S00 heavier seed amp=1.2 — mass dependence",
			Args: func() []string {
				a := base("kx=0", "q_detune=1.2")
				for i := range a {
					if a[i] == "amp=0.8" {
						a[i] = "amp=1.2"
					}
				}
				return a
			}(), Snap: true},
	}
}
