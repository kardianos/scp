package main

import "strconv"

// buildMonoArms: monodispersity / rad fixed-point scan on quiet unitary blob.
// Varies k_rad × q_detune; short horizon (caller sets T, default 2000).
func buildMonoArms(seed int, T float64) []Arm {
	ts := strconv.FormatFloat(T, 'f', -1, 64)
	sd := strconv.Itoa(seed)
	base := func(extra ...string) []string {
		a := []string{
			"L=16", "dt=0.02", "bath=1",
			"T=" + ts, "seed=" + sd,
			"diag_every=1000",
			"hop_order=1",
			"exp=blob", "amp=0.8", "sigma=2.5",
			"amp_nat=2", "amp_logate=1",
			"wf_on=1", "noise_amp=0", "kx=0",
			"p_rad=4",
		}
		return append(a, extra...)
	}
	// 2×4 = 8 arms: law trap vs deep × k_rad ladder
	var arms []Arm
	for _, qd := range []struct {
		tag string
		v   string
	}{
		{"qd12", "1.2"},
		{"qd36", "3.6"},
	} {
		for _, kr := range []struct {
			tag string
			v   string
		}{
			{"kr00", "0"},
			{"kr02", "0.02"},
			{"kr05", "0.05"},
			{"kr10", "0.10"},
		} {
			name := "M_" + qd.tag + "_" + kr.tag
			arms = append(arms, Arm{
				Name:   name,
				Family: "mono",
				Desc:   "mono scan q_detune=" + qd.v + " k_rad=" + kr.v,
				Args:   base("q_detune="+qd.v, "k_rad="+kr.v),
				Snap:   true,
			})
		}
	}
	return arms
}
