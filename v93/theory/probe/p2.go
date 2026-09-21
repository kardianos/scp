package main

import "strconv"

// buildP2Arms: dual-intake thermometer test (qgraph SIM_PREDICTIONS P2).
//
// Claim: with rad_dual=1, hoard load x* tracks bath xb with slope ~1.
// Incumbent (rad_dual=0): x* is bath-independent (slope ~0).
//
// Apparatus: quiet unitary blob + k_rad=0.05 p_rad=4 + bath_xb hold on vacuum.
// Sweep xb ∈ {0.3, 0.6, 1.0, 1.5} × rad_dual ∈ {0, 1}.
func buildP2Arms(seed int, T float64) []Arm {
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
			"q_detune=1.2",
			"k_rad=0.05", "p_rad=4", "rad_clock=0",
			"bath_em=0.05",
		}
		return append(a, extra...)
	}
	var arms []Arm
	for _, xb := range []struct {
		tag string
		v   string
	}{
		{"xb03", "0.3"},
		{"xb06", "0.6"},
		{"xb10", "1.0"},
		{"xb15", "1.5"},
	} {
		for _, dual := range []struct {
			tag string
			v   string
		}{
			{"d0", "0"},
			{"d1", "1"},
		} {
			name := "P2_" + xb.tag + "_" + dual.tag
			arms = append(arms, Arm{
				Name:   name,
				Family: "p2",
				Desc:   "P2 dual intake xb=" + xb.v + " rad_dual=" + dual.v,
				Args:   base("bath_xb="+xb.v, "rad_dual="+dual.v),
				Snap:   true,
			})
		}
	}
	// Control: no radiance, dual off (pure hold baseline)
	arms = append(arms, Arm{
		Name:   "P2_ctrl_kr0",
		Family: "p2",
		Desc:   "P2 control k_rad=0 (no tax, no dual)",
		Args:   base("bath_xb=0.6", "rad_dual=0", "k_rad=0"),
		Snap:   true,
	})
	return arms
}
