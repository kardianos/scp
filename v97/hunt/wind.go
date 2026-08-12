package hunt

import "math"

// Winding of th2 on a motif, sorted by azimuth about the motif
// centroid (same construction as report/analyze_ring.py).
func Winding(fr Frame, m Motif) (W, Rm float64, n int) {
	type row struct{ phi, th, em float64 }
	rs := make([]row, 0, m.N)
	L := fr.L
	for _, i := range m.Cells {
		c := fr.Cells[i]
		dx := wrap(c.X-m.Cx, L)
		dy := wrap(c.Y-m.Cy, L)
		rs = append(rs, row{math.Atan2(dy, dx), c.Th2, c.Em})
	}
	n = len(rs)
	if n < 3 {
		return 0, 0, n
	}
	// insertion sort by phi
	for i := 1; i < n; i++ {
		j := i
		for j > 0 && rs[j].phi < rs[j-1].phi {
			rs[j], rs[j-1] = rs[j-1], rs[j]
			j--
		}
	}
	var sum float64
	for i := 0; i < n; i++ {
		d := rs[(i+1)%n].th - rs[i].th
		sum += wrapPi(d)
	}
	W = sum / (2 * math.Pi)
	// R_m at m = round(W), fallback m=2
	mm := int(math.Round(W))
	if mm == 0 {
		mm = 2
	}
	var re, im float64
	for _, r := range rs {
		ang := r.th - float64(mm)*r.phi
		re += math.Cos(ang)
		im += math.Sin(ang)
	}
	Rm = math.Hypot(re, im) / float64(n)
	return W, Rm, n
}
