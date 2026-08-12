package hunt

import "math"

// Currents are reconstructed from the snapshot (no live hop ledger).
// J_ij = 2 Im(ψ_i* ψ_j) with ψ = √Em e^{i th2}.
// Space flow rate uses the pass-S pressure law (not *dt, not outflow-clamped).

type Currents struct {
	// Net dense current Σ J û  (fd-like).
	JnetX, JnetY, JnetZ float64
	// Circulation proxy: Σ J on internal motif links (signed, no direction).
	Circ float64
	// Monopole mass: ½ Σ_cells |net outgoing J|. Zero on a pure loop.
	Jmono float64
	// Dipole of cell monopoles about the motif centroid.
	JdipX, JdipY, JdipZ float64
	// Space flux through the motif surface (links with one end inside).
	// Thru = ½ Σ |f|; Net > 0 leaves the motif.
	EsThru float64
	EsNet  float64
	// Internal space slosh (both ends inside).
	EsInt float64
}

func psi(c Cell) (m1, m2 float64) {
	if c.Em <= 0 {
		return 0, 0
	}
	a := math.Sqrt(c.Em)
	return a * math.Cos(c.Th2), a * math.Sin(c.Th2)
}

// linkJ is 2 Im(ψ_i* ψ_j).
func linkJ(a, b Cell) float64 {
	m1i, m2i := psi(a)
	m1j, m2j := psi(b)
	return 2 * (m1i*m2j - m2i*m1j)
}

func pressure(c Cell, sdisp float64) float64 {
	return c.Es + sdisp*(c.Em+c.Ee)
}

// GlobalJ is Σ J û over every link — the snapshot analog of fd.
func GlobalJ(fr Frame) (x, y, z float64) {
	return sumJ(fr, nil)
}

// TaggedJ is Σ J û on links with both ends tagged (the painted blob/ring).
func TaggedJ(fr Frame) (x, y, z float64) {
	tag := make([]bool, len(fr.Cells))
	n := 0
	for i, c := range fr.Cells {
		if c.Tag > 0.5 {
			tag[i] = true
			n++
		}
	}
	if n < 2 {
		return 0, 0, 0
	}
	return sumJ(fr, tag)
}

func sumJ(fr Frame, mask []bool) (x, y, z float64) {
	for _, ln := range fr.Links {
		i, j := ln.I, ln.J
		if i < 0 || j < 0 || i >= len(fr.Cells) || j >= len(fr.Cells) {
			continue
		}
		if mask != nil && !(mask[i] && mask[j]) {
			continue
		}
		ci, cj := fr.Cells[i], fr.Cells[j]
		dx := wrap(cj.X-ci.X, fr.L)
		dy := wrap(cj.Y-ci.Y, fr.L)
		dz := wrap(cj.Z-ci.Z, fr.L)
		d := ln.D
		if d <= 1e-12 {
			d = hypot3(dx, dy, dz)
		}
		if d <= 1e-12 {
			continue
		}
		J := linkJ(ci, cj)
		x += J * dx / d
		y += J * dy / d
		z += J * dz / d
	}
	return
}

// CoreFlux is pass-S surface flux of the high-Em core (Em > floor).
func CoreFlux(fr Frame, floor float64, s *Stream) (thru, net, em float64, n int) {
	if floor <= 0 {
		floor = 0.5
	}
	inside := make([]bool, len(fr.Cells))
	for i, c := range fr.Cells {
		if c.Em > floor {
			inside[i] = true
			em += c.Em
			n++
		}
	}
	sk, sdisp, aref, dref := defaultSK, defaultSDisp, math.Pi*defaultR0*defaultR0, 2*defaultR0
	if s != nil {
		sk, sdisp, aref, dref = s.SK, s.SDisp, s.Aref, s.Dref
	}
	for _, ln := range fr.Links {
		i, j := ln.I, ln.J
		if i < 0 || j < 0 || i >= len(fr.Cells) || j >= len(fr.Cells) {
			continue
		}
		inI, inJ := inside[i], inside[j]
		if inI == inJ {
			continue
		}
		var src, dst Cell
		sign := 1.0
		if inI {
			src, dst = fr.Cells[i], fr.Cells[j]
		} else {
			src, dst = fr.Cells[j], fr.Cells[i]
			sign = -1
		}
		d := ln.D
		if d <= 1e-12 {
			d = dist3(src.X, src.Y, src.Z, dst.X, dst.Y, dst.Z, fr.L)
		}
		if d <= 1e-12 {
			continue
		}
		w := 1.0
		if ln.A > 0 && aref > 0 && dref > 0 {
			w = (ln.A / aref) * (dref / d)
		}
		f := sk * w * (pressure(src, sdisp) - pressure(dst, sdisp)) * sign
		net += f
		thru += 0.5 * math.Abs(f)
	}
	return
}

func MeasureCurrents(fr Frame, m Motif, s *Stream) Currents {
	var out Currents
	if m.N == 0 || len(fr.Cells) == 0 {
		return out
	}
	inside := make([]bool, len(fr.Cells))
	for _, i := range m.Cells {
		if i >= 0 && i < len(inside) {
			inside[i] = true
		}
	}
	q := make([]float64, len(fr.Cells)) // net outgoing J
	sk, sdisp, aref, dref := defaultSK, defaultSDisp, math.Pi*defaultR0*defaultR0, 2*defaultR0
	if s != nil {
		sk, sdisp, aref, dref = s.SK, s.SDisp, s.Aref, s.Dref
		if aref <= 0 {
			aref = math.Pi * s.R0 * s.R0
		}
		if dref <= 0 {
			dref = 2 * s.R0
		}
	}
	L := fr.L
	for _, ln := range fr.Links {
		i, j := ln.I, ln.J
		if i < 0 || j < 0 || i >= len(fr.Cells) || j >= len(fr.Cells) {
			continue
		}
		inI, inJ := inside[i], inside[j]
		if !inI && !inJ {
			continue
		}
		ci, cj := fr.Cells[i], fr.Cells[j]
		dx := wrap(cj.X-ci.X, L)
		dy := wrap(cj.Y-ci.Y, L)
		dz := wrap(cj.Z-ci.Z, L)
		d := ln.D
		if d <= 1e-12 {
			d = hypot3(dx, dy, dz)
		}
		if d <= 1e-12 {
			continue
		}
		ux, uy, uz := dx/d, dy/d, dz/d

		if inI && inJ {
			J := linkJ(ci, cj)
			out.Circ += J
			out.JnetX += J * ux
			out.JnetY += J * uy
			out.JnetZ += J * uz
			q[i] += J
			q[j] -= J
			// internal space slosh
			dp := pressure(ci, sdisp) - pressure(cj, sdisp)
			w := 1.0
			if ln.A > 0 && aref > 0 {
				w = (ln.A / aref) * (dref / d)
			}
			out.EsInt += 0.5 * math.Abs(sk*w*dp)
		} else {
			// surface: one end in. Orient so + leaves the motif.
			var src, dst Cell
			var sign float64
			if inI {
				src, dst, sign = ci, cj, 1
			} else {
				src, dst, sign = cj, ci, -1
				ux, uy, uz = -ux, -uy, -uz
			}
			_ = uy
			_ = uz
			dp := pressure(src, sdisp) - pressure(dst, sdisp)
			w := 1.0
			if ln.A > 0 && aref > 0 {
				w = (ln.A / aref) * (dref / d)
			}
			f := sk * w * dp * sign
			// f > 0: pressure pushes out of the motif
			if f > 0 {
				out.EsNet += f
			} else {
				out.EsNet += f
			}
			out.EsThru += 0.5 * math.Abs(f)
		}
	}
	for _, i := range m.Cells {
		qi := q[i]
		out.Jmono += math.Abs(qi)
		dx := wrap(fr.Cells[i].X-m.Cx, L)
		dy := wrap(fr.Cells[i].Y-m.Cy, L)
		dz := wrap(fr.Cells[i].Z-m.Cz, L)
		out.JdipX += dx * qi
		out.JdipY += dy * qi
		out.JdipZ += dz * qi
	}
	out.Jmono *= 0.5
	return out
}

func (c Currents) JnetAbs() float64 {
	return hypot3(c.JnetX, c.JnetY, c.JnetZ)
}

func (c Currents) JdipAbs() float64 {
	return hypot3(c.JdipX, c.JdipY, c.JdipZ)
}
