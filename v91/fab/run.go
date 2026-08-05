package fab

// init + seeds + main loop + report — port of v89/freecell.c main()
// (lines 1393-2252). Log-line formats match the C instrument except the
// first line, which identifies the Go port (A/B comparators skip it).

import (
	"fmt"
	"io"
	"math"
	"os"
)

func fprintf(w io.Writer, format string, args ...any) {
	fmt.Fprintf(w, format, args...)
}

func (s *Sim) buildBath(bx, by, bz []float64, nmax int) int {
	P := &s.P
	gn := int(math.Ceil(P.L / P.Dmin))
	if gn < 1 {
		gn = 1
	}
	nbins := gn * gn * gn
	head := make([]int, nbins)
	nxt := make([]int, nmax)
	for b := 0; b < nbins; b++ {
		head[b] = -1
	}
	n, fails := 0, 0
	d2min := P.Dmin * P.Dmin
	for fails < 30000 && n < nmax {
		x := s.frand() * P.L
		y := s.frand() * P.L
		z := s.frand() * P.L
		bx2, by2, bz2 := int(x/P.L*float64(gn)), int(y/P.L*float64(gn)), int(z/P.L*float64(gn))
		if bx2 >= gn {
			bx2 = gn - 1
		}
		if by2 >= gn {
			by2 = gn - 1
		}
		if bz2 >= gn {
			bz2 = gn - 1
		}
		ok := true
		for ax := bx2 - 1; ax <= bx2+1 && ok; ax++ {
			for ay := by2 - 1; ay <= by2+1 && ok; ay++ {
				for az := bz2 - 1; az <= bz2+1 && ok; az++ {
					wx, wy2, wz := (ax+gn)%gn, (ay+gn)%gn, (az+gn)%gn
					for q := head[(wx*gn+wy2)*gn+wz]; q >= 0; q = nxt[q] {
						dx := s.wr(bx[q] - x)
						dy := s.wr(by[q] - y)
						dz := s.wr(bz[q] - z)
						if dx*dx+dy*dy+dz*dz < d2min {
							ok = false
							break
						}
					}
				}
			}
		}
		if !ok {
			fails++
			continue
		}
		bx[n], by[n], bz[n] = x, y, z
		b := (bx2*gn+by2)*gn + bz2
		nxt[n] = head[b]
		head[b] = n
		n++
		fails = 0
	}
	if P.BathFrac < 1.0 {
		m := 0
		for i := 0; i < n; i++ {
			if s.frand() < P.BathFrac {
				bx[m], by[m], bz[m] = bx[i], by[i], bz[i]
				m++
			}
		}
		n = m
	}
	return n
}

// 2D dart throw in the z = L/2 plane (DS medium; coplanar dynamics is
// exact — every force the kernel applies to coplanar cells is in-plane)
func (s *Sim) buildBath2d(bx, by, bz []float64, nmax int) int {
	P := &s.P
	gn := int(math.Ceil(P.L / P.Dmin))
	if gn < 1 {
		gn = 1
	}
	nbins := gn * gn
	head := make([]int, nbins)
	nxt := make([]int, nmax)
	for b := 0; b < nbins; b++ {
		head[b] = -1
	}
	n, fails := 0, 0
	d2min := P.Dmin * P.Dmin
	for fails < 30000 && n < nmax {
		x := s.frand() * P.L
		y := s.frand() * P.L
		bx2, by2 := int(x/P.L*float64(gn)), int(y/P.L*float64(gn))
		if bx2 >= gn {
			bx2 = gn - 1
		}
		if by2 >= gn {
			by2 = gn - 1
		}
		ok := true
		for ax := bx2 - 1; ax <= bx2+1 && ok; ax++ {
			for ay := by2 - 1; ay <= by2+1 && ok; ay++ {
				wx, wy2 := (ax+gn)%gn, (ay+gn)%gn
				for q := head[wx*gn+wy2]; q >= 0; q = nxt[q] {
					dx := s.wr(bx[q] - x)
					dy := s.wr(by[q] - y)
					if dx*dx+dy*dy < d2min {
						ok = false
						break
					}
				}
			}
		}
		if !ok {
			fails++
			continue
		}
		bx[n], by[n], bz[n] = x, y, 0.5*P.L
		b := bx2*gn + by2
		nxt[n] = head[b]
		head[b] = n
		n++
		fails = 0
	}
	if P.BathFrac < 1.0 {
		m := 0
		for i := 0; i < n; i++ {
			if s.frand() < P.BathFrac {
				bx[m], by[m], bz[m] = bx[i], by[i], bz[i]
				m++
			}
		}
		n = m
	}
	return n
}

// load a voice to occupancy x: Em = x*cap exactly, space pulled
func (s *Sim) loadVoice(i int, x float64) {
	P := &s.P
	add := x * P.Cap / (1.0 + P.SPull)
	s.Em[i] += add
	pull := P.SPull * add
	avail := s.Es[i] - P.EsFloor
	if pull > avail {
		if avail > 0 {
			pull = avail
		} else {
			pull = 0
		}
	}
	s.Es[i] -= pull
	s.Em[i] += pull
}

func (s *Sim) dstarOf(x float64, pp, qq, m int) float64 {
	w := s.P.W2 / (1.0 + s.P.QDetune*x)
	return TwoPi * float64(m) * s.P.C / (float64(qq+pp) * w)
}

func (s *Sim) seedCellDefaults(i int) {
	s.cr0[i] = s.P.R0
	s.cr[i] = s.cr0[i]
	s.n1x[i], s.n1y[i], s.n1z[i] = s.randUnit()
	s.n2x[i], s.n2y[i], s.n2z[i] = s.randUnit()
	s.th2[i] = s.frand() * TwoPi
	s.cbeta[i] = s.frand() * TwoPi
	s.Es[i] = s.P.Es0
}

func (s *Sim) normalsTransverse(i int, ux, uy, uz float64) {
	var ax, ay float64
	az := 0.0
	if math.Abs(ux) < 0.9 {
		ax, ay = 1, 0
	} else {
		ax, ay = 0, 1
	}
	v1x := uy*az - uz*ay
	v1y := uz*ax - ux*az
	v1z := ux*ay - uy*ax
	n := math.Sqrt(v1x*v1x + v1y*v1y + v1z*v1z)
	v1x /= n
	v1y /= n
	v1z /= n
	v2x := uy*v1z - uz*v1y
	v2y := uz*v1x - ux*v1z
	v2z := ux*v1y - uy*v1x
	s.n1x[i], s.n1y[i], s.n1z[i] = v1x, v1y, v1z
	s.n2x[i], s.n2y[i], s.n2z[i] = v2x, v2y, v2z
}

// Run executes the configured experiment. Mirrors C main().
func (s *Sim) Run() {
	if s.Out == nil {
		s.Out = os.Stdout
	}
	if s.Errw == nil {
		s.Errw = os.Stderr
	}
	P := &s.P

	if P.Exp == "p2lc" {
		// P2 local-clock prototype (Go-only; localclock.go)
		s.RunP2LC()
		return
	}

	s.combBuild()
	s.rng = P.Seed
	if s.rng == 0 {
		s.rng = 1
	}
	for w := 0; w < 8; w++ {
		s.xrand() // cellfab warm-up
	}

	if P.GammaResM < 0 {
		s.gm = P.GammaRes
	} else {
		s.gm = P.GammaResM
	}
	s.lockf = P.LockFloor * P.Cap
	s.Aref = math.Pi * P.R0 * P.R0
	s.dref = 2.0 * P.R0
	if P.QuantA0 < 0 {
		s.A0eff = P.Es0 * 1.5053 / P.C
	} else {
		s.A0eff = P.QuantA0
	}
	s.repPairI, s.repPairJ = -1, -1

	fprintf(s.Out, "# freecell-go — v90 Go port of v89/freecell.c (exp=%s)\n", P.Exp)
	fprintf(s.Out, "# laws = laws_V2g VERBATIM (defaults): C=%.6g w1=%.6g w2=%.6g q_detune=%.6g\n",
		P.C, P.W1, P.W2, P.QDetune)
	fprintf(s.Out, "# gamma_res=%.6g gamma_res_m=%.6g p_gate=%.6g lock_floor=%.6g k_dep=%.6g k_dep_m=%.6g cap=%.6g\n",
		P.GammaRes, P.GammaResM, P.PGate, P.LockFloor, P.KDep, P.KDepM, P.Cap)
	fprintf(s.Out, "# e_s0=%.6g es_floor=%.6g e_cond=%.6g f_conv=%.6g f_evap=%.6g s_pull=%.6g\n",
		P.Es0, P.EsFloor, P.ECond, P.FConv, P.FEvap, P.SPull)
	fprintf(s.Out, "# kappa_lock=%.6g kappa_align=%.6g s_k=%.6g s_disp=%.6g sigma_tumble=%.6g\n",
		P.KappaLock, P.KappaAlign, P.SK, P.SDisp, P.SigmaTumble)
	fprintf(s.Out, "# comb_limit=%d rough_k=%.6g gamma_rough=%.6g mob_sym=%d mob_floor=%.6g field_J=%.6g\n",
		P.CombLimit, P.RoughK, P.GammaRough, P.MobSym, P.MobFloor, P.FieldJ)
	fprintf(s.Out, "# quant_A0=%.6g quant_mode=%d (A0eff=%.6g)\n", P.QuantA0, P.QuantMode, s.A0eff)
	fprintf(s.Out, "# v91 radiance (laws_V3r candidate A): k_rad=%.6g p_rad=%.6g rad_clock=%d\n",
		P.KRad, P.PRad, P.RadClock)
	fprintf(s.Out, "# v91 cantus (coherent-channel candidate B): k_cant=%.6g k_tune=%.6g cant_tau=%.6g cant_seed=%d cant_grow=%d\n",
		P.KCant, P.KTune, P.CantTau, P.CantSeed, P.CantGrow)
	fprintf(s.Out, "# v91 registry (exchange-registry lane, REGISTRY.md): reg_tau=%.6g reg_gate=%d reg_f0=%.6g\n",
		P.RegTau, P.RegGate, P.RegF0)
	fprintf(s.Out, "# GEOMETRY (apparatus): cfac=%.6g k_rep=%.6g mob_geo=%.6g kappa_bond=%.6g freeze_geo=%d\n",
		P.Cfac, P.KRep, P.MobGeo, P.KappaBond, P.FreezeGeo)
	fprintf(s.Out, "# bath=%d bath_frac=%.6g jam_sweeps=%d jam_k=%.6g L=%.6g dt=%.6g T=%.6g seed=%d diag_every=%d\n",
		P.Bath, P.BathFrac, P.JamSweeps, P.JamK, P.L, P.Dt, P.T, P.Seed, P.DiagEvery)

	// ---------------- build ----------------
	V := P.L * P.L * P.L
	nmax := int(0.8*V/(P.Dmin*P.Dmin*P.Dmin)) + 4096
	bx := make([]float64, nmax)
	by := make([]float64, nmax)
	bz := make([]float64, nmax)
	nb := 0
	if P.Bath != 0 {
		if P.Exp == "slit" {
			nb = s.buildBath2d(bx, by, bz, nmax)
		} else {
			nb = s.buildBath(bx, by, bz, nmax)
		}
	}
	if P.GradR0 > 0 && P.GradFrac < 1.0 && nb > 0 {
		// FORGE: thin the dart bath inside r0 of the box centre
		gc := 0.5 * P.L
		m, nin := 0, 0
		for i := 0; i < nb; i++ {
			dx := s.wr(bx[i] - gc)
			dy := s.wr(by[i] - gc)
			dz := s.wr(bz[i] - gc)
			inside := dx*dx+dy*dy+dz*dz < P.GradR0*P.GradR0
			if inside {
				nin++
			}
			if !inside || s.frand() < P.GradFrac {
				bx[m], by[m], bz[m] = bx[i], by[i], bz[i]
				m++
			}
		}
		fprintf(s.Out, "# GRAD carved: removed %d of %d cells inside r=%.6g (frac=%.6g)\n",
			nb-m, nin, P.GradR0, P.GradFrac)
		nb = m
	}

	extra := 0
	switch {
	case P.Exp == "pair" && P.Bath == 0:
		extra = 2
	case P.Exp == "ring":
		extra = P.RingN
	case P.Exp == "oct":
		extra = 6
	case P.Exp == "tri":
		extra = 3
	case P.Exp == "tri2":
		extra = 6
	case P.Exp == "rings":
		g := 3
		if P.RingsKind == 0 {
			g = 2
		} else if P.RingsKind == 3 {
			g = 6
		}
		extra = g * P.RingsNv
	}

	s.allocAll(nb + extra)
	for i := 0; i < nb; i++ {
		s.px[i], s.py[i], s.pz[i] = bx[i], by[i], bz[i]
		// cellfab per-cell draw order: cr0 jitter, n1, n2, th (th1
		// skipped: field phase lives in fa), th2, cbeta
		s.cr0[i] = P.R0 * (1.0 + P.Rjit*(2.0*s.frand()-1.0))
		s.cr[i] = s.cr0[i]
		s.n1x[i], s.n1y[i], s.n1z[i] = s.randUnit()
		s.n2x[i], s.n2y[i], s.n2z[i] = s.randUnit()
		_ = s.frand() // th1 draw consumed, stream aligned
		s.th2[i] = s.frand() * TwoPi
		s.cbeta[i] = s.frand() * TwoPi
		s.Es[i] = P.Es0
	}
	s.NC = nb + extra

	cx0 := 0.5 * P.L
	if P.Px >= 0 {
		cx0 = P.Px
	}
	cy0 := 0.5 * P.L
	if P.Py >= 0 {
		cy0 = P.Py
	}
	cz0 := 0.5 * P.L
	if P.Pz >= 0 {
		cz0 = P.Pz
	}
	s.cenx, s.ceny, s.cenz = cx0, cy0, cz0

	// jam-settle the bath BEFORE seeding
	if nb > 0 {
		ncFull := s.NC
		s.NC = nb
		s.jamSettle()
		s.NC = ncFull
	} else {
		fprintf(s.Out, "# JAM skipped (no bath)\n")
	}

	// ---------------- seeds ----------------
	switch P.Exp {
	case "pair":
		var i, j int
		x := P.PairX0
		x1 := P.PairX1
		if x1 < 0 {
			x1 = x
		}
		wA := P.W2 / (1.0 + P.QDetune*x)
		wB := P.W2 / (1.0 + P.QDetune*x1)
		ds := TwoPi * float64(P.PairM) * P.C / (float64(P.PairQQ)*wA + float64(P.PairPP)*wB)
		d0 := ds + P.PairDoff
		if P.Bath == 0 {
			i, j = nb, nb+1
			s.seedCellDefaults(i)
			s.seedCellDefaults(j)
			s.px[i], s.py[i], s.pz[i] = cx0-0.5*d0, cy0, cz0
			s.px[j], s.py[j], s.pz[j] = cx0+0.5*d0, cy0, cz0
		} else {
			bi, bj := -1, -1
			bcost := 1e300
			for a := 0; a < nb; a++ {
				for b2 := a + 1; b2 < nb; b2++ {
					dx := s.wr(s.px[b2] - s.px[a])
					dy := s.wr(s.py[b2] - s.py[a])
					dz := s.wr(s.pz[b2] - s.pz[a])
					d := math.Sqrt(dx*dx + dy*dy + dz*dz)
					mx := s.wr(s.px[a]-cx0) + 0.5*dx
					my := s.wr(s.py[a]-cy0) + 0.5*dy
					mz := s.wr(s.pz[a]-cz0) + 0.5*dz
					cost := math.Abs(d-d0)*10 + math.Sqrt(mx*mx+my*my+mz*mz)
					if d > 0.9*d0 && d < 1.3*d0 && cost < bcost {
						bcost = cost
						bi, bj = a, b2
					}
				}
			}
			i, j = bi, bj
			if i < 0 {
				s.fatal(1, "no bath pair found near d0\n")
			}
		}
		dx := s.wr(s.px[j] - s.px[i])
		dy := s.wr(s.py[j] - s.py[i])
		dz := s.wr(s.pz[j] - s.pz[i])
		d := math.Sqrt(dx*dx + dy*dy + dz*dz)
		s.normalsTransverse(i, dx/d, dy/d, dz/d)
		s.normalsTransverse(j, dx/d, dy/d, dz/d)
		if x > 0 {
			s.loadVoice(i, x)
		}
		if x1 > 0 {
			s.loadVoice(j, x1)
		}
		if P.Seedlock != 0 {
			s.th2[i] = s.frand() * TwoPi
			s.th2[j] = math.Mod((float64(P.PairQQ)*(s.th2[i]-wA*d/P.C))/float64(P.PairPP)+
				8.0*TwoPi, TwoPi)
		}
		s.tag[i] = 1
		s.tag[j] = 1
		s.repPairI, s.repPairJ = i, j
		fprintf(s.Out, "# SEED pair: i=%d j=%d x=%.4f/%.4f p:q=%d:%d m=%d d*=%.6f d0=%.6f (doff %+0.4f)\n",
			i, j, x, x1, P.PairPP, P.PairQQ, P.PairM, ds, d, P.PairDoff)
		fprintf(s.Out, "# SEED pair radii r=%.4f contact=%.4f  Em/voice=%.4f\n",
			s.cr0[i]*math.Cbrt(s.Es[i]), s.cr0[i]*math.Cbrt(s.Es[i])+s.cr0[j]*math.Cbrt(s.Es[j]), s.Em[i])

	case "ring", "oct":
		nvo := P.RingN
		x := P.RingX
		doff := P.RingDoff
		if P.Exp == "oct" {
			nvo = 6
			x = P.OctX
			doff = P.OctDoff
		}
		ds := s.dstarOf(x, 1, 1, P.PairM)
		de := ds + doff
		// bath embedding: carve a cavity
		if P.Bath != 0 && nb > 0 {
			den := nvo
			if den <= 2 {
				den = 3
			}
			R := de / (2.0 * math.Sin(math.Pi/float64(den)))
			clear := R + P.Dmin + 0.2
			m := 0
			for i := 0; i < nb; i++ {
				dx := s.wr(s.px[i] - cx0)
				dy := s.wr(s.py[i] - cy0)
				dz := s.wr(s.pz[i] - cz0)
				keepR2 := dx*dx + dy*dy + dz*dz
				keep := keepR2 >= clear*clear
				if keep {
					if m != i {
						s.px[m], s.py[m], s.pz[m] = s.px[i], s.py[i], s.pz[i]
						s.cr0[m], s.cr[m] = s.cr0[i], s.cr[i]
						s.n1x[m], s.n1y[m], s.n1z[m] = s.n1x[i], s.n1y[i], s.n1z[i]
						s.n2x[m], s.n2y[m], s.n2z[m] = s.n2x[i], s.n2y[i], s.n2z[i]
						s.th2[m], s.cbeta[m] = s.th2[i], s.cbeta[i]
						s.Es[m], s.Em[m], s.Ee[m] = s.Es[i], s.Em[i], s.Ee[i]
						s.fa1[m], s.fa2[m] = s.fa1[i], s.fa2[i]
					}
					m++
				}
			}
			fprintf(s.Out, "# CARVE removed %d bath cells (clear=%.2f)\n", nb-m, clear)
			nb = m
			s.NC = nb + extra
		}
		baseI := nb
		s.trussN = 0
		if P.Exp == "oct" {
			a := de / math.Sqrt(2.0)
			vx := [6]float64{a, -a, 0, 0, 0, 0}
			vy := [6]float64{0, 0, a, -a, 0, 0}
			vz := [6]float64{0, 0, 0, 0, a, -a}
			for k := 0; k < 6; k++ {
				i := baseI + k
				s.seedCellDefaults(i)
				s.px[i], s.py[i], s.pz[i] = cx0+vx[k], cy0+vy[k], cz0+vz[k]
				// global-aligned normals: dense plane n2 = z, n1 = x
				s.n1x[i], s.n1y[i], s.n1z[i] = 1, 0, 0
				s.n2x[i], s.n2y[i], s.n2z[i] = 0, 0, 1
				s.loadVoice(i, x)
				s.tag[i] = 1
			}
			for a2 := 0; a2 < 6; a2++ {
				for b2 := a2 + 1; b2 < 6; b2++ {
					dx := s.px[baseI+b2] - s.px[baseI+a2]
					dy := s.py[baseI+b2] - s.py[baseI+a2]
					dz := s.pz[baseI+b2] - s.pz[baseI+a2]
					d := math.Sqrt(dx*dx + dy*dy + dz*dz)
					if d < 1.3*de {
						s.trussE[s.trussN][0] = baseI + a2
						s.trussE[s.trussN][1] = baseI + b2
						s.trussDstar[s.trussN] = ds
						s.trussN++
					}
				}
			}
		} else {
			R := de / (2.0 * math.Sin(math.Pi/float64(nvo)))
			for k := 0; k < nvo; k++ {
				i := baseI + k
				s.seedCellDefaults(i)
				a := TwoPi * float64(k) / float64(nvo)
				s.px[i] = cx0 + R*math.Cos(a)
				s.py[i] = cy0 + R*math.Sin(a)
				s.pz[i] = cz0
				// COPLANAR normals: both planes = the ring plane (gpl=1)
				s.n1x[i], s.n1y[i], s.n1z[i] = 0, 0, 1
				s.n2x[i], s.n2y[i], s.n2z[i] = 0, 0, 1
				s.loadVoice(i, x)
				s.tag[i] = 1
			}
			for k := 0; k < nvo; k++ {
				s.trussE[s.trussN][0] = baseI + k
				s.trussE[s.trussN][1] = baseI + (k+1)%nvo
				s.trussDstar[s.trussN] = ds
				s.trussN++
			}
		}
		// phases: lock around the edge list (unison chain)
		if P.Seedlock != 0 {
			we := P.W2 / (1.0 + P.QDetune*x)
			s.th2[baseI] = s.frand() * TwoPi
			for e := 0; e < s.trussN; e++ {
				i, j := s.trussE[e][0], s.trussE[e][1]
				dx := s.px[j] - s.px[i]
				dy := s.py[j] - s.py[i]
				dz := s.pz[j] - s.pz[i]
				d := math.Sqrt(dx*dx + dy*dy + dz*dz)
				s.th2[j] = math.Mod(s.th2[i]-we*d/P.C+8.0*TwoPi, TwoPi)
			}
		}
		fprintf(s.Out, "# SEED %s: n=%d x=%.4f d*=%.6f d_edge=%.6f edges=%d\n",
			P.Exp, nvo, x, ds, de, s.trussN)

	case "blob":
		// e3a-class Gaussian dense blob (cellfab.c:1754 recipe)
		nload := 0
		for i := 0; i < s.NC; i++ {
			dx := s.wr(s.px[i] - cx0)
			dy := s.wr(s.py[i] - cy0)
			dz := s.wr(s.pz[i] - cz0)
			r2 := dx*dx + dy*dy + dz*dz
			add := P.Amp * math.Exp(-r2/(2.0*P.Sigma*P.Sigma))
			if add < 1e-4 {
				continue
			}
			if s.Em[i]+add > 0.95*P.Cap {
				add = 0.95*P.Cap - s.Em[i]
			}
			if add <= 0 {
				continue
			}
			s.Em[i] += add
			pull := P.SPull * add
			avail := s.Es[i] - P.EsFloor
			if pull > avail {
				if avail > 0 {
					pull = avail
				} else {
					pull = 0
				}
			}
			s.Es[i] -= pull
			s.Em[i] += pull
			// e3b tilt: phase ramp along x
			if P.Kx != 0 {
				s.th2[i] = math.Mod(-(P.Kx*dx)+8.0*TwoPi, TwoPi)
			} else {
				s.th2[i] = 0
			}
			if add > 0.05*P.Amp {
				s.tag[i] = 1
				nload++
			}
		}
		if P.Kx != 0 {
			for i := 0; i < s.NC; i++ {
				s.normalsTransverse(i, 1, 0, 0)
			}
		}
		fprintf(s.Out, "# SEED blob: amp=%.6g sigma=%.6g kx=%.6g tagged=%d (add>5%% of amp)\n",
			P.Amp, P.Sigma, P.Kx, nload)

	case "tri", "tri2":
		// FCQ — the fifth-triangle: the incomplete harmonic as binding
		xU := P.TriXU
		xD := P.TriXD
		if xD < 0 {
			xD = (0.5 + 1.8*xU) / 1.2
		}
		wU := P.W2 / (1.0 + P.QDetune*xU)
		wD := P.W2 / (1.0 + P.QDetune*xD)
		dUU := math.Pi / wU
		dDD := math.Pi / wD
		dUD := 2.0 * TwoPi / (2.0*wU + 3.0*wD)
		s.ntri = 1
		if P.Exp == "tri2" {
			s.ntri = 2
		}
		s.triOn = true
		baseI := nb
		s.trussN = 0
		for t := 0; t < s.ntri; t++ {
			ox := 0.0
			kbr := P.TriBranch
			if t != 0 {
				ox = P.Tri2Sep
				kbr = P.Tri2K2
			}
			v0 := baseI + 3*t
			for kk := 0; kk < 3; kk++ {
				i := v0 + kk
				s.seedCellDefaults(i)
				s.n1x[i], s.n1y[i], s.n1z[i] = 0, 0, 1 // coplanar, gpl=1
				s.n2x[i], s.n2y[i], s.n2z[i] = 0, 0, 1
				s.tag[i] = 1
				s.triV[t][kk] = i
			}
			if P.TriKind == 0 {
				// UUD equilateral, edge dUU (= dUD to FP)
				de := dUU + P.TriDoff
				R := de / math.Sqrt(3.0)
				for kk := 0; kk < 3; kk++ {
					a := math.Pi/2.0 + TwoPi*float64(kk)/3.0
					s.px[v0+kk] = s.fold(cx0 + ox + R*math.Cos(a))
					s.py[v0+kk] = s.fold(cy0 + R*math.Sin(a))
					s.pz[v0+kk] = cz0
				}
				s.loadVoice(v0+0, xU)
				s.loadVoice(v0+1, xU)
				s.loadVoice(v0+2, xD)
				s.triP[t][0], s.triQ[t][0] = 1, 1
				s.triP[t][1], s.triQ[t][1] = 3, 2
				s.triP[t][2], s.triQ[t][2] = 2, 3
				if P.Seedlock != 0 {
					s.th2[v0] = s.frand() * TwoPi
					s.th2[v0+1] = math.Mod(s.th2[v0]-wU*de+8.0*TwoPi, TwoPi)
					s.th2[v0+2] = math.Mod((2.0*(s.th2[v0+1]-wU*dUD)+
						TwoPi*float64(kbr))/3.0+8.0*TwoPi, TwoPi)
				}
				s.trussE[s.trussN][0] = v0
				s.trussE[s.trussN][1] = v0 + 1
				s.trussDstar[s.trussN] = dUU
				s.trussN++
				s.trussE[s.trussN][0] = v0 + 1
				s.trussE[s.trussN][1] = v0 + 2
				s.trussDstar[s.trussN] = dUD
				s.trussN++
				s.trussE[s.trussN][0] = v0 + 2
				s.trussE[s.trussN][1] = v0
				s.trussDstar[s.trussN] = dUD
				s.trussN++
			} else {
				// UDD isoceles: U apex, D-D base at dDD, sides dUD
				hb := 0.5 * dDD
				hh := math.Sqrt(dUD*dUD - hb*hb)
				s.px[v0+0] = s.fold(cx0 + ox)
				s.py[v0+0] = s.fold(cy0 + hh)
				s.pz[v0+0] = cz0
				s.px[v0+1] = s.fold(cx0 + ox - hb)
				s.py[v0+1] = cy0
				s.pz[v0+1] = cz0
				s.px[v0+2] = s.fold(cx0 + ox + hb)
				s.py[v0+2] = cy0
				s.pz[v0+2] = cz0
				s.loadVoice(v0+0, xU)
				s.loadVoice(v0+1, xD)
				s.loadVoice(v0+2, xD)
				s.triP[t][0], s.triQ[t][0] = 3, 2 // U->D fifth
				s.triP[t][1], s.triQ[t][1] = 1, 1 // D->D unison
				s.triP[t][2], s.triQ[t][2] = 2, 3 // D->U inverse
				if P.Seedlock != 0 {
					s.th2[v0] = s.frand() * TwoPi
					s.th2[v0+1] = math.Mod((2.0*(s.th2[v0]-wU*dUD)+
						TwoPi*float64(kbr))/3.0+8.0*TwoPi, TwoPi)
					s.th2[v0+2] = math.Mod(s.th2[v0+1]-wD*dDD+8.0*TwoPi, TwoPi)
				}
				s.trussE[s.trussN][0] = v0
				s.trussE[s.trussN][1] = v0 + 1
				s.trussDstar[s.trussN] = dUD
				s.trussN++
				s.trussE[s.trussN][0] = v0 + 1
				s.trussE[s.trussN][1] = v0 + 2
				s.trussDstar[s.trussN] = dDD
				s.trussN++
				s.trussE[s.trussN][0] = v0 + 2
				s.trussE[s.trussN][1] = v0
				s.trussDstar[s.trussN] = dUD
				s.trussN++
			}
		}
		kindStr := "UUD"
		if P.TriKind != 0 {
			kindStr = "UDD"
		}
		fprintf(s.Out, "# SEED tri: kind=%s ntri=%d xU=%.4f xD=%.4f wU=%.4f wD=%.4f (wD/wU=%.6f)\n",
			kindStr, s.ntri, xU, xD, wU, wD, wD/wU)
		fprintf(s.Out, "# SEED tri: d*_UU=%.6f d*_UD(m=2)=%.6f d*_DD=%.6f branch_k=%d doff=%+.3f\n",
			dUU, dUD, dDD, P.TriBranch, P.TriDoff)
		ctD := math.Cbrt(1.0 - P.SPull*xD*P.Cap/(1+P.SPull))
		ddOpen := "closed"
		if dDD > 2*P.R0*ctD {
			ddOpen = "OPEN (no contact)"
		}
		fprintf(s.Out, "# SEED tri: contacts UU=%.4f UD=%.4f DD=%.4f (DD edge %s)\n",
			2*P.R0*math.Cbrt(s.Es[baseI]),
			P.R0*(math.Cbrt(s.Es[baseI])+ctD),
			2*P.R0*ctD, ddOpen)

	case "pulse":
		// Tier-A e2 analog: field packet with phase tilt, prealigned
		for i := 0; i < s.NC; i++ {
			dx := s.wr(s.px[i] - cx0)
			dy := s.wr(s.py[i] - cy0)
			dz := s.wr(s.pz[i] - cz0)
			g := math.Exp(-(dx*dx + dy*dy + dz*dz) / (2.0 * P.Sigma * P.Sigma))
			if g < 1e-8 {
				continue
			}
			tilt := -(P.Kx * dx)
			s.fa1[i] += math.Sqrt(P.Amp*g) * math.Cos(tilt)
			s.fa2[i] += math.Sqrt(P.Amp*g) * math.Sin(tilt)
			s.Ee[i] = s.fa1[i]*s.fa1[i] + s.fa2[i]*s.fa2[i]
		}
		for i := 0; i < s.NC; i++ {
			s.normalsTransverse(i, 1, 0, 0)
		}
		fprintf(s.Out, "# SEED pulse: amp=%.6g sigma=%.6g kx=%.6g (prealigned transverse)\n",
			P.Amp, P.Sigma, P.Kx)

	case "blob2":
		// MOTION #28: two dense blobs with opposing tilts (approach),
		// per-blob tags for COM/exchange meters
		cxA := s.fold(cx0 - 0.5*P.Blob2Sep)
		cxB := s.fold(cx0 + 0.5*P.Blob2Sep)
		nA, nB := 0, 0
		for i := 0; i < s.NC; i++ {
			dxa := s.wr(s.px[i] - cxA)
			dya := s.wr(s.py[i] - cy0)
			dza := s.wr(s.pz[i] - cz0)
			dxb := s.wr(s.px[i] - cxB)
			dyb := s.wr(s.py[i] - cy0)
			dzb := s.wr(s.pz[i] - cz0)
			ga := P.Amp * math.Exp(-(dxa*dxa+dya*dya+dza*dza)/(2.0*P.Sigma*P.Sigma))
			gb := P.Amp * math.Exp(-(dxb*dxb+dyb*dyb+dzb*dzb)/(2.0*P.Sigma*P.Sigma))
			add := ga + gb
			if add < 1e-4 {
				continue
			}
			if s.Em[i]+add > 0.95*P.Cap {
				add = 0.95*P.Cap - s.Em[i]
			}
			if add <= 0 {
				continue
			}
			s.Em[i] += add
			pull := P.SPull * add
			avail := s.Es[i] - P.EsFloor
			if pull > avail {
				if avail > 0 {
					pull = avail
				} else {
					pull = 0
				}
			}
			s.Es[i] -= pull
			s.Em[i] += pull
			// tilt toward each other; phase ramp of the DOMINANT blob
			if P.Blob2Kx != 0 {
				if ga >= gb {
					s.th2[i] = math.Mod(-(P.Blob2Kx*dxa)+8.0*TwoPi, TwoPi)
				} else {
					s.th2[i] = math.Mod(+(P.Blob2Kx*dxb)+8.0*TwoPi, TwoPi)
				}
			}
			if ga >= gb && ga > 0.05*P.Amp {
				s.tag[i] = 1
				nA++
			} else if gb > 0.05*P.Amp {
				s.tag[i] = 2
				nB++
			}
		}
		if P.Blob2Kx != 0 {
			for i := 0; i < s.NC; i++ {
				s.normalsTransverse(i, 1, 0, 0)
			}
		}
		fprintf(s.Out, "# SEED blob2: amp=%.6g sigma=%.6g sep=%.6g kx=%.6g tagged A=%d B=%d\n",
			P.Amp, P.Sigma, P.Blob2Sep, P.Blob2Kx, nA, nB)

	case "rings":
		// COMPOSITE.md CQ0/CQ3b/CQ9: groups of rings as many-celled
		// quarks. kind 0 pair, 1 UUD, 2 UDD, 3 = TWO UUD nucleons.
		nv := P.RingsNv
		ngrp := 1
		if P.RingsKind == 3 {
			ngrp = 2
		}
		npg := 3
		if P.RingsKind == 0 {
			npg = 2
		}
		nring := ngrp * npg
		xU := P.RingsXU
		xD := P.RingsXD
		if xD < 0 {
			xD = (0.5 + 1.8*xU) / 1.2
		}
		wU := P.W2 / (1.0 + P.QDetune*xU)
		wD := P.W2 / (1.0 + P.QDetune*xD)
		dUD := 2.0 * TwoPi / (2.0*wU + 3.0*wD)
		var xs, ws [6]float64
		for k := 0; k < nring; k++ {
			lk := k % npg
			x := xU
			if P.RingsKind == 2 && lk >= 1 {
				x = xD
			}
			if (P.RingsKind == 1 || P.RingsKind == 3) && lk == 2 {
				x = xD
			}
			xs[k] = x
			ws[k] = P.W2 / (1.0 + P.QDetune*x)
		}
		var rload, dint, R [6]float64
		var mode [6]int
		for k := 0; k < nring; k++ {
			pull := P.SPull * (xs[k] * P.Cap / (1.0 + P.SPull))
			cap2 := P.Es0 - P.EsFloor
			if pull > cap2 {
				pull = cap2
			}
			es := P.Es0 - pull
			rload[k] = P.R0 * math.Cbrt(es/P.Es0)
			rung := math.Pi / ws[k]
			contact := 2.0 * rload[k]
			if rung < contact*P.Cfac {
				dint[k] = rung
				mode[k] = 1
			} else {
				dint[k] = 0.98 * contact
				mode[k] = 0
			}
			R[k] = dint[k] / (2.0 * math.Sin(math.Pi/float64(nv)))
		}
		// group-local centers (SSS), then group offsets
		var cxr, cyr [6]float64
		for g := 0; g < ngrp; g++ {
			b0 := g * npg
			gx := 0.0
			if ngrp == 2 {
				if g == 0 {
					gx = -0.5 * P.RingsSep
				} else {
					gx = 0.5 * P.RingsSep
				}
			}
			if npg == 2 {
				bg := math.Pi/ws[b0] + P.RingsGapoff
				cdist := R[b0] + R[b0+1] + bg
				cxr[b0], cyr[b0] = gx-0.5*cdist, 0
				cxr[b0+1], cyr[b0+1] = gx+0.5*cdist, 0
			} else {
				bgOf := func(a, b int) float64 {
					if math.Abs(ws[a]-ws[b]) < 1e-12 {
						return math.Pi/ws[a] + P.RingsGapoff
					}
					return dUD + P.RingsGapoff
				}
				bg01 := bgOf(b0, b0+1)
				bg12 := bgOf(b0+1, b0+2)
				bg02 := bgOf(b0, b0+2)
				c01 := R[b0] + R[b0+1] + bg01
				c12 := R[b0+1] + R[b0+2] + bg12
				c02 := R[b0] + R[b0+2] + bg02
				x0, y0, x1, y1 := 0.0, 0.0, c01, 0.0
				x2 := (c01*c01 + c02*c02 - c12*c12) / (2.0 * c01)
				y2s := c02*c02 - x2*x2
				y2 := 0.0
				if y2s > 0 {
					y2 = math.Sqrt(y2s)
				}
				mx := (x0 + x1 + x2) / 3.0
				my := (y0 + y1 + y2) / 3.0
				cxr[b0], cyr[b0] = gx+x0-mx, y0-my
				cxr[b0+1], cyr[b0+1] = gx+x1-mx, y1-my
				cxr[b0+2], cyr[b0+2] = gx+x2-mx, y2-my
			}
		}
		// bath: spherical carve per ring
		if P.Bath != 0 && nb > 0 {
			m, rem := 0, 0
			for i := 0; i < nb; i++ {
				keep := true
				for k := 0; k < nring && keep; k++ {
					dx := s.wr(s.px[i] - (cx0 + cxr[k]))
					dy := s.wr(s.py[i] - (cy0 + cyr[k]))
					dz := s.wr(s.pz[i] - cz0)
					clear := R[k] + 0.8
					if dx*dx+dy*dy+dz*dz < clear*clear {
						keep = false
					}
				}
				if keep {
					if m != i {
						s.px[m], s.py[m], s.pz[m] = s.px[i], s.py[i], s.pz[i]
						s.cr0[m], s.cr[m] = s.cr0[i], s.cr[i]
						s.n1x[m], s.n1y[m], s.n1z[m] = s.n1x[i], s.n1y[i], s.n1z[i]
						s.n2x[m], s.n2y[m], s.n2z[m] = s.n2x[i], s.n2y[i], s.n2z[i]
						s.th2[m], s.cbeta[m] = s.th2[i], s.cbeta[i]
						s.Es[m], s.Em[m], s.Ee[m] = s.Es[i], s.Em[i], s.Ee[i]
						s.fa1[m], s.fa2[m] = s.fa1[i], s.fa2[i]
					}
					m++
				} else {
					rem++
				}
			}
			fprintf(s.Out, "# CARVE rings: removed %d bath cells\n", rem)
			nb = m
			s.NC = nb + extra
		}
		s.ringsNring = nring
		baseI := nb
		s.trussN = 0
		for k := 0; k < nring; k++ {
			s.ringsV0[k] = baseI + k*nv
			s.ringsMode[k] = mode[k]
			s.ringsW[k] = ws[k]
			s.ringsGrp[k] = k / npg
			g, lk := k/npg, k%npg
			nbr := g*npg + (lk+1)%npg
			az := math.Atan2(cyr[nbr]-cyr[k], cxr[nbr]-cxr[k])
			for q := 0; q < nv; q++ {
				i := s.ringsV0[k] + q
				s.seedCellDefaults(i)
				a := az + TwoPi*float64(q)/float64(nv)
				s.px[i] = s.fold(cx0 + cxr[k] + R[k]*math.Cos(a))
				s.py[i] = s.fold(cy0 + cyr[k] + R[k]*math.Sin(a))
				s.pz[i] = cz0
				s.n1x[i], s.n1y[i], s.n1z[i] = 0, 0, 1
				s.n2x[i], s.n2y[i], s.n2z[i] = 0, 0, 1
				s.tag[i] = 1
			}
			for q := 0; q < nv; q++ {
				s.trussE[s.trussN][0] = s.ringsV0[k] + q
				s.trussE[s.trussN][1] = s.ringsV0[k] + (q+1)%nv
				s.trussDstar[s.trussN] = dint[k]
				s.trussN++
			}
		}
		// facing vertices (within group)
		var face [6][6]int
		for a := 0; a < nring; a++ {
			for b2 := 0; b2 < nring; b2++ {
				face[a][b2] = -1
				if a == b2 || s.ringsGrp[a] != s.ringsGrp[b2] {
					continue
				}
				az := math.Atan2(cyr[b2]-cyr[a], cxr[b2]-cxr[a])
				g, la := s.ringsGrp[a], a%npg
				nbr0 := g*npg + (la+1)%npg
				az0 := math.Atan2(cyr[nbr0]-cyr[a], cxr[nbr0]-cxr[a])
				rel := az - az0
				q := int(math.Floor(rel/(TwoPi/float64(nv)) + 0.5))
				q = ((q % nv) + nv) % nv
				face[a][b2] = s.ringsV0[a] + q
			}
		}
		// loads: per-flavor WHOLE-RING flight compensation (CQ3b v2)
		for k := 0; k < nring; k++ {
			xl := xs[k]
			if math.Abs(ws[k]-wU) < 1e-12 {
				xl -= P.RingsXcomp
			} else {
				xl -= P.RingsXcompD
			}
			if xl < 0 {
				xl = 0
			}
			for q := 0; q < nv; q++ {
				s.loadVoice(s.ringsV0[k]+q, xl)
			}
		}
		// inter-ring boundary edges (within groups)
		s.rintN = 0
		for g := 0; g < ngrp; g++ {
			nedge := 3
			if npg == 2 {
				nedge = 1
			}
			for e := 0; e < nedge; e++ {
				a := g*npg + e
				b2 := g*npg + (e+1)%npg
				ia, ib := face[a][b2], face[b2][a]
				s.rintE[s.rintN][0] = ia
				s.rintE[s.rintN][1] = ib
				if math.Abs(ws[a]-ws[b2]) < 1e-12 {
					s.rintP[s.rintN] = 1
					s.rintQ[s.rintN] = 1
					s.rintDstar[s.rintN] = math.Pi / ws[a]
				} else if ws[a] > ws[b2] {
					s.rintP[s.rintN] = 3
					s.rintQ[s.rintN] = 2
					s.rintDstar[s.rintN] = dUD
				} else {
					s.rintP[s.rintN] = 2
					s.rintQ[s.rintN] = 3
					s.rintDstar[s.rintN] = dUD
				}
				s.rintN++
			}
		}
		// phases: per-group chains; branch on the first U->D edge
		if P.Seedlock != 0 {
			for g := 0; g < ngrp; g++ {
				b0 := g * npg
				s.th2[s.ringsV0[b0]] = s.frand() * TwoPi
				nedge := 3
				if npg == 2 {
					nedge = 1
				}
				for k := b0; k < b0+npg; k++ {
					if k > b0 {
						e := g*nedge + (k - b0 - 1)
						ia, ib := s.rintE[e][0], s.rintE[e][1]
						dx := s.wr(s.px[ib] - s.px[ia])
						dy := s.wr(s.py[ib] - s.py[ia])
						d := math.Sqrt(dx*dx + dy*dy)
						if s.rintP[e] == 1 && s.rintQ[e] == 1 {
							s.th2[ib] = math.Mod(s.th2[ia]-ws[k-1]*d/P.C+8.0*TwoPi, TwoPi)
						} else if s.rintP[e] == 3 {
							s.th2[ib] = math.Mod((2.0*(s.th2[ia]-ws[k-1]*d/P.C)+
								TwoPi*float64(P.RingsBranch))/3.0+8.0*TwoPi, TwoPi)
						} else {
							s.th2[ib] = math.Mod((3.0*(s.th2[ia]-ws[k-1]*d/P.C))/2.0+
								8.0*TwoPi, TwoPi)
						}
					}
					start := s.ringsV0[b0]
					if k != b0 {
						start = s.rintE[g*nedge+(k-b0-1)][1]
					}
					sq := start - s.ringsV0[k]
					for st := 1; st < nv; st++ {
						qa, qb := (sq+st-1)%nv, (sq+st)%nv
						ia, ib := s.ringsV0[k]+qa, s.ringsV0[k]+qb
						dx := s.wr(s.px[ib] - s.px[ia])
						dy := s.wr(s.py[ib] - s.py[ia])
						d := math.Sqrt(dx*dx + dy*dy)
						s.th2[ib] = math.Mod(s.th2[ia]-ws[k]*d/P.C+8.0*TwoPi, TwoPi)
					}
				}
			}
		}
		fprintf(s.Out, "# SEED rings: kind=%d ngrp=%d nv=%d xU=%.4f xD=%.4f wU=%.4f wD=%.4f xcomp=%.4f sep=%.2f\n",
			P.RingsKind, ngrp, nv, xU, xD, wU, wD, P.RingsXcomp, P.RingsSep)
		for k := 0; k < nring; k++ {
			ms := "DROPLET"
			if mode[k] != 0 {
				ms = "MOLECULE"
			}
			fprintf(s.Out, "# SEED ring%d: grp=%d mode=%s x=%.4f w=%.4f d_int=%.4f R=%.4f\n",
				k, s.ringsGrp[k], ms, xs[k], ws[k], dint[k], R[k])
		}
		for e := 0; e < s.rintN; e++ {
			fprintf(s.Out, "# SEED rint%d: %d-%d p:q=%d:%d d*=%.6f\n",
				e, s.rintE[e][0], s.rintE[e][1], s.rintP[e], s.rintQ[e], s.rintDstar[e])
		}

	case "slit":
		// DS tier-0 (DS.md): carve the vacuum wall out of the jammed 2D
		// bath (except the slit bridges), then seed a quasi-plane packet
		yA := cy0 - 0.5*P.SlitSep
		yB := cy0 + 0.5*P.SlitSep
		m := 0
		for i := 0; i < nb; i++ {
			inwall := math.Abs(s.wr(s.px[i]-P.SlitWallx)) < 0.5*P.SlitTh
			inA := math.Abs(s.py[i]-yA) < P.SlitHw
			inB := math.Abs(s.py[i]-yB) < P.SlitHw
			keep := true
			if P.SlitMask != 3 && inwall {
				keep = false
				if P.SlitMask == 0 && (inA || inB) {
					keep = true
				}
				if P.SlitMask == 1 && inA {
					keep = true
				}
				if P.SlitMask == 2 && inB {
					keep = true
				}
			}
			if keep {
				if m != i {
					s.px[m], s.py[m], s.pz[m] = s.px[i], s.py[i], s.pz[i]
					s.cr0[m], s.cr[m] = s.cr0[i], s.cr[i]
					s.n1x[m], s.n1y[m], s.n1z[m] = s.n1x[i], s.n1y[i], s.n1z[i]
					s.n2x[m], s.n2y[m], s.n2z[m] = s.n2x[i], s.n2y[i], s.n2z[i]
					s.th2[m], s.cbeta[m] = s.th2[i], s.cbeta[i]
					s.Es[m], s.Em[m], s.Ee[m] = s.Es[i], s.Em[i], s.Ee[i]
					s.fa1[m], s.fa2[m] = s.fa1[i], s.fa2[i]
				}
				m++
			}
		}
		fprintf(s.Out, "# CARVE slit wall: removed %d cells (mask=%d wall_x=%.6g th=%.6g slits y=%.2f/%.2f hw=%.6g)\n",
			nb-m, P.SlitMask, P.SlitWallx, P.SlitTh, yA, yB, P.SlitHw)
		nb = m
		s.NC = nb
		// THE WALL IS A FIXTURE: the face cells are pinned (pass D skips
		// them); a carved vacuum gap in the LIVE foam heals in ~15 t.u.
		npin := 0
		if P.SlitMask != 3 {
			for i := 0; i < s.NC; i++ {
				if math.Abs(s.wr(s.px[i]-P.SlitWallx)) < P.SlitPinw {
					s.pin[i] = 1
					npin++
				}
			}
		}
		fprintf(s.Out, "# PIN wall fixture: %d cells within +-%.6g of x=%.6g\n",
			npin, P.SlitPinw, P.SlitWallx)
		// quasi-plane packet: Gaussian sigma_x = P.Sigma, sigma_y = slit_sy
		for i := 0; i < s.NC; i++ {
			dx := s.wr(s.px[i] - P.SlitSrcx)
			dy := s.wr(s.py[i] - cy0)
			g := math.Exp(-dx*dx/(2.0*P.Sigma*P.Sigma) -
				dy*dy/(2.0*P.SlitSy*P.SlitSy))
			if g < 1e-8 {
				continue
			}
			tilt := -(P.Kx * dx)
			s.fa1[i] += math.Sqrt(P.Amp*g) * math.Cos(tilt)
			s.fa2[i] += math.Sqrt(P.Amp*g) * math.Sin(tilt)
			s.Ee[i] = s.fa1[i]*s.fa1[i] + s.fa2[i]*s.fa2[i]
		}
		for i := 0; i < s.NC; i++ {
			s.normalsTransverse(i, 1, 0, 0)
		}
		if P.SlitObj != 0 {
			// MOTION #29/XSEC: dense object in the beam; obj_y (<0 = box
			// centre) sets the impact parameter b = obj_y - cy0
			nobj := 0
			oy := cy0
			if P.ObjY >= 0 {
				oy = P.ObjY
			}
			for i := 0; i < s.NC; i++ {
				dx := s.wr(s.px[i] - cx0)
				dy := s.wr(s.py[i] - oy)
				add := P.ObjAmp * math.Exp(-(dx*dx+dy*dy)/(2.0*P.ObjSigma*P.ObjSigma))
				if add < 1e-4 {
					continue
				}
				if s.Em[i]+add > 0.95*P.Cap {
					add = 0.95*P.Cap - s.Em[i]
				}
				if add <= 0 {
					continue
				}
				s.Em[i] += add
				pull := P.SPull * add
				avail := s.Es[i] - P.EsFloor
				if pull > avail {
					if avail > 0 {
						pull = avail
					} else {
						pull = 0
					}
				}
				s.Es[i] -= pull
				s.Em[i] += pull
				if add > 0.05*P.ObjAmp {
					s.tag[i] = 1
					nobj++
				}
			}
			fprintf(s.Out, "# SEED slit_obj: amp=%.6g sigma=%.6g y=%.2f tagged=%d (dense occulter)\n",
				P.ObjAmp, P.ObjSigma, oy, nobj)
		}
		if P.SlitClicks != 0 {
			nsc := 0
			for i := 0; i < s.NC; i++ {
				if math.Abs(s.wr(s.px[i]-P.SlitScreenx)) <= 0.75 {
					s.scond[i] = 1
					nsc++
				}
			}
			fprintf(s.Out, "# SEED slit_clicks: %d screen cells condensation-active (e_cond=0 there)\n", nsc)
		}
		lam := 0.0
		if P.Kx != 0 {
			lam = TwoPi / P.Kx
		}
		D := P.SlitScreenx - P.SlitWallx
		fprintf(s.Out, "# SEED slit: src_x=%.6g sigma_x=%.6g sigma_y=%.6g kx=%.6g amp=%.6g -> lam_seed=%.4f\n",
			P.SlitSrcx, P.Sigma, P.SlitSy, P.Kx, P.Amp, lam)
		fprintf(s.Out, "# SEED slit loci: D=%.2f s=%.2f dy_smallangle=%.3f screen_x=%.6g gate=[%.6g,%.6g]\n",
			D, P.SlitSep, lam*D/P.SlitSep, P.SlitScreenx, P.SlitT0, P.SlitT1)
	}
	// exp=bath: nothing to seed

	if P.NoiseAmp > 0 {
		for i := 0; i < s.NC; i++ {
			e0 := s.frand() * P.NoiseAmp
			ph := s.frand() * TwoPi
			s.fa1[i] += math.Sqrt(e0) * math.Cos(ph)
			s.fa2[i] += math.Sqrt(e0) * math.Sin(ph)
			s.Ee[i] = s.fa1[i]*s.fa1[i] + s.fa2[i]*s.fa2[i]
		}
		fprintf(s.Out, "# SEED noise: amp=%.6g (e1-class churn)\n", P.NoiseAmp)
	}

	if P.TagR > 0 {
		// FORGE: region tag WITHOUT matter (convtag region ledger)
		nt := 0
		for i := 0; i < s.NC; i++ {
			dx := s.wr(s.px[i] - cx0)
			dy := s.wr(s.py[i] - cy0)
			dz := s.wr(s.pz[i] - cz0)
			if dx*dx+dy*dy+dz*dz < P.TagR*P.TagR {
				s.tag[i] = 1
				nt++
			}
		}
		fprintf(s.Out, "# TAG region: %d cells within r=%.6g of centre (no load added)\n",
			nt, P.TagR)
	}

	// v91 cantus init: holdings memory starts AT the seeded part;
	// cant_seed=1 arms tagged object voices at full amplitude
	for i := 0; i < s.NC; i++ {
		s.cxl[i] = s.Em[i] / P.Cap
	}

	// initial radii + topology + ledger
	for i := 0; i < s.NC; i++ {
		ratio := s.Es[i] / P.Es0
		if ratio < 0 {
			ratio = 0
		}
		s.cr[i] = s.cr0[i] * math.Cbrt(ratio)
	}
	s.topoRefresh()
	// v91 cantus v1.1: cant_seed=1 arms the BONDS between tagged
	// voices at full gauge memory (apparatus; default = self-growth)
	if P.CantSeed != 0 {
		nseed := 0
		for sl := 0; sl < s.NSLOT; sl++ {
			if s.sst[sl] == sFree {
				continue
			}
			if s.tag[s.sli[sl]] != 0 && s.tag[s.slj[sl]] != 0 {
				s.sgg[sl] = 1.0
				nseed++
			}
		}
		fprintf(s.Out, "# SEED cantus: sgg=1 on %d tagged-pair slots\n", nseed)
	}
	for i := 0; i < s.NC; i++ {
		s.fxb[i], s.fyb[i], s.fzb[i] = 0, 0, 0
	}
	s.E0 = s.totalEnergy()
	s.births, s.deaths, s.betaReturns = 0, 0, 0
	s.betaEnergy = 0
	{
		phi, zl, nla, _, dbar, sig, _ := s.geoStatsV()
		fprintf(s.Out, "# INIT NC=%d slots=%d live=%d phi_nom=%.4f z_live=%.2f dbar=%.4f sigma_d=%.2f%% E0=%.9f\n",
			s.NC, s.NSLOT, nla, phi, zl, dbar, 100*sig, s.E0)
	}
	if s.trussN > 0 {
		s.gyrationEigs(&s.trussShape0)
	}

	if P.SectMeter != 0 {
		s.sectCx = cx0
		if P.SectX >= 0 {
			s.sectCx = P.SectX
		}
		s.sectCy = cy0
		if P.SectY >= 0 {
			s.sectCy = P.SectY
		}
		if P.SectN > sectMax {
			P.SectN = sectMax
		}
		if P.SectN < 1 {
			P.SectN = 1
		}
		fprintf(s.Out, "# SECT meter: centre=(%.2f,%.2f) r=[%.6g,%.6g) n=%d gate=[%.6g,%.6g]\n",
			s.sectCx, s.sectCy, P.SectR0, P.SectR1, P.SectN, P.SectT0, P.SectT1)
	}

	var stream *fcsW
	if P.SnapFile != "" {
		var err error
		stream, err = fcsOpen(P.SnapFile, P.SnapComp != 0)
		if err != nil {
			fprintf(s.Errw, "# WARN snap_file open: %v\n", err)
		} else {
			s.fcsBegin(stream)
		}
	}

	// ---------------- run ----------------
	steps := int(P.T/P.Dt + 0.5)
	sheared := false
	fprintf(s.Out, "#\n# t | drift_rel | phi_nom z_live NLa NLd births deaths | dbar sigma_d%% | maxstep maxbond | Em_tag ret com_dev rms conn z_tag | pair: d delta gf gb | core/far\n")
	emt0 := -1.0
	snapIdx := 0
	for st := 0; st <= steps; st++ {
		s.simT = float64(st) * P.Dt
		// instantaneous deviatoric strain test (truss shear leg)
		if P.ShearEps > 0 && !sheared && s.simT >= P.ShearT && P.ShearT > 0 {
			sx := 1.0 + P.ShearEps
			sy := 1.0 / math.Sqrt(1.0+P.ShearEps)
			sz := sy
			for i := 0; i < s.NC; i++ {
				if s.tag[i] == 0 {
					continue
				}
				s.px[i] = s.fold(cx0 + s.wr(s.px[i]-cx0)*sx)
				s.py[i] = s.fold(cy0 + s.wr(s.py[i]-cy0)*sy)
				s.pz[i] = s.fold(cz0 + s.wr(s.pz[i]-cz0)*sz)
			}
			sheared = true
			fprintf(s.Out, "# SHEAR applied: eps=%.6g at t=%.2f (volume-preserving deviatoric)\n",
				P.ShearEps, s.simT)
		}
		if st%P.DiagEvery == 0 || st == steps {
			Etot := s.totalEnergy()
			den := s.E0
			if den == 0 {
				den = 1
			}
			drift := (Etot - s.E0) / den
			phi, zl, nla, nld, dbar, sig, mld := s.geoStatsV()
			emt, comx, comy, comz, rms, conn, zb := s.tagStats()
			if emt0 < 0 && emt > 0 {
				emt0 = emt
			}
			ret := 0.0
			if emt0 > 0 {
				ret = emt / emt0
			}
			comdev := 0.0
			if emt > 1e-12 {
				comdev = math.Sqrt((comx-s.cenx)*(comx-s.cenx) + (comy-s.ceny)*(comy-s.ceny) +
					(comz-s.cenz)*(comz-s.cenz))
			}
			mstep := s.maxNetForce()
			fprintf(s.Out, "t=%8.2f %+9.2e | %.4f %5.2f %6d %4d %6d %6d | %.4f %5.2f | %8.2e %8.2e | %8.4f %6.4f %7.4f %7.4f %5.3f %5.2f",
				s.simT, drift, phi, zl, nla, nld, s.births, s.deaths,
				dbar, 100*sig, mstep, mld, emt, ret, comdev, rms, conn, zb)
			if s.repPairI >= 0 {
				i, j := s.repPairI, s.repPairJ
				lo, hi := i, j
				if j < i {
					lo, hi = j, i
				}
				sl := s.hfind(lo, hi)
				if sl >= 0 && s.sst[sl] != sFree {
					delta := wrapPi((float64(s.slq[sl])*s.w2e[i] + float64(s.slp[sl])*s.w2e[j]) * s.sd[sl] / P.C)
					gf := s.gateOf(float64(s.slq[sl])*s.th2[i] - float64(s.slq[sl])*s.w2e[i]*s.sd[sl]/P.C - float64(s.slp[sl])*s.th2[j])
					gb := s.gateOf(float64(s.slp[sl])*s.th2[j] - float64(s.slp[sl])*s.w2e[j]*s.sd[sl]/P.C - float64(s.slq[sl])*s.th2[i])
					fprintf(s.Out, " | d=%.6f dl=%+.4f gf=%.4f gb=%.4f x=%.4f",
						s.sd[sl], delta, gf, gb, (s.Em[i]+s.flload[i])/P.Cap)
				} else {
					fprintf(s.Out, " | PAIR-CHANNEL-DEAD")
				}
			}
			if P.Exp == "blob" {
				var sh [8]float64
				s.esShells(&sh)
				far := (sh[5] + sh[6] + sh[7]) / 3.0
				cf := 0.0
				if far > 0 {
					cf = sh[0] / far
				}
				fprintf(s.Out, " | es_core/far=%.4f", cf)
				if s.nfsamp < 512 {
					s.fsT[s.nfsamp] = s.simT
					s.fsR[s.nfsamp] = comx
					s.fsY[s.nfsamp] = comy
					s.fsZ[s.nfsamp] = comz
					s.nfsamp++
				}
			}
			if P.Exp == "blob2" {
				ax, aw, bx2, bw, tx, tw := 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
				for i2 := 0; i2 < s.NC; i2++ {
					w := s.Em[i2] + s.flload[i2]
					xx := s.cenx + s.wr(s.px[i2]-s.cenx)
					if w > 1e-12 {
						tx += w * xx
						tw += w
					}
					if s.tag[i2] == 1 && w > 1e-12 {
						ax += w * xx
						aw += w
					}
					if s.tag[i2] == 2 && w > 1e-12 {
						bx2 += w * xx
						bw += w
					}
				}
				comA, comB := 0.0, 0.0
				if aw > 0 {
					comA = ax / aw
				}
				if bw > 0 {
					comB = bx2 / bw
				}
				tcom := 0.0
				if tw > 0 {
					tcom = tx / tw
				}
				fprintf(s.Out, " | A=%.4f B=%.4f sepx=%.4f tot=%.5f EmA=%.3f EmB=%.3f",
					comA, comB, comB-comA, tcom, aw, bw)
				if s.nfsamp < 512 {
					s.fsT[s.nfsamp] = s.simT
					s.b2ax[s.nfsamp] = comA
					s.b2bx[s.nfsamp] = comB
					s.b2tx[s.nfsamp] = tcom
					s.nfsamp++
				}
			}
			if P.Exp == "pulse" {
				// field front = max radius over cells with Ee > 5% amp
				rf := 0.0
				thr := 0.05
				if P.Amp > 0 {
					thr = 0.05 * P.Amp
				}
				for i2 := 0; i2 < s.NC; i2++ {
					if s.Ee[i2] <= thr {
						continue
					}
					dx := s.wr(s.px[i2] - s.cenx)
					dy := s.wr(s.py[i2] - s.ceny)
					dz := s.wr(s.pz[i2] - s.cenz)
					rr := math.Sqrt(dx*dx + dy*dy + dz*dz)
					if rr > rf {
						rf = rr
					}
				}
				fprintf(s.Out, " | Ee_front=%.3f", rf)
				if s.nfsamp < 512 {
					s.fsT[s.nfsamp] = s.simT
					s.fsR[s.nfsamp] = rf
					s.nfsamp++
				}
			}
			if s.trussN > 0 {
				mx, sm, gmean := 0.0, 0.0, 0.0
				for e := 0; e < s.trussN; e++ {
					i, j := s.trussE[e][0], s.trussE[e][1]
					dx := s.px[j] - s.px[i]
					dy := s.py[j] - s.py[i]
					dz := s.pz[j] - s.pz[i]
					d := math.Sqrt(dx*dx + dy*dy + dz*dz)
					dev := math.Abs(d - s.trussDstar[e])
					if dev > mx {
						mx = dev
					}
					sm += dev
					lo, hi := i, j
					if j < i {
						lo, hi = j, i
					}
					sl := s.hfind(lo, hi)
					if sl >= 0 && s.sst[sl] != sFree {
						gf := s.gateOf(float64(s.slq[sl])*s.th2[i] - float64(s.slq[sl])*s.w2e[i]*s.sd[sl]/P.C - float64(s.slp[sl])*s.th2[j])
						gb := s.gateOf(float64(s.slp[sl])*s.th2[j] - float64(s.slp[sl])*s.w2e[j]*s.sd[sl]/P.C - float64(s.slq[sl])*s.th2[i])
						gmean += gf * gb
					}
				}
				var eg [3]float64
				s.gyrationEigs(&eg)
				sdev := math.Abs(eg[0]-s.trussShape0[0]) + math.Abs(eg[1]-s.trussShape0[1]) +
					math.Abs(eg[2]-s.trussShape0[2])
				xm := 0.0
				nxm := 0
				for i2 := 0; i2 < s.NC; i2++ {
					if s.tag[i2] != 0 {
						xm += (s.Em[i2] + s.flload[i2]) / P.Cap
						nxm++
					}
				}
				xbar := 0.0
				if nxm > 0 {
					xbar = xm / float64(nxm)
				}
				fprintf(s.Out, " | edge_dev mean=%.4f max=%.4f gg=%.4f shape_dev=%.4f xbar=%.4f",
					sm/float64(s.trussN), mx, gmean/float64(s.trussN), sdev, xbar)
			}
			if s.rintN > 0 {
				fprintf(s.Out, " | IR")
				for e := 0; e < s.rintN; e++ {
					i, j := s.rintE[e][0], s.rintE[e][1]
					lo, hi := i, j
					if j < i {
						lo, hi = j, i
					}
					sl := s.hfind(lo, hi)
					if sl < 0 || s.sst[sl] == sFree {
						fprintf(s.Out, " e%d=DEAD", e)
						continue
					}
					p := float64(s.rintP[e])
					q := float64(s.rintQ[e])
					psi := wrapPi(q*s.th2[i] - q*s.w2e[i]*s.sd[sl]/P.C - p*s.th2[j])
					gf := s.gateOf(q*s.th2[i] - q*s.w2e[i]*s.sd[sl]/P.C - p*s.th2[j])
					gb := s.gateOf(p*s.th2[j] - p*s.w2e[j]*s.sd[sl]/P.C - q*s.th2[i])
					fprintf(s.Out, " e%d:psi=%+.3f gg=%.3f d=%.4f", e, psi, gf*gb, s.sd[sl])
				}
			}
			if s.triOn {
				for t2 := 0; t2 < s.ntri; t2++ {
					var psi [3]float64
					var ggm, circ float64
					s.triMeters(t2, &psi, &ggm, &circ)
					fprintf(s.Out, " | T%d psi=(%+.3f,%+.3f,%+.3f) ggm=%.3f circ=%+.4f xUDD=(%.3f,%.3f,%.3f)",
						t2, psi[0], psi[1], psi[2], ggm, circ,
						(s.Em[s.triV[t2][0]]+s.flload[s.triV[t2][0]])/P.Cap,
						(s.Em[s.triV[t2][1]]+s.flload[s.triV[t2][1]])/P.Cap,
						(s.Em[s.triV[t2][2]]+s.flload[s.triV[t2][2]])/P.Cap)
				}
				if s.ntri == 2 {
					c1x, c1y, c2x, c2y := 0.0, 0.0, 0.0, 0.0
					for kk := 0; kk < 3; kk++ {
						c1x += s.cenx + s.wr(s.px[s.triV[0][kk]]-s.cenx)
						c1y += s.ceny + s.wr(s.py[s.triV[0][kk]]-s.ceny)
						c2x += s.cenx + s.wr(s.px[s.triV[1][kk]]-s.cenx)
						c2y += s.ceny + s.wr(s.py[s.triV[1][kk]]-s.ceny)
					}
					ddx := (c2x - c1x) / 3
					ddy := (c2y - c1y) / 3
					fprintf(s.Out, " | sep=%.4f", math.Sqrt(ddx*ddx+ddy*ddy))
				}
			}
			fprintf(s.Out, "\n")
			if P.P1Meter != 0 {
				tx := s.p1sp[0] + s.p1fl[0] + s.p1fd[0] + s.p1gm[0]
				ty := s.p1sp[1] + s.p1fl[1] + s.p1fd[1] + s.p1gm[1]
				tz := s.p1sp[2] + s.p1fl[2] + s.p1fd[2] + s.p1gm[2]
				fprintf(s.Out, "# P1 t=%.2f tot=(%+.6e,%+.6e,%+.6e) sp=%+.3e fl=%+.3e fd=%+.3e gm=%+.3e (x)\n",
					s.simT, tx, ty, tz, s.p1sp[0], s.p1fl[0], s.p1fd[0], s.p1gm[0])
			}
			if P.SlitObj != 0 || P.Convtag != 0 {
				fprintf(s.Out, "# CONVTAG t=%.2f rough=%.6f cond=%.6f evap=%.6f backs=%.6f\n",
					s.simT, s.ctRough, s.ctCond, s.ctEvap, s.ctBacks)
			}
			if P.KCant > 0 || P.KTune > 0 {
				at, am, xt := 0.0, 0.0, 0.0
				ntg, nl := 0, 0
				for i2 := 0; i2 < s.NC; i2++ {
					av := s.cantAmpOf(i2)
					if av > am {
						am = av
					}
					if av > 0.5 {
						nl++
					}
					if s.tag[i2] != 0 {
						at += av
						xt += s.cxl[i2]
						ntg++
					}
				}
				atv, xtv := 0.0, 0.0
				if ntg > 0 {
					atv = at / float64(ntg)
					xtv = xt / float64(ntg)
				}
				fprintf(s.Out, "# CANT t=%.2f a_tag=%.4f a_max=%.4f xl_tag=%.4f nlock=%d tune=%.6f\n",
					s.simT, atv, am, xtv, nl, s.tuneTotal)
			}
			if P.RegTau > 0 {
				ntp, a25, a50, a75, a90, agr, afl := s.regStats(0)
				nba, b25, b50, b75, b90, bgr, bfl := s.regStats(1)
				fprintf(s.Out, "# REG t=%.2f tp n=%d rho=[%.3f %.3f %.3f %.3f] gr=%.5f fl=%.3f | ba n=%d rho=[%.3f %.3f %.3f %.3f] gr=%.5f fl=%.3f\n",
					s.simT, ntp, a25, a50, a75, a90, agr, afl,
					nba, b25, b50, b75, b90, bgr, bfl)
			}
		}
		if P.SnapEvery > 0 && s.P.SnapDir != "" && st%P.SnapEvery == 0 {
			s.writeFCS(snapIdx)
			snapIdx++
		}
		if stream != nil && P.SnapEvery > 0 && st%P.SnapEvery == 0 {
			s.fcsCellFrame(stream)
			s.fcsInstrument(stream)
		}
		if st == steps {
			break
		}
		s.step()
		if P.Exp == "slit" {
			tn := float64(st+1) * P.Dt
			if tn >= P.SlitT0 && tn <= P.SlitT1 {
				if s.dsCellI == nil {
					s.dsCellI = make([]float64, s.NC)
				}
				for i := 0; i < s.NC; i++ {
					if math.Abs(s.wr(s.px[i]-P.SlitScreenx)) > 0.75 {
						continue
					}
					b := int(s.py[i] / P.L * dsNBin)
					if b < 0 {
						b = 0
					}
					if b >= dsNBin {
						b = dsNBin - 1
					}
					s.dsI[b] += s.Ee[i] * P.Dt
					s.dsCellI[i] += s.Ee[i] * P.Dt
					s.dsExpo += s.Ee[i] * P.Dt
				}
			}
		}
		if P.SectMeter != 0 {
			tn := float64(st+1) * P.Dt
			if tn >= P.SectT0 && tn <= P.SectT1 {
				NS := P.SectN
				for i := 0; i < s.NC; i++ {
					dx := s.wr(s.px[i] - s.sectCx)
					dy := s.wr(s.py[i] - s.sectCy)
					r2 := dx*dx + dy*dy
					if r2 < P.SectR0*P.SectR0 || r2 >= P.SectR1*P.SectR1 {
						continue
					}
					// half-bin rotation: sector k is CENTRED at k*2pi/NS
					u := (math.Atan2(dy, dx) + math.Pi/float64(NS)) / TwoPi
					u -= math.Floor(u)
					k := int(u * float64(NS))
					if k >= NS {
						k = NS - 1
					}
					s.sectE[k] += s.Ee[i] * P.Dt
					s.sectN2[k] += P.Dt
				}
			}
		}
	}

	// ---------------- final report ----------------
	{
		Etot := s.totalEnergy()
		if stream != nil {
			den := s.E0
			if den == 0 {
				den = 1
			}
			s.fcsAnlzText(stream, fmt.Sprintf(
				"RESULT drift_rel %.3e\nRESULT births %d deaths %d beta_returns %d\n"+
					"RESULT conv rough=%.6f cond=%.6f evap=%.6f backs=%.6f rad=%.6f\n",
				(Etot-s.E0)/den, s.births, s.deaths, s.betaReturns,
				s.roughTotal, s.condTotal, s.evapTotal, s.backsTotal, s.radTotal))
			stream.close()
			stream = nil
		}
		den := s.E0
		if den == 0 {
			den = 1
		}
		fprintf(s.Out, "#\n# RESULT drift_rel %.3e\n", (Etot-s.E0)/den)
		fprintf(s.Out, "# RESULT births %d deaths %d beta_returns %d beta_energy %.6e\n",
			s.births, s.deaths, s.betaReturns, s.betaEnergy)
		fprintf(s.Out, "# RESULT conv rough=%.6f cond=%.6f evap=%.6f backs=%.6f rad=%.6f\n",
			s.roughTotal, s.condTotal, s.evapTotal, s.backsTotal, s.radTotal)
		{ // pad-31 apparatus: the standing sub-atom credit inventory
			cd, cf := 0.0, 0.0
			for ii := 0; ii < s.NC; ii++ {
				cd += s.qcnvD[ii]
				cf += s.qcnvF[ii]
			}
			fprintf(s.Out, "# RESULT credit qcnvD=%.6f qcnvF=%.6f\n", cd, cf)
		}
		if P.KCant > 0 || P.KTune > 0 {
			at, am, xt := 0.0, 0.0, 0.0
			ntg, nl := 0, 0
			for ii := 0; ii < s.NC; ii++ {
				av := s.cantAmpOf(ii)
				if av > am {
					am = av
				}
				if av > 0.5 {
					nl++
				}
				if s.tag[ii] != 0 {
					at += av
					xt += s.cxl[ii]
					ntg++
				}
			}
			atv, xtv := 0.0, 0.0
			if ntg > 0 {
				atv = at / float64(ntg)
				xtv = xt / float64(ntg)
			}
			fprintf(s.Out, "# RESULT cantus a_tag=%.4f a_max=%.4f xl_tag=%.4f nlock=%d tune_total=%.6f\n",
				atv, am, xtv, nl, s.tuneTotal)
		}
		if P.RegTau > 0 {
			ntp, a25, a50, a75, a90, agr, afl := s.regStats(0)
			nba, b25, b50, b75, b90, bgr, bfl := s.regStats(1)
			fprintf(s.Out, "# RESULT reg tp n=%d rho=[%.4f %.4f %.4f %.4f] gr=%.6f fl=%.4f | ba n=%d rho=[%.4f %.4f %.4f %.4f] gr=%.6f fl=%.4f\n",
				ntp, a25, a50, a75, a90, agr, afl,
				nba, b25, b50, b75, b90, bgr, bfl)
		}
		if P.SlitObj != 0 || P.Convtag != 0 {
			// net field capture at the occulter = cond - evap - rough + backs
			fprintf(s.Out, "# RESULT convtag rough=%.6f cond=%.6f evap=%.6f backs=%.6f net=%.6f\n",
				s.ctRough, s.ctCond, s.ctEvap, s.ctBacks,
				s.ctCond-s.ctEvap-s.ctRough+s.ctBacks)
		}
		if P.P1Meter != 0 {
			tx := s.p1sp[0] + s.p1fl[0] + s.p1fd[0] + s.p1gm[0]
			ty := s.p1sp[1] + s.p1fl[1] + s.p1fd[1] + s.p1gm[1]
			tz := s.p1sp[2] + s.p1fl[2] + s.p1fd[2] + s.p1gm[2]
			fprintf(s.Out, "# RESULT p1 tot=(%+.6e,%+.6e,%+.6e) rate=(%+.6e,%+.6e,%+.6e)\n",
				tx, ty, tz, tx/P.T, ty/P.T, tz/P.T)
			fprintf(s.Out, "# RESULT p1x sp=%+.6e fl=%+.6e fd=%+.6e gm=%+.6e\n",
				s.p1sp[0], s.p1fl[0], s.p1fd[0], s.p1gm[0])
		}
		phi, zl, _, _, dbar, sig, _ := s.geoStatsV()
		fprintf(s.Out, "# RESULT geometry phi_nom=%.4f z_live=%.2f dbar=%.4f sigma_d=%.4f\n",
			phi, zl, dbar, sig)
		if s.repPairI >= 0 {
			i, j := s.repPairI, s.repPairJ
			lo, hi := i, j
			if j < i {
				lo, hi = j, i
			}
			sl := s.hfind(lo, hi)
			if sl >= 0 && s.sst[sl] != sFree {
				x1r := P.PairX1
				if x1r < 0 {
					x1r = P.PairX0
				}
				wA := P.W2 / (1.0 + P.QDetune*P.PairX0)
				wB := P.W2 / (1.0 + P.QDetune*x1r)
				ds := TwoPi * float64(P.PairM) * P.C / (float64(P.PairQQ)*wA + float64(P.PairPP)*wB)
				dsl := TwoPi * float64(P.PairM) * P.C / (float64(s.slq[sl])*s.w2e[i] + float64(s.slp[sl])*s.w2e[j])
				fprintf(s.Out, "# RESULT pair d_final=%.6f d_star_seed=%.6f d_star_live=%.6f off_live=%+.6f Em_i=%.5f Em_j=%.5f\n",
					s.sd[sl], ds, dsl, s.sd[sl]-dsl, s.Em[i], s.Em[j])
				// directed ledger (PAULI-0 rate control) — C kernel parity
				fprintf(s.Out, "# RESULT pairflux dep_ij=%.6e dep_ji=%.6e lem_ij=%.6e lem_ji=%.6e\n",
					s.sfluxd[2*sl], s.sfluxd[2*sl+1], s.slem[2*sl], s.slem[2*sl+1])
			} else {
				fprintf(s.Out, "# RESULT pair CHANNEL DEAD\n")
			}
		}
		if s.trussN > 0 {
			mx, sm := 0.0, 0.0
			for e := 0; e < s.trussN; e++ {
				i, j := s.trussE[e][0], s.trussE[e][1]
				dx := s.px[j] - s.px[i]
				dy := s.py[j] - s.py[i]
				dz := s.pz[j] - s.pz[i]
				d := math.Sqrt(dx*dx + dy*dy + dz*dz)
				dev := math.Abs(d - s.trussDstar[e])
				if dev > mx {
					mx = dev
				}
				sm += dev
			}
			var eg [3]float64
			s.gyrationEigs(&eg)
			fprintf(s.Out, "# RESULT truss edge_dev_mean=%.6f edge_dev_max=%.6f shape=(%.4f,%.4f,%.4f) shape0=(%.4f,%.4f,%.4f)\n",
				sm/float64(s.trussN), mx, eg[0], eg[1], eg[2],
				s.trussShape0[0], s.trussShape0[1], s.trussShape0[2])
			for e := 0; e < s.trussN; e++ {
				i, j := s.trussE[e][0], s.trussE[e][1]
				lo, hi := i, j
				if j < i {
					lo, hi = j, i
				}
				sl := s.hfind(lo, hi)
				if sl < 0 || s.sst[sl] == sFree {
					fprintf(s.Out, "# EDGE %d %d-%d CHANNEL DEAD\n", e, i, j)
					continue
				}
				gf := s.gateOf(float64(s.slq[sl])*s.th2[i] - float64(s.slq[sl])*s.w2e[i]*s.sd[sl]/P.C - float64(s.slp[sl])*s.th2[j])
				gb := s.gateOf(float64(s.slp[sl])*s.th2[j] - float64(s.slp[sl])*s.w2e[j]*s.sd[sl]/P.C - float64(s.slq[sl])*s.th2[i])
				d1i := s.n1x[i]*s.sux[sl] + s.n1y[i]*s.suy[sl] + s.n1z[i]*s.suz[sl]
				d2i := s.n2x[i]*s.sux[sl] + s.n2y[i]*s.suy[sl] + s.n2z[i]*s.suz[sl]
				d1j := s.n1x[j]*s.sux[sl] + s.n1y[j]*s.suy[sl] + s.n1z[j]*s.suz[sl]
				d2j := s.n2x[j]*s.sux[sl] + s.n2y[j]*s.suy[sl] + s.n2z[j]*s.suz[sl]
				axi := (1.0 - d1i*d1i) * (1.0 - d2i*d2i)
				if axi < 0 {
					axi = 0
				}
				axj := (1.0 - d1j*d1j) * (1.0 - d2j*d2j)
				if axj < 0 {
					axj = 0
				}
				dnn := s.n2x[i]*s.n2x[j] + s.n2y[i]*s.n2y[j] + s.n2z[i]*s.n2z[j]
				fprintf(s.Out, "# EDGE %d %d-%d d=%.5f dev=%+.5f gg=%.4f gpl=%.4f\n",
					e, i, j, s.sd[sl], s.sd[sl]-s.trussDstar[e], gf*gb,
					axi*axj*dnn*dnn)
			}
		}
		if s.rintN > 0 {
			for e := 0; e < s.rintN; e++ {
				i, j := s.rintE[e][0], s.rintE[e][1]
				lo, hi := i, j
				if j < i {
					lo, hi = j, i
				}
				sl := s.hfind(lo, hi)
				if sl < 0 || s.sst[sl] == sFree {
					fprintf(s.Out, "# RESULT rint%d %d-%d CHANNEL DEAD\n", e, i, j)
					continue
				}
				p := float64(s.rintP[e])
				q := float64(s.rintQ[e])
				psi := wrapPi(q*s.th2[i] - q*s.w2e[i]*s.sd[sl]/P.C - p*s.th2[j])
				gf := s.gateOf(q*s.th2[i] - q*s.w2e[i]*s.sd[sl]/P.C - p*s.th2[j])
				gb := s.gateOf(p*s.th2[j] - p*s.w2e[j]*s.sd[sl]/P.C - q*s.th2[i])
				fprintf(s.Out, "# RESULT rint%d %d-%d p:q=%d:%d d=%.6f dev=%+.6f psi=%+.4f gg=%.4f lem=%.6e\n",
					e, i, j, s.rintP[e], s.rintQ[e], s.sd[sl], s.sd[sl]-s.rintDstar[e],
					psi, gf*gb, s.slem[2*sl]+s.slem[2*sl+1])
			}
			for k := 0; k < s.ringsNring; k++ {
				ggs, xm, ems := 0.0, 0.0, 0.0
				ne := 0
				for e := 0; e < s.trussN; e++ {
					i := s.trussE[e][0]
					if i < s.ringsV0[k] || i >= s.ringsV0[k]+P.RingsNv {
						continue
					}
					j := s.trussE[e][1]
					lo, hi := i, j
					if j < i {
						lo, hi = j, i
					}
					sl := s.hfind(lo, hi)
					if sl >= 0 && s.sst[sl] != sFree {
						gf := s.gateOf(float64(s.slq[sl])*s.th2[i] - float64(s.slq[sl])*s.w2e[i]*s.sd[sl]/P.C - float64(s.slp[sl])*s.th2[j])
						gb := s.gateOf(float64(s.slp[sl])*s.th2[j] - float64(s.slp[sl])*s.w2e[j]*s.sd[sl]/P.C - float64(s.slq[sl])*s.th2[i])
						ggs += gf * gb
						ne++
					}
				}
				for q := 0; q < P.RingsNv; q++ {
					i := s.ringsV0[k] + q
					xm += (s.Em[i] + s.flload[i]) / P.Cap
					ems += s.Em[i]
				}
				ggi := 0.0
				if ne > 0 {
					ggi = ggs / float64(ne)
				}
				ms := "DROPLET"
				if s.ringsMode[k] != 0 {
					ms = "MOLECULE"
				}
				fprintf(s.Out, "# RESULT ringq%d mode=%s w=%.4f gg_int=%.4f xbar=%.4f Em=%.5f\n",
					k, ms, s.ringsW[k], ggi, xm/float64(P.RingsNv), ems)
			}
		}
		if s.triOn {
			for t2 := 0; t2 < s.ntri; t2++ {
				var psi [3]float64
				var ggm, circ float64
				s.triMeters(t2, &psi, &ggm, &circ)
				fprintf(s.Out, "# RESULT tri%d psi=(%+.4f,%+.4f,%+.4f) ggm=%.4f circ=%+.5f Em=(%.4f,%.4f,%.4f)\n",
					t2, psi[0], psi[1], psi[2], ggm, circ,
					s.Em[s.triV[t2][0]], s.Em[s.triV[t2][1]], s.Em[s.triV[t2][2]])
			}
		}
		if P.Exp == "pulse" && s.nfsamp >= 6 {
			h := s.nfsamp / 2
			n := s.nfsamp - h
			sx, sy, sxx, sxy := 0.0, 0.0, 0.0, 0.0
			for k := h; k < s.nfsamp; k++ {
				sx += s.fsT[k]
				sy += s.fsR[k]
				sxx += s.fsT[k] * s.fsT[k]
				sxy += s.fsT[k] * s.fsR[k]
			}
			den2 := float64(n)*sxx - sx*sx
			v := 0.0
			if den2 != 0 {
				v = (float64(n)*sxy - sx*sy) / den2
			}
			fprintf(s.Out, "# RESULT front_speed v=%.4f v_over_C=%.4f (second-half fit, n=%d)\n",
				v, v/P.C, n)
		}
		if P.Exp == "blob2" && s.nfsamp >= 6 {
			h := s.nfsamp / 2
			// approach speeds from the FIRST half (before contact)
			sx, sxx := 0.0, 0.0
			sya, syb, syt := 0.0, 0.0, 0.0
			sxa, sxb, sxt := 0.0, 0.0, 0.0
			for k2 := 0; k2 < h; k2++ {
				sx += s.fsT[k2]
				sxx += s.fsT[k2] * s.fsT[k2]
				sya += s.b2ax[k2]
				sxa += s.fsT[k2] * s.b2ax[k2]
				syb += s.b2bx[k2]
				sxb += s.fsT[k2] * s.b2bx[k2]
				syt += s.b2tx[k2]
				sxt += s.fsT[k2] * s.b2tx[k2]
			}
			den2 := float64(h)*sxx - sx*sx
			vA, vB, vT := 0.0, 0.0, 0.0
			if den2 != 0 {
				vA = (float64(h)*sxa - sx*sya) / den2
				vB = (float64(h)*sxb - sx*syb) / den2
				vT = (float64(h)*sxt - sx*syt) / den2
			}
			fprintf(s.Out, "# RESULT blob2 vA=%.6f vB=%.6f closing=%.6f vTotalCOM=%.2e sep_final=%.4f\n",
				vA, vB, vA-vB, vT, s.b2bx[s.nfsamp-1]-s.b2ax[s.nfsamp-1])
		}
		if P.SlitClicks != 0 {
			esum := 0.0
			for k2 := 0; k2 < s.nclick; k2++ {
				esum += s.clickE[k2]
			}
			fprintf(s.Out, "# RESULT clicks n=%d e_sum=%.4f\n", s.nclick, esum)
			for k2 := 0; k2 < s.nclick; k2++ {
				fprintf(s.Out, "# CLICK t=%.2f y=%.4f e=%.5f\n",
					s.clickT[k2], s.clickY[k2], s.clickE[k2])
			}
		}
		if P.Exp == "blob" {
			var sh [8]float64
			s.esShells(&sh)
			fprintf(s.Out, "# RESULT es_shells %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",
				sh[0], sh[1], sh[2], sh[3], sh[4], sh[5], sh[6], sh[7])
			emt, comx, comy, comz, rms, conn, zb := s.tagStats()
			fprintf(s.Out, "# RESULT blob Em_tag=%.5f com=(%.3f,%.3f,%.3f) rms=%.4f conn=%.4f z_tag=%.2f\n",
				emt, comx, comy, comz, rms, conn, zb)
			if s.nfsamp >= 6 {
				h := s.nfsamp / 2
				n := s.nfsamp - h
				sx, sxx := 0.0, 0.0
				syx, syy, syz := 0.0, 0.0, 0.0
				sxyx, sxyy, sxyz := 0.0, 0.0, 0.0
				for k := h; k < s.nfsamp; k++ {
					sx += s.fsT[k]
					sxx += s.fsT[k] * s.fsT[k]
					syx += s.fsR[k]
					sxyx += s.fsT[k] * s.fsR[k]
					syy += s.fsY[k]
					sxyy += s.fsT[k] * s.fsY[k]
					syz += s.fsZ[k]
					sxyz += s.fsT[k] * s.fsZ[k]
				}
				den2 := float64(n)*sxx - sx*sx
				vx, vy, vz := 0.0, 0.0, 0.0
				if den2 != 0 {
					vx = (float64(n)*sxyx - sx*syx) / den2
					vy = (float64(n)*sxyy - sx*syy) / den2
					vz = (float64(n)*sxyz - sx*syz) / den2
				}
				sp := math.Sqrt(vx*vx + vy*vy + vz*vz)
				cosk := 0.0
				if sp > 0 {
					cosk = vx / sp
				}
				fprintf(s.Out, "# RESULT blob_drift speed=%.6f cos_to_kdir=%.4f v=(%.2e,%.2e,%.2e)\n",
					sp, cosk, vx, vy, vz)
			}
		}
		if P.SectMeter != 0 {
			stot := 0.0
			for k := 0; k < P.SectN; k++ {
				stot += s.sectE[k]
			}
			fprintf(s.Out, "# RESULT sect Etot=%.6f n=%d r0=%.6g r1=%.6g centre=(%.2f,%.2f) gate=[%.6g,%.6g]\n",
				stot, P.SectN, P.SectR0, P.SectR1, s.sectCx, s.sectCy,
				P.SectT0, P.SectT1)
			for k := 0; k < P.SectN; k++ {
				thc := float64(k) * TwoPi / float64(P.SectN)
				if thc > math.Pi {
					thc -= TwoPi
				}
				fprintf(s.Out, "# SECT k=%d th=%+.4f E=%.8f n=%.2f\n", k, thc, s.sectE[k], s.sectN2[k])
			}
		}
		if P.Exp == "slit" {
			fprintf(s.Out, "# RESULT ds exposure=%.6f gate=[%.6g,%.6g] screen_x=%.6g nbin=%d\n",
				s.dsExpo, P.SlitT0, P.SlitT1, P.SlitScreenx, dsNBin)
			for b := 0; b < dsNBin; b++ {
				fprintf(s.Out, "# SCREEN y=%.4f I=%.8f\n",
					(float64(b)+0.5)*P.L/dsNBin, s.dsI[b])
			}
			// per-cell record: the analyzer smooths these
			if s.dsCellI != nil {
				for i := 0; i < s.NC; i++ {
					if s.dsCellI[i] > 0 {
						fprintf(s.Out, "# SCREENCELL y=%.4f x=%.4f I=%.8f\n",
							s.py[i], s.px[i], s.dsCellI[i])
					}
				}
			}
		}
	}
}

// geoStatsV is geoStats with the multi-return spelled for run.go use.
func (s *Sim) geoStatsV() (phi, zl float64, nla, nld int, dbar, sigD, maxldd float64) {
	return s.geoStats()
}
