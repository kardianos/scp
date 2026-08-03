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
	fprintf(s.Out, "# GEOMETRY (apparatus): cfac=%.6g k_rep=%.6g mob_geo=%.6g kappa_bond=%.6g freeze_geo=%d\n",
		P.Cfac, P.KRep, P.MobGeo, P.KappaBond, P.FreezeGeo)
	fprintf(s.Out, "# bath=%d bath_frac=%.6g jam_sweeps=%d jam_k=%.6g L=%.6g dt=%.6g T=%.6g seed=%d\n",
		P.Bath, P.BathFrac, P.JamSweeps, P.JamK, P.L, P.Dt, P.T, P.Seed)

	// ---------------- build ----------------
	V := P.L * P.L * P.L
	nmax := int(0.8*V/(P.Dmin*P.Dmin*P.Dmin)) + 4096
	bx := make([]float64, nmax)
	by := make([]float64, nmax)
	bz := make([]float64, nmax)
	nb := 0
	if P.Bath != 0 {
		nb = s.buildBath(bx, by, bz, nmax)
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

	// initial radii + topology + ledger
	for i := 0; i < s.NC; i++ {
		ratio := s.Es[i] / P.Es0
		if ratio < 0 {
			ratio = 0
		}
		s.cr[i] = s.cr0[i] * math.Cbrt(ratio)
	}
	s.topoRefresh()
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
		}
		if P.SnapEvery > 0 && s.P.SnapDir != "" && st%P.SnapEvery == 0 {
			s.writeFCS(snapIdx)
			snapIdx++
		}
		if st == steps {
			break
		}
		s.step()
	}

	// ---------------- final report ----------------
	{
		Etot := s.totalEnergy()
		den := s.E0
		if den == 0 {
			den = 1
		}
		fprintf(s.Out, "#\n# RESULT drift_rel %.3e\n", (Etot-s.E0)/den)
		fprintf(s.Out, "# RESULT births %d deaths %d beta_returns %d beta_energy %.6e\n",
			s.births, s.deaths, s.betaReturns, s.betaEnergy)
		fprintf(s.Out, "# RESULT conv rough=%.6f cond=%.6f evap=%.6f backs=%.6f\n",
			s.roughTotal, s.condTotal, s.evapTotal, s.backsTotal)
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
	}
}

// geoStatsV is geoStats with the multi-return spelled for run.go use.
func (s *Sim) geoStatsV() (phi, zl float64, nla, nld int, dbar, sigD, maxldd float64) {
	return s.geoStats()
}
