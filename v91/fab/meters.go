package fab

// totals + diagnostics — port of v89/freecell.c lines 1252-1593.

import "math"

func (s *Sim) totalEnergy() float64 {
	sum, comp := 0.0, 0.0
	for i := 0; i < s.NC; i++ {
		v := s.Es[i] + s.Em[i] + s.Ee[i]
		y := v - comp
		t := sum + y
		comp = (t - sum) - y
		sum = t
	}
	for sl := 0; sl < s.NSLOT; sl++ {
		if s.sst[sl] == sFree {
			continue
		}
		v := s.slem[2*sl] + s.slem[2*sl+1]
		y := v - comp
		t := sum + y
		comp = (t - sum) - y
		sum = t
	}
	return sum
}

func (s *Sim) geoStats() (phi, zl float64, nla, nld int, dbar, sigD, maxldd float64) {
	vol := 0.0
	for i := 0; i < s.NC; i++ {
		vol += (4.0 / 3.0) * math.Pi * s.cr[i] * s.cr[i] * s.cr[i]
	}
	phi = vol / (s.P.L * s.P.L * s.P.L)
	var na, nd int64
	sm, sm2 := 0.0, 0.0
	var nl int64
	ml := 0.0
	for sl := 0; sl < s.NSLOT; sl++ {
		if s.sst[sl] == sFree {
			continue
		}
		if s.sst[sl] == sDying {
			nd++
			continue
		}
		if s.sA[sl] > 0 {
			na++
			sm += s.sd[sl]
			sm2 += s.sd[sl] * s.sd[sl]
			nl++
		}
		if math.Abs(s.sldd[sl]) > ml {
			ml = math.Abs(s.sldd[sl])
		}
	}
	den := s.NC
	if den <= 0 {
		den = 1
	}
	zl = 2.0 * float64(na) / float64(den)
	nla, nld = int(na), int(nd)
	m := 0.0
	if nl > 0 {
		m = sm / float64(nl)
	}
	vr := 0.0
	if nl > 0 {
		vr = sm2/float64(nl) - m*m
	}
	if vr < 0 {
		vr = 0
	}
	dbar = m
	if m > 0 {
		sigD = math.Sqrt(vr) / m
	}
	maxldd = ml
	return
}

func (s *Sim) maxNetForce() float64 {
	mx := 0.0
	for i := 0; i < s.NC; i++ {
		f := math.Sqrt(s.fxb[i]*s.fxb[i] + s.fyb[i]*s.fyb[i] + s.fzb[i]*s.fzb[i])
		if f > mx {
			mx = f
		}
	}
	return mx
}

// tagged-object metrics: Em retention, COM, RMS radius, connectivity
func (s *Sim) tagStats() (emt, comx, comy, comz, rms, conn, zblob float64) {
	em, cx2, cy2, cz2, wt := 0.0, 0.0, 0.0, 0.0, 0.0
	nt := 0
	for i := 0; i < s.NC; i++ {
		if s.tag[i] == 0 {
			continue
		}
		nt++
		// dense inventory NOW: store + this cell's live half of in-flight
		fl := 0.0
		for q := s.cls[i]; q < s.cls[i+1]; q++ {
			sl := s.clidx[q]
			fl += s.slem[2*sl] + s.slem[2*sl+1]
		}
		w := s.Em[i] + 0.5*fl
		em += w
		cx2 += w * (s.cenx + s.wr(s.px[i]-s.cenx))
		cy2 += w * (s.ceny + s.wr(s.py[i]-s.ceny))
		cz2 += w * (s.cenz + s.wr(s.pz[i]-s.cenz))
		wt += w
	}
	emt = em
	if wt > 1e-12 {
		cx2 /= wt
		cy2 /= wt
		cz2 /= wt
	}
	comx, comy, comz = cx2, cy2, cz2
	r2 := 0.0
	if wt > 1e-12 {
		for i := 0; i < s.NC; i++ {
			if s.tag[i] == 0 {
				continue
			}
			w := s.Em[i] + s.flload[i]
			dx := s.wr(s.px[i] - cx2)
			dy := s.wr(s.py[i] - cy2)
			dz := s.wr(s.pz[i] - cz2)
			r2 += w * (dx*dx + dy*dy + dz*dz)
		}
	}
	if wt > 1e-12 {
		rms = math.Sqrt(r2 / wt)
	}
	// largest connected component among tagged via live (A>0) links
	if s.comp == nil {
		s.comp = make([]int, s.NC)
		s.stack = make([]int, s.NC)
	}
	for i := 0; i < s.NC; i++ {
		s.comp[i] = -1
	}
	bigc := 0
	var degsum int64
	for i := 0; i < s.NC; i++ {
		if s.tag[i] == 0 {
			continue
		}
		var dg int64
		for q := s.cls[i]; q < s.cls[i+1]; q++ {
			if s.sA[s.clidx[q]] > 0 {
				dg++
			}
		}
		degsum += dg
		if s.comp[i] >= 0 {
			continue
		}
		sz, sp := 0, 0
		s.stack[sp] = i
		sp++
		s.comp[i] = i
		for sp > 0 {
			sp--
			u := s.stack[sp]
			sz++
			for q := s.cls[u]; q < s.cls[u+1]; q++ {
				sl := s.clidx[q]
				if s.sA[sl] <= 0 {
					continue
				}
				v := s.slj[sl]
				if s.sli[sl] != u {
					v = s.sli[sl]
				}
				if s.tag[v] != 0 && s.comp[v] < 0 {
					s.comp[v] = i
					s.stack[sp] = v
					sp++
				}
			}
		}
		if sz > bigc {
			bigc = sz
		}
	}
	if nt > 0 {
		conn = float64(bigc) / float64(nt)
		zblob = float64(degsum) / float64(nt)
	}
	return
}

func (s *Sim) esShells(sh *[8]float64) {
	var cnt [8]float64
	for k := 0; k < 8; k++ {
		sh[k] = 0
	}
	drs := 0.5 * s.P.L / 8.0
	for i := 0; i < s.NC; i++ {
		dx := s.wr(s.px[i] - s.cenx)
		dy := s.wr(s.py[i] - s.ceny)
		dz := s.wr(s.pz[i] - s.cenz)
		r := math.Sqrt(dx*dx + dy*dy + dz*dz)
		k := int(r / drs)
		if k < 8 {
			sh[k] += s.Es[i]
			cnt[k] += 1
		}
	}
	for k := 0; k < 8; k++ {
		if cnt[k] > 0 {
			sh[k] /= cnt[k]
		}
	}
}

// the fixed-chart edge defect psi_e and cumulative circulation around a
// declared 3-cycle
func (s *Sim) triMeters(t int, psi *[3]float64, ggm, circ *float64) {
	*ggm = 1
	*circ = 0
	for e := 0; e < 3; e++ {
		a, b := s.triV[t][e], s.triV[t][(e+1)%3]
		lo, hi := a, b
		if b < a {
			lo, hi = b, a
		}
		sl := s.hfind(lo, hi)
		if sl < 0 || s.sst[sl] == sFree {
			psi[e] = 99
			*ggm = 0
			continue
		}
		p, q := float64(s.triP[t][e]), float64(s.triQ[t][e])
		psi[e] = wrapPi(q*s.th2[a] - q*s.w2e[a]*s.sd[sl]/s.P.C - p*s.th2[b])
		gf := s.gateOf(q*s.th2[a] - q*s.w2e[a]*s.sd[sl]/s.P.C - p*s.th2[b])
		gb := s.gateOf(p*s.th2[b] - p*s.w2e[b]*s.sd[sl]/s.P.C - q*s.th2[a])
		*ggm *= math.Sqrt(gf * gb)
		// cycle-direction flux minus counter-flux: a->b is dir0 iff a==sli
		fwd := 1
		if s.sli[sl] == a {
			fwd = 0
		}
		*circ += s.sfluxd[2*sl+fwd] - s.sfluxd[2*sl+(1-fwd)]
	}
}

func (s *Sim) gyrationEigs(out *[3]float64) {
	// sorted eigenvalues of the tagged gyration tensor (shape metric)
	cx2, cy2, cz2 := 0.0, 0.0, 0.0
	nt := 0
	for i := 0; i < s.NC; i++ {
		if s.tag[i] != 0 {
			cx2 += s.cenx + s.wr(s.px[i]-s.cenx)
			cy2 += s.ceny + s.wr(s.py[i]-s.ceny)
			cz2 += s.cenz + s.wr(s.pz[i]-s.cenz)
			nt++
		}
	}
	if nt == 0 {
		out[0], out[1], out[2] = 0, 0, 0
		return
	}
	cx2 /= float64(nt)
	cy2 /= float64(nt)
	cz2 /= float64(nt)
	var m [6]float64
	for i := 0; i < s.NC; i++ {
		if s.tag[i] == 0 {
			continue
		}
		dx := s.wr(s.px[i] - cx2)
		dy := s.wr(s.py[i] - cy2)
		dz := s.wr(s.pz[i] - cz2)
		m[0] += dx * dx
		m[1] += dy * dy
		m[2] += dz * dz
		m[3] += dx * dy
		m[4] += dx * dz
		m[5] += dy * dz
	}
	for k := 0; k < 6; k++ {
		m[k] /= float64(nt)
	}
	// eigenvalues of symmetric 3x3 (analytic, Smith)
	p1 := m[3]*m[3] + m[4]*m[4] + m[5]*m[5]
	q := (m[0] + m[1] + m[2]) / 3.0
	p2 := (m[0]-q)*(m[0]-q) + (m[1]-q)*(m[1]-q) + (m[2]-q)*(m[2]-q) + 2*p1
	p := math.Sqrt(p2 / 6.0)
	if p < 1e-15 {
		out[0], out[1], out[2] = q, q, q
		return
	}
	var B [6]float64
	B[0] = (m[0] - q) / p
	B[1] = (m[1] - q) / p
	B[2] = (m[2] - q) / p
	B[3] = m[3] / p
	B[4] = m[4] / p
	B[5] = m[5] / p
	detB := B[0]*(B[1]*B[2]-B[5]*B[5]) - B[3]*(B[3]*B[2]-B[5]*B[4]) + B[4]*(B[3]*B[5]-B[1]*B[4])
	r := detB / 2.0
	if r < -1 {
		r = -1
	}
	if r > 1 {
		r = 1
	}
	phi := math.Acos(r) / 3.0
	e1 := q + 2*p*math.Cos(phi)
	e3 := q + 2*p*math.Cos(phi+2.0*math.Pi/3.0)
	e2 := 3*q - e1 - e3
	out[0], out[1], out[2] = e3, e2, e1 // ascending
}
