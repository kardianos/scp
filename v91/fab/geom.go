package fab

// geometry: spatial hash + candidate refresh + birth/death (rule α)
// + jam settle — port of v89/freecell.c lines 594-745, 1433-1481.

import "math"

func (s *Sim) binsBuild() {
	rmax := 0.0
	for i := 0; i < s.NC; i++ {
		if s.cr[i] > rmax {
			rmax = s.cr[i]
		}
	}
	cut := s.P.Cfac * 2.0 * rmax
	if cut < 1e-6 {
		cut = 1e-6
	}
	g := int(s.P.L / cut)
	if g < 1 {
		g = 1
	}
	if g > 128 {
		g = 128
	}
	if g != s.binG {
		s.binHead = make([]int, g*g*g)
		s.binG = g
	}
	s.binSz = s.P.L / float64(g)
	for b := 0; b < g*g*g; b++ {
		s.binHead[b] = -1
	}
	for i := 0; i < s.NC; i++ {
		bx, by, bz := int(s.px[i]/s.binSz), int(s.py[i]/s.binSz), int(s.pz[i]/s.binSz)
		if bx >= g {
			bx = g - 1
		}
		if by >= g {
			by = g - 1
		}
		if bz >= g {
			bz = g - 1
		}
		if bx < 0 {
			bx = 0
		}
		if by < 0 {
			by = 0
		}
		if bz < 0 {
			bz = 0
		}
		b := (bx*g+by)*g + bz
		s.binNext[i] = s.binHead[b]
		s.binHead[b] = i
	}
}

func pkey(i, j int) int64 { return int64(i)<<20 | int64(j) }

func (s *Sim) hfind(i, j int) int {
	k := pkey(i, j)
	h := int((uint64(k)*0x9E3779B97F4A7C15)>>40) & (s.hsize - 1)
	for {
		kk := s.hkey[h]
		if kk == k {
			return s.hval[h]
		}
		if kk == hEmpty {
			return -1
		}
		h = (h + 1) & (s.hsize - 1)
	}
}

func (s *Sim) hput(i, j, slot int) {
	k := pkey(i, j)
	h := int((uint64(k)*0x9E3779B97F4A7C15)>>40) & (s.hsize - 1)
	for {
		kk := s.hkey[h]
		if kk == hEmpty || kk == hTomb || kk == k {
			s.hkey[h] = k
			s.hval[h] = slot
			return
		}
		h = (h + 1) & (s.hsize - 1)
	}
}

func (s *Sim) hdel(i, j int) {
	k := pkey(i, j)
	h := int((uint64(k)*0x9E3779B97F4A7C15)>>40) & (s.hsize - 1)
	for {
		kk := s.hkey[h]
		if kk == k {
			s.hkey[h] = hTomb
			return
		}
		if kk == hEmpty {
			return
		}
		h = (h + 1) & (s.hsize - 1)
	}
}

func (s *Sim) slotNew(i, j int) int {
	if s.nfree == 0 {
		s.fatal(2, "FATAL: slot table full at pair (%d,%d) NSLOT=%d NC=%d "+
			"pos_i=(%.2f,%.2f,%.2f) pos_j=(%.2f,%.2f,%.2f)\n",
			i, j, s.NSLOT, s.NC, s.px[i], s.py[i], s.pz[i], s.px[j], s.py[j], s.pz[j])
	}
	s.nfree--
	sl := s.freelist[s.nfree]
	s.sli[sl] = i
	s.slj[sl] = j
	s.sst[sl] = sAlive
	s.slem[2*sl] = 0
	s.slem[2*sl+1] = 0
	s.slph[2*sl] = 0
	s.slph[2*sl+1] = 0
	s.slp[sl] = 1
	s.slq[sl] = 1
	s.sA[sl] = 0
	s.sldd[sl] = 0
	s.swl[sl] = 0
	s.sgg[sl] = 0 // cantus: a reborn bond starts mute
	s.rfp[sl] = 0 // registry: a reborn pair has no past
	s.rfm[sl] = 0
	s.rdel[2*sl] = 0
	s.rdel[2*sl+1] = 0
	s.swant[2*sl] = 0
	s.swant[2*sl+1] = 0
	s.sflux[sl] = 0
	s.hput(i, j, sl)
	if sl >= s.NSLOT {
		s.NSLOT = sl + 1
	}
	s.births++
	return sl
}

// candidate refresh: births, dying transitions, rule-α frees, live d/u.
// Canonical order everywhere (cells ascending, neighbours sorted).
func (s *Sim) topoRefresh() {
	s.binsBuild()
	var ncand int64
	g := s.binG
	var nbuf [512]int
	// mark: any ALIVE slot no longer in the candidate set -> DYING
	for sl := 0; sl < s.NSLOT; sl++ {
		if s.sst[sl] == sAlive {
			s.sst[sl] = sDying // provisional
		}
	}
	for i := 0; i < s.NC; i++ {
		bx, by, bz := int(s.px[i]/s.binSz), int(s.py[i]/s.binSz), int(s.pz[i]/s.binSz)
		if bx >= g {
			bx = g - 1
		}
		if by >= g {
			by = g - 1
		}
		if bz >= g {
			bz = g - 1
		}
		if bx < 0 {
			bx = 0
		}
		if by < 0 {
			by = 0
		}
		if bz < 0 {
			bz = 0
		}
		nn := 0
		for ax := bx - 1; ax <= bx+1; ax++ {
			for ay := by - 1; ay <= by+1; ay++ {
				for az := bz - 1; az <= bz+1; az++ {
					wx, wy2, wz := (ax+g)%g, (ay+g)%g, (az+g)%g
					for q := s.binHead[(wx*g+wy2)*g+wz]; q >= 0; q = s.binNext[q] {
						if q <= i {
							continue
						}
						dx := s.wr(s.px[q] - s.px[i])
						dy := s.wr(s.py[q] - s.py[i])
						dz := s.wr(s.pz[q] - s.pz[i])
						d2 := dx*dx + dy*dy + dz*dz
						cut := s.P.Cfac * (s.cr[i] + s.cr[q])
						if d2 >= cut*cut || d2 <= 1e-12 {
							continue
						}
						dup := false
						for a2 := 0; a2 < nn; a2++ {
							if nbuf[a2] == q {
								dup = true
							}
						}
						if !dup && nn < 512 {
							nbuf[nn] = q
							nn++
						}
					}
				}
			}
		}
		// canonical: sort the neighbour list ascending
		for a := 1; a < nn; a++ {
			v := nbuf[a]
			b := a - 1
			for b >= 0 && nbuf[b] > v {
				nbuf[b+1] = nbuf[b]
				b--
			}
			nbuf[b+1] = v
		}
		for a := 0; a < nn; a++ {
			j := nbuf[a]
			ncand++
			sl := s.hfind(i, j)
			if sl < 0 {
				s.slotNew(i, j)
			} else {
				s.sst[sl] = sAlive // re-affirm / resurrect
			}
		}
	}
	s.ncandLast = ncand
	// rule α: a DYING slot with no in-flight energy is freed
	for sl := 0; sl < s.NSLOT; sl++ {
		if s.sst[sl] == sDying && s.slem[2*sl] == 0 && s.slem[2*sl+1] == 0 {
			s.hdel(s.sli[sl], s.slj[sl])
			s.sst[sl] = sFree
			s.freelist[s.nfree] = sl
			s.nfree++
			s.deaths++
		}
	}
	// live geometry for every non-free slot
	for sl := 0; sl < s.NSLOT; sl++ {
		if s.sst[sl] == sFree {
			continue
		}
		i, j := s.sli[sl], s.slj[sl]
		dx := s.wr(s.px[j] - s.px[i])
		dy := s.wr(s.py[j] - s.py[i])
		dz := s.wr(s.pz[j] - s.pz[i])
		d := math.Sqrt(dx*dx + dy*dy + dz*dz)
		if d < 1e-9 {
			d = 1e-9
		}
		s.sd[sl] = d
		s.sux[sl] = dx / d
		s.suy[sl] = dy / d
		s.suz[sl] = dz / d
		ri, rj := s.cr[i], s.cr[j]
		A := 0.0
		if d < ri+rj {
			t := d*d - rj*rj + ri*ri
			a2 := (4.0*d*d*ri*ri - t*t) / (4.0 * d * d)
			if a2 > 0 {
				A = math.Pi * a2
			} else if d < math.Abs(ri-rj) {
				rm := ri
				if rj < ri {
					rm = rj
				}
				A = math.Pi * rm * rm
			}
		}
		s.sA[sl] = A
	}
	// CSR, canonical
	for i := 0; i <= s.NC; i++ {
		s.cls[i] = 0
	}
	for sl := 0; sl < s.NSLOT; sl++ {
		if s.sst[sl] == sFree {
			continue
		}
		s.cls[s.sli[sl]+1]++
		s.cls[s.slj[sl]+1]++
	}
	for i := 0; i < s.NC; i++ {
		s.cls[i+1] += s.cls[i]
	}
	fill := make([]int, s.NC)
	for sl := 0; sl < s.NSLOT; sl++ {
		if s.sst[sl] == sFree {
			continue
		}
		s.clidx[s.cls[s.sli[sl]]+fill[s.sli[sl]]] = sl
		fill[s.sli[sl]]++
		s.clidx[s.cls[s.slj[sl]]+fill[s.slj[sl]]] = sl
		fill[s.slj[sl]]++
	}
	// sort each cell's incident list by (other endpoint) for canonical order
	for i := 0; i < s.NC; i++ {
		lo, hi := s.cls[i], s.cls[i+1]
		for a := lo + 1; a < hi; a++ {
			v := s.clidx[a]
			vk := s.slj[v]
			if s.sli[v] != i {
				vk = s.sli[v]
			}
			b := a - 1
			for b >= lo {
				u := s.clidx[b]
				uk := s.slj[u]
				if s.sli[u] != i {
					uk = s.sli[u]
				}
				if uk <= vk {
					break
				}
				s.clidx[b+1] = u
				b--
			}
			s.clidx[b+1] = v
		}
	}
}

func (s *Sim) jamSettle() {
	// pure repulsion to jamming (LIVEFAB relax), geometry only
	var sw int
	mf, mfPrev, jk := 0.0, 1e300, s.P.JamK
	for sw = 0; sw < s.P.JamSweeps; sw++ {
		s.binsBuild()
		for i := 0; i < s.NC; i++ {
			s.fxb[i], s.fyb[i], s.fzb[i] = 0, 0, 0
		}
		g := s.binG
		for i := 0; i < s.NC; i++ {
			bx, by, bz := int(s.px[i]/s.binSz), int(s.py[i]/s.binSz), int(s.pz[i]/s.binSz)
			if bx >= g {
				bx = g - 1
			}
			if by >= g {
				by = g - 1
			}
			if bz >= g {
				bz = g - 1
			}
			if bx < 0 {
				bx = 0
			}
			if by < 0 {
				by = 0
			}
			if bz < 0 {
				bz = 0
			}
			for ax := bx - 1; ax <= bx+1; ax++ {
				for ay := by - 1; ay <= by+1; ay++ {
					for az := bz - 1; az <= bz+1; az++ {
						wx, wy2, wz := (ax+g)%g, (ay+g)%g, (az+g)%g
						for q := s.binHead[(wx*g+wy2)*g+wz]; q >= 0; q = s.binNext[q] {
							if q <= i {
								continue
							}
							dx := s.wr(s.px[q] - s.px[i])
							dy := s.wr(s.py[q] - s.py[i])
							dz := s.wr(s.pz[q] - s.pz[i])
							d2 := dx*dx + dy*dy + dz*dz
							contact := s.cr[i] + s.cr[q]
							if d2 >= contact*contact || d2 <= 0 {
								continue
							}
							d := math.Sqrt(d2)
							push := jk * (contact - d)
							if push > 0.3 {
								push = 0.3 // bounded shove
							}
							f := 0.0
							if d > 1e-6 {
								f = push / d
							}
							s.fxb[i] -= f * dx
							s.fyb[i] -= f * dy
							s.fzb[i] -= f * dz
							s.fxb[q] += f * dx
							s.fyb[q] += f * dy
							s.fzb[q] += f * dz
						}
					}
				}
			}
		}
		mf = 0
		for i := 0; i < s.NC; i++ {
			s.px[i] = s.fold(s.px[i] + s.fxb[i])
			s.py[i] = s.fold(s.py[i] + s.fyb[i])
			s.pz[i] = s.fold(s.pz[i] + s.fzb[i])
			f := math.Sqrt(s.fxb[i]*s.fxb[i] + s.fyb[i]*s.fyb[i] + s.fzb[i]*s.fzb[i])
			if f > mf {
				mf = f
			}
		}
		// adaptive damping: a rising witness means overshoot — soften
		if mf > mfPrev && jk > 1e-3 {
			jk *= 0.7
		}
		mfPrev = mf
		if mf < s.P.JamTol {
			sw++
			break
		}
	}
	fprintf(s.Out, "# JAM sweeps=%d max_step=%.3e (tol %.1e, k_end=%.4f) — witness for the settle\n",
		sw, mf, s.P.JamTol, jk)
}
