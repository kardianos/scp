package fab

// the step — passes ported from v89/freecell.c step() (itself ported
// from cellfab.c step_field()) in order: 0, S, 1, G0, 2, 3, 4, 5, D, F,
// 6, 7. Serial semantics preserved exactly.

import "math"

func (s *Sim) step() {
	dt := s.P.Dt
	P := &s.P

	// pass 0: pitch share of in-flight dense energy (cellfab.c:2795)
	for i := 0; i < s.NC; i++ {
		fl := 0.0
		for q := s.cls[i]; q < s.cls[i+1]; q++ {
			sl := s.clidx[q]
			fl += s.slem[2*sl] + s.slem[2*sl+1]
		}
		s.flload[i] = 0.5 * fl
	}

	// pass S: space-mode transport (cellfab.c:2811) — stale-by-one lA
	if P.SK > 0 {
		for sl := 0; sl < s.NSLOT; sl++ {
			s.swl[sl] = 0.0
			if s.sst[sl] == sFree || s.sA[sl] <= 0 {
				continue
			}
			i, j := s.sli[sl], s.slj[sl]
			pi := s.Es[i] + P.SDisp*(s.Em[i]+s.Ee[i])
			pj := s.Es[j] + P.SDisp*(s.Em[j]+s.Ee[j])
			dp := pi - pj
			if dp == 0 {
				continue
			}
			w := (s.sA[sl] / s.Aref) * (s.dref / s.sd[sl])
			s.swl[sl] = P.SK * dt * w * dp
			if P.BedK > 0 {
				s.swl[sl] *= s.sbed[sl] // FLOW: the bed
			}
		}
		for i := 0; i < s.NC; i++ {
			rq := 0.0
			for q := s.cls[i]; q < s.cls[i+1]; q++ {
				sl := s.clidx[q]
				f := s.swl[sl]
				if (f > 0 && s.sli[sl] == i) || (f < 0 && s.slj[sl] == i) {
					rq += math.Abs(f)
				}
			}
			s.sprq[i] = rq
			avail := s.Es[i] - P.EsFloor
			if avail <= 0 {
				s.sscl[i] = 0.0
			} else if rq > avail {
				s.sscl[i] = avail / rq
			} else {
				s.sscl[i] = 1.0
			}
		}
		if P.P1Meter != 0 {
			// P1 site 1: space transport src->dst = sgn(f) * d * u_hat
			for sl := 0; sl < s.NSLOT; sl++ {
				f := s.swl[sl]
				if f == 0 {
					continue
				}
				src := s.sli[sl]
				if f <= 0 {
					src = s.slj[sl]
				}
				mag := math.Abs(f) * s.sscl[src]
				if mag <= 0 {
					continue
				}
				m := mag * s.sd[sl]
				if f <= 0 {
					m = -m
				}
				s.p1sp[0] += m * s.sux[sl]
				s.p1sp[1] += m * s.suy[sl]
				s.p1sp[2] += m * s.suz[sl]
			}
		}
		for i := 0; i < s.NC; i++ {
			de := 0.0
			for q := s.cls[i]; q < s.cls[i+1]; q++ {
				sl := s.clidx[q]
				f := s.swl[sl]
				if f == 0 {
					continue
				}
				src := s.sli[sl]
				if f <= 0 {
					src = s.slj[sl]
				}
				mag := math.Abs(f) * s.sscl[src]
				if src == i {
					de -= mag
				} else {
					de += mag
				}
			}
			s.Es[i] += de
		}
		// FLOW (FLOW.md): the bed-digging channel law. Net memory from
		// the ACTUAL moved flow; growth on |net|; ZERO-SUM per-cell
		// renorm (the anti-ignition structure); clamp.
		if P.BedK > 0 {
			kb := 1.0
			if P.BedTau > dt {
				kb = dt / P.BedTau
			}
			for sl := 0; sl < s.NSLOT; sl++ {
				if s.sst[sl] == sFree {
					continue
				}
				f := s.swl[sl]
				flow := 0.0
				if f != 0 {
					src := s.sli[sl]
					if f <= 0 {
						src = s.slj[sl]
					}
					flow = math.Abs(f) * s.sscl[src]
					if f <= 0 {
						flow = -flow
					}
				}
				s.bednet[sl] += kb * (flow/dt - s.bednet[sl])
				s.sbed[sl] *= 1.0 + P.BedK*dt*math.Abs(s.bednet[sl])
			}
			for i := 0; i < s.NC; i++ {
				sum := 0.0
				n := 0
				for q := s.cls[i]; q < s.cls[i+1]; q++ {
					sl := s.clidx[q]
					if s.sst[sl] == sFree || s.sA[sl] <= 0 {
						continue
					}
					sum += s.sbed[sl]
					n++
				}
				if n > 0 && sum > 0 {
					s.bedf[i] = float64(n) / sum
				} else {
					s.bedf[i] = 1.0
				}
			}
			for sl := 0; sl < s.NSLOT; sl++ {
				if s.sst[sl] == sFree || s.sA[sl] <= 0 {
					continue
				}
				s.sbed[sl] *= math.Sqrt(s.bedf[s.sli[sl]] * s.bedf[s.slj[sl]])
				if s.sbed[sl] < 0.2 {
					s.sbed[sl] = 0.2
				}
				if s.sbed[sl] > 5.0 {
					s.sbed[sl] = 5.0
				}
			}
		}
	}

	// pass 1: live radii, effective frequencies (cellfab.c:2890)
	for i := 0; i < s.NC; i++ {
		ratio := s.Es[i] / P.Es0
		if ratio < 0 {
			ratio = 0
		}
		s.cr[i] = s.cr0[i] * math.Cbrt(ratio)
		x := (s.Em[i] + s.flload[i]) / P.Cap
		det := 1.0 + P.QDetune*x
		s.w1e[i] = P.W1 / det
		s.w2e[i] = P.W2 / det
	}

	// G0: topology from live geometry
	if P.FreezeGeo == 0 {
		s.topoRefresh()
	}

	// pass 2: channel lens area + dense wants + bond misfit buffer
	if P.AmpDrv > 0 && P.AmpTau > 0 {
		for i := 0; i < s.NC; i++ { // L-1 zero-sum bookkeeping
			s.ampre[i] = 0
		}
	}
	for sl := 0; sl < s.NSLOT; sl++ {
		if s.sst[sl] == sFree {
			s.sA[sl] = 0
			continue
		}
		i, j := s.sli[sl], s.slj[sl]
		d, ri, rj := s.sd[sl], s.cr[i], s.cr[j]
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
		s.swant[2*sl] = 0
		s.swant[2*sl+1] = 0
		s.shau[sl] = 0
		s.sflux[sl] = 0
		s.sldd[sl] = 0
		if A <= 0 {
			continue
		}

		geo := (A / s.Aref) * (s.dref / d)
		occI := (s.Em[i] + s.Ee[i]) / P.Cap
		occJ := (s.Em[j] + s.Ee[j]) / P.Cap
		headI := 1.0 - occI
		if headI < 0 {
			headI = 0
		}
		if headI > 1 {
			headI = 1
		}
		headJ := 1.0 - occJ
		if headJ < 0 {
			headJ = 0
		}
		if headJ > 1 {
			headJ = 1
		}

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
		axij := axi * axj
		if axij < 1e-10 {
			continue
		}

		mobI, mobJ := s.Em[i], s.Em[j]
		if mobI < 1e-15 && mobJ < 1e-15 {
			continue
		}
		dnn := s.n2x[i]*s.n2x[j] + s.n2y[i]*s.n2y[j] + s.n2z[i]*s.n2z[j]
		gpl := axij * dnn * dnn
		if gpl < 1e-8 {
			continue
		}

		wi, wj := s.w2e[i], s.w2e[j]
		thi, thj := s.th2[i], s.th2[j]

		// C1 sympathetic joining through coincident partials
		best := -1.0
		bp, bq := 1, 1
		for tt := 0; tt < s.ncomb; tt++ {
			pp := float64(s.combp[tt])
			qq := float64(s.combq[tt])
			gw := s.gm / (pp * qq)
			det := qq*wi - pp*wj
			sc := (gw * gw / (gw*gw + det*det)) / (pp * qq)
			if sc > best {
				best = sc
				bp = s.combp[tt]
				bq = s.combq[tt]
			}
		}
		res := best
		s.slp[sl] = int8(bp)
		s.slq[sl] = int8(bq)
		fbp, fbq := float64(bp), float64(bq)
		gIJ := s.gateOf(fbq*thi - fbq*wi*d/P.C - fbp*thj)
		gJI := s.gateOf(fbp*thj - fbp*wj*d/P.C - fbq*thi)
		miEff, mjEff := mobI, mobJ
		if P.MobSym != 0 {
			fl0 := P.MobFloor * P.Cap
			mjr := mobJ
			if mjr < fl0 {
				mjr = fl0
			}
			mir := mobI
			if mir < fl0 {
				mir = fl0
			}
			miEff = math.Sqrt(mobI * mjr)
			mjEff = math.Sqrt(mobJ * mir)
		}
		kd := P.KDep * P.KDepM
		base := kd * dt * geo * gpl * res
		wIJ := base * gIJ * headJ * miEff
		wJI := base * gJI * headI * mjEff

		if P.KappaReac > 0 {
			// S2 — the choir's correction, DERIVED (cellfab.c:3042 verbatim)
			psIJ := wrapPi(fbq*thi - fbq*wi*d/P.C - fbp*thj)
			psJI := wrapPi(fbp*thj - fbp*wj*d/P.C - fbq*thi)
			Sm2 := math.Sqrt(miEff * mjEff)
			hh := math.Sqrt(headI * headJ)
			reac := P.KappaReac * 0.5 * base * hh * Sm2
			wIJ -= reac * gIJ * lutSin(psIJ)
			wJI -= reac * gJI * lutSin(psJI)
			if wIJ < 0 {
				wIJ = 0
			}
			if wJI < 0 {
				wJI = 0
			}
		}

		// THE BOND — P15 (cellfab.c:3068) on the geometric conjugate,
		// buffered Jacobi, applied to POSITIONS in pass D
		if P.KappaBond > 0 {
			psF := wrapPi(fbq*thi - fbq*wi*d/P.C - fbp*thj)
			psB := wrapPi(fbp*thj - fbp*wj*d/P.C - fbq*thi)
			Smp := math.Sqrt(miEff * mjEff)
			if Smp > 0 {
				hf := 0.5 * (1.0 + lutCos(psF))
				hb := 0.5 * (1.0 + lutCos(psB))
				Gf := math.Pow(hf, P.PGate)
				Gb := math.Pow(hb, P.PGate)
				Gfp := 0.0
				if hf > 1e-12 {
					Gfp = -0.5 * P.PGate * (Gf / hf) * lutSin(psF)
				}
				Gbp := 0.0
				if hb > 1e-12 {
					Gbp = -0.5 * P.PGate * (Gb / hb) * lutSin(psB)
				}
				dVdd := (Gfp*Gb*fbq*wi + Gf*Gbp*fbp*wj) / P.C
				dd := -P.KappaBond * base * Smp * dVdd
				// integrator guard: cap the per-step kick
				if dd > 0.05 {
					dd = 0.05
				}
				if dd < -0.05 {
					dd = -0.05
				}
				s.sldd[sl] = dd
			}
		}

		// v92 AMPLITUDE Phase L lane L-1 (L0_DESIGN.md §3, rev phase-current
		// per L1_FINDINGS.md finding 3): shadow promoted from meter to
		// driver via its IN-PHASE current — project each direction's shadow
		// onto the receiver's transport frame (closure phase), so the bias
		// responds to a phase gradient (the e3b tilt). Signed, saturating,
		// slot-borne. amp_mmin: 1 = unison (e3b/p1); 2 = chords. (Mirror of
		// freecell.c pass-2 L-1 block.)
		if P.AmpDrv > 0 && P.AmpTau > 0 {
			ms := int(s.slp[sl])
			if int(s.slq[sl]) > ms {
				ms = int(s.slq[sl])
			}
			if float64(ms) >= P.AmpMmin {
				// record pre-bias outflow for the zero-sum renorm
				s.ampre[i] += wIJ
				s.ampre[j] += wJI
				phr0 := float64(s.slp[sl]) * s.th2[j]
				J0 := s.sreA[2*sl]*lutCos(phr0) + s.simA[2*sl]*lutSin(phr0)
				phr1 := float64(s.slq[sl]) * s.th2[i]
				J1 := s.sreA[2*sl+1]*lutCos(phr1) + s.simA[2*sl+1]*lutSin(phr1)
				Jnet := J0 - J1
				Jabs := math.Abs(J0) + math.Abs(J1) + 1e-15
				bias := P.AmpDrv * Jnet / Jabs
				wIJ *= (1.0 + bias)
				wJI *= (1.0 - bias)
			}
		}

		if P.AmpNat > 0 {
			// v93 UNITARY DENSE CHANNEL: fold the want ingredients into a
			// pairwise hop angle tau_s (mobility carried by the amplitude
			// magnitude itself, as in pass F; door enforces cap; closure
			// gate survives as the angle envelope). Byte-inert at amp_nat=0.
			gsym := math.Sqrt(gIJ * gJI)
			tau := P.AmpNat * base * gsym * math.Sqrt(headI*headJ)
			if tau > 0.5 {
				tau = 0.5
			}
			s.shau[sl] = tau
		} else {
			if wIJ > 0 {
				s.swant[2*sl] = wIJ
			}
			if wJI > 0 {
				s.swant[2*sl+1] = wJI
			}
		}
	}

	// L-1 zero-sum renorm (FLOW architecture, L1_FINDINGS/grok fix): hold
	// each cell's total outflow at its pre-bias level — the bias only
	// REDISTRIBUTES direction (anti-ignition + coherence-runaway bound).
	if P.AmpDrv > 0 && P.AmpTau > 0 {
		for i := 0; i < s.NC; i++ {
			if s.ampre[i] <= 0 {
				continue
			}
			post := 0.0
			for q := s.cls[i]; q < s.cls[i+1]; q++ {
				sl := s.clidx[q]
				dir := 0
				if s.sli[sl] != i {
					dir = 1
				}
				post += s.swant[2*sl+dir]
			}
			if post > 1e-15 {
				fac := s.ampre[i] / post
				for q := s.cls[i]; q < s.cls[i+1]; q++ {
					sl := s.clidx[q]
					dir := 0
					if s.sli[sl] != i {
						dir = 1
					}
					s.swant[2*sl+dir] *= fac
				}
			}
		}
	}

	// pass U: v93 UNITARY DENSE CHANNEL — within-mode dense transport as a
	// product of UNITARY PAIRWISE PLANE ROTATIONS (pass F's cousin,
	// v93/README.md §II.3). Engages only when amp_nat>0; in that regime the
	// additive magnitude want above is bypassed (swant stays 0, passes 3-5
	// are no-ops). The dense amplitude psi_m = sqrt(Em) e^{i th2} hops
	// between link endpoints by a Givens rotation of angle shau[sl]; each
	// hop conserves the two-cell norm EXACTLY (conservation is a theorem of
	// the update). The cross term 2 Im(psi_i* psi_j) IS the link current =
	// dense momentum. The door (pass 6) is never unitarized. Hops applied in
	// canonical link order (from the i side). Byte-inert at amp_nat=0.
	if P.AmpNat > 0 {
		for i := 0; i < s.NC; i++ {
			e := s.Em[i]
			if e > 0 {
				r := math.Sqrt(e)
				s.dm1[i] = r * lutCos(s.th2[i])
				s.dm2[i] = r * lutSin(s.th2[i])
			} else {
				s.dm1[i] = 0
				s.dm2[i] = 0
			}
		}
		for i := 0; i < s.NC; i++ {
			for q := s.cls[i]; q < s.cls[i+1]; q++ {
				sl := s.clidx[q]
				if s.sli[sl] != i {
					continue // canonical: apply from the i side
				}
				if s.sst[sl] == sFree || s.sA[sl] <= 0 {
					continue
				}
				tau := s.shau[sl]
				if tau <= 0 {
					continue
				}
				j := s.slj[sl]
				ss, cc := math.Sincos(tau)
				m1i, m2i, m1j, m2j := s.dm1[i], s.dm2[i], s.dm1[j], s.dm2[j]
				s.dm1[i] = cc*m1i + ss*m2j
				s.dm2[i] = cc*m2i - ss*m1j
				s.dm1[j] = cc*m1j + ss*m2i
				s.dm2[j] = cc*m2j - ss*m1i
				if P.P1Meter != 0 {
					// P1 site (dense): energy the rotation moved into j
					mj := (s.dm1[j]*s.dm1[j] + s.dm2[j]*s.dm2[j]) - (m1j*m1j + m2j*m2j)
					m := mj * s.sd[sl]
					s.p1fd[0] += m * s.sux[sl]
					s.p1fd[1] += m * s.suy[sl]
					s.p1fd[2] += m * s.suz[sl]
				}
			}
		}
		for i := 0; i < s.NC; i++ {
			a, b := s.dm1[i], s.dm2[i]
			s.Em[i] = a*a + b*b
			ph := math.Atan2(b, a)
			if ph < 0 {
				ph += TwoPi
			}
			s.th2[i] = ph
		}
	}

	// pass 3: outflow limiter (cellfab.c:3109), dense only
	for i := 0; i < s.NC; i++ {
		r1 := 0.0
		for q := s.cls[i]; q < s.cls[i+1]; q++ {
			sl := s.clidx[q]
			dir := 1
			if s.sli[sl] == i {
				dir = 0
			}
			r1 += s.swant[2*sl+dir]
		}
		s.req1[i] = r1
		a1 := 0.98 * s.Em[i]
		if r1 > a1 && r1 > 0 {
			s.scl1[i] = a1 / r1
		} else {
			s.scl1[i] = 1.0
		}
	}

	// pass 4: resolve deposits, debit sources, entrain receivers
	for i := 0; i < s.NC; i++ {
		s.th2s[i] = s.th2[i]
	}
	for sl := 0; sl < s.NSLOT; sl++ {
		if s.sst[sl] == sFree || s.sA[sl] <= 0 {
			continue
		}
		for dir := 0; dir < 2; dir++ {
			f := s.swant[2*sl+dir]
			if f <= 0 {
				continue
			}
			src := s.sli[sl]
			if dir != 0 {
				src = s.slj[sl]
			}
			f *= s.scl1[src]
			s.swant[2*sl+dir] = f
			if s.slem[2*sl+dir] <= 0 {
				s.slph[2*sl+dir] = 0
			}
			s.slem[2*sl+dir] += f
			// identity lane (IDENTITY.md §1.3.2): the departing parcel
			// carries the depositor's episode gid (I-D1 last-depositor)
			if P.ParTau > 0 {
				s.slgid[2*sl+dir] = s.cgid[src]
			}
			// AMPLITUDE Phase M: shadow deposit sqrt(f) e^{i m th}
			if P.AmpTau > 0 {
				m := float64(s.slq[sl])
				if dir != 0 {
					m = float64(s.slp[sl])
				}
				ph := m * s.th2[src]
				s.sreA[2*sl+dir] += math.Sqrt(f) * lutCos(ph)
				s.simA[2*sl+dir] += math.Sqrt(f) * lutSin(ph)
			}
			s.sflux[sl] += f
			s.sfluxd[2*sl+dir] += f // directed ledger for circulation
			if P.P1Meter != 0 && f > 0 {
				// P1 site 2: dense deposit src -> link midpoint
				m := 0.5 * f * s.sd[sl]
				if dir != 0 {
					m = -m
				}
				s.p1fl[0] += m * s.sux[sl]
				s.p1fl[1] += m * s.suy[sl]
				s.p1fl[2] += m * s.suz[sl]
			}
		}
	}
	for i := 0; i < s.NC; i++ {
		d1 := 0.0
		for q := s.cls[i]; q < s.cls[i+1]; q++ {
			sl := s.clidx[q]
			if s.sA[sl] <= 0 {
				continue
			}
			dir := 1
			if s.sli[sl] == i {
				dir = 0
			}
			f := s.swant[2*sl+dir]
			if f > 0 {
				d1 += f
			}
		}
		s.Em[i] -= d1
	}
	for i := 0; i < s.NC; i++ {
		for q := s.cls[i]; q < s.cls[i+1]; q++ {
			sl := s.clidx[q]
			if s.sA[sl] <= 0 {
				continue
			}
			dir := 0 // the direction ARRIVING at i
			if s.sli[sl] == i {
				dir = 1
			}
			src := s.sli[sl]
			if dir != 0 {
				src = s.slj[sl]
			}
			f := s.swant[2*sl+dir]
			if f <= 0 {
				continue
			}
			mobr := s.Em[i]
			wsrc := s.w2e[src]
			thsrc := s.th2s[src]
			var ms, mr float64
			if dir == 0 {
				ms, mr = float64(s.slq[sl]), float64(s.slp[sl])
			} else {
				ms, mr = float64(s.slp[sl]), float64(s.slq[sl])
			}
			err := wrapPi(ms*(thsrc-wsrc*s.sd[sl]/P.C)-mr*s.th2[i]) / mr
			mix := f / (f + mobr + s.lockf)
			s.th2[i] += P.KappaLock * mix * err
		}
	}

	// pass 5: advance cycles; deliver on completion (cellfab.c:3198).
	// ALIVE pinched channels hold flight; DYING channels advance and
	// flush (rule α), returning unreceivable residue to their source.
	for i2 := 0; i2 < s.NC; i2++ {
		s.th2s[i2] = s.th2[i2]
	}
	for recv := 0; recv < s.NC; recv++ {
		for q := s.cls[recv]; q < s.cls[recv+1]; q++ {
			sl := s.clidx[q]
			dying := s.sst[sl] == sDying
			if s.sA[sl] <= 0 && !dying {
				continue
			}
			i, j := s.sli[sl], s.slj[sl]
			dir := 0
			if s.sli[sl] == recv {
				dir = 1
			}
			send := i
			if dir != 0 {
				send = j
			}
			adv := dt * P.C / s.sd[sl]
			sslot := 2*sl + dir
			if s.slem[sslot] <= 0 {
				continue
			}
			s.slph[sslot] += adv
			// AMPLITUDE Phase M: transit rotation -m w dt
			if P.AmpTau > 0 && (s.sreA[sslot] != 0 || s.simA[sslot] != 0) {
				m := float64(s.slq[sl])
				if dir != 0 {
					m = float64(s.slp[sl])
				}
				dph := -m * s.w2e[send] * dt
				c := lutCos(dph)
				sn := lutSin(dph)
				re0 := s.sreA[sslot]
				s.sreA[sslot] = re0*c - s.simA[sslot]*sn
				s.simA[sslot] = re0*sn + s.simA[sslot]*c
			}
			if s.slph[sslot] < 1.0 {
				continue
			}
			freec := P.Cap - (s.Em[recv] + s.Ee[recv])
			take := s.slem[sslot]
			if take > freec {
				if freec > 0 {
					take = freec
				} else {
					take = 0
				}
			}
			if take > 0 {
				mobprev := s.Em[recv]
				// v91 registry stamp (REGISTRY.md §1.3): a delivered
				// parcel between two continuous identities. Deliveries
				// only — the rule-alpha flush is a RETURN, not stamped.
				if P.RegTau > 0 {
					s.rdel[sslot] += take
				}
				// identity lane (IDENTITY.md §1.3.3): identity-carried
				// delivery — counts only when the parcel's label is the
				// source endpoint's LIVING episode. Deliveries only.
				if P.ParTau > 0 {
					s.parDelTot += take
					if s.slgid[sslot] != 0 && s.slgid[sslot] == s.cgid[send] {
						s.pdel[sslot] += take
						s.parDelID += take
					}
				}
				// AMPLITUDE Phase M: deliver the amplitude fraction,
				// composed into the dst chart frame
				if P.AmpTau > 0 {
					efl := s.slem[sslot]
					a2 := s.sreA[sslot]*s.sreA[sslot] + s.simA[sslot]*s.simA[sslot]
					if efl > 1e-30 && a2 > 0 {
						frac := take / efl
						if frac > 1 {
							frac = 1
						}
						sf := math.Sqrt(frac)
						dre := s.sreA[sslot] * sf
						dim := s.simA[sslot] * sf
						keep := math.Sqrt(1.0 - frac)
						s.sreA[sslot] *= keep
						s.simA[sslot] *= keep
						n := float64(s.slp[sl])
						if dir != 0 {
							n = float64(s.slq[sl])
						}
						ph := -n * s.th2[recv]
						c := lutCos(ph)
						sn := lutSin(ph)
						s.amdre[sl] += dre*c - dim*sn
						s.amdim[sl] += dre*sn + dim*c
						s.amdc[sl] += math.Sqrt(dre*dre + dim*dim)
					}
				}
				s.slem[sslot] -= take
				if P.P1Meter != 0 {
					// P1 site 3: flight arrival, link midpoint -> recv
					m := 0.5 * take * s.sd[sl]
					if dir != 0 {
						m = -m
					}
					s.p1fl[0] += m * s.sux[sl]
					s.p1fl[1] += m * s.suy[sl]
					s.p1fl[2] += m * s.suz[sl]
				}
				// C2: dissonance radiates (cellfab.c:3228)
				det := float64(s.slq[sl])*s.w2e[i] - float64(s.slp[sl])*s.w2e[j]
				R := 2.0 * math.Abs(det) * P.GammaRough /
					(det*det + P.GammaRough*P.GammaRough)
				rough := take * P.RoughK * R
				reF := s.A0eff * s.w2e[recv] / TwoPi
				reD := s.A0eff * s.w1e[recv] / TwoPi
				rough = s.atomsFire(rough, reF, reD, &s.qcnvF[recv])
				rough = s.atomsClamp(rough, take, reF, reD, &s.qcnvF[recv])
				if rough > 0 {
					s.roughq[recv] += rough
				}
				backS := rough * P.SPull / (1.0 + P.SPull)
				s.backsTotal += backS
				s.Em[recv] += take - rough
				s.fieldInject(recv, rough-backS)
				s.Es[recv] += backS
				s.roughTotal += rough
				if s.tag[recv] != 0 {
					s.ctRough += rough
					s.ctBacks += backS
				}
				wsend := s.w2e[send]
				thsend := s.th2s[send]
				var ms, mr float64
				if dir == 0 {
					ms, mr = float64(s.slq[sl]), float64(s.slp[sl])
				} else {
					ms, mr = float64(s.slp[sl]), float64(s.slq[sl])
				}
				err := wrapPi(ms*(thsend-wsend*s.sd[sl]/P.C)-mr*s.th2[recv]) / mr
				mix := take / (take + mobprev + s.lockf)
				s.th2[recv] += P.KappaLock * mix * err
			}
			if dying && s.slem[sslot] > 1e-17 {
				// rule α flush: unreceivable residue returns to source
				if P.P1Meter != 0 {
					// P1 site 4: flush, link midpoint -> send (return)
					m := -0.5 * s.slem[sslot] * s.sd[sl]
					if dir != 0 {
						m = -m
					}
					s.p1fl[0] += m * s.sux[sl]
					s.p1fl[1] += m * s.suy[sl]
					s.p1fl[2] += m * s.suz[sl]
				}
				s.Em[send] += s.slem[sslot]
				s.betaEnergy += s.slem[sslot]
				s.betaReturns++
				s.slem[sslot] = 0
				s.slph[sslot] = 0
				s.sreA[sslot] = 0 // flush = return, not delivery
				s.simA[sslot] = 0
			}
			if s.slem[sslot] <= 1e-17 {
				s.slem[sslot] = 0
				s.slph[sslot] = 0
				s.sreA[sslot] = 0
				s.simA[sslot] = 0
			} else if take <= 0 {
				s.slph[sslot] = 0
			} else {
				s.slph[sslot] -= 1.0
			}
		}
	}

	// pass H0: v91 EXCHANGE REGISTRY ledger (REGISTRY.md §1.3 item 2)
	// — mirror of the C kernel's pass H0, operation for operation.
	// All non-free slots (incl. pinched), sl ascending; reg_tau=0 =>
	// does not execute.
	if P.RegTau > 0 {
		kreg := 1.0
		if P.RegTau > dt {
			kreg = dt / P.RegTau
		}
		for sl := 0; sl < s.NSLOT; sl++ {
			if s.sst[sl] == sFree {
				continue
			}
			s.rfp[sl] += kreg * (s.rdel[2*sl]/dt - s.rfp[sl])
			s.rfm[sl] += kreg * (s.rdel[2*sl+1]/dt - s.rfm[sl])
			s.rdel[2*sl] = 0
			s.rdel[2*sl+1] = 0
		}
	}

	// pass H0b: v91 IDENTITY lane (IDENTITY.md §1.3 items 1+4) —
	// mirror of the C kernel's pass H0b, operation for operation:
	// episode gids by hysteresis on x=(Em+flload)/cap, then the
	// identity-carried ledgers and bond stamps. par_tau=0 => skip.
	if P.ParTau > 0 {
		for i := 0; i < s.NC; i++ {
			xep := (s.Em[i] + s.flload[i]) / P.Cap
			if s.cgid[i] == 0 {
				if xep >= P.ParHi {
					s.cgid[i] = s.gidNext
					s.gidNext++
					s.cbirth[i] = s.simT
					s.parMints++
				}
			} else if xep < P.ParLo {
				s.cgid[i] = 0
				s.parRetires++
			}
		}
		kpar := 1.0
		if P.ParTau > dt {
			kpar = dt / P.ParTau
		}
		for sl := 0; sl < s.NSLOT; sl++ {
			if s.sst[sl] == sFree {
				continue
			}
			s.pdp[sl] += kpar * (s.pdel[2*sl]/dt - s.pdp[sl])
			s.pdm[sl] += kpar * (s.pdel[2*sl+1]/dt - s.pdm[sl])
			s.pdel[2*sl] = 0
			s.pdel[2*sl+1] = 0
			ga := s.cgid[s.sli[sl]]
			gb := s.cgid[s.slj[sl]]
			if s.pstampa[sl] == 0 {
				if ga != 0 && gb != 0 && s.pdp[sl] > 0 && s.pdm[sl] > 0 {
					s.pstampa[sl] = ga
					s.pstampb[sl] = gb
					s.pstampt[sl] = s.simT
				}
			} else if ga != s.pstampa[sl] || gb != s.pstampb[sl] {
				s.pstampa[sl] = 0
				s.pstampb[sl] = 0
				s.pstampt[sl] = 0
			}
		}
	}

	// pass H0c: v91 AMPLITUDE Phase M (AMPLITUDE.md §1) — low-pass the
	// per-step delivered shadow. amp_tau=0 => skip.
	if P.AmpTau > 0 {
		kam := 1.0
		if P.AmpTau > dt {
			kam = dt / P.AmpTau
		}
		for sl := 0; sl < s.NSLOT; sl++ {
			if s.sst[sl] == sFree {
				continue
			}
			s.amre[sl] += kam * (s.amdre[sl]/dt - s.amre[sl])
			s.amim[sl] += kam * (s.amdim[sl]/dt - s.amim[sl])
			s.amc[sl] += kam * (s.amdc[sl]/dt - s.amc[sl])
			s.amdre[sl] = 0
			s.amdim[sl] = 0
			s.amdc[sl] = 0
		}
	}

	// pass H: v91 CANTUS (coherent-channel candidate B, CANTUS.md) —
	// mirror of the C kernel's pass H, operation for operation. See
	// kernel/freecell.c for the full rationale comment.
	if P.KCant > 0 || P.KTune > 0 {
		ktau := 1.0
		if P.CantTau > dt {
			ktau = dt / P.CantTau
		}
		for i := 0; i < s.NC; i++ {
			s.th2s[i] = s.th2[i]
			s.dthH[i] = 0
		}
		for i := 0; i < s.NC; i++ {
			for q2 := s.cls[i]; q2 < s.cls[i+1]; q2++ {
				sl := s.clidx[q2]
				if s.sli[sl] != i {
					continue // visit once, from i
				}
				if s.sst[sl] == sFree || s.sA[sl] <= 0 {
					continue
				}
				j := s.slj[sl]
				if s.Em[i] <= 1e-15 || s.Em[j] <= 1e-15 {
					continue
				}
				pp := float64(s.slp[sl])
				qq := float64(s.slq[sl])
				d := s.sd[sl]
				psF := wrapPi(qq*s.th2s[i] - qq*s.w2e[i]*d/P.C - pp*s.th2s[j])
				psB := wrapPi(pp*s.th2s[j] - pp*s.w2e[j]*d/P.C - qq*s.th2s[i])
				gg := s.gateOf(psF) * s.gateOf(psB)
				// v1.1: the LINK's own gauge memory (holds through
				// lens blinks — non-eligible steps skip this update).
				// reg_gate=1 (REGISTRY.md R-G3): growth target gated
				// by the registry match (form F-A) — coherence may
				// only grow on identity-continuous reciprocal exchange.
				if P.CantGrow != 0 || s.sgg[sl] > 0 {
					tgt := gg
					if P.RegGate != 0 {
						gross := s.rfp[sl] + s.rfm[sl]
						mn := s.rfp[sl]
						if s.rfm[sl] < mn {
							mn = s.rfm[sl]
						}
						mult := 0.0
						if P.RegGate == 1 {
							// F-B: balance x flow saturation (f0=0 => F-A)
							if gross > 0 {
								mult = (2.0 * mn / gross) * (gross / (gross + P.RegF0))
							}
						} else {
							// F-D: reciprocal-flow saturation, s = 2*min
							s2 := 2.0 * mn
							if s2 > 0 {
								mult = s2 / (s2 + P.RegF0)
							}
						}
						tgt *= mult
					}
					if P.ParGate != 0 {
						// IDENTITY.md §1.3.5-6: maturity-clocked
						// identity continuity — stamp AGE is the
						// gate variable the lock cannot manufacture.
						rid := 0.0
						if s.pstampa[sl] != 0 &&
							s.cgid[i] == s.pstampa[sl] &&
							s.cgid[j] == s.pstampb[sl] &&
							s.pdp[sl] > 0 && s.pdm[sl] > 0 {
							page := s.simT - s.pstampt[sl]
							mat := 1.0
							if P.ParMature > 0 && page < P.ParMature {
								mat = page / P.ParMature
							}
							if P.ParForm == 1 {
								gr2 := s.pdp[sl] + s.pdm[sl]
								mn2 := s.pdp[sl]
								if s.pdm[sl] < mn2 {
									mn2 = s.pdm[sl]
								}
								if gr2 > 0 {
									rid = mat * 2.0 * mn2 / gr2
								}
							} else {
								rid = mat
							}
						}
						tgt *= rid
					}
					s.sgg[sl] += ktau * (tgt - s.sgg[sl])
				}
				amp := s.sgg[sl]
				if amp <= 0 {
					continue
				}
				if P.KCant > 0 {
					e := 0.5 * wrapPi(psF-psB)
					wl := P.KCant * dt * amp
					n2q := pp*pp + qq*qq
					s.dthH[i] -= wl * (qq / n2q) * e
					s.dthH[j] += wl * (pp / n2q) * e
				}
				if P.KTune > 0 {
					ui := s.Em[i]/P.Cap - s.cxl[i]
					uj := s.Em[j]/P.Cap - s.cxl[j]
					J := P.KTune * dt * amp * P.Cap * (ui - uj)
					src, dst := i, j
					if J <= 0 {
						src, dst = j, i
					}
					mag := math.Abs(J)
					av := 0.98 * s.Em[src]
					if mag > av {
						mag = av
					}
					freec := P.Cap - (s.Em[dst] + s.Ee[dst])
					if freec < 0 {
						freec = 0
					}
					if mag > freec {
						mag = freec
					}
					if mag > 0 {
						s.Em[src] -= mag
						s.Em[dst] += mag
						s.tuneTotal += mag
						if P.P1Meter != 0 {
							// P1 site H: within-mode current src->dst over
							// the full link (deposit+arrival telescoped)
							m := mag * d
							if src != i {
								m = -m
							}
							s.p1fl[0] += m * s.sux[sl]
							s.p1fl[1] += m * s.suy[sl]
							s.p1fl[2] += m * s.suz[sl]
						}
					}
				}
			}
		}
		for i := 0; i < s.NC; i++ {
			if s.dthH[i] != 0 {
				s.th2[i] += s.dthH[i]
			}
			s.cxl[i] += ktau * (s.Em[i]/P.Cap - s.cxl[i])
		}
	}

	// pass D: MOTION — the free-cell core. Bond displacements (Jacobi
	// buffered, ±u_hat/2) + contact repulsion, overdamped.
	if P.FreezeGeo == 0 {
		for i := 0; i < s.NC; i++ {
			s.fxb[i], s.fyb[i], s.fzb[i] = 0, 0, 0
		}
		for i := 0; i < s.NC; i++ {
			dxx, dyy, dzz := 0.0, 0.0, 0.0
			if s.pin[i] != 0 {
				s.fxb[i], s.fyb[i], s.fzb[i] = 0, 0, 0
				continue
			}
			for q := s.cls[i]; q < s.cls[i+1]; q++ {
				sl := s.clidx[q]
				sgn := 1.0
				if s.sli[sl] == i {
					sgn = -1.0
				}
				// bond: d changes by sldd -> each endpoint ±sldd/2
				if s.sldd[sl] != 0 {
					h := 0.5 * s.sldd[sl] * sgn
					dxx += h * s.sux[sl]
					dyy += h * s.suy[sl]
					dzz += h * s.suz[sl]
				}
				// repulsion below contact
				contact := s.cr[s.sli[sl]] + s.cr[s.slj[sl]]
				if s.sd[sl] < contact {
					f := P.KRep * (contact - s.sd[sl]) * sgn
					dxx += P.MobGeo * P.Dt * f * s.sux[sl]
					dyy += P.MobGeo * P.Dt * f * s.suy[sl]
					dzz += P.MobGeo * P.Dt * f * s.suz[sl]
				}
			}
			// per-step displacement cap (quasi-static motion guard)
			dm := math.Sqrt(dxx*dxx + dyy*dyy + dzz*dzz)
			if dm > 0.1 {
				sc := 0.1 / dm
				dxx *= sc
				dyy *= sc
				dzz *= sc
			}
			s.fxb[i], s.fyb[i], s.fzb[i] = dxx, dyy, dzz
		}
		if P.P1Meter != 0 {
			// P1 site 6: geometry — cells carry their energy when they
			// move; link flight (midpoint) moves with the endpoint mean
			for i := 0; i < s.NC; i++ {
				ei := s.Es[i] + s.Em[i] + s.Ee[i]
				s.p1gm[0] += ei * s.fxb[i]
				s.p1gm[1] += ei * s.fyb[i]
				s.p1gm[2] += ei * s.fzb[i]
			}
			for sl := 0; sl < s.NSLOT; sl++ {
				if s.sst[sl] == sFree {
					continue
				}
				le := s.slem[2*sl] + s.slem[2*sl+1]
				if le <= 0 {
					continue
				}
				i2, j2 := s.sli[sl], s.slj[sl]
				s.p1gm[0] += le * 0.5 * (s.fxb[i2] + s.fxb[j2])
				s.p1gm[1] += le * 0.5 * (s.fyb[i2] + s.fyb[j2])
				s.p1gm[2] += le * 0.5 * (s.fzb[i2] + s.fzb[j2])
			}
		}
		for i := 0; i < s.NC; i++ {
			s.px[i] = s.fold(s.px[i] + s.fxb[i])
			s.py[i] = s.fold(s.py[i] + s.fyb[i])
			s.pz[i] = s.fold(s.pz[i] + s.fzb[i])
			if math.IsNaN(s.px[i]) || math.IsNaN(s.py[i]) || math.IsNaN(s.pz[i]) {
				fprintf(s.Errw, "FATAL NaN position cell %d at t=%.3f: f=(%g,%g,%g) "+
					"Em=%g Es=%g cr=%g\n", i, s.simT, s.fxb[i], s.fyb[i], s.fzb[i],
					s.Em[i], s.Es[i], s.cr[i])
				for q := s.cls[i]; q < s.cls[i+1]; q++ {
					sl := s.clidx[q]
					fprintf(s.Errw, "  slot %d (%d-%d) st=%d d=%g A=%g ldd=%g lem=%g/%g\n",
						sl, s.sli[sl], s.slj[sl], s.sst[sl], s.sd[sl], s.sA[sl], s.sldd[sl],
						s.slem[2*sl], s.slem[2*sl+1])
				}
				s.fatal(4, "")
			}
		}
	}

	// pass F: unitary field hops (cellfab.c:3289), canonical link order
	for i := 0; i < s.NC; i++ {
		ang := s.w1e[i] * dt
		ss, cc := math.Sincos(ang)
		a1, a2 := s.fa1[i], s.fa2[i]
		s.fa1[i] = cc*a1 + ss*a2
		s.fa2[i] = -ss*a1 + cc*a2
	}
	for i := 0; i < s.NC; i++ {
		sacc := 0.0
		for q := s.cls[i]; q < s.cls[i+1]; q++ {
			sl := s.clidx[q]
			if s.sA[sl] <= 0 {
				continue
			}
			sacc += (s.sA[sl] / s.Aref) * (s.dref / s.sd[sl])
		}
		s.fsum[i] = sacc
	}
	for i := 0; i < s.NC; i++ {
		for q := s.cls[i]; q < s.cls[i+1]; q++ {
			sl := s.clidx[q]
			if s.sli[sl] != i { // canonical: apply from the i side
				continue
			}
			if s.sA[sl] <= 0 {
				continue
			}
			j := s.slj[sl]
			si, sj := s.fsum[i], s.fsum[j]
			if si <= 1e-12 || sj <= 1e-12 {
				continue
			}
			w := (s.sA[sl] / s.Aref) * (s.dref / s.sd[sl])
			tau := P.FieldJ * w / math.Sqrt(si*sj) * dt
			ss, cc := math.Sincos(tau)
			a1i, a2i, a1j, a2j := s.fa1[i], s.fa2[i], s.fa1[j], s.fa2[j]
			s.fa1[i] = cc*a1i + ss*a2j
			s.fa2[i] = cc*a2i - ss*a1j
			s.fa1[j] = cc*a1j + ss*a2i
			s.fa2[j] = cc*a2j - ss*a1i
			if P.P1Meter != 0 {
				// P1 site 5: field hop — energy the rotation moved into j
				mj := (s.fa1[j]*s.fa1[j] + s.fa2[j]*s.fa2[j]) - (a1j*a1j + a2j*a2j)
				m := mj * s.sd[sl]
				s.p1fd[0] += m * s.sux[sl]
				s.p1fd[1] += m * s.suy[sl]
				s.p1fd[2] += m * s.suz[sl]
			}
		}
	}
	for i := 0; i < s.NC; i++ {
		s.Ee[i] = s.fa1[i]*s.fa1[i] + s.fa2[i]*s.fa2[i]
	}

	// pass 6: dense clock + beat-gated conversion (cellfab.c:3448), serial
	for i := 0; i < s.NC; i++ {
		if s.roughq[i] > 0 {
			s.qatomDiag(0, s.atomsW(s.w2e[i], s.w1e[i]), s.roughq[i], i, s.Em[i])
			s.roughq[i] = 0
		}
	}
	s.bhNin = 0
	for i := 0; i < s.NC; i++ {
		if P.BhR > 0 {
			hdx := s.wr(s.px[i] - 0.5*P.L)
			hdy := s.wr(s.py[i] - 0.5*P.L)
			hdz := s.wr(s.pz[i] - 0.5*P.L)
			if P.BhSep > 0 {
				// FLOW: two centres at cx +- sep/2; use the nearer
				dxa := s.wr(s.px[i] - (0.5*P.L - 0.5*P.BhSep))
				dxb := s.wr(s.px[i] - (0.5*P.L + 0.5*P.BhSep))
				if math.Abs(dxa) < math.Abs(dxb) {
					hdx = dxa
				} else {
					hdx = dxb
				}
			}
			if hdx*hdx+hdy*hdy+hdz*hdz < P.BhR*P.BhR {
				s.bhNin++
				if s.Ee[i] > 0 {
					s.ehAdd(s.Ee[i])
					s.bhEatF += s.Ee[i]
					s.Ee[i] = 0
					s.fa1[i] = 0
					s.fa2[i] = 0
				}
				if s.Em[i] > 0 {
					s.ehAdd(s.Em[i])
					s.bhEatM += s.Em[i]
					s.Em[i] = 0
				}
				avail := s.Es[i] - P.EsFloor
				if avail > 0 {
					d := P.BhK * avail * dt
					if d > avail {
						d = avail
					}
					s.Es[i] -= d
					s.ehAdd(d)
					s.bhEatS += d
				}
				continue // pitchless: no clock, no beat, no door
			}
		}
		s.th2[i] = math.Mod(s.th2[i]+s.w2e[i]*dt, TwoPi)
		s.cbeta[i] += (s.w1e[i] - s.w2e[i]) * dt
		beatFire := false
		if s.cbeta[i] >= TwoPi {
			s.cbeta[i] -= TwoPi
			beatFire = true
		} else if s.cbeta[i] <= -TwoPi {
			s.cbeta[i] += TwoPi
			beatFire = true
		}
		if beatFire {
			econdI := P.ECond
			if P.WfOn != 0 {
				// WORKFN W-A: presence threshold (WORKFN.md §1)
				if s.Em[i] >= P.WfFloor {
					econdI = 0
				} else {
					econdI = P.WfFar
				}
			}
			if s.scond[i] != 0 {
				econdI = 0.0
			}
			if s.Ee[i] > econdI {
				d1 := P.FConv * (s.Ee[i] - econdI)
				eF := s.A0eff * s.w1e[i] / TwoPi
				eD := s.A0eff * s.w2e[i] / TwoPi
				d1 = s.atomsFire(d1, eF, eD, &s.qcnvD[i])
				d1 = s.atomsClamp(d1, 0.98*s.Ee[i], eF, eD, &s.qcnvD[i])
				if d1 > 0 {
					s.condTotal += d1
					if s.tag[i] != 0 {
						s.ctCond += d1
					}
					if s.scond[i] != 0 && s.nclick < clickMax {
						s.clickT[s.nclick] = s.simT
						s.clickY[s.nclick] = s.py[i]
						s.clickE[s.nclick] = d1
						s.nclick++
					}
					s.qatomDiag(1, s.atomsW(s.w1e[i], s.w2e[i]), d1, i, s.Em[i])
					dsp := P.SPull * d1
					avail := s.Es[i] - P.EsFloor
					if avail < 0 {
						avail = 0
					}
					if dsp > avail {
						dsp = avail
					}
					fac := 0.0
					if s.Ee[i] > 0 {
						fac = math.Sqrt((s.Ee[i] - d1) / s.Ee[i])
					}
					s.fa1[i] *= fac
					s.fa2[i] *= fac
					s.Ee[i] -= d1
					s.Es[i] -= dsp
					s.Em[i] += d1 + dsp
					if P.QpPhase > 0 {
						// QUENCH-3: phase crosses the door (mix <= 1
						// by construction: d1 <= Em after the add)
						aph := math.Atan2(s.fa2[i], s.fa1[i])
						mix := P.QpPhase * d1 / s.Em[i]
						s.th2[i] = math.Mod(s.th2[i]+mix*wrapPi(aph-s.th2[i])+
							8.0*TwoPi, TwoPi)
					}
				}
			}
			tot := s.Em[i] + s.Ee[i]
			if tot > P.Cap {
				d2 := P.FEvap * (tot - P.Cap)
				eF2 := s.A0eff * s.w2e[i] / TwoPi
				eD2 := s.A0eff * s.w1e[i] / TwoPi
				d2 = s.atomsFire(d2, eF2, eD2, &s.qcnvF[i])
				d2 = s.atomsClamp(d2, s.Em[i], eF2, eD2, &s.qcnvF[i])
				if d2 > 0 {
					s.evapTotal += d2
					s.qatomDiag(0, s.atomsW(s.w2e[i], s.w1e[i]), d2, i, s.Em[i])
					bs := d2 * P.SPull / (1.0 + P.SPull)
					s.backsTotal += bs
					if s.tag[i] != 0 {
						s.ctEvap += d2
						s.ctBacks += bs
					}
					s.Em[i] -= d2
					s.fieldInject(i, d2-bs)
					s.Es[i] += bs
				}
			}
			// v91 RADIANCE (laws_V3r candidate A) — mirrors the C
			// kernel's pass-6 block exactly; k_rad=0 => V2g.
			if P.KRad > 0 && s.Em[i] > 0 {
				xr := s.Em[i] / P.Cap
				comp := 1.0
				if P.RadClock == 0 {
					comp = P.W2 / s.w2e[i]
				}
				dr := P.KRad * P.Cap * math.Pow(xr, P.PRad) * comp
				eFr := s.A0eff * s.w2e[i] / TwoPi
				eDr := s.A0eff * s.w1e[i] / TwoPi
				dr = s.atomsFire(dr, eFr, eDr, &s.qcnvF[i])
				dr = s.atomsClamp(dr, s.Em[i], eFr, eDr, &s.qcnvF[i])
				if dr > 0 {
					s.radTotal += dr
					s.qatomDiag(0, s.atomsW(s.w2e[i], s.w1e[i]), dr, i, s.Em[i])
					bsr := dr * P.SPull / (1.0 + P.SPull)
					s.backsTotal += bsr
					if s.tag[i] != 0 {
						s.ctEvap += dr
						s.ctBacks += bsr
					}
					s.Em[i] -= dr
					s.fieldInject(i, dr-bsr)
					s.Es[i] += bsr
				}
			}
		}
	}

	// pass 7: plane re-alignment + tumble (cellfab.c:3515)
	sq := P.SigmaTumble * math.Sqrt(dt)
	if sq > 0 {
		for i := 0; i < s.NC; i++ {
			for k := 0; k < 6; k++ {
				s.rngbuf[6*i+k] = s.grand()
			}
		}
	}
	for i := 0; i < s.NC; i++ {
		s.nsnap[6*i+0] = s.n1x[i]
		s.nsnap[6*i+1] = s.n1y[i]
		s.nsnap[6*i+2] = s.n1z[i]
		s.nsnap[6*i+3] = s.n2x[i]
		s.nsnap[6*i+4] = s.n2y[i]
		s.nsnap[6*i+5] = s.n2z[i]
	}
	for i := 0; i < s.NC; i++ {
		for c := 0; c < 2; c++ {
			var mx, my, mz float64
			if c == 0 {
				mx, my, mz = s.n1x[i], s.n1y[i], s.n1z[i]
			} else {
				mx, my, mz = s.n2x[i], s.n2y[i], s.n2z[i]
			}
			ax, ay, az, fs := 0.0, 0.0, 0.0, 0.0
			if c == 1 {
				for q := s.cls[i]; q < s.cls[i+1]; q++ {
					sl := s.clidx[q]
					fl := s.sflux[sl]
					if fl <= 0 {
						continue
					}
					o := s.slj[sl]
					if s.sli[sl] != i {
						o = s.sli[sl]
					}
					no := s.nsnap[6*o+3 : 6*o+6]
					sgn := 1.0
					if mx*no[0]+my*no[1]+mz*no[2] < 0 {
						sgn = -1.0
					}
					du := mx*s.sux[sl] + my*s.suy[sl] + mz*s.suz[sl]
					ax += fl * (sgn*no[0] - du*s.sux[sl])
					ay += fl * (sgn*no[1] - du*s.suy[sl])
					az += fl * (sgn*no[2] - du*s.suz[sl])
					fs += fl
				}
			}
			vx, vy, vz := mx, my, mz
			if fs > 0 {
				// stable form: |ax/fs| <= 2 by construction, but
				// kappa_align*dt/fs overflows to inf on subnormal-dust
				// flux (kernel/freecell.c pass 7) — divide first
				w := P.KappaAlign * dt
				vx += w * (ax / fs)
				vy += w * (ay / fs)
				vz += w * (az / fs)
			}
			if sq > 0 {
				g := s.rngbuf[6*i+3*c : 6*i+3*c+3]
				vx += sq * g[0]
				vy += sq * g[1]
				vz += sq * g[2]
			}
			nn := math.Sqrt(vx*vx + vy*vy + vz*vz)
			if nn > 1e-12 {
				if c == 0 {
					s.n1x[i], s.n1y[i], s.n1z[i] = vx/nn, vy/nn, vz/nn
				} else {
					s.n2x[i], s.n2y[i], s.n2z[i] = vx/nn, vy/nn, vz/nn
				}
			}
		}
	}
}
