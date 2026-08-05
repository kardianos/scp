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

		if wIJ > 0 {
			s.swant[2*sl] = wIJ
		}
		if wJI > 0 {
			s.swant[2*sl+1] = wJI
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
			}
			if s.slem[sslot] <= 1e-17 {
				s.slem[sslot] = 0
				s.slph[sslot] = 0
			} else if take <= 0 {
				s.slph[sslot] = 0
			} else {
				s.slph[sslot] -= 1.0
			}
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
			s.supH[i] = 0
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
				if gg > s.supH[i] {
					s.supH[i] = gg
				}
				if gg > s.supH[j] {
					s.supH[j] = gg
				}
				amp := math.Sqrt(s.ca[i] * s.ca[j])
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
			a := s.ca[i] + ktau*(s.supH[i]-s.ca[i])
			if a < 0 {
				a = 0
			}
			if a > 1 {
				a = 1
			}
			s.ca[i] = a
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
	for i := 0; i < s.NC; i++ {
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
