package par

import "scp/pfc/core"

func passSpaceWantsLive(sh *Shard, law core.Law) {
	dt := law.Dt
	for _, si := range sh.Live {
		s := int(si)
		sl := &sh.Slots[s]
		i, j := sl.I, sl.J
		pi := sh.Lo[i].Es + law.SDisp*(sh.Lo[i].Em+sh.Lo[i].Ee)
		pj := sh.Lo[j].Es + law.SDisp*(sh.Lo[j].Em+sh.Lo[j].Ee)
		dp := pi - pj
		if dp == 0 {
			sl.Swl = 0
			continue
		}
		w := (sl.A / law.Aref) * (law.Dref / sl.D)
		sl.Swl = law.SK * dt * w * dp
	}
}

func passSpaceScaleActive(sh *Shard, law core.Law) {
	for i := 0; i < sh.NOwn; i++ {
		sh.sscl[i] = 1
	}
	for _, ai := range sh.Active {
		i := int(ai)
		rq := 0.0
		for q := sh.Cls[i]; q < sh.Cls[i+1]; q++ {
			sl := &sh.Slots[sh.Inc[q]]
			f := sl.Swl
			if f == 0 {
				continue
			}
			if (f > 0 && sl.I == int32(i)) || (f < 0 && sl.J == int32(i)) {
				if f < 0 {
					rq -= f
				} else {
					rq += f
				}
			}
		}
		avail := sh.Lo[i].Es - law.EsFloor
		if avail <= 0 {
			sh.sscl[i] = 0
		} else if rq > avail {
			sh.sscl[i] = avail / rq
		}
	}
}

func passSpaceApplyActive(sh *Shard, law core.Law) {
	for _, ai := range sh.Active {
		i := int(ai)
		de := 0.0
		for q := sh.Cls[i]; q < sh.Cls[i+1]; q++ {
			sl := &sh.Slots[sh.Inc[q]]
			f := sl.Swl
			if f == 0 {
				continue
			}
			var src int32
			if f > 0 {
				src = sl.I
			} else {
				src = sl.J
			}
			mag := abs(f) * sh.sscl[src]
			if src == int32(i) {
				de -= mag
			} else {
				de += mag
			}
		}
		sh.Lo[i].Es += de
		if sh.Lo[i].Es < law.EsFloor {
			sh.Lo[i].Es = law.EsFloor
		}
	}
}

func passFieldSumLive(sh *Shard, law core.Law) {
	for i := range sh.fsum {
		sh.fsum[i] = 0
	}
	for _, si := range sh.Live {
		sl := &sh.Slots[si]
		w := (sl.A / law.Aref) * (law.Dref / sl.D)
		sh.fsum[sl.I] += w
		sh.fsum[sl.J] += w
	}
}

func passHopWeightLive(sh *Shard, law core.Law) {
	for _, si := range sh.Live {
		s := int(si)
		sl := &sh.Slots[s]
		i, j := int(sl.I), int(sl.J)
		si_, sj := sh.fsum[i], sh.fsum[j]
		if si_ <= 1e-12 || sj <= 1e-12 {
			sl.WField = 0
			continue
		}
		w := (sl.A / law.Aref) * (law.Dref / sl.D)
		sl.WField = law.FieldJ * w / sqrt(si_*sj) * law.Dt
	}
}

func hopFa(fa1, fa2 []float64, i, j int, tau float64) {
	cc, ss := cos(tau), sin(tau)
	a1i, a2i := fa1[i], fa2[i]
	a1j, a2j := fa1[j], fa2[j]
	fa1[i] = cc*a1i + ss*a2j
	fa2[i] = cc*a2i - ss*a1j
	fa1[j] = cc*a1j + ss*a2i
	fa2[j] = cc*a2j - ss*a1i
}

func forceBoundaryOwned(sh *Shard, s int, law core.Law) {
	sl := &sh.Slots[s]
	if sl.Alive == 0 {
		return
	}
	d0 := sh.Lo[sl.I].Cr + sh.Lo[sl.J].Cr
	if sl.D >= d0 || sl.D <= 1e-12 {
		return
	}
	mag := law.KRep * (d0 - sl.D)
	fx, fy, fz := mag*sl.Ux, mag*sl.Uy, mag*sl.Uz
	if int(sl.I) < sh.NOwn {
		sh.Hi[sl.I].Fx -= fx
		sh.Hi[sl.I].Fy -= fy
		sh.Hi[sl.I].Fz -= fz
	}
	if int(sl.J) < sh.NOwn {
		sh.Hi[sl.J].Fx += fx
		sh.Hi[sl.J].Fy += fy
		sh.Hi[sl.J].Fz += fz
	}
}

func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}
