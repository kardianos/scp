package par

// Phase ids: Step sets wphase and wakes each shard's goroutine.
const (
	phPullGeom = iota + 1
	phExchSsclPub
	phExchSsclPull
	phApplyCrPub
	phFieldHops
	phFieldOnly // SK==0 path after geom
	phFinish
)

// workerPhase runs one phase on shard id. Must not allocate.
func (c *Cluster) workerPhase(id int, phase int) {
	sh := &c.shards[id]
	law := c.Law
	switch phase {
	case phPullGeom:
		sh.pullGhostsFromPub(c.pubLo)
		corePassPitch(sh, law)
		sh.refreshGeom(law.L)
		if law.SK > 0 {
			passSpaceWantsLive(sh, law)
			passSpaceScaleActive(sh, law)
		}

	case phExchSsclPub:
		for _, li := range sh.Publish {
			i := int(li)
			c.pubS[sh.G0+i] = sh.sscl[i]
		}

	case phExchSsclPull:
		for k, gid := range sh.GhostGID {
			sh.sscl[sh.NOwn+k] = c.pubS[gid]
		}

	case phApplyCrPub:
		passSpaceApplyActive(sh, law)
		corePassPitch(sh, law)
		sh.publishFaceCr(c.pubS)

	case phFieldHops:
		sh.pullGhostCr(c.pubS)
		corePassFieldPrecess(sh, law)
		passFieldSumLive(sh, law)
		passHopWeightLive(sh, law)
		for col := 0; col < len(sh.buckets); col++ {
			for _, si := range sh.buckets[col] {
				s := int(si)
				if !sh.slotLive(s) {
					continue
				}
				coreApplyHopOne(sh, s)
			}
		}
		for _, li := range sh.BndTouch {
			i := int(li)
			c.pubFa1[sh.G0+i] = sh.Hi[i].Fa1Next
			c.pubFa2[sh.G0+i] = sh.Hi[i].Fa2Next
		}

	case phFieldOnly:
		corePassFieldPrecess(sh, law)
		passFieldSumLive(sh, law)
		passHopWeightLive(sh, law)
		for col := 0; col < len(sh.buckets); col++ {
			for _, si := range sh.buckets[col] {
				s := int(si)
				if !sh.slotLive(s) {
					continue
				}
				coreApplyHopOne(sh, s)
			}
		}
		for _, li := range sh.BndTouch {
			i := int(li)
			c.pubFa1[sh.G0+i] = sh.Hi[i].Fa1Next
			c.pubFa2[sh.G0+i] = sh.Hi[i].Fa2Next
		}

	case phFinish:
		for _, li := range sh.BndTouch {
			i := int(li)
			sh.Hi[i].Fa1Next = c.pubFa1[sh.G0+i]
			sh.Hi[i].Fa2Next = c.pubFa2[sh.G0+i]
		}
		corePassCommitField(sh)
		for i := 0; i < sh.NOwn; i++ {
			sh.Hi[i].Fx, sh.Hi[i].Fy, sh.Hi[i].Fz = 0, 0, 0
		}
		for col := 0; col < len(sh.buckets); col++ {
			for _, si := range sh.buckets[col] {
				s := int(si)
				if !sh.slotLive(s) {
					continue
				}
				coreForceOne(sh, s, law)
			}
		}
		for _, bi := range sh.Boundary {
			s := int(bi)
			if sh.slotLive(s) {
				forceBoundaryOwned(sh, s, law)
			}
		}
		corePassMotion(sh, law)
		corePassClock(sh, law)
		sh.publishFace(c.pubLo)
	}
}
