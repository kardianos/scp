package par

import (
	"time"
)

// Step — zero heap allocations after FromWorld (see alloc_test.go).
//
// N shard goroutines were started in FromWorld; Step only sets a phase
// id and barriers (see workers.go). No pool abstraction, no go per phase.
func (c *Cluster) Step() {
	law := c.Law
	timed := phaseOn.Load()
	mark := func(id int, t0 time.Time) {
		if timed {
			phaseAdd(id, time.Since(t0))
		}
	}

	t0 := time.Now()
	c.runPhase(phPullGeom)
	mark(PhPullGeom, t0)

	if law.SK > 0 {
		t0 = time.Now()
		c.runPhase(phExchSsclPub)
		c.runPhase(phExchSsclPull)
		mark(PhExch, t0)

		t0 = time.Now()
		c.runPhase(phApplyCrPub)
		c.runPhase(phFieldHops)
		mark(PhWork, t0)
	} else {
		t0 = time.Now()
		c.runPhase(phFieldOnly)
		mark(PhWork, t0)
	}

	t0 = time.Now()
	for e := 0; e < len(c.bndGi); e++ {
		sh := &c.shards[c.bndShard[e]]
		s := int(c.bndLocal[e])
		sl := &sh.Slots[s]
		if sl.Alive == 0 || sl.WField == 0 || !sh.slotLive(s) {
			continue
		}
		tau := sl.WField
		if tau > 0.5 {
			tau = 0.5
		}
		hopFa(c.pubFa1, c.pubFa2, int(c.bndGi[e]), int(c.bndGj[e]), tau)
	}
	mark(PhBoundary, t0)

	t0 = time.Now()
	c.runPhase(phFinish)
	mark(PhFinish, t0)
}
