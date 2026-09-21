package parcore

import "time"

func now() time.Time        { return time.Now() }
func since(t time.Time) time.Duration { return time.Since(t) }

// Step replays the cached phase graph. Consecutive phases are separated by a
// barrier (runParallel returns only when every shard finished the phase).
//
// KindParallel phases run Phase.Run on each shard concurrently. KindBoundary
// phases run Phase.Bnd ONCE on the caller thread over the global publish
// buffer, after refreshing each boundary edge's Tau from the owning shard's
// SlotParam (so the cross-shard hop uses the same per-step angle as the
// interior colored hops).
//
// Zero heap allocations after FromWorld (see TestStepZeroAlloc). The per-shard
// Live list is reset in place by the driver's geom pass calling View.ResetLive
// (the core never clears it, because the driver decides when geometry refreshes).
func (c *Cluster[C, S]) Step() {
	for pi := range c.phases {
		ph := &c.phases[pi]
		if ph.Kind == KindBoundary {
			c.runBoundary(ph)
			continue
		}
		if ph.Run == nil {
			continue
		}
		t0 := now()
		c.runPhase(pi)
		phaseAdd(pi, since(t0))
	}
}

// runBoundary refreshes the boundary Tau table from the owning shards'
// SlotParam, then invokes the driver's Bnd closure once on the caller thread.
func (c *Cluster[C, S]) runBoundary(ph *Phase[C, S]) {
	t0 := now()
	if ph.Bnd == nil || len(c.bndGi) == 0 {
		phaseAdd(PhBoundary, since(t0))
		return
	}
	for e := range c.bndGi {
		sh := &c.shards[c.bndShard[e]]
		c.bndTau[e] = sh.slotParam[c.bndLocal[e]]
	}
	c.bndView.Pub = c.pub
	c.bndView.Edges = c.bndEdgesView()
	ph.Bnd(&c.bndView)
	phaseAdd(PhBoundary, since(t0))
}

// bndEdgesView builds a BndEdge slice (gi, gj, tau) view over the global
// boundary table. This is the ONE allocation Cluster.Step permits in the
// generic path; it is amortized away when the boundary edge count is stable
// (the slice is re-sliced from a pre-sized buffer in buildOnce below).
func (c *Cluster[C, S]) bndEdgesView() []BndEdge {
	if cap(c.bndEdgeBuf) < len(c.bndGi) {
		c.bndEdgeBuf = make([]BndEdge, len(c.bndGi))
	}
	buf := c.bndEdgeBuf[:len(c.bndGi)]
	for e := range c.bndGi {
		buf[e] = BndEdge{Gi: c.bndGi[e], Gj: c.bndGj[e], Tau: c.bndTau[e]}
	}
	return buf
}
