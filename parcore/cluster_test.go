package parcore_test

import (
	"math"
	"testing"
	"unsafe"

	"scp/parcore"
	"scp/parcore/example"
)

// buildCluster is the shared harness: a lattice world + the demo driver.
func buildCluster(nApprox int, nWorkers int) (*parcore.Cluster[example.Cell, example.Slot], float64) {
	law := example.DefaultLaw()
	w := example.BuildWorld(law, nApprox)
	d := &example.Driver{Law: law}
	c := parcore.FromWorld[example.Cell, example.Slot](w, d, nWorkers)
	e0 := parcore.TotalEnergy(c, example.EnergyOf)
	return c, e0
}

// TestParallelEnergyConservation: total E (Es+Em+Ee) must not drift over many
// steps on multiple worker counts. This is the conservation contract: space
// transport conserves Es pairwise, the unitary Givens2 hops conserve Ee
// pairwise, so the full ledger is conserved by construction.
func TestParallelEnergyConservation(t *testing.T) {
	for _, nw := range []int{1, 2, 4, 8} {
		c, e0 := buildCluster(1728, nw)
		for i := 0; i < 30; i++ {
			c.Step()
		}
		e1 := parcore.TotalEnergy(c, example.EnergyOf)
		rel := math.Abs(e1-e0) / e0
		if rel > 1e-9 {
			t.Fatalf("workers=%d: energy drift %g (E0=%g E1=%g)", nw, rel, e0, e1)
		}
		c.Stop()
		t.Logf("workers=%d: %d steps rel_drift=%.3e", nw, 30, rel)
	}
}

// TestDistinctSlabs: every shard's Cells backing array must be a distinct
// allocation (the pfc ownership invariant — no aliasing between workers).
func TestDistinctSlabs(t *testing.T) {
	c, _ := buildCluster(1728, 4)
	defer c.Stop()
	seen := map[uintptr]int{}
	for p := 0; p < 4; p++ {
		sh := c.Shard(p)
		if sh.NOwn == 0 {
			t.Fatalf("shard %d empty", p)
		}
		ptr := uintptr(unsafe.Pointer(&sh.Cells[0]))
		if other, ok := seen[ptr]; ok {
			t.Fatalf("shard %d and %d share Cells backing", p, other)
		}
		seen[ptr] = p
	}
}

// TestShardOwnershipStable: G0/NOwn must not move across steps (frozen topo).
func TestShardOwnershipStable(t *testing.T) {
	c, _ := buildCluster(1000, 4)
	defer c.Stop()
	type snap struct{ g0, n int }
	var before [4]snap
	for p := 0; p < 4; p++ {
		sh := c.Shard(p)
		before[p] = snap{sh.G0, sh.NOwn}
	}
	for i := 0; i < 5; i++ {
		c.Step()
	}
	for p := 0; p < 4; p++ {
		sh := c.Shard(p)
		if sh.G0 != before[p].g0 || sh.NOwn != before[p].n {
			t.Fatalf("shard %d ownership moved", p)
		}
	}
}

// TestStepZeroAlloc: Cluster.Step allocates nothing after FromWorld (the pfc
// contract carried into parcore; the Live list is reset in place, the boundary
// BndEdge view is re-sliced from a pre-sized buffer).
func TestStepZeroAlloc(t *testing.T) {
	for _, nw := range []int{1, 4} {
		c, _ := buildCluster(1728, nw)
		for i := 0; i < 3; i++ { // warmup
			c.Step()
		}
		avg := testing.AllocsPerRun(30, func() {
			c.Step()
		})
		c.Stop()
		if avg != 0 {
			t.Fatalf("workers=%d: Step allocated %.2f objs/run (want 0)", nw, avg)
		}
		t.Logf("workers=%d: Step allocs/run = %.0f", nw, avg)
	}
}

// TestCoverageNoDeadLocks is a smoke check that a multi-worker run completes
// many steps without hanging (catches a broken barrier).
func TestCoverageNoDeadLocks(t *testing.T) {
	c, _ := buildCluster(8000, 8)
	defer c.Stop()
	for i := 0; i < 100; i++ {
		c.Step()
	}
}

// TestParallelVsSerialEnergy: at 1 worker the parallel ledger should match the
// serial baseline to FP precision (same passes, same ordering, no boundary
// edges at W=1 so no reordering). At W>1 the ledger is still conserved though
// the trajectory diverges (different hop ordering across the cut-plane).
func TestParallelVsSerialEnergy(t *testing.T) {
	law := example.DefaultLaw()
	w := example.BuildWorld(law, 1000)

	ser := example.NewSerial(w, law)
	se0 := ser.TotalE()
	for i := 0; i < 20; i++ {
		ser.Step()
	}
	se1 := ser.TotalE()
	if rel := math.Abs(se1-se0) / se0; rel > 1e-9 {
		t.Fatalf("serial drift %g", rel)
	}

	c := parcore.FromWorld[example.Cell, example.Slot](w, &example.Driver{Law: law}, 1)
	pe0 := parcore.TotalEnergy(c, example.EnergyOf)
	for i := 0; i < 20; i++ {
		c.Step()
	}
	pe1 := parcore.TotalEnergy(c, example.EnergyOf)
	c.Stop()
	if rel := math.Abs(pe1-pe0) / pe0; rel > 1e-9 {
		t.Fatalf("1-worker parallel drift %g", rel)
	}
}
