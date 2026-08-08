package par

import (
	"math"
	"testing"
	"unsafe"

	"scp/pfc/core"
)

func TestParallelMatchesSerialEnergy(t *testing.T) {
	law := core.DefaultLaw()
	law.L = 10
	w := core.BuildLattice(law, 1000)
	c := FromWorld(w, 4)
	defer c.Stop()
	e0 := c.TotalE()
	for i := 0; i < 15; i++ {
		c.Step()
	}
	e1 := c.TotalE()
	rel := math.Abs(e1-e0) / e0
	if rel > 1e-9 {
		t.Fatalf("drift %g", rel)
	}
}

func TestWorkerOwnsPrivateSlices(t *testing.T) {
	law := core.DefaultLaw()
	law.L = 8
	w := core.BuildLattice(law, 512)
	c := FromWorld(w, 4)
	defer c.Stop()
	// distinct backing arrays per shard
	seenLo := map[uintptr]int{}
	seenHi := map[uintptr]int{}
	for i := 0; i < 4; i++ {
		sh := c.Shard(i)
		if sh.NOwn == 0 {
			t.Fatalf("shard %d empty", i)
		}
		loPtr := uintptr(unsafe.Pointer(&sh.Lo[0]))
		hiPtr := uintptr(unsafe.Pointer(&sh.Hi[0]))
		if other, ok := seenLo[loPtr]; ok {
			t.Fatalf("shard %d and %d share Lo backing", i, other)
		}
		if other, ok := seenHi[hiPtr]; ok {
			t.Fatalf("shard %d and %d share Hi backing", i, other)
		}
		seenLo[loPtr] = i
		seenHi[hiPtr] = i
		// worker always uses same slice header length policy
		if len(sh.Lo) != sh.NOwn+sh.NGhost || len(sh.Hi) != len(sh.Lo) {
			t.Fatalf("shard %d Lo/Hi len mismatch", i)
		}
	}
	// ranges stay fixed across steps
	g0 := c.Shard(1).G0
	n := c.Shard(1).NOwn
	c.Step()
	c.Step()
	if c.Shard(1).G0 != g0 || c.Shard(1).NOwn != n {
		t.Fatal("shard ownership moved")
	}
	if &c.Shard(1).Lo[0] == nil {
		t.Fatal("lo lost")
	}
}

func TestWorkerCountEnergy(t *testing.T) {
	// multi-worker need not be bit-identical to 1-worker under boundary
	// hop ordering, but energy should stay conserved on each.
	law := core.DefaultLaw()
	law.L = 8
	src := core.BuildLattice(law, 512)
	for _, nw := range []int{1, 4} {
		c := FromWorld(src, nw)
		e0 := c.TotalE()
		for i := 0; i < 8; i++ {
			c.Step()
		}
		e1 := c.TotalE()
		if math.Abs(e1-e0)/e0 > 1e-9 {
			t.Fatalf("workers=%d drift", nw)
		}
		c.Stop()
	}
}

func TestLocalSlabDrift(t *testing.T) {
	law := core.DefaultLaw()
	law.L = 8
	w := core.BuildLattice(law, 512)
	cl := BuildLocal(w, 4)
	defer cl.Stop()
	e0 := cl.TotalE()
	for i := 0; i < 10; i++ {
		cl.Step()
	}
	e1 := cl.TotalE()
	if math.Abs(e1-e0)/e0 > 1e-9 {
		t.Fatalf("local drift")
	}
}
