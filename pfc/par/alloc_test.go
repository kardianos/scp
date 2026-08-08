package par

import (
	"testing"

	"scp/pfc/core"
)

// TestStepZeroAlloc confirms Cluster.Step allocates nothing after init.
func TestStepZeroAlloc(t *testing.T) {
	law := core.DefaultLaw()
	law.L = 10
	// modest size so this stays light if other work is on the machine
	w := core.BuildLattice(law, 1000)
	for _, nw := range []int{1, 4} {
		c := FromWorld(w, nw)
		// warmup: first-touch / worker settle
		for i := 0; i < 5; i++ {
			c.Step()
		}
		avg := testing.AllocsPerRun(50, func() {
			c.Step()
		})
		c.Stop()
		if avg != 0 {
			t.Fatalf("workers=%d: Step allocated %.2f objects/run (want 0)", nw, avg)
		}
		t.Logf("workers=%d: Step allocs/run = %.0f", nw, avg)
	}
}

func TestSerialWorldStepZeroAlloc(t *testing.T) {
	law := core.DefaultLaw()
	law.L = 8
	w := core.BuildLattice(law, 512)
	for i := 0; i < 3; i++ {
		w.Step() // first steps may size sscl/fsum
	}
	avg := testing.AllocsPerRun(50, func() {
		w.Step()
	})
	if avg != 0 {
		t.Fatalf("World.Step allocated %.2f objects/run (want 0)", avg)
	}
}
