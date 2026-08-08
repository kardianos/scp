// Smoke: serial energy, owned-shard layout, parallel drift.
package main

import (
	"fmt"
	"math"
	"os"
	"unsafe"

	"scp/pfc/core"
	"scp/pfc/par"
)

func main() {
	var lo core.CellLo
	var hi core.CellHi
	fmt.Printf("sizeof CellLo=%d CellHi=%d Slot=%d\n",
		unsafe.Sizeof(lo), unsafe.Sizeof(hi), unsafe.Sizeof(core.Slot{}))

	law := core.DefaultLaw()
	law.L = 12
	w := core.BuildLattice(law, 1728)
	fmt.Printf("world NC=%d slots=%d\n", w.NC(), len(w.Slots))

	e0 := w.TotalE()
	for i := 0; i < 10; i++ {
		w.Step()
	}
	e1 := w.TotalE()
	rel := math.Abs(e1-e0) / e0
	fmt.Printf("serial 10 steps: E0=%.6f E1=%.6f rel_drift=%.3e\n", e0, e1, rel)
	if rel > 1e-9 {
		fmt.Fprintf(os.Stderr, "FAIL serial energy drift\n")
		os.Exit(1)
	}

	w2 := core.BuildLattice(law, 1728)
	c4 := par.FromWorld(w2, 4)
	defer c4.Stop()
	// prove distinct allocations
	for p := 0; p < 4; p++ {
		sh := c4.Shard(p)
		fmt.Printf("  shard %d G0=%d own=%d ghost=%d boundary_slots=%d Lo=%p Hi=%p\n",
			p, sh.G0, sh.NOwn, sh.NGhost, len(sh.Boundary), &sh.Lo[0], &sh.Hi[0])
	}
	fmt.Printf("  cross-shard edges (global table)=%d\n", c4.NBoundary())
	eL0 := c4.TotalE()
	for i := 0; i < 10; i++ {
		c4.Step()
	}
	eL1 := c4.TotalE()
	relL := math.Abs(eL1-eL0) / eL0
	fmt.Printf("owned 4-shard 10 steps: E0=%.6f E1=%.6f rel_drift=%.3e\n", eL0, eL1, relL)
	if relL > 1e-9 {
		fmt.Fprintf(os.Stderr, "FAIL shard energy drift\n")
		os.Exit(1)
	}

	fmt.Println("smoke OK")
}
