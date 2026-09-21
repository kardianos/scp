package parcore

import (
	"math"
	"testing"
)

// TestEdgeColorDisjoint confirms each color bucket has endpoint-disjoint slots
// (the race-freedom invariant for dual-endpoint writes).
func TestEdgeColorDisjoint(t *testing.T) {
	si := []int32{0, 1, 2, 0, 1, 0}
	sj := []int32{1, 2, 0, 2, 3, 3}
	alive := make([]bool, len(si))
	for i := range alive {
		alive[i] = true
	}
	interior := make([]int32, len(si))
	for i := range interior {
		interior[i] = int32(i)
	}
	buckets := EdgeColorPublic(si, sj, alive, interior, 4)
	for c, bucket := range buckets {
		seen := map[int32]bool{}
		for _, s := range bucket {
			if seen[si[s]] {
				t.Fatalf("color %d: cell %d in two slots", c, si[s])
			}
			if seen[sj[s]] {
				t.Fatalf("color %d: cell %d in two slots", c, sj[s])
			}
			seen[si[s]] = true
			seen[sj[s]] = true
		}
	}
}

// TestPartitionZSlab confirms z-slab partitioning when Side is set.
func TestPartitionZSlab(t *testing.T) {
	const side = 4
	nc := side * side * side // 64
	starts, counts := partitionCells(nc, 4, side)
	total := 0
	for p := 0; p < 4; p++ {
		// each slab must be whole xy-planes
		if counts[p]%(side*side) != 0 {
			t.Fatalf("slab %d not whole planes: %d", p, counts[p])
		}
		if starts[p]%(side*side) != 0 {
			t.Fatalf("slab %d start not plane-aligned: %d", p, starts[p])
		}
		total += counts[p]
	}
	if total != nc {
		t.Fatalf("partition lost cells: %d != %d", total, nc)
	}
}

// TestCSRSymmetric confirms every alive slot is incident on both endpoints.
func TestCSRSymmetric(t *testing.T) {
	si := []int32{0, 1, 2, 3}
	sj := []int32{1, 2, 3, 0}
	alive := []bool{true, true, true, true}
	cls, inc := BuildCSRPublic(si, sj, alive, 4)
	// each cell should be incident on exactly 2 slots
	for i := 0; i < 4; i++ {
		deg := cls[i+1] - cls[i]
		if deg != 2 {
			t.Fatalf("cell %d degree %d != 2", i, deg)
		}
	}
	_ = inc
	// every slot appears in exactly two cells' incidence lists
	count := make(map[int32]int)
	for i := 0; i < 4; i++ {
		for q := cls[i]; q < cls[i+1]; q++ {
			count[inc[q]]++
		}
	}
	for s := int32(0); s < 4; s++ {
		if count[s] != 2 {
			t.Fatalf("slot %d incident %d != 2", s, count[s])
		}
	}
}

// TestGivens2NormConservation: |a'i|²+|a'j|² == |ai|²+|aj|² to FP precision.
func TestGivens2NormConservation(t *testing.T) {
	cases := []struct{ ai1, ai2, aj1, aj2, tau float64 }{
		{0.3, 0.4, 0.1, 0.2, 0.0},
		{0.3, 0.4, 0.1, 0.2, 0.1},
		{0.3, 0.4, 0.1, 0.2, 0.5},
		{1.0, 0.0, 0.0, 1.0, 0.7853981633974483}, // pi/4
		{0.0, 0.0, 0.0, 0.0, 0.3},                // zero pair
	}
	for i, c := range cases {
		n1i, n2i, n1j, n2j := Givens2(c.ai1, c.ai2, c.aj1, c.aj2, c.tau)
		before := c.ai1*c.ai1 + c.ai2*c.ai2 + c.aj1*c.aj1 + c.aj2*c.aj2
		after := n1i*n1i + n2i*n2i + n1j*n1j + n2j*n2j
		if math.Abs(before-after) > 1e-15*math.Max(1, before) {
			t.Fatalf("case %d: norm not conserved %g -> %g", i, before, after)
		}
	}
}

// TestGivens2Tau0Identity: at tau=0 the rotation is the identity exactly.
func TestGivens2Tau0Identity(t *testing.T) {
	n1i, n2i, n1j, n2j := Givens2(0.3, 0.4, 0.1, 0.2, 0.0)
	if n1i != 0.3 || n2i != 0.4 || n1j != 0.1 || n2j != 0.2 {
		t.Fatal("tau=0 not identity")
	}
}
