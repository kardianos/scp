package example

import (
	"math"
	"testing"
)

// TestSerialConservation: the serial baseline conserves total E (Es+Em+Ee).
func TestSerialConservation(t *testing.T) {
	law := DefaultLaw()
	w := BuildWorld(law, 1000)
	s := NewSerial(w, law)
	e0 := s.TotalE()
	for i := 0; i < 30; i++ {
		s.Step()
	}
	e1 := s.TotalE()
	if rel := math.Abs(e1-e0) / e0; rel > 1e-9 {
		t.Fatalf("serial drift %g", rel)
	}
	t.Logf("serial 30 steps rel_drift=%.3e", math.Abs(e1-e0)/e0)
}

// TestSerialZeroAlloc: SerialSim.Step allocates nothing after warmup.
func TestSerialZeroAlloc(t *testing.T) {
	law := DefaultLaw()
	w := BuildWorld(law, 1000)
	s := NewSerial(w, law)
	for i := 0; i < 3; i++ {
		s.Step()
	}
	avg := testing.AllocsPerRun(30, func() {
		s.Step()
	})
	if avg != 0 {
		t.Fatalf("SerialSim.Step allocated %.2f objs/run (want 0)", avg)
	}
}

// TestPitchMonotoneRadius: higher Es → larger Cr (the load-flattens-pitch law).
func TestPitchMonotoneRadius(t *testing.T) {
	law := DefaultLaw()
	c := []Cell{{Es: 0.5}, {Es: 1.0}, {Es: 2.0}}
	PassPitch(law, c, 0, 3)
	if !(c[0].Cr < c[1].Cr && c[1].Cr < c[2].Cr) {
		t.Fatalf("radius not monotone in Es: %g %g %g", c[0].Cr, c[1].Cr, c[2].Cr)
	}
}
