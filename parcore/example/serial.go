package example

import "scp/parcore"

// SerialSim is the serial baseline: the same passes as Driver, run on one
// global array with no sharding, no ghosts, no barriers. It exists so the
// bench can measure parcore's overhead vs a truly serial run of identical
// physics. The pass functions are shared with the parallel driver — this is
// the proof that parcore imposes no physics-side overhead and that the driver
// is reusable without the core.
type SerialSim struct {
	Law      Law
	Cells    []Cell
	Slots    []Slot
	SI, SJ   []int32
	Cls, Inc []int32
	Colors   [][]int32
	Live     []int32
	Fsum     []float64
	Param    []float64
}

// NewSerial builds a SerialSim from a World (copies state, builds global CSR
// + edge coloring over ALL slots — there are no boundary edges in serial).
func NewSerial(w *parcore.World[Cell, Slot], law Law) *SerialSim {
	nc := len(w.Cells)
	s := &SerialSim{
		Law:   law,
		Cells: append([]Cell(nil), w.Cells...),
		Slots: append([]Slot(nil), w.Slots...),
		SI:    append([]int32(nil), w.SI...),
		SJ:    append([]int32(nil), w.SJ...),
		Live:  make([]int32, 0, len(w.SI)),
		Fsum:  make([]float64, nc),
		Param: make([]float64, len(w.SI)),
	}
	alive := make([]bool, len(w.SI))
	for i := range alive {
		alive[i] = true
	}
	s.Cls, s.Inc = parcore.BuildCSRPublic(s.SI, s.SJ, alive, nc)
	all := make([]int32, len(w.SI))
	for i := range all {
		all[i] = int32(i)
	}
	s.Colors = parcore.EdgeColorPublic(s.SI, s.SJ, alive, all, nc)
	return s
}

// Step runs one serial timestep. Same pass ordering as Driver.Step but with
// no halo exchanges and no boundary round (every slot is interior).
func (s *SerialSim) Step() {
	law := s.Law
	nc := len(s.Cells)
	PassPitch(law, s.Cells, 0, nc)
	PassGeom(law, s.Cells, s.Slots, s.SI, s.SJ, &s.Live)
	PassSpaceWants(law, s.Cells, s.Slots, s.SI, s.SJ, s.Live)
	PassSpaceScale(law, s.Cells, s.Slots, s.SI, s.SJ, s.Cls, s.Inc, nc)
	PassSpaceApply(law, s.Cells, s.Slots, s.SI, s.SJ, s.Cls, s.Inc, nc)
	PassPitch(law, s.Cells, 0, nc)
	PassFieldPrecess(law, s.Cells, 0, nc)
	PassFieldSum(s.Slots, s.SI, s.SJ, s.Fsum, s.Live)
	PassHopWeight(law, s.Slots, s.SI, s.SJ, s.Fsum, s.Param, s.Live)
	for _, bucket := range s.Colors {
		for _, si32 := range bucket {
			ApplyHop(s.Cells, s.SI, s.SJ, s.Param, int(si32))
		}
	}
	PassCommit(s.Cells, 0, nc)
}

// TotalE sums Es+Em+Ee over all cells.
func (s *SerialSim) TotalE() float64 {
	var e float64
	for i := range s.Cells {
		e += s.Cells[i].Es + s.Cells[i].Em + s.Cells[i].Ee
	}
	return e
}
