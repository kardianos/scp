package example

import (
	"math"

	"scp/parcore"
)

// Driver implements parcore.Driver[Cell, Slot]. It is a value struct holding
// the immutable Law; Step() returns the cached phase graph. Per-shard mutable
// state lives on the View (fsum scratch is allocated lazily and cached on the
// driver's per-shard scratch map — but since this driver keeps fsum in the
// View via a ScratchDriver hook would be cleaner; for v1 we store fsum in a
// small per-shard slice created in InitShard).
type Driver struct {
	Law Law
}

// fsumScratch is the per-shard field-sum buffer (sized to nOwn+nGhost).
type fsumScratch struct {
	fsum []float64
}

// InitShard sizes the field-sum buffer once per shard (zero-alloc on Step).
func (d *Driver) InitShard(nOwn, nGhost, nSlot int) any {
	return &fsumScratch{fsum: make([]float64, nOwn+nGhost)}
}

// Step returns the 4-phase timestep graph (3 barriers/step — pfc's fused-
// wave target). Each KindParallel phase runs on every shard concurrently;
// consecutive phases are separated by a barrier. Publishes are FUSED into the
// compute phases (no standalone publish phase): each compute phase ends by
// publishing its face into the global buffer, so the barrier that closes it
// also settles the halo for the next phase's Pull.
//
// Phase layout:
//
//	0  start     Pull + Pitch + Geom + Wants + Scale + Publish(Sscl)  (parallel)
//	1  mid       Pull + Apply + Pitch + Precess + Sum + HopWeight +
//	             InteriorHops + Publish(Fa)                          (parallel)
//	2  boundary  BoundaryHops over pub                                (serial)
//	3  finish    MergeFaBack + Commit + Publish(end)                  (parallel)
//
// The Pull at the top of phase 0 reads the buffer that phase 3 of the previous
// step published (settled by the step-boundary barrier). The Pull at the top
// of phase 1 reads the Sscl that phase 0 just published (settled by the
// inter-phase barrier). Both are correct because a full-cell Publish carries
// every field, including the scratch the next phase needs.
func (d *Driver) Step() []parcore.Phase[Cell, Slot] {
	law := d.Law
	return []parcore.Phase[Cell, Slot]{
		{
			Name: "start",
			Kind: parcore.KindParallel,
			Run: func(v *parcore.View[Cell, Slot]) {
				v.Pull() // ghosts fresh from prev step's end publish
				PassPitch(law, v.Cells, 0, v.NOwn)
				PassGeom(law, v.Cells, v.Slots, v.SI, v.SJ, &v.Live)
				PassSpaceWants(law, v.Cells, v.Slots, v.SI, v.SJ, v.Live)
				PassSpaceScale(law, v.Cells, v.Slots, v.SI, v.SJ, v.Cls, v.Inc, v.NOwn)
				v.Publish() // publish face (carries Sscl) for the mid pull
			},
		},
		{
			Name: "mid",
			Kind: parcore.KindParallel,
			Run: func(v *parcore.View[Cell, Slot]) {
				v.Pull() // ghosts now carry fresh Sscl
				PassSpaceApply(law, v.Cells, v.Slots, v.SI, v.SJ, v.Cls, v.Inc, v.NOwn)
				PassPitch(law, v.Cells, 0, v.NOwn)
				PassFieldPrecess(law, v.Cells, 0, v.NOwn)
				fs := v.Scratch.(*fsumScratch).fsum
				PassFieldSum(v.Slots, v.SI, v.SJ, fs, v.Live)
				PassHopWeight(law, v.Slots, v.SI, v.SJ, fs, v.SlotParam, v.Live)
				for _, bucket := range v.Colors {
					for _, si32 := range bucket {
						ApplyHop(v.Cells, v.SI, v.SJ, v.SlotParam, int(si32))
					}
				}
				v.Publish() // publish face (carries Fa) for the boundary round
			},
		},
		parcore.BoundaryHops[Cell, Slot](ApplyHopPub),
		{
			Name: "finish",
			Kind: parcore.KindParallel,
			Run: func(v *parcore.View[Cell, Slot]) {
				v.MergeBndTouch()
				PassCommit(v.Cells, 0, v.NOwn)
				v.Publish() // end-of-step publish (for next step's start Pull)
			},
		},
	}
}

// PhaseNames returns the phase labels for the timer report.
func (d *Driver) PhaseNames() []string {
	return []string{"start", "mid", "<boundary>", "finish"}
}

// BuildWorld constructs a World[C,S] on a cubic lattice with a smooth field
// amplitude seed (so the unitary hops have something to transport).
func BuildWorld(law Law, nApprox int) *parcore.World[Cell, Slot] {
	return parcore.BuildLattice[Cell, Slot](nApprox, law.L, law.R0,
		func(ix, iy, iz int, dx float64) Cell {
			x := (float64(ix) + 0.5) * dx
			y := (float64(iy) + 0.5) * dx
			z := (float64(iz) + 0.5) * dx
			fa1 := 0.05 * math.Sin(0.3*x+0.1*y)
			fa2 := 0.05 * math.Cos(0.2*z+0.1*x)
			var c Cell
			c.X, c.Y, c.Z = x, y, z
			c.Cr = law.R0
			c.Es = law.Es0
			c.W1e = law.W1
			c.Fa1, c.Fa2 = fa1, fa2
			c.Ee = fa1*fa1 + fa2*fa2
			return c
		})
}

// EnergyOf is the per-cell energy measure (Es + Em + Ee) for conservation
// checks and the bench.
func EnergyOf(c Cell) float64 { return c.Es + c.Em + c.Ee }
