package parcore

// View is the per-shard handle a Phase.Run closure operates on. It carries
// the shard's owned+ghost cell slice, the slot channel slice, the core-owned
// endpoint/CSR/color topology, the live-slot list, the per-slot hop-angle
// scratch (SlotParam), and the driver's per-shard Scratch.
//
// All indices into Cells/Slots/SI/SJ/Cls/Inc/Colors/Live/SlotParam are LOCAL
// to this shard. Cells layout: [0:NOwn) owned, [NOwn:NOwn+NGhost) ghosts.
// Ghosts are refreshed by Pull (from the global publish buffer, written by
// other shards' Publish).
//
// A Phase.Run closure must:
//   - write owned cells [0:NOwn) freely (that is this shard's exclusive state);
//   - treat ghosts [NOwn:) as READ-ONLY (they are a frozen copy until the next
//     Pull; the core never mutates them mid-phase);
//   - iterate Colors bucket-by-bucket for dual-endpoint slot writes (each
//     bucket's slots are endpoint-disjoint, so a pairwise rotation is safe);
//   - write SlotParam[s] for any slot whose hop angle it computed, so the
//     cross-shard boundary pass can apply the same angle to boundary edges.
type View[C any, S any] struct {
	G0   int   // global id of owned cell 0
	NOwn int   // count of owned cells
	Cells []C  // [0:NOwn) owned, [NOwn:NOwn+NGhost) ghosts
	Slots []S  // channel data; indices parallel SI/SJ

	// Topology (core-owned; rebuilt only if the graph is dynamic — v1 frozen).
	SI, SJ   []int32 // slot endpoints (LOCAL cell indices)
	Cls, Inc []int32 // CSR: cell -> incident local slot indices
	Colors   [][]int32 // edge-color buckets over INTERIOR slots (both ends owned)
	Live     []int32   // live slot indices (driver marks via SetLive during geom)

	// SlotParam is the per-slot hop angle. The driver's hop-weight pass writes
	// it; the core reads it for cross-shard boundary edges. Kept core-owned so
	// the boundary pass does not reach into the driver's slot type.
	SlotParam []float64

	// BndTouch: owned local cells that are an endpoint of at least one
	// cross-shard edge. After a KindBoundary hop over the publish buffer, a
	// Phase merges the published result back into these owned cells
	// (v.Cells[li] = v.pub[v.G0+int(li)]). Empty on interior-only shards.
	BndTouch []int32

	// Scratch is the driver's per-shard scratch (from ScratchDriver.InitShard,
	// or nil if the driver carries scratch inside C).
	Scratch any

	// ---- exchange internals (core-managed; Phase closures call Publish/Pull) ----
	pub        []C       // global publish buffer, indexed by global cell id
	publishSet []int32   // owned local indices this shard must publish (its face)
	ghostGID   []int32   // global id of ghost cell at local index NOwn+k
	NGhost     int
	box        float64
}

// Publish copies this shard's face cells (publishSet) into the global publish
// buffer. Called as a KindParallel Phase (PublishFaces). Every shard must run
// its Publish before any shard runs Pull — the inter-Phase barrier enforces it.
func (v *View[C, S]) Publish() {
	for _, li := range v.publishSet {
		v.pub[v.G0+int(li)] = v.Cells[li]
	}
}

// Pull refreshes this shard's ghosts from the global publish buffer. Called as
// a KindParallel Phase (PullGhosts) AFTER every shard's Publish (the barrier
// between them guarantees the buffer is settled).
func (v *View[C, S]) Pull() {
	for k, gid := range v.ghostGID {
		v.Cells[v.NOwn+k] = v.pub[gid]
	}
}

// PublishedCell returns the published copy of owned cell li (local index) from
// the global publish buffer. Drivers use it after a KindBoundary hop to merge
// the boundary-updated state back into the owned face cells (v.BndTouch).
func (v *View[C, S]) PublishedCell(li int) C {
	return v.pub[v.G0+li]
}

// MergeBndTouch copies the published copy of every BndTouch owned cell back
// into local storage. Standard post-boundary-hop merge: after the serial
// boundary round over the publish buffer, the owned face cells pick up their
// rotated field amplitude (and identical copies of everything else, which is
// harmless).
func (v *View[C, S]) MergeBndTouch() {
	for _, li := range v.BndTouch {
		v.Cells[li] = v.pub[v.G0+int(li)]
	}
}

// SetLive marks slot s live and appends it to the Live list. The driver's
// geometry pass calls this for every slot with a positive channel (A>0). The
// Live list is reset (zero-alloc) at the start of each geom refresh by the
// driver calling ResetLive; the core never clears it.
func (v *View[C, S]) SetLive(s int) {
	v.Live = append(v.Live, int32(s))
}

// ResetLive truncates the Live list in place (zero heap) so the driver's geom
// pass can repopulate it. Call at the top of each geometry refresh.
func (v *View[C, S]) ResetLive() {
	v.Live = v.Live[:0]
}

// Box returns the periodic box length (for minimum-image geometry).
func (v *View[C, S]) Box() float64 { return v.box }

// BoundaryView is the handle for a KindBoundary Phase. It carries the global
// publish buffer and the list of cross-shard dual-endpoint edges (global ids
// + the per-step hop angle Tau, refreshed by the core from SlotParam). The
// driver's Bnd closure loops Edges and applies the SAME pairwise update used
// inside the colored interior rounds to Pub[gi], Pub[gj].
type BoundaryView[C any] struct {
	Pub   []C
	Edges []BndEdge
}

// BndEdge is one cross-shard dual-endpoint edge.
type BndEdge struct {
	Gi, Gj int32   // global cell ids (endpoints)
	Tau    float64 // per-step hop angle, refreshed from the owning shard's SlotParam
}

// ShardView is a read-only snapshot of a shard returned by Cluster.Shard for
// diagnostics and tests. It is a copy of the public fields; the Cells/Slots
// slices share the shard's backing arrays (so pointers are stable and the
// zero-alloc property holds).
type ShardView[C any, S any] struct {
	ID     int
	G0     int
	NOwn   int
	NGhost int
	Cells  []C
	Slots  []S
}
