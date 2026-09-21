// Package parcore is a generic, pluggable, accelerated parallel core for
// free-cell-shaped simulations.
//
// # The split
//
// The DRIVER owns the physics: a per-cell value type [C], a per-slot channel
// type [S], and an ordered list of Phase funcs that compose one timestep.
// The driver is a plain value struct (law + config; no shared mutable state).
//
// The CORE owns the machinery: sharding into contiguous AoS shards, the
// two-phase publish/compute/apply barriers, ghost-halo exchange, edge-coloring
// for dual-endpoint updates, the CSR topology, the procedural worker pool, and
// timers. The driver never sees a lock, a goroutine, or another shard.
//
// This design generalizes the structurally-isolated pfc experiment
// (scp/pfc/par + scp/pfc/core) into a reusable core. pfc proved the structure
// (AoS cells with no pointers, int32 index topology, two-phase barriers,
// edge-color rounds, z-slab ownership, per-shard ghost halos). parcore makes
// the cell/slot types and the pass list generic so a new simulation plugs in
// without touching the harness.
//
// # What a driver supplies
//
//  1. Concrete cell and slot value types (no pointers inside either).
//  2. A driver value implementing Driver[C,S]: its Step() returns the ordered
//     []Phase that defines one timestep. Phase closures capture the (immutable)
//     driver and operate on a per-shard *View.
//  3. A World[C,S]: the serial global state (one contiguous []C, one []S, and
//     the endpoint/geometry tables) that FromWorld partitions into shards.
//
// # The parallel model
//
// FromWorld partitions cells into Workers z-slabs (or contiguous id ranges).
// Each shard gets its own make()'d []C and []S backing arrays (distinct, see
// TestDistinctSlabs). Each shard also carries ghost cells [NOwn:NOwn+NGhost)
// that mirror remote owned cells. At Cluster.Step, the core runs every Phase
// on every shard in parallel and barriers between consecutive phases. The
// driver triggers a halo refresh by inserting PublishFaces / PullGhosts phases
// (full-cell copy in this version; selective exchange is a documented future
// hook). Dual-endpoint slot updates (e.g. unitary hops) run over edge-color
// buckets so no two updates in one round touch the same cell; cross-shard
// dual-endpoint edges run serially on the caller thread in a KindBoundary
// phase over the global publish buffer.
//
// # Conservation contract
//
// The core makes the free-cell invariants HARD to break by accident:
//
//   - Contiguous value cells (no *C fields): all mutation is &Cells[i] on a
//     slice backing array; there is no background of persistent re-valued
//     pointers.
//   - Edge-coloring of dual-endpoint slot ops: a KindSlotColor phase iterates
//     color buckets whose slot sets are endpoint-disjoint, so a pairwise
//     rotation is race-free and its two-cell norm preservation is exact.
//   - Givens2 / Givens2Buf: the blessed pairwise unitary primitive (the
//     v93/freecell pass-F/pass-U rotation). Drivers SHOULD compose dual-
//     endpoint updates from it so conservation is a theorem of the update,
//     not a hope.
//   - Full Es/Ee/Em ledger check (ConservationCheck) in the test suite.
//
// # What is NOT here (v1 seams, see DESIGN.md)
//
//   - Dynamic topology: the roster/graph is frozen at FromWorld. The seam is
//     named (the Topology interface + a Refresh step the core would call); a
//     future DynamicTopology implements birth/death + spatial re-neighboring.
//   - Selective (per-field) ghost exchange: v1 publishes whole cells. pfc's
//     Fa1Next-only and Cr-only exchanges are the bandwidth optimization.
//   - Full v93 law parity / C≡Go: the bundled driver is a freecell-shaped
//     demo (pitch + space transport + field precession + unitary hops), NOT a
//     port of the law table. It exists to exercise the machinery.
//   - Bit-identical serial Givens ordering across worker counts: parallel
//     colored application is physically related to the serial canonical sweep
//     but is not float-identical; the ledger is conserved either way.
package parcore

import "sync"

// PhaseKind classifies how the core schedules a Phase.
type PhaseKind int

const (
	// KindParallel runs Phase.Run once per shard, concurrently across shards,
	// and barriers when every shard has returned. This is the common case:
	// per-cell passes, CSR gathers, colored interior hops, and the
	// Publish/Pull exchange phases.
	KindParallel PhaseKind = iota

	// KindBoundary runs Phase.Bnd exactly ONCE on the caller thread, after all
	// shards have published the endpoints of cross-shard edges into the global
	// publish buffer. It applies the same dual-endpoint update used inside the
	// colored interior rounds to the boundary edges (gi, gj, tau). Running it
	// single-threaded avoids a per-edge region mutex on the dense contact
	// graph; boundary edges are a thin set (the cut-plane).
	KindBoundary
)

// Phase is one step in a driver's timestep graph. The driver returns a slice
// of these from Driver.Step(); the core caches it at build and replays it
// every Cluster.Step().
//
// For KindParallel: set Run. It executes on each shard's *View in parallel.
// For KindBoundary: set Bnd. It executes once on the caller thread over the
// global publish buffer and the cross-shard edge list (see BoundaryView).
//
// Every consecutive Phase pair is separated by a barrier, so a halo refresh is
// just a PublishFaces phase followed (after the inter-phase barrier) by a
// PullGhosts phase. The driver composes the graph; the core enforces the sync.
type Phase[C any, S any] struct {
	Name string
	Kind PhaseKind
	// Run executes on one shard. KindParallel only.
	Run func(v *View[C, S])
	// Bnd executes once on the caller thread. KindBoundary only.
	Bnd func(bv *BoundaryView[C])
}

// Driver is the contract between user physics and the parallel core. The
// concrete driver is a value struct (law + config); its Step() is called once
// at build and the returned []Phase is replayed every step.
//
// A driver MAY also implement ScratchDriver to allocate per-shard scratch
// (e.g. Jacobi "next" buffers kept out of the cell value). If it does not,
// View.Scratch is nil and the driver is expected to carry any scratch inside
// the cell value C itself.
type Driver[C any, S any] interface {
	Step() []Phase[C, S]
}

// ScratchDriver is an optional extension of Driver. InitShard is called once
// per shard at build with that shard's owned/ghost/slot counts; the returned
// value is stored on the shard's View.Scratch and the driver's Phase closures
// cast it back to the concrete scratch type. Pre-sizing scratch here is what
// keeps Cluster.Step at zero heap allocations after build.
type ScratchDriver[C any, S any] interface {
	Driver[C, S]
	InitShard(nOwn, nGhost, nSlot int) any
}

// World is the serial global state that FromWorld partitions. It is the
// generic analog of pfc/core.World. The core reads positions/topology from it
// at build time and never mutates it afterwards (the shards own their copies).
//
// SI/SJ are the slot endpoints as GLOBAL cell indices; Alive marks live
// channels. The core rebuilds a local-index CSR + edge coloring per shard.
type World[C any, S any] struct {
	Cells []C
	Slots []S
	SI, SJ []int32
	Alive  []bool
	// Side is the cubic lattice edge when built by BuildLattice (0 otherwise).
	// z-slab ownership uses it to cut whole xy-planes (minimal boundary).
	Side int
	// Box is the periodic box length (for minimum-image geometry on shards).
	Box float64
}

// NC returns len(Cells).
func (w *World[C, S]) NC() int { return len(w.Cells) }

// Cluster is a domain-decomposed simulation: one shard per worker, plus the
// shared global publish buffer and cross-shard boundary edge table. Build with
// FromWorld; drive with Step; park the worker goroutines with Stop.
type Cluster[C any, S any] struct {
	driver  Driver[C, S]
	Workers int
	shards  []shard[C, S]

	// global publish buffer, indexed by global cell id. PublishFaces writes
	// owned face cells here; PullGhosts copies into ghost slots; KindBoundary
	// hops operate on it directly.
	pub []C

	// cross-shard dual-endpoint edges (global ids) + the shard/slot that owns
	// the per-step hop angle (Tau, refreshed each step from SlotParam).
	bndGi      []int32
	bndGj      []int32
	bndShard   []int32
	bndLocal   []int32
	bndTau     []float64
	bndEdgeBuf []BndEdge      // pre-sized BndEdge view (zero-alloc per step)
	bndView    BoundaryView[C] // reused each step so &bv does not escape

	globalNC int
	side     int
	box      float64

	// cached phase graph (from driver.Step(), built once).
	phases []Phase[C, S]

	// current-generation phase index for the worker runtime (workers.go).
	// Set under wmu by runPhase; read by workers after wake. No closure is
	// allocated per step — workers index into c.phases — so Step is zero-alloc.
	curPhase int

	// procedural worker runtime (see workers.go).
	wmu   sync.Mutex
	wcond *sync.Cond
	wgen  uint64
	wdone int
	wquit bool
}

// Shard returns worker id's permanent shard. Stored value (not pointer) so
// callers get a stable snapshot; cell slices share the shard's backing arrays.
func (c *Cluster[C, S]) Shard(id int) ShardView[C, S] {
	return c.shards[id].public()
}

// Workers count, GlobalNC total cells, NBoundary cross-shard edges.
func (c *Cluster[C, S]) GlobalNC() int     { return c.globalNC }
func (c *Cluster[C, S]) NBoundary() int    { return len(c.bndGi) }
func (c *Cluster[C, S]) Box() float64      { return c.box }
