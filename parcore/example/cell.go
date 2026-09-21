// Package example is a minimal but real free-cell-shaped driver plugged into
// parcore. It is NOT a port of the v93 law table — it exists to exercise every
// part of the parallel machinery (per-cell passes, CSR gather, two-phase
// ghost exchange, edge-colored dual-endpoint unitary hops, cross-shard
// boundary hops, and the publish/compute/apply barrier cycle).
//
// Physics: pitch (radius + load-detuned frequency) → geometry (lens area) →
// space transport (Jacobi wants/scale/apply with a mid-step halo exchange) →
// field precession → unitary field hops (Givens2, conserved pairwise) on
// interior colored rounds + cross-shard boundary rounds. Total E = Es + Em +
// Ee is conserved (Es pairwise by space transport, Ee = Fa1²+Fa2² pairwise by
// the unitary hops).
//
// The pass functions are PURE: they take plain slices and operate on explicit
// [a,b) ranges. The same functions are called (a) directly on flat arrays by
// the serial baseline (SerialStep) and (b) from the Phase closures that
// Driver.Step returns. This proves the parcore abstraction adds no physics
// overhead and that the driver is reusable without the core.
package example

// Cell is the per-site value (AoS, no pointers). Field order is roughly
// hot-first (positions/radius/energies/field) with the space-scale scratch
// (Sscl) at the tail — a single-slab layout chosen for v1 simplicity over the
// pfc CellLo/CellHi split. The split is a documented future optimization.
type Cell struct {
	X, Y, Z float64 // position (geometric)
	Cr      float64 // live radius (from Es)
	Es      float64 // space energy
	Em      float64 // matter energy
	Ee      float64 // field energy (== Fa1²+Fa2², kept in sync)
	Fa1, Fa2 float64 // field amplitude components (unitary hops act here)
	W1e     float64 // effective field frequency (load-detuned)
	Sscl    float64 // space-scale scratch (published mid-step to ghosts)
}

// Slot is the per-channel value (no pointers). Endpoints are NOT here — the
// parcore core owns SI/SJ; the slot carries only channel geometry + scratch.
type Slot struct {
	A          float64 // lens overlap area
	D          float64 // separation
	Ux, Uy, Uz float64 // unit vector i->j
	Swl        float64 // space want (signed magnitude)
}

// Law is the driver's constant table (a freecell-shaped subset, not the V2g
// table verbatim).
type Law struct {
	Dt, L, R0 float64
	Es0, EsFloor float64
	Cap, W1, QDetune float64
	SK, SDisp float64
	Aref, Dref float64
	FieldJ float64
}

// DefaultLaw returns a working subset. R0=0.58 gives positive lens area on a
// unit lattice (R0=0.5 sits exactly at the contact edge → A=0 everywhere).
func DefaultLaw() Law {
	return Law{
		Dt: 0.02, L: 16, R0: 0.58,
		Es0: 1, EsFloor: 0.05,
		Cap: 2.5, W1: 1.65, QDetune: 1.2,
		SK: 0.06, SDisp: 0.3,
		Aref: 1.0, Dref: 1.0,
		FieldJ: 1.8,
	}
}
