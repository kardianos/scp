// Package core is the freecell-shaped physics and data layout.
// No concurrency, no locks, no proc awareness.
//
// Cell state is split into two value types living in parallel slices that
// share an index i:
//
//	Lo[i]  — hot fields for pitch / space / field / clocks (pass-1 class)
//	Hi[i]  — planes, Jacobi scratch, forces (pass-2 class)
//
// A pass that only needs Lo streams a denser packing (more cells per
// cache line) without dragging Hi. Neither type contains pointers.
package core

// CellLo is the hot half of a free-cell site.
type CellLo struct {
	X, Y, Z float64
	Cr, Cr0 float64

	Es, Em, Ee float64
	Fa1, Fa2   float64
	Th2, Cbeta float64
	W1e, W2e   float64
}

// CellHi is the cold / scratch half (same index as CellLo).
type CellHi struct {
	N1x, N1y, N1z float64
	N2x, N2y, N2z float64

	EsNext  float64
	Fa1Next float64
	Fa2Next float64
	Fx, Fy, Fz float64
}

// Slot is a contact channel. Endpoints are indices into a Lo/Hi pair
// of slices (local to a World or a worker shard after ghost packing).
type Slot struct {
	I, J  int32
	Alive uint8
	D, Ux, Uy, Uz float64
	A             float64
	Swl           float64
	Sldd          float64
	WField        float64
}

// Law is a minimal constant table (subset of laws_V2g-shaped numbers).
type Law struct {
	C, Dt float64
	L     float64

	R0, Es0, EsFloor float64
	Cap              float64
	W1, W2, QDetune  float64

	SK, SDisp float64
	FieldJ    float64
	KRep      float64
	MobGeo    float64
	Cfac      float64

	Aref, Dref float64
}

// DefaultLaw returns V2g-ish defaults for the experimental subset.
func DefaultLaw() Law {
	return Law{
		C: 1, Dt: 0.02, L: 16,
		// R0=0.58 ⇒ on unit lattice spacing, d=1 < 2*R0 so lens A>0
		// (R0=0.5 sits exactly at contact edge and yields A=0 everywhere).
		R0: 0.58, Es0: 1, EsFloor: 0.05, Cap: 2.5,
		W1: 1.65, W2: 2.9, QDetune: 1.2,
		SK: 0.06, SDisp: 0.3,
		FieldJ: 1.8,
		KRep:   1.0, MobGeo: 1.0, Cfac: 1.15,
		Aref: 1.0, Dref: 1.0,
	}
}

// SyncEe writes Ee from Fa on the low half.
func (c *CellLo) SyncEe() { c.Ee = c.Fa1*c.Fa1 + c.Fa2*c.Fa2 }
