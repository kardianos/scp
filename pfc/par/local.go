package par

import "scp/pfc/core"

// BuildLocal is the owned-slab constructor (same as FromWorld).
// Kept so older call sites / bench labels still work.
func BuildLocal(w *core.World, nProc int) *Cluster {
	return FromWorld(w, nProc)
}
