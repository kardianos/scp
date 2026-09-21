package parcore

// PublishFaces returns a KindParallel Phase that copies this shard's face cells
// (the publish set) into the global publish buffer. Pair it with PullGhosts
// across an inter-phase barrier to refresh remote halos.
func PublishFaces[C any, S any]() Phase[C, S] {
	return Phase[C, S]{
		Name: "publish-faces",
		Kind: KindParallel,
		Run:  func(v *View[C, S]) { v.Publish() },
	}
}

// PullGhosts returns a KindParallel Phase that refreshes this shard's ghost
// cells from the global publish buffer. Must run AFTER every shard's
// PublishFaces; the inter-phase barrier enforces it.
func PullGhosts[C any, S any]() Phase[C, S] {
	return Phase[C, S]{
		Name: "pull-ghosts",
		Kind: KindParallel,
		Run:  func(v *View[C, S]) { v.Pull() },
	}
}

// BoundaryHops returns a KindBoundary Phase that runs the driver's dual-
// endpoint hop function once on the caller thread over every cross-shard edge
// in the global publish buffer. hop receives the buffer and one edge at a time
// (gi, gj global ids + the per-step tau). Use the SAME hop the interior
// colored rounds use (typically built from Givens2) so the ledger is
// conserved identically on interior and boundary edges.
func BoundaryHops[C any, S any](hop func(pub []C, gi, gj int, tau float64)) Phase[C, S] {
	return Phase[C, S]{
		Name: "boundary-hops",
		Kind: KindBoundary,
		Bnd: func(bv *BoundaryView[C]) {
			for _, e := range bv.Edges {
				hop(bv.Pub, int(e.Gi), int(e.Gj), e.Tau)
			}
		},
	}
}
