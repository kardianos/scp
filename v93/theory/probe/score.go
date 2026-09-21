// optional: go run ./theory/probe - but scoring is in main.
// This file holds shared sustain classification docs only if needed.
package main

// Sustain classes used in classifySustain:
//
//	HOLD       — Em substantial, last-half slope small, not melting
//	GROW       — still accumulating (open I > O)
//	DRAIN      — systematic loss
//	MELT       — Em collapsed
//	empty      — never formed
//	alive      — Em present but slope not plateau-clean
//	quasi-hold — Em ok but envelope drifted
