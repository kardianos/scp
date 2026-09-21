package qgraph

import "math"

// Nu is the per-mode occupation ceiling: the cap-refusal quantum.
//
// Nu = 1 is the programme's PAULI-0 seed made exact — a graph harmonic either
// carries its atom or does not. Nu > 1 recovers a boson-like ladder. Nothing
// in the law below imports fermion statistics; Nu is the cap, and the cap is
// already in laws_V3a.
const Nu = 1

// Shell is a degeneracy class of the Laplacian: all modes sharing one
// eigenvalue. Shells are NOT chosen — they are the eigenspaces the graph's own
// automorphisms force, which is why their capacities are integers.
type Shell struct {
	K      int     // shell index, ordered by energy
	Modes  []int   // indices into the mode slice
	Lambda float64 // shared eigenvalue
	Omega  float64 // shared pitch
	Cap    int     // capacity = len(Modes) * Nu
}

// Shells groups a sorted mode slice into degeneracy classes.
//
// tol is the absolute eigenvalue tolerance for calling two modes degenerate.
func Shells(modes []Mode, tol float64) []Shell {
	var out []Shell
	for i := 0; i < len(modes); {
		j := i
		for j < len(modes) && math.Abs(modes[j].Lambda-modes[i].Lambda) <= tol {
			j++
		}
		idx := make([]int, 0, j-i)
		for k := i; k < j; k++ {
			idx = append(idx, k)
		}
		out = append(out, Shell{
			K:      len(out),
			Modes:  idx,
			Lambda: modes[i].Lambda,
			Omega:  modes[i].Omega,
			Cap:    len(idx) * Nu,
		})
		i = j
	}
	return out
}

// Well is a bound harmonic stack: the shells of a graph, truncated at the
// binding threshold.
//
// The threshold is the programme's own Es depression depth: harmonics with
// pitch below it are bound to the object, harmonics above it belong to the
// exterior field. NBound is therefore geometry, not a dial.
type Well struct {
	Shells []Shell
	NBound int // number of shells that are bound (0..len(Shells))
}

// NewWell truncates a shell list at a pitch threshold.
func NewWell(shells []Shell, omegaBind float64) Well {
	nb := 0
	for _, s := range shells {
		if s.Omega <= omegaBind {
			nb++
		}
	}
	return Well{Shells: shells, NBound: nb}
}

// Z is the well's total bound capacity — the spectral invariant that plays the
// role a nuclear charge plays in chemistry, derived from the contact graph.
func (w Well) Z() int {
	z := 0
	for k := 0; k < w.NBound && k < len(w.Shells); k++ {
		z += w.Shells[k].Cap
	}
	return z
}

// MagicNumbers returns the closed-shell totals Z_j = sum_{k<=j} Cap[k] over
// the bound shells. These are the ONLY neutral occupancies, and therefore the
// only masses a species may have (UNIFIED_LAW §7).
func (w Well) MagicNumbers() []int {
	var out []int
	z := 0
	for k := 0; k < w.NBound && k < len(w.Shells); k++ {
		z += w.Shells[k].Cap
		out = append(out, z)
	}
	return out
}

// IsMagic reports whether n closes a bound shell exactly.
func (w Well) IsMagic(n int) bool {
	for _, m := range w.MagicNumbers() {
		if m == n {
			return true
		}
	}
	return false
}

// Aufbau fills n quanta into the shells from the bottom up and returns the
// per-shell occupation. Overflow beyond the bound stack lands in the first
// unbound shells (a captured-but-unbound quantum: the anion branch).
func (w Well) Aufbau(n int) []int {
	occ := make([]int, len(w.Shells))
	for k := range w.Shells {
		if n <= 0 {
			break
		}
		take := w.Shells[k].Cap
		if take > n {
			take = n
		}
		occ[k] = take
		n -= take
	}
	return occ
}

// Energy returns the total mode energy of an occupation, sum_k occ_k * eps_k,
// with eps the action atom A0*omega/(2*pi).
func (w Well) Energy(occ []int, a0 float64) float64 {
	e := 0.0
	for k, n := range occ {
		if k < len(w.Shells) {
			e += float64(n) * a0 * w.Shells[k].Omega / (2 * math.Pi)
		}
	}
	return e
}

// Gap returns the shell gap Delta_j = eps_(j+1) - eps_j in action-atom units.
//
// Delta_j / T is the ONLY source of monodispersity in this law: the residence
// time at a closed configuration is ~ exp(Delta_j/2T) (monodispersity.mac).
func (w Well) Gap(j int, a0 float64) float64 {
	if j < 0 || j+1 >= len(w.Shells) {
		return 0
	}
	return a0 * (w.Shells[j+1].Omega - w.Shells[j].Omega) / (2 * math.Pi)
}

// ValenceShell returns the index of the highest shell with any occupancy,
// and -1 for an empty configuration.
func ValenceShell(occ []int) int {
	v := -1
	for k, n := range occ {
		if n > 0 {
			v = k
		}
	}
	return v
}

// Incomplete reports whether the valence shell is partially filled — the
// structural condition the user's postulate names as the origin of charge.
func (w Well) Incomplete(occ []int) bool {
	v := ValenceShell(occ)
	if v < 0 || v >= len(w.Shells) {
		return false
	}
	return occ[v] > 0 && occ[v] < w.Shells[v].Cap
}
