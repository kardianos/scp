package qgraph

import "math"

// ---------------------------------------------------------------------------
// Charge as an occupancy defect of the graph's own harmonic shells.
//
// Two integers describe an object:
//
//	Z  = the well's bound capacity        — a SPECTRAL invariant of the graph
//	N  = the number of quanta it carries  — a CONSERVED count
//
//	Q  = Z - N
//
// Q is not a renamed energy: Z counts eigenspace dimensions (a topological /
// symmetry fact about the contact graph) while N counts atoms. Q > 0 is an
// unclosed shell (a deficit that will pull); Q < 0 is a quantum held above
// closure (an excess that will push). Q = 0 exactly at the magic numbers.
//
// Because emission, absorption and internal promotion all leave N invariant
// (they move a quantum between modes and pay the difference in field quanta),
// Q changes ONLY by transfer between objects or by graph surgery. That is the
// conservation law, and it is exact — see ChargeIsConserved.
// ---------------------------------------------------------------------------

// Object is a bound harmonic stack with an occupation.
type Object struct {
	Well Well
	Occ  []int // per-shell occupancy
}

// NewObject fills n quanta into a well by aufbau.
func NewObject(w Well, n int) *Object { return &Object{Well: w, Occ: w.Aufbau(n)} }

// N is the total quantum count.
func (o *Object) N() int {
	n := 0
	for _, v := range o.Occ {
		n += v
	}
	return n
}

// Charge returns Q = Z - N.
func (o *Object) Charge() int { return o.Well.Z() - o.N() }

// IsNeutral reports Q == 0, equivalently: every bound shell is closed.
func (o *Object) IsNeutral() bool { return o.Charge() == 0 }

// Energy is the total action-atom energy of the occupation.
func (o *Object) Energy(a0 float64) float64 { return o.Well.Energy(o.Occ, a0) }

// Promote moves one quantum from shell a to shell b, respecting capacity.
// It returns the emitted field quantum's pitch (positive when the object
// falls and radiates, negative when it climbs and must absorb), and false if
// the move is blocked.
//
// This single operation IS both emission and absorption: the sign of the
// returned pitch selects the direction. There is no second law for intake.
func (o *Object) Promote(a, b int) (omega float64, ok bool) {
	if a < 0 || b < 0 || a >= len(o.Occ) || b >= len(o.Occ) || a == b {
		return 0, false
	}
	if o.Occ[a] == 0 || o.Occ[b] >= o.Well.Shells[b].Cap {
		return 0, false // empty source, or cap refusal at the destination
	}
	o.Occ[a]--
	o.Occ[b]++
	return o.Well.Shells[a].Omega - o.Well.Shells[b].Omega, true
}

// Transfer moves one quantum from src to dst (ionization / capture). This is
// the only single-object operation that changes charge, and it changes the two
// charges oppositely — total charge is untouched.
func Transfer(src, dst *Object) bool {
	v := ValenceShell(src.Occ)
	if v < 0 {
		return false
	}
	for k := range dst.Occ {
		if dst.Occ[k] < dst.Well.Shells[k].Cap {
			src.Occ[v]--
			dst.Occ[k]++
			return true
		}
	}
	return false
}

// TotalCharge sums Q over a set of objects.
func TotalCharge(objs ...*Object) int {
	q := 0
	for _, o := range objs {
		q += o.Charge()
	}
	return q
}

// ---------------------------------------------------------------------------
// Circulation — the dynamical face of the same defect.
//
// On a ring C_N the harmonics are chi_m(j) = exp(2*pi*i*m*j/N)/sqrt(N) with
// winding label m and link current
//
//	circ_m ~ sin(2*pi*m/N).
//
// A full shell {+m,-m} cancels: sum circ = 0 (shells.mac proves this
// identically). Both singlets (m=0 and, for even N, m=N/2) are real standing
// modes with circ = 0. So NET CIRCULATION IS CARRIED EXACTLY BY THE UNCLOSED
// SHELL — the occupancy defect and the residual current are the same object.
// ---------------------------------------------------------------------------

// RingMode is an analytic C_N harmonic carrying its winding label.
type RingMode struct {
	M      int     // winding label in (-N/2, N/2]
	Lambda float64 // 2(1-cos(2*pi*m/N))
	Omega  float64
	Circ   float64 // sin(2*pi*m/N): the physical link current of the mode
}

// RingModes returns the N analytic harmonics of C_N, sorted by eigenvalue.
func RingModes(n int) []RingMode {
	out := make([]RingMode, 0, n)
	for m := 0; m < n; m++ {
		mm := m
		if 2*m > n {
			mm = m - n
		}
		th := 2 * math.Pi * float64(m) / float64(n)
		lam := 2 * (1 - math.Cos(th))
		out = append(out, RingMode{M: mm, Lambda: lam, Omega: math.Sqrt(lam), Circ: math.Sin(th)})
	}
	// stable sort by lambda then |m|
	for i := 1; i < len(out); i++ {
		for j := i; j > 0 && (out[j].Lambda < out[j-1].Lambda-1e-12); j-- {
			out[j], out[j-1] = out[j-1], out[j]
		}
	}
	return out
}

// RingShells groups analytic ring modes into +-m doublets and singlets.
func RingShells(n int) [][]RingMode {
	modes := RingModes(n)
	var out [][]RingMode
	for i := 0; i < len(modes); {
		j := i
		for j < len(modes) && math.Abs(modes[j].Lambda-modes[i].Lambda) < 1e-9 {
			j++
		}
		out = append(out, modes[i:j])
		i = j
	}
	return out
}

// ShellCirculation sums the link current of the OCCUPIED members of a shell.
// occ[i] is the occupancy of member i of the shell.
func ShellCirculation(shell []RingMode, occ []int) float64 {
	c := 0.0
	for i, m := range shell {
		if i < len(occ) {
			c += float64(occ[i]) * m.Circ
		}
	}
	return c
}

// ShellWinding sums the integer winding labels of the occupied members.
// This is the quantized charge label; ShellCirculation is its current.
func ShellWinding(shell []RingMode, occ []int) int {
	w := 0
	for i, m := range shell {
		if i < len(occ) {
			w += occ[i] * m.M
		}
	}
	return w
}
