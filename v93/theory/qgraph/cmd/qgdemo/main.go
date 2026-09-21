// Command qgdemo prints the quantized-harmonic-graph tables the theory note
// quotes, so every number in docs/UNIFIED_LAW.md is reproducible:
//
//	go run ./cmd/qgdemo
package main

import (
	"fmt"
	"math"

	"scp/v93/theory/qgraph"
)

const a0 = 2 * math.Pi // eps = omega

func main() {
	fmt.Println("== 1. SHELLS ARE GRAPH SYMMETRY (Laplacian degeneracy classes) ==")
	for _, c := range []struct {
		name string
		g    *qgraph.Graph
	}{
		{"C6   ring", qgraph.NewCycle(6)},
		{"C10  ring", qgraph.NewCycle(10)},
		{"Q3   cube", qgraph.NewHypercube(3)},
		{"Q4   tesseract", qgraph.NewHypercube(4)},
		{"icosahedron", qgraph.NewIcosahedron()},
		{"P8   path (tree)", qgraph.NewPath(8)},
	} {
		sh := qgraph.Shells(qgraph.Spectrum(c.g), 1e-6)
		w := qgraph.NewWell(sh, math.Inf(1))
		caps := make([]int, len(sh))
		for i, s := range sh {
			caps[i] = s.Cap
		}
		fmt.Printf("  %-16s b1=%2d  caps=%v  magic=%v\n",
			c.name, c.g.CycleRank(), caps, w.MagicNumbers())
	}

	fmt.Println("\n== 2. CHARGE = SHELL DEFICIT (icosahedral well bound at shell 2, Z=9) ==")
	sh := qgraph.Shells(qgraph.Spectrum(qgraph.NewIcosahedron()), 1e-6)
	w := qgraph.NewWell(sh, sh[2].Omega+1e-9)
	fmt.Printf("  Z = %d      N:  ", w.Z())
	for n := 6; n <= 11; n++ {
		fmt.Printf("%2d ", n)
	}
	fmt.Printf("\n              Q:  ")
	for n := 6; n <= 11; n++ {
		fmt.Printf("%2d ", qgraph.NewObject(w, n).Charge())
	}
	fmt.Println("\n  (Q = 0 only at closure; both signs occur; N is conserved so Q is)")

	fmt.Println("\n== 3. RING SHELL NEUTRALITY (circulation of a full vs half shell) ==")
	for _, n := range []int{6, 10} {
		fmt.Printf("  C%-3d", n)
		for _, s := range qgraph.RingShells(n) {
			full := make([]int, len(s))
			for i := range full {
				full[i] = 1
			}
			c := qgraph.ShellCirculation(s, full)
			half := math.NaN()
			if len(s) == 2 {
				half = qgraph.ShellCirculation(s, []int{1, 0})
			}
			fmt.Printf(" [m=%+d cap=%d full=%.0e half=%+.3f]", s[0].M, len(s), math.Abs(c), half)
		}
		fmt.Println()
	}

	fmt.Println("\n== 4. DOORS ARE SELECTION RULES (ring C12) ==")
	for _, d := range []struct {
		name string
		p    qgraph.DoorProfile
	}{
		{"uniform  (neutral light)", qgraph.UniformDoor(12)},
		{"wound k=2 (organized field)", qgraph.WoundDoor(12, 2)},
		{"local     (point contact)", qgraph.LocalDoor(12, 0)},
	} {
		allowed := qgraph.AllowedTransfers(12, d.p, 1e-9)
		mean, variance := qgraph.CirculationTransfer(12, d.p)
		fmt.Printf("  %-28s allowed Dm=%v  mean=%+.3f  sd=%.3f\n",
			d.name, allowed, mean, math.Sqrt(variance))
	}

	fmt.Println("\n== 5. DUAL I/O FIXED POINT (k=0.05 cap=2.5 p=4, bath xb=0.62) ==")
	for _, g := range []float64{0, 0.5, 1, 2, 3} {
		x := qgraph.FixedPointLoad(0.05, 2.5, 0.62, 4, g)
		phi := qgraph.Throughput(0.05, 2.5, 0.62, 4, g)
		fmt.Printf("  coherence g=%.1f  x*=%.4f  x*/xb=%.4f  Phi*=%.5f  dJ=%+.4f\n",
			g, x, x/0.62, phi, qgraph.Restoring(0.05, 2.5, 0.62, 4, g))
	}
	fmt.Println("  g=0 (random-phase bath): x* = xb EXACTLY -> box-fill, never a particle")

	fmt.Println("\n== 6. STEFAN EXPONENT = SPECTRAL DIMENSION + 1 ==")
	for _, d := range []int{1, 2, 3} {
		L := map[int]int{1: 60, 2: 30, 3: 14}[d]
		g := qgraph.NewTorus(d, L)
		rw := qgraph.ReturnDimension(g, 128, []int{0, g.N / 3, 2 * g.N / 3})
		eig := qgraph.SpectralDimensionFromEigenvalues(qgraph.TorusSpectrum(d, L), 0.10)
		fmt.Printf("  %dD torus L=%-3d  d_s(walk)=%.3f  d_s(count)=%.3f  =>  p_rad = %.0f\n",
			d, L, rw, eig, qgraph.StefanExponent(float64(d)))
	}
	fmt.Printf("  Stefan constant at d_s=3: %.6f  (pi^4/15 = %.6f)\n",
		qgraph.StefanConstant(3), math.Pow(math.Pi, 4)/15)
	fmt.Println("  low-load (frozen stack) exponents p_eff = 4/(dm+1):")
	for dm := 1.0; dm <= 3; dm++ {
		fmt.Printf("    dm=%.0f -> p_eff=%.3f\n", dm, qgraph.EffectiveRadianceExponent(3, dm))
	}

	fmt.Println("\n== 7. MONODISPERSITY IS Delta/T (C10 well, mu=1.7) ==")
	wr := qgraph.NewWell(qgraph.Shells(qgraph.Spectrum(qgraph.NewCycle(10)), 1e-8), math.Inf(1))
	fmt.Printf("  magic numbers %v\n", wr.MagicNumbers())
	for _, T := range []float64{0.05, 0.3, 1.5} {
		p := qgraph.MassDistribution(wr, 10, 1.7, T, a0)
		n, pk := qgraph.PeakOccupancy(p)
		tag := "hierarchy"
		if wr.IsMagic(n) && pk > 0.5 {
			tag = "SPECIES"
		}
		fmt.Printf("  T=%-5.2f peak n=%-3d P=%.3f  %s\n", T, n, pk, tag)
	}
	fmt.Println("  residence time exp(Delta/2T)/2:")
	for _, r := range []float64{1, 4, 8, 16} {
		fmt.Printf("    Delta/T=%-3.0f tau=%.1f\n", r, qgraph.ResidenceTime(r, 1))
	}
}
