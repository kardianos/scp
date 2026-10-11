// unsigned.go — v103 lane c1, post-hoc (not registered): is there an
// unsigned long-range void-void part beside the signed Coulomb part?
//
// From the four-void (E4) runs, the pair potential of each type (same /
// opposite sublattice) is fitted jointly on shells r >= 2 as
//
//	V = alpha_type + beta_s * s * G(r) + beta_u * G(r),   s = +1 same, -1 opp
//
// beta_s is the signed (charge) coupling, beta_u the unsigned one. Four
// charges in a periodic box screen one another, so a pure two-sign Coulomb
// gas can show a nonzero beta_u by itself. The null is that gas: two A and
// two B charges on the same lattice, weight exp(-kappa * sum q_i q_j G),
// kappa = -beta from the two-void (E2) fit, analysed identically.
//
//	go run unsigned.go            (from cas/unsigned, reads ../out)
package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"regexp"
	"sort"
	"strconv"
	"strings"
)

var rules = []string{"none", "a", "b", "c", "ab", "ac", "bc", "abc"}

type Bin struct {
	R     float64
	Opp   bool
	Ideal float64
	Cnt   []float64
	G0    float64
}

type Run struct {
	Lat         string
	N, Dim, Bip int
	RB          bool
	D           int
	Bins        []*Bin
}

func atof(s string) float64 { v, _ := strconv.ParseFloat(s, 64); return v }

var hdr = regexp.MustCompile(`lat=(\S+) .* N=(\d+) .* dim=(\d+) bipartite=(\d+) rules=\S+ \(a=\d b=(\d) c=\d\) CAP=\d q=\d D=(\d+)`)

func parse(path string) *Run {
	f, err := os.Open(path)
	if err != nil {
		return nil
	}
	defer f.Close()
	r := &Run{}
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 1<<20), 1<<20)
	done := false
	for sc.Scan() {
		l := sc.Text()
		switch {
		case strings.HasPrefix(l, "# chain"):
			m := hdr.FindStringSubmatch(l)
			if m == nil {
				return nil
			}
			r.Lat = m[1]
			r.N, _ = strconv.Atoi(m[2])
			r.Dim, _ = strconv.Atoi(m[3])
			r.Bip, _ = strconv.Atoi(m[4])
			r.RB = m[5] == "1"
			r.D, _ = strconv.Atoi(m[6])
		case strings.HasPrefix(l, "SHELL"):
			fs := strings.Fields(l)
			b := &Bin{}
			b.R = atof(strings.TrimPrefix(fs[2], "r="))
			b.Opp = fs[3] == "type=opp"
			b.Ideal = atof(strings.TrimPrefix(fs[4], "ideal="))
			i := 6
			for ; fs[i] != "G"; i++ {
				b.Cnt = append(b.Cnt, atof(fs[i]))
			}
			b.G0 = atof(fs[i+1])
			r.Bins = append(r.Bins, b)
		case strings.HasPrefix(l, "CHECK"):
			done = true
		}
	}
	if !done {
		return nil
	}
	return r
}

func meanSE(v []float64) (float64, float64) {
	n := float64(len(v))
	m := 0.0
	for _, x := range v {
		m += x
	}
	m /= n
	s := 0.0
	for _, x := range v {
		s += (x - m) * (x - m)
	}
	return m, math.Sqrt(s / (n - 1) / n)
}

type Pt struct {
	R, V, SE, G float64
	Opp         bool
}

// pair potential per type, normalised within the type, per seed
func potential(bins []*Bin, opp bool, rmin float64) []Pt {
	var bs []*Bin
	for _, b := range bins {
		if b.Opp == opp && b.Ideal > 0 {
			bs = append(bs, b)
		}
	}
	if len(bs) == 0 {
		return nil
	}
	ns := len(bs[0].Cnt)
	idealTot := 0.0
	tot := make([]float64, ns)
	for _, b := range bs {
		idealTot += b.Ideal
		for s := 0; s < ns; s++ {
			tot[s] += b.Cnt[s]
		}
	}
	var out []Pt
	for _, b := range bs {
		if b.R < rmin-1e-9 {
			continue
		}
		vs := make([]float64, 0, ns)
		ok := true
		for s := 0; s < ns; s++ {
			if b.Cnt[s] <= 0 || tot[s] <= 0 {
				ok = false
				break
			}
			vs = append(vs, -math.Log((b.Cnt[s]/tot[s])/(b.Ideal/idealTot)))
		}
		if !ok {
			continue
		}
		m, se := meanSE(vs)
		if se <= 0 {
			continue
		}
		out = append(out, Pt{R: b.R, V: m, SE: se, G: b.G0, Opp: opp})
	}
	sort.Slice(out, func(i, j int) bool { return out[i].R < out[j].R })
	return out
}

// weighted linear least squares; returns params, their se, chi2
func lsq(rows [][]float64, y, se []float64) ([]float64, []float64, float64) {
	p := len(rows[0])
	A := make([][]float64, p)
	for i := range A {
		A[i] = make([]float64, p+1)
	}
	for k := range rows {
		w := 1 / (se[k] * se[k])
		for i := 0; i < p; i++ {
			for j := 0; j < p; j++ {
				A[i][j] += w * rows[k][i] * rows[k][j]
			}
			A[i][p] += w * rows[k][i] * y[k]
		}
	}
	// invert the normal matrix (Gauss-Jordan) for parameters and covariance
	n := p
	M := make([][]float64, n)
	for i := range M {
		M[i] = make([]float64, 2*n+1)
		copy(M[i], A[i][:n])
		M[i][n+i] = 1
		M[i][2*n] = A[i][n]
	}
	for c := 0; c < n; c++ {
		piv := c
		for r := c + 1; r < n; r++ {
			if math.Abs(M[r][c]) > math.Abs(M[piv][c]) {
				piv = r
			}
		}
		M[c], M[piv] = M[piv], M[c]
		d := M[c][c]
		for j := range M[c] {
			M[c][j] /= d
		}
		for r := 0; r < n; r++ {
			if r != c {
				f := M[r][c]
				for j := range M[r] {
					M[r][j] -= f * M[c][j]
				}
			}
		}
	}
	par := make([]float64, n)
	pse := make([]float64, n)
	for i := 0; i < n; i++ {
		par[i] = M[i][2*n]
		pse[i] = math.Sqrt(M[i][n+i])
	}
	chi := 0.0
	for k := range rows {
		f := 0.0
		for i := 0; i < p; i++ {
			f += par[i] * rows[k][i]
		}
		chi += (y[k] - f) * (y[k] - f) / (se[k] * se[k])
	}
	return par, pse, chi
}

type Dec struct {
	Bs, BsSE, Bu, BuSE, Chi float64
	Dof                     int
	Two                     bool
}

// joint fit of both types (or one type when only one exists)
func decompose(bins []*Bin, rmin float64) Dec {
	ps := potential(bins, false, rmin)
	po := potential(bins, true, rmin)
	var rows [][]float64
	var y, se []float64
	two := len(po) >= 3 && len(ps) >= 3
	if two {
		for _, p := range ps {
			rows = append(rows, []float64{1, 0, p.G, p.G})
			y, se = append(y, p.V), append(se, p.SE)
		}
		for _, p := range po {
			rows = append(rows, []float64{0, 1, -p.G, p.G})
			y, se = append(y, p.V), append(se, p.SE)
		}
		par, pse, chi := lsq(rows, y, se)
		return Dec{par[2], pse[2], par[3], pse[3], chi, len(rows) - 4, true}
	}
	all := append(ps, po...)
	if len(all) < 3 {
		return Dec{}
	}
	for _, p := range all {
		rows = append(rows, []float64{1, p.G})
		y, se = append(y, p.V), append(se, p.SE)
	}
	par, pse, chi := lsq(rows, y, se)
	return Dec{0, 0, par[1], pse[1], chi, len(rows) - 2, false}
}

// E2 coupling: V = alpha + beta G on r >= 2 (as fit.go)
func betaE2(r *Run) (float64, float64) {
	pts := potential(r.Bins, true, 2)
	var rows [][]float64
	var y, se []float64
	for _, p := range pts {
		rows = append(rows, []float64{1, p.G})
		y, se = append(y, p.V), append(se, p.SE)
	}
	par, pse, _ := lsq(rows, y, se)
	return par[1], pse[1]
}

// ---------------- null: four-charge lattice Coulomb gas ----------------

type Lat struct {
	L, Dim, N int
	G         []float64 // by displacement index
	shellOf   []int     // displacement index -> shell
	shr       []float64
}

func newLat(L, dim int) *Lat {
	lt := &Lat{L: L, Dim: dim}
	N := 1
	for d := 0; d < dim; d++ {
		N *= L
	}
	lt.N = N
	coords := func(i int) []int {
		c := make([]int, dim)
		for d := 0; d < dim; d++ {
			c[d] = i % L
			i /= L
		}
		return c
	}
	// G(r) = (1/N) sum_{k != 0} cos(k.r) / sum_mu 2(1 - cos k_mu): the zero-mean
	// solution of L G = delta - 1/N, as chain.c builds it by CG
	lt.G = make([]float64, N)
	for r := 0; r < N; r++ {
		rc := coords(r)
		s := 0.0
		for k := 1; k < N; k++ {
			kc := coords(k)
			den, ph := 0.0, 0
			for d := 0; d < dim; d++ {
				den += 2 * (1 - math.Cos(2*math.Pi*float64(kc[d])/float64(L)))
				ph += kc[d] * rc[d]
			}
			s += math.Cos(2*math.Pi*float64(ph%L)/float64(L)) / den
		}
		lt.G[r] = s / float64(N)
	}
	// shells by minimum-image distance
	lt.shellOf = make([]int, N)
	var rs []float64
	dist := make([]float64, N)
	for r := 0; r < N; r++ {
		rc := coords(r)
		d2 := 0.0
		for d := 0; d < dim; d++ {
			x := rc[d]
			if x > L/2 {
				x = L - x
			}
			d2 += float64(x * x)
		}
		dist[r] = math.Sqrt(d2)
		found := false
		for _, v := range rs {
			if math.Abs(v-dist[r]) < 1e-6 {
				found = true
				break
			}
		}
		if !found {
			rs = append(rs, dist[r])
		}
	}
	sort.Float64s(rs)
	lt.shr = rs
	for r := 0; r < N; r++ {
		for s, v := range rs {
			if math.Abs(v-dist[r]) < 1e-6 {
				lt.shellOf[r] = s
				break
			}
		}
	}
	return lt
}

func (lt *Lat) disp(i, j int) int {
	r, mul := 0, 1
	for d := 0; d < lt.Dim; d++ {
		a, b := i%lt.L, j%lt.L
		i /= lt.L
		j /= lt.L
		r += ((b - a + lt.L) % lt.L) * mul
		mul *= lt.L
	}
	return r
}

func (lt *Lat) parity(i int) int {
	s := 0
	for d := 0; d < lt.Dim; d++ {
		s += i % lt.L
		i /= lt.L
	}
	return s & 1
}

type rng struct{ s uint64 }

func (r *rng) next() uint64 {
	r.s += 0x9e3779b97f4a7c15
	z := r.s
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb
	return z ^ (z >> 31)
}
func (r *rng) f() float64 { return float64(r.next()>>11) / (1 << 53) }

// four charges (+,+ on parity 0; -,- on parity 1), kappa > 0 attracts unlike
func null(lt *Lat, kappa float64, seeds int, sweeps int) []*Bin {
	NS := len(lt.shr)
	sites := [2][]int{}
	for i := 0; i < lt.N; i++ {
		sites[lt.parity(i)] = append(sites[lt.parity(i)], i)
	}
	q := []float64{1, 1, -1, -1}
	cls := []int{0, 0, 1, 1}
	cnt := [2][][]float64{}
	ideal := [2][]float64{make([]float64, NS), make([]float64, NS)}
	for t := 0; t < 2; t++ {
		cnt[t] = make([][]float64, NS)
		for s := range cnt[t] {
			cnt[t][s] = make([]float64, seeds)
		}
	}
	// ideal: pairs (i fixed in class, j over its class / the other class)
	for _, i := range []int{sites[0][0]} {
		for j := 0; j < lt.N; j++ {
			t := lt.parity(j) != lt.parity(i)
			ti := 0
			if t {
				ti = 1
			}
			ideal[ti][lt.shellOf[lt.disp(i, j)]] += 1
		}
	}
	for sd := 0; sd < seeds; sd++ {
		g := &rng{uint64(0x51ed270b) + 1000003*uint64(sd+1)}
		pos := make([]int, 4)
		for k := 0; k < 4; k++ {
			c := sites[cls[k]]
			pos[k] = c[int(g.next()%uint64(len(c)))]
		}
		energy := func(k, p int) float64 {
			e := 0.0
			for m := 0; m < 4; m++ {
				if m != k {
					e += kappa * q[k] * q[m] * lt.G[lt.disp(p, pos[m])]
				}
			}
			return e
		}
		for sw := 0; sw < sweeps+sweeps/50; sw++ {
			for k := 0; k < 4; k++ {
				c := sites[cls[k]]
				np := c[int(g.next()%uint64(len(c)))]
				dE := energy(k, np) - energy(k, pos[k])
				if dE <= 0 || g.f() < math.Exp(-dE) {
					pos[k] = np
				}
			}
			if sw < sweeps/50 {
				continue
			}
			for a := 0; a < 4; a++ {
				for b := a + 1; b < 4; b++ {
					t := 0
					if cls[a] != cls[b] {
						t = 1
					}
					cnt[t][lt.shellOf[lt.disp(pos[a], pos[b])]][sd] += 1
				}
			}
		}
	}
	var bins []*Bin
	for s := 0; s < NS; s++ {
		for t := 0; t < 2; t++ {
			if ideal[t][s] == 0 {
				continue
			}
			// shell-averaged G over the displacements in this shell and type
			gs, w := 0.0, 0.0
			for r := 0; r < lt.N; r++ {
				if lt.shellOf[r] == s && (lt.parity(r) == 1) == (t == 1) {
					gs += lt.G[r]
					w++
				}
			}
			bins = append(bins, &Bin{R: lt.shr[s], Opp: t == 1, Ideal: ideal[t][s], Cnt: cnt[t][s], G0: gs / w})
		}
	}
	return bins
}

func main() {
	fmt.Printf("# c1 post-hoc: signed and unsigned parts of the four-void potential (shells r >= 2)\n\n")
	fmt.Printf("| rules | d | beta_s (signed) | beta_u (unsigned) | chi2/dof | beta_E2 |\n|---|---|---|---|---|---|\n")
	type key struct{ r, d string }
	data := map[key]Dec{}
	runs := map[key]*Run{}
	for _, d := range []string{"2d", "3d"} {
		for _, r := range rules {
			run := parse(fmt.Sprintf("../out/c1_%s_%s_E4.out", d, r))
			if run == nil {
				fmt.Printf("| %s | %s | (missing) | | | |\n", r, d)
				continue
			}
			dc := decompose(run.Bins, 2)
			data[key{r, d}] = dc
			runs[key{r, d}] = run
			e2 := "—"
			if r2 := parse(fmt.Sprintf("../out/c1_%s_%s_E2.out", d, r)); r2 != nil && r2.Bip == 1 && r2.RB {
				b, bse := betaE2(r2)
				e2 = fmt.Sprintf("%.3f ± %.3f", b, bse)
			}
			bs := "—"
			if dc.Two {
				bs = fmt.Sprintf("%+.3f ± %.3f", dc.Bs, dc.BsSE)
			}
			fmt.Printf("| %s | %s | %s | %+.3f ± %.3f | %.2f | %s |\n", r, d, bs, dc.Bu, dc.BuSE, dc.Chi/float64(dc.Dof), e2)
		}
	}
	fmt.Printf("\n## null: four-charge lattice Coulomb gas, kappa = -beta_E2 (no unsigned part by construction)\n\n")
	fmt.Printf("| rules | d | kappa | null beta_s | null beta_u | data beta_u - null beta_u | z |\n|---|---|---|---|---|---|---|\n")
	for _, d := range []string{"2d", "3d"} {
		var lt *Lat
		if d == "2d" {
			lt = newLat(32, 2)
		} else {
			lt = newLat(12, 3)
		}
		for _, r := range []string{"bc", "abc"} {
			r2 := parse(fmt.Sprintf("../out/c1_%s_%s_E2.out", d, r))
			if r2 == nil {
				continue
			}
			b, _ := betaE2(r2)
			kappa := -b
			bins := null(lt, kappa, 8, 2000000)
			// the null's shell-averaged G must be the tool's
			gmax := 0.0
			for _, nb := range bins {
				for _, db := range r2.Bins {
					if math.Abs(nb.R-db.R) < 1e-6 && nb.Opp == db.Opp {
						gmax = math.Max(gmax, math.Abs(nb.G0-db.G0))
					}
				}
			}
			fmt.Fprintf(os.Stderr, "%s %s: max |G_null - G_tool| over shells = %.2e\n", r, d, gmax)
			nd := decompose(bins, 2)
			dd := data[key{r, d}]
			diff := dd.Bu - nd.Bu
			se := math.Hypot(dd.BuSE, nd.BuSE)
			fmt.Printf("| %s | %s | %.3f | %+.3f ± %.3f | %+.3f ± %.3f | %+.3f ± %.3f | %.1f |\n", r, d, kappa, nd.Bs, nd.BsSE, nd.Bu, nd.BuSE, diff, se, diff/se)
		}
	}
}
