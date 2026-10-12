// tab.go — v103 lane c, step c3: the registered tables and the ending
// (../PREREG.md §7). Reads out/c3_*.out, written by c3pair.c.
//
//	go run tab.go > out/tab.out
package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

type rec map[string]string

func (r rec) f(k string) float64 {
	x, err := strconv.ParseFloat(r[k], 64)
	if err != nil {
		return math.NaN()
	}
	return x
}

func (r rec) i(k string) int {
	x, err := strconv.Atoi(r[k])
	if err != nil {
		return -1
	}
	return x
}

// one process: one (set, L, placement k)
type proc struct {
	set     string
	L, k    int
	done    bool
	med     rec
	made    map[string]rec // kind/R
	quench  map[string]rec
	windows map[string]rec
	pairs   map[string]rec
	checks  []rec
	fails   []rec
	skips   []rec
}

func key(kind string, R int) string { return fmt.Sprintf("%s/%d", kind, R) }

func parse(path string) (*proc, error) {
	fh, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer fh.Close()
	p := &proc{made: map[string]rec{}, quench: map[string]rec{}, windows: map[string]rec{}, pairs: map[string]rec{}}
	sc := bufio.NewScanner(fh)
	sc.Buffer(make([]byte, 1<<20), 1<<20)
	for sc.Scan() {
		fs := strings.Fields(sc.Text())
		if len(fs) == 0 {
			continue
		}
		r := rec{}
		for _, x := range fs[1:] {
			if kv := strings.SplitN(x, "=", 2); len(kv) == 2 {
				r[kv[0]] = kv[1]
			}
		}
		switch fs[0] {
		case "RUN":
			p.set, p.L, p.k = r["set"], r.i("L"), r.i("k")
		case "MED":
			p.med = r
		case "MADE":
			p.made[key(r["kind"], r.i("R"))] = r
		case "QUENCH":
			p.quench[key(r["kind"], r.i("R"))] = r
		case "WINDOW":
			p.windows[key(r["kind"], r.i("R"))] = r
		case "PAIR":
			p.pairs[key(r["kind"], r.i("R"))] = r
		case "CHECK":
			p.checks = append(p.checks, r)
		case "FAIL":
			p.fails = append(p.fails, r)
		case "SKIP":
			p.skips = append(p.skips, r)
		case "DONE":
			p.done = true
		}
	}
	return p, sc.Err()
}

var (
	sets  = []string{"none", "a", "b", "c", "ab", "ac", "bc", "abc"}
	Ls    = []int{8, 12, 16}
	barLs = []int{12, 16}
	Rs    = []int{1, 2, 3, 4, 6, 8}
	kinds = []string{"GG", "XX", "GX"}
	data  = map[string]map[int][]*proc{} // set -> L -> processes
)

func clean(x []float64) []float64 {
	y := []float64{}
	for _, v := range x {
		if !math.IsNaN(v) {
			y = append(y, v)
		}
	}
	return y
}

func median(x []float64) float64 {
	y := clean(x)
	if len(y) == 0 {
		return math.NaN()
	}
	sort.Float64s(y)
	n := len(y)
	if n%2 == 1 {
		return y[n/2]
	}
	return 0.5 * (y[n/2-1] + y[n/2])
}

func mean(x []float64) float64 {
	y := clean(x)
	if len(y) == 0 {
		return math.NaN()
	}
	s := 0.0
	for _, v := range y {
		s += v
	}
	return s / float64(len(y))
}

// a field of the PAIR lines of one kind and separation, over the placements
func collect(s string, L int, kind string, R int, field string) []float64 {
	out := []float64{}
	for _, p := range data[s][L] {
		if r, ok := p.pairs[key(kind, R)]; ok {
			out = append(out, r.f(field))
		}
	}
	return out
}

// a field of the QUENCH or WINDOW line of one single (G1, X1 at R = 0)
func single(s string, L int, kind string, which string, field string) []float64 {
	out := []float64{}
	for _, p := range data[s][L] {
		m := p.quench
		if which == "window" {
			m = p.windows
		} else if which == "made" {
			m = p.made
		}
		if r, ok := m[key(kind, 0)]; ok {
			out = append(out, r.f(field))
		}
	}
	return out
}

func healed(s string, L int, kind string, R int, field string) (h, n int) {
	for _, p := range data[s][L] {
		if r, ok := p.pairs[key(kind, R)]; ok {
			n++
			if r[field] == "1" {
				h++
			}
		}
	}
	return
}

func yes(b bool) string {
	if b {
		return "✔"
	}
	return "✘"
}

func num(x float64, f string) string {
	if math.IsNaN(x) {
		return "(missing)"
	}
	return fmt.Sprintf(f, x)
}

// PREREG §7: the class of one kind at contact (R = 1) at one size
func classOf(s, kind string, L int) string {
	vals := clean(collect(s, L, kind, 1, "Eint_w"))
	if len(vals) < 3 {
		return "unreadable"
	}
	h, n := healed(s, L, kind, 1, "healed_end")
	if 4*h >= 3*n {
		return "annihilate"
	}
	m := median(vals)
	switch {
	case m <= -2:
		return "attract"
	case m >= 2:
		return "repel"
	}
	return "none"
}

func letter(c map[string]string) string {
	att := func(x string) bool { return x == "attract" || x == "annihilate" }
	for _, x := range c {
		if x == "unreadable" || x == "size-dependent" {
			return "M"
		}
	}
	gg, xx, gx := c["GG"], c["XX"], c["GX"]
	switch {
	case att(gg) && att(xx) && att(gx):
		return "U"
	case att(gg) && att(xx):
		return "L"
	case att(gx) && !att(gg) && !att(xx):
		return "Q"
	case gg == "none" && xx == "none" && gx == "none":
		return "N"
	}
	return "M"
}

// the largest R with |median E_int| >= 2t for any kind; and whether every kind is below 2t at R = 8 (L = 16 only)
func reach(s string, L int) (int, bool) {
	RE, far := 0, true
	for _, kind := range kinds {
		for _, R := range Rs {
			if 2*R > L {
				continue
			}
			m := median(collect(s, L, kind, R, "Eint_w"))
			if !math.IsNaN(m) && math.Abs(m) >= 2 {
				if R > RE {
					RE = R
				}
				if R == 8 {
					far = false
				}
			}
		}
	}
	return RE, far
}

func main() {
	files, _ := filepath.Glob("out/c3_*.out")
	nfile, ndone, nbad, nfail, nskip, nchk := 0, 0, 0, 0, 0, 0
	secs := 0.0
	for _, fn := range files {
		p, err := parse(fn)
		if err != nil || p.set == "" {
			fmt.Printf("unreadable: %s\n", fn)
			continue
		}
		nfile++
		if p.done {
			ndone++
		}
		for _, c := range p.checks {
			nchk++
			if c["ok"] != "1" {
				nbad++
			}
		}
		for _, w := range p.windows {
			if v := w.f("secs"); !math.IsNaN(v) {
				secs += v
			}
		}
		nfail += len(p.fails)
		nskip += len(p.skips)
		if data[p.set] == nil {
			data[p.set] = map[int][]*proc{}
		}
		data[p.set][p.L] = append(data[p.set][p.L], p)
	}
	fmt.Printf("# v103 c3 — registered tables (PREREG §7)\n\n")
	fmt.Printf("## 0. Ledger\n\n- processes read: %d of 96; finished (DONE): %d\n- CHECK lines: %d, failed: %d\n- FAIL lines: %d; SKIP lines: %d\n- CPU in media: %.1f h\n", nfile, ndone, nchk, nbad, nfail, nskip, secs/3600)
	for _, s := range sets {
		for _, L := range Ls {
			for _, p := range data[s][L] {
				for _, f := range p.fails {
					fmt.Printf("- FAIL set=%s L=%d k=%d kind=%s R=%s what=%s\n", s, L, p.k, f["kind"], f["R"], f["what"])
				}
				for _, f := range p.skips {
					fmt.Printf("- SKIP set=%s L=%d k=%d R=%s what=%s\n", s, L, p.k, f["R"], f["what"])
				}
			}
		}
	}

	// 1. the medium, and which sets are scored
	scored := map[string]bool{}
	fmt.Printf("\n## 1. The medium s0 (every placement), and the scored sets\n\n| set | L=8 | L=12 | L=16 | scored |\n|---|---|---|---|---|\n")
	for _, s := range sets {
		cells := []string{}
		ok := true
		for _, L := range Ls {
			ps := data[s][L]
			if len(ps) == 0 {
				cells = append(cells, "(missing)")
				if L != 8 {
					ok = false
				}
				continue
			}
			cryst, nm := true, 0.0
			for _, p := range ps {
				if p.med == nil || p.med.i("changed_vs_flat") != 0 || p.med.i("nmis") != 0 {
					cryst = false
				}
				if p.med != nil {
					nm += p.med.f("nmis")
				}
			}
			if cryst {
				cells = append(cells, "crystal (unchanged)")
			} else {
				cells = append(cells, fmt.Sprintf("re-formed (mismatch %.0f)", nm/float64(len(ps))))
				if L != 8 {
					ok = false
				}
			}
		}
		scored[s] = ok
		fmt.Printf("| %s | %s | %s |\n", s, strings.Join(cells, " | "), yes(ok))
	}

	// 2. one defect alone
	fmt.Printf("\n## 2. One defect alone (median over placements): energy after making → after the quench, and its clusters\n\n| set | L | gap E (t) | gap quench flips | gap clusters | excess E (t) | excess quench flips | excess clusters |\n|---|---|---|---|---|---|---|---|\n")
	for _, s := range sets {
		for _, L := range Ls {
			if len(data[s][L]) == 0 {
				continue
			}
			fmt.Printf("| %s | %d | %s → %s | %s | %s | %s → %s | %s | %s |\n", s, L,
				num(median(single(s, L, "G1", "made", "E_after")), "%.0f"), num(median(single(s, L, "G1", "quench", "E")), "%.0f"),
				num(median(single(s, L, "G1", "made", "quench_steps")), "%.0f"), num(median(single(s, L, "G1", "quench", "ncl")), "%.0f"),
				num(median(single(s, L, "X1", "made", "E_after")), "%.0f"), num(median(single(s, L, "X1", "quench", "E")), "%.0f"),
				num(median(single(s, L, "X1", "made", "quench_steps")), "%.0f"), num(median(single(s, L, "X1", "quench", "ncl")), "%.0f"))
		}
	}

	// 3, 4. the interaction energy
	head := "| set | kind | L | R=1 | 2 | 3 | 4 | 6 | 8 |\n|---|---|---|---|---|---|---|---|---|\n"
	eTable := func(title, field string, showHealed bool) {
		fmt.Printf("\n%s\n\n%s", title, head)
		for _, s := range sets {
			for _, kind := range kinds {
				for _, L := range Ls {
					if len(data[s][L]) == 0 {
						continue
					}
					cells := []string{}
					for _, R := range Rs {
						if 2*R > L {
							cells = append(cells, "—")
							continue
						}
						v := collect(s, L, kind, R, field)
						c := num(median(v), "%+.2f")
						if showHealed {
							if h, n := healed(s, L, kind, R, "healed_end"); h > 0 {
								c += fmt.Sprintf(" [healed %d/%d]", h, n)
							}
						}
						cells = append(cells, c)
					}
					fmt.Printf("| %s | %s | %d | %s |\n", s, kind, L, strings.Join(cells, " | "))
				}
			}
		}
	}
	eTable("## 3. Interaction energy E_int = E_pair − E_1 − E_2 + E_0 (window at t, median over placements; t)", "Eint_w", true)
	eTable("## 4. The same, exact, after the quench (median over placements; t)", "Eint_q", false)

	// 5. the sign at contact and the ending
	fmt.Printf("\n## 5. Sign at contact (R = 1; L = 12 \\| 16) and the ending\n\n| set | GG | XX | GX | ending | energy reach R_E (12 \\| 16) | at R = 8 (L = 16) |\n|---|---|---|---|---|---|---|\n")
	letters := map[string]string{}
	for _, s := range sets {
		c := map[string]string{}
		cells := []string{}
		for _, kind := range kinds {
			a, b := classOf(s, kind, 12), classOf(s, kind, 16)
			if a == b {
				c[kind] = a
			} else {
				c[kind] = "size-dependent"
			}
			cells = append(cells, fmt.Sprintf("%s \\| %s", a, b))
		}
		r12, _ := reach(s, 12)
		r16, far := reach(s, 16)
		reachTxt := "contact (below 2t)"
		if !far {
			reachTxt = "extended"
		}
		lt := letter(c)
		if scored[s] {
			letters[s] = lt
		} else {
			lt = "not scored (" + lt + ")"
		}
		fmt.Printf("| %s | %s | %s | %d \\| %d | %s |\n", s, strings.Join(cells, " | "), lt, r12, r16, reachTxt)
	}
	agg, same, found := "", true, false
	for _, s := range sets {
		if l, ok := letters[s]; ok {
			if !found {
				agg, found = l, true
			} else if l != agg {
				same = false
			}
		}
	}
	switch {
	case !found:
		fmt.Printf("\n**Ending: no set is scored** (no medium held at L = 12 and 16).\n")
	case same:
		fmt.Printf("\n**Ending: %s** in every scored set.\n", agg)
	default:
		fmt.Printf("\n**Ending: mixed by set** (table above).\n")
	}

	// 6. the pair's less-where against the sum of the singles'
	fmt.Printf("\n## 6. The pair's less-where S_P = ⟨d_P − d_0⟩ over common quanta, and ρ_W = S_P/(S_1+S_2) (means over placements)\n\n%s", head)
	for _, s := range sets {
		for _, kind := range kinds {
			for _, L := range Ls {
				if len(data[s][L]) == 0 {
					continue
				}
				cells := []string{}
				for _, R := range Rs {
					if 2*R > L {
						cells = append(cells, "—")
						continue
					}
					sp, s1, s2 := mean(collect(s, L, kind, R, "S_P")), mean(collect(s, L, kind, R, "S_1")), mean(collect(s, L, kind, R, "S_2"))
					cells = append(cells, fmt.Sprintf("%s (ρ %s)", num(sp, "%+.4f"), num(sp/(s1+s2), "%.3f")))
				}
				fmt.Printf("| %s | %s | %d | %s |\n", s, kind, L, strings.Join(cells, " | "))
			}
		}
	}

	// 7. between the two defects
	fmt.Printf("\n## 7. Between the two defects: B_P = ⟨d_P − d_0⟩ over quanta within 2 of one site and within 2 of the other, and B_int = B_P − B_1 − B_2 (means)\n\n%s", head)
	for _, s := range sets {
		for _, kind := range kinds {
			for _, L := range Ls {
				if len(data[s][L]) == 0 {
					continue
				}
				cells := []string{}
				for _, R := range Rs {
					if 2*R > L {
						cells = append(cells, "—")
						continue
					}
					cells = append(cells, fmt.Sprintf("%s (int %s; d0 %s)", num(mean(collect(s, L, kind, R, "B_P")), "%+.3f"), num(mean(collect(s, L, kind, R, "B_int")), "%+.3f"),
						num(mean(collect(s, L, kind, R, "B_0")), "%.2f")))
				}
				fmt.Printf("| %s | %s | %d | %s |\n", s, kind, L, strings.Join(cells, " | "))
			}
		}
	}

	// 8. collapse on R/L
	fmt.Printf("\n## 8. Collapse on R/L (scored sets): I = S_P − S_1 − S_2 and B_int at R/L = 1/4 and 1/2, L = 8, 12, 16\n\n| set | kind | I at 1/4 (8, 12, 16) | I at 1/2 | B_int at 1/4 | B_int at 1/2 |\n|---|---|---|---|---|---|\n")
	quarter := map[int]int{8: 2, 12: 3, 16: 4}
	half := map[int]int{8: 4, 12: 6, 16: 8}
	for _, s := range sets {
		if !scored[s] {
			continue
		}
		for _, kind := range kinds {
			row := func(at map[int]int, field string) string {
				c := []string{}
				for _, L := range Ls {
					c = append(c, num(mean(collect(s, L, kind, at[L], field)), "%+.4f"))
				}
				return strings.Join(c, ", ")
			}
			fmt.Printf("| %s | %s | %s | %s | %s | %s |\n", s, kind, row(quarter, "I"), row(half, "I"), row(quarter, "B_int"), row(half, "B_int"))
		}
	}

	// 9. where: sign and additivity (reported with the ending)
	fmt.Printf("\n## 9. Where, reported with the ending (scored sets; L = 12, 16): sign of S_P (healed pairs excluded) and additivity over R ≥ 4\n\n| set | kind | pairs | S_P ≥ 0 | healed | mean ρ_W (R ≥ 4) | class |\n|---|---|---|---|---|---|---|\n")
	for _, s := range sets {
		if !scored[s] {
			continue
		}
		for _, kind := range kinds {
			tot, viol, nh := 0, 0, 0
			rhos := []float64{}
			for _, L := range barLs {
				for _, R := range Rs {
					if 2*R > L {
						continue
					}
					for _, p := range data[s][L] {
						r, ok := p.pairs[key(kind, R)]
						if !ok {
							continue
						}
						if r["healed_q"] == "1" {
							nh++
							continue
						}
						tot++
						if r.f("S_P") >= 0 {
							viol++
						}
					}
					if R >= 4 {
						sp, s1, s2 := mean(collect(s, L, kind, R, "S_P")), mean(collect(s, L, kind, R, "S_1")), mean(collect(s, L, kind, R, "S_2"))
						rhos = append(rhos, sp/(s1+s2))
					}
				}
			}
			rho := mean(rhos)
			class := "additive"
			switch {
			case math.IsNaN(rho):
				class = "—"
			case rho < 0.95:
				class = "sub-additive"
			case rho > 1.05:
				class = "super-additive"
			}
			fmt.Printf("| %s | %s | %d | %d | %d | %s | %s |\n", s, kind, tot, viol, nh, num(rho, "%.3f"), class)
		}
	}

	// 10. what the pairs become
	fmt.Printf("\n## 10. What the pair becomes (L = 16; clusters at the quench → at the end, and the centre distance when two; per placement)\n\n%s", head)
	for _, s := range sets {
		for _, kind := range kinds {
			L := 16
			if len(data[s][L]) == 0 {
				continue
			}
			cells := []string{}
			for _, R := range Rs {
				c := []string{}
				for _, p := range data[s][L] {
					if r, ok := p.pairs[key(kind, R)]; ok {
						t := fmt.Sprintf("%s→%s", r["ncl_q"], r["ncl_end"])
						if r.i("dcc") >= 0 {
							t += fmt.Sprintf(" d%s", r["dcc"])
						}
						c = append(c, t)
					}
				}
				cells = append(cells, strings.Join(c, ", "))
			}
			fmt.Printf("| %s | %s | %d | %s |\n", s, kind, L, strings.Join(cells, " | "))
		}
	}

	// 11. persistence and exact locality
	fmt.Printf("\n## 11. Persistence in the window, and exact locality after the quench\n\n| set | media | largest contacts changed | media with any change | pairs with E_int^q ≠ 0 at R ≥ 6 |\n|---|---|---|---|---|\n")
	for _, s := range sets {
		nm, mx, nch, nloc := 0, 0, 0, 0
		for _, L := range Ls {
			for _, p := range data[s][L] {
				for _, w := range p.windows {
					nm++
					c := w.i("changed")
					if c > mx {
						mx = c
					}
					if c > 0 {
						nch++
					}
				}
				for _, r := range p.pairs {
					if r.i("R") >= 6 && r.f("Eint_q") != 0 {
						nloc++
					}
				}
			}
		}
		fmt.Printf("| %s | %d | %d | %d | %d |\n", s, nm, mx, nch, nloc)
	}
}
