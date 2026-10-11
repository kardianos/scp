// tab.go — v103 lane a1: the survival ledger (../PREREG.md §5) from the
// bx outputs in out/.
//
//	go run tab.go
package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"regexp"
	"strconv"
	"strings"
)

var rules = []string{"none", "a", "b", "c", "ab", "ac", "bc", "abc"}

type Run struct {
	M                           int
	Gap, S2G, W, Zmin, S2V, DAB float64
	HasEx, HasBr                bool
	ZminEx, BrE, BrW, Thr, R    float64
	Pdouble, Bind, EvoidMin     float64
	Check                       string
	Done                        bool
}

func num(s string) float64 {
	if s == "inf" {
		return math.Inf(1)
	}
	v, err := strconv.ParseFloat(s, 64)
	if err != nil {
		return math.NaN()
	}
	return v
}

func rx(s, pat string) float64 {
	m := regexp.MustCompile(pat).FindStringSubmatch(s)
	if m == nil {
		return math.NaN()
	}
	return num(m[1])
}

const F = `(-?[0-9.eE+-]+|inf|nan)`

func parse(path string) *Run {
	f, err := os.Open(path)
	if err != nil {
		return nil
	}
	defer f.Close()
	r := &Run{}
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 1<<20), 1<<20)
	for sc.Scan() {
		l := sc.Text()
		switch {
		case strings.HasPrefix(l, "A1 model"):
			r.M = int(rx(l, `M=(\d+)`))
		case strings.HasPrefix(l, "MEDIUM"):
			r.Gap = rx(l, `gap=`+F)
			r.S2G = rx(l, `S2=`+F)
		case strings.HasPrefix(l, "VOID W="):
			r.W = rx(l, `W=`+F)
			r.Zmin = rx(l, `Zmin=`+F)
		case strings.HasPrefix(l, "VOID band"):
			mn := math.Inf(1)
			for _, m := range regexp.MustCompile(`E=`+F).FindAllStringSubmatch(l, -1) {
				if v := num(m[1]); v < mn {
					mn = v
				}
			}
			r.EvoidMin = mn
		case strings.HasPrefix(l, "VOID S2="):
			r.S2V = rx(l, `VOID S2=`+F)
		case strings.HasPrefix(l, "AB "):
			r.DAB = rx(l, `dAB=`+F)
		case strings.HasPrefix(l, "EXCESS Zmin"):
			r.HasEx = true
			r.ZminEx = rx(l, `Zmin=`+F)
		case strings.HasPrefix(l, "BREATH E="):
			r.HasBr = true
			r.BrE = rx(l, `E=`+F)
			r.BrW = rx(l, `weight=`+F)
			r.Thr = rx(l, `threshold=`+F)
			r.R = rx(l, `R=`+F)
		case strings.HasPrefix(l, "TWOVOID"):
			r.Pdouble = rx(l, `Pdouble=`+F)
			r.Bind = rx(l, `bind=`+F)
		case strings.HasPrefix(l, "CHECK"):
			r.Check = l
			r.Done = true
		}
	}
	return r
}

type Bars struct{ A1, A2, A4, A5, A6 bool }

func bars(r *Run) Bars {
	var b Bars
	if r == nil || !r.Done {
		return b
	}
	b.A1 = r.Gap >= 1e-6 && r.W >= 1e-3 && r.Zmin >= 0.5 && math.Abs(r.DAB) >= 1e-8
	b.A2 = r.Pdouble <= 0.01
	b.A4 = r.M == 2 && r.Gap >= 1e-6 && r.S2G <= 1e-6 && math.Abs(r.S2V-0.75) <= 1e-6
	b.A5 = r.HasEx && r.ZminEx >= 0.5
	b.A6 = r.HasBr && r.R < 1 && r.BrW >= 0.1
	return b
}

func yn(v bool) string {
	if v {
		return "✔"
	}
	return "✘"
}

func main() {
	fmt.Printf("## a1 measurements\n\n| rules | L | medium gap | medium S² | void W | Z(min) | void S² | dAB | excess Z(min) | breath E / thr (R), weight | P_double | two-void bind |\n|---|---|---|---|---|---|---|---|---|---|---|---|\n")
	all := map[string][2]Bars{}
	for _, r := range rules {
		var bb [2]Bars
		for i, L := range []int{6, 8} {
			run := parse(fmt.Sprintf("out/a1_%s_L%d.out", r, L))
			if run == nil || !run.Done {
				fmt.Printf("| %s | %d | (missing) | | | | | | | | | |\n", r, L)
				continue
			}
			ex, br := "none", "none (O|G⟩ = 0)"
			if run.HasEx {
				ex = fmt.Sprintf("%.4f", run.ZminEx)
			}
			if run.HasBr {
				br = fmt.Sprintf("%.4f / %.4f (%.3f), %.3f", run.BrE, run.Thr, run.R, run.BrW)
			}
			s2v := "—"
			if run.M == 2 {
				s2v = fmt.Sprintf("%.6f", run.S2V)
			}
			fmt.Printf("| %s | %d | %.3g | %.6f | %.4f | %.4f | %s | %.3e | %s | %s | %.4f | %.4f |\n", r, L, run.Gap, run.S2G, run.W, run.Zmin, s2v, run.DAB, ex, br, run.Pdouble, run.Bind)
			bb[i] = bars(run)
		}
		all[r] = bb
	}
	fmt.Printf("\n## a1 ledger (survives = bar holds at L = 6 and L = 8)\n\n| rules | A1 wave | A2 exclusion | A4 spin ½ | A5 antiparticle | A6 breathing line | count |\n|---|---|---|---|---|---|---|\n")
	best, bestR := -1, ""
	for _, r := range rules {
		b := all[r]
		s := Bars{b[0].A1 && b[1].A1, b[0].A2 && b[1].A2, b[0].A4 && b[1].A4, b[0].A5 && b[1].A5, b[0].A6 && b[1].A6}
		n := 0
		for _, v := range []bool{s.A1, s.A2, s.A4, s.A5, s.A6} {
			if v {
				n++
			}
		}
		if n > best {
			best, bestR = n, r
		} else if n == best {
			bestR += ", " + r
		}
		fmt.Printf("| %s | %s | %s | %s | %s | %s | %d |\n", r, yn(s.A1), yn(s.A2), yn(s.A4), yn(s.A5), yn(s.A6), n)
	}
	ending := fmt.Sprintf("P(%d)", best)
	if best == 5 {
		ending = "F"
	}
	fmt.Printf("\n**Ending: %s** (best: %s)\n", ending, bestR)
	fmt.Printf("\n## checks\n\n")
	for _, r := range rules {
		for _, L := range []int{6, 8} {
			if run := parse(fmt.Sprintf("out/a1_%s_L%d.out", r, L)); run != nil {
				fmt.Printf("- %s L=%d: %s\n", r, L, run.Check)
			}
		}
	}
}
