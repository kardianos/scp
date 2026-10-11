// tab.go — v103 lane c2: ledger, candidate roles and the registered ending
// (../PREREG.md §5) from the gap1d outputs in out/.
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

type Pair struct {
	Wmin, Rho float64
	Prof      []string
}

type Run struct {
	Bgg, Bee    float64
	P           map[string]Pair
	Mach        string
	Check       string
	Sectors     []string
	Done        bool
	MaxRes      float64
	Trunc, Viol int
}

func num(s string) float64 {
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

const F = `(-?[0-9.eE+-]+|nan|inf|-nan)`

func parse(path string) *Run {
	f, err := os.Open(path)
	if err != nil {
		return nil
	}
	defer f.Close()
	r := &Run{P: map[string]Pair{}}
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 1<<20), 1<<20)
	for sc.Scan() {
		l := sc.Text()
		switch {
		case strings.HasPrefix(l, "SECTOR"):
			r.Sectors = append(r.Sectors, l)
		case strings.HasPrefix(l, "BIND"):
			r.Bgg = rx(l, `gg=`+F)
			r.Bee = rx(l, `ee=`+F)
		case strings.HasPrefix(l, "PAIR"):
			fs := strings.Fields(l)
			p := r.P[fs[1]]
			p.Wmin = rx(l, `Wmin=`+F)
			p.Rho = rx(l, `rho=`+F)
			r.P[fs[1]] = p
		case strings.HasPrefix(l, "PROFILE"):
			fs := strings.Fields(l)
			p := r.P[fs[1]]
			for _, x := range fs[2:] {
				// D=k:p=..,pref=..,V=..
				d := strings.SplitN(x, ":", 2)
				v := rx(x, `V=`+F)
				p.Prof = append(p.Prof, fmt.Sprintf("%s:%s", strings.TrimPrefix(d[0], "D="), fv(v)))
			}
			r.P[fs[1]] = p
		case strings.HasPrefix(l, "MACH"):
			r.Mach = l
		case strings.HasPrefix(l, "CHECK"):
			r.Check = l
			r.MaxRes = rx(l, `maxresid=`+F)
			r.Trunc = int(rx(l, `trunc=(\d+)`))
			r.Viol = int(rx(l, `moment_viol=(\d+)`))
			r.Done = true
		}
	}
	return r
}

func fv(v float64) string {
	switch {
	case math.IsInf(v, 1):
		return "∞"
	case math.IsNaN(v):
		return "—"
	}
	return fmt.Sprintf("%.2f", v)
}

func f3(v float64) string {
	if math.IsNaN(v) {
		return "—"
	}
	return fmt.Sprintf("%+.3f", v)
}

func name(a, b, c, d int) string {
	s := ""
	for i, on := range []int{a, b, c, d} {
		if on == 1 {
			s += string("abcd"[i])
		}
	}
	if s == "" {
		s = "none"
	}
	return s
}

// attraction of a like pair at one size: B <= -1e-6 and rho <= 0.95
func attracts(B, rho float64) bool {
	return !math.IsNaN(B) && !math.IsNaN(rho) && B <= -1e-6 && rho <= 0.95
}
func repels(rho float64) bool { return !math.IsNaN(rho) && rho >= 1.05 }

func verdict(B, rho float64) string {
	switch {
	case math.IsNaN(rho):
		return "—"
	case attracts(B, rho):
		return "attract"
	case repels(rho):
		return "repel"
	}
	return "neither"
}

func main() {
	Ls := []int{12, 16}
	type key struct{ a, b, c, d int }
	runs := map[key][2]*Run{}
	missing := 0
	for m := 0; m < 16; m++ {
		k := key{m >> 3 & 1, m >> 2 & 1, m >> 1 & 1, m & 1}
		var rr [2]*Run
		for i, L := range Ls {
			rr[i] = parse(fmt.Sprintf("out/c2_a%db%dc%dd%d_L%d.out", k.a, k.b, k.c, k.d, L))
			if rr[i] == nil || !rr[i].Done {
				missing++
			}
		}
		runs[k] = rr
	}
	fmt.Printf("## c2 measurements (L = 12 | L = 16)\n\n")
	fmt.Printf("| set | B_gg | ρ_gg (W_min) | B_ee | ρ_ee (W_min) | ρ_ge, reported | gg | ee |\n|---|---|---|---|---|---|---|---|\n")
	type res struct{ gg, ee bool }
	R := map[key]res{}
	order := []key{}
	for _, c := range []int{0, 1} {
		for _, a := range []int{1, 0} {
			for _, b := range []int{0, 1} {
				for _, d := range []int{0, 1} {
					order = append(order, key{a, b, c, d})
				}
			}
		}
	}
	for _, k := range order {
		rr := runs[k]
		cellB := func(get func(*Run) float64) string {
			s := []string{}
			for _, r := range rr {
				if r == nil || !r.Done {
					s = append(s, "(missing)")
					continue
				}
				s = append(s, f3(get(r)))
			}
			return strings.Join(s, " \\| ")
		}
		cellP := func(t string) string {
			s := []string{}
			for _, r := range rr {
				if r == nil || !r.Done {
					s = append(s, "(missing)")
					continue
				}
				p := r.P[t]
				if math.IsNaN(p.Rho) {
					s = append(s, "—")
					continue
				}
				s = append(s, fmt.Sprintf("%.3f (%.3f)", p.Rho, p.Wmin))
			}
			return strings.Join(s, " \\| ")
		}
		gg, ee := true, true
		vg, ve := []string{}, []string{}
		for _, r := range rr {
			if r == nil || !r.Done {
				gg, ee = false, false
				continue
			}
			gg = gg && attracts(r.Bgg, r.P["gg"].Rho)
			ee = ee && attracts(r.Bee, r.P["ee"].Rho)
			vg = append(vg, verdict(r.Bgg, r.P["gg"].Rho))
			ve = append(ve, verdict(r.Bee, r.P["ee"].Rho))
		}
		R[k] = res{gg, ee}
		fmt.Printf("| %s | %s | %s | %s | %s | %s | %s | %s |\n", name(k.a, k.b, k.c, k.d),
			cellB(func(r *Run) float64 { return r.Bgg }), cellP("gg"),
			cellB(func(r *Run) float64 { return r.Bee }), cellP("ee"), cellP("ge"),
			strings.Join(vg, ", "), strings.Join(ve, ", "))
	}

	fmt.Printf("\n## like-pair attraction at both sizes\n\n| set | gaps attract | excesses attract |\n|---|---|---|\n")
	anyU, anyS, anyA := []string{}, []string{}, []string{}
	for _, k := range order {
		r := R[k]
		yn := func(v bool) string {
			if v {
				return "✔"
			}
			return "✘"
		}
		ee := yn(r.ee)
		if k.a == 1 {
			ee = "— (no excess)"
		}
		fmt.Printf("| %s | %s | %s |\n", name(k.a, k.b, k.c, k.d), yn(r.gg), ee)
		n := name(k.a, k.b, k.c, k.d)
		switch {
		case k.a == 0 && r.gg && r.ee:
			anyU = append(anyU, n)
		case k.a == 0 && (r.gg != r.ee):
			anyS = append(anyS, n)
		case k.a == 1 && r.gg:
			anyA = append(anyA, n)
		}
	}
	ending := "N"
	switch {
	case len(anyU) > 0:
		ending = "U (" + strings.Join(anyU, ", ") + ")"
	case len(anyS) > 0:
		ending = "S (" + strings.Join(anyS, ", ") + ")"
	case len(anyA) > 0:
		ending = "A (" + strings.Join(anyA, ", ") + ")"
	}
	fmt.Printf("\nsign-dependent sets (one like pair only): %v; gaps-only sets (a on): %v\n", anyS, anyA)
	fmt.Printf("\n**Ending: %s**\n", ending)

	// candidate roles over the sets where gap pairs attract
	fmt.Printf("\n## candidate roles (gap pairs attracting at both sizes)\n\n")
	att := []key{}
	for _, k := range order {
		if R[k].gg {
			att = append(att, k)
		}
	}
	for i, c := range []string{"a", "b", "c", "d"} {
		in, out := 0, 0
		for _, k := range att {
			v := []int{k.a, k.b, k.c, k.d}[i]
			if v == 1 {
				in++
			} else {
				out++
			}
		}
		role := "mixed"
		switch {
		case len(att) == 0:
			role = "—"
		case out == 0:
			role = "present in every attracting set (necessary)"
		case in == 0:
			role = "absent from every attracting set"
		}
		fmt.Printf("- %s: on in %d, off in %d of %d attracting sets — %s\n", c, in, out, len(att), role)
	}

	fmt.Printf("\n## Mach structure (depends on a and b only)\n\n")
	for _, k := range order {
		if k.c == 1 && k.d == 0 {
			for i, r := range runs[k] {
				if r != nil {
					fmt.Printf("- a=%d b=%d L=%d: %s\n", k.a, k.b, Ls[i], strings.TrimPrefix(r.Mach, "MACH "))
				}
			}
		}
	}

	fmt.Printf("\n## profiles V(D), D in quanta (reported, not scored)\n\n")
	for _, k := range order {
		for i, r := range runs[k] {
			if r == nil {
				continue
			}
			for _, t := range []string{"gg", "ee"} {
				if p := r.P[t]; !math.IsNaN(p.Rho) {
					fmt.Printf("- %s L=%d %s: %s\n", name(k.a, k.b, k.c, k.d), Ls[i], t, strings.Join(p.Prof, " "))
				}
			}
		}
	}

	// c-off check: b on and b off agree exactly (PREREG §2)
	fmt.Printf("\n## checks\n\n")
	cok := true
	for _, k := range order {
		if k.c != 0 || k.b != 0 {
			continue
		}
		k2 := key{k.a, 1, 0, k.d}
		for i := range Ls {
			r1, r2 := runs[k][i], runs[k2][i]
			if r1 == nil || r2 == nil {
				cok = false
				continue
			}
			same := fmt.Sprint(r1.Bgg, r1.Bee, r1.P["gg"], r1.P["ee"], r1.P["ge"]) == fmt.Sprint(r2.Bgg, r2.Bee, r2.P["gg"], r2.P["ee"], r2.P["ge"])
			if !same {
				cok = false
				fmt.Printf("- c-off mismatch a=%d d=%d L=%d\n", k.a, k.d, Ls[i])
			}
		}
	}
	fmt.Printf("- c-off b-on/b-off agreement (all four pairs, both sizes): %v\n", cok)
	worst, trunc, viol := 0.0, 0, 0
	for _, rr := range runs {
		for _, r := range rr {
			if r == nil {
				continue
			}
			if r.MaxRes > worst {
				worst = r.MaxRes
			}
			trunc += r.Trunc
			viol += r.Viol
		}
	}
	fmt.Printf("- worst residual %.1e; truncated manifolds %d; moment violations %d; missing runs %d\n", worst, trunc, viol, missing)
}
