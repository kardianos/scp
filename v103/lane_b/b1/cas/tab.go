// tab.go — v103 lane b1: where classical begins (../PREREG.md §5) from the
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
	"sort"
	"strconv"
	"strings"
)

var rules = []string{"none", "a", "b", "c", "ab", "ac", "bc", "abc"}

type Pt struct {
	NB               int
	G1, Sig, Beta, Z float64
	Var0, Var1, P    float64
	Check            string
}

func num(s string) float64 { v, _ := strconv.ParseFloat(s, 64); return v }
func rx(s, pat string) float64 {
	m := regexp.MustCompile(pat).FindStringSubmatch(s)
	if m == nil {
		return math.NaN()
	}
	return num(m[1])
}

const F = `(-?[0-9.eE+-]+)`

func parse(path string) *Pt {
	f, err := os.Open(path)
	if err != nil {
		return nil
	}
	defer f.Close()
	p := &Pt{}
	done := false
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		l := sc.Text()
		switch {
		case strings.HasPrefix(l, "B1 model"):
			p.NB = int(rx(l, `nb=(\d+)`))
		case strings.HasPrefix(l, "MEDIUM"):
			p.G1 = rx(l, `g1\(L/2\)=`+F)
			p.Sig = rx(l, `sigmaN/nb=`+F)
		case strings.HasPrefix(l, "VOID"):
			p.Beta = rx(l, `beta=`+F)
			p.Z = rx(l, `Z=`+F)
			p.Var0 = rx(l, `var0=`+F)
			p.Var1 = rx(l, `var1=`+F)
			p.P = rx(l, `p=`+F)
		case strings.HasPrefix(l, "CHECK"):
			p.Check = l
			done = true
		}
	}
	if !done {
		return nil
	}
	return p
}

func ladder(r string, L int) []*Pt {
	var out []*Pt
	for _, nb := range []int{1, 2, 3, 4, 6, 8} {
		if p := parse(fmt.Sprintf("out/b1_%s_L%d_n%d.out", r, L, nb)); p != nil {
			out = append(out, p)
		}
	}
	sort.Slice(out, func(i, j int) bool { return out[i].NB < out[j].NB })
	return out
}

func main() {
	fmt.Printf("## b1 ladders\n\n| rules | L | n̄ | g1(L/2) | σ_N/n̄ | β (void discreteness) | Z | regime |\n|---|---|---|---|---|---|---|---|\n")
	type key struct {
		r string
		L int
	}
	top := map[key]bool{}
	anyTwo, anyTop := false, false
	for _, r := range rules {
		for _, L := range []int{4, 6} {
			lad := ladder(r, L)
			for i, p := range lad {
				cb, dv := p.G1 >= 0.9, p.Beta >= 0.5
				reg := "neither"
				switch {
				case cb && dv:
					reg = "**two-regime**"
					anyTwo = true
				case cb:
					reg = "classical bulk"
				case dv:
					reg = "discrete void"
				}
				if i == len(lad)-1 {
					top[key{r, L}] = cb && dv
				}
				fmt.Printf("| %s | %d | %d | %.4f | %.4f | %.4f | %.4f | %s |\n", r, L, p.NB, p.G1, p.Sig, p.Beta, p.Z, reg)
			}
		}
	}
	// U: two-regime at the top of the ladder, at L = 4 and at L = 6 where run
	unified := ""
	for _, r := range rules {
		ok4 := top[key{r, 4}]
		_, has6 := top[key{r, 6}]
		ok := ok4 && (!has6 || top[key{r, 6}])
		if len(ladder(r, 6)) == 0 {
			ok = ok4
		}
		if ok {
			anyTop = true
			unified += r + " "
		}
	}
	ending := "C"
	switch {
	case anyTop:
		ending = "U (" + strings.TrimSpace(unified) + ")"
	case anyTwo:
		ending = "P"
	}
	fmt.Printf("\n**Ending: %s**\n\n## checks\n\n", ending)
	for _, r := range rules {
		for _, L := range []int{4, 6} {
			for _, p := range ladder(r, L) {
				fmt.Printf("- %s L=%d n̄=%d: %s\n", r, L, p.NB, p.Check)
			}
		}
	}
}
