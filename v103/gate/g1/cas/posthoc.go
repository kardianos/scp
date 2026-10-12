// posthoc.go — v103 gate g1, POST-HOC readings (not registered; ../DECIDE.md
// labels them so). The registered tables are tab.go's.
//
// 1. How the mean-distance change of one gap (and one excess) scales with the
//    size of the medium: δ⟨d⟩ at every L.
// 2. Whether the gap's distance footprint Δm(r) collapses on r/L: Δm read at
//    r/L = 1/8, 1/4, 3/8, 1/2 for every L. A footprint that collapses on r/L
//    is set by the size of the closed medium, not by a decay of its own.
//
//	go run posthoc.go > out/posthoc.out
package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
)

type acc struct{ n, dm map[int]float64 }

func read(path string) (mM, mG, mX float64, prof map[string]*acc, ok bool) {
	f, err := os.Open(path)
	if err != nil {
		return
	}
	defer f.Close()
	prof = map[string]*acc{"G": {map[int]float64{}, map[int]float64{}}, "X": {map[int]float64{}, map[int]float64{}}}
	sum := map[string]float64{}
	cnt := map[string]float64{}
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 1<<20), 1<<20)
	for sc.Scan() {
		fs := strings.Fields(sc.Text())
		if len(fs) == 0 {
			continue
		}
		v := map[string]string{}
		for _, x := range fs[1:] {
			if p := strings.SplitN(x, "=", 2); len(p) == 2 {
				v[p[0]] = p[1]
			}
		}
		num := func(k string) float64 { x, _ := strconv.ParseFloat(v[k], 64); return x }
		switch fs[0] {
		case "SAMPLE":
			sum[v["sec"]] += num("meand")
			cnt[v["sec"]]++
		case "PROF":
			if a := prof[v["sec"]]; a != nil {
				r := int(num("r"))
				a.n[r] += num("n")
				a.dm[r] += num("dm") * num("n")
			}
		case "DONE":
			ok = true
		}
	}
	mM, mG, mX = sum["M"]/cnt["M"], sum["G"]/cnt["G"], sum["X"]/cnt["X"]
	return
}

func main() {
	sets := []struct {
		name    string
		a, b, c int
	}{{"a", 1, 0, 0}, {"ac", 1, 0, 1}, {"b", 0, 1, 0}, {"bc", 0, 1, 1}, {"ab", 1, 1, 0}, {"abc", 1, 1, 1}}
	Ls := []int{6, 8, 10, 12, 16}
	fmt.Printf("## POST-HOC 1: mean-distance change δ⟨d⟩ = ⟨d⟩_sector − ⟨d⟩_M, by size\n\n| set | sector | L=6 | L=8 | L=10 | L=12 | L=16 |\n|---|---|---|---|---|---|---|\n")
	type rec struct {
		mM  float64
		pG  *acc
		ok  bool
		row map[string]float64
	}
	data := map[string]map[int]rec{}
	for _, s := range sets {
		data[s.name] = map[int]rec{}
		for _, L := range Ls {
			mM, mG, mX, prof, ok := read(fmt.Sprintf("out/g1_a%db%dc%d_L%d.out", s.a, s.b, s.c, L))
			r := rec{mM: mM, ok: ok, row: map[string]float64{"G": mG - mM, "X": mX - mM}}
			if prof != nil {
				r.pG = prof["G"]
			}
			data[s.name][L] = r
		}
		for _, sec := range []string{"G", "X"} {
			cells := []string{}
			for _, L := range Ls {
				r := data[s.name][L]
				if !r.ok {
					cells = append(cells, "(missing)")
					continue
				}
				cells = append(cells, fmt.Sprintf("%+.4f", r.row[sec]))
			}
			fmt.Printf("| %s | %s | %s |\n", s.name, sec, strings.Join(cells, " | "))
		}
	}
	fracs := []float64{0.125, 0.25, 0.375, 0.5}
	fmt.Printf("\n## POST-HOC 2: the gap's distance footprint Δm against r/L (Δm(r) = m̄(r) − m_ref, interpolated)\n\n| set | L | r/L = 1/8 | 1/4 | 3/8 | 1/2 |\n|---|---|---|---|---|---|\n")
	for _, s := range sets {
		for _, L := range Ls {
			r := data[s.name][L]
			if !r.ok || r.pG == nil {
				continue
			}
			dm := func(x float64) float64 {
				lo := int(math.Floor(x))
				hi := lo + 1
				get := func(k int) float64 {
					if r.pG.n[k] == 0 {
						return math.NaN()
					}
					return r.pG.dm[k]/r.pG.n[k] - r.mM
				}
				a, b := get(lo), get(hi)
				return a + (b-a)*(x-float64(lo))
			}
			cells := []string{}
			for _, f := range fracs {
				cells = append(cells, fmt.Sprintf("%+.4f", dm(f*float64(L))))
			}
			fmt.Printf("| %s | %d | %s |\n", s.name, L, strings.Join(cells, " | "))
		}
	}
}
