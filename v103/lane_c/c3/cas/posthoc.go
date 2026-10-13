// v103 lane c step c3: post-hoc readings (not registered; written after the grid, labelled as such in DECIDE.md).
//
//	go run posthoc.go > out/posthoc.out
//
// 1. Persistence: window contact changes by set and kind, and for the pairs by separation.
// 2. The where-interaction I = S_P - S_1 - S_2 relative to |S_1 + S_2|, by kind, L and R (scored sets, unhealed pairs).
package main

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

var scored = map[string]bool{"a": true, "b": true, "ab": true, "ac": true, "bc": true, "abc": true}
var setOrder = []string{"a", "b", "ab", "ac", "bc", "abc"}
var Rs = []int{1, 2, 3, 4, 6, 8}

func fields(line string) map[string]string {
	m := map[string]string{}
	for _, f := range strings.Fields(line)[1:] {
		if i := strings.IndexByte(f, '='); i > 0 {
			m[f[:i]] = f[i+1:]
		}
	}
	return m
}

func num(s string) float64 {
	v, err := strconv.ParseFloat(s, 64)
	if err != nil {
		panic(err)
	}
	return v
}

type acc struct {
	n, changed, max int
}

type wsum struct {
	n      int
	I, S12 float64
}

func main() {
	files, _ := filepath.Glob("out/c3_*.out")
	sort.Strings(files)
	win := map[string]*acc{}    // set|kind
	winR := map[string]*acc{}   // kind|R (pairs, scored sets)
	where := map[string]*wsum{} // kind|L|R
	for _, fn := range files {
		f, err := os.Open(fn)
		if err != nil {
			panic(err)
		}
		set, L := "", 0
		sc := bufio.NewScanner(f)
		sc.Buffer(make([]byte, 1<<20), 1<<20)
		for sc.Scan() {
			line := sc.Text()
			switch {
			case strings.HasPrefix(line, "RUN "):
				m := fields(line)
				set = m["set"]
				L = int(num(m["L"]))
			case strings.HasPrefix(line, "WINDOW "):
				if !scored[set] {
					continue
				}
				m := fields(line)
				ch := int(num(m["changed"]))
				k := set + "|" + m["kind"]
				if win[k] == nil {
					win[k] = &acc{}
				}
				a := win[k]
				a.n++
				if ch > 0 {
					a.changed++
				}
				if ch > a.max {
					a.max = ch
				}
				kind := m["kind"]
				if kind == "GG" || kind == "XX" || kind == "GX" {
					k2 := kind + "|" + m["R"]
					if winR[k2] == nil {
						winR[k2] = &acc{}
					}
					winR[k2].n++
					if ch > 0 {
						winR[k2].changed++
					}
				}
			case strings.HasPrefix(line, "PAIR "):
				if !scored[set] {
					continue
				}
				m := fields(line)
				if m["healed_end"] != "0" {
					continue
				}
				k := fmt.Sprintf("%s|%d|%s", m["kind"], L, m["R"])
				if where[k] == nil {
					where[k] = &wsum{}
				}
				w := where[k]
				w.n++
				w.I += num(m["I"])
				w.S12 += num(m["S_1"]) + num(m["S_2"])
			}
		}
		f.Close()
	}

	fmt.Println("# v103 c3 — post-hoc readings (not registered)")
	fmt.Println()
	fmt.Println("## 1. Persistence in the window, scored sets: media with any contact changed (of all), and the largest change")
	fmt.Println()
	kinds := []string{"G1", "X1", "G2", "X2", "GG", "XX", "GX"}
	fmt.Printf("| set | %s |\n", strings.Join(kinds, " | "))
	fmt.Printf("|---|%s\n", strings.Repeat("---|", len(kinds)))
	for _, s := range setOrder {
		var cells []string
		for _, kd := range kinds {
			a := win[s+"|"+kd]
			if a == nil {
				cells = append(cells, "(missing)")
				continue
			}
			cells = append(cells, fmt.Sprintf("%d/%d (max %d)", a.changed, a.n, a.max))
		}
		fmt.Printf("| %s | %s |\n", s, strings.Join(cells, " | "))
	}
	fmt.Println()
	fmt.Println("Pairs by separation (all scored sets and sizes): media with any contact changed, of all")
	fmt.Println()
	fmt.Println("| kind | R=1 | 2 | 3 | 4 | 6 | 8 |")
	fmt.Println("|---|---|---|---|---|---|---|")
	for _, kd := range []string{"GG", "XX", "GX"} {
		var cells []string
		for _, R := range Rs {
			a := winR[fmt.Sprintf("%s|%d", kd, R)]
			if a == nil {
				cells = append(cells, "—")
				continue
			}
			cells = append(cells, fmt.Sprintf("%d/%d", a.changed, a.n))
		}
		fmt.Printf("| %s | %s |\n", kd, strings.Join(cells, " | "))
	}
	fmt.Println()
	fmt.Println("## 2. The where-interaction relative to the singles: I/|S_1 + S_2| (means over scored sets and placements; unhealed pairs)")
	fmt.Println()
	fmt.Println("| kind | L | R=1 | 2 | 3 | 4 | 6 | 8 |")
	fmt.Println("|---|---|---|---|---|---|---|---|")
	for _, kd := range []string{"GG", "XX", "GX"} {
		for _, L := range []int{8, 12, 16} {
			var cells []string
			for _, R := range Rs {
				w := where[fmt.Sprintf("%s|%d|%d", kd, L, R)]
				if w == nil || w.n == 0 {
					cells = append(cells, "—")
					continue
				}
				cells = append(cells, fmt.Sprintf("%.3f", -w.I/w.S12))
			}
			fmt.Printf("| %s | %d | %s |\n", kd, L, strings.Join(cells, " | "))
		}
	}
}
