// p1ds — P1 offline: spectral dimension of free-cell contact graph.
//
//	go run ./theory/qgraph/cmd/p1ds path/to/file.fcs [more.fcs...]
//
// Builds contact graph from fcsdump -mode cells last frame (r_i+r_j contact),
// runs qgraph.ReturnDimension, prints d_s (claim: 3.00 ± 0.15 for p_rad=4).
package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"

	"scp/v93/theory/qgraph"
)

type cell struct {
	x, y, z, r float64
	em         float64
}

func main() {
	if len(os.Args) < 2 {
		fmt.Println("usage: go run ./theory/qgraph/cmd/p1ds file.fcs ...")
		os.Exit(2)
	}
	fcsdump := findFcsdump()
	fmt.Println("# P1 spectral dimension  claim ds=3±0.15 for p_rad=4")
	fmt.Println("file n_cells n_edges b1 ds_return")
	for _, fcs := range os.Args[1:] {
		cells := lastFrame(fcsdump, fcs)
		if len(cells) == 0 {
			fmt.Printf("%s EMPTY\n", filepath.Base(fcs))
			continue
		}
		// restrict to Em>0.05 support for "matter foam"
		var use []cell
		for _, c := range cells {
			if c.em > 0.05 {
				use = append(use, c)
			}
		}
		if len(use) < 8 {
			use = cells // fall back to all
		}
		g := contactGraph(use, 16)
		// sources: a few vertices of high degree
		src := []int{0}
		if g.N > 3 {
			src = []int{0, g.N / 3, 2 * g.N / 3}
		}
		ds := qgraph.ReturnDimension(g, 128, src)
		b1 := len(g.Edges) - g.N + 1 // single component approx
		fmt.Printf("%s %d %d %d %.4f\n", filepath.Base(fcs), g.N, len(g.Edges), b1, ds)
	}
}

func findFcsdump() string {
	// walk up from cwd
	wd, _ := os.Getwd()
	for d := wd; d != "/" && d != "."; d = filepath.Dir(d) {
		p := filepath.Join(d, "fcsdump")
		if _, err := os.Stat(p); err == nil {
			return p
		}
	}
	return "fcsdump"
}

func lastFrame(fcsdump, fcs string) []cell {
	cmd := exec.Command(fcsdump, "-mode", "cells", fcs)
	out, err := cmd.Output()
	if err != nil {
		return nil
	}
	var lastT float64 = -1
	var cur, last []cell
	sc := bufio.NewScanner(strings.NewReader(string(out)))
	sc.Buffer(make([]byte, 0, 1<<20), 32<<20)
	if sc.Scan() { // header
	}
	for sc.Scan() {
		f := strings.Fields(sc.Text())
		if len(f) < 14 {
			continue
		}
		t, _ := strconv.ParseFloat(f[0], 64)
		c := cell{
			x: atof(f[2]), y: atof(f[3]), z: atof(f[4]), r: atof(f[5]),
			em: atof(f[7]),
		}
		if t != lastT {
			if lastT >= 0 {
				last = cur
			}
			cur = cur[:0]
			lastT = t
		}
		cur = append(cur, c)
	}
	if len(cur) > 0 {
		last = cur
	}
	return last
}

func atof(s string) float64 { v, _ := strconv.ParseFloat(s, 64); return v }

func wr(dx, L float64) float64 { return math.Mod(dx+0.5*L, L) - 0.5*L }

func contactGraph(cells []cell, L float64) *qgraph.Graph {
	n := len(cells)
	// Build via NewPath then clear — better: empty graph
	g := &qgraph.Graph{N: n, Adj: make([][]int, n), W: make([][]float64, n)}
	for i := 0; i < n; i++ {
		g.W[i] = make([]float64, n)
	}
	// spatial bins
	maxR := 0.5
	for _, c := range cells {
		if c.r > maxR {
			maxR = c.r
		}
	}
	bin := 2.1 * maxR
	if bin < 0.5 {
		bin = 0.5
	}
	gsize := int(math.Floor(L / bin))
	if gsize < 1 {
		gsize = 1
	}
	bin = L / float64(gsize)
	type key struct{ x, y, z int }
	heads := map[key][]int{}
	for i, c := range cells {
		kx := int(c.x/bin) % gsize
		if kx < 0 {
			kx += gsize
		}
		ky := int(c.y/bin) % gsize
		if ky < 0 {
			ky += gsize
		}
		kz := int(c.z/bin) % gsize
		if kz < 0 {
			kz += gsize
		}
		k := key{kx, ky, kz}
		heads[k] = append(heads[k], i)
	}
	add := func(i, j int) {
		if i > j {
			i, j = j, i
		}
		for _, e := range g.Edges {
			if e[0] == i && e[1] == j {
				return
			}
		}
		g.Edges = append(g.Edges, [2]int{i, j})
		g.Adj[i] = append(g.Adj[i], j)
		g.Adj[j] = append(g.Adj[j], i)
		g.W[i][j] = 1
		g.W[j][i] = 1
	}
	for k, mem := range heads {
		for dx := -1; dx <= 1; dx++ {
			for dy := -1; dy <= 1; dy++ {
				for dz := -1; dz <= 1; dz++ {
					ok := key{(k.x + dx + gsize) % gsize, (k.y + dy + gsize) % gsize, (k.z + dz + gsize) % gsize}
					oth := heads[ok]
					for _, i := range mem {
						for _, j := range oth {
							if j <= i {
								continue
							}
							ddx := wr(cells[j].x-cells[i].x, L)
							ddy := wr(cells[j].y-cells[i].y, L)
							ddz := wr(cells[j].z-cells[i].z, L)
							d2 := ddx*ddx + ddy*ddy + ddz*ddz
							ct := cells[i].r + cells[j].r
							if d2 < ct*ct && d2 > 1e-20 {
								add(i, j)
							}
						}
					}
				}
			}
		}
	}
	return g
}
