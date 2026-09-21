package main

import (
	"bufio"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
)

// ProgMeters: progressive / Bohmian-method offline meters from last fcs frame.
type ProgMeters struct {
	EmTot, EeTot, EmMax float64
	PR, FTag            float64
	JAbs, Jx, R2f       float64 // field
	JdAbs               float64 // dense |J|
	NCells, NLinks      int
	OK                  bool
	Note                string
}

type cell struct {
	x, y, z, r float64
	em, ee     float64
	tag        float64
	fa1, fa2   float64
	th2        float64
}

// measureFCS runs fcsdump -mode cells and computes progressive meters on the last frame.
func measureFCS(fcsdump, fcs string, box float64) ProgMeters {
	var m ProgMeters
	if fcs == "" {
		m.Note = "no-fcs"
		return m
	}
	if _, err := os.Stat(fcs); err != nil {
		m.Note = "missing-fcs"
		return m
	}
	cmd := exec.Command(fcsdump, "-mode", "cells", fcs)
	out, err := cmd.StdoutPipe()
	if err != nil {
		m.Note = err.Error()
		return m
	}
	if err := cmd.Start(); err != nil {
		m.Note = err.Error()
		return m
	}
	// keep only last frame
	var lastT float64 = -1
	var cur []cell
	var last []cell
	sc := bufio.NewScanner(out)
	sc.Buffer(make([]byte, 0, 256*1024), 16*1024*1024)
	// header
	if sc.Scan() {
		// t i x y z r es em ee xload tag fa1 fa2 th2
	}
	for sc.Scan() {
		p := strings.Fields(sc.Text())
		if len(p) < 14 {
			continue
		}
		t, _ := strconv.ParseFloat(p[0], 64)
		c := cell{
			x: atof(p[2]), y: atof(p[3]), z: atof(p[4]), r: atof(p[5]),
			em: atof(p[7]), ee: atof(p[8]), tag: atof(p[10]),
			fa1: atof(p[11]), fa2: atof(p[12]), th2: atof(p[13]),
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
	_ = cmd.Wait()
	if len(cur) > 0 {
		last = cur
	}
	if len(last) == 0 {
		m.Note = "empty"
		return m
	}
	return computeMeters(last, box)
}

func atof(s string) float64 {
	v, _ := strconv.ParseFloat(s, 64)
	return v
}

func wr(dx, box float64) float64 {
	return math.Mod(dx+0.5*box, box) - 0.5*box
}

func computeMeters(cells []cell, box float64) ProgMeters {
	var m ProgMeters
	m.NCells = len(cells)
	m.OK = true
	var em2, tagEm float64
	for i := range cells {
		m.EmTot += cells[i].em
		m.EeTot += cells[i].ee
		if cells[i].em > m.EmMax {
			m.EmMax = cells[i].em
		}
		em2 += cells[i].em * cells[i].em
		if cells[i].tag > 0.5 {
			tagEm += cells[i].em
		}
	}
	if m.EmTot > 0 {
		m.PR = m.EmTot * m.EmTot / em2
		m.FTag = tagEm / m.EmTot
	}
	// spatial bins for contacts
	maxR := 0.5
	for i := range cells {
		if cells[i].r > maxR {
			maxR = cells[i].r
		}
	}
	binSz := 2.1 * maxR
	if binSz < 0.5 {
		binSz = 0.5
	}
	g := int(math.Floor(box / binSz))
	if g < 1 {
		g = 1
	}
	binSz = box / float64(g)
	type key struct{ x, y, z int }
	heads := map[key][]int{}
	for i := range cells {
		kx := int(cells[i].x/binSz) % g
		ky := int(cells[i].y/binSz) % g
		kz := int(cells[i].z/binSz) % g
		if kx < 0 {
			kx += g
		}
		if ky < 0 {
			ky += g
		}
		if kz < 0 {
			kz += g
		}
		k := key{kx, ky, kz}
		heads[k] = append(heads[k], i)
	}
	var jAbs, jX, wF, jdAbs, wD float64
	var cR, sR, wR float64 // field global R2
	for i := range cells {
		ei := cells[i].fa1*cells[i].fa1 + cells[i].fa2*cells[i].fa2
		if ei > 1e-14 {
			ph := math.Atan2(cells[i].fa2, cells[i].fa1)
			cR += ei * math.Cos(ph)
			sR += ei * math.Sin(ph)
			wR += ei
		}
	}
	if wR > 0 {
		m.R2f = math.Hypot(cR, sR) / wR
	}
	for k, mem := range heads {
		for dx := -1; dx <= 1; dx++ {
			for dy := -1; dy <= 1; dy++ {
				for dz := -1; dz <= 1; dz++ {
					ok := key{(k.x + dx + g) % g, (k.y + dy + g) % g, (k.z + dz + g) % g}
					oth := heads[ok]
					for _, i := range mem {
						for _, j := range oth {
							if j <= i {
								continue
							}
							ddx := wr(cells[j].x-cells[i].x, box)
							ddy := wr(cells[j].y-cells[i].y, box)
							ddz := wr(cells[j].z-cells[i].z, box)
							d2 := ddx*ddx + ddy*ddy + ddz*ddz
							contact := cells[i].r + cells[j].r
							if d2 >= contact*contact || d2 <= 1e-20 {
								continue
							}
							d := math.Sqrt(d2)
							m.NLinks++
							// field J
							ei := cells[i].fa1*cells[i].fa1 + cells[i].fa2*cells[i].fa2
							ej := cells[j].fa1*cells[j].fa1 + cells[j].fa2*cells[j].fa2
							if ei > 1e-14 && ej > 1e-14 {
								Jf := 2 * (cells[i].fa1*cells[j].fa2 - cells[i].fa2*cells[j].fa1)
								wf := math.Sqrt(ei * ej)
								jAbs += math.Abs(Jf) * wf
								jX += Jf * (ddx / d) * wf
								wF += wf
							}
							// dense J
							ri, rj := cells[i].em, cells[j].em
							if ri > 1e-12 && rj > 1e-12 {
								dphi := math.Atan2(math.Sin(cells[j].th2-cells[i].th2),
									math.Cos(cells[j].th2-cells[i].th2))
								Jd := 2 * math.Sqrt(ri*rj) * math.Sin(dphi)
								wd := math.Sqrt(ri * rj)
								jdAbs += math.Abs(Jd) * wd
								wD += wd
							}
						}
					}
				}
			}
		}
	}
	if wF > 0 {
		m.JAbs = jAbs / wF
		m.Jx = jX / wF
	}
	if wD > 0 {
		m.JdAbs = jdAbs / wD
	}
	return m
}

func fcsPath(runDir, name string) string {
	return filepath.Join(runDir, name+".fcs")
}
