// tab.go — v102 gather step: tables and the registered verdicts
// (../PREREG.md §5) from the gx outputs in out/.
//
//	go run tab.go            # all tables + bars + ending
//	go run tab.go s7         # the 3+2 n = 7 supplement (labelled, not scored)
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

type Run struct {
	Tag                         string
	WLock, WLockSE, WFinal      float64
	F, Nsat, Nlight             float64
	P09, PF09, FRel             float64
	NFRel                       int
	PBoth, PSpace, PLeft, PDark float64
	T, Arr, TFirst              float64
	Ev                          map[string]float64
}

type File struct {
	Runs   map[string]*Run
	Order  []string
	Static []string
	Check  string
	Start  []string
}

var rules = []string{"none", "a", "b", "c", "ab", "ac", "bc", "abc"}

func num(s string) float64 { v, _ := strconv.ParseFloat(s, 64); return v }

func rx(s, pat string) float64 {
	m := regexp.MustCompile(pat).FindStringSubmatch(s)
	if m == nil {
		return math.NaN()
	}
	return num(m[1])
}

var kv = regexp.MustCompile(`([A-Za-z_0-9/]+)=(-?[0-9.eE+-]+)`)

func parse(path string) *File {
	f, err := os.Open(path)
	if err != nil {
		return nil
	}
	defer f.Close()
	out := &File{Runs: map[string]*Run{}}
	var cur *Run
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 1<<20), 1<<20)
	for sc.Scan() {
		l := sc.Text()
		switch {
		case strings.HasPrefix(l, "RUN "):
			cur = &Run{Tag: strings.Fields(l)[1], Ev: map[string]float64{}}
			out.Runs[cur.Tag] = cur
			out.Order = append(out.Order, cur.Tag)
		case strings.HasPrefix(l, "STATIC") || strings.HasPrefix(l, "  k=") || strings.HasPrefix(l, "  maxF") || strings.HasPrefix(l, "  L0:"):
			out.Static = append(out.Static, l)
		case strings.HasPrefix(l, "start_restricted"):
			out.Start = append(out.Start, l)
		case strings.HasPrefix(l, "CHECK"):
			out.Check = l
		case cur != nil:
			t := strings.TrimSpace(l)
			fs := strings.Fields(t)
			if len(fs) == 0 {
				continue
			}
			m := map[string]float64{}
			for _, x := range kv.FindAllStringSubmatch(t, -1) {
				m[x[1]] = num(x[2])
			}
			switch fs[0] {
			case "w_lock":
				cur.WLock, cur.WLockSE, cur.WFinal = m["mean"], m["se"], m["final"]
			case "F":
				cur.F = m["mean"]
			case "nsat":
				cur.Nsat = m["mean"]
			case "Nlight":
				cur.Nlight = m["mean"]
			case "events/traj:":
				for k, v := range m {
					cur.Ev[k] = v
				}
			default:
				if strings.HasPrefix(t, "P(w_lock") {
					cur.P09 = rx(t, `^P\(w_lock>=0\.9\)=([-0-9.eE+]+)`)
					cur.PF09 = rx(t, `P\(final w_lock>=0\.9\)=([-0-9.eE+]+)`)
					cur.FRel = rx(t, `F/w_lock\(w>=0\.5\)=([-0-9.eE+]+)`)
					cur.NFRel = int(rx(t, `\(n=([0-9]+)\)`))
				} else if strings.HasPrefix(t, "P(both") {
					cur.PBoth = rx(t, `P\(both light out\)=([-0-9.eE+]+)`)
					cur.PSpace = rx(t, `P\(any space out\)=([-0-9.eE+]+)`)
					cur.PLeft = rx(t, `P\(quanta left in chain\)=([-0-9.eE+]+)`)
					cur.PDark = rx(t, `P\(ended dark-stable\)=([-0-9.eE+]+)`)
				} else if strings.HasPrefix(t, "mean t") {
					cur.T = rx(t, `mean t=([-0-9.eE+]+)`)
					cur.Arr = rx(t, `arrivals=([-0-9.eE+]+)`)
					cur.TFirst = rx(t, `t_first\(w_lock>=0\.9\)=([-0-9.eE+]+)`)
				}
			}
		}
	}
	return out
}

// gpre and spre name the grid and side-variant outputs (g/gs, or s7/s7s).
var gpre, spre = "g", "gs"

func get(fab, r string) *File { return parse(fmt.Sprintf("out/%s_%s_%s.out", gpre, fab, r)) }

type Bars struct {
	G1, G2, G3, G4, G6 bool
	G1where            string
}

func bars(f *File, r string) (b Bars) {
	if f == nil {
		return
	}
	// G1: some X1 run with mean end w_lock >= 0.5
	frel, nfr := 0.0, 0
	for _, t := range f.Order {
		if !strings.HasPrefix(t, "X1_") {
			continue
		}
		x := f.Runs[t]
		if x.WLock >= 0.5 {
			b.G1 = true
			b.G1where += strings.TrimPrefix(t, "X1_") + " "
			if x.NFRel > 0 {
				frel += x.FRel * float64(x.NFRel)
				nfr += x.NFRel
			}
		}
	}
	d, s := f.Runs["X2_dark"], f.Runs["X2_S_dense"]
	if d != nil && s != nil {
		b.G2 = d.PF09 >= 0.9 && s.WLock >= 0.5
		if d.NFRel > 0 {
			frel += d.FRel * float64(d.NFRel)
			nfr += d.NFRel
		}
	}
	if nfr > 0 {
		b.G3 = frel/float64(nfr) >= 0.9
	}
	if x := f.Runs["X3_dark"]; x != nil {
		b.G4 = x.PBoth >= 0.9
	}
	b.G6 = strings.Contains(r, "c")
	return
}

func yn(v bool) string {
	if v {
		return "✔"
	}
	return "✘"
}

func main() {
	fabs := []string{"M2", "M4"}
	names := map[string]string{"M2": "3+2 (n = 6, Λ = 4)", "M4": "W (n = 5, Λ = 3)"}
	if len(os.Args) > 1 && os.Args[1] == "s7" {
		gpre, spre = "s7", "s7s"
		fabs = []string{"M2"}
		names["M2"] = "3+2 supplement (n = 7, Λ = 5; labelled, not scored)"
	}
	for _, fab := range fabs {
		fmt.Printf("\n## %s\n", names[fab])
		fmt.Printf("\n### X4 static (lock sector)\n\n| rules | L0 residual | max-F state: F, w_lock, rate out | lowest state: w_lock, rate out |\n|---|---|---|---|\n")
		for _, r := range rules {
			f := get(fab, r)
			if f == nil {
				fmt.Printf("| %s | (missing) | | |\n", r)
				continue
			}
			var res, mf, lo string
			for _, l := range f.Static {
				m := map[string]string{}
				for _, x := range kv.FindAllStringSubmatch(l, -1) {
					m[x[1]] = x[2]
				}
				switch {
				case strings.HasPrefix(l, "  L0:"):
					res = m["residual"]
				case strings.HasPrefix(l, "  maxF"):
					mf = fmt.Sprintf("%s, %s, %s", m["F"], m["w_lock"], m["rate_out"])
				case strings.HasPrefix(l, "  k=0 "):
					lo = fmt.Sprintf("%s, %s", m["w_lock"], m["rate_out"])
				}
			}
			fmt.Printf("| %s | %s | %s | %s |\n", r, res, mf, lo)
		}
		tabs := []struct {
			title string
			tags  []string
		}{
			{"X1 gather (from the empty chain)", []string{"X1_S_dilute", "X1_S_dense", "X1_L_dilute", "X1_L_dense", "X1_SL_dilute", "X1_SL_dense"}},
			{"X2 hold (from |L0⟩)", []string{"X2_dark", "X2_S_dense", "X2_SL_dense"}},
			{"X3 light", []string{"X3_dark"}},
			{"X5 vacancy fate (1-D)", []string{"X5e_S_dilute", "X5e_S_dense", "X5b_S_dilute", "X5b_S_dense"}},
		}
		for _, tb := range tabs {
			fmt.Printf("\n### %s\n\n| rules | run | w_lock (se) | final | P(w≥0.9 end) | F/w | nsat | N_light | P(both light) | P(space out) | P(left) | arrivals | events: L out/in, S out/in |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|\n", tb.title)
			for _, r := range rules {
				f := get(fab, r)
				if f == nil {
					continue
				}
				for _, t := range tb.tags {
					x := f.Runs[t]
					if x == nil {
						continue
					}
					fmt.Printf("| %s | %s | %.4f (%.4f) | %.4f | %.3f | %.4f | %.3f | %.3f | %.3f | %.3f | %.3f | %.1f | %.1f / %.1f, %.2f / %.2f |\n",
						r, t, x.WLock, x.WLockSE, x.WFinal, x.PF09, x.FRel, x.Nsat, x.Nlight, x.PBoth, x.PSpace, x.PLeft, x.Arr,
						x.Ev["light_out"], x.Ev["light_in"], x.Ev["space_out"], x.Ev["space_in"])
				}
			}
		}
		fmt.Printf("\n### X5 side variant (c off)\n\n| rules | run | w_lock | final | P(refilled at end) | t_first(w≥0.9) | side a | side form | S in/out |\n|---|---|---|---|---|---|---|---|---|\n")
		for _, r := range []string{"none", "a", "b", "ab"} {
			f := parse(fmt.Sprintf("out/%s_%s_%s.out", spre, fab, r))
			if f == nil {
				continue
			}
			for _, t := range f.Order {
				x := f.Runs[t]
				fmt.Printf("| %s | %s | %.4f | %.4f | %.3f | %.4g | %.2f | %.2f | %.2f / %.2f |\n", r, t, x.WLock, x.WFinal, x.PF09, x.TFirst, x.Ev["side_a"], x.Ev["side_form"], x.Ev["space_in"], x.Ev["space_out"])
			}
		}
		fmt.Printf("\n### checks\n\n")
		for _, r := range rules {
			if f := get(fab, r); f != nil {
				fmt.Printf("- %s: %s\n", r, f.Check)
			}
		}
	}
	if gpre == "s7" {
		fmt.Printf("\n## Bars, 3+2 n = 7 (not scored)\n\n| rules | G1 G2 G3 G4 | G6 |\n|---|---|---|\n")
		for _, r := range rules {
			b := bars(get("M2", r), r)
			fmt.Printf("| %s | %s %s %s %s | %s |\n", r, yn(b.G1), yn(b.G2), yn(b.G3), yn(b.G4), yn(b.G6))
		}
		return
	}
	fmt.Printf("\n## Bars (PREREG §5)\n\n| rules | 3+2: G1 G2 G3 G4 | W: G1 G2 G3 G4 | G6 | gathers in both |\n|---|---|---|---|---|\n")
	ending := "N"
	anyC, anyOff, anyHold := false, false, false
	for _, r := range rules {
		b2, b4 := bars(get("M2", r), r), bars(get("M4", r), r)
		all := b2.G1 && b2.G2 && b2.G3 && b2.G4 && b4.G1 && b4.G2 && b4.G3 && b4.G4
		hold := b2.G2 && b2.G3 && b2.G4 && b4.G2 && b4.G3 && b4.G4 && !b2.G1 && !b4.G1
		fmt.Printf("| %s | %s %s %s %s | %s %s %s %s | %s | %s |\n", r, yn(b2.G1), yn(b2.G2), yn(b2.G3), yn(b2.G4), yn(b4.G1), yn(b4.G2), yn(b4.G3), yn(b4.G4), yn(b2.G6), yn(all))
		if all && b2.G6 {
			anyC = true
		}
		if all && !b2.G6 {
			anyOff = true
		}
		if hold && b2.G6 {
			anyHold = true
		}
	}
	anyG1 := false
	for _, r := range rules {
		if bars(get("M2", r), r).G1 || bars(get("M4", r), r).G1 {
			anyG1 = true
		}
	}
	switch {
	case anyC:
		ending = "Gc"
	case anyOff:
		ending = "Gl"
	case !anyG1 && anyHold:
		ending = "H"
	}
	fmt.Printf("\n**Ending: %s**\n", ending)
}
