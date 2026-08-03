// volview — volumetric viewer for free-cell .fcs snapshots (headless).
//
// Basis: the sfa/volview GL ray marcher (emission-absorption compositing,
// per-channel transfer function), reduced to a software orthographic
// splatter for the free-cell point clouds. The interactive GL front-end
// is the P1+ enhancement; this renders committed evidence images now.
//
// Usage:
//   volview -view em -axis z -px 768 -out frame.png snap_00000_t0.000.fcs
//   volview -view ee -all -outdir imgs/ rundir/          (every .fcs in dir)
//
// Channels: es (space), em (dense), ee (field), x (load), r (radius).
// Tagged cells are tinted cyan so seeded objects stand out.
package main

import (
	"encoding/binary"
	"flag"
	"fmt"
	"image"
	"image/color"
	"image/png"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"
)

type cellRec struct {
	x, y, z, r, es, em, ee, xl, tag float32
}

type frame struct {
	t     float64
	L     float64
	cells []cellRec
}

func readFCS(path string) (*frame, error) {
	b, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	if len(b) < 28 || string(b[:4]) != "FCS1" {
		return nil, fmt.Errorf("%s: not an FCS1 file", path)
	}
	ver := binary.LittleEndian.Uint32(b[4:])
	if ver != 1 {
		return nil, fmt.Errorf("%s: FCS version %d unsupported", path, ver)
	}
	n := int(binary.LittleEndian.Uint32(b[8:]))
	f := &frame{
		t: math.Float64frombits(binary.LittleEndian.Uint64(b[12:])),
		L: math.Float64frombits(binary.LittleEndian.Uint64(b[20:])),
	}
	off := 28
	need := off + n*9*4
	if len(b) < need {
		return nil, fmt.Errorf("%s: truncated (%d < %d)", path, len(b), need)
	}
	f.cells = make([]cellRec, n)
	for i := 0; i < n; i++ {
		g := func(k int) float32 {
			return math.Float32frombits(binary.LittleEndian.Uint32(b[off+4*k:]))
		}
		f.cells[i] = cellRec{g(0), g(1), g(2), g(3), g(4), g(5), g(6), g(7), g(8)}
		off += 9 * 4
	}
	return f, nil
}

func channel(c *cellRec, view string) float64 {
	switch view {
	case "es":
		return float64(c.es)
	case "em":
		return float64(c.em)
	case "ee":
		return float64(c.ee)
	case "x":
		return float64(c.xl)
	case "r":
		return float64(c.r)
	}
	return float64(c.em)
}

// heat ramp black -> deep orange -> white (volview transfer-function
// spirit: monotone luminance, alpha follows value)
func ramp(v float64) (r, g, b float64) {
	if v < 0 {
		v = 0
	}
	if v > 1 {
		v = 1
	}
	r = math.Min(1, 2.2*v)
	g = math.Max(0, math.Min(1, 1.6*v-0.25))
	b = math.Max(0, math.Min(1, 2.5*v-1.5))
	return
}

func render(f *frame, view, axis string, px int, gamma float64) *image.RGBA {
	img := image.NewRGBA(image.Rect(0, 0, px, px))
	scale := float64(px) / f.L

	// project: u,v = in-plane coords, w = depth
	proj := func(c *cellRec) (u, v, w float64) {
		switch axis {
		case "x":
			return float64(c.y), float64(c.z), float64(c.x)
		case "y":
			return float64(c.x), float64(c.z), float64(c.y)
		}
		return float64(c.x), float64(c.y), float64(c.z)
	}

	// normalize channel to its 99th percentile (robust vmax)
	vals := make([]float64, 0, len(f.cells))
	for i := range f.cells {
		vals = append(vals, channel(&f.cells[i], view))
	}
	sorted := append([]float64(nil), vals...)
	sort.Float64s(sorted)
	vmax := 1e-12
	if len(sorted) > 0 {
		vmax = math.Max(vmax, sorted[len(sorted)*99/100])
	}

	// back-to-front by depth, emission-absorption composite
	order := make([]int, len(f.cells))
	for i := range order {
		order[i] = i
	}
	sort.Slice(order, func(a, b int) bool {
		_, _, wa := proj(&f.cells[order[a]])
		_, _, wb := proj(&f.cells[order[b]])
		return wa < wb
	})

	buf := make([]float64, px*px*3)
	for _, idx := range order {
		c := &f.cells[idx]
		u, v, _ := proj(c)
		val := channel(c, view) / vmax
		if val <= 0 {
			continue
		}
		cr, cg, cb := ramp(math.Pow(val, gamma))
		if c.tag > 0 { // tinted so seeded objects stand out
			cr, cg, cb = 0.2*cr, 0.5+0.5*cg, 0.5+0.5*cb
		}
		rad := float64(c.r) * scale
		cu, cv := u*scale, float64(px)-v*scale // y up
		alpha := math.Min(0.85, 0.15+0.7*val)
		r2 := rad * rad
		lo := func(x float64) int { return int(math.Max(0, math.Floor(x))) }
		hi := func(x float64) int { return int(math.Min(float64(px-1), math.Ceil(x))) }
		for yy := lo(cv - rad); yy <= hi(cv+rad); yy++ {
			for xx := lo(cu - rad); xx <= hi(cu+rad); xx++ {
				du, dv := float64(xx)+0.5-cu, float64(yy)+0.5-cv
				d2 := du*du + dv*dv
				if d2 > r2 {
					continue
				}
				fall := 1.0 - d2/r2 // soft sphere profile
				a := alpha * fall
				o := (yy*px + xx) * 3
				buf[o] = buf[o]*(1-a) + cr*a
				buf[o+1] = buf[o+1]*(1-a) + cg*a
				buf[o+2] = buf[o+2]*(1-a) + cb*a
			}
		}
	}
	for y := 0; y < px; y++ {
		for x := 0; x < px; x++ {
			o := (y*px + x) * 3
			img.SetRGBA(x, y, color.RGBA{
				uint8(255 * math.Min(1, buf[o])),
				uint8(255 * math.Min(1, buf[o+1])),
				uint8(255 * math.Min(1, buf[o+2])), 255})
		}
	}
	return img
}

func main() {
	view := flag.String("view", "em", "channel: es|em|ee|x|r")
	axis := flag.String("axis", "z", "projection axis: x|y|z")
	px := flag.Int("px", 768, "image size")
	gamma := flag.Float64("gamma", 0.7, "transfer gamma")
	out := flag.String("out", "", "output png (single input)")
	outdir := flag.String("outdir", ".", "output dir (directory input)")
	avg := flag.Bool("avg", false, "time-average the channel over all frames (per cell index, positions from the last frame) — the shell/response-structure view")
	flag.Parse()
	if flag.NArg() < 1 {
		fmt.Fprintln(os.Stderr, "usage: volview [flags] file.fcs|dir")
		os.Exit(2)
	}
	arg := flag.Arg(0)
	st, err := os.Stat(arg)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	var files []string
	if st.IsDir() {
		ents, _ := os.ReadDir(arg)
		for _, e := range ents {
			if strings.HasSuffix(e.Name(), ".fcs") {
				files = append(files, filepath.Join(arg, e.Name()))
			}
		}
		sort.Strings(files)
	} else {
		files = []string{arg}
	}
	if *avg && len(files) > 1 {
		var acc *frame
		nf := 0
		for _, fp := range files {
			f, err := readFCS(fp)
			if err != nil {
				fmt.Fprintln(os.Stderr, err)
				os.Exit(1)
			}
			if acc == nil {
				acc = f
			} else if len(f.cells) == len(acc.cells) {
				for i := range f.cells {
					f.cells[i].es += acc.cells[i].es
					f.cells[i].em += acc.cells[i].em
					f.cells[i].ee += acc.cells[i].ee
					f.cells[i].xl += acc.cells[i].xl
				}
				acc = f // positions/r/tag ride the newest frame
			} else {
				fmt.Fprintf(os.Stderr, "skip %s: cell count changed\n", fp)
				continue
			}
			nf++
		}
		inv := float32(1.0 / float64(nf))
		for i := range acc.cells {
			acc.cells[i].es *= inv
			acc.cells[i].em *= inv
			acc.cells[i].ee *= inv
			acc.cells[i].xl *= inv
		}
		img := render(acc, *view, *axis, *px, *gamma)
		op := *out
		if op == "" {
			op = filepath.Join(*outdir, fmt.Sprintf("avg_%dframes_%s.png", nf, *view))
		}
		w, err := os.Create(op)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
		if err := png.Encode(w, img); err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
		w.Close()
		fmt.Printf("# volview avg of %d frames -> %s\n", nf, op)
		return
	}
	for _, fp := range files {
		f, err := readFCS(fp)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
		img := render(f, *view, *axis, *px, *gamma)
		op := *out
		if op == "" || len(files) > 1 {
			base := strings.TrimSuffix(filepath.Base(fp), ".fcs")
			op = filepath.Join(*outdir, fmt.Sprintf("%s_%s.png", base, *view))
		}
		w, err := os.Create(op)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
		if err := png.Encode(w, img); err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
		w.Close()
		fmt.Printf("# volview %s t=%.3f cells=%d -> %s\n", fp, f.t, len(f.cells), op)
	}
}
