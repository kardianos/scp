// volview — viewer for free-cell .fcs snapshots (FCS.md: v1/v2 legacy
// single-frame files, v3 chunked streams with interleaved ANLZ frames).
//
// Basis: the sfa/volview GL ray marcher (emission-absorption compositing,
// per-channel transfer function) — headless software splatter plus the
// interactive GL mode (interactive.go).
//
// Usage:
//   volview -view em -axis z -out frame.png run.fcs        headless render
//   volview -view phase -outdir imgs/ rundir/              every frame
//   volview -avg -view ee -out avg.png run.fcs             time-average
//   volview -info run.fcs                                  chunk inventory + ANLZ tables
//   volview -i run.fcs                                     interactive (orbit/scrub/channels)
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
	fa1, fa2, th2                   float32 // FCS v2+
}

type linkRec struct {
	i, j          uint32
	d, a, lem, gg float32
}

type frame struct {
	t     float64
	L     float64
	cells []cellRec
	links []linkRec
	extra map[string][]float32 // schema columns beyond the base set
}

type anlzTable struct {
	name string
	t    float64
	cols []string
	rows []float64 // nrows*ncols, row-major
}

// ------------------------------------------------------------------
// readers
// ------------------------------------------------------------------

func le32(b []byte) uint32  { return binary.LittleEndian.Uint32(b) }
func le64(b []byte) uint64  { return binary.LittleEndian.Uint64(b) }
func lef32(b []byte) float32 { return math.Float32frombits(le32(b)) }
func lef64(b []byte) float64 { return math.Float64frombits(le64(b)) }

var baseCellCols = map[string]int{"x": 0, "y": 1, "z": 2, "r": 3, "es": 4,
	"em": 5, "ee": 6, "xload": 7, "tag": 8, "fa1": 9, "fa2": 10, "th2": 11}

func parseV12(b []byte, path string, ver uint32) (*frame, error) {
	nf := 9
	if ver == 2 {
		nf = 12
	}
	n := int(le32(b[8:]))
	f := &frame{t: lef64(b[12:]), L: lef64(b[20:])}
	off := 28
	if len(b) < off+n*nf*4 {
		return nil, fmt.Errorf("%s: truncated", path)
	}
	f.cells = make([]cellRec, n)
	for i := 0; i < n; i++ {
		g := func(k int) float32 { return lef32(b[off+4*k:]) }
		c := cellRec{x: g(0), y: g(1), z: g(2), r: g(3), es: g(4), em: g(5),
			ee: g(6), xl: g(7), tag: g(8)}
		if ver == 2 {
			c.fa1, c.fa2, c.th2 = g(9), g(10), g(11)
		}
		f.cells[i] = c
		off += nf * 4
	}
	if ver == 2 && off+4 <= len(b) {
		nl := int(le32(b[off:]))
		off += 4
		if off+nl*24 <= len(b) {
			f.links = make([]linkRec, nl)
			for k := 0; k < nl; k++ {
				f.links[k] = linkRec{i: le32(b[off:]), j: le32(b[off+4:]),
					d: lef32(b[off+8:]), a: lef32(b[off+12:]),
					lem: lef32(b[off+16:]), gg: lef32(b[off+20:])}
				off += 24
			}
		}
	}
	return f, nil
}

func parseNames(b []byte, off int) ([]string, int) {
	n := int(le32(b[off:]))
	off += 4
	names := make([]string, 0, n)
	for k := 0; k < n; k++ {
		l := int(b[off])
		off++
		names = append(names, string(b[off:off+l]))
		off += l
	}
	return names, off
}

func parseV3(b []byte, path string) ([]*frame, []anlzTable, string, error) {
	var frames []*frame
	var tables []anlzTable
	var cfg strings.Builder
	var cellCols, linkCols []string
	off := 8
	for off+12 <= len(b) {
		cc := string(b[off : off+4])
		plen := int(le64(b[off+4:]))
		off += 12
		if off+plen > len(b) {
			break // truncated final chunk: stop cleanly
		}
		p := b[off : off+plen]
		off += plen
		switch cc {
		case "CFG ":
			cfg.Write(p)
		case "SCHM":
			var o int
			cellCols, o = parseNames(p, 0)
			linkCols, _ = parseNames(p, o)
		case "CELL":
			if cellCols == nil {
				return nil, nil, "", fmt.Errorf("%s: CELL before SCHM", path)
			}
			nc := len(cellCols)
			t, L := lef64(p), lef64(p[8:])
			n := int(le32(p[16:]))
			f := &frame{t: t, L: L, cells: make([]cellRec, n)}
			var extraNames []string
			for _, name := range cellCols {
				if _, known := baseCellCols[name]; !known {
					extraNames = append(extraNames, name)
				}
			}
			if len(extraNames) > 0 {
				f.extra = map[string][]float32{}
				for _, name := range extraNames {
					f.extra[name] = make([]float32, n)
				}
			}
			o := 20
			for i := 0; i < n; i++ {
				c := &f.cells[i]
				for k, name := range cellCols {
					v := lef32(p[o+4*k:])
					switch name {
					case "x":
						c.x = v
					case "y":
						c.y = v
					case "z":
						c.z = v
					case "r":
						c.r = v
					case "es":
						c.es = v
					case "em":
						c.em = v
					case "ee":
						c.ee = v
					case "xload":
						c.xl = v
					case "tag":
						c.tag = v
					case "fa1":
						c.fa1 = v
					case "fa2":
						c.fa2 = v
					case "th2":
						c.th2 = v
					default:
						f.extra[name][i] = v
					}
				}
				o += 4 * nc
			}
			frames = append(frames, f)
		case "LINK":
			if len(frames) == 0 {
				continue
			}
			f := frames[len(frames)-1]
			nl := int(le32(p[8:]))
			nfc := len(linkCols)
			rec := 8 + 4*nfc
			o := 12
			f.links = make([]linkRec, 0, nl)
			for k := 0; k < nl && o+rec <= len(p); k++ {
				lr := linkRec{i: le32(p[o:]), j: le32(p[o+4:])}
				for c, name := range linkCols {
					v := lef32(p[o+8+4*c:])
					switch name {
					case "d":
						lr.d = v
					case "A":
						lr.a = v
					case "lem":
						lr.lem = v
					case "gg":
						lr.gg = v
					}
				}
				f.links = append(f.links, lr)
				o += rec
			}
		case "ANLZ":
			if len(p) < 1 {
				continue
			}
			if p[0] == 0 { // text
				cfg.WriteString("--- ANLZ text ---\n")
				cfg.Write(p[1:])
			} else { // table
				o := 1
				nl := int(p[o])
				o++
				name := string(p[o : o+nl])
				o += nl
				t := lef64(p[o:])
				o += 8
				cols, o2 := parseNames(p, o)
				o = o2
				nrows := int(le32(p[o:]))
				o += 4
				rows := make([]float64, nrows*len(cols))
				for k := range rows {
					rows[k] = lef64(p[o+8*k:])
				}
				tables = append(tables, anlzTable{name, t, cols, rows})
			}
		default:
			// unknown chunk: skipped by length — the format's extension point
		}
	}
	return frames, tables, cfg.String(), nil
}

// loadPath reads one .fcs file of any version.
func loadPath(path string) ([]*frame, []anlzTable, string, error) {
	b, err := os.ReadFile(path)
	if err != nil {
		return nil, nil, "", err
	}
	if len(b) < 8 || string(b[:4]) != "FCS1" {
		return nil, nil, "", fmt.Errorf("%s: not an FCS file", path)
	}
	ver := le32(b[4:])
	switch ver {
	case 1, 2:
		f, err := parseV12(b, path, ver)
		if err != nil {
			return nil, nil, "", err
		}
		return []*frame{f}, nil, "", nil
	case 3:
		return parseV3(b, path)
	}
	return nil, nil, "", fmt.Errorf("%s: FCS version %d unsupported", path, ver)
}

// ------------------------------------------------------------------
// channels
// ------------------------------------------------------------------

var channelNames = []string{"es", "em", "ee", "x", "r", "phase", "fa1", "thd"}

func channelCyclic(view string) bool { return view == "phase" || view == "thd" }

func (f *frame) chanVal(i int, view string) float64 {
	c := &f.cells[i]
	switch view {
	case "es":
		return float64(c.es)
	case "em":
		return float64(c.em)
	case "ee":
		return float64(c.ee)
	case "x", "xload":
		return float64(c.xl)
	case "r":
		return float64(c.r)
	case "phase": // field phase arg(psi) in [0,1); -1 = no field (dark)
		if c.fa1*c.fa1+c.fa2*c.fa2 < 1e-12 {
			return -1
		}
		return (math.Atan2(float64(c.fa2), float64(c.fa1)) + math.Pi) / (2 * math.Pi)
	case "fa1":
		return float64(c.fa1)
	case "thd":
		return float64(c.th2) / (2 * math.Pi)
	}
	if f.extra != nil {
		if col, ok := f.extra[view]; ok {
			return float64(col[i])
		}
	}
	return float64(c.em)
}

// runChannels = base channels + any schema extras present.
func runChannels(f *frame) []string {
	out := append([]string(nil), channelNames...)
	if f.extra != nil {
		var ex []string
		for name := range f.extra {
			ex = append(ex, name)
		}
		sort.Strings(ex)
		out = append(out, ex...)
	}
	return out
}

// ------------------------------------------------------------------
// headless render
// ------------------------------------------------------------------

// heat ramp black -> deep orange -> white
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

	proj := func(c *cellRec) (u, v, w float64) {
		switch axis {
		case "x":
			return float64(c.y), float64(c.z), float64(c.x)
		case "y":
			return float64(c.x), float64(c.z), float64(c.y)
		}
		return float64(c.x), float64(c.y), float64(c.z)
	}

	vals := make([]float64, len(f.cells))
	for i := range f.cells {
		vals[i] = f.chanVal(i, view)
	}
	sorted := append([]float64(nil), vals...)
	sort.Float64s(sorted)
	vmax := 1e-12
	if len(sorted) > 0 {
		vmax = math.Max(vmax, sorted[len(sorted)*99/100])
	}
	if channelCyclic(view) {
		vmax = 1 // already normalized
	}

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
		val := vals[idx] / vmax
		if val <= 0 {
			continue
		}
		cr, cg, cb := ramp(math.Pow(val, gamma))
		if c.tag > 0 {
			cr, cg, cb = 0.2*cr, 0.5+0.5*cg, 0.5+0.5*cb
		}
		rad := float64(c.r) * scale
		cu, cv := u*scale, float64(px)-v*scale
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
				fall := 1.0 - d2/r2
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

func savePNG(img *image.RGBA, path string) {
	w, err := os.Create(path)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	if err := png.Encode(w, img); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	w.Close()
}

func buildAvgFrame(frames []*frame) *frame {
	last := frames[len(frames)-1]
	av := &frame{t: last.t, L: last.L, cells: append([]cellRec(nil), last.cells...),
		links: last.links, extra: last.extra}
	n := 0
	for i := range av.cells {
		var es, em, ee, xl, f1, f2 float64
		n = 0
		for _, f := range frames {
			if len(f.cells) != len(av.cells) {
				continue
			}
			c := &f.cells[i]
			es += float64(c.es)
			em += float64(c.em)
			ee += float64(c.ee)
			xl += float64(c.xl)
			f1 += float64(c.fa1)
			f2 += float64(c.fa2)
			n++
		}
		if n == 0 {
			break
		}
		inv := 1.0 / float64(n)
		av.cells[i].es = float32(es * inv)
		av.cells[i].em = float32(em * inv)
		av.cells[i].ee = float32(ee * inv)
		av.cells[i].xl = float32(xl * inv)
		av.cells[i].fa1 = float32(f1 * inv)
		av.cells[i].fa2 = float32(f2 * inv)
	}
	return av
}

// ------------------------------------------------------------------

func main() {
	view := flag.String("view", "em", "channel: es|em|ee|x|r|phase|fa1|thd|<schema name>")
	axis := flag.String("axis", "z", "projection axis: x|y|z")
	px := flag.Int("px", 768, "image size")
	gamma := flag.Float64("gamma", 0.7, "transfer gamma")
	out := flag.String("out", "", "output png (single frame)")
	outdir := flag.String("outdir", ".", "output dir (multi-frame)")
	avg := flag.Bool("avg", false, "time-average the channel over all frames")
	interactive := flag.Bool("i", false, "interactive viewer: orbit, time scrub, channels, links (needs a display)")
	selfcheck := flag.Bool("selfcheck", false, "with -i: render a few frames, screenshot, exit")
	info := flag.Bool("info", false, "print stream inventory, CFG provenance, and ANLZ tables")
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

	if *info {
		for _, fp := range files {
			frames, tables, cfg, err := loadPath(fp)
			if err != nil {
				fmt.Fprintln(os.Stderr, err)
				os.Exit(1)
			}
			fmt.Printf("== %s: %d frame(s)\n", fp, len(frames))
			for k, f := range frames {
				fmt.Printf("   frame %d t=%.3f L=%g cells=%d links=%d channels=%s\n",
					k, f.t, f.L, len(f.cells), len(f.links), strings.Join(runChannels(f), ","))
			}
			if cfg != "" {
				fmt.Println("-- CFG/provenance:")
				fmt.Print(cfg)
				if !strings.HasSuffix(cfg, "\n") {
					fmt.Println()
				}
			}
			for _, tb := range tables {
				fmt.Printf("-- ANLZ %s t=%.3f (%d rows)\n", tb.name, tb.t, len(tb.rows)/len(tb.cols))
				fmt.Println("   " + strings.Join(tb.cols, "\t"))
				nr := len(tb.rows) / len(tb.cols)
				for r := 0; r < nr; r++ {
					var sb strings.Builder
					for c := range tb.cols {
						if c > 0 {
							sb.WriteByte('\t')
						}
						fmt.Fprintf(&sb, "%.6g", tb.rows[r*len(tb.cols)+c])
					}
					fmt.Println("   " + sb.String())
				}
			}
		}
		return
	}

	// load all frames (a dir of files and/or multi-frame streams)
	var frames []*frame
	for _, fp := range files {
		fr, _, _, err := loadPath(fp)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
		frames = append(frames, fr...)
	}
	if len(frames) == 0 {
		fmt.Fprintln(os.Stderr, "no frames")
		os.Exit(1)
	}
	sort.SliceStable(frames, func(a, b int) bool { return frames[a].t < frames[b].t })

	if *interactive {
		runInteractive(frames, *selfcheck)
		return
	}
	if *avg && len(frames) > 1 {
		av := buildAvgFrame(frames)
		img := render(av, *view, *axis, *px, *gamma)
		op := *out
		if op == "" {
			op = filepath.Join(*outdir, fmt.Sprintf("avg_%dframes_%s.png", len(frames), *view))
		}
		savePNG(img, op)
		fmt.Printf("# volview avg of %d frames -> %s\n", len(frames), op)
		return
	}
	if len(frames) == 1 || *out != "" {
		f := frames[0]
		img := render(f, *view, *axis, *px, *gamma)
		op := *out
		if op == "" {
			op = filepath.Join(*outdir, fmt.Sprintf("frame_t%.3f_%s.png", f.t, *view))
		}
		savePNG(img, op)
		fmt.Printf("# volview t=%.3f cells=%d -> %s\n", f.t, len(f.cells), op)
		return
	}
	for _, f := range frames {
		img := render(f, *view, *axis, *px, *gamma)
		op := filepath.Join(*outdir, fmt.Sprintf("frame_t%08.3f_%s.png", f.t, *view))
		savePNG(img, op)
		fmt.Printf("# volview t=%.3f cells=%d -> %s\n", f.t, len(f.cells), op)
	}
}
