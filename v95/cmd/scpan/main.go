// scpan — SCP analysis: structured parsing of fabrun/freecell simulation output.
//
// Replaces fragile awk/sed/rg pipelines (which silently break when column
// indices shift) by parsing the diag/FCS HEADERS to build a column-name →
// value map, then answering queries by NAME.
//
// Subcommands:
//
//	scpan diag    <logfile>     extract named diag columns over time
//	scpan result  <logfile>     extract RESULT key=value pairs
//	scpan ring    <fcsdump.txt> ring winding analysis (port of analyze_ring.py)
//	scpan cor     <fcsdump.txt> spatial correlation between two FCS columns
//	scpan profile <fcsdump.txt> spatial profile of one FCS column along an axis
//
// The FCS subcommands consume the TEXT output of `fcsdump -mode cells`
// (pass a text file, or `-` to read from stdin).
package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
)

func main() {
	if len(os.Args) < 2 {
		usage()
		os.Exit(2)
	}
	cmd := os.Args[1]
	args := os.Args[2:]
	var err error
	switch cmd {
	case "diag":
		err = runDiag(args)
	case "result":
		err = runResult(args)
	case "ring":
		err = runRing(args)
	case "cor":
		err = runCor(args)
	case "profile":
		err = runProfile(args)
	case "-h", "--help", "help":
		usage()
		return
	default:
		fmt.Fprintf(os.Stderr, "unknown subcommand %q\n", cmd)
		usage()
		os.Exit(2)
	}
	if err != nil {
		fmt.Fprintln(os.Stderr, "error:", err)
		os.Exit(1)
	}
}

func usage() {
	fmt.Fprint(os.Stderr, "scpan — structured SCP simulation analysis\n\n")
	fmt.Fprintln(os.Stderr, "Usage:")
	fmt.Fprintln(os.Stderr, "  scpan diag    <logfile>      [flags]   named diag columns over time")
	fmt.Fprintln(os.Stderr, "  scpan result  <logfile>      [flags]   RESULT key=value pairs")
	fmt.Fprintln(os.Stderr, "  scpan ring    <fcsdump.txt>  [flags]   ring winding (port of analyze_ring.py)")
	fmt.Fprintln(os.Stderr, "  scpan cor     <fcsdump.txt>  [flags]   spatial correlation of two columns")
	fmt.Fprintln(os.Stderr, "  scpan profile <fcsdump.txt>  [flags]   spatial profile of one column")
	fmt.Fprintln(os.Stderr)
	fmt.Fprintln(os.Stderr, "ring/cor/profile read the TEXT output of `fcsdump -mode cells`;")
	fmt.Fprintln(os.Stderr, "pass a text file or `-` for stdin.")
	fmt.Fprintln(os.Stderr)
	fmt.Fprintln(os.Stderr, "Each subcommand has a --help flag.")
}

// ---------------------------------------------------------------------------
// diag: named columns from the pipe-delimited diag rows
// ---------------------------------------------------------------------------

// parseDiagHeader turns a header line like
//
//	# t | drift_rel | phi_nom z_live NLa NLd births deaths | dbar sigma_d% |
//	  maxstep maxbond | Em_tag ret com_dev rms conn z_tag | pair: d delta gf gb |
//	  core/far | Em_tot
//
// into an ordered list of positional column names. Unit suffixes (the '%' in
// sigma_d%) and the 'pair:' label are stripped.
func parseDiagHeader(line string) []string {
	s := strings.TrimLeft(line, "# \t")
	var cols []string
	for _, sec := range strings.Split(s, "|") {
		sec = strings.TrimSpace(sec)
		sec = strings.TrimPrefix(sec, "pair:")
		sec = strings.TrimSpace(sec)
		if sec == "" {
			continue
		}
		for _, name := range strings.Fields(sec) {
			name = strings.TrimSuffix(name, "%")
			cols = append(cols, name)
		}
	}
	return cols
}

// parseDiagRow turns a data row like
//
//	t=    0.00 +0.00e+00 | 1.6733 15.29 20271 0 0 0 | ... | Em_tot=4.8750
//
// into a column-name → value map. Numeric bare tokens map positionally onto
// the header columns (in order); key=value tokens map by key; other bare
// tokens (section labels like 'edge_dev') are skipped. Returns the row map,
// the t value, and whether the line was a diag row.
func parseDiagRow(line string, cols []string) (map[string]string, float64, bool) {
	line = strings.TrimSpace(line)
	if !strings.HasPrefix(line, "t=") {
		return nil, 0, false
	}
	line = strings.TrimPrefix(line, "t=")
	line = strings.ReplaceAll(line, "|", " ")
	m := make(map[string]string)
	var positional []string
	for _, tok := range strings.Fields(line) {
		if eq := strings.IndexByte(tok, '='); eq >= 0 {
			m[tok[:eq]] = tok[eq+1:]
			continue
		}
		if _, err := strconv.ParseFloat(tok, 64); err == nil {
			positional = append(positional, tok)
		}
		// non-numeric bare token → section label (e.g. edge_dev); skip
	}
	for i, v := range positional {
		if i < len(cols) {
			m[cols[i]] = v
		}
	}
	var t float64
	if len(positional) > 0 {
		t, _ = strconv.ParseFloat(positional[0], 64)
	}
	return m, t, true
}

func runDiag(args []string) error {
	fs := flag.NewFlagSet("diag", flag.ExitOnError)
	fs.Usage = func() {
		fmt.Fprintln(os.Stderr, "usage: scpan diag <logfile> [--col a,b] [--from T] [--to T] [--last N] [--stats]")
		fs.PrintDefaults()
	}
	col := fs.String("col", "", "comma-separated column names (default: all present)")
	from := fs.Float64("from", 0, "start time")
	to := fs.Float64("to", math.Inf(1), "end time")
	last := fs.Int("last", 0, "last N rows only (0 = all)")
	stats := fs.Bool("stats", false, "print summary stats instead of time series")
	if err := fs.Parse(args); err != nil {
		return err
	}
	logfile := fs.Arg(0)
	if logfile == "" {
		return fmt.Errorf("diag: logfile required")
	}
	f, err := os.Open(logfile)
	if err != nil {
		return err
	}
	defer f.Close()

	type row struct {
		t float64
		m map[string]string
	}
	var cols []string
	var rows []row
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 1<<20), 1<<20)
	for sc.Scan() {
		line := sc.Text()
		trim := strings.TrimSpace(line)
		if cols == nil && strings.HasPrefix(trim, "#") && strings.Contains(line, "drift_rel") && strings.Contains(line, "|") {
			cols = parseDiagHeader(line)
			continue
		}
		if !strings.HasPrefix(trim, "t=") {
			continue
		}
		if cols == nil {
			cols = []string{"t", "drift_rel"} // header missing: minimal skeleton
		}
		m, t, ok := parseDiagRow(line, cols)
		if !ok {
			continue
		}
		if t < *from || t > *to {
			continue
		}
		rows = append(rows, row{t, m})
	}
	if err := sc.Err(); err != nil {
		return err
	}
	if cols == nil {
		return fmt.Errorf("diag: no header or t= rows found in %s", logfile)
	}

	// determine output columns
	var outCols []string
	if *col != "" {
		for _, c := range strings.Split(*col, ",") {
			outCols = append(outCols, strings.TrimSpace(c))
		}
	} else {
		// default: positional header columns (in order) then any extra kv keys,
		// dropping any that never appear in a selected row.
		seen := map[string]bool{}
		for _, c := range cols {
			seen[c] = true
			outCols = append(outCols, c)
		}
		var extra []string
		for _, r := range rows {
			keys := make([]string, 0, len(r.m))
			for k := range r.m {
				keys = append(keys, k)
			}
			sort.Strings(keys)
			for _, k := range keys {
				if !seen[k] {
					seen[k] = true
					extra = append(extra, k)
				}
			}
		}
		outCols = append(outCols, extra...)
		present := map[string]bool{}
		for _, r := range rows {
			for k := range r.m {
				present[k] = true
			}
		}
		filt := outCols[:0]
		for _, c := range outCols {
			if present[c] {
				filt = append(filt, c)
			}
		}
		outCols = filt
	}

	if *last > 0 && len(rows) > *last {
		rows = rows[len(rows)-*last:]
	}

	// warn on explicitly-requested columns never seen
	if *col != "" {
		present := map[string]bool{}
		for _, r := range rows {
			for k := range r.m {
				present[k] = true
			}
		}
		for _, c := range outCols {
			if !present[c] {
				fmt.Fprintf(os.Stderr, "warning: column %q not present in any row\n", c)
			}
		}
	}

	w := bufio.NewWriter(os.Stdout)
	defer w.Flush()

	if *stats {
		fmt.Fprintln(w, "col\tmin\tmax\tmean\tmedian\tn")
		for _, c := range outCols {
			vals := []float64{}
			for _, r := range rows {
				if s, ok := r.m[c]; ok {
					if v, err := strconv.ParseFloat(s, 64); err == nil {
						vals = append(vals, v)
					}
				}
			}
			if len(vals) == 0 {
				fmt.Fprintf(w, "%s\t\t\t\t\t0\n", c)
				continue
			}
			sorted := append([]float64(nil), vals...)
			sort.Float64s(sorted)
			var sum float64
			for _, v := range vals {
				sum += v
			}
			mean := sum / float64(len(vals))
			n := len(sorted)
			med := sorted[n/2]
			if n%2 == 0 {
				med = (sorted[n/2-1] + sorted[n/2]) / 2
			}
			fmt.Fprintf(w, "%s\t%g\t%g\t%g\t%g\t%d\n", c, sorted[0], sorted[len(sorted)-1], mean, med, n)
		}
		return nil
	}

	header := append([]string{"t"}, outCols...)
	fmt.Fprintln(w, strings.Join(header, "\t"))
	for _, r := range rows {
		parts := []string{strconv.FormatFloat(r.t, 'f', -1, 64)}
		for _, c := range outCols {
			parts = append(parts, r.m[c])
		}
		fmt.Fprintln(w, strings.Join(parts, "\t"))
	}
	return nil
}

// ---------------------------------------------------------------------------
// result: RESULT key=value pairs
// ---------------------------------------------------------------------------

func runResult(args []string) error {
	fs := flag.NewFlagSet("result", flag.ExitOnError)
	fs.Usage = func() {
		fmt.Fprintln(os.Stderr, "usage: scpan result <logfile> [--key a,b]")
		fs.PrintDefaults()
	}
	key := fs.String("key", "", "comma-separated keys to extract (default: all)")
	if err := fs.Parse(args); err != nil {
		return err
	}
	logfile := fs.Arg(0)
	if logfile == "" {
		return fmt.Errorf("result: logfile required")
	}
	f, err := os.Open(logfile)
	if err != nil {
		return err
	}
	defer f.Close()

	var order []string
	kv := map[string]string{}
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 1<<20), 1<<20)
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if !strings.HasPrefix(line, "# RESULT") {
			continue
		}
		toks := strings.Fields(strings.TrimSpace(strings.TrimPrefix(line, "# RESULT")))
		i := 0
		for i < len(toks) {
			tok := toks[i]
			if eq := strings.IndexByte(tok, '='); eq >= 0 {
				k, v := tok[:eq], tok[eq+1:]
				if _, ok := kv[k]; !ok {
					order = append(order, k)
				}
				kv[k] = v
				i++
				continue
			}
			// bare token followed by another bare token → key/value pair
			if i+1 < len(toks) && !strings.Contains(toks[i+1], "=") {
				k, v := tok, toks[i+1]
				if _, ok := kv[k]; !ok {
					order = append(order, k)
				}
				kv[k] = v
				i += 2
				continue
			}
			// standalone label/kind (e.g. 'conv', 'truss'); skip
			i++
		}
	}
	if err := sc.Err(); err != nil {
		return err
	}

	want := order
	if *key != "" {
		want = nil
		for _, k := range strings.Split(*key, ",") {
			want = append(want, strings.TrimSpace(k))
		}
	}
	w := bufio.NewWriter(os.Stdout)
	defer w.Flush()
	for _, k := range want {
		if v, ok := kv[k]; ok {
			fmt.Fprintf(w, "%s=%s\n", k, v)
		} else {
			fmt.Fprintf(os.Stderr, "warning: key %q not found\n", k)
		}
	}
	return nil
}

// ---------------------------------------------------------------------------
// FCS text dump reading (shared by ring/cor/profile)
// ---------------------------------------------------------------------------

type fcsFrame struct {
	t  float64
	rw [][]float64
}

// readFCS reads the text output of `fcsdump -mode cells`, parsing the header
// for column names, and groups rows by frame time. Returns frames sorted by t
// and the column-name → index map.
func readFCS(r io.Reader) ([]fcsFrame, map[string]int, error) {
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 1<<20), 1<<24)
	if !sc.Scan() {
		return nil, nil, fmt.Errorf("fcs: empty input")
	}
	hdr := strings.Fields(sc.Text())
	idx := make(map[string]int, len(hdr))
	for i, name := range hdr {
		idx[name] = i
	}
	tIdx, ok := idx["t"]
	if !ok {
		return nil, nil, fmt.Errorf("fcs: header has no 't' column: %q", sc.Text())
	}
	framesMap := map[float64]*fcsFrame{}
	var order []float64
	for sc.Scan() {
		toks := strings.Fields(sc.Text())
		if len(toks) <= tIdx {
			continue
		}
		t, err := strconv.ParseFloat(toks[tIdx], 64)
		if err != nil {
			continue
		}
		row := make([]float64, len(toks))
		ok := true
		for i, tk := range toks {
			v, err := strconv.ParseFloat(tk, 64)
			if err != nil {
				ok = false
				break
			}
			row[i] = v
		}
		if !ok {
			continue
		}
		fr, exists := framesMap[t]
		if !exists {
			fr = &fcsFrame{t: t}
			framesMap[t] = fr
			order = append(order, t)
		}
		fr.rw = append(fr.rw, row)
	}
	if err := sc.Err(); err != nil {
		return nil, nil, err
	}
	sort.Float64s(order)
	frames := make([]fcsFrame, 0, len(order))
	for _, t := range order {
		frames = append(frames, *framesMap[t])
	}
	return frames, idx, nil
}

func openFCSArg(arg string) (io.ReadCloser, error) {
	if arg == "-" || arg == "" {
		return io.NopCloser(bufio.NewReader(os.Stdin)), nil
	}
	return os.Open(arg)
}

// resolveFrame returns the frames to process given --frame = "last" | "all" | <t>
func resolveFrame(frames []fcsFrame, spec string) ([]fcsFrame, error) {
	switch spec {
	case "", "last":
		if len(frames) == 0 {
			return nil, fmt.Errorf("no frames")
		}
		return frames[len(frames)-1:], nil
	case "all":
		return frames, nil
	default:
		tv, err := strconv.ParseFloat(spec, 64)
		if err != nil {
			return nil, fmt.Errorf("invalid --frame %q (use last|all|<time>)", spec)
		}
		var best *fcsFrame
		bestD := math.Inf(1)
		for i := range frames {
			if d := math.Abs(frames[i].t - tv); d < bestD {
				bestD = d
				best = &frames[i]
			}
		}
		if best == nil {
			return nil, fmt.Errorf("no frame near t=%s", spec)
		}
		return []fcsFrame{*best}, nil
	}
}

func colVal(idx map[string]int, row []float64, name string) (float64, error) {
	i, ok := idx[name]
	if !ok {
		return 0, fmt.Errorf("fcs: column %q not in header", name)
	}
	if i >= len(row) {
		return 0, fmt.Errorf("fcs: column %q out of range", name)
	}
	return row[i], nil
}

// ---------------------------------------------------------------------------
// ring: winding analysis (port of report/analyze_ring.py)
// ---------------------------------------------------------------------------

func runRing(args []string) error {
	fs := flag.NewFlagSet("ring", flag.ExitOnError)
	fs.Usage = func() {
		fmt.Fprintln(os.Stderr, "usage: scpan ring <fcsdump.txt> [--m 2] [--cx 8] [--cy 8] [--frame last]")
		fs.PrintDefaults()
	}
	m := fs.Int("m", 2, "seeded winding number")
	cx := fs.Float64("cx", 8, "box center x")
	cy := fs.Float64("cy", 8, "box center y")
	frame := fs.String("frame", "last", "frame: last | all | <time>")
	if err := fs.Parse(args); err != nil {
		return err
	}
	in, err := openFCSArg(fs.Arg(0))
	if err != nil {
		return err
	}
	defer in.Close()
	frames, idx, err := readFCS(in)
	if err != nil {
		return err
	}
	sel, err := resolveFrame(frames, *frame)
	if err != nil {
		return err
	}
	for _, n := range []string{"x", "y", "em", "tag", "th2"} {
		if _, ok := idx[n]; !ok {
			return fmt.Errorf("fcs: required column %q missing", n)
		}
	}

	w := bufio.NewWriter(os.Stdout)
	defer w.Flush()
	fmt.Fprintln(w, "t\tW\tR_m\tEm_min\tEm_mean\tkeep")
	for _, fr := range sel {
		var rx, ry, rem, rth []float64
		var emAll, emRing float64
		for _, row := range fr.rw {
			x, _ := colVal(idx, row, "x")
			y, _ := colVal(idx, row, "y")
			em, _ := colVal(idx, row, "em")
			tg, _ := colVal(idx, row, "tag")
			th, _ := colVal(idx, row, "th2")
			emAll += em
			if tg > 0.5 {
				rx = append(rx, x)
				ry = append(ry, y)
				rem = append(rem, em)
				rth = append(rth, th)
				emRing += em
			}
		}
		n := len(rx)
		if n < 3 {
			fmt.Fprintf(w, "%.2f\tring lost (n=%d)\n", fr.t, n)
			continue
		}
		phi := make([]float64, n)
		for i := range rx {
			phi[i] = math.Atan2(ry[i]-*cy, rx[i]-*cx)
		}
		ord := argsortFloats(phi)
		tr := make([]float64, n)
		for i, j := range ord {
			tr[i] = rth[j]
		}
		var dsum float64
		for i := 0; i < n; i++ {
			dsum += wrapPi(tr[(i+1)%n] - tr[i])
		}
		W := dsum / (2 * math.Pi)

		var re, im float64
		for i := 0; i < n; i++ {
			a := rth[i] - float64(*m)*phi[i]
			re += math.Cos(a)
			im += math.Sin(a)
		}
		re /= float64(n)
		im /= float64(n)
		Rm := math.Sqrt(re*re + im*im)

		emMin := rem[0]
		var emMean float64
		for _, v := range rem {
			if v < emMin {
				emMin = v
			}
			emMean += v
		}
		emMean /= float64(n)

		var keep float64
		if emAll > 0 {
			keep = emRing / emAll
		}
		fmt.Fprintf(w, "%.2f\t%+.3f\t%.3f\t%.4f\t%.4f\t%.4f\n",
			fr.t, W, Rm, emMin, emMean, keep)
	}
	return nil
}

// wrapPi wraps an angle difference into [-pi, pi) with positive-divisor modulo
// semantics, matching numpy's `(d + pi) % (2pi) - pi`.
func wrapPi(d float64) float64 {
	d = math.Mod(d+math.Pi, 2*math.Pi)
	if d < 0 {
		d += 2 * math.Pi
	}
	return d - math.Pi
}

func argsortFloats(a []float64) []int {
	idx := make([]int, len(a))
	for i := range idx {
		idx[i] = i
	}
	sort.SliceStable(idx, func(i, j int) bool { return a[idx[i]] < a[idx[j]] })
	return idx
}

// ---------------------------------------------------------------------------
// cor: spatial correlation between two FCS columns
// ---------------------------------------------------------------------------

func runCor(args []string) error {
	fs := flag.NewFlagSet("cor", flag.ExitOnError)
	fs.Usage = func() {
		fmt.Fprintln(os.Stderr, "usage: scpan cor <fcsdump.txt> [--colA em] [--colB ee] [--x 8] [--z 8] [--dx 2] [--dz 2] [--axis y] [--bw 1]")
		fs.PrintDefaults()
	}
	colA := fs.String("colA", "em", "first column")
	colB := fs.String("colB", "ee", "second column")
	cx := fs.Float64("x", 8, "region center x")
	cz := fs.Float64("z", 8, "region center z")
	dx := fs.Float64("dx", 2, "region half-width in x")
	dz := fs.Float64("dz", 2, "region half-width in z")
	axis := fs.String("axis", "y", "axis to bin along (x|y|z)")
	bw := fs.Float64("bw", 1, "bin width along the axis")
	frame := fs.String("frame", "last", "frame: last | all | <time>")
	if err := fs.Parse(args); err != nil {
		return err
	}
	in, err := openFCSArg(fs.Arg(0))
	if err != nil {
		return err
	}
	defer in.Close()
	frames, idx, err := readFCS(in)
	if err != nil {
		return err
	}
	sel, err := resolveFrame(frames, *frame)
	if err != nil {
		return err
	}
	for _, n := range []string{"x", "z", *axis, *colA, *colB} {
		if _, ok := idx[n]; !ok {
			return fmt.Errorf("fcs: column %q missing", n)
		}
	}
	axI := idx[*axis]

	w := bufio.NewWriter(os.Stdout)
	defer w.Flush()
	fmt.Fprintf(w, "t\tpearson_r\tn_bins\tmean_%s\tmean_%s\n", *colA, *colB)
	for _, fr := range sel {
		type acc struct {
			sa, sb float64
			n      int
		}
		bins := map[int]*acc{}
		for _, row := range fr.rw {
			x, _ := colVal(idx, row, "x")
			z, _ := colVal(idx, row, "z")
			if math.Abs(x-*cx) > *dx || math.Abs(z-*cz) > *dz {
				continue
			}
			bk := int(math.Round(row[axI] / *bw))
			a, _ := colVal(idx, row, *colA)
			b, _ := colVal(idx, row, *colB)
			bn := bins[bk]
			if bn == nil {
				bn = &acc{}
				bins[bk] = bn
			}
			bn.sa += a
			bn.sb += b
			bn.n++
		}
		keys := make([]int, 0, len(bins))
		for k := range bins {
			keys = append(keys, k)
		}
		sort.Ints(keys)
		if len(keys) < 2 {
			fmt.Fprintf(w, "%.2f\tnan\t%d\t\t\n", fr.t, len(keys))
			continue
		}
		xv := make([]float64, len(keys))
		yv := make([]float64, len(keys))
		var sumA, sumB float64
		var cnt int
		for i, k := range keys {
			bn := bins[k]
			xv[i] = bn.sa / float64(bn.n)
			yv[i] = bn.sb / float64(bn.n)
			sumA += bn.sa
			sumB += bn.sb
			cnt += bn.n
		}
		r := pearson(xv, yv)
		var mA, mB float64
		if cnt > 0 {
			mA = sumA / float64(cnt)
			mB = sumB / float64(cnt)
		}
		fmt.Fprintf(w, "%.2f\t%.4f\t%d\t%.4f\t%.4f\n", fr.t, r, len(keys), mA, mB)
	}
	return nil
}

func pearson(x, y []float64) float64 {
	n := len(x)
	if n == 0 {
		return math.NaN()
	}
	var sx, sy float64
	for i := 0; i < n; i++ {
		sx += x[i]
		sy += y[i]
	}
	mx, my := sx/float64(n), sy/float64(n)
	var num, dx2, dy2 float64
	for i := 0; i < n; i++ {
		a := x[i] - mx
		b := y[i] - my
		num += a * b
		dx2 += a * a
		dy2 += b * b
	}
	den := math.Sqrt(dx2 * dy2)
	if den == 0 {
		return math.NaN()
	}
	return num / den
}

// ---------------------------------------------------------------------------
// profile: spatial profile of one FCS column along an axis
// ---------------------------------------------------------------------------

func runProfile(args []string) error {
	fs := flag.NewFlagSet("profile", flag.ExitOnError)
	fs.Usage = func() {
		fmt.Fprintln(os.Stderr, "usage: scpan profile <fcsdump.txt> [--col em] [--axis y] [--x 8] [--z 8] [--dx 2] [--dz 2] [--bw 1]")
		fs.PrintDefaults()
	}
	c := fs.String("col", "em", "column to profile")
	cx := fs.Float64("x", 8, "region center x")
	cz := fs.Float64("z", 8, "region center z")
	dx := fs.Float64("dx", 2, "region half-width in x")
	dz := fs.Float64("dz", 2, "region half-width in z")
	axis := fs.String("axis", "y", "axis to bin along (x|y|z)")
	bw := fs.Float64("bw", 1, "bin width along the axis")
	frame := fs.String("frame", "last", "frame: last | all | <time>")
	if err := fs.Parse(args); err != nil {
		return err
	}
	in, err := openFCSArg(fs.Arg(0))
	if err != nil {
		return err
	}
	defer in.Close()
	frames, idx, err := readFCS(in)
	if err != nil {
		return err
	}
	sel, err := resolveFrame(frames, *frame)
	if err != nil {
		return err
	}
	for _, n := range []string{"x", "z", *axis, *c} {
		if _, ok := idx[n]; !ok {
			return fmt.Errorf("fcs: column %q missing", n)
		}
	}
	axI := idx[*axis]

	w := bufio.NewWriter(os.Stdout)
	defer w.Flush()
	fmt.Fprintf(w, "t\t%s\t%s_mean\tn\n", *axis, *c)
	for _, fr := range sel {
		type acc struct {
			s float64
			n int
		}
		bins := map[int]*acc{}
		for _, row := range fr.rw {
			x, _ := colVal(idx, row, "x")
			z, _ := colVal(idx, row, "z")
			if math.Abs(x-*cx) > *dx || math.Abs(z-*cz) > *dz {
				continue
			}
			bk := int(math.Round(row[axI] / *bw))
			v, _ := colVal(idx, row, *c)
			bn := bins[bk]
			if bn == nil {
				bn = &acc{}
				bins[bk] = bn
			}
			bn.s += v
			bn.n++
		}
		keys := make([]int, 0, len(bins))
		for k := range bins {
			keys = append(keys, k)
		}
		sort.Ints(keys)
		for _, k := range keys {
			bn := bins[k]
			fmt.Fprintf(w, "%.2f\t%.4f\t%.6f\t%d\n", fr.t, float64(k)**bw, bn.s/float64(bn.n), bn.n)
		}
	}
	return nil
}
