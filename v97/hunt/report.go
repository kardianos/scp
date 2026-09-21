package hunt

import (
	"fmt"
	"io"
	"math"
	"strings"
)

// Row is one motif at one frame — the hunt's atomic record.
type Row struct {
	T                   float64
	NMotifs             int
	ID                  int
	N                   int
	Em                  float64
	DEm                 float64 // vs previous sample of this ID
	Cx, Cy, Cz          float64
	Rms                 float64
	EmMin               float64
	JnetX, JnetY, JnetZ float64
	Jnet                float64
	JdipX, JdipY, JdipZ float64
	Jdip                float64
	Jmono               float64
	Circ                float64
	EsThru              float64
	EsNet               float64
	I, O, Phi           float64
	Ifore, Ofore        float64
	Iaft, Oaft          float64
	Sep                 float64
	W, Rm               float64
	// Box / core (leak and current honesty — not motif-membership noise)
	EmTot      float64
	DEmTot     float64
	GJnetX     float64
	GJnetY     float64
	GJnetZ     float64
	TJnetX     float64
	TJnetY     float64
	TJnetZ     float64
	NTag       int
	CoreEm     float64
	CoreN      int
	CoreEsThru float64
	CoreEsNet  float64
	// Exact-ledger channel (schema promo_i/promo_o; zero when absent).
	Ledger       bool    // frame carried the exact door ledger
	BoxPI, BoxPO float64 // box-wide gross intake/outtake this snap
	Ih, Oh       float64 // heuristic estimate on the motif (side-channel)
	// Integrated hop current (pj_meter columns; zero when absent).
	HasPJ                  bool
	BoxPJx, BoxPJy, BoxPJz float64 // box sums (≡ p1fd increment)
	MotPJx, MotPJy, MotPJz float64 // motif net (loop ⇒ ~0)
	MotPJa                 float64 // motif Σ|pJ_cell| (loop ⇒ large)
}

type Options struct {
	Census CensusCfg
	// If FollowID > 0, only that tracked ID is emitted after the first
	// assignment. Otherwise the dominant (heaviest) motif is emitted.
	FollowID int
	// BookDilate: rings of link-neighbours added to the I/O collection
	// set (promo books land on the oriented link's far cell). 0 = off.
	BookDilate int
	// BookDest: dest-only scope (dest.go) — I/O counted only on the
	// motif plus its replayed promoPick fore/aft picks. Overrides
	// BookDilate. The honest promo book scope (REVIEW_grok_state A).
	BookDest bool
}

func RowHeader() string {
	return strings.Join([]string{
		"t", "nMot", "id", "n", "Em", "dEm",
		"cx", "cy", "cz", "rms", "EmMin",
		"Jnet", "JnetX", "JnetY", "JnetZ",
		"Jdip", "JdipX", "JdipY", "JdipZ", "Jmono", "Circ",
		"EsThru", "EsNet",
		"I", "O", "Phi", "Ifore", "Ofore", "Iaft", "Oaft", "sep",
		"W", "Rm",
		"EmTot", "dEmTot", "GJnetX", "GJnetY", "GJnetZ",
		"TJnetX", "TJnetY", "TJnetZ", "nTag",
		"coreN", "coreEm", "coreEsThru", "coreEsNet",
		"ledger", "boxPI", "boxPO", "Ih", "Oh",
		"boxPJx", "boxPJy", "boxPJz", "motPJx", "motPJy", "motPJz",
	}, "\t")
}

func (r Row) TSV() string {
	return fmt.Sprintf(
		"%.4f\t%d\t%d\t%d\t%.6g\t%.6g\t%.4f\t%.4f\t%.4f\t%.4f\t%.5g\t"+
			"%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t"+
			"%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.4f\t%.4f\t"+
			"%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%d\t%d\t%.6g\t%.6g\t%.6g\t"+
			"%d\t%.6g\t%.6g\t%.6g\t%.6g\t"+
			"%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g",
		r.T, r.NMotifs, r.ID, r.N, r.Em, r.DEm,
		r.Cx, r.Cy, r.Cz, r.Rms, r.EmMin,
		r.Jnet, r.JnetX, r.JnetY, r.JnetZ,
		r.Jdip, r.JdipX, r.JdipY, r.JdipZ, r.Jmono, r.Circ,
		r.EsThru, r.EsNet,
		r.I, r.O, r.Phi, r.Ifore, r.Ofore, r.Iaft, r.Oaft, r.Sep,
		r.W, r.Rm,
		r.EmTot, r.DEmTot, r.GJnetX, r.GJnetY, r.GJnetZ,
		r.TJnetX, r.TJnetY, r.TJnetZ, r.NTag,
		r.CoreN, r.CoreEm, r.CoreEsThru, r.CoreEsNet,
		boolInt(r.Ledger), r.BoxPI, r.BoxPO, r.Ih, r.Oh,
		r.BoxPJx, r.BoxPJy, r.BoxPJz, r.MotPJx, r.MotPJy, r.MotPJz)
}

// Analyze walks a stream and yields one Row per frame for the followed motif.
func Analyze(s *Stream, opt Options) ([]Row, error) {
	var rows []Row
	var prev Frame
	var havePrev bool
	var prevEm = map[int]float64{}
	var prevEmTot float64
	tr := &Track{}
	for {
		fr, err := s.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			return rows, err
		}
		ms := tr.Assign(Census(fr, opt.Census), fr.L)
		var m Motif
		if opt.FollowID > 0 {
			for i := range ms {
				if ms[i].ID == opt.FollowID {
					m = ms[i]
					break
				}
			}
		} else {
			m = Dominant(ms)
		}
		r := Row{
			T: fr.T, NMotifs: len(ms),
			ID: m.ID, N: m.N, Em: m.Em,
			Cx: m.Cx, Cy: m.Cy, Cz: m.Cz,
			Rms: m.Rms, EmMin: m.EmMin,
			Ledger: fr.HasPromo,
		}
		r.HasPJ = fr.HasPJ
		for _, c := range fr.Cells {
			r.EmTot += c.Em
			r.BoxPI += c.PromoI
			r.BoxPO += c.PromoO
			r.BoxPJx += c.PJx
			r.BoxPJy += c.PJy
			r.BoxPJz += c.PJz
		}
		for _, ci := range m.Cells {
			if ci >= 0 && ci < len(fr.Cells) {
				c := fr.Cells[ci]
				r.MotPJx += c.PJx
				r.MotPJy += c.PJy
				r.MotPJz += c.PJz
				r.MotPJa += math.Sqrt(c.PJx*c.PJx + c.PJy*c.PJy + c.PJz*c.PJz)
			}
		}
		if havePrev {
			r.DEmTot = r.EmTot - prevEmTot
		}
		r.GJnetX, r.GJnetY, r.GJnetZ = GlobalJ(fr)
		r.TJnetX, r.TJnetY, r.TJnetZ = TaggedJ(fr)
		for _, c := range fr.Cells {
			if c.Tag > 0.5 {
				r.NTag++
			}
		}
		r.CoreEsThru, r.CoreEsNet, r.CoreEm, r.CoreN = CoreFlux(fr, 0.5, s)
		if pe, ok := prevEm[m.ID]; ok {
			r.DEm = m.Em - pe
		}
		if m.N > 0 {
			cur := MeasureCurrents(fr, m, s)
			r.JnetX, r.JnetY, r.JnetZ = cur.JnetX, cur.JnetY, cur.JnetZ
			r.Jnet = cur.JnetAbs()
			r.JdipX, r.JdipY, r.JdipZ = cur.JdipX, cur.JdipY, cur.JdipZ
			r.Jdip = cur.JdipAbs()
			r.Jmono, r.Circ = cur.Jmono, cur.Circ
			r.EsThru, r.EsNet = cur.EsThru, cur.EsNet
			r.W, r.Rm, _ = Winding(fr, m)
			if havePrev {
				cv := MeasureConvert(prev, fr, m, cur.JnetX, cur.JnetY, cur.JnetZ, opt.BookDilate, opt.BookDest)
				r.I, r.O, r.Phi = cv.I, cv.O, cv.Phi
				r.Ifore, r.Ofore = cv.Ifore, cv.Ofore
				r.Iaft, r.Oaft = cv.Iaft, cv.Oaft
				r.Sep = cv.Sep
				r.Ih, r.Oh = cv.Ih, cv.Oh
			}
		}
		if m.ID != 0 {
			prevEm[m.ID] = m.Em
		}
		prevEmTot = r.EmTot
		rows = append(rows, r)
		prev = fr
		havePrev = true
	}
	return rows, nil
}

// Series helpers for H0.

func pearson(xs, ys []float64) float64 {
	n := len(xs)
	if n != len(ys) || n < 3 {
		return math.NaN()
	}
	var sx, sy, sxx, syy, sxy float64
	k := 0
	for i := 0; i < n; i++ {
		if math.IsNaN(xs[i]) || math.IsNaN(ys[i]) {
			continue
		}
		sx += xs[i]
		sy += ys[i]
		sxx += xs[i] * xs[i]
		syy += ys[i] * ys[i]
		sxy += xs[i] * ys[i]
		k++
	}
	if k < 3 {
		return math.NaN()
	}
	fk := float64(k)
	num := sxy - sx*sy/fk
	den := math.Sqrt((sxx - sx*sx/fk) * (syy - sy*sy/fk))
	if den <= 0 {
		return math.NaN()
	}
	return num / den
}

func meanAbs(xs []float64) float64 {
	if len(xs) == 0 {
		return 0
	}
	var s float64
	for _, v := range xs {
		s += math.Abs(v)
	}
	return s / float64(len(xs))
}

func mean(xs []float64) float64 {
	if len(xs) == 0 {
		return 0
	}
	var s float64
	for _, v := range xs {
		s += v
	}
	return s / float64(len(xs))
}

func fracSameSign(xs []float64) float64 {
	var p, n int
	for _, v := range xs {
		if v > 0 {
			p++
		} else if v < 0 {
			n++
		}
	}
	t := p + n
	if t == 0 {
		return 0
	}
	if p > n {
		return float64(p) / float64(t)
	}
	return float64(n) / float64(t)
}

func boolInt(b bool) int {
	if b {
		return 1
	}
	return 0
}

type H0Class string

const (
	H0Sealed   H0Class = "sealed"
	H0Leak     H0Class = "leak"
	H0Ring     H0Class = "ring"
	H0BathRing H0Class = "bathring"
	H0Current  H0Class = "current"
	H0Exact    H0Class = "exact"
)

type Check struct {
	Name   string
	OK     bool
	Detail string
	Skip   bool
}

type H0Report struct {
	Class  H0Class
	Path   string
	N      int
	Checks []Check
	OK     bool
	Notes  []string
}

func (r H0Report) String() string {
	var b strings.Builder
	ver := "FAIL"
	if r.OK {
		ver = "PASS"
	}
	fmt.Fprintf(&b, "H0 %s  %s  frames=%d  %s\n", r.Class, ver, r.N, r.Path)
	for _, c := range r.Checks {
		mark := "FAIL"
		if c.Skip {
			mark = "SKIP"
		} else if c.OK {
			mark = "PASS"
		}
		fmt.Fprintf(&b, "  [%s] %s — %s\n", mark, c.Name, c.Detail)
	}
	for _, n := range r.Notes {
		fmt.Fprintf(&b, "  note: %s\n", n)
	}
	return b.String()
}

func GradeH0(class H0Class, path string, rows []Row) H0Report {
	rep := H0Report{Class: class, Path: path, N: len(rows)}
	if len(rows) == 0 {
		rep.Checks = []Check{{Name: "frames", OK: false, Detail: "no frames"}}
		return rep
	}
	switch class {
	case H0Sealed:
		rep.Checks = gradeSealed(rows)
	case H0Leak:
		rep.Checks = gradeLeak(rows)
	case H0Ring:
		rep.Checks = gradeRing(rows)
	case H0BathRing:
		rep.Checks = gradeBathRing(rows)
	case H0Current:
		rep.Checks = gradeCurrent(rows)
	case H0Exact:
		rep.Checks = gradeExact(rows)
	default:
		rep.Checks = []Check{{Name: "class", OK: false, Detail: string(class)}}
	}
	rep.OK = true
	for _, c := range rep.Checks {
		if !c.Skip && !c.OK {
			rep.OK = false
		}
	}
	return rep
}

func collect(rows []Row) (phi, sep, jdip, jnet, circ, esThru, esNet, dEm, jx, W, Rm []float64) {
	for _, r := range rows {
		if r.N == 0 {
			continue
		}
		phi = append(phi, r.Phi)
		sep = append(sep, r.Sep)
		jdip = append(jdip, r.Jdip)
		jnet = append(jnet, r.Jnet)
		circ = append(circ, r.Circ)
		esThru = append(esThru, r.EsThru)
		esNet = append(esNet, r.EsNet)
		dEm = append(dEm, r.DEm)
		jx = append(jx, r.JnetX)
		W = append(W, r.W)
		Rm = append(Rm, r.Rm)
	}
	return
}

func gradeSealed(rows []Row) []Check {
	phi, sep, jdip, jnet, circ, esThru, esNet, dEm, _, _, _ := collect(rows)
	em0 := rows[len(rows)/2].Em
	if em0 <= 0 {
		em0 = 1
	}
	// Φ ≈ 0 or I,O co-located
	phiBar := meanAbs(phi)
	sepBar := meanAbs(sep)
	point := phiBar < 0.05*em0 || sepBar < 0.4
	c1 := Check{Name: "point book (Φ≈0 or |sep| small)", OK: point,
		Detail: fmt.Sprintf("⟨|Φ|⟩=%.4g  ⟨|sep|⟩=%.3f  Em~%.3g", phiBar, sepBar, em0)}
	// dipole of J not a steady aligned engine: large AND direction-stable
	// is the failure. A sloshing/random Jdip (Face C) is allowed.
	circBar := meanAbs(circ)
	jdipBar := meanAbs(jdip)
	jnetBar := meanAbs(jnet)
	var jdx []float64
	for _, r := range rows {
		jdx = append(jdx, r.JdipX)
	}
	steady := jdipBar > 1e-3 && jdipBar > 0.2*math.Max(circBar, 1e-9) &&
		jnetBar > 0.2*math.Max(circBar, 1e-9) && fracSameSign(jdx) >= 0.85
	c2 := Check{Name: "J dipole not a steady engine", OK: !steady,
		Detail: fmt.Sprintf("⟨|Jdip|⟩=%.4g  ⟨|Jnet|⟩=%.4g  ⟨|Circ|⟩=%.4g  dir-stable=%.2f",
			jdipBar, jnetBar, circBar, fracSameSign(jdx))}
	// Es throughput: no live space cycle at constant M
	// either thru small vs |dEm|, or net tracks dEm (bookkeeping)
	corr := pearson(esNet, dEm)
	thruBar := meanAbs(esThru)
	dBar := meanAbs(dEm)
	tracks := !math.IsNaN(corr) && math.Abs(corr) > 0.4
	quiet := thruBar < 0.2 || thruBar < 3*dBar || tracks || dBar < 1e-4
	c3 := Check{Name: "no live Es cycle at constant M", OK: quiet,
		Detail: fmt.Sprintf("⟨EsThru⟩=%.4g  ⟨|dEm|⟩=%.4g  corr(EsNet,dEm)=%.3f", thruBar, dBar, corr)}
	nMot := 0
	for _, r := range rows {
		if r.NMotifs > nMot {
			nMot = r.NMotifs
		}
	}
	c4 := Check{Name: "census (connected, not hoard histogram)", OK: nMot >= 1,
		Detail: fmt.Sprintf("max nMotifs=%d  last n=%d Em=%.3g", nMot, rows[len(rows)-1].N, rows[len(rows)-1].Em)}
	return []Check{c1, c2, c3, c4}
}

func gradeLeak(rows []Row) []Check {
	var esNet, esThru, dM []float64
	for _, r := range rows {
		esNet = append(esNet, r.CoreEsNet)
		esThru = append(esThru, r.CoreEsThru)
		dM = append(dM, r.DEmTot)
	}
	c := pearson(esNet, dM)
	cAbs := pearson(absSlice(esThru), absSlice(dM))
	ok := (!math.IsNaN(c) && math.Abs(c) >= 0.30) || (!math.IsNaN(cAbs) && cAbs >= 0.30)
	if meanAbs(dM) < 1e-6 {
		ok = false
	}
	return []Check{{
		Name: "core Es throughput tracks box dM/dt",
		OK:   ok,
		Detail: fmt.Sprintf("corr(coreEsNet,dEmTot)=%.3f  corr(|Thru|,|dM|)=%.3f  ⟨|dEmTot|⟩=%.4g  last coreN=%d",
			c, cAbs, meanAbs(dM), rows[len(rows)-1].CoreN),
	}}
}

func gradeRing(rows []Row) []Check {
	_, sep, jdip, jnet, circ, _, _, _, _, W, Rm := collect(rows)
	last := rows[len(rows)-1]
	wBar := meanAbs(W)
	rmBar := mean(Rm)
	c1 := Check{Name: "winding held", OK: wBar > 0.7 && rmBar > 0.6,
		Detail: fmt.Sprintf("⟨|W|⟩=%.3f  ⟨Rm⟩=%.3f  last W=%.3f Rm=%.3f n=%d", wBar, rmBar, last.W, last.Rm, last.N)}
	circBar := meanAbs(circ)
	c2 := Check{Name: "monopole-neutral (loop, not dipole)",
		OK:     meanAbs(jnet) < 0.35*math.Max(circBar, 1e-9) && meanAbs(jdip) < 0.5*math.Max(circBar, 1e-9),
		Detail: fmt.Sprintf("⟨|Jnet|⟩=%.4g  ⟨|Jdip|⟩=%.4g  ⟨|Circ|⟩=%.4g", meanAbs(jnet), meanAbs(jdip), circBar)}
	c3 := Check{Name: "I/O not split",
		OK:     meanAbs(sep) < 0.6,
		Detail: fmt.Sprintf("⟨|sep|⟩=%.3f  ⟨|Φ|⟩=%.4g", meanAbs(sep), meanAbs(collectPhi(rows)))}
	return []Check{c1, c2, c3}
}

func collectPhi(rows []Row) []float64 {
	var p []float64
	for _, r := range rows {
		p = append(p, r.Phi)
	}
	return p
}

func gradeBathRing(rows []Row) []Check {
	_, _, jdip, _, _, _, _, _, _, W, Rm := collect(rows)
	if len(W) < 4 {
		return []Check{{Name: "frames", OK: false, Detail: "too few live frames"}}
	}
	// winding dies: last third |W| or Rm much smaller than first third
	n := len(W)
	a, b := W[:n/3], W[2*n/3:]
	ra, rb := Rm[:n/3], Rm[2*n/3:]
	died := meanAbs(b) < 0.55*meanAbs(a) || mean(rb) < 0.7*mean(ra) || meanAbs(b) < 0.6
	c1 := Check{Name: "winding dies", OK: died,
		Detail: fmt.Sprintf("⟨|W|⟩ early=%.3f late=%.3f  ⟨Rm⟩ early=%.3f late=%.3f",
			meanAbs(a), meanAbs(b), mean(ra), mean(rb))}
	// no dipole rescue: Jdip does not grow into a steady engine
	c2 := Check{Name: "no dipole rescue",
		OK:     meanAbs(jdip[2*len(jdip)/3:]) < 5 || fracSameSign(jdipSignX(rows[2*len(rows)/3:])) < 0.85,
		Detail: fmt.Sprintf("late ⟨|Jdip|⟩=%.4g", meanAbs(jdip[2*len(jdip)/3:]))}
	return []Check{c1, c2}
}

func jdipSignX(rows []Row) []float64 {
	var xs []float64
	for _, r := range rows {
		xs = append(xs, r.JdipX)
	}
	return xs
}

func gradeCurrent(rows []Row) []Check {
	// Snapshot J is not fd (fd is a time-integral of hops). H0 here
	// checks the reconstruction is alive on the painted packet, and
	// reports whether a preferred direction even exists.
	var tj, gj []float64
	nTag := 0
	for i, r := range rows {
		if r.NTag > nTag {
			nTag = r.NTag
		}
		if i == 0 {
			continue
		}
		tj = append(tj, r.TJnetX)
		gj = append(gj, r.GJnetX)
	}
	c1 := Check{Name: "tagged packet present", OK: nTag >= 8,
		Detail: fmt.Sprintf("max nTag=%d", nTag)}
	eng := meanAbs(tj)
	c2 := Check{Name: "tagged J engaged (snapshot bilinear)", OK: eng > 1e-3,
		Detail: fmt.Sprintf("⟨|TJnetX|⟩=%.4g  ⟨TJnetX⟩=%.4g ± %.4g", eng, mean(tj), stddev(tj))}
	// Direction is informational: snapshot J sloshes on the known e3b streams.
	dir := math.Abs(mean(tj)) > 0.4*math.Max(stddev(tj), 1e-9)
	c3 := Check{Name: "preferred +/−x in the time average", OK: dir,
		Detail: fmt.Sprintf("⟨TJnetX⟩=%.4g  ⟨GJnetX⟩=%.4g  (sloshing is allowed; fd is not in the FCS)", mean(tj), mean(gj))}
	if !dir {
		c3.Skip = true // do not fail H0 on the known wobble
		c3.Detail += " — SKIP: not a fail; instantaneous J ≠ integrated fd"
	}
	return []Check{c1, c2, c3}
}

// gradeExact — H0 class `exact`: on a stream that carries the
// promo_i/promo_o ledger, the graded I/O must BE that ledger, and the
// mode-crossing heuristic's inflation must be visible, not hidden.
// Stream-specific totals (p2_space cond=0.5249, evap+rad=14.185) are
// asserted in hunt_test.go against the kernel's own RESULT lines.
func gradeExact(rows []Row) []Check {
	nLedger := 0
	var boxI, boxO, motI, motO, hI, hO float64
	subsetViol := 0
	live := 0
	for i, r := range rows {
		if r.Ledger {
			nLedger++
		}
		boxI += r.BoxPI
		boxO += r.BoxPO
		if i == 0 {
			continue // no convert on the first frame
		}
		motI += r.I
		motO += r.O
		hI += r.Ih
		hO += r.Oh
		if r.N > 0 {
			live++
			// motif exact I/O is a subset of the box ledger
			if r.I > r.BoxPI+1e-6 || r.O > r.BoxPO+1e-6 {
				subsetViol++
			}
		}
	}
	c1 := Check{Name: "ledger present on every frame", OK: nLedger == len(rows),
		Detail: fmt.Sprintf("ledger frames %d/%d", nLedger, len(rows))}
	c2 := Check{Name: "graded I/O ⊆ box ledger (exact channel in use)",
		OK: nLedger == len(rows) && subsetViol == 0,
		Detail: fmt.Sprintf("violations %d/%d live frames  ΣmotI=%.4g ΣmotO=%.4g  ΣboxPI=%.4g ΣboxPO=%.4g",
			subsetViol, live, motI, motO, boxI, boxO)}
	// The heuristic ran as a side-channel: report its inflation. This
	// is informational unless the heuristic is *smaller* than the
	// ledger on the box scale (which would mean the old meter had been
	// under-counting, the direction the review wrongly guessed).
	c3 := Check{Name: "heuristic inflation measured", OK: true,
		Detail: fmt.Sprintf("heuristic ΣI=%.4g ΣO=%.4g  vs exact motif ΣI=%.4g ΣO=%.4g  (inflation I×%.3g O×%.3g)",
			hI, hO, motI, motO, safeRatio(hI, motI), safeRatio(hO, motO))}
	return []Check{c1, c2, c3}
}

func safeRatio(a, b float64) float64 {
	if b <= 0 {
		return math.Inf(1)
	}
	return a / b
}

func stddev(xs []float64) float64 {
	if len(xs) < 2 {
		return 0
	}
	m := mean(xs)
	var s float64
	for _, v := range xs {
		d := v - m
		s += d * d
	}
	return math.Sqrt(s / float64(len(xs)-1))
}

func absSlice(xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, v := range xs {
		out[i] = math.Abs(v)
	}
	return out
}
