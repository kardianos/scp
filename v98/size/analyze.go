package size

import (
	"fmt"
	"io"
	"sort"

	"scp/v97/hunt"
)

// Options for a size-spectrum pass. Defaults match v97 hunt census.
type Options struct {
	Census hunt.CensusCfg
	// ChartNMax: skip Laplacian if motif.N is above this (default 24).
	ChartNMax int
}

func (o Options) norm() Options {
	if o.Census.EmFloor <= 0 {
		o.Census.EmFloor = 0.05
	}
	if o.Census.MinN <= 0 {
		o.Census.MinN = 3
	}
	if o.ChartNMax <= 0 {
		o.ChartNMax = chartNMax
	}
	return o
}

// MotifSnap is one connected live-Em component on one frame.
type MotifSnap struct {
	T      float64
	N      int
	NEff   float64
	Em     float64
	EmMax  float64
	Chart  LapSpec
	Tagged bool
}

// FrameSnap is the size book of one FCS frame.
type FrameSnap struct {
	T       float64
	NMotif  int
	DomN    int
	DomNEff float64
	DomEm   float64
	EmTot   float64
	Motifs  []MotifSnap
}

// Report is the size spectrum of a whole stream.
type Report struct {
	Path    string
	NFrame  int
	// HistN[n] = motif-counts with that N (all frames).
	HistN map[int]int
	// HistDom[n] = frames whose dominant motif has that N.
	HistDom map[int]int
	Snaps   []FrameSnap
	EmFirst float64
	EmLast  float64
	NFirst  float64 // mean dominant N, first third of frames
	NLast   float64
}

func neff(fr hunt.Frame, cells []int) float64 {
	var s1, s2 float64
	for _, i := range cells {
		e := fr.Cells[i].Em
		if e < 0 {
			continue
		}
		s1 += e
		s2 += e * e
	}
	if s2 <= 0 {
		return 0
	}
	return s1 * s1 / s2
}

// Analyze reads every frame and builds the size spectrum.
func Analyze(path string, opt Options) (*Report, error) {
	opt = opt.norm()
	s, err := hunt.Open(path)
	if err != nil {
		return nil, err
	}
	defer s.Close()
	rep := &Report{
		Path:    path,
		HistN:   map[int]int{},
		HistDom: map[int]int{},
	}
	for {
		fr, err := s.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		mots := hunt.Census(fr, opt.Census)
		fs := FrameSnap{T: fr.T, NMotif: len(mots)}
		var emTot float64
		for _, c := range fr.Cells {
			emTot += c.Em
		}
		fs.EmTot = emTot
		for i, m := range mots {
			ms := MotifSnap{
				T:      fr.T,
				N:      m.N,
				NEff:   neff(fr, m.Cells),
				Em:     m.Em,
				EmMax:  m.EmMax,
				Tagged: m.NTag > 0 && m.NTag*2 >= m.N,
			}
			if m.N <= opt.ChartNMax {
				ms.Chart = analyzeChart(fr, m.Cells)
			}
			fs.Motifs = append(fs.Motifs, ms)
			rep.HistN[m.N]++
			if i == 0 {
				fs.DomN = m.N
				fs.DomNEff = ms.NEff
				fs.DomEm = m.Em
				rep.HistDom[m.N]++
			}
		}
		rep.Snaps = append(rep.Snaps, fs)
	}
	rep.NFrame = len(rep.Snaps)
	if rep.NFrame > 0 {
		rep.EmFirst = rep.Snaps[0].EmTot
		rep.EmLast = rep.Snaps[rep.NFrame-1].EmTot
		a, b := thirdMeans(rep.Snaps)
		rep.NFirst, rep.NLast = a, b
	}
	return rep, nil
}

func thirdMeans(sn []FrameSnap) (first, last float64) {
	if len(sn) == 0 {
		return
	}
	k := len(sn) / 3
	if k < 1 {
		k = 1
	}
	var sf, sl float64
	var nf, nl int
	for i := 0; i < k; i++ {
		if sn[i].DomN > 0 {
			sf += float64(sn[i].DomN)
			nf++
		}
	}
	for i := len(sn) - k; i < len(sn); i++ {
		if sn[i].DomN > 0 {
			sl += float64(sn[i].DomN)
			nl++
		}
	}
	if nf > 0 {
		first = sf / float64(nf)
	}
	if nl > 0 {
		last = sl / float64(nl)
	}
	return
}

func (r *Report) distinctN() int {
	return len(r.HistN)
}

func (r *Report) topBinFrac() float64 {
	var tot, top int
	for _, c := range r.HistN {
		tot += c
		if c > top {
			top = c
		}
	}
	if tot == 0 {
		return 0
	}
	return float64(top) / float64(tot)
}

func (r *Report) topDomFrac() float64 {
	var tot, top int
	for _, c := range r.HistDom {
		tot += c
		if c > top {
			top = c
		}
	}
	if tot == 0 {
		return 0
	}
	return float64(top) / float64(tot)
}

func (r *Report) fracDomN(n int) float64 {
	if r.NFrame == 0 {
		return 0
	}
	return float64(r.HistDom[n]) / float64(r.NFrame)
}

func (r *Report) fracMatch6() float64 {
	var ok, n6 int
	for _, fs := range r.Snaps {
		if len(fs.Motifs) == 0 {
			continue
		}
		m := fs.Motifs[0]
		if m.N != 6 || !m.Chart.OK {
			continue
		}
		n6++
		if m.Chart.Match6 {
			ok++
		}
	}
	if n6 == 0 {
		return 0
	}
	return float64(ok) / float64(n6)
}

func (r *Report) meanNeffOverN() float64 {
	var s float64
	var n int
	for _, fs := range r.Snaps {
		for _, m := range fs.Motifs {
			if m.N > 0 {
				s += m.NEff / float64(m.N)
				n++
			}
		}
	}
	if n == 0 {
		return 0
	}
	return s / float64(n)
}

func (r *Report) histLine() string {
	keys := make([]int, 0, len(r.HistN))
	for k := range r.HistN {
		keys = append(keys, k)
	}
	sort.Ints(keys)
	s := ""
	for i, k := range keys {
		if i > 0 {
			s += " "
		}
		s += fmt.Sprintf("%d:%d", k, r.HistN[k])
	}
	if s == "" {
		s = "(empty)"
	}
	return s
}

func (r *Report) dumpSummary() string {
	return fmt.Sprintf(
		"frames=%d motifs_hist=%s distinct=%d top=%.3f dom_top=%.3f N_eff/N=%.3f Em %.4g→%.4g domN %.2f→%.2f",
		r.NFrame, r.histLine(), r.distinctN(), r.topBinFrac(), r.topDomFrac(),
		r.meanNeffOverN(), r.EmFirst, r.EmLast, r.NFirst, r.NLast)
}

// Dump is the one-line stream summary (stderr / H0).
func (r *Report) Dump() string { return r.dumpSummary() }

// SummaryTSVHeader is the per-frame TSV header for sizespec.
func (r *Report) SummaryTSVHeader() string {
	return "t\tn_motif\tdom_n\tdom_neff\tdom_em\tem_tot\tmatch6"
}

// TSV is one frame, dominant motif.
func (fs FrameSnap) TSV() string {
	m6 := 0
	if len(fs.Motifs) > 0 && fs.Motifs[0].Chart.Match6 {
		m6 = 1
	}
	return fmt.Sprintf("%.6g\t%d\t%d\t%.4f\t%.6g\t%.6g\t%d",
		fs.T, fs.NMotif, fs.DomN, fs.DomNEff, fs.DomEm, fs.EmTot, m6)
}
