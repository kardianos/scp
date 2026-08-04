package fab

// FCS v3 chunked stream writer (see v90/FCS.md) — same bytes as the C
// kernel: CFG provenance, SCHM column schema, CELL+LINK frames, ANLZ
// instrumentation frames interleaved. snap_dir writes per-frame v3
// files; snap_file appends one stream. Pure output — consumes no RNG.

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"math"
	"os"
	"path/filepath"

	"github.com/klauspost/compress/zstd"
)

// codec 1: columnar transpose (4-byte units per record) + byte shuffle +
// zstd (per-chunk self-contained; measured ~0.46-0.55 at ~175 MB/s).
var fcsZenc, _ = zstd.NewWriter(nil, zstd.WithEncoderLevel(zstd.SpeedFastest))

func fcsColT(p []byte, hdr, stride int) []byte {
	out := make([]byte, len(p))
	copy(out, p[:hdr])
	body := p[hdr:]
	n := len(body) / (4 * stride)
	o := out[hdr:]
	for c := 0; c < stride; c++ {
		for i := 0; i < n; i++ {
			copy(o[(c*n+i)*4:(c*n+i)*4+4], body[(i*stride+c)*4:])
		}
	}
	return out
}

func fcsShuffle4(p []byte) []byte {
	m := len(p) / 4
	out := make([]byte, len(p))
	for i := 0; i < m; i++ {
		out[i] = p[4*i]
		out[m+i] = p[4*i+1]
		out[2*m+i] = p[4*i+2]
		out[3*m+i] = p[4*i+3]
	}
	copy(out[4*m:], p[4*m:])
	return out
}

var fcsCellCols = []string{"x", "y", "z", "r", "es", "em", "ee", "xload", "tag", "fa1", "fa2", "th2"}
var fcsLinkCols = []string{"d", "A", "lem", "gg"}

type fcsW struct {
	f    *os.File
	w    *bufio.Writer
	comp bool
}

// write one chunk; CELL/LINK optionally wrapped in CMPD (codec 1)
func (o *fcsW) chunkOut(cc string, p []byte, hdr, stride int) {
	if o.comp && stride > 0 {
		t := fcsShuffle4(fcsColT(p, hdr, stride))
		c := fcsZenc.EncodeAll(t, nil)
		if len(c)+13 < len(p) {
			o.w.WriteString("CMPD")
			o.u64(uint64(4 + 8 + 1 + len(c)))
			o.w.WriteString(cc)
			o.u64(uint64(len(p)))
			o.w.WriteByte(1)
			o.w.Write(c)
			o.w.Flush()
			return
		}
	}
	o.w.WriteString(cc)
	o.u64(uint64(len(p)))
	o.w.Write(p)
	o.w.Flush()
}

func (o *fcsW) u32(v uint32) {
	var b [4]byte
	binary.LittleEndian.PutUint32(b[:], v)
	o.w.Write(b[:])
}

func (o *fcsW) u64(v uint64) {
	var b [8]byte
	binary.LittleEndian.PutUint64(b[:], v)
	o.w.Write(b[:])
}

func (o *fcsW) f64(v float64) { o.u64(math.Float64bits(v)) }
func (o *fcsW) f32(v float32) { o.u32(math.Float32bits(v)) }

func (s *Sim) fcsBegin(o *fcsW) {
	P := &s.P
	o.w.WriteString("FCS1")
	o.u32(3)
	cfg := fmt.Sprintf(
		"exp=%s seed=%d L=%.6g dt=%.6g T=%.6g\n"+
			"laws_V2g: C=%.6g w1=%.6g w2=%.6g q_detune=%.6g gamma_res=%.6g gamma_res_m=%.6g "+
			"p_gate=%.6g lock_floor=%.6g k_dep=%.6g k_dep_m=%.6g cap=%.6g e_s0=%.6g es_floor=%.6g "+
			"e_cond=%.6g f_conv=%.6g f_evap=%.6g s_pull=%.6g kappa_lock=%.6g kappa_align=%.6g "+
			"kappa_reac=%.6g s_k=%.6g s_disp=%.6g sigma_tumble=%.6g comb_limit=%d rough_k=%.6g "+
			"gamma_rough=%.6g mob_sym=%d mob_floor=%.6g field_J=%.6g quant_A0=%.6g quant_mode=%d\n"+
			"geometry: cfac=%.6g k_rep=%.6g mob_geo=%.6g kappa_bond=%.6g freeze_geo=%d "+
			"bath=%d bath_frac=%.6g\n",
		P.Exp, P.Seed, P.L, P.Dt, P.T,
		P.C, P.W1, P.W2, P.QDetune, P.GammaRes, P.GammaResM,
		P.PGate, P.LockFloor, P.KDep, P.KDepM, P.Cap, P.Es0, P.EsFloor,
		P.ECond, P.FConv, P.FEvap, P.SPull, P.KappaLock, P.KappaAlign,
		P.KappaReac, P.SK, P.SDisp, P.SigmaTumble, P.CombLimit, P.RoughK,
		P.GammaRough, P.MobSym, P.MobFloor, P.FieldJ, P.QuantA0, P.QuantMode,
		P.Cfac, P.KRep, P.MobGeo, P.KappaBond, P.FreezeGeo,
		P.Bath, P.BathFrac)
	o.w.WriteString("CFG ")
	o.u64(uint64(len(cfg)))
	o.w.WriteString(cfg)
	// SCHM
	var sb []byte
	app32 := func(v uint32) {
		var b [4]byte
		binary.LittleEndian.PutUint32(b[:], v)
		sb = append(sb, b[:]...)
	}
	app32(uint32(len(fcsCellCols)))
	for _, c := range fcsCellCols {
		sb = append(sb, byte(len(c)))
		sb = append(sb, c...)
	}
	app32(uint32(len(fcsLinkCols)))
	for _, c := range fcsLinkCols {
		sb = append(sb, byte(len(c)))
		sb = append(sb, c...)
	}
	o.w.WriteString("SCHM")
	o.u64(uint64(len(sb)))
	o.w.Write(sb)
}

func (s *Sim) fcsCellFrame(o *fcsW) {
	nc := len(fcsCellCols)
	cp := make([]byte, 0, 20+s.NC*nc*4)
	put64 := func(v float64) {
		var b [8]byte
		binary.LittleEndian.PutUint64(b[:], math.Float64bits(v))
		cp = append(cp, b[:]...)
	}
	put32 := func(v uint32) {
		var b [4]byte
		binary.LittleEndian.PutUint32(b[:], v)
		cp = append(cp, b[:]...)
	}
	putf := func(v float32) { put32(math.Float32bits(v)) }
	put64(s.simT)
	put64(s.P.L)
	put32(uint32(s.NC))
	for i := 0; i < s.NC; i++ {
		xload := (s.Em[i] + s.flload[i]) / s.P.Cap
		putf(float32(s.px[i]))
		putf(float32(s.py[i]))
		putf(float32(s.pz[i]))
		putf(float32(s.cr[i]))
		putf(float32(s.Es[i]))
		putf(float32(s.Em[i]))
		putf(float32(s.Ee[i]))
		putf(float32(xload))
		putf(float32(s.tag[i]))
		putf(float32(s.fa1[i]))
		putf(float32(s.fa2[i]))
		putf(float32(s.th2[i]))
	}
	o.chunkOut("CELL", cp, 20, nc)

	var nl uint32
	for sl := 0; sl < s.NSLOT; sl++ {
		if s.sst[sl] != sFree {
			nl++
		}
	}
	lp := make([]byte, 0, 12+int(nl)*(8+len(fcsLinkCols)*4))
	cp = lp // reuse the append helpers on lp via cp alias
	put64(s.simT)
	put32(nl)
	for sl := 0; sl < s.NSLOT; sl++ {
		if s.sst[sl] == sFree {
			continue
		}
		i, j := s.sli[sl], s.slj[sl]
		gf := s.gateOf(float64(s.slq[sl])*s.th2[i] - float64(s.slq[sl])*s.w2e[i]*s.sd[sl]/s.P.C - float64(s.slp[sl])*s.th2[j])
		gb := s.gateOf(float64(s.slp[sl])*s.th2[j] - float64(s.slp[sl])*s.w2e[j]*s.sd[sl]/s.P.C - float64(s.slq[sl])*s.th2[i])
		put32(uint32(i))
		put32(uint32(j))
		putf(float32(s.sd[sl]))
		putf(float32(s.sA[sl]))
		putf(float32(s.slem[2*sl] + s.slem[2*sl+1]))
		putf(float32(gf * gb))
	}
	o.chunkOut("LINK", cp, 12, 2+len(fcsLinkCols))
}

func (s *Sim) fcsAnlzTable(o *fcsW, name string, cols []string, rows []float64, nrows int) {
	plen := uint64(1 + 1 + len(name) + 8 + 4 + 4 + nrows*len(cols)*8)
	for _, c := range cols {
		plen += uint64(1 + len(c))
	}
	o.w.WriteString("ANLZ")
	o.u64(plen)
	o.w.WriteByte(1)
	o.w.WriteByte(byte(len(name)))
	o.w.WriteString(name)
	o.f64(s.simT)
	o.u32(uint32(len(cols)))
	for _, c := range cols {
		o.w.WriteByte(byte(len(c)))
		o.w.WriteString(c)
	}
	o.u32(uint32(nrows))
	for _, v := range rows[:nrows*len(cols)] {
		o.f64(v)
	}
	o.w.Flush()
}

func (s *Sim) fcsAnlzText(o *fcsW, txt string) {
	o.w.WriteString("ANLZ")
	o.u64(uint64(1 + len(txt)))
	o.w.WriteByte(0)
	o.w.WriteString(txt)
	o.w.Flush()
}

var fcsMeterCols = []string{"drift", "phi", "z_live", "nla", "nld", "births", "deaths"}

func (s *Sim) fcsInstrument(o *fcsW) {
	Etot := s.totalEnergy()
	den := s.E0
	if den == 0 {
		den = 1
	}
	phi, zl, nla, nld, _, _, _ := s.geoStats()
	rows := []float64{(Etot - s.E0) / den, phi, zl, float64(nla), float64(nld),
		float64(s.births), float64(s.deaths)}
	s.fcsAnlzTable(o, "meters", fcsMeterCols, rows, 1)
	if s.P.Exp == "slit" {
		var dsRows [dsNBin * 2]float64
		for b := 0; b < dsNBin; b++ {
			dsRows[2*b] = (float64(b) + 0.5) * s.P.L / dsNBin
			dsRows[2*b+1] = s.dsI[b]
		}
		s.fcsAnlzTable(o, "ds_screen", []string{"y", "I"}, dsRows[:], dsNBin)
	}
	if s.P.SectMeter != 0 {
		var xRows [sectMax * 3]float64
		for k := 0; k < s.P.SectN; k++ {
			thc := float64(k) * TwoPi / float64(s.P.SectN)
			if thc > math.Pi {
				thc -= TwoPi
			}
			xRows[3*k] = thc
			xRows[3*k+1] = s.sectE[k]
			xRows[3*k+2] = s.sectN2[k]
		}
		s.fcsAnlzTable(o, "sector", []string{"th", "E", "n"}, xRows[:], s.P.SectN)
	}
}

func fcsOpen(path string, comp bool) (*fcsW, error) {
	f, err := os.Create(path)
	if err != nil {
		return nil, err
	}
	return &fcsW{f: f, w: bufio.NewWriterSize(f, 1<<20), comp: comp}, nil
}

func (o *fcsW) close() {
	o.w.Flush()
	o.f.Close()
}

// per-frame v3 file under snap_dir (a v3 stream with one frame)
func (s *Sim) writeFCS(idx int) {
	if err := os.MkdirAll(s.P.SnapDir, 0o755); err != nil {
		fprintf(s.Errw, "# WARN snap mkdir: %v\n", err)
		return
	}
	path := filepath.Join(s.P.SnapDir, fmt.Sprintf("snap_%05d_t%.3f.fcs", idx, s.simT))
	o, err := fcsOpen(path, s.P.SnapComp != 0)
	if err != nil {
		fprintf(s.Errw, "# WARN snap create: %v\n", err)
		return
	}
	s.fcsBegin(o)
	s.fcsCellFrame(o)
	o.close()
}
