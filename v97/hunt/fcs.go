package hunt

// FCS v3 stream reader. Maps CELL/LINK columns by SCHM name when
// present; falls back to the v90/v97 writer schema. Skips ANLZ and
// unknown chunks. CMPD codec 1 is the same shuffle+colT+zstd used
// by fab/fcs.go and cmd/fcsdump.

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"os"
	"strconv"
	"strings"

	"github.com/klauspost/compress/zstd"
)

const (
	defaultL     = 16.0
	defaultDt    = 0.02
	defaultSK    = 0.06
	defaultSDisp = 0.3
	defaultR0    = 0.85
	defaultCap   = 2.5
)

var (
	defCellCols = []string{"x", "y", "z", "r", "es", "em", "ee", "xload", "tag", "fa1", "fa2", "th2"}
	defLinkCols = []string{"d", "A", "lem", "gg"}
)

type Stream struct {
	Cfg      string
	Keys     map[string]float64
	CellCols []string
	LinkCols []string
	L        float64
	Dt       float64
	SK       float64
	SDisp    float64
	R0       float64
	Cap      float64
	Aref     float64
	Dref     float64

	r       *bufio.Reader
	zdec    *zstd.Decoder
	f       *os.File
	pending *Frame // CELL already parsed, waiting for its LINK or a return
}

type Cell struct {
	X, Y, Z, R     float64
	Es, Em, Ee     float64
	Xload, Tag     float64
	Fa1, Fa2, Th2  float64
	PromoI, PromoO float64
}

type Link struct {
	I, J      int
	D, A, Lem float64
	GG        float64
}

type Frame struct {
	T     float64
	L     float64
	Cells []Cell
	Links []Link
}

func Open(path string) (*Stream, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	s := &Stream{
		Keys:     map[string]float64{},
		CellCols: append([]string(nil), defCellCols...),
		LinkCols: append([]string(nil), defLinkCols...),
		L:        defaultL,
		Dt:       defaultDt,
		SK:       defaultSK,
		SDisp:    defaultSDisp,
		R0:       defaultR0,
		Cap:      defaultCap,
		f:        f,
		r:        bufio.NewReaderSize(f, 1<<20),
	}
	s.zdec, _ = zstd.NewReader(nil)
	magic := make([]byte, 8)
	if _, err := io.ReadFull(s.r, magic); err != nil || string(magic[:4]) != "FCS1" {
		f.Close()
		return nil, fmt.Errorf("not an FCS v3 stream: %s", path)
	}
	s.Aref = math.Pi * s.R0 * s.R0
	s.Dref = 2 * s.R0
	return s, nil
}

func (s *Stream) Close() error {
	if s.f != nil {
		return s.f.Close()
	}
	return nil
}

func (s *Stream) applyKeys() {
	if v, ok := s.Keys["L"]; ok && v > 0 {
		s.L = v
	}
	if v, ok := s.Keys["dt"]; ok && v > 0 {
		s.Dt = v
	}
	if v, ok := s.Keys["s_k"]; ok {
		s.SK = v
	}
	if v, ok := s.Keys["s_disp"]; ok {
		s.SDisp = v
	}
	if v, ok := s.Keys["r0"]; ok && v > 0 {
		s.R0 = v
	}
	if v, ok := s.Keys["cap"]; ok && v > 0 {
		s.Cap = v
	}
	s.Aref = math.Pi * s.R0 * s.R0
	s.Dref = 2 * s.R0
}

func parseKeys(cfg string) map[string]float64 {
	out := map[string]float64{}
	for _, tok := range strings.FieldsFunc(cfg, func(r rune) bool {
		return r == ' ' || r == '\n' || r == '\t' || r == '\r'
	}) {
		k, v, ok := strings.Cut(tok, "=")
		if !ok || k == "" {
			continue
		}
		// strip trailing punctuation from prose dumps
		v = strings.TrimRight(v, ",;")
		f, err := strconv.ParseFloat(v, 64)
		if err != nil {
			continue
		}
		out[k] = f
	}
	return out
}

func (s *Stream) readChunk() (string, []byte, error) {
	tag := make([]byte, 4)
	if _, err := io.ReadFull(s.r, tag); err != nil {
		return "", nil, err
	}
	lb := make([]byte, 8)
	if _, err := io.ReadFull(s.r, lb); err != nil {
		return "", nil, err
	}
	n := binary.LittleEndian.Uint64(lb)
	if n > 1<<32 {
		return "", nil, fmt.Errorf("chunk %q implausibly large (%d)", tag, n)
	}
	p := make([]byte, n)
	if _, err := io.ReadFull(s.r, p); err != nil {
		return "", nil, err
	}
	ct := string(tag)
	if ct != "CMPD" {
		return ct, p, nil
	}
	inner := string(p[:4])
	rawlen := binary.LittleEndian.Uint64(p[4:12])
	codec := p[12]
	if codec != 1 {
		return inner, nil, fmt.Errorf("unknown CMPD codec %d", codec)
	}
	raw, err := s.zdec.DecodeAll(p[13:], nil)
	if err != nil || uint64(len(raw)) != rawlen {
		return inner, nil, fmt.Errorf("CMPD decode: %v", err)
	}
	hdr, stride := 20, len(s.CellCols)
	if inner == "LINK" {
		hdr, stride = 12, 2+len(s.LinkCols)
	}
	return inner, unColT(unshuffle4(raw), hdr, stride), nil
}

func unshuffle4(p []byte) []byte {
	m := len(p) / 4
	out := make([]byte, len(p))
	for i := 0; i < m; i++ {
		out[4*i] = p[i]
		out[4*i+1] = p[m+i]
		out[4*i+2] = p[2*m+i]
		out[4*i+3] = p[3*m+i]
	}
	copy(out[4*m:], p[4*m:])
	return out
}

func unColT(p []byte, hdr, stride int) []byte {
	out := make([]byte, len(p))
	copy(out, p[:hdr])
	body := p[hdr:]
	if stride <= 0 || len(body)%(4*stride) != 0 {
		return out
	}
	n := len(body) / (4 * stride)
	o := out[hdr:]
	for c := 0; c < stride; c++ {
		for i := 0; i < n; i++ {
			copy(o[(i*stride+c)*4:(i*stride+c)*4+4], body[(c*n+i)*4:])
		}
	}
	return out
}

func f32(b []byte) float64 { return float64(math.Float32frombits(binary.LittleEndian.Uint32(b))) }
func f64(b []byte) float64 { return math.Float64frombits(binary.LittleEndian.Uint64(b)) }
func u32(b []byte) uint32  { return binary.LittleEndian.Uint32(b) }

func (s *Stream) parseSCHM(p []byte) {
	if len(p) < 4 {
		return
	}
	off := 0
	readNames := func() []string {
		if off+4 > len(p) {
			return nil
		}
		n := int(u32(p[off:]))
		off += 4
		names := make([]string, 0, n)
		for i := 0; i < n; i++ {
			if off >= len(p) {
				break
			}
			ln := int(p[off])
			off++
			if off+ln > len(p) {
				break
			}
			names = append(names, string(p[off:off+ln]))
			off += ln
		}
		return names
	}
	if cc := readNames(); len(cc) > 0 {
		s.CellCols = cc
	}
	if lc := readNames(); len(lc) > 0 {
		s.LinkCols = lc
	}
}

func colIdx(cols []string, name string) int {
	for i, c := range cols {
		if c == name {
			return i
		}
	}
	return -1
}

func (s *Stream) parseCELL(p []byte) (Frame, error) {
	if len(p) < 20 {
		return Frame{}, fmt.Errorf("short CELL")
	}
	fr := Frame{T: f64(p[0:8]), L: f64(p[8:16])}
	if fr.L > 0 {
		s.L = fr.L
	} else {
		fr.L = s.L
	}
	N := int(u32(p[16:20]))
	nc := len(s.CellCols)
	rec := p[20:]
	if nc <= 0 || len(rec) < N*nc*4 {
		return Frame{}, fmt.Errorf("CELL body short: N=%d cols=%d bytes=%d", N, nc, len(rec))
	}
	ix := colIdx(s.CellCols, "x")
	iy := colIdx(s.CellCols, "y")
	iz := colIdx(s.CellCols, "z")
	ir := colIdx(s.CellCols, "r")
	ies := colIdx(s.CellCols, "es")
	iem := colIdx(s.CellCols, "em")
	iee := colIdx(s.CellCols, "ee")
	ixl := colIdx(s.CellCols, "xload")
	itg := colIdx(s.CellCols, "tag")
	ifa1 := colIdx(s.CellCols, "fa1")
	ifa2 := colIdx(s.CellCols, "fa2")
	ith := colIdx(s.CellCols, "th2")
	ipi := colIdx(s.CellCols, "promo_i")
	ipo := colIdx(s.CellCols, "promo_o")
	get := func(b []byte, i int) float64 {
		if i < 0 {
			return 0
		}
		return f32(b[i*4:])
	}
	fr.Cells = make([]Cell, N)
	for i := 0; i < N; i++ {
		b := rec[i*nc*4:]
		fr.Cells[i] = Cell{
			X: get(b, ix), Y: get(b, iy), Z: get(b, iz), R: get(b, ir),
			Es: get(b, ies), Em: get(b, iem), Ee: get(b, iee),
			Xload: get(b, ixl), Tag: get(b, itg),
			Fa1: get(b, ifa1), Fa2: get(b, ifa2), Th2: get(b, ith),
			PromoI: get(b, ipi), PromoO: get(b, ipo),
		}
	}
	return fr, nil
}

func (s *Stream) parseLINK(p []byte) ([]Link, float64, error) {
	if len(p) < 12 {
		return nil, 0, fmt.Errorf("short LINK")
	}
	t := f64(p[0:8])
	nl := int(u32(p[8:12]))
	nf := len(s.LinkCols)
	rec := p[12:]
	stride := 8 + nf*4
	if nf < 0 || len(rec) < nl*stride {
		return nil, t, fmt.Errorf("LINK body short: nl=%d cols=%d bytes=%d", nl, nf, len(rec))
	}
	id := colIdx(s.LinkCols, "d")
	iA := colIdx(s.LinkCols, "A")
	il := colIdx(s.LinkCols, "lem")
	ig := colIdx(s.LinkCols, "gg")
	get := func(b []byte, i int) float64 {
		if i < 0 {
			return 0
		}
		return f32(b[8+i*4:])
	}
	out := make([]Link, nl)
	for k := 0; k < nl; k++ {
		b := rec[k*stride:]
		out[k] = Link{
			I: int(u32(b[0:4])), J: int(u32(b[4:8])),
			D: get(b, id), A: get(b, iA), Lem: get(b, il), GG: get(b, ig),
		}
	}
	return out, t, nil
}

// Next returns the next CELL frame. A LINK at the same t is attached
// when present. A following CELL is held in pending for the next call.
func (s *Stream) Next() (Frame, error) {
	var fr Frame
	have := false
	if s.pending != nil {
		fr = *s.pending
		s.pending = nil
		have = true
	}
	for {
		ct, p, err := s.readChunk()
		if err == io.EOF {
			if have {
				return fr, nil
			}
			return Frame{}, io.EOF
		}
		if err != nil {
			if have {
				return fr, nil
			}
			return Frame{}, err
		}
		switch ct {
		case "CFG ":
			s.Cfg = string(p)
			s.Keys = parseKeys(s.Cfg)
			s.applyKeys()
		case "SCHM":
			s.parseSCHM(p)
		case "CELL":
			nf, err := s.parseCELL(p)
			if err != nil {
				return Frame{}, err
			}
			if have {
				s.pending = &nf
				return fr, nil
			}
			fr = nf
			have = true
		case "LINK":
			ls, t, err := s.parseLINK(p)
			if err != nil {
				return Frame{}, err
			}
			if have && almostEq(t, fr.T) {
				fr.Links = ls
				return fr, nil
			}
		}
	}
}

func almostEq(a, b float64) bool {
	d := a - b
	if d < 0 {
		d = -d
	}
	return d < 1e-9
}
