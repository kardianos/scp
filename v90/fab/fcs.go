package fab

// .fcs snapshot writer — FCS version 2 (see v90/FCS.md): per-cell state
// including the two-plane field amplitude + dense phase, followed by the
// live channel graph. Byte-compatible with the C kernel's write_fcs().
// Pure output — consumes no RNG, cannot perturb determinism.

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"math"
	"os"
	"path/filepath"
)

func (s *Sim) writeFCS(idx int) {
	if err := os.MkdirAll(s.P.SnapDir, 0o755); err != nil {
		fprintf(s.Errw, "# WARN snap mkdir: %v\n", err)
		return
	}
	path := filepath.Join(s.P.SnapDir, fmt.Sprintf("snap_%05d_t%.3f.fcs", idx, s.simT))
	f, err := os.Create(path)
	if err != nil {
		fprintf(s.Errw, "# WARN snap create: %v\n", err)
		return
	}
	w := bufio.NewWriterSize(f, 1<<20)
	pu32 := func(v uint32) {
		var b [4]byte
		binary.LittleEndian.PutUint32(b[:], v)
		w.Write(b[:])
	}
	pf64 := func(v float64) {
		var b [8]byte
		binary.LittleEndian.PutUint64(b[:], math.Float64bits(v))
		w.Write(b[:])
	}
	w.WriteString("FCS1")
	pu32(2) // version
	pu32(uint32(s.NC))
	pf64(s.simT)
	pf64(s.P.L)
	var rec [12 * 4]byte
	for i := 0; i < s.NC; i++ {
		xload := (s.Em[i] + s.flload[i]) / s.P.Cap
		vals := [12]float32{
			float32(s.px[i]), float32(s.py[i]), float32(s.pz[i]),
			float32(s.cr[i]), float32(s.Es[i]), float32(s.Em[i]),
			float32(s.Ee[i]), float32(xload), float32(s.tag[i]),
			float32(s.fa1[i]), float32(s.fa2[i]), float32(s.th2[i]),
		}
		for k, v := range vals {
			binary.LittleEndian.PutUint32(rec[4*k:], math.Float32bits(v))
		}
		w.Write(rec[:])
	}
	// link block: every non-free slot (dying slots have A=0, lem>0)
	var nl uint32
	for sl := 0; sl < s.NSLOT; sl++ {
		if s.sst[sl] != sFree {
			nl++
		}
	}
	pu32(nl)
	var lrec [6 * 4]byte
	for sl := 0; sl < s.NSLOT; sl++ {
		if s.sst[sl] == sFree {
			continue
		}
		i, j := s.sli[sl], s.slj[sl]
		gf := s.gateOf(float64(s.slq[sl])*s.th2[i] - float64(s.slq[sl])*s.w2e[i]*s.sd[sl]/s.P.C - float64(s.slp[sl])*s.th2[j])
		gb := s.gateOf(float64(s.slp[sl])*s.th2[j] - float64(s.slp[sl])*s.w2e[j]*s.sd[sl]/s.P.C - float64(s.slq[sl])*s.th2[i])
		binary.LittleEndian.PutUint32(lrec[0:], uint32(i))
		binary.LittleEndian.PutUint32(lrec[4:], uint32(j))
		binary.LittleEndian.PutUint32(lrec[8:], math.Float32bits(float32(s.sd[sl])))
		binary.LittleEndian.PutUint32(lrec[12:], math.Float32bits(float32(s.sA[sl])))
		binary.LittleEndian.PutUint32(lrec[16:], math.Float32bits(float32(s.slem[2*sl]+s.slem[2*sl+1])))
		binary.LittleEndian.PutUint32(lrec[20:], math.Float32bits(float32(gf*gb)))
		w.Write(lrec[:])
	}
	if err := w.Flush(); err != nil {
		fprintf(s.Errw, "# WARN snap flush: %v\n", err)
	}
	f.Close()
}
