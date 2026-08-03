package fab

// .fcs snapshot writer — v90 extension (the C parser carries
// snap_every/snap_dir; the writer lands here and in the v90 C kernel).
// Format FCS1, little-endian:
//   magic  [4]byte "FCS1"
//   ver    uint32 (=1)
//   ncells uint32
//   t      float64
//   L      float64
//   per cell, 9 x float32: x y z r Es Em Ee xload tag
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
	w.WriteString("FCS1")
	var u32 [4]byte
	binary.LittleEndian.PutUint32(u32[:], 1)
	w.Write(u32[:])
	binary.LittleEndian.PutUint32(u32[:], uint32(s.NC))
	w.Write(u32[:])
	var u64 [8]byte
	binary.LittleEndian.PutUint64(u64[:], math.Float64bits(s.simT))
	w.Write(u64[:])
	binary.LittleEndian.PutUint64(u64[:], math.Float64bits(s.P.L))
	w.Write(u64[:])
	var rec [9 * 4]byte
	for i := 0; i < s.NC; i++ {
		xload := (s.Em[i] + s.flload[i]) / s.P.Cap
		vals := [9]float32{
			float32(s.px[i]), float32(s.py[i]), float32(s.pz[i]),
			float32(s.cr[i]), float32(s.Es[i]), float32(s.Em[i]),
			float32(s.Ee[i]), float32(xload), float32(s.tag[i]),
		}
		for k, v := range vals {
			binary.LittleEndian.PutUint32(rec[4*k:], math.Float32bits(v))
		}
		w.Write(rec[:])
	}
	if err := w.Flush(); err != nil {
		fprintf(s.Errw, "# WARN snap flush: %v\n", err)
	}
	f.Close()
}
