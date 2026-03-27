// volview3 — GPU-accelerated volume ray marcher for SFA files (OpenGL 4.3 + GLFW)
//
// Architecture:
//   1. Load SFA file, decompress frames into float32 arrays
//   2. Compute visualization channels (R=|P|, G=phi^2, B=theta^2)
//   3. Upload as 3D RGBA32F texture to GPU via glTexImage3D
//   4. Full-screen triangle with fragment shader does ray marching
//   5. Hardware trilinear sampling, alpha compositing, transfer function
//
// Build: cd sfa/viewer/volview3 && go build -o volview3 .
// Usage: ./volview3 archive.sfa
package main

import (
	"encoding/binary"
	"flag"
	"fmt"
	"image"
	"io"
	"log"
	"math"
	"net/http"
	_ "net/http/pprof"
	"os"
	"runtime"
	"runtime/pprof"
	"strings"
	"sync"
	"sync/atomic"
	"syscall"
	"time"
	"unsafe"

	"github.com/HugoSmits86/nativewebp"
	"github.com/go-gl/gl/v4.3-core/gl"
	"github.com/go-gl/glfw/v3.3/glfw"
	"github.com/klauspost/compress/zstd"
)

func init() {
	// GLFW event handling must run on the main thread
	runtime.LockOSThread()
}

// ============================================================
// SFA Reader (minimal, pure Go, no CGO)
// ============================================================

const (
	sfaCodecRaw     = 0
	sfaCodecZstd    = 1
	sfaCodecBSS     = 2
	sfaCodecF32BSS  = 3
	sfaCodecF16BSS  = 4
	sfaCodecColZstd = 7
	sfaCodecBQ8    = 5

	sfaDtypeF16 = 0
	sfaDtypeF32 = 1
	sfaDtypeF64 = 2
)

var sfaDtypeSize = [12]int{2, 4, 8, 16, 1, 2, 4, 8, 1, 2, 4, 8}

type sfaColumn struct {
	name      string
	dtype     uint8
	semantic  uint8
	component uint8
	flags     uint8
	scale     float64
}

type sfaL1Entry struct {
	jmpfOffset uint64
	firstFrame uint32
	frameCount uint32
}

type sfaL2Entry struct {
	time           float64
	offset         uint64
	compressedSize uint64
	checksum       uint32
	reserved       uint32
}

type sfaFile struct {
	fp *os.File

	// Header
	version        uint32
	flags          uint32
	Nx, Ny, Nz     uint32
	Lx, Ly, Lz     float64
	dt             float64
	nColumns       uint32
	totalFrames    uint32
	firstJtopOff   uint64
	cdefOffset     uint64
	jtopMax        uint32
	jmpfMax        uint32
	nTotal         uint64
	frameBytes     uint64

	columns []sfaColumn

	// mmap'd file data (nil if mmap failed — falls back to read)
	mmapData []byte

	// Pre-allocated buffers (allocated once at open, reused for every frame)
	compBuf    []byte      // buffer for compressed frame data (only used if no mmap)
	rawBuf     []byte      // buffer for decompressed frame data
	colF32Bufs [][]float32 // per-column float32 buffers, indexed by column
	colBSSBufs [][]byte    // per-column BSS temp buffers for decompression

	// Intermediate buffers for single-pass computeFieldView
	absPBuf  []float32 // |P| per voxel (cached between passes)
	phi2Buf  []float32 // Σφ² per voxel
	theta2Buf []float32 // Σθ² per voxel
}

func sfaOpen(path string) (*sfaFile, error) {
	fp, err := os.Open(path)
	if err != nil {
		return nil, err
	}

	s := &sfaFile{fp: fp}

	// Read chunk header: type[4] + size[8]
	var chunkType [4]byte
	var chunkSize uint64
	binary.Read(fp, binary.LittleEndian, &chunkType)
	binary.Read(fp, binary.LittleEndian, &chunkSize)
	if string(chunkType[:]) != "SFAH" {
		fp.Close()
		return nil, fmt.Errorf("not an SFA file (got %q)", chunkType)
	}

	// SFAH payload
	binary.Read(fp, binary.LittleEndian, &s.version)
	binary.Read(fp, binary.LittleEndian, &s.flags)
	binary.Read(fp, binary.LittleEndian, &s.Nx)
	binary.Read(fp, binary.LittleEndian, &s.Ny)
	binary.Read(fp, binary.LittleEndian, &s.Nz)
	binary.Read(fp, binary.LittleEndian, &s.Lx)
	binary.Read(fp, binary.LittleEndian, &s.Ly)
	binary.Read(fp, binary.LittleEndian, &s.Lz)
	binary.Read(fp, binary.LittleEndian, &s.dt)
	binary.Read(fp, binary.LittleEndian, &s.nColumns)
	binary.Read(fp, binary.LittleEndian, &s.totalFrames)
	binary.Read(fp, binary.LittleEndian, &s.firstJtopOff)
	binary.Read(fp, binary.LittleEndian, &s.cdefOffset)
	binary.Read(fp, binary.LittleEndian, &s.jtopMax)
	binary.Read(fp, binary.LittleEndian, &s.jmpfMax)

	// Read CDEF
	fp.Seek(int64(s.cdefOffset), io.SeekStart)
	binary.Read(fp, binary.LittleEndian, &chunkType)
	binary.Read(fp, binary.LittleEndian, &chunkSize)

	s.columns = make([]sfaColumn, s.nColumns)
	for c := uint32(0); c < s.nColumns; c++ {
		var nameBuf [12]byte
		fp.Read(nameBuf[:])
		// Trim null bytes
		nameEnd := 12
		for i := 0; i < 12; i++ {
			if nameBuf[i] == 0 {
				nameEnd = i
				break
			}
		}
		s.columns[c].name = string(nameBuf[:nameEnd])
		binary.Read(fp, binary.LittleEndian, &s.columns[c].dtype)
		binary.Read(fp, binary.LittleEndian, &s.columns[c].semantic)
		binary.Read(fp, binary.LittleEndian, &s.columns[c].component)
		binary.Read(fp, binary.LittleEndian, &s.columns[c].flags)
		binary.Read(fp, binary.LittleEndian, &s.columns[c].scale)
	}

	// Compute frame bytes
	s.nTotal = uint64(s.Nx) * uint64(s.Ny) * uint64(s.Nz)
	s.frameBytes = 0
	for c := uint32(0); c < s.nColumns; c++ {
		s.frameBytes += s.nTotal * uint64(sfaDtypeSize[s.columns[c].dtype])
	}

	// mmap the file for zero-copy access to compressed frame data
	fi, err := fp.Stat()
	if err == nil {
		data, err := syscall.Mmap(int(fp.Fd()), 0, int(fi.Size()),
			syscall.PROT_READ, syscall.MAP_PRIVATE)
		if err == nil {
			s.mmapData = data
			fmt.Printf("  mmap: %.1f GB mapped\n", float64(len(data))/1e9)
		}
	}

	// Pre-allocate reusable buffers
	s.rawBuf = make([]byte, s.frameBytes)
	s.colF32Bufs = make([][]float32, s.nColumns)
	s.colBSSBufs = make([][]byte, s.nColumns)
	for c := uint32(0); c < s.nColumns; c++ {
		s.colF32Bufs[c] = make([]float32, s.nTotal)
		colBytes := s.nTotal * uint64(sfaDtypeSize[s.columns[c].dtype])
		s.colBSSBufs[c] = make([]byte, colBytes)
	}

	// Intermediate buffers for single-pass field view computation
	s.absPBuf = make([]float32, s.nTotal)
	s.phi2Buf = make([]float32, s.nTotal)
	s.theta2Buf = make([]float32, s.nTotal)

	return s, nil
}

func (s *sfaFile) close() {
	if s.mmapData != nil {
		syscall.Munmap(s.mmapData)
		s.mmapData = nil
	}
	if s.fp != nil {
		s.fp.Close()
	}
}

// findFrame locates a frame's L2 entry in the two-level index
func (s *sfaFile) findFrame(frameIdx uint32) (*sfaL2Entry, error) {
	jtopOff := s.firstJtopOff
	fidx := frameIdx

	for jtopOff != 0 {
		// Skip chunk header (12 bytes)
		s.fp.Seek(int64(jtopOff)+12, io.SeekStart)

		var maxEnt, curEnt uint32
		var nextJtop uint64
		binary.Read(s.fp, binary.LittleEndian, &maxEnt)
		binary.Read(s.fp, binary.LittleEndian, &curEnt)
		binary.Read(s.fp, binary.LittleEndian, &nextJtop)

		l1Idx := fidx / s.jmpfMax
		if l1Idx < curEnt {
			// Read JTOP entry: jtop_offset + 12 + 16 + l1_idx * 16
			s.fp.Seek(int64(jtopOff)+12+16+int64(l1Idx)*16, io.SeekStart)
			var l1 sfaL1Entry
			binary.Read(s.fp, binary.LittleEndian, &l1.jmpfOffset)
			binary.Read(s.fp, binary.LittleEndian, &l1.firstFrame)
			binary.Read(s.fp, binary.LittleEndian, &l1.frameCount)

			l2Idx := fidx % s.jmpfMax
			if l2Idx >= l1.frameCount {
				return nil, fmt.Errorf("frame %d: l2_idx %d >= frame_count %d", frameIdx, l2Idx, l1.frameCount)
			}

			// Read JMPF entry: jmpf_offset + 12 + 8 + l2_idx * 32
			s.fp.Seek(int64(l1.jmpfOffset)+12+8+int64(l2Idx)*32, io.SeekStart)
			var entry sfaL2Entry
			binary.Read(s.fp, binary.LittleEndian, &entry.time)
			binary.Read(s.fp, binary.LittleEndian, &entry.offset)
			binary.Read(s.fp, binary.LittleEndian, &entry.compressedSize)
			binary.Read(s.fp, binary.LittleEndian, &entry.checksum)
			binary.Read(s.fp, binary.LittleEndian, &entry.reserved)
			return &entry, nil
		}

		fidx -= curEnt * s.jmpfMax
		jtopOff = nextJtop
	}
	return nil, fmt.Errorf("frame %d not found in index", frameIdx)
}

// bssDecode reverses byte-stream split for a column
func bssDecode(src []byte, nValues uint64, elemSize int) []byte {
	dst := make([]byte, nValues*uint64(elemSize))
	for b := 0; b < elemSize; b++ {
		for i := uint64(0); i < nValues; i++ {
			dst[i*uint64(elemSize)+uint64(b)] = src[uint64(b)*nValues+i]
		}
	}
	return dst
}

// bssDecodeFrame reverses BSS for an entire frame (all columns)
func (s *sfaFile) bssDecodeFrame(bss []byte) []byte {
	raw := make([]byte, s.frameBytes)
	off := uint64(0)
	for c := uint32(0); c < s.nColumns; c++ {
		es := sfaDtypeSize[s.columns[c].dtype]
		colBytes := s.nTotal * uint64(es)
		decoded := bssDecode(bss[off:off+colBytes], s.nTotal, es)
		copy(raw[off:], decoded)
		off += colBytes
	}
	return raw
}

// readFrame reads and decompresses a frame into s.rawBuf (zero-alloc after init).
func (s *sfaFile) readFrame(frameIdx uint32) (float64, error) {
	entry, err := s.findFrame(frameIdx)
	if err != nil {
		return 0, err
	}

	// Get compressed data: mmap (zero-copy) or fallback to read
	var compressed []byte
	dataOff := entry.offset + 12 // skip FRMD chunk header
	dataEnd := dataOff + entry.compressedSize

	if s.mmapData != nil && dataEnd <= uint64(len(s.mmapData)) {
		// Zero-copy: slice directly from mmap'd region
		compressed = s.mmapData[dataOff:dataEnd]
	} else {
		// Fallback: read into compBuf
		if uint64(cap(s.compBuf)) < entry.compressedSize {
			s.compBuf = make([]byte, entry.compressedSize)
		}
		compressed = s.compBuf[:entry.compressedSize]
		s.fp.Seek(int64(dataOff), io.SeekStart)
		if _, err := io.ReadFull(s.fp, compressed); err != nil {
			return 0, fmt.Errorf("read compressed data: %w", err)
		}
	}

	codec := s.flags & 0xF

	switch codec {
	case sfaCodecRaw:
		copy(s.rawBuf, compressed)
	case sfaCodecZstd:
		if err := zstdDecompressInto(compressed, s.rawBuf); err != nil {
			return 0, fmt.Errorf("zstd decompress: %w", err)
		}
	case sfaCodecColZstd:
		if err := s.readFrameColZstdInto(compressed); err != nil {
			return 0, fmt.Errorf("colzstd decompress: %w", err)
		}
	case sfaCodecBSS, sfaCodecF32BSS, sfaCodecF16BSS:
		if err := zstdDecompressInto(compressed, s.rawBuf); err != nil {
			return 0, fmt.Errorf("zstd decompress: %w", err)
		}
		// BSS decode in-place (need a temp buffer)
		if uint64(cap(s.compBuf)) < s.frameBytes {
			s.compBuf = make([]byte, s.frameBytes)
		}
		bssTmp := s.compBuf[:s.frameBytes]
		copy(bssTmp, s.rawBuf)
		off := uint64(0)
		for c := uint32(0); c < s.nColumns; c++ {
			es := sfaDtypeSize[s.columns[c].dtype]
			colBytes := s.nTotal * uint64(es)
			for b := 0; b < es; b++ {
				for v := uint64(0); v < s.nTotal; v++ {
					s.rawBuf[off+v*uint64(es)+uint64(b)] = bssTmp[off+uint64(b)*s.nTotal+v]
				}
			}
			off += colBytes
		}
	default:
		return 0, fmt.Errorf("unsupported codec %d", codec)
	}

	return entry.time, nil
}

// readFrameColZstdInto decompresses COLZSTD with parallel goroutines.
// For f32 columns: fuses BSS decode + float32 extraction directly into colF32Bufs (skips rawBuf).
// For other dtypes: BSS decodes into rawBuf using pre-allocated colBSSBufs.
func (s *sfaFile) readFrameColZstdInto(data []byte) error {
	if len(data) < 4 {
		return fmt.Errorf("colzstd data too short")
	}
	nc := int(binary.LittleEndian.Uint32(data[:4]))
	if nc != int(s.nColumns) {
		return fmt.Errorf("colzstd n_cols mismatch (%d vs %d)", nc, s.nColumns)
	}

	headerSize := 4 + nc*8

	// Parse per-column compressed offsets (stack-allocated for nc<=12)
	type colInfo struct {
		compOff, compSize int
		rawOff            uint64
	}
	var colsBuf [16]colInfo
	cols := colsBuf[:nc]

	compOff := headerSize
	var rawOff uint64
	for c := 0; c < nc; c++ {
		sz := int(binary.LittleEndian.Uint64(data[4+c*8 : 4+c*8+8]))
		es := int(sfaDtypeSize[s.columns[c].dtype])
		cols[c] = colInfo{compOff, sz, rawOff}
		compOff += sz
		rawOff += s.nTotal * uint64(es)
	}

	var wg sync.WaitGroup
	var errs [16]error

	nVals := int(s.nTotal)
	for c := 0; c < nc; c++ {
		wg.Add(1)
		go func(c int) {
			defer wg.Done()
			ci := cols[c]
			src := data[ci.compOff : ci.compOff+ci.compSize]
			es := int(sfaDtypeSize[s.columns[c].dtype])
			bssBuf := s.colBSSBufs[c] // pre-allocated, no alloc

			// Decompress into pre-allocated BSS buffer
			if err := zstdDecompressInto(src, bssBuf); err != nil {
				errs[c] = fmt.Errorf("col %d: %w", c, err)
				return
			}

			if es == 4 && s.columns[c].dtype == sfaDtypeF32 {
				// FAST PATH: fused BSS decode + f32 extraction
				// Combine 4 byte streams directly into float32 output
				// bssBuf layout: [b0_v0..b0_vN, b1_v0..b1_vN, b2_v0..b2_vN, b3_v0..b3_vN]
				out := s.colF32Bufs[c]
				b0 := bssBuf[0*nVals : 1*nVals]
				b1 := bssBuf[1*nVals : 2*nVals]
				b2 := bssBuf[2*nVals : 3*nVals]
				b3 := bssBuf[3*nVals : 4*nVals]

				// Process in blocks for cache efficiency
				const blockSize = 64
				fullBlocks := nVals / blockSize
				for blk := 0; blk < fullBlocks; blk++ {
					base := blk * blockSize
					for v := 0; v < blockSize; v++ {
						idx := base + v
						bits := uint32(b0[idx]) | uint32(b1[idx])<<8 |
							uint32(b2[idx])<<16 | uint32(b3[idx])<<24
						out[idx] = math.Float32frombits(bits)
					}
				}
				// Handle remainder
				for v := fullBlocks * blockSize; v < nVals; v++ {
					bits := uint32(b0[v]) | uint32(b1[v])<<8 |
						uint32(b2[v])<<16 | uint32(b3[v])<<24
					out[v] = math.Float32frombits(bits)
				}
			} else {
				// Generic path: BSS decode into rawBuf
				dst := s.rawBuf[ci.rawOff : ci.rawOff+uint64(nVals*es)]
				for b := 0; b < es; b++ {
					bOff := b * nVals
					for v := 0; v < nVals; v++ {
						dst[v*es+b] = bssBuf[bOff+v]
					}
				}
			}
		}(c)
	}
	wg.Wait()

	for i := 0; i < nc; i++ {
		if errs[i] != nil {
			return errs[i]
		}
	}
	return nil
}

// readFrameColZstd decompresses a COLZSTD frame with parallel per-column zstd.
// Layout: [n_cols(4)] [comp_size(8)]×n_cols [col0_data...] [col1_data...]
// extractColumnF32 extracts a column into s.colF32Bufs[colIdx] (zero-alloc).
// For COLZSTD with f32 columns, the fused decode already filled colF32Bufs — this is a no-op.
func (s *sfaFile) extractColumnF32(colIdx int) []float32 {
	// If COLZSTD with f32, the fused BSS+f32 path already wrote directly to colF32Bufs
	codec := s.flags & 0xF
	if codec == sfaCodecColZstd && s.columns[colIdx].dtype == sfaDtypeF32 {
		return s.colF32Bufs[colIdx] // already populated by readFrameColZstdInto
	}

	result := s.colF32Bufs[colIdx]

	// Find offset to column in rawBuf
	off := uint64(0)
	for c := 0; c < colIdx; c++ {
		off += s.nTotal * uint64(sfaDtypeSize[s.columns[c].dtype])
	}

	dtype := s.columns[colIdx].dtype
	switch dtype {
	case sfaDtypeF64:
		for i := uint64(0); i < s.nTotal; i++ {
			bits := binary.LittleEndian.Uint64(s.rawBuf[off+i*8 : off+i*8+8])
			result[i] = float32(math.Float64frombits(bits))
		}
	case sfaDtypeF32:
		// Unsafe reinterpret: rawBuf bytes → float32 slice (zero-copy)
		src := (*[1 << 30]float32)(unsafe.Pointer(&s.rawBuf[off]))[:s.nTotal:s.nTotal]
		copy(result, src)
	case sfaDtypeF16:
		for i := uint64(0); i < s.nTotal; i++ {
			h := binary.LittleEndian.Uint16(s.rawBuf[off+i*2 : off+i*2+2])
			result[i] = f16ToF32(h)
		}
	}
	return result
}

func f16ToF32(h uint16) float32 {
	sign := h & 0x8000
	exp := (h >> 10) & 0x1F
	mant := h & 0x3FF
	if exp == 0 {
		return 0
	}
	if exp == 31 {
		if sign != 0 {
			return float32(math.Inf(-1))
		}
		return float32(math.Inf(1))
	}
	x := (uint32(sign) << 16) | (uint32(exp-15+127) << 23) | (uint32(mant) << 13)
	return math.Float32frombits(x)
}

// ============================================================
// Zstd decompressor (using klauspost/compress, pure Go)
// ============================================================

var zstdDecoder *zstd.Decoder

func initZstd() {
	var err error
	zstdDecoder, err = zstd.NewReader(nil, zstd.WithDecoderConcurrency(1))
	if err != nil {
		log.Fatalf("zstd init: %v", err)
	}
}

func zstdDecompress(src []byte, expectedSize uint64) ([]byte, error) {
	dst := make([]byte, 0, expectedSize)
	dst, err := zstdDecoder.DecodeAll(src, dst)
	if err != nil {
		return nil, err
	}
	return dst, nil
}

// zstdDecompressInto decompresses into an existing buffer (zero-alloc).
func zstdDecompressInto(src []byte, dst []byte) error {
	result, err := zstdDecoder.DecodeAll(src, dst[:0])
	if err != nil {
		return err
	}
	if len(result) != len(dst) {
		return fmt.Errorf("decompressed %d bytes, expected %d", len(result), len(dst))
	}
	return nil
}

// ============================================================
// Volume data
// ============================================================

type volumeData struct {
	n         int
	rgba      []float32 // N*N*N*4 floats (RGBA)
	loaded    bool
	absPBuf   []float32 // intermediate: |P| per voxel
	phi2Buf   []float32 // intermediate: Σφ² per voxel
	theta2Buf []float32 // intermediate: Σθ² per voxel
}

// makeTestVolume generates a test volume (two overlapping Gaussians)
func makeTestVolume(n int) *volumeData {
	vol := &volumeData{n: n, loaded: true}
	total := n * n * n
	vol.rgba = make([]float32, total*4)

	center := float32(n) / 2
	for z := 0; z < n; z++ {
		for y := 0; y < n; y++ {
			for x := 0; x < n; x++ {
				// First Gaussian: red/orange core
				dx1 := float32(x) - center*0.8
				dy1 := float32(y) - center
				dz1 := float32(z) - center
				r1 := float32(math.Sqrt(float64(dx1*dx1 + dy1*dy1 + dz1*dz1)))
				d1 := float32(math.Exp(-float64(r1*r1) / float64(center*center/6)))

				// Second Gaussian: blue torsion shell
				dx2 := float32(x) - center*1.2
				dy2 := float32(y) - center
				dz2 := float32(z) - center
				r2 := float32(math.Sqrt(float64(dx2*dx2 + dy2*dy2 + dz2*dz2)))
				shell := float32(math.Exp(-float64((r2-center*0.3)*(r2-center*0.3)) / float64(center*center/20)))

				// Binding product (red)
				rv := d1 * d1 * d1
				// Fabric phi^2 (green) - everywhere with field
				gv := d1 * d1 * (1 - rv)
				// Torsion theta^2 (blue) - shell
				bv := shell * 0.8

				idx := (z*n*n + y*n + x) * 4
				vol.rgba[idx+0] = rv
				vol.rgba[idx+1] = gv
				vol.rgba[idx+2] = bv
				vol.rgba[idx+3] = 1.0
			}
		}
	}
	return vol
}

// f32abs returns absolute value via bit manipulation (branchless, SIMD-friendly).
func f32abs(x float32) float32 {
	return math.Float32frombits(math.Float32bits(x) & 0x7FFFFFFF)
}

// f32sqrt returns sqrt of a float32.
func f32sqrt(x float32) float32 {
	return float32(math.Sqrt(float64(x)))
}

// f32min returns the smaller of two float32 values (branchless via compiler).
func f32min(a, b float32) float32 {
	if a < b {
		return a
	}
	return b
}

// computeFieldView fills RGBA volume from SFA column data.
// Uses pre-allocated intermediate buffers from sfa.absPBuf/phi2Buf/theta2Buf
// to avoid reading the input arrays twice.
func computeFieldView(n int, phi0, phi1, phi2, theta0, theta1, theta2 []float32, vol *volumeData) {
	total := n * n * n
	vol.n = n
	vol.loaded = true
	if len(vol.rgba) != total*4 {
		vol.rgba = make([]float32, total*4)
	}

	hasTheta := theta0 != nil && theta1 != nil && theta2 != nil

	nWorkers := runtime.GOMAXPROCS(0)
	if nWorkers < 1 {
		nWorkers = 1
	}
	chunkSize := (total + nWorkers - 1) / nWorkers

	// Use caller's pre-allocated intermediate buffers if available, else allocate
	absPBuf := vol.absPBuf
	phi2Buf := vol.phi2Buf
	theta2Buf := vol.theta2Buf
	if len(absPBuf) < total {
		absPBuf = make([]float32, total)
		phi2Buf = make([]float32, total)
		theta2Buf = make([]float32, total)
	}

	// Pass 1: compute intermediates AND find max (reads input arrays once)
	type maxVals struct {
		maxP, maxPhi2, maxTheta2 float32
		_pad [4]float32 // avoid false sharing between cache lines
	}
	var resultsBuf [64]maxVals
	results := resultsBuf[:nWorkers]

	var wg sync.WaitGroup
	for w := 0; w < nWorkers; w++ {
		wg.Add(1)
		go func(w int) {
			defer wg.Done()
			start := w * chunkSize
			end := start + chunkSize
			if end > total {
				end = total
			}
			var mp, mph2, mt2 float32
			ap := absPBuf[start:end]
			p2 := phi2Buf[start:end]
			t2 := theta2Buf[start:end]
			p0 := phi0[start:end]
			p1 := phi1[start:end]
			p2in := phi2[start:end]

			for i := range ap {
				v := f32abs(p0[i] * p1[i] * p2in[i])
				ap[i] = v
				if v > mp {
					mp = v
				}
				pv := p0[i]*p0[i] + p1[i]*p1[i] + p2in[i]*p2in[i]
				p2[i] = pv
				if pv > mph2 {
					mph2 = pv
				}
			}
			if hasTheta {
				t0 := theta0[start:end]
				t1 := theta1[start:end]
				t2in := theta2[start:end]
				for i := range t2 {
					tv := t0[i]*t0[i] + t1[i]*t1[i] + t2in[i]*t2in[i]
					t2[i] = tv
					if tv > mt2 {
						mt2 = tv
					}
				}
			}
			results[w] = maxVals{maxP: mp, maxPhi2: mph2, maxTheta2: mt2}
		}(w)
	}
	wg.Wait()

	var maxP, maxPhi2, maxTheta2 float32
	for _, r := range results {
		if r.maxP > maxP {
			maxP = r.maxP
		}
		if r.maxPhi2 > maxPhi2 {
			maxPhi2 = r.maxPhi2
		}
		if r.maxTheta2 > maxTheta2 {
			maxTheta2 = r.maxTheta2
		}
	}

	invNormP := float32(1)
	if maxP*0.3 > 1e-20 {
		invNormP = 1.0 / (maxP * 0.3)
	}
	invNormPhi2 := float32(1)
	if maxPhi2*0.15 > 1e-20 {
		invNormPhi2 = 1.0 / (maxPhi2 * 0.15)
	}
	invNormTheta2 := float32(1)
	if maxTheta2*0.3 > 1e-20 {
		invNormTheta2 = 1.0 / (maxTheta2 * 0.3)
	}

	// Pass 2: normalize cached intermediates → RGBA (reads only 3 small arrays, not 6 input arrays)
	rgba := vol.rgba
	for w := 0; w < nWorkers; w++ {
		wg.Add(1)
		go func(w int) {
			defer wg.Done()
			start := w * chunkSize
			end := start + chunkSize
			if end > total {
				end = total
			}
			ap := absPBuf[start:end]
			p2 := phi2Buf[start:end]
			t2 := theta2Buf[start:end]
			out := rgba[start*4 : end*4]

			for i := range ap {
				pn := ap[i] * invNormP
				// Clamp pn for the sqrt — values > 1 would make gv negative
				sqrtPn := pn
				if sqrtPn > 1 {
					sqrtPn = 1
				}
				gv := p2[i] * (1 - f32sqrt(sqrtPn)) * invNormPhi2
				bv := t2[i] * invNormTheta2

				j := i * 4
				out[j+0] = pn  // unclamped — shader handles display clamp
				out[j+1] = gv  // non-negative now (sqrtPn ≤ 1)
				out[j+2] = bv
				out[j+3] = 1.0
			}
		}(w)
	}
	wg.Wait()
}

// ============================================================
// GLSL 430 ray marching shader
// ============================================================

const vertexShaderSource = `#version 430

out vec2 uv;

void main() {
    // Full-screen triangle: 3 vertices cover [-1,1]^2
    vec2 positions[3] = vec2[3](
        vec2(-1.0, -1.0),
        vec2( 3.0, -1.0),
        vec2(-1.0,  3.0)
    );
    gl_Position = vec4(positions[gl_VertexID], 0.0, 1.0);
    uv = positions[gl_VertexID] * 0.5 + 0.5;
}
` + "\x00"

const fragmentShaderSource = `#version 430

uniform sampler3D volume;
uniform vec3 camPos;
uniform vec3 camTarget;
uniform float brightness;
uniform float opacity;
uniform uint bgWhite;
uniform vec2 winSize;
uniform float fov;
uniform vec3 channelMask;

in vec2 uv;
out vec4 fragColor;

// Ray-box intersection for unit cube [0,1]^3
vec2 intersectBox(vec3 origin, vec3 dirInv) {
    vec3 t0 = (vec3(0.0) - origin) * dirInv;
    vec3 t1 = (vec3(1.0) - origin) * dirInv;
    vec3 tminV = min(t0, t1);
    vec3 tmaxV = max(t0, t1);
    float tmin = max(max(tminV.x, tminV.y), tminV.z);
    float tmax = min(min(tmaxV.x, tmaxV.y), tmaxV.z);
    return vec2(tmin, tmax);
}

void main() {
    float aspect = winSize.x / winSize.y;

    // Camera setup
    vec3 fwd = normalize(camTarget - camPos);

    // Construct camera basis
    vec3 worldUp = vec3(0.0, 0.0, 1.0);
    if (abs(dot(fwd, worldUp)) > 0.999) {
        worldUp = vec3(0.0, 1.0, 0.0);
    }
    vec3 right = normalize(cross(fwd, worldUp));
    vec3 up = cross(right, fwd);

    // Pixel to ray direction
    float px = (uv.x * 2.0 - 1.0) * fov * aspect;
    float py = (uv.y * 2.0 - 1.0) * fov;
    vec3 rayDir = normalize(fwd + px * right + py * up);

    // Ray-box intersection
    vec3 dirInv = 1.0 / rayDir;
    vec2 tRange = intersectBox(camPos, dirInv);

    // Background color
    vec3 bg;
    if (bgWhite == 1u) {
        bg = vec3(1.0, 1.0, 1.0);
    } else {
        bg = vec3(0.05, 0.05, 0.07);
    }

    if (tRange.x >= tRange.y || tRange.y <= 0.0) {
        fragColor = vec4(bg, 1.0);
        return;
    }

    // Ray marching
    float tmin = max(tRange.x, 0.0);
    float tmax = tRange.y;
    int numSteps = 200;
    float step = (tmax - tmin) / float(numSteps);

    vec3 accumColor = vec3(0.0);
    float accumAlpha = 0.0;

    for (int i = 0; i < numSteps; i++) {
        if (accumAlpha > 0.98) {
            break;
        }

        float t = tmin + (float(i) + 0.5) * step;
        vec3 pos = camPos + t * rayDir;

        // Sample volume (hardware trilinear)
        vec4 s = texture(volume, pos);

        // Clamp to [0,1] and apply channel toggles
        float sr = clamp(s.r, 0.0, 1.0) * channelMask.r;
        float sg = clamp(s.g, 0.0, 1.0) * channelMask.g;
        float sb = clamp(s.b, 0.0, 1.0) * channelMask.b;

        // Transfer function
        float density = (sr + sg + sb) * brightness;
        float alpha = 1.0 - exp(-density * step * opacity);

        if (alpha < 0.001) {
            continue;
        }

        // Gamma-corrected color
        float cr = pow(sr, 0.6) * brightness;
        float cg = pow(sg, 0.6) * brightness;
        float cb = pow(sb, 0.6) * brightness;
        vec3 color = vec3(cr, cg, cb);

        // Front-to-back compositing
        accumColor += (1.0 - accumAlpha) * color * alpha;
        accumAlpha += (1.0 - accumAlpha) * alpha;
    }

    // Final color with background
    vec3 finalColor;
    if (bgWhite == 1u) {
        // White background: energy darkens
        finalColor = bg - accumColor;
        finalColor = max(finalColor, vec3(0.0));
    } else {
        // Dark background: energy brightens
        finalColor = accumColor + (1.0 - accumAlpha) * bg;
    }

    fragColor = vec4(clamp(finalColor, vec3(0.0), vec3(1.0)), 1.0);
}
` + "\x00"

// ============================================================
// Camera / view state
// ============================================================

type viewParams struct {
	// Spherical camera
	azimuth   float32 // phi (horizontal rotation)
	elevation float32 // theta (vertical angle from pole)
	distance  float32

	camPos    [3]float32
	camTarget [3]float32

	// Rendering
	brightness float32
	opacity    float32
	fov        float32
	bgWhite    bool
	showR      bool
	showG      bool
	showB      bool

	// State
	viewMode   int
	curFrame   int
	playing    bool
	dragging   bool
	lastMouseX float64
	lastMouseY float64
}

func newViewParams() *viewParams {
	p := &viewParams{
		azimuth:   0.3,
		elevation: 0.7,
		distance:  2.5,
		camTarget: [3]float32{0.5, 0.5, 0.5},
		brightness: 3.0,
		opacity:    2.0,
		fov:        0.8,
		showR:      true,
		showG:      true,
		showB:      true,
	}
	p.updateCamera()
	return p
}

func (p *viewParams) updateCamera() {
	// Spherical to Cartesian (Z-up)
	sinEl := float32(math.Sin(float64(p.elevation)))
	cosEl := float32(math.Cos(float64(p.elevation)))
	sinAz := float32(math.Sin(float64(p.azimuth)))
	cosAz := float32(math.Cos(float64(p.azimuth)))

	p.camPos[0] = p.camTarget[0] + p.distance*sinEl*cosAz
	p.camPos[1] = p.camTarget[1] + p.distance*sinEl*sinAz
	p.camPos[2] = p.camTarget[2] + p.distance*cosEl
}

func (p *viewParams) reset() {
	p.azimuth = 0.3
	p.elevation = 0.7
	p.distance = 2.5
	p.brightness = 3.0
	p.opacity = 2.0
	p.updateCamera()
}

// ============================================================
// OpenGL helpers
// ============================================================

func compileShader(source string, shaderType uint32) (uint32, error) {
	shader := gl.CreateShader(shaderType)

	csources, free := gl.Strs(source)
	gl.ShaderSource(shader, 1, csources, nil)
	free()
	gl.CompileShader(shader)

	var status int32
	gl.GetShaderiv(shader, gl.COMPILE_STATUS, &status)
	if status == gl.FALSE {
		var logLength int32
		gl.GetShaderiv(shader, gl.INFO_LOG_LENGTH, &logLength)
		logStr := strings.Repeat("\x00", int(logLength+1))
		gl.GetShaderInfoLog(shader, logLength, nil, gl.Str(logStr))
		return 0, fmt.Errorf("shader compile error: %v", logStr)
	}

	return shader, nil
}

func linkProgram(vertexShader, fragmentShader uint32) (uint32, error) {
	program := gl.CreateProgram()
	gl.AttachShader(program, vertexShader)
	gl.AttachShader(program, fragmentShader)
	gl.LinkProgram(program)

	var status int32
	gl.GetProgramiv(program, gl.LINK_STATUS, &status)
	if status == gl.FALSE {
		var logLength int32
		gl.GetProgramiv(program, gl.INFO_LOG_LENGTH, &logLength)
		logStr := strings.Repeat("\x00", int(logLength+1))
		gl.GetProgramInfoLog(program, logLength, nil, gl.Str(logStr))
		return 0, fmt.Errorf("program link error: %v", logStr)
	}

	return program, nil
}

// ============================================================
// GPU resources
// ============================================================

type gpuResources struct {
	program uint32
	vao     uint32

	// Volume texture
	volTexture uint32
	volSize    int

	// Uniform locations
	locCamPos      int32
	locCamTarget   int32
	locBrightness  int32
	locOpacity     int32
	locBgWhite     int32
	locWinSize     int32
	locFov         int32
	locChannelMask int32
	locVolume      int32
}

func initGPUResources() (*gpuResources, error) {
	g := &gpuResources{}

	// Compile shaders
	vertShader, err := compileShader(vertexShaderSource, gl.VERTEX_SHADER)
	if err != nil {
		return nil, fmt.Errorf("vertex shader: %w", err)
	}
	fragShader, err := compileShader(fragmentShaderSource, gl.FRAGMENT_SHADER)
	if err != nil {
		return nil, fmt.Errorf("fragment shader: %w", err)
	}

	// Link program
	g.program, err = linkProgram(vertShader, fragShader)
	if err != nil {
		return nil, fmt.Errorf("link: %w", err)
	}

	// Clean up shader objects
	gl.DeleteShader(vertShader)
	gl.DeleteShader(fragShader)

	// Get uniform locations
	gl.UseProgram(g.program)
	g.locCamPos = gl.GetUniformLocation(g.program, gl.Str("camPos\x00"))
	g.locCamTarget = gl.GetUniformLocation(g.program, gl.Str("camTarget\x00"))
	g.locBrightness = gl.GetUniformLocation(g.program, gl.Str("brightness\x00"))
	g.locOpacity = gl.GetUniformLocation(g.program, gl.Str("opacity\x00"))
	g.locBgWhite = gl.GetUniformLocation(g.program, gl.Str("bgWhite\x00"))
	g.locWinSize = gl.GetUniformLocation(g.program, gl.Str("winSize\x00"))
	g.locFov = gl.GetUniformLocation(g.program, gl.Str("fov\x00"))
	g.locChannelMask = gl.GetUniformLocation(g.program, gl.Str("channelMask\x00"))
	g.locVolume = gl.GetUniformLocation(g.program, gl.Str("volume\x00"))

	// Set texture unit 0 for volume sampler
	gl.Uniform1i(g.locVolume, 0)

	// Create VAO (required for core profile, even with no vertex attribs)
	gl.GenVertexArrays(1, &g.vao)

	return g, nil
}

func (g *gpuResources) createVolumeTexture(n int) {
	if g.volTexture != 0 && g.volSize == n {
		return
	}
	if g.volTexture != 0 {
		gl.DeleteTextures(1, &g.volTexture)
	}

	g.volSize = n
	gl.GenTextures(1, &g.volTexture)
	gl.ActiveTexture(gl.TEXTURE0)
	gl.BindTexture(gl.TEXTURE_3D, g.volTexture)

	// Trilinear sampling, clamp to edge
	gl.TexParameteri(gl.TEXTURE_3D, gl.TEXTURE_MIN_FILTER, gl.LINEAR)
	gl.TexParameteri(gl.TEXTURE_3D, gl.TEXTURE_MAG_FILTER, gl.LINEAR)
	gl.TexParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE)
	gl.TexParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE)
	gl.TexParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_R, gl.CLAMP_TO_EDGE)

	// Allocate texture storage
	gl.TexImage3D(gl.TEXTURE_3D, 0, gl.RGBA32F,
		int32(n), int32(n), int32(n),
		0, gl.RGBA, gl.FLOAT, nil)
}

func (g *gpuResources) uploadVolume(vol *volumeData) {
	n := int32(vol.n)
	gl.ActiveTexture(gl.TEXTURE0)
	gl.BindTexture(gl.TEXTURE_3D, g.volTexture)
	gl.TexSubImage3D(gl.TEXTURE_3D, 0,
		0, 0, 0, n, n, n,
		gl.RGBA, gl.FLOAT, unsafe.Pointer(&vol.rgba[0]))
}

func (g *gpuResources) updateUniforms(p *viewParams, winW, winH float32) {
	gl.UseProgram(g.program)

	gl.Uniform3f(g.locCamPos, p.camPos[0], p.camPos[1], p.camPos[2])
	gl.Uniform3f(g.locCamTarget, p.camTarget[0], p.camTarget[1], p.camTarget[2])
	gl.Uniform1f(g.locBrightness, p.brightness)
	gl.Uniform1f(g.locOpacity, p.opacity)

	var bgWhite uint32
	if p.bgWhite {
		bgWhite = 1
	}
	gl.Uniform1ui(g.locBgWhite, bgWhite)
	gl.Uniform2f(g.locWinSize, winW, winH)
	gl.Uniform1f(g.locFov, p.fov)

	var sr, sg, sb float32
	if p.showR {
		sr = 1
	}
	if p.showG {
		sg = 1
	}
	if p.showB {
		sb = 1
	}
	gl.Uniform3f(g.locChannelMask, sr, sg, sb)
}

func (g *gpuResources) render() {
	gl.BindVertexArray(g.vao)
	gl.DrawArrays(gl.TRIANGLES, 0, 3)
}

// ============================================================
// Async frame loader (copied from volview2)
// ============================================================

type frameLoader struct {
	sfa    *sfaFile
	volN   int
	mu     sync.Mutex
	ready  *volumeData // non-nil when a new frame is ready for upload
	readyT float64     // time of the ready frame
	readyF int         // frame index of the ready frame
	busy   atomic.Bool // true while loading
	workVol volumeData // pre-allocated work buffer for computing volume
}

func (fl *frameLoader) isLoading() bool {
	return fl.busy.Load()
}

// startLoad kicks off async frame loading. Non-blocking.
func (fl *frameLoader) startLoad(frameIdx int) {
	if fl.busy.Load() {
		return // already loading
	}
	fl.busy.Store(true)
	fmt.Printf("Loading frame %d/%d...\n", frameIdx+1, fl.sfa.totalFrames)
	go func() {
		defer fl.busy.Store(false)

		frameTime, err := fl.sfa.readFrame(uint32(frameIdx))
		if err != nil {
			fmt.Printf("Failed to load frame %d: %v\n", frameIdx, err)
			return
		}

		// Extract only the columns we need (lazy: skip velocity/acceleration)
		phi0 := fl.sfa.extractColumnF32(0)
		phi1 := fl.sfa.extractColumnF32(1)
		phi2 := fl.sfa.extractColumnF32(2)
		var theta0, theta1, theta2 []float32
		if fl.sfa.nColumns >= 6 {
			theta0 = fl.sfa.extractColumnF32(3)
			theta1 = fl.sfa.extractColumnF32(4)
			theta2 = fl.sfa.extractColumnF32(5)
		}

		computeFieldView(fl.volN, phi0, phi1, phi2, theta0, theta1, theta2, &fl.workVol)

		fl.mu.Lock()
		fl.ready = &fl.workVol
		fl.readyT = frameTime
		fl.readyF = frameIdx
		fl.mu.Unlock()
	}()
}

// consume returns the loaded volume if ready, or nil.
func (fl *frameLoader) consume() (*volumeData, float64, int, bool) {
	fl.mu.Lock()
	defer fl.mu.Unlock()
	if fl.ready == nil {
		return nil, 0, 0, false
	}
	v := fl.ready
	t := fl.readyT
	f := fl.readyF
	fl.ready = nil
	return v, t, f, true
}

func navigateFrame(fl *frameLoader, params *viewParams, delta int) {
	if fl.isLoading() {
		return
	}
	next := params.curFrame + delta
	if next < 0 {
		next = 0
	}
	if next >= int(fl.sfa.totalFrames) {
		next = int(fl.sfa.totalFrames) - 1
	}
	if next == params.curFrame {
		return
	}
	params.curFrame = next
	fl.startLoad(next)
}

func loadSFAFrame(fl *frameLoader, params *viewParams) {
	fl.startLoad(params.curFrame)
}

// ============================================================
// Screenshot + Animation
// ============================================================

// readFramebuffer reads the current OpenGL framebuffer as an NRGBA image, flipped vertically.
func readFramebuffer(x, y, w, h int) *image.NRGBA {
	pixels := make([]uint8, w*h*4)
	gl.ReadPixels(int32(x), int32(y), int32(w), int32(h), gl.RGBA, gl.UNSIGNED_BYTE, unsafe.Pointer(&pixels[0]))

	img := image.NewNRGBA(image.Rect(0, 0, w, h))
	stride := w * 4
	for row := 0; row < h; row++ {
		srcRow := h - 1 - row // flip vertically (OpenGL origin is bottom-left)
		copy(img.Pix[row*stride:(row+1)*stride], pixels[srcRow*stride:(srcRow+1)*stride])
	}
	return img
}

// saveScreenshot captures the current framebuffer and writes a WebP file.
func saveScreenshot(w, h int) {
	img := readFramebuffer(0, 0, w, h)
	ts := time.Now().Format("20060102_150405")
	fname := fmt.Sprintf("screenshot_%s.webp", ts)

	f, err := os.Create(fname)
	if err != nil {
		fmt.Printf("Screenshot error: %v\n", err)
		return
	}
	defer f.Close()

	if err := nativewebp.Encode(f, img, nil); err != nil {
		fmt.Printf("WebP encode error: %v\n", err)
		return
	}
	fmt.Printf("Screenshot saved: %s (%dx%d)\n", fname, w, h)
}

// renderQuadrant renders the scene into a specific quadrant of an FBO.
// quadrant: 0=top-left, 1=top-right, 2=bottom-left, 3=bottom-right
// channelMask: RGB channel visibility
// bgWhite: use white background
func renderQuadrant(gpu *gpuResources, params *viewParams, qx, qy, qw, qh int, maskR, maskG, maskB float32, bgW bool) {
	gl.Viewport(int32(qx), int32(qy), int32(qw), int32(qh))
	gl.UseProgram(gpu.program)

	gl.Uniform3f(gpu.locCamPos, params.camPos[0], params.camPos[1], params.camPos[2])
	gl.Uniform3f(gpu.locCamTarget, params.camTarget[0], params.camTarget[1], params.camTarget[2])
	gl.Uniform1f(gpu.locBrightness, params.brightness)
	gl.Uniform1f(gpu.locOpacity, params.opacity)
	gl.Uniform1f(gpu.locFov, params.fov)
	gl.Uniform2f(gpu.locWinSize, float32(qw), float32(qh))
	gl.Uniform3f(gpu.locChannelMask, maskR, maskG, maskB)

	var bg uint32
	if bgW {
		bg = 1
	}
	gl.Uniform1ui(gpu.locBgWhite, bg)

	gl.BindVertexArray(gpu.vao)
	gl.DrawArrays(gl.TRIANGLES, 0, 3)
}

// exportAnimation renders a 2x2 composite for each SFA frame and saves as animated WebP.
func exportAnimation(gpu *gpuResources, params *viewParams, fl *frameLoader, sfa *sfaFile, vol *volumeData, winW, winH int) {
	startFrame := params.curFrame
	endFrame := int(sfa.totalFrames)
	numFrames := endFrame - startFrame
	if numFrames <= 0 {
		fmt.Println("No frames to export (already at last frame)")
		return
	}

	// Composite is 2x the window size (2x2 grid of half-res quadrants = same total as 2x2*win)
	compW := winW * 2
	compH := winH * 2
	halfW := winW
	halfH := winH

	// Create offscreen FBO with color texture
	var fbo, fboTex uint32
	gl.GenFramebuffers(1, &fbo)
	gl.BindFramebuffer(gl.FRAMEBUFFER, fbo)

	gl.GenTextures(1, &fboTex)
	gl.BindTexture(gl.TEXTURE_2D, fboTex)
	gl.TexImage2D(gl.TEXTURE_2D, 0, gl.RGBA8, int32(compW), int32(compH), 0, gl.RGBA, gl.UNSIGNED_BYTE, nil)
	gl.TexParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR)
	gl.TexParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR)
	gl.FramebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, fboTex, 0)

	status := gl.CheckFramebufferStatus(gl.FRAMEBUFFER)
	if status != gl.FRAMEBUFFER_COMPLETE {
		fmt.Printf("FBO not complete: 0x%x\n", status)
		gl.BindFramebuffer(gl.FRAMEBUFFER, 0)
		gl.DeleteTextures(1, &fboTex)
		gl.DeleteFramebuffers(1, &fbo)
		return
	}

	ts := time.Now().Format("20060102_150405")
	fname := fmt.Sprintf("animation_%s.webp", ts)

	images := make([]image.Image, 0, numFrames)
	durations := make([]uint, 0, numFrames)
	disposals := make([]uint, 0, numFrames)

	fmt.Printf("Recording... %d frames to %s\n", numFrames, fname)

	// Bind the volume texture (stays the same throughout)
	gl.ActiveTexture(gl.TEXTURE0)
	gl.BindTexture(gl.TEXTURE_3D, gpu.volTexture)

	for fi := startFrame; fi < endFrame; fi++ {
		// Load frame synchronously
		frameTime, err := sfa.readFrame(uint32(fi))
		if err != nil {
			fmt.Printf("  Frame %d load error: %v\n", fi, err)
			break
		}

		// Extract columns and compute field view
		phi0 := sfa.extractColumnF32(0)
		phi1 := sfa.extractColumnF32(1)
		phi2 := sfa.extractColumnF32(2)
		var theta0, theta1, theta2 []float32
		if sfa.nColumns >= 6 {
			theta0 = sfa.extractColumnF32(3)
			theta1 = sfa.extractColumnF32(4)
			theta2 = sfa.extractColumnF32(5)
		}
		computeFieldView(int(sfa.Nx), phi0, phi1, phi2, theta0, theta1, theta2, vol)

		// Upload new volume data
		gpu.uploadVolume(vol)

		// Render 4 quadrants into the FBO
		gl.BindFramebuffer(gl.FRAMEBUFFER, fbo)
		gl.Clear(gl.COLOR_BUFFER_BIT)

		// Top-left: Field view (R=|P|, G=phi^2, B=theta^2) -- full color
		renderQuadrant(gpu, params, 0, halfH, halfW, halfH, 1, 1, 1, false)

		// Top-right: RB view (R=|P|, B=theta^2, no green) -- binding + EM
		renderQuadrant(gpu, params, halfW, halfH, halfW, halfH, 1, 0, 1, false)

		// Bottom-left: White background contrast view
		renderQuadrant(gpu, params, 0, 0, halfW, halfH, 1, 1, 1, true)

		// Bottom-right: Green only (phi^2 -- the fabric)
		renderQuadrant(gpu, params, halfW, 0, halfW, halfH, 0, 1, 0, false)

		// Read back the composite
		img := readFramebuffer(0, 0, compW, compH)

		images = append(images, img)
		durations = append(durations, 200)
		disposals = append(disposals, 1)

		fmt.Printf("  Frame %d/%d (t=%.3f)\n", fi-startFrame+1, numFrames, frameTime)
	}

	// Restore default framebuffer
	gl.BindFramebuffer(gl.FRAMEBUFFER, 0)
	gl.DeleteTextures(1, &fboTex)
	gl.DeleteFramebuffers(1, &fbo)

	if len(images) == 0 {
		fmt.Println("No frames captured")
		return
	}

	// Encode animated WebP
	fmt.Printf("Encoding %d frames as animated WebP...\n", len(images))

	f, err := os.Create(fname)
	if err != nil {
		fmt.Printf("Animation file error: %v\n", err)
		return
	}
	defer f.Close()

	ani := &nativewebp.Animation{
		Images:          images,
		Durations:       durations,
		Disposals:       disposals,
		LoopCount:       0, // infinite loop
		BackgroundColor: 0xFF0D0D12, // dark bg in BGRA
	}
	if err := nativewebp.EncodeAll(f, ani, nil); err != nil {
		fmt.Printf("Animated WebP encode error: %v\n", err)
		return
	}

	fmt.Printf("Animation saved: %s (%dx%d, %d frames)\n", fname, compW, compH, len(images))
}

// ============================================================
// Main
// ============================================================

const (
	defaultWinW = 900
	defaultWinH = 700
)

func printHelp() {
	fmt.Println()
	fmt.Println("=== volview3 Controls ===")
	fmt.Println("  Mouse drag     Rotate camera")
	fmt.Println("  Scroll         Zoom in/out")
	fmt.Println("  Left/Right,A/D Previous/next frame")
	fmt.Println("  Home / End     First/last frame")
	fmt.Println("  Space          Play/pause animation")
	fmt.Println("  1 / 2 / 3     Toggle Red/Green/Blue channel")
	fmt.Println("  + / -          Brightness up/down")
	fmt.Println("  O / P          Opacity down/up")
	fmt.Println("  B              Toggle white/dark background")
	fmt.Println("  S              Save screenshot (WebP)")
	fmt.Println("  Shift+A        Export animated WebP (2x2 composite)")
	fmt.Println("  R              Reset view")
	fmt.Println("  Escape / Q     Quit")
	fmt.Println()
	fmt.Println("  RED=|P| (binding), GREEN=phi^2 (fabric), BLUE=theta^2 (torsion)")
	fmt.Println()
}

func main() {
	fmt.Println("volview3: GPU volume ray marcher (OpenGL 4.3)")

	// Parse flags (before positional args)
	cpuprofile := flag.String("cpuprofile", "", "write CPU profile to file")
	flag.Parse()

	// CPU profiling to file
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatalf("cpuprofile: %v", err)
		}
		defer f.Close()
		if err := pprof.StartCPUProfile(f); err != nil {
			log.Fatalf("cpuprofile start: %v", err)
		}
		defer pprof.StopCPUProfile()
		fmt.Printf("CPU profile: writing to %s\n", *cpuprofile)
	}

	// HTTP pprof server
	go func() {
		fmt.Println("Profile: http://localhost:6060/debug/pprof/")
		http.ListenAndServe("localhost:6060", nil)
	}()

	initZstd()

	// Parse positional arg (SFA path) from remaining args after flags
	args := flag.Args()
	var sfaPath string
	if len(args) > 0 {
		sfaPath = args[0]
	}

	// Load SFA if provided
	var sfa *sfaFile
	if sfaPath != "" {
		var err error
		sfa, err = sfaOpen(sfaPath)
		if err != nil {
			log.Fatalf("Failed to open SFA: %v", err)
		}
		defer sfa.close()
		fmt.Printf("SFA: %dx%dx%d, %d columns, %d frames\n",
			sfa.Nx, sfa.Ny, sfa.Nz, sfa.nColumns, sfa.totalFrames)
		for i, col := range sfa.columns {
			fmt.Printf("  [%d] %s (dtype=%d)\n", i, col.name, col.dtype)
		}
	}

	if sfa == nil {
		fmt.Fprintf(os.Stderr, "Usage: volview3 [-cpuprofile file] <file.sfa>\n")
		os.Exit(1)
	}

	volN := int(sfa.Nx)
	vol := &volumeData{n: volN, rgba: make([]float32, volN*volN*volN*4)}

	// Initialize GLFW
	if err := glfw.Init(); err != nil {
		log.Fatalf("GLFW init: %v", err)
	}
	defer glfw.Terminate()

	glfw.WindowHint(glfw.Resizable, glfw.True)
	glfw.WindowHint(glfw.ContextVersionMajor, 4)
	glfw.WindowHint(glfw.ContextVersionMinor, 3)
	glfw.WindowHint(glfw.OpenGLProfile, glfw.OpenGLCoreProfile)
	glfw.WindowHint(glfw.OpenGLForwardCompatible, glfw.True)

	window, err := glfw.CreateWindow(defaultWinW, defaultWinH, "volview3 -- GPU Volume Viewer", nil, nil)
	if err != nil {
		log.Fatalf("Create window: %v", err)
	}
	window.MakeContextCurrent()
	glfw.SwapInterval(1) // vsync

	// Initialize OpenGL
	if err := gl.Init(); err != nil {
		log.Fatalf("OpenGL init: %v", err)
	}
	fmt.Printf("OpenGL: %s\n", gl.GoStr(gl.GetString(gl.VERSION)))
	fmt.Printf("Renderer: %s\n", gl.GoStr(gl.GetString(gl.RENDERER)))

	// Initialize GPU resources
	gpu, err := initGPUResources()
	if err != nil {
		log.Fatalf("Init GPU: %v", err)
	}

	// Create and upload volume texture
	gpu.createVolumeTexture(volN)
	gpu.uploadVolume(vol)

	fmt.Printf("Volume: %dx%dx%d loaded\n", volN, volN, volN)
	printHelp()

	params := newViewParams()

	// Start async frame loader
	var fl *frameLoader
	if sfa != nil {
		total := volN * volN * volN
		fl = &frameLoader{sfa: sfa, volN: volN}
		fl.workVol.absPBuf = make([]float32, total)
		fl.workVol.phi2Buf = make([]float32, total)
		fl.workVol.theta2Buf = make([]float32, total)
		fl.workVol.rgba = make([]float32, total*4)
		loadSFAFrame(fl, params)
	}

	lastAnimTime := time.Now()
	animInterval := 200 * time.Millisecond

	// -- GLFW callbacks --

	// Scroll callback
	window.SetScrollCallback(func(w *glfw.Window, xoff, yoff float64) {
		if yoff > 0 {
			params.distance *= 0.9
		} else if yoff < 0 {
			params.distance *= 1.1
		}
		if params.distance < 0.5 {
			params.distance = 0.5
		}
		if params.distance > 10 {
			params.distance = 10
		}
		params.updateCamera()
	})

	// Mouse button callback
	window.SetMouseButtonCallback(func(w *glfw.Window, button glfw.MouseButton, action glfw.Action, mods glfw.ModifierKey) {
		if button == glfw.MouseButtonLeft {
			if action == glfw.Press {
				params.dragging = true
				params.lastMouseX, params.lastMouseY = w.GetCursorPos()
			} else if action == glfw.Release {
				params.dragging = false
			}
		}
	})

	// Cursor position callback (mouse drag)
	window.SetCursorPosCallback(func(w *glfw.Window, xpos, ypos float64) {
		if !params.dragging {
			return
		}
		dx := xpos - params.lastMouseX
		dy := ypos - params.lastMouseY
		if dx != 0 || dy != 0 {
			params.azimuth += float32(dx) * 0.01
			params.elevation -= float32(dy) * 0.01
			if params.elevation < 0.1 {
				params.elevation = 0.1
			}
			if params.elevation > 3.04 {
				params.elevation = 3.04
			}
			params.updateCamera()
		}
		params.lastMouseX = xpos
		params.lastMouseY = ypos
	})

	// Key callback
	window.SetKeyCallback(func(w *glfw.Window, key glfw.Key, scancode int, action glfw.Action, mods glfw.ModifierKey) {
		if action != glfw.Press {
			return
		}

		switch key {
		case glfw.KeyEscape, glfw.KeyQ:
			w.SetShouldClose(true)
		case glfw.KeyR:
			params.reset()
		case glfw.KeyB:
			params.bgWhite = !params.bgWhite
			if params.bgWhite {
				fmt.Println("Background: white")
			} else {
				fmt.Println("Background: dark")
			}
		case glfw.KeySpace:
			params.playing = !params.playing
			if params.playing {
				fmt.Println("[PLAY]")
			} else {
				fmt.Println("[PAUSE]")
			}
		case glfw.Key1:
			params.showR = !params.showR
		case glfw.Key2:
			params.showG = !params.showG
		case glfw.Key3:
			params.showB = !params.showB
		case glfw.KeyEqual: // + key
			params.brightness *= 1.3
		case glfw.KeyMinus:
			params.brightness /= 1.3
		case glfw.KeyO:
			params.opacity *= 0.7
			if params.opacity < 0.1 {
				params.opacity = 0.1
			}
		case glfw.KeyP:
			params.opacity *= 1.4
			if params.opacity > 20 {
				params.opacity = 20
			}

		// Screenshot
		case glfw.KeyS:
			fbW, fbH := w.GetFramebufferSize()
			saveScreenshot(fbW, fbH)

		// Frame navigation
		case glfw.KeyLeft:
			if fl != nil {
				navigateFrame(fl, params, -1)
			}
		case glfw.KeyA:
			if mods&glfw.ModShift != 0 {
				// Shift+A: export animated WebP
				if fl != nil {
					fbW, fbH := w.GetFramebufferSize()
					exportAnimation(gpu, params, fl, sfa, vol, fbW, fbH)
					// Reload current frame to restore display
					loadSFAFrame(fl, params)
				}
			} else {
				// Plain A: previous frame
				if fl != nil {
					navigateFrame(fl, params, -1)
				}
			}
		case glfw.KeyRight, glfw.KeyD:
			if fl != nil {
				navigateFrame(fl, params, 1)
			}
		case glfw.KeyHome:
			if fl != nil {
				params.curFrame = 0
				loadSFAFrame(fl, params)
			}
		case glfw.KeyEnd:
			if fl != nil {
				params.curFrame = int(sfa.totalFrames) - 1
				loadSFAFrame(fl, params)
			}
		}
	})

	// -- Main loop --
	for !window.ShouldClose() {
		glfw.PollEvents()

		// Check async frame loader
		if fl != nil {
			if newVol, frameTime, frameIdx, ok := fl.consume(); ok {
				// Recreate texture if size changed
				if newVol.n != gpu.volSize {
					gpu.createVolumeTexture(newVol.n)
				}
				gpu.uploadVolume(newVol)
				fmt.Printf("Frame %d/%d t=%.3f (ready)\n", frameIdx+1, sfa.totalFrames, frameTime)
			}
		}

		// Animation
		if params.playing && fl != nil && !fl.isLoading() {
			now := time.Now()
			if now.Sub(lastAnimTime) > animInterval {
				nextFrame := params.curFrame + 1
				if nextFrame >= int(sfa.totalFrames) {
					nextFrame = 0
				}
				params.curFrame = nextFrame
				loadSFAFrame(fl, params)
				lastAnimTime = now
			}
		}

		// Render
		w, h := window.GetFramebufferSize()
		gl.Viewport(0, 0, int32(w), int32(h))
		gl.Clear(gl.COLOR_BUFFER_BIT)

		gpu.updateUniforms(params, float32(w), float32(h))
		gpu.render()

		window.SwapBuffers()
	}

	fmt.Println("Exit")
}
