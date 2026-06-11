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
	frameType      uint32 // 0=voxel(FRMD), 1=vec_I(FRVD), 2=vec_P(FRVD)
}

const (
	sfaFrameVoxel = 0
	sfaFrameVecI  = 1
	sfaFrameVecP  = 2
	sfaFrameVecK  = 3
	sfaFrameMesh  = 4 // Voronoi mesh frame (FMSH chunk)
	sfaFrameCell  = 5 // per-cell field values (FCEL chunk)
	sfaFrameCellP = 6 // sparse delta on cells (FCEP chunk)
	sfaFrameTemporalModel = 7 // standalone Fourier model (FMTL chunk)
	sfaChunkFRVD  = 0x44565246 // "FRVD" little-endian
	sfaChunkFMSH  = 0x48534D46 // "FMSH"
	sfaChunkFCEL  = 0x4C454346 // "FCEL"
	sfaChunkFCEP  = 0x50454346 // "FCEP"
	sfaChunkFMTL  = 0x4C544D46 // "FMTL"
)

// Default rendering resolution when displaying a cell-native (FMSH/FCEL) SFA.
// Can be overridden with the -cell-voxel-N CLI flag.
const defaultCellVoxelN = 192

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

	// Scratch buffers for complex-modulus / gauge derived arrays (v66+ files)
	effBufs [9][]float32

	// Vector frame state
	vecCoeffs   []float32 // n_patches × nCoeffs: current patch coefficients
	vecOrigins  [][3]int16 // n_patches × (ox,oy,oz)
	vecNPatches uint32
	vecNCoeffs  uint16
	vecBS       int
	vecOrder    int
	// Temporal model for P-frame prediction
	vecTempMean     []float32
	vecTempAmp      []float32
	vecTempPhase    []float32
	vecTempOmega    float32
	vecTempValid    bool
	vecTempIFrameIdx int32 // which I-frame loaded the current model (-1 = none)

	// Cell-native (FMSH + FCEL[+FCEP]) state
	cellNative      bool      // file contains FMSH/FCEL/FCEP frames (auto-detected)
	cellVoxelN      int       // rendering-grid resolution
	meshNCells      uint32    // current mesh cell count
	meshL           float64   // current mesh box half-extent
	meshCellPos     []float64 // 3 × N_cells (interleaved x,y,z)
	meshCellVol     []float64 // N_cells
	voxelToCell     []int32   // cellVoxelN^3 → cell index (precomputed nearest-neighbour)
	cellDataBuf     []byte    // per-cell raw payload buffer (parsed FCEL)
	cellLoadedMesh  int32     // index of mesh frame whose voxelToCell is current (-1 = none)
	// Per-cell state buffer (column-major, n_columns × N_cells f32). Updated
	// by FCEL frames (full state) and FCEP frames (sparse delta from model).
	cellState       []float32
	// Temporal model loaded from an FCEL v2 I-frame; reused by FCEP frames.
	cellTOmega      float32
	cellTMean       []float32 // n_columns × N_cells
	cellTAmp        []float32
	cellTPhase      []float32
	cellTValid      bool
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
	fileSize := uint64(0)
	if err == nil {
		fileSize = uint64(fi.Size())
		data, err := syscall.Mmap(int(fp.Fd()), 0, int(fi.Size()),
			syscall.PROT_READ, syscall.MAP_PRIVATE)
		if err == nil {
			s.mmapData = data
			fmt.Printf("  mmap: %.1f GB mapped\n", float64(len(data))/1e9)
		}
	}

	// Handle streaming files: totalFrames may be 0 or stale.
	// Scan the JMPF to count actually-written frames.
	if s.totalFrames == 0 || (s.flags&0x20) != 0 { // SFA_FLAG_STREAMING = 0x20
		counted := s.countValidFrames(fileSize)
		if counted > 0 {
			fmt.Printf("  streaming: found %d valid frames (header said %d)\n", counted, s.totalFrames)
			s.totalFrames = uint32(counted)
		}
	}

	// Cell-native auto-detection: scan JMPF for any FMSH frame.
	s.cellLoadedMesh = -1
	if s.detectCellNative() {
		s.cellNative = true
		s.cellVoxelN = defaultCellVoxelN
		// Override the SFA's "voxel grid" with our rendering resolution.
		// The recorded Nx/Ny/Nz (= N_cells, 1, 1) was a placeholder.
		s.Nx = uint32(s.cellVoxelN)
		s.Ny = uint32(s.cellVoxelN)
		s.Nz = uint32(s.cellVoxelN)
		s.nTotal = uint64(s.Nx) * uint64(s.Ny) * uint64(s.Nz)
		s.frameBytes = 0
		for c := uint32(0); c < s.nColumns; c++ {
			s.frameBytes += s.nTotal * uint64(sfaDtypeSize[s.columns[c].dtype])
		}
		fmt.Printf("  cell-native (FMSH+FCEL) — rendering at %d^3 voxels\n", s.cellVoxelN)
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

// countValidFrames scans the JMPF index to count frames with valid offsets
// that fit within the file. Handles streaming files where totalFrames=0.
func (s *sfaFile) countValidFrames(fileSize uint64) int {
	// Navigate: JTOP → first L1 entry → JMPF offset
	s.fp.Seek(int64(s.firstJtopOff)+12+4+4+8, io.SeekStart) // skip chunk hdr + max + cur + next
	var jmpfOff uint64
	binary.Read(s.fp, binary.LittleEndian, &jmpfOff)

	// Read JMPF header
	s.fp.Seek(int64(jmpfOff)+12, io.SeekStart) // skip chunk header
	var jmpfMax, jmpfCur uint32
	binary.Read(s.fp, binary.LittleEndian, &jmpfMax)
	binary.Read(s.fp, binary.LittleEndian, &jmpfCur)

	// Scan entries
	count := 0
	limit := jmpfMax
	if jmpfCur > 0 && jmpfCur < limit {
		limit = jmpfCur
	}
	for i := uint32(0); i < limit; i++ {
		var ftime float64
		var foffset, fcompSize uint64
		var fcrc, freserved uint32
		binary.Read(s.fp, binary.LittleEndian, &ftime)
		binary.Read(s.fp, binary.LittleEndian, &foffset)
		binary.Read(s.fp, binary.LittleEndian, &fcompSize)
		binary.Read(s.fp, binary.LittleEndian, &fcrc)
		binary.Read(s.fp, binary.LittleEndian, &freserved)

		if foffset == 0 && fcompSize == 0 && i > 0 {
			break // empty slot
		}
		if fileSize > 0 && foffset+12+fcompSize > fileSize {
			break // truncated
		}
		count++
	}
	return count
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

// detectCellNative scans the JMPF index for any FMSH frame; returns true
// if found. The file is then handled in cell-native mode.
func (s *sfaFile) detectCellNative() bool {
	// JTOP → first L1 entry → JMPF offset
	s.fp.Seek(int64(s.firstJtopOff)+12+4+4+8, io.SeekStart)
	var jmpfOff uint64
	if err := binary.Read(s.fp, binary.LittleEndian, &jmpfOff); err != nil {
		return false
	}
	s.fp.Seek(int64(jmpfOff)+12, io.SeekStart)
	var jmpfMax, jmpfCur uint32
	binary.Read(s.fp, binary.LittleEndian, &jmpfMax)
	binary.Read(s.fp, binary.LittleEndian, &jmpfCur)
	limit := jmpfMax
	if jmpfCur > 0 && jmpfCur < limit {
		limit = jmpfCur
	}
	for i := uint32(0); i < limit; i++ {
		var ftime float64
		var foffset, fcompSize uint64
		var fcrc, ftype uint32
		binary.Read(s.fp, binary.LittleEndian, &ftime)
		binary.Read(s.fp, binary.LittleEndian, &foffset)
		binary.Read(s.fp, binary.LittleEndian, &fcompSize)
		binary.Read(s.fp, binary.LittleEndian, &fcrc)
		binary.Read(s.fp, binary.LittleEndian, &ftype)
		if foffset == 0 && fcompSize == 0 && i > 0 {
			break
		}
		if ftype == sfaFrameMesh || ftype == sfaFrameCell ||
			ftype == sfaFrameCellP || ftype == sfaFrameTemporalModel {
			return true
		}
	}
	return false
}

// readMeshPayload decompresses a FMSH chunk and parses cell positions/volumes.
func (s *sfaFile) readMeshPayload(entry sfaL2Entry) error {
	dataOff := entry.offset + 12 + 8 // skip chunk header (12) + time (8)
	var compressed []byte
	if s.mmapData != nil && dataOff+entry.compressedSize <= uint64(len(s.mmapData)) {
		compressed = s.mmapData[dataOff : dataOff+entry.compressedSize]
	} else {
		if uint64(cap(s.compBuf)) < entry.compressedSize {
			s.compBuf = make([]byte, entry.compressedSize)
		}
		compressed = s.compBuf[:entry.compressedSize]
		s.fp.Seek(int64(dataOff), io.SeekStart)
		if _, err := io.ReadFull(s.fp, compressed); err != nil {
			return fmt.Errorf("read FMSH: %w", err)
		}
	}
	// Decompress to a fresh buffer (mesh payload may be large; allocate once)
	dec, err := zstdDecompress(compressed, uint64(len(compressed)*4))
	if err != nil {
		return fmt.Errorf("decompress FMSH: %w", err)
	}
	if len(dec) < 28 {
		return fmt.Errorf("FMSH payload too short")
	}
	magic := binary.LittleEndian.Uint32(dec[0:4])
	if magic != sfaChunkFMSH {
		return fmt.Errorf("FMSH bad magic 0x%X", magic)
	}
	// version dec[4:8]
	nCells := binary.LittleEndian.Uint32(dec[8:12])
	nFaces := binary.LittleEndian.Uint32(dec[12:16])
	flags := dec[16]
	// reserved dec[17:20]
	L := math.Float64frombits(binary.LittleEndian.Uint64(dec[20:28]))
	posF64 := (flags & 0x04) != 0
	posSize := 4
	if posF64 {
		posSize = 8
	}
	off := 28
	posBytes := int(nCells) * 3 * posSize
	volBytes := int(nCells) * posSize
	if len(dec) < off+posBytes+volBytes {
		return fmt.Errorf("FMSH payload truncated (cells=%d)", nCells)
	}
	pos := make([]float64, 3*int(nCells))
	vol := make([]float64, int(nCells))
	if posF64 {
		for i := 0; i < 3*int(nCells); i++ {
			pos[i] = math.Float64frombits(binary.LittleEndian.Uint64(dec[off : off+8]))
			off += 8
		}
		for i := 0; i < int(nCells); i++ {
			vol[i] = math.Float64frombits(binary.LittleEndian.Uint64(dec[off : off+8]))
			off += 8
		}
	} else {
		for i := 0; i < 3*int(nCells); i++ {
			pos[i] = float64(math.Float32frombits(binary.LittleEndian.Uint32(dec[off : off+4])))
			off += 4
		}
		for i := 0; i < int(nCells); i++ {
			vol[i] = float64(math.Float32frombits(binary.LittleEndian.Uint32(dec[off : off+4])))
			off += 4
		}
	}
	_ = nFaces // we ignore faces/CSR for visualization

	s.meshNCells = nCells
	s.meshL = L
	s.meshCellPos = pos
	s.meshCellVol = vol
	// Reset state buffer; size computed from the file's actual column count.
	N := int(nCells) * int(s.nColumns)
	if cap(s.cellState) < N {
		s.cellState = make([]float32, N)
	} else {
		s.cellState = s.cellState[:N]
		for i := range s.cellState {
			s.cellState[i] = 0
		}
	}
	s.cellTValid = false
	s.cellTMean = nil
	s.cellTAmp = nil
	s.cellTPhase = nil
	return nil
}

// resampleCellStateToRawBuf takes the current per-cell state buffer
// (column-major, cellState[col*nCells + cell]) and writes the voxel-grid
// equivalent into rawBuf via the precomputed voxelToCell map.
func (s *sfaFile) resampleCellStateToRawBuf() {
	N3 := int(s.nTotal)
	nCells := int(s.meshNCells)
	for col := 0; col < int(s.nColumns); col++ {
		colSrc := s.cellState[col*nCells : (col+1)*nCells]
		dstOff := col * N3 * 4
		for v := 0; v < N3; v++ {
			c := s.voxelToCell[v]
			var val float32
			if c >= 0 {
				val = colSrc[c]
			}
			binary.LittleEndian.PutUint32(s.rawBuf[dstOff+v*4:dstOff+v*4+4], math.Float32bits(val))
		}
	}
}

// buildVoxelToCell precomputes the nearest-cell map for the rendering grid.
// Uses a uniform spatial bin index, identical pattern to foam_to_voxel.
func (s *sfaFile) buildVoxelToCell() {
	N := int(s.cellVoxelN)
	N3 := int64(N) * int64(N) * int64(N)
	if int64(len(s.voxelToCell)) != N3 {
		s.voxelToCell = make([]int32, N3)
	}
	nCells := int(s.meshNCells)
	L := s.meshL
	boxVol := math.Pow(2.0*L, 3.0)
	cellVolAvg := boxVol / float64(nCells)
	binSize := math.Pow(cellVolAvg, 1.0/3.0) * 1.5
	Ng := int(math.Ceil(2.0 * L / binSize))
	binSize = 2.0 * L / float64(Ng)
	Ng3 := Ng * Ng * Ng

	binCount := make([]int32, Ng3)
	binOff := make([]int32, Ng3)
	cellsByBin := make([]int32, nCells)
	for i := 0; i < nCells; i++ {
		bx := int((s.meshCellPos[3*i+0] + L) / binSize)
		by := int((s.meshCellPos[3*i+1] + L) / binSize)
		bz := int((s.meshCellPos[3*i+2] + L) / binSize)
		if bx < 0 { bx = 0 }; if bx >= Ng { bx = Ng - 1 }
		if by < 0 { by = 0 }; if by >= Ng { by = Ng - 1 }
		if bz < 0 { bz = 0 }; if bz >= Ng { bz = Ng - 1 }
		binCount[bx*Ng*Ng+by*Ng+bz]++
	}
	total := int32(0)
	for b := 0; b < Ng3; b++ {
		binOff[b] = total
		total += binCount[b]
	}
	cursor := make([]int32, Ng3)
	copy(cursor, binOff)
	for i := 0; i < nCells; i++ {
		bx := int((s.meshCellPos[3*i+0] + L) / binSize)
		by := int((s.meshCellPos[3*i+1] + L) / binSize)
		bz := int((s.meshCellPos[3*i+2] + L) / binSize)
		if bx < 0 { bx = 0 }; if bx >= Ng { bx = Ng - 1 }
		if by < 0 { by = 0 }; if by >= Ng { by = Ng - 1 }
		if bz < 0 { bz = 0 }; if bz >= Ng { bz = Ng - 1 }
		cellsByBin[cursor[bx*Ng*Ng+by*Ng+bz]] = int32(i)
		cursor[bx*Ng*Ng+by*Ng+bz]++
	}

	dx := 2.0 * L / float64(N-1)
	S := 2.0 * L
	wrap := func(v int) int {
		v = v % Ng
		if v < 0 {
			v += Ng
		}
		return v
	}

	// Parallelise over i (outer voxel index) with goroutines
	work := make(chan int, 8)
	done := make(chan struct{})
	worker := func() {
		for i := range work {
			x := -L + float64(i)*dx
			bx0 := wrap(int((x + L) / binSize))
			for j := 0; j < N; j++ {
				y := -L + float64(j)*dx
				by0 := wrap(int((y + L) / binSize))
				for k := 0; k < N; k++ {
					z := -L + float64(k)*dx
					bz0 := wrap(int((z + L) / binSize))
					best := int32(-1)
					bestD2 := 1e30
					for di := -1; di <= 1; di++ {
						gi := wrap(bx0 + di)
						for dj := -1; dj <= 1; dj++ {
							gj := wrap(by0 + dj)
							for dk := -1; dk <= 1; dk++ {
								gk := wrap(bz0 + dk)
								b := gi*Ng*Ng + gj*Ng + gk
								n := int(binCount[b])
								o := int(binOff[b])
								for q := 0; q < n; q++ {
									c := int(cellsByBin[o+q])
									cdx := s.meshCellPos[3*c+0] - x
									cdy := s.meshCellPos[3*c+1] - y
									cdz := s.meshCellPos[3*c+2] - z
									if cdx > 0.5*S { cdx -= S }
									if cdx < -0.5*S { cdx += S }
									if cdy > 0.5*S { cdy -= S }
									if cdy < -0.5*S { cdy += S }
									if cdz > 0.5*S { cdz -= S }
									if cdz < -0.5*S { cdz += S }
									d2 := cdx*cdx + cdy*cdy + cdz*cdz
									if d2 < bestD2 {
										bestD2 = d2
										best = int32(c)
									}
								}
							}
						}
					}
					s.voxelToCell[int64(i)*int64(N)*int64(N)+int64(j)*int64(N)+int64(k)] = best
				}
			}
		}
		done <- struct{}{}
	}
	const nWorkers = 8
	for w := 0; w < nWorkers; w++ {
		go worker()
	}
	for i := 0; i < N; i++ {
		work <- i
	}
	close(work)
	for w := 0; w < nWorkers; w++ {
		<-done
	}
}

// readCellPayload decompresses an FCEL chunk, copies cell values into the
// per-cell state buffer, and (for v2 with bit2 flag) loads the temporal
// model. Then resamples state to rawBuf for the rendering pipeline.
func (s *sfaFile) readCellPayload(entry sfaL2Entry) error {
	dataOff := entry.offset + 12 + 8
	var compressed []byte
	if s.mmapData != nil && dataOff+entry.compressedSize <= uint64(len(s.mmapData)) {
		compressed = s.mmapData[dataOff : dataOff+entry.compressedSize]
	} else {
		if uint64(cap(s.compBuf)) < entry.compressedSize {
			s.compBuf = make([]byte, entry.compressedSize)
		}
		compressed = s.compBuf[:entry.compressedSize]
		s.fp.Seek(int64(dataOff), io.SeekStart)
		if _, err := io.ReadFull(s.fp, compressed); err != nil {
			return fmt.Errorf("read FCEL: %w", err)
		}
	}
	dec, err := zstdDecompress(compressed, uint64(len(compressed)*4))
	if err != nil {
		return fmt.Errorf("decompress FCEL: %w", err)
	}
	if len(dec) < 20 {
		return fmt.Errorf("FCEL payload too short")
	}
	magic := binary.LittleEndian.Uint32(dec[0:4])
	if magic != sfaChunkFCEL {
		return fmt.Errorf("FCEL bad magic")
	}
	version := binary.LittleEndian.Uint32(dec[4:8])
	nCells := binary.LittleEndian.Uint32(dec[8:12])
	nCols := binary.LittleEndian.Uint32(dec[12:16])
	dtype := dec[16]
	flags := dec[17]
	if nCells != s.meshNCells {
		return fmt.Errorf("FCEL N_cells=%d != mesh N_cells=%d", nCells, s.meshNCells)
	}
	es := sfaDtypeSize[dtype]
	off := 20
	colBytes := int(nCells) * es
	if len(dec) < off+int(nCols)*colBytes {
		return fmt.Errorf("FCEL payload truncated")
	}

	// Copy each column into cellState (column-major)
	stateNeeded := int(s.nColumns) * int(nCells)
	if cap(s.cellState) < stateNeeded {
		s.cellState = make([]float32, stateNeeded)
	} else {
		s.cellState = s.cellState[:stateNeeded]
	}
	for col := uint32(0); col < nCols && col < s.nColumns; col++ {
		colStart := off + int(col)*colBytes
		dst := s.cellState[int(col)*int(nCells) : int(col+1)*int(nCells)]
		for c := uint32(0); c < nCells; c++ {
			switch dtype {
			case sfaDtypeF32:
				dst[c] = math.Float32frombits(binary.LittleEndian.Uint32(dec[colStart+int(c)*4 : colStart+int(c)*4+4]))
			case sfaDtypeF64:
				dst[c] = float32(math.Float64frombits(binary.LittleEndian.Uint64(dec[colStart+int(c)*8 : colStart+int(c)*8+8])))
			case sfaDtypeF16:
				dst[c] = f16ToF32(binary.LittleEndian.Uint16(dec[colStart+int(c)*2 : colStart+int(c)*2+2]))
			}
		}
	}

	// Optional v2: temporal model
	if version >= 2 && (flags&0x04) != 0 {
		modelOff := off + int(nCols)*colBytes
		if modelOff+4 > len(dec) {
			return fmt.Errorf("FCEL v2 missing model header")
		}
		s.cellTOmega = math.Float32frombits(binary.LittleEndian.Uint32(dec[modelOff : modelOff+4]))
		modelOff += 4
		mb := int(nCols) * int(nCells) * 4
		if modelOff+3*mb > len(dec) {
			return fmt.Errorf("FCEL v2 model truncated")
		}
		// Read mean / amp / phase
		s.cellTMean = make([]float32, int(nCols)*int(nCells))
		s.cellTAmp = make([]float32, int(nCols)*int(nCells))
		s.cellTPhase = make([]float32, int(nCols)*int(nCells))
		for i := 0; i < int(nCols)*int(nCells); i++ {
			s.cellTMean[i] = math.Float32frombits(binary.LittleEndian.Uint32(dec[modelOff+i*4 : modelOff+i*4+4]))
			s.cellTAmp[i] = math.Float32frombits(binary.LittleEndian.Uint32(dec[modelOff+mb+i*4 : modelOff+mb+i*4+4]))
			s.cellTPhase[i] = math.Float32frombits(binary.LittleEndian.Uint32(dec[modelOff+2*mb+i*4 : modelOff+2*mb+i*4+4]))
		}
		s.cellTValid = true
	}

	if len(s.rawBuf) < int(s.nColumns)*int(s.nTotal)*4 {
		s.rawBuf = make([]byte, int(s.nColumns)*int(s.nTotal)*4)
	}
	s.resampleCellStateToRawBuf()
	return nil
}

// readTemporalModelPayload decompresses a FMTL chunk and stores the
// per-cell Fourier model parameters. Subsequent FCEP frames use this
// model for prediction.
func (s *sfaFile) readTemporalModelPayload(entry sfaL2Entry) error {
	dataOff := entry.offset + 12 + 8
	var compressed []byte
	if s.mmapData != nil && dataOff+entry.compressedSize <= uint64(len(s.mmapData)) {
		compressed = s.mmapData[dataOff : dataOff+entry.compressedSize]
	} else {
		if uint64(cap(s.compBuf)) < entry.compressedSize {
			s.compBuf = make([]byte, entry.compressedSize)
		}
		compressed = s.compBuf[:entry.compressedSize]
		s.fp.Seek(int64(dataOff), io.SeekStart)
		if _, err := io.ReadFull(s.fp, compressed); err != nil {
			return fmt.Errorf("read FMTL: %w", err)
		}
	}
	dec, err := zstdDecompress(compressed, uint64(len(compressed)*4))
	if err != nil {
		return fmt.Errorf("decompress FMTL: %w", err)
	}
	if len(dec) < 28 {
		return fmt.Errorf("FMTL payload too short")
	}
	magic := binary.LittleEndian.Uint32(dec[0:4])
	if magic != sfaChunkFMTL {
		return fmt.Errorf("FMTL bad magic")
	}
	nCells := binary.LittleEndian.Uint32(dec[8:12])
	nCols := binary.LittleEndian.Uint32(dec[12:16])
	// flags dec[16], reserved dec[17:20]
	s.cellTOmega = math.Float32frombits(binary.LittleEndian.Uint32(dec[20:24]))
	// pad dec[24:28]
	mb := int(nCols) * int(nCells) * 4
	if 28+3*mb > len(dec) {
		return fmt.Errorf("FMTL truncated")
	}
	s.cellTMean = make([]float32, int(nCols)*int(nCells))
	s.cellTAmp = make([]float32, int(nCols)*int(nCells))
	s.cellTPhase = make([]float32, int(nCols)*int(nCells))
	for i := 0; i < int(nCols)*int(nCells); i++ {
		s.cellTMean[i] = math.Float32frombits(binary.LittleEndian.Uint32(dec[28+i*4 : 28+i*4+4]))
		s.cellTAmp[i] = math.Float32frombits(binary.LittleEndian.Uint32(dec[28+mb+i*4 : 28+mb+i*4+4]))
		s.cellTPhase[i] = math.Float32frombits(binary.LittleEndian.Uint32(dec[28+2*mb+i*4 : 28+2*mb+i*4+4]))
	}
	s.cellTValid = true
	return nil
}

// readCellPPayload decompresses an FCEP chunk: applies sparse delta on
// top of model prediction at this frame's time t. Updates cellState
// and rawBuf.
func (s *sfaFile) readCellPPayload(entry sfaL2Entry) error {
	if !s.cellTValid {
		return fmt.Errorf("FCEP frame without valid temporal model")
	}
	dataOff := entry.offset + 12 + 8
	var compressed []byte
	if s.mmapData != nil && dataOff+entry.compressedSize <= uint64(len(s.mmapData)) {
		compressed = s.mmapData[dataOff : dataOff+entry.compressedSize]
	} else {
		if uint64(cap(s.compBuf)) < entry.compressedSize {
			s.compBuf = make([]byte, entry.compressedSize)
		}
		compressed = s.compBuf[:entry.compressedSize]
		s.fp.Seek(int64(dataOff), io.SeekStart)
		if _, err := io.ReadFull(s.fp, compressed); err != nil {
			return fmt.Errorf("read FCEP: %w", err)
		}
	}
	dec, err := zstdDecompress(compressed, uint64(len(compressed)*4))
	if err != nil {
		return fmt.Errorf("decompress FCEP: %w", err)
	}
	if len(dec) < 24 {
		return fmt.Errorf("FCEP payload too short")
	}
	magic := binary.LittleEndian.Uint32(dec[0:4])
	if magic != sfaChunkFCEP {
		return fmt.Errorf("FCEP bad magic")
	}
	nCells := binary.LittleEndian.Uint32(dec[8:12])
	nCols := binary.LittleEndian.Uint32(dec[12:16])
	nChanged := binary.LittleEndian.Uint32(dec[16:20])
	if nCells != s.meshNCells {
		return fmt.Errorf("FCEP N_cells mismatch")
	}
	off := 24

	// Reconstruct state from temporal model.
	// Model storage (writer side) is CELL-major: t_mean[c*6 + col].
	// State storage is COLUMN-major: cellState[col*N_cells + c].
	t := float32(entry.time)
	for c := 0; c < int(nCells); c++ {
		mbase := c * int(nCols)        // model index base
		for col := 0; col < int(nCols); col++ {
			pred := s.cellTMean[mbase+col] +
				s.cellTAmp[mbase+col]*float32(math.Cos(float64(s.cellTOmega*t+s.cellTPhase[mbase+col])))
			s.cellState[col*int(nCells)+c] = pred
		}
	}
	// Apply sparse deltas. Delta layout in payload: per cell record,
	// n_columns floats column-major within that record.
	idsOff := off
	deltasOff := off + int(nChanged)*4
	for i := uint32(0); i < nChanged; i++ {
		c := binary.LittleEndian.Uint32(dec[idsOff+int(i)*4 : idsOff+int(i)*4+4])
		for col := uint32(0); col < nCols; col++ {
			d := math.Float32frombits(binary.LittleEndian.Uint32(
				dec[deltasOff+int(i)*int(nCols)*4+int(col)*4 : deltasOff+int(i)*int(nCols)*4+int(col)*4+4]))
			s.cellState[int(col)*int(nCells)+int(c)] += d
		}
	}

	if len(s.rawBuf) < int(s.nColumns)*int(s.nTotal)*4 {
		s.rawBuf = make([]byte, int(s.nColumns)*int(s.nTotal)*4)
	}
	s.resampleCellStateToRawBuf()
	return nil
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
			binary.Read(s.fp, binary.LittleEndian, &entry.frameType)
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

	// Vector frame: decompress FRVD payload, evaluate patches to voxels
	if entry.frameType == sfaFrameVecI || entry.frameType == sfaFrameVecP {
		return s.readVecFrame(frameIdx, *entry)
	}

	// Mesh frame (FMSH): parse cell positions, build voxel-to-cell map.
	// No voxel data is produced — caller should advance to the next frame
	// (a subsequent FCEL).
	if entry.frameType == sfaFrameMesh {
		if err := s.readMeshPayload(*entry); err != nil {
			return 0, err
		}
		s.buildVoxelToCell()
		s.cellLoadedMesh = int32(frameIdx)
		// Zero rawBuf so any consumer that reads voxels gets a clean slate.
		for i := range s.rawBuf {
			s.rawBuf[i] = 0
		}
		return entry.time, nil
	}

	// Standalone temporal-model frame (FMTL): just load the model.
	if entry.frameType == sfaFrameTemporalModel {
		if err := s.readTemporalModelPayload(*entry); err != nil {
			return 0, err
		}
		return entry.time, nil
	}

	// Cell P-frame (FCEP): needs the most recent FMSH AND a valid
	// temporal model (from either an FCEL v2 I-frame or an FMTL frame).
	if entry.frameType == sfaFrameCellP {
		if !s.cellTValid || s.cellLoadedMesh < 0 {
			// Locate prior FMSH
			meshIdx := int32(-1)
			for back := int32(frameIdx) - 1; back >= 0; back-- {
				e2, err := s.findFrame(uint32(back))
				if err != nil {
					continue
				}
				if e2.frameType == sfaFrameMesh && meshIdx < 0 {
					if err := s.readMeshPayload(*e2); err == nil {
						s.buildVoxelToCell()
						s.cellLoadedMesh = back
						meshIdx = back
						break
					}
				}
			}
			if meshIdx < 0 {
				return 0, fmt.Errorf("FCEP %d: no FMSH found", frameIdx)
			}
			// Locate most recent model source (FMTL or FCEL v2)
			// at or before frameIdx, plus the most recent FCEL (data).
			modelLoaded := false
			cellLoaded := false
			for back := int32(frameIdx) - 1; back >= meshIdx; back-- {
				e2, err := s.findFrame(uint32(back))
				if err != nil {
					continue
				}
				if !modelLoaded && e2.frameType == sfaFrameTemporalModel {
					if err := s.readTemporalModelPayload(*e2); err == nil {
						modelLoaded = true
					}
				}
				if !cellLoaded && e2.frameType == sfaFrameCell {
					if err := s.readCellPayload(*e2); err == nil {
						cellLoaded = true
						if s.cellTValid {
							modelLoaded = true
						}
					}
				}
				if modelLoaded && cellLoaded {
					break
				}
			}
			if !modelLoaded {
				return 0, fmt.Errorf("FCEP %d: no temporal model found (FMTL or FCEL v2)", frameIdx)
			}
		}
		if err := s.readCellPPayload(*entry); err != nil {
			return 0, err
		}
		return entry.time, nil
	}

	// Cell-data frame (FCEL): resample current mesh's cell values to the
	// voxel rendering grid.
	if entry.frameType == sfaFrameCell {
		if s.cellLoadedMesh < 0 || s.voxelToCell == nil {
			// Find the most-recent FMSH frame at or before frameIdx and load it.
			loaded := false
			for back := int32(frameIdx); back >= 0; back-- {
				e2, err := s.findFrame(uint32(back))
				if err != nil {
					continue
				}
				if e2.frameType == sfaFrameMesh {
					if err := s.readMeshPayload(*e2); err == nil {
						s.buildVoxelToCell()
						s.cellLoadedMesh = back
						loaded = true
					}
					break
				}
			}
			if !loaded {
				return 0, fmt.Errorf("FCEL frame %d with no preceding mesh", frameIdx)
			}
		}
		if err := s.readCellPayload(*entry); err != nil {
			return 0, err
		}
		return entry.time, nil
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

// readVecFrame reads an FRVD vector frame (I or P), maintains patch state,
// evaluates patches to voxels, and fills rawBuf for the rendering pipeline.
func (s *sfaFile) readVecFrame(frameIdx uint32, entry sfaL2Entry) (float64, error) {
	// For P-frames: ensure the correct temporal model is loaded.
	// Find the preceding I-frame and load its model if needed.
	if entry.frameType == sfaFrameVecP {
		// Search backward for the I-frame that covers this P-frame
		needIFrame := int32(-1)
		for search := int(frameIdx) - 1; search >= 0; search-- {
			e, err := s.findFrame(uint32(search))
			if err != nil {
				continue
			}
			if e.frameType == sfaFrameVecI {
				needIFrame = int32(search)
				break
			}
			if e.frameType == sfaFrameVoxel {
				break // crossed a voxel frame boundary — no I-frame before us
			}
		}
		if needIFrame >= 0 && needIFrame != s.vecTempIFrameIdx {
			// Load the correct I-frame to get its temporal model
			iEntry, err := s.findFrame(uint32(needIFrame))
			if err == nil && iEntry.frameType == sfaFrameVecI {
				s.readVecFrame(uint32(needIFrame), *iEntry)
				// Now the temporal model is loaded from the correct I-frame
			}
		}
	}
	// Read and decompress payload
	// FRVD layout: chunk_header(12) + time(8) + compressed_data
	dataOff := entry.offset + 12 + 8 // skip chunk header + embedded time
	var compressed []byte
	if s.mmapData != nil && dataOff+entry.compressedSize <= uint64(len(s.mmapData)) {
		compressed = s.mmapData[dataOff : dataOff+entry.compressedSize]
	} else {
		compressed = make([]byte, entry.compressedSize)
		s.fp.Seek(int64(dataOff), io.SeekStart)
		if _, err := io.ReadFull(s.fp, compressed); err != nil {
			return 0, fmt.Errorf("read FRVD: %w", err)
		}
	}
	payload, err := zstdDecoder.DecodeAll(compressed, nil)
	if err != nil {
		return 0, fmt.Errorf("decompress FRVD: %w", err)
	}

	if entry.frameType == sfaFrameVecI {
		// I-frame: [n_patches(4)][bs(1)][nc(2)][flags(1)][origins(n*6)][coeffs(n*nc*4)]
		// If flags & 1: [omega(4)][mean(n*nc*4)][amp(n*nc*4)][phase(n*nc*4)]
		if len(payload) < 8 {
			return 0, fmt.Errorf("FRVD I-frame too short")
		}
		np := binary.LittleEndian.Uint32(payload[0:4])
		bs := int(payload[4])
		nc := int(binary.LittleEndian.Uint16(payload[5:7]))
		flags := payload[7]

		originOff := 8
		coeffOff := originOff + int(np)*6
		if len(payload) < coeffOff+int(np)*nc*4 {
			return 0, fmt.Errorf("FRVD I-frame truncated")
		}

		s.vecNPatches = np
		s.vecNCoeffs = uint16(nc)
		s.vecBS = bs
		s.vecOrder = 3
		if nc <= 27 { s.vecOrder = 2 }
		if nc <= 8 { s.vecOrder = 1 }

		s.vecOrigins = make([][3]int16, np)
		for p := uint32(0); p < np; p++ {
			off := originOff + int(p)*6
			s.vecOrigins[p][0] = int16(binary.LittleEndian.Uint16(payload[off:]))
			s.vecOrigins[p][1] = int16(binary.LittleEndian.Uint16(payload[off+2:]))
			s.vecOrigins[p][2] = int16(binary.LittleEndian.Uint16(payload[off+4:]))
		}

		nTotal := int(np) * nc
		s.vecCoeffs = make([]float32, nTotal)
		for i := range s.vecCoeffs {
			s.vecCoeffs[i] = math.Float32frombits(binary.LittleEndian.Uint32(payload[coeffOff+i*4:]))
		}

		// Parse temporal model if present
		modelOff := coeffOff + nTotal*4
		if flags&0x01 != 0 && len(payload) >= modelOff+4+nTotal*4*3 {
			s.vecTempOmega = math.Float32frombits(binary.LittleEndian.Uint32(payload[modelOff:]))
			modelOff += 4
			s.vecTempMean = make([]float32, nTotal)
			s.vecTempAmp = make([]float32, nTotal)
			s.vecTempPhase = make([]float32, nTotal)
			for i := 0; i < nTotal; i++ {
				s.vecTempMean[i] = math.Float32frombits(binary.LittleEndian.Uint32(payload[modelOff+i*4:]))
			}
			modelOff += nTotal * 4
			for i := 0; i < nTotal; i++ {
				s.vecTempAmp[i] = math.Float32frombits(binary.LittleEndian.Uint32(payload[modelOff+i*4:]))
			}
			modelOff += nTotal * 4
			for i := 0; i < nTotal; i++ {
				s.vecTempPhase[i] = math.Float32frombits(binary.LittleEndian.Uint32(payload[modelOff+i*4:]))
			}
			s.vecTempValid = true
			s.vecTempIFrameIdx = int32(frameIdx)
		} else {
			s.vecTempValid = false
			s.vecTempIFrameIdx = int32(frameIdx)
		}

	} else if entry.frameType == sfaFrameVecP {
		// P-frame: [n_deltas(4)][n_coeffs(2)][pad(2)][indices(n*4)][delta_coeffs(n*nc*4)]
		if len(payload) < 8 {
			return 0, fmt.Errorf("FRVD P-frame too short")
		}
		nd := binary.LittleEndian.Uint32(payload[0:4])
		nc := int(binary.LittleEndian.Uint16(payload[4:6]))

		if s.vecCoeffs == nil {
			return entry.time, nil
		}

		nTotal := len(s.vecCoeffs)

		// Start from temporal prediction or previous coefficients
		if s.vecTempValid && len(s.vecTempMean) == nTotal {
			t := float32(entry.time)
			for i := 0; i < nTotal; i++ {
				s.vecCoeffs[i] = s.vecTempMean[i] +
					s.vecTempAmp[i]*float32(math.Cos(float64(s.vecTempOmega*t+s.vecTempPhase[i])))
			}
		}
		// If no temporal model, vecCoeffs already holds previous frame state

		idxOff := 8
		coeffOff := idxOff + int(nd)*4

		for d := uint32(0); d < nd; d++ {
			patchIdx := binary.LittleEndian.Uint32(payload[idxOff+int(d)*4:])
			if int(patchIdx) >= nTotal/nc {
				continue
			}
			base := int(patchIdx) * nc
			dOff := coeffOff + int(d)*nc*4
			for c := 0; c < nc; c++ {
				delta := math.Float32frombits(binary.LittleEndian.Uint32(payload[dOff+c*4:]))
				s.vecCoeffs[base+c] += delta
			}
		}
	}

	// Evaluate current patch state to voxels in rawBuf
	s.vecPatchesToRawBuf()
	return entry.time, nil
}

// vecPatchesToRawBuf evaluates all stored patches into rawBuf for rendering.
func (s *sfaFile) vecPatchesToRawBuf() {
	if s.vecCoeffs == nil || s.vecNPatches == 0 {
		return
	}

	// Clear rawBuf
	for i := range s.rawBuf {
		s.rawBuf[i] = 0
	}

	nx := int(s.Nx)
	ny := int(s.Ny)
	nz := int(s.Nz)
	nc := int(s.vecNCoeffs)
	bs := s.vecBS
	o1 := s.vecOrder + 1

	// Detect multi-field layout
	pfNC := nc    // per-field coefficients
	nFields := 1
	if nc > 64 {
		if nc%64 == 0 { nFields = nc / 64; pfNC = 64 } else if nc%27 == 0 { nFields = nc / 27; pfNC = 27 } else if nc%8 == 0 { nFields = nc / 8; pfNC = 8 }
	}
	pfOrder := 3
	if pfNC <= 27 { pfOrder = 2 }
	if pfNC <= 8 { pfOrder = 1 }
	pfO1 := pfOrder + 1
	_ = o1 // use pfO1 instead

	// Precompute column byte offsets
	colOffsets := make([]uint64, s.nColumns)
	{ off := uint64(0)
	  for cc := uint32(0); cc < s.nColumns; cc++ {
		  colOffsets[cc] = off
		  off += s.nTotal * uint64(sfaDtypeSize[s.columns[cc].dtype])
	  }
	}

	for pi := uint32(0); pi < s.vecNPatches; pi++ {
		ox := int(s.vecOrigins[pi][0])
		oy := int(s.vecOrigins[pi][1])
		oz := int(s.vecOrigins[pi][2])
		cBase := int(pi) * nc

		for di := 0; di < bs; di++ {
			gi := ox + di
			if gi < 0 || gi >= nx { continue }
			var tx float64
			if bs > 1 { tx = float64(di) / float64(bs-1) }
			sx := 2*tx - 1
			txa := [4]float64{1, sx, 2*sx*sx - 1, 4*sx*sx*sx - 3*sx}

			for dj := 0; dj < bs; dj++ {
				gj := oy + dj
				if gj < 0 || gj >= ny { continue }
				var ty float64
				if bs > 1 { ty = float64(dj) / float64(bs-1) }
				sy := 2*ty - 1
				tya := [4]float64{1, sy, 2*sy*sy - 1, 4*sy*sy*sy - 3*sy}

				for dk := 0; dk < bs; dk++ {
					gk := oz + dk
					if gk < 0 || gk >= nz { continue }
					var tz float64
					if bs > 1 { tz = float64(dk) / float64(bs-1) }
					sz := 2*tz - 1
					tza := [4]float64{1, sz, 2*sz*sz - 1, 4*sz*sz*sz - 3*sz}

					idx := int64(gi)*int64(ny)*int64(nz) + int64(gj)*int64(nz) + int64(gk)

					// Evaluate each field separately
					for field := 0; field < nFields && field < int(s.nColumns); field++ {
						fcBase := cBase + field*pfNC
						var val float64
						for a := 0; a < pfO1; a++ {
							for b := 0; b < pfO1; b++ {
								ab := txa[a] * tya[b]
								for c := 0; c < pfO1; c++ {
									ci := a*pfO1*pfO1 + b*pfO1 + c
									val += float64(s.vecCoeffs[fcBase+ci]) * ab * tza[c]
								}
							}
						}

						col := uint32(field)
						es := sfaDtypeSize[s.columns[col].dtype]
						f := float32(val)
						if s.columns[col].dtype == sfaDtypeF16 {
							bits := math.Float32bits(f)
							sign := uint16((bits >> 16) & 0x8000)
							exp := int((bits>>23)&0xFF) - 127 + 15
							mant := uint16((bits >> 13) & 0x3FF)
							var h uint16
							if exp <= 0 { h = sign } else if exp >= 31 { h = sign | 0x7C00 } else { h = sign | uint16(exp)<<10 | mant }
							off := colOffsets[col] + uint64(idx)*uint64(es)
							if off+1 < uint64(len(s.rawBuf)) {
								s.rawBuf[off] = byte(h)
								s.rawBuf[off+1] = byte(h >> 8)
							}
						} else if s.columns[col].dtype == sfaDtypeF32 {
							bits := math.Float32bits(f)
							off := colOffsets[col] + uint64(idx)*4
							if off+3 < uint64(len(s.rawBuf)) {
								s.rawBuf[off] = byte(bits)
								s.rawBuf[off+1] = byte(bits >> 8)
								s.rawBuf[off+2] = byte(bits >> 16)
								s.rawBuf[off+3] = byte(bits >> 24)
							}
						}
					}
				}
			}
		}
	}
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
	var colsBuf [32]colInfo
	var cols []colInfo
	if nc <= len(colsBuf) {
		cols = colsBuf[:nc]
	} else {
		cols = make([]colInfo, nc)
	}

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
	errs := make([]error, nc)

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

			if es == 1 && s.columns[c].dtype == 8 { // SFA_U8
				// FAST PATH: uint8 → float32 (0-1 range)
				out := s.colF32Bufs[c]
				for v := 0; v < nVals; v++ {
					out[v] = float32(bssBuf[v]) / 255.0
				}
			} else if es == 4 && s.columns[c].dtype == sfaDtypeF32 {
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
	// If COLZSTD with f32 or u8, the fused decode already wrote to colF32Bufs
	codec := s.flags & 0xF
	if codec == sfaCodecColZstd && (s.columns[colIdx].dtype == sfaDtypeF32 || s.columns[colIdx].dtype == 8) {
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
	case 8: // SFA_U8
		for i := uint64(0); i < s.nTotal; i++ {
			result[i] = float32(s.rawBuf[off+i]) / 255.0
		}
	case 9: // SFA_U16
		for i := uint64(0); i < s.nTotal; i++ {
			v := binary.LittleEndian.Uint16(s.rawBuf[off+i*2 : off+i*2+2])
			result[i] = float32(v) / 65535.0
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
// Vec3D file reader
// ============================================================

type vec3dPatch struct {
	fitType   uint8
	fieldTerm uint8
	nCoeffs   uint16
	bx, by, bz int16
	sx, sy, sz int16
	order     uint8
	maxError  float32
	rmsError  float32
	coeffs    []float32
}

type vec3dField struct {
	fieldTerm uint8
	patches   []vec3dPatch
}

type vec3dFile struct {
	version        uint32
	nFields        uint32
	Nx, Ny, Nz     uint32
	fields         []vec3dField
}

// Field term names: 0=phi_x, 1=phi_y, 2=phi_z, 3=theta_x, 4=theta_y, 5=theta_z
var fieldTermColors = [][3]float32{
	{1.0, 0.0, 0.0}, // phi_x = red
	{0.0, 1.0, 0.0}, // phi_y = green
	{0.0, 0.0, 1.0}, // phi_z = blue
	{0.0, 1.0, 1.0}, // theta_x = cyan
	{1.0, 0.0, 1.0}, // theta_y = magenta
	{1.0, 1.0, 0.0}, // theta_z = yellow
}

func vec3dOpen(path string) (*vec3dFile, error) {
	fp, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer fp.Close()

	// Read magic
	var magic [8]byte
	if _, err := io.ReadFull(fp, magic[:]); err != nil {
		return nil, fmt.Errorf("read magic: %w", err)
	}
	if string(magic[:5]) != "VEC3D" {
		return nil, fmt.Errorf("not a VEC3D file (got %q)", magic)
	}

	v := &vec3dFile{}
	binary.Read(fp, binary.LittleEndian, &v.version)
	binary.Read(fp, binary.LittleEndian, &v.nFields)
	binary.Read(fp, binary.LittleEndian, &v.Nx)
	binary.Read(fp, binary.LittleEndian, &v.Ny)
	binary.Read(fp, binary.LittleEndian, &v.Nz)

	v.fields = make([]vec3dField, v.nFields)
	totalPatches := 0
	for f := uint32(0); f < v.nFields; f++ {
		// VecFrameHeader: magic(4) + n_patches(4) + field_term(1) + reserved(3) = 12 bytes
		var frameMagic uint32
		var nPatches uint32
		var fieldTerm uint8
		var reserved [3]byte
		binary.Read(fp, binary.LittleEndian, &frameMagic)
		binary.Read(fp, binary.LittleEndian, &nPatches)
		binary.Read(fp, binary.LittleEndian, &fieldTerm)
		fp.Read(reserved[:])

		if frameMagic != 0x56454333 { // "VEC3"
			return nil, fmt.Errorf("field %d: bad VEC3 magic 0x%08x", f, frameMagic)
		}

		v.fields[f].fieldTerm = fieldTerm
		v.fields[f].patches = make([]vec3dPatch, nPatches)
		totalPatches += int(nPatches)

		for i := uint32(0); i < nPatches; i++ {
			p := &v.fields[f].patches[i]
			// VecPatchHeader: 28 bytes total (C struct with float alignment padding)
			// Layout: fit_type(1) field_term(1) n_coeffs(2) bx(2) by(2) bz(2)
			//         sx(2) sy(2) sz(2) order(1) reserved(1) pad(2) max_error(4) rms_error(4)
			binary.Read(fp, binary.LittleEndian, &p.fitType)
			binary.Read(fp, binary.LittleEndian, &p.fieldTerm)
			binary.Read(fp, binary.LittleEndian, &p.nCoeffs)
			binary.Read(fp, binary.LittleEndian, &p.bx)
			binary.Read(fp, binary.LittleEndian, &p.by)
			binary.Read(fp, binary.LittleEndian, &p.bz)
			binary.Read(fp, binary.LittleEndian, &p.sx)
			binary.Read(fp, binary.LittleEndian, &p.sy)
			binary.Read(fp, binary.LittleEndian, &p.sz)
			binary.Read(fp, binary.LittleEndian, &p.order)
			var pad [3]byte // 1 byte reserved + 2 bytes alignment padding before float
			fp.Read(pad[:])
			binary.Read(fp, binary.LittleEndian, &p.maxError)
			binary.Read(fp, binary.LittleEndian, &p.rmsError)

			p.coeffs = make([]float32, p.nCoeffs)
			for c := uint16(0); c < p.nCoeffs; c++ {
				binary.Read(fp, binary.LittleEndian, &p.coeffs[c])
			}
		}
	}

	fmt.Printf("VEC3D: %dx%dx%d, %d fields, %d total patches\n",
		v.Nx, v.Ny, v.Nz, v.nFields, totalPatches)
	return v, nil
}

// buildWireframeData generates line vertices for all patches as wireframe boxes.
// Returns interleaved [x,y,z, r,g,b,a] float32 arrays and vertex count.
// Positions are normalized to [0,1]^3.
func (v *vec3dFile) buildWireframeData() ([]float32, int) {
	// Count total patches
	totalPatches := 0
	for _, f := range v.fields {
		totalPatches += len(f.patches)
	}

	// 12 edges per box, 2 vertices per edge = 24 vertices per box
	// 7 floats per vertex (x,y,z, r,g,b,a)
	data := make([]float32, 0, totalPatches*24*7)

	invNx := float32(1.0) / float32(v.Nx)
	invNy := float32(1.0) / float32(v.Ny)
	invNz := float32(1.0) / float32(v.Nz)

	for _, field := range v.fields {
		// Get color for this field term
		ft := int(field.fieldTerm)
		var cr, cg, cb float32
		if ft < len(fieldTermColors) {
			cr = fieldTermColors[ft][0]
			cg = fieldTermColors[ft][1]
			cb = fieldTermColors[ft][2]
		} else {
			cr, cg, cb = 1.0, 1.0, 1.0
		}

		for _, p := range field.patches {
			// Find max absolute coefficient for opacity
			var maxCoeff float32
			for _, c := range p.coeffs {
				ac := f32abs(c)
				if ac > maxCoeff {
					maxCoeff = ac
				}
			}

			// Box corners in normalized [0,1] space
			x0 := float32(p.bx) * invNx
			y0 := float32(p.by) * invNy
			z0 := float32(p.bz) * invNz
			x1 := float32(int(p.bx)+int(p.sx)) * invNx
			y1 := float32(int(p.by)+int(p.sy)) * invNy
			z1 := float32(int(p.bz)+int(p.sz)) * invNz

			// Alpha proportional to max coefficient magnitude
			// Use sqrt for perceptual scaling, clamp to [0.05, 1.0]
			alpha := float32(math.Sqrt(float64(maxCoeff)))
			if alpha > 1.0 {
				alpha = 1.0
			}
			if alpha < 0.05 {
				alpha = 0.05
			}

			// 8 corners of the box
			corners := [8][3]float32{
				{x0, y0, z0}, {x1, y0, z0}, {x1, y1, z0}, {x0, y1, z0},
				{x0, y0, z1}, {x1, y0, z1}, {x1, y1, z1}, {x0, y1, z1},
			}

			// 12 edges: bottom face(4), top face(4), verticals(4)
			edges := [12][2]int{
				{0, 1}, {1, 2}, {2, 3}, {3, 0}, // bottom
				{4, 5}, {5, 6}, {6, 7}, {7, 4}, // top
				{0, 4}, {1, 5}, {2, 6}, {3, 7}, // verticals
			}

			for _, e := range edges {
				c0 := corners[e[0]]
				c1 := corners[e[1]]
				data = append(data, c0[0], c0[1], c0[2], cr, cg, cb, alpha)
				data = append(data, c1[0], c1[1], c1[2], cr, cg, cb, alpha)
			}
		}
	}

	return data, totalPatches * 24
}

// ============================================================
// VecStream file reader (pure Go, no CGO)
// ============================================================

const (
	vsMaxOrder = 3
	vsNCoeffs  = 64 // (vsMaxOrder+1)^3

	vsFrameI = 0
	vsFrameP = 1
	vsFrameK = 2
)

type vsPatch struct {
	originX, originY, originZ int16
	sizeX, sizeY, sizeZ       uint8
	order                      uint8
	nCoeffs                    uint16
	coeffs                     [vsNCoeffs]float32
}

type vsIndexEntry struct {
	time      float64
	offset    uint64
	frameType uint8
	fieldIdx  uint8
}

type vecstreamFile struct {
	fp *os.File

	// Header
	version              uint32
	Nx, Ny, Nz           uint32
	Lx, Ly, Lz, dt      float64
	nFields              uint16
	blockSize            uint16
	nFrames              uint32
	flags                uint32

	// Derived
	blocksX, blocksY, blocksZ uint32
	totalPatches              uint32

	// Index
	index []vsIndexEntry

	// Reconstruction state: cached patches per field
	// fieldPatches[fieldIdx] holds the current accumulated patch state
	fieldPatches [][]vsPatch
	// fieldBase[fieldIdx] = index entry of last loaded I/K-frame
	fieldBase []int
	// fieldApplied[fieldIdx] = index entry of last applied P-frame
	fieldApplied []int

	// Pre-allocated voxel buffer for reconstruction
	voxelBuf []float32
}

func vecstreamOpen(path string) (*vecstreamFile, error) {
	fp, err := os.Open(path)
	if err != nil {
		return nil, err
	}

	vs := &vecstreamFile{fp: fp}

	// Read 64-byte file header
	var hdr [64]byte
	if _, err := io.ReadFull(fp, hdr[:]); err != nil {
		fp.Close()
		return nil, fmt.Errorf("read header: %w", err)
	}
	if string(hdr[0:4]) != "VECS" {
		fp.Close()
		return nil, fmt.Errorf("not a vecstream file (got %q)", hdr[0:4])
	}

	vs.version = binary.LittleEndian.Uint32(hdr[4:8])
	vs.Nx = binary.LittleEndian.Uint32(hdr[8:12])
	vs.Ny = binary.LittleEndian.Uint32(hdr[12:16])
	vs.Nz = binary.LittleEndian.Uint32(hdr[16:20])
	vs.Lx = math.Float64frombits(binary.LittleEndian.Uint64(hdr[20:28]))
	vs.Ly = math.Float64frombits(binary.LittleEndian.Uint64(hdr[28:36]))
	vs.Lz = math.Float64frombits(binary.LittleEndian.Uint64(hdr[36:44]))
	vs.dt = math.Float64frombits(binary.LittleEndian.Uint64(hdr[44:52]))
	vs.nFields = binary.LittleEndian.Uint16(hdr[52:54])
	vs.blockSize = binary.LittleEndian.Uint16(hdr[54:56])
	vs.nFrames = binary.LittleEndian.Uint32(hdr[56:60])
	vs.flags = binary.LittleEndian.Uint32(hdr[60:64])

	bs := uint32(vs.blockSize)
	vs.blocksX = (vs.Nx + bs - 1) / bs
	vs.blocksY = (vs.Ny + bs - 1) / bs
	vs.blocksZ = (vs.Nz + bs - 1) / bs
	vs.totalPatches = vs.blocksX * vs.blocksY * vs.blocksZ

	// Read footer (last 16 bytes) to find index
	if _, err := fp.Seek(-16, io.SeekEnd); err != nil {
		fp.Close()
		return nil, fmt.Errorf("seek footer: %w", err)
	}
	var footer [16]byte
	if _, err := io.ReadFull(fp, footer[:]); err != nil {
		fp.Close()
		return nil, fmt.Errorf("read footer: %w", err)
	}
	if string(footer[12:16]) != "VSND" {
		fp.Close()
		return nil, fmt.Errorf("bad footer magic (got %q)", footer[12:16])
	}

	indexOffset := binary.LittleEndian.Uint64(footer[0:8])

	// Read index header (16 bytes)
	if _, err := fp.Seek(int64(indexOffset), io.SeekStart); err != nil {
		fp.Close()
		return nil, fmt.Errorf("seek index: %w", err)
	}
	var ixHdr [16]byte
	if _, err := io.ReadFull(fp, ixHdr[:]); err != nil {
		fp.Close()
		return nil, fmt.Errorf("read index header: %w", err)
	}
	if string(ixHdr[0:4]) != "IXVS" {
		fp.Close()
		return nil, fmt.Errorf("bad index magic (got %q)", ixHdr[0:4])
	}

	nEntries := binary.LittleEndian.Uint32(ixHdr[4:8])
	vs.index = make([]vsIndexEntry, nEntries)

	// Read index entries (24 bytes each)
	for i := uint32(0); i < nEntries; i++ {
		var entry [24]byte
		if _, err := io.ReadFull(fp, entry[:]); err != nil {
			break
		}
		vs.index[i].time = math.Float64frombits(binary.LittleEndian.Uint64(entry[0:8]))
		vs.index[i].offset = binary.LittleEndian.Uint64(entry[8:16])
		vs.index[i].frameType = entry[16]
		vs.index[i].fieldIdx = entry[17]
	}

	// Initialize per-field reconstruction state
	nf := int(vs.nFields)
	vs.fieldPatches = make([][]vsPatch, nf)
	vs.fieldBase = make([]int, nf)
	vs.fieldApplied = make([]int, nf)
	for i := 0; i < nf; i++ {
		vs.fieldPatches[i] = make([]vsPatch, vs.totalPatches)
		vs.fieldBase[i] = -1
		vs.fieldApplied[i] = -1
	}

	// Pre-allocate voxel buffer
	vs.voxelBuf = make([]float32, uint64(vs.Nx)*uint64(vs.Ny)*uint64(vs.Nz))

	fmt.Printf("VecStream: %dx%dx%d, %d fields, block=%d, %d index entries, %d total patches\n",
		vs.Nx, vs.Ny, vs.Nz, vs.nFields, vs.blockSize, nEntries, vs.totalPatches)
	return vs, nil
}

func (vs *vecstreamFile) close() {
	if vs.fp != nil {
		vs.fp.Close()
	}
}

// readFramePayload reads and decompresses a frame payload at the given index entry.
func (vs *vecstreamFile) readFramePayload(idx int) ([]byte, error) {
	entry := vs.index[idx]

	// Seek to frame header
	if _, err := vs.fp.Seek(int64(entry.offset), io.SeekStart); err != nil {
		return nil, fmt.Errorf("seek frame %d: %w", idx, err)
	}

	// Read 32-byte frame header
	var hdr [32]byte
	if _, err := io.ReadFull(vs.fp, hdr[:]); err != nil {
		return nil, fmt.Errorf("read frame header %d: %w", idx, err)
	}
	if string(hdr[0:4]) != "FRVS" {
		return nil, fmt.Errorf("bad frame magic at index %d (got %q)", idx, hdr[0:4])
	}

	compSize := binary.LittleEndian.Uint64(hdr[16:24])
	// uncompSize := binary.LittleEndian.Uint64(hdr[24:32])

	// Read compressed data
	comp := make([]byte, compSize)
	if _, err := io.ReadFull(vs.fp, comp); err != nil {
		return nil, fmt.Errorf("read frame payload %d: %w", idx, err)
	}

	// Decompress
	data, err := zstdDecoder.DecodeAll(comp, nil)
	if err != nil {
		return nil, fmt.Errorf("zstd decompress frame %d: %w", idx, err)
	}
	return data, nil
}

// parseIframe parses decompressed I-frame payload into patches.
func parseIframe(data []byte, patches []vsPatch) (uint32, error) {
	if len(data) < 4 {
		return 0, fmt.Errorf("I-frame data too short")
	}
	np := binary.LittleEndian.Uint32(data[0:4])
	p := data[4:]

	for i := uint32(0); i < np; i++ {
		if len(p) < 16 {
			return 0, fmt.Errorf("I-frame truncated at patch %d", i)
		}
		if int(i) >= len(patches) {
			// skip extra patches beyond our capacity
			nc := binary.LittleEndian.Uint16(p[10:12])
			p = p[16+int(nc)*4:]
			continue
		}
		vp := &patches[i]
		vp.originX = int16(binary.LittleEndian.Uint16(p[0:2]))
		vp.originY = int16(binary.LittleEndian.Uint16(p[2:4]))
		vp.originZ = int16(binary.LittleEndian.Uint16(p[4:6]))
		vp.sizeX = p[6]
		vp.sizeY = p[7]
		vp.sizeZ = p[8]
		vp.order = p[9]
		vp.nCoeffs = binary.LittleEndian.Uint16(p[10:12])
		// skip max_err_f16 and rms_err_f16 (bytes 12-15) -- not needed for reconstruction
		p = p[16:]

		nc := int(vp.nCoeffs)
		if nc > vsNCoeffs {
			nc = vsNCoeffs
		}
		// Clear all coefficients first
		vp.coeffs = [vsNCoeffs]float32{}
		for c := 0; c < nc; c++ {
			vp.coeffs[c] = math.Float32frombits(binary.LittleEndian.Uint32(p[c*4 : c*4+4]))
		}
		p = p[int(vp.nCoeffs)*4:]
	}
	return np, nil
}

// applyPframe applies P-frame deltas to existing patches.
func applyPframe(data []byte, patches []vsPatch) error {
	if len(data) < 8 {
		return fmt.Errorf("P-frame data too short")
	}
	nd := binary.LittleEndian.Uint32(data[0:4])
	cpp := int(binary.LittleEndian.Uint16(data[4:6]))
	// 2 bytes padding
	p := data[8:]

	for i := uint32(0); i < nd; i++ {
		if len(p) < 4+cpp*4 {
			return fmt.Errorf("P-frame truncated at delta %d", i)
		}
		pidx := binary.LittleEndian.Uint32(p[0:4])
		p = p[4:]

		if int(pidx) < len(patches) {
			nc := int(patches[pidx].nCoeffs)
			if nc > cpp {
				nc = cpp
			}
			for c := 0; c < nc; c++ {
				delta := math.Float32frombits(binary.LittleEndian.Uint32(p[c*4 : c*4+4]))
				patches[pidx].coeffs[c] += delta
			}
		}
		p = p[cpp*4:]
	}
	return nil
}

// parseKframe parses K-frame raw voxels into the output buffer.
func (vs *vecstreamFile) parseKframe(data []byte, output []float32) error {
	n3 := uint64(vs.Nx) * uint64(vs.Ny) * uint64(vs.Nz)
	expected := n3 * 4
	if uint64(len(data)) < expected {
		return fmt.Errorf("K-frame data too small (%d < %d)", len(data), expected)
	}
	for i := uint64(0); i < n3; i++ {
		output[i] = math.Float32frombits(binary.LittleEndian.Uint32(data[i*4 : i*4+4]))
	}
	return nil
}

// evalPatch evaluates a tensor-product polynomial patch at local grid point (di,dj,dk).
func evalPatch(p *vsPatch, di, dj, dk int) float32 {
	o1 := int(p.order) + 1

	bsx := int(p.sizeX)
	bsy := int(p.sizeY)
	bsz := int(p.sizeZ)

	var tx, ty, tz float64
	if bsx > 1 {
		tx = float64(di) / float64(bsx-1)
	}
	if bsy > 1 {
		ty = float64(dj) / float64(bsy-1)
	}
	if bsz > 1 {
		tz = float64(dk) / float64(bsz-1)
	}

	// Precompute powers (up to order 3)
	txa := [4]float64{1, tx, tx * tx, tx * tx * tx}
	tya := [4]float64{1, ty, ty * ty, ty * ty * ty}
	tza := [4]float64{1, tz, tz * tz, tz * tz * tz}

	var val float64
	for a := 0; a < o1; a++ {
		for b := 0; b < o1; b++ {
			ab := txa[a] * tya[b]
			base := a*o1*o1 + b*o1
			for c := 0; c < o1; c++ {
				val += float64(p.coeffs[base+c]) * ab * tza[c]
			}
		}
	}
	return float32(val)
}

// patchesToVoxels reconstructs a full voxel grid from patches.
func (vs *vecstreamFile) patchesToVoxels(patches []vsPatch, output []float32) {
	n3 := int(vs.Nx) * int(vs.Ny) * int(vs.Nz)
	for i := 0; i < n3; i++ {
		output[i] = 0
	}

	ny := int(vs.Ny)
	nz := int(vs.Nz)
	nx := int(vs.Nx)

	nWorkers := runtime.GOMAXPROCS(0)
	if nWorkers < 1 {
		nWorkers = 1
	}
	np := len(patches)
	chunkSize := (np + nWorkers - 1) / nWorkers

	var wg sync.WaitGroup
	for w := 0; w < nWorkers; w++ {
		wg.Add(1)
		go func(w int) {
			defer wg.Done()
			start := w * chunkSize
			end := start + chunkSize
			if end > np {
				end = np
			}
			for pi := start; pi < end; pi++ {
				p := &patches[pi]
				ox := int(p.originX)
				oy := int(p.originY)
				oz := int(p.originZ)
				bsx := int(p.sizeX)
				bsy := int(p.sizeY)
				bsz := int(p.sizeZ)

				for di := 0; di < bsx; di++ {
					gi := ox + di
					if gi < 0 || gi >= nx {
						continue
					}
					for dj := 0; dj < bsy; dj++ {
						gj := oy + dj
						if gj < 0 || gj >= ny {
							continue
						}
						for dk := 0; dk < bsz; dk++ {
							gk := oz + dk
							if gk < 0 || gk >= nz {
								continue
							}
							idx := gi*ny*nz + gj*nz + gk
							output[idx] = evalPatch(p, di, dj, dk)
						}
					}
				}
			}
		}(w)
	}
	wg.Wait()
}

// vsTimestampFrames returns the list of unique simulation times (deduplicated across fields).
// It also returns for each unique time the list of index entries.
func (vs *vecstreamFile) uniqueTimes() []float64 {
	seen := make(map[float64]bool)
	var times []float64
	for _, e := range vs.index {
		if !seen[e.time] {
			seen[e.time] = true
			times = append(times, e.time)
		}
	}
	return times
}

// reconstruct reconstructs field fieldIdx at the given index entry position.
// It uses cached patches when possible (incremental P-frame application).
func (vs *vecstreamFile) reconstruct(targetIdx int, fieldIdx int) ([]float32, error) {
	fi := fieldIdx

	// Find the most recent I-frame or K-frame for this field at or before targetIdx
	baseIdx := -1
	for i := targetIdx; i >= 0; i-- {
		e := vs.index[i]
		if int(e.fieldIdx) != fi {
			continue
		}
		if e.frameType == vsFrameI || e.frameType == vsFrameK {
			baseIdx = i
			break
		}
	}
	if baseIdx < 0 {
		return nil, fmt.Errorf("no I/K-frame found for field %d at or before index %d", fi, targetIdx)
	}

	// Check if we can reuse cached patches
	patches := vs.fieldPatches[fi]
	needBase := true
	if vs.fieldBase[fi] == baseIdx {
		// Same base -- can we just apply new P-frames?
		needBase = false
	}

	if needBase {
		entry := vs.index[baseIdx]
		if entry.frameType == vsFrameK {
			// K-frame: decompress raw voxels directly
			data, err := vs.readFramePayload(baseIdx)
			if err != nil {
				return nil, fmt.Errorf("read K-frame: %w", err)
			}
			if err := vs.parseKframe(data, vs.voxelBuf); err != nil {
				return nil, fmt.Errorf("parse K-frame: %w", err)
			}
			// If target is the K-frame itself, return voxels directly
			if baseIdx == targetIdx {
				vs.fieldBase[fi] = baseIdx
				vs.fieldApplied[fi] = baseIdx
				return vs.voxelBuf, nil
			}
			// Need to re-vectorize the K-frame into patches for P-frame application
			// For simplicity, just read the I-frame and P-frames starting fresh
			// Actually, for K-frames followed by P-frames, we need patches.
			// Re-fit from raw voxels (expensive). Better: just recompute from voxels directly.
			// For now: return the K-frame voxels as-is and warn if P-frames follow.
			// In practice, K-frames are usually standalone or at the end.
			vs.fieldBase[fi] = baseIdx
			vs.fieldApplied[fi] = baseIdx
			// We'll skip P-frames after K-frames for now -- reconstruct from voxels directly
			return vs.voxelBuf, nil
		}

		// I-frame: load patches
		data, err := vs.readFramePayload(baseIdx)
		if err != nil {
			return nil, fmt.Errorf("read I-frame: %w", err)
		}
		if _, err := parseIframe(data, patches); err != nil {
			return nil, fmt.Errorf("parse I-frame: %w", err)
		}
		vs.fieldBase[fi] = baseIdx
		vs.fieldApplied[fi] = baseIdx
	}

	// Apply P-frames from after the last applied frame to targetIdx
	startFrom := vs.fieldApplied[fi] + 1
	for i := startFrom; i <= targetIdx; i++ {
		e := vs.index[i]
		if int(e.fieldIdx) != fi {
			continue
		}
		if e.frameType != vsFrameP {
			continue
		}
		data, err := vs.readFramePayload(i)
		if err != nil {
			return nil, fmt.Errorf("read P-frame %d: %w", i, err)
		}
		if err := applyPframe(data, patches); err != nil {
			return nil, fmt.Errorf("apply P-frame %d: %w", i, err)
		}
	}
	vs.fieldApplied[fi] = targetIdx

	// Evaluate patches to voxels
	vs.patchesToVoxels(patches, vs.voxelBuf)
	return vs.voxelBuf, nil
}

// findDisplayFrames returns the list of logical "display frames" for a vecstream file.
// A display frame is a time step where we can show reconstructed data.
// We group index entries by time and return the unique times in order.
// For each logical frame, we pick the first index entry at that time for field 0.
type vsDisplayFrame struct {
	time     float64
	indexIdx int   // index of the entry for the primary field (field 0) at this time
	entryIdx []int // all index entries at this time
}

func (vs *vecstreamFile) buildDisplayFrames() []vsDisplayFrame {
	// Group entries by time
	type timeGroup struct {
		time    float64
		entries []int
	}
	var groups []timeGroup
	timeMap := make(map[float64]int) // time -> index in groups

	for i, e := range vs.index {
		gi, ok := timeMap[e.time]
		if !ok {
			gi = len(groups)
			timeMap[e.time] = gi
			groups = append(groups, timeGroup{time: e.time})
		}
		groups[gi].entries = append(groups[gi].entries, i)
	}

	// Build display frames
	frames := make([]vsDisplayFrame, len(groups))
	for i, g := range groups {
		frames[i].time = g.time
		frames[i].entryIdx = g.entries
		// Pick the entry for field 0 (or the first entry if field 0 not present)
		frames[i].indexIdx = g.entries[0]
		for _, ei := range g.entries {
			if vs.index[ei].fieldIdx == 0 {
				frames[i].indexIdx = ei
				break
			}
		}
	}
	return frames
}

// frameTypeName returns "I", "P", or "K" for display.
func vsFrameTypeName(ft uint8) string {
	switch ft {
	case vsFrameI:
		return "I"
	case vsFrameP:
		return "P"
	case vsFrameK:
		return "K"
	default:
		return "?"
	}
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

// f32cbrt returns the cube root of a float32.
func f32cbrt(x float32) float32 {
	return float32(math.Cbrt(float64(x)))
}

// PGA component → hue mapping. 8 components evenly spaced on the
// color wheel so that even a "soup" of all 8 active still produces a
// distinguishable RGB. Rotor (M0..M3) skews warm; spatial (M4..M7)
// skews cool, so the bulk-Higgs vacuum (M[4]≈v, others≈0) reads as
// cyan and rotor activity reads as warm hues.
var pgaHueColors = [8][3]float32{
	{1.00, 0.00, 0.00}, // M0 e₄₁₀  red
	{1.00, 0.40, 0.00}, // M1 e₄₂₀  orange
	{1.00, 1.00, 0.00}, // M2 e₄₃₀  yellow
	{0.30, 1.00, 0.00}, // M3 𝟙     green
	{0.00, 1.00, 0.80}, // M4 e₁    cyan
	{0.00, 0.40, 1.00}, // M5 e₂    blue
	{0.55, 0.10, 1.00}, // M6 e₃    violet
	{1.00, 0.00, 0.70}, // M7 e₃₂₁  magenta
}

// computePGAView fills RGBA from the 8 PGA components.
//
//	mode 0 (spectrum):  every component contributes |M[k] − mean_k|
//	                    weighted by its hue color; sum → RGB. Shows
//	                    DEVIATION from per-frame mean — vacuum-like
//	                    constant fields contribute nothing, structure
//	                    in any of the 8 components becomes visible.
//	mode 1 (rotation):  3-of-8 components mapped 1:1 to RGB, starting
//	                    at `phase` (cycles M[phase], M[phase+1], M[phase+2]).
//	                    Also shows |M[k] − mean_k|.
//
// Per-component normalisation by max |deviation| so rotor (range ~2)
// and spatial (range ~v≈0.5) sectors share the same display scale.
func computePGAView(n int, cols [8][]float32, vol *volumeData, mode int, phase int) {
	total := n * n * n
	if len(vol.rgba) != total*4 {
		vol.rgba = make([]float32, total*4)
	}

	// Per-component mean and max |deviation|.
	var mean [8]float32
	for k := 0; k < 8; k++ {
		var s float64
		for i := 0; i < total; i++ {
			s += float64(cols[k][i])
		}
		mean[k] = float32(s / float64(total))
	}
	var maxDev [8]float32
	for k := 0; k < 8; k++ {
		var m float32
		mk := mean[k]
		for i := 0; i < total; i++ {
			v := f32abs(cols[k][i] - mk)
			if v > m {
				m = v
			}
		}
		maxDev[k] = m
	}
	var invNorm [8]float32
	for k := 0; k < 8; k++ {
		if maxDev[k] > 1e-12 {
			invNorm[k] = 1.0 / maxDev[k]
		}
	}

	if mode == 0 {
		// Spectrum: weighted sum over all 8 component hues. Factor 0.4
		// keeps a soup of all 8 active components from blowing out the
		// shader clamp; single dominant components stay punchy.
		const scale = float32(0.40)
		for i := 0; i < total; i++ {
			var r, g, b float32
			for k := 0; k < 8; k++ {
				a := f32abs(cols[k][i]-mean[k]) * invNorm[k]
				r += a * pgaHueColors[k][0]
				g += a * pgaHueColors[k][1]
				b += a * pgaHueColors[k][2]
			}
			vol.rgba[i*4+0] = r * scale
			vol.rgba[i*4+1] = g * scale
			vol.rgba[i*4+2] = b * scale
			vol.rgba[i*4+3] = 1.0
		}
	} else {
		// Rotation: select 3 consecutive components (mod 8) → RGB.
		kr := ((phase % 8) + 8) % 8
		kg := (kr + 1) % 8
		kb := (kr + 2) % 8
		for i := 0; i < total; i++ {
			vol.rgba[i*4+0] = f32abs(cols[kr][i]-mean[kr]) * invNorm[kr]
			vol.rgba[i*4+1] = f32abs(cols[kg][i]-mean[kg]) * invNorm[kg]
			vol.rgba[i*4+2] = f32abs(cols[kb][i]-mean[kb]) * invNorm[kb]
			vol.rgba[i*4+3] = 1.0
		}
	}
	vol.n = n
	vol.loaded = true
}

// computeFieldView fills RGBA volume from SFA column data.
// Uses pre-allocated intermediate buffers from sfa.absPBuf/phi2Buf/theta2Buf
// to avoid reading the input arrays twice.
func computeFieldView(n int, phi0, phi1, phi2, theta0, theta1, theta2 []float32, vol *volumeData, fixedMax float32) {
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
	invNormPhi2 := float32(1)
	invNormTheta2 := float32(1)
	if fixedMax > 0 {
		// Fixed scale: all channels use the same user-specified max
		invNormP = 1.0 / fixedMax
		invNormPhi2 = 1.0 / fixedMax
		invNormTheta2 = 1.0 / fixedMax
	} else {
		// Auto-scale from data max
		if maxP*0.3 > 1e-20 {
			invNormP = 1.0 / (maxP * 0.3)
		}
		if maxPhi2*0.15 > 1e-20 {
			invNormPhi2 = 1.0 / (maxPhi2 * 0.15)
		}
		if maxTheta2*0.3 > 1e-20 {
			invNormTheta2 = 1.0 / (maxTheta2 * 0.3)
		}
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

// computeVelocityView fills RGBA volume from velocity columns.
// R = |v| = sqrt(vphi0^2 + vphi1^2 + vphi2^2)
// G = |vphi0| (x-velocity component)
// B = |vphi1| (y-velocity component)
// Each channel normalized independently.
func computeVelocityView(n int, phi0, phi1, phi2, vphi0, vphi1, vphi2 []float32, vol *volumeData) {
	total := n * n * n
	vol.n = n
	vol.loaded = true
	if len(vol.rgba) != total*4 {
		vol.rgba = make([]float32, total*4)
	}

	nWorkers := runtime.GOMAXPROCS(0)
	if nWorkers < 1 {
		nWorkers = 1
	}
	chunkSize := (total + nWorkers - 1) / nWorkers

	// Reuse intermediate buffers: absPBuf=|v|, phi2Buf=|vphi0|, theta2Buf=|vphi1|
	absPBuf := vol.absPBuf
	phi2Buf := vol.phi2Buf
	theta2Buf := vol.theta2Buf
	if len(absPBuf) < total {
		absPBuf = make([]float32, total)
		phi2Buf = make([]float32, total)
		theta2Buf = make([]float32, total)
	}

	// Pass 1: compute intermediates and find max
	type maxVals struct {
		maxVmag, maxVx, maxVy float32
		_pad                  [4]float32
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
			var mvm, mvx, mvy float32
			vm := absPBuf[start:end]
			vx := phi2Buf[start:end]
			vy := theta2Buf[start:end]
			v0 := vphi0[start:end]
			v1 := vphi1[start:end]
			v2 := vphi2[start:end]

			for i := range vm {
				mag := f32sqrt(v0[i]*v0[i] + v1[i]*v1[i] + v2[i]*v2[i])
				vm[i] = mag
				if mag > mvm {
					mvm = mag
				}
				ax := f32abs(v0[i])
				vx[i] = ax
				if ax > mvx {
					mvx = ax
				}
				ay := f32abs(v1[i])
				vy[i] = ay
				if ay > mvy {
					mvy = ay
				}
			}
			results[w] = maxVals{maxVmag: mvm, maxVx: mvx, maxVy: mvy}
		}(w)
	}
	wg.Wait()

	var maxVmag, maxVx, maxVy float32
	for _, r := range results {
		if r.maxVmag > maxVmag {
			maxVmag = r.maxVmag
		}
		if r.maxVx > maxVx {
			maxVx = r.maxVx
		}
		if r.maxVy > maxVy {
			maxVy = r.maxVy
		}
	}

	invVmag := float32(1)
	if maxVmag*0.3 > 1e-20 {
		invVmag = 1.0 / (maxVmag * 0.3)
	}
	invVx := float32(1)
	if maxVx*0.3 > 1e-20 {
		invVx = 1.0 / (maxVx * 0.3)
	}
	invVy := float32(1)
	if maxVy*0.3 > 1e-20 {
		invVy = 1.0 / (maxVy * 0.3)
	}

	// Pass 2: normalize to RGBA
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
			vm := absPBuf[start:end]
			vx := phi2Buf[start:end]
			vy := theta2Buf[start:end]
			out := rgba[start*4 : end*4]

			for i := range vm {
				j := i * 4
				out[j+0] = vm[i] * invVmag
				out[j+1] = vx[i] * invVx
				out[j+2] = vy[i] * invVy
				out[j+3] = 1.0
			}
		}(w)
	}
	wg.Wait()
}

// computeAccelView fills RGBA volume with a combined force/energy view.
// R = |P|^(1/3) — cube root enhances weak binding regions
// G = theta_rms^2 with enhanced contrast
// B = |v| velocity magnitude
func computeAccelView(n int, phi0, phi1, phi2, theta0, theta1, theta2, vphi0, vphi1, vphi2 []float32, vol *volumeData) {
	total := n * n * n
	vol.n = n
	vol.loaded = true
	if len(vol.rgba) != total*4 {
		vol.rgba = make([]float32, total*4)
	}

	hasTheta := theta0 != nil && theta1 != nil && theta2 != nil
	hasVel := vphi0 != nil && vphi1 != nil && vphi2 != nil

	nWorkers := runtime.GOMAXPROCS(0)
	if nWorkers < 1 {
		nWorkers = 1
	}
	chunkSize := (total + nWorkers - 1) / nWorkers

	// Reuse intermediate buffers: absPBuf=|P|^(1/3), phi2Buf=theta_rms^2, theta2Buf=|v|
	absPBuf := vol.absPBuf
	phi2Buf := vol.phi2Buf
	theta2Buf := vol.theta2Buf
	if len(absPBuf) < total {
		absPBuf = make([]float32, total)
		phi2Buf = make([]float32, total)
		theta2Buf = make([]float32, total)
	}

	// Pass 1: compute intermediates and find max
	type maxVals struct {
		maxCbrtP, maxTheta2, maxVmag float32
		_pad                         [4]float32
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
			var mcp, mt2, mvm float32
			cp := absPBuf[start:end]
			t2 := phi2Buf[start:end]
			vm := theta2Buf[start:end]
			p0 := phi0[start:end]
			p1 := phi1[start:end]
			p2 := phi2[start:end]

			for i := range cp {
				absP := f32abs(p0[i] * p1[i] * p2[i])
				cbrtP := f32cbrt(absP)
				cp[i] = cbrtP
				if cbrtP > mcp {
					mcp = cbrtP
				}
			}
			if hasTheta {
				t0s := theta0[start:end]
				t1s := theta1[start:end]
				t2s := theta2[start:end]
				for i := range t2 {
					tv := t0s[i]*t0s[i] + t1s[i]*t1s[i] + t2s[i]*t2s[i]
					t2[i] = tv
					if tv > mt2 {
						mt2 = tv
					}
				}
			}
			if hasVel {
				v0 := vphi0[start:end]
				v1 := vphi1[start:end]
				v2 := vphi2[start:end]
				for i := range vm {
					mag := f32sqrt(v0[i]*v0[i] + v1[i]*v1[i] + v2[i]*v2[i])
					vm[i] = mag
					if mag > mvm {
						mvm = mag
					}
				}
			}
			results[w] = maxVals{maxCbrtP: mcp, maxTheta2: mt2, maxVmag: mvm}
		}(w)
	}
	wg.Wait()

	var maxCbrtP, maxTheta2, maxVmag float32
	for _, r := range results {
		if r.maxCbrtP > maxCbrtP {
			maxCbrtP = r.maxCbrtP
		}
		if r.maxTheta2 > maxTheta2 {
			maxTheta2 = r.maxTheta2
		}
		if r.maxVmag > maxVmag {
			maxVmag = r.maxVmag
		}
	}

	// Stronger normalization for |P|^(1/3) to enhance weak binding
	invCbrtP := float32(1)
	if maxCbrtP*0.1 > 1e-20 {
		invCbrtP = 1.0 / (maxCbrtP * 0.1)
	}
	invTheta2 := float32(1)
	if maxTheta2*0.15 > 1e-20 {
		invTheta2 = 1.0 / (maxTheta2 * 0.15)
	}
	invVmag := float32(1)
	if maxVmag*0.3 > 1e-20 {
		invVmag = 1.0 / (maxVmag * 0.3)
	}

	// Pass 2: normalize to RGBA
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
			cp := absPBuf[start:end]
			t2 := phi2Buf[start:end]
			vm := theta2Buf[start:end]
			out := rgba[start*4 : end*4]

			for i := range cp {
				j := i * 4
				out[j+0] = cp[i] * invCbrtP
				out[j+1] = t2[i] * invTheta2
				out[j+2] = vm[i] * invVmag
				out[j+3] = 1.0
			}
		}(w)
	}
	wg.Wait()
}

// ============================================================
// Complex / U(1)-gauge aware views (v66+ kernels, 24/30-column files)
// ============================================================

// findColumn returns the index of the named column, or -1.
func (s *sfaFile) findColumn(name string) int {
	for i := 0; i < int(s.nColumns); i++ {
		if strings.TrimRight(string(s.columns[i].name[:]), "\x00") == name {
			return i
		}
	}
	return -1
}

// extCols holds name-resolved indices for the complexified (phiim_*) and
// gauged (th_*, E_*) column sets. Index -1 = absent.
type extCols struct {
	im, thim, vre, vim, tvre, tvim, link, ef [3]int
	hasIm, hasThIm, hasVel, hasThVel, hasGauge bool
}

func resolveExtCols(s *sfaFile) extCols {
	var e extCols
	names := func(dst *[3]int, has *bool, n0, n1, n2 string) {
		*has = true
		for a, nm := range [3]string{n0, n1, n2} {
			dst[a] = s.findColumn(nm)
			if dst[a] < 0 {
				*has = false
			}
		}
	}
	names(&e.im, &e.hasIm, "phiim_x", "phiim_y", "phiim_z")
	names(&e.thim, &e.hasThIm, "thetaim_x", "thetaim_y", "thetaim_z")
	names(&e.vre, &e.hasVel, "phi_vx", "phi_vy", "phi_vz")
	var hasVim bool
	names(&e.vim, &hasVim, "phiim_vx", "phiim_vy", "phiim_vz")
	e.hasVel = e.hasVel && hasVim
	names(&e.tvre, &e.hasThVel, "theta_vx", "theta_vy", "theta_vz")
	var hasTvim bool
	names(&e.tvim, &hasTvim, "thetaim_vx", "thetaim_vy", "thetaim_vz")
	e.hasThVel = e.hasThVel && hasTvim
	names(&e.link, &e.hasGauge, "th_x", "th_y", "th_z")
	var hasEf bool
	names(&e.ef, &hasEf, "E_x", "E_y", "E_z")
	e.hasGauge = e.hasGauge && hasEf
	return e
}

// effBuf returns a cached scratch buffer for derived per-voxel arrays.
func (s *sfaFile) effBuf(i, n int) []float32 {
	if len(s.effBufs[i]) < n {
		s.effBufs[i] = make([]float32, n)
	}
	return s.effBufs[i][:n]
}

// parallelChunks runs fn(start, end, worker) over [0,total) in GOMAXPROCS chunks.
func parallelChunks(total int, fn func(start, end, w int)) int {
	nWorkers := runtime.GOMAXPROCS(0)
	if nWorkers < 1 {
		nWorkers = 1
	}
	if nWorkers > 64 {
		nWorkers = 64
	}
	chunk := (total + nWorkers - 1) / nWorkers
	var wg sync.WaitGroup
	for w := 0; w < nWorkers; w++ {
		wg.Add(1)
		go func(w int) {
			defer wg.Done()
			start := w * chunk
			end := start + chunk
			if end > total {
				end = total
			}
			if start < end {
				fn(start, end, w)
			}
		}(w)
	}
	wg.Wait()
	return nWorkers
}

// complexModulus writes sqrt(re^2+im^2) into dst (parallel).
func complexModulus(dst, re, im []float32) {
	parallelChunks(len(dst), func(start, end, _ int) {
		d, r, m := dst[start:end], re[start:end], im[start:end]
		for i := range d {
			d[i] = float32(math.Sqrt(float64(r[i]*r[i] + m[i]*m[i])))
		}
	})
}

// computeGaugeView: U(1) sector. R = |E| (electric field / Coulomb halo),
// G = |Phi|^2 (matter density, dimmed context), B = |A| (link angles).
func computeGaugeView(n int, phi, ef, link [3][]float32, vol *volumeData) {
	total := n * n * n
	vol.n = n
	vol.loaded = true
	if len(vol.rgba) != total*4 {
		vol.rgba = make([]float32, total*4)
	}
	var maxE, maxR, maxA [64]float32
	nw := parallelChunks(total, func(start, end, w int) {
		var me, mr, ma float32
		for i := start; i < end; i++ {
			e := float32(math.Sqrt(float64(ef[0][i]*ef[0][i] + ef[1][i]*ef[1][i] + ef[2][i]*ef[2][i])))
			r := phi[0][i]*phi[0][i] + phi[1][i]*phi[1][i] + phi[2][i]*phi[2][i]
			a := float32(math.Sqrt(float64(link[0][i]*link[0][i] + link[1][i]*link[1][i] + link[2][i]*link[2][i])))
			vol.rgba[i*4+0] = e
			vol.rgba[i*4+1] = r
			vol.rgba[i*4+2] = a
			vol.rgba[i*4+3] = 1.0
			if e > me {
				me = e
			}
			if r > mr {
				mr = r
			}
			if a > ma {
				ma = a
			}
		}
		maxE[w], maxR[w], maxA[w] = me, mr, ma
	})
	var mE, mR, mA float32
	for w := 0; w < nw; w++ {
		if maxE[w] > mE {
			mE = maxE[w]
		}
		if maxR[w] > mR {
			mR = maxR[w]
		}
		if maxA[w] > mA {
			mA = maxA[w]
		}
	}
	invE := float32(0)
	if mE > 1e-30 {
		invE = 1.0 / mE
	}
	invR := float32(0)
	if mR > 1e-30 {
		invR = 0.5 / mR // dimmed: matter is context here
	}
	invA := float32(0)
	if mA > 1e-30 {
		invA = 1.0 / mA
	}
	parallelChunks(total, func(start, end, _ int) {
		for i := start; i < end; i++ {
			vol.rgba[i*4+0] *= invE
			vol.rgba[i*4+1] *= invR
			vol.rgba[i*4+2] *= invA
		}
	})
	fmt.Printf("  gauge view: max|E|=%.3e max|Phi|^2=%.3e max|A|=%.3e\n", mE, mR, mA)
}

// computeChargeView: Noether charge density rhoQ = sum_a (u vdot - v udot)
// (+ theta terms when present). R = +rhoQ, B = -rhoQ, G = |Phi|^2 context.
func computeChargeView(n int, re, im, vre, vim [3][]float32,
	tre, tim, tvre, tvim [3][]float32, hasThetaQ bool, vol *volumeData) {
	total := n * n * n
	vol.n = n
	vol.loaded = true
	if len(vol.rgba) != total*4 {
		vol.rgba = make([]float32, total*4)
	}
	var maxQ, maxR [64]float32
	nw := parallelChunks(total, func(start, end, w int) {
		var mq, mr float32
		for i := start; i < end; i++ {
			var q, r float32
			for a := 0; a < 3; a++ {
				q += re[a][i]*vim[a][i] - im[a][i]*vre[a][i]
				r += re[a][i]*re[a][i] + im[a][i]*im[a][i]
			}
			if hasThetaQ {
				for a := 0; a < 3; a++ {
					q += tre[a][i]*tvim[a][i] - tim[a][i]*tvre[a][i]
				}
			}
			vol.rgba[i*4+0] = q // signed for now; split after normalize
			vol.rgba[i*4+1] = r
			vol.rgba[i*4+3] = 1.0
			aq := q
			if aq < 0 {
				aq = -aq
			}
			if aq > mq {
				mq = aq
			}
			if r > mr {
				mr = r
			}
		}
		maxQ[w], maxR[w] = mq, mr
	})
	var mQ, mR float32
	for w := 0; w < nw; w++ {
		if maxQ[w] > mQ {
			mQ = maxQ[w]
		}
		if maxR[w] > mR {
			mR = maxR[w]
		}
	}
	invQ := float32(0)
	if mQ > 1e-30 {
		invQ = 1.0 / mQ
	}
	invR := float32(0)
	if mR > 1e-30 {
		invR = 0.3 / mR // faint context
	}
	parallelChunks(total, func(start, end, _ int) {
		for i := start; i < end; i++ {
			q := vol.rgba[i*4+0] * invQ
			if q >= 0 {
				vol.rgba[i*4+0] = q
				vol.rgba[i*4+2] = 0
			} else {
				vol.rgba[i*4+0] = 0
				vol.rgba[i*4+2] = -q
			}
			vol.rgba[i*4+1] *= invR
		}
	})
	fmt.Printf("  charge view: max|rhoQ|=%.3e (R=+Q, B=-Q)\n", mQ)
}

// computeStandardView is the shared column-resolution + view dispatch for
// non-PGA, non-preview files. For complexified files (phiim_* present) the
// field/velocity/accel views receive per-component complex moduli |Phi_a|,
// |Theta_a| — phase-invariant, so Q-balls render as static objects instead
// of flickering at their internal clock frequency.
// Modes: 0=field, 1=velocity, 2=accel, 3=U(1) gauge, 4=charge.
func computeStandardView(s *sfaFile, n int, vol *volumeData, mode int, fixedMax float32) {
	total := n * n * n
	ec := resolveExtCols(s)

	phi := [3][]float32{s.extractColumnF32(0), s.extractColumnF32(1), s.extractColumnF32(2)}
	var theta [3][]float32
	hasTheta := s.nColumns >= 6
	if hasTheta {
		theta = [3][]float32{s.extractColumnF32(3), s.extractColumnF32(4), s.extractColumnF32(5)}
	}

	// charge view needs the RAW re/im parts — capture before modulus replacement
	if mode == 4 {
		if !(ec.hasIm && ec.hasVel) {
			fmt.Println("Warning: charge view needs complex columns (phiim_*, *_v*), falling back to field view")
			mode = 0
		} else {
			im := [3][]float32{s.extractColumnF32(ec.im[0]), s.extractColumnF32(ec.im[1]), s.extractColumnF32(ec.im[2])}
			vre := [3][]float32{s.extractColumnF32(ec.vre[0]), s.extractColumnF32(ec.vre[1]), s.extractColumnF32(ec.vre[2])}
			vim := [3][]float32{s.extractColumnF32(ec.vim[0]), s.extractColumnF32(ec.vim[1]), s.extractColumnF32(ec.vim[2])}
			var tre, tim, tvre, tvim [3][]float32
			hasThetaQ := hasTheta && ec.hasThIm && ec.hasThVel
			if hasThetaQ {
				tre = theta
				tim = [3][]float32{s.extractColumnF32(ec.thim[0]), s.extractColumnF32(ec.thim[1]), s.extractColumnF32(ec.thim[2])}
				tvre = [3][]float32{s.extractColumnF32(ec.tvre[0]), s.extractColumnF32(ec.tvre[1]), s.extractColumnF32(ec.tvre[2])}
				tvim = [3][]float32{s.extractColumnF32(ec.tvim[0]), s.extractColumnF32(ec.tvim[1]), s.extractColumnF32(ec.tvim[2])}
			}
			computeChargeView(n, phi, im, vre, vim, tre, tim, tvre, tvim, hasThetaQ, vol)
			return
		}
	}

	// complexified file: replace components by their moduli (phase-invariant)
	if ec.hasIm {
		for a := 0; a < 3; a++ {
			eff := s.effBuf(a, total)
			complexModulus(eff, phi[a], s.extractColumnF32(ec.im[a]))
			phi[a] = eff
		}
		if hasTheta && ec.hasThIm {
			for a := 0; a < 3; a++ {
				eff := s.effBuf(3+a, total)
				complexModulus(eff, theta[a], s.extractColumnF32(ec.thim[a]))
				theta[a] = eff
			}
		}
	}

	switch mode {
	case 1: // velocity
		var v [3][]float32
		if ec.hasIm && ec.hasVel {
			for a := 0; a < 3; a++ {
				eff := s.effBuf(6+a, total)
				complexModulus(eff, s.extractColumnF32(ec.vre[a]), s.extractColumnF32(ec.vim[a]))
				v[a] = eff
			}
		} else if s.nColumns >= 9 {
			v = [3][]float32{s.extractColumnF32(6), s.extractColumnF32(7), s.extractColumnF32(8)}
		}
		if v[0] != nil {
			computeVelocityView(n, phi[0], phi[1], phi[2], v[0], v[1], v[2], vol)
		} else {
			fmt.Println("Warning: no velocity columns (need >=9 columns), falling back to field view")
			computeFieldView(n, phi[0], phi[1], phi[2], theta[0], theta[1], theta[2], vol, fixedMax)
		}
	case 2: // accel (legacy view; complex files get modulus inputs)
		var v0, v1, v2 []float32
		if s.nColumns >= 9 {
			v0, v1, v2 = s.extractColumnF32(6), s.extractColumnF32(7), s.extractColumnF32(8)
		}
		computeAccelView(n, phi[0], phi[1], phi[2], theta[0], theta[1], theta[2], v0, v1, v2, vol)
	case 3: // U(1) gauge
		if !ec.hasGauge {
			fmt.Println("Warning: no gauge columns (th_*, E_*) — not a complex_gauge=1 file; falling back to field view")
			computeFieldView(n, phi[0], phi[1], phi[2], theta[0], theta[1], theta[2], vol, fixedMax)
			return
		}
		ef := [3][]float32{s.extractColumnF32(ec.ef[0]), s.extractColumnF32(ec.ef[1]), s.extractColumnF32(ec.ef[2])}
		link := [3][]float32{s.extractColumnF32(ec.link[0]), s.extractColumnF32(ec.link[1]), s.extractColumnF32(ec.link[2])}
		computeGaugeView(n, phi, ef, link, vol)
	default:
		computeFieldView(n, phi[0], phi[1], phi[2], theta[0], theta[1], theta[2], vol, fixedMax)
	}
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
// GLSL 430 line shaders (for wireframe overlay)
// ============================================================

const lineVertexShaderSource = `#version 430

layout(location = 0) in vec3 inPos;    // position in [0,1]^3
layout(location = 1) in vec4 inColor;  // RGBA

uniform vec3 camPos;
uniform vec3 camTarget;
uniform vec2 winSize;
uniform float fov;

out vec4 vColor;

void main() {
    float aspect = winSize.x / winSize.y;

    // Camera setup (same as volume shader)
    vec3 fwd = normalize(camTarget - camPos);
    vec3 worldUp = vec3(0.0, 0.0, 1.0);
    if (abs(dot(fwd, worldUp)) > 0.999) {
        worldUp = vec3(0.0, 1.0, 0.0);
    }
    vec3 right = normalize(cross(fwd, worldUp));
    vec3 up = cross(right, fwd);

    // View-space position
    vec3 p = inPos - camPos;
    float z = dot(p, fwd);
    float x = dot(p, right);
    float y = dot(p, up);

    // Perspective projection
    float invZ = 1.0 / max(z, 0.001);
    float px = x * invZ / (fov * aspect);
    float py = y * invZ / (fov);

    gl_Position = vec4(px, py, -1.0/max(z, 0.001), 1.0);
    vColor = inColor;
}
` + "\x00"

const lineFragmentShaderSource = `#version 430

in vec4 vColor;
out vec4 fragColor;

void main() {
    fragColor = vColor;
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

	// Scale mode
	fixedScale    bool    // true = use fixedScaleVal instead of auto-ranging
	fixedScaleVal float32 // the fixed max value for color mapping

	// Vec3d overlay
	showVec bool // toggle wireframe overlay

	// State
	viewMode   int
	pgaMode    int  // 0=spectrum, 1=rotation (used when v56 PGA 8-comp file detected)
	pgaPhase   int  // rotation slot (which 3-of-8 components are R/G/B)
	curFrame   int
	playing    bool
	dragging   bool
	lastMouseX float64
	lastMouseY float64

	// Vecstream field cycling
	vsFieldIdx int
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
		fixedScale:    false,
		fixedScaleVal: 0.01,
		showVec:       true,
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
// Line rendering GPU resources (wireframe overlay)
// ============================================================

type lineGPUResources struct {
	program uint32
	vao     uint32
	vbo     uint32
	nVerts  int32

	// Uniform locations
	locCamPos    int32
	locCamTarget int32
	locWinSize   int32
	locFov       int32
}

func initLineGPUResources() (*lineGPUResources, error) {
	lg := &lineGPUResources{}

	vertShader, err := compileShader(lineVertexShaderSource, gl.VERTEX_SHADER)
	if err != nil {
		return nil, fmt.Errorf("line vertex shader: %w", err)
	}
	fragShader, err := compileShader(lineFragmentShaderSource, gl.FRAGMENT_SHADER)
	if err != nil {
		return nil, fmt.Errorf("line fragment shader: %w", err)
	}
	lg.program, err = linkProgram(vertShader, fragShader)
	if err != nil {
		return nil, fmt.Errorf("line program link: %w", err)
	}
	gl.DeleteShader(vertShader)
	gl.DeleteShader(fragShader)

	gl.UseProgram(lg.program)
	lg.locCamPos = gl.GetUniformLocation(lg.program, gl.Str("camPos\x00"))
	lg.locCamTarget = gl.GetUniformLocation(lg.program, gl.Str("camTarget\x00"))
	lg.locWinSize = gl.GetUniformLocation(lg.program, gl.Str("winSize\x00"))
	lg.locFov = gl.GetUniformLocation(lg.program, gl.Str("fov\x00"))

	gl.GenVertexArrays(1, &lg.vao)
	gl.GenBuffers(1, &lg.vbo)

	return lg, nil
}

func (lg *lineGPUResources) uploadWireframe(data []float32, nVerts int) {
	lg.nVerts = int32(nVerts)

	gl.BindVertexArray(lg.vao)
	gl.BindBuffer(gl.ARRAY_BUFFER, lg.vbo)
	gl.BufferData(gl.ARRAY_BUFFER, len(data)*4, unsafe.Pointer(&data[0]), gl.STATIC_DRAW)

	// Stride = 7 floats (pos xyz + color rgba) = 28 bytes
	stride := int32(7 * 4)

	// Attribute 0: position (vec3)
	gl.EnableVertexAttribArray(0)
	gl.VertexAttribPointerWithOffset(0, 3, gl.FLOAT, false, stride, 0)

	// Attribute 1: color (vec4)
	gl.EnableVertexAttribArray(1)
	gl.VertexAttribPointerWithOffset(1, 4, gl.FLOAT, false, stride, uintptr(3*4))

	gl.BindVertexArray(0)
}

func (lg *lineGPUResources) render(p *viewParams, winW, winH float32) {
	if lg.nVerts == 0 {
		return
	}

	gl.UseProgram(lg.program)
	gl.Uniform3f(lg.locCamPos, p.camPos[0], p.camPos[1], p.camPos[2])
	gl.Uniform3f(lg.locCamTarget, p.camTarget[0], p.camTarget[1], p.camTarget[2])
	gl.Uniform2f(lg.locWinSize, winW, winH)
	gl.Uniform1f(lg.locFov, p.fov)

	gl.Enable(gl.BLEND)
	gl.BlendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA)
	gl.Enable(gl.LINE_SMOOTH)
	gl.LineWidth(1.0)

	gl.BindVertexArray(lg.vao)
	gl.DrawArrays(gl.LINES, 0, lg.nVerts)
	gl.BindVertexArray(0)

	gl.Disable(gl.BLEND)
	gl.Disable(gl.LINE_SMOOTH)
}

// ============================================================
// Async frame loader (copied from volview2)
// ============================================================

type frameLoader struct {
	sfa      *sfaFile
	volN     int
	viewMode int        // 0=field, 1=velocity, 2=acceleration
	pgaMode  int        // 0=spectrum, 1=rotation (used for v56 8-component PGA files)
	pgaPhase int        // rotation slot (M[phase..phase+2] → R,G,B)
	fixedMax float32    // >0: fixed color scale; 0: auto-scale from data
	mu       sync.Mutex
	ready    *volumeData // non-nil when a new frame is ready for upload
	readyT   float64     // time of the ready frame
	readyF   int         // frame index of the ready frame
	busy     atomic.Bool // true while loading
	workVol  volumeData  // pre-allocated work buffer for computing volume
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
	mode := fl.viewMode
	fmt.Printf("Loading frame %d/%d (mode %d)...\n", frameIdx+1, fl.sfa.totalFrames, mode)
	go func() {
		defer fl.busy.Store(false)

		frameTime, err := fl.sfa.readFrame(uint32(frameIdx))
		if err != nil {
			fmt.Printf("Failed to load frame %d: %v\n", frameIdx, err)
			return
		}

		// Detect preview mode: 3 columns named P_abs, phi_sq, theta_sq
		colName0 := strings.TrimRight(string(fl.sfa.columns[0].name[:]), "\x00")
		isPreview := fl.sfa.nColumns == 3 && colName0 == "P_abs"
		if isPreview {
			fmt.Printf("  Preview mode: 3 uint8 channels (P_abs, phi_sq, theta_sq)\n")
		}
		// Detect v56 PGA file: 8 columns whose first column is named "M0_…".
		isPGA := fl.sfa.nColumns == 8 && strings.HasPrefix(colName0, "M0")

		if isPGA {
			var pgaCols [8][]float32
			for k := 0; k < 8; k++ {
				pgaCols[k] = fl.sfa.extractColumnF32(k)
			}
			computePGAView(fl.volN, pgaCols, &fl.workVol, fl.pgaMode, fl.pgaPhase)
		} else if isPreview {
			// Preview mode: columns ARE the pre-computed derived quantities.
			// Column 0 = |P| (0-1 normalized), Column 1 = |φ|² (0-1), Column 2 = |θ|² (0-1)
			pAbs := fl.sfa.extractColumnF32(0)
			phiSq := fl.sfa.extractColumnF32(1)
			thetaSq := fl.sfa.extractColumnF32(2)
			n := fl.volN
			vol := &fl.workVol
			vol.n = n
			nn := n * n * n
			if len(vol.rgba) != nn*4 {
				vol.rgba = make([]float32, nn*4)
			}
			// Find per-frame max for normalization (same as computeFieldView)
			var maxP, maxPhi, maxTh float32
			for i := 0; i < nn; i++ {
				if pAbs[i] > maxP { maxP = pAbs[i] }
				if phiSq[i] > maxPhi { maxPhi = phiSq[i] }
				if thetaSq[i] > maxTh { maxTh = thetaSq[i] }
			}
			invP := float32(1); if maxP*0.3 > 1e-20 { invP = 1.0 / (maxP * 0.3) }
			invPhi := float32(1); if maxPhi*0.15 > 1e-20 { invPhi = 1.0 / (maxPhi * 0.15) }
			invTh := float32(1); if maxTh*0.3 > 1e-20 { invTh = 1.0 / (maxTh * 0.3) }

			for i := 0; i < nn; i++ {
				pn := pAbs[i] * invP
				sqP := pn; if sqP > 1 { sqP = 1 }
				// Match computeFieldView: G = phi² * (1-sqrt(P_norm))
				gv := phiSq[i] * (1 - float32(math.Sqrt(float64(sqP)))) * invPhi
				bv := thetaSq[i] * invTh
				vol.rgba[i*4+0] = pn
				vol.rgba[i*4+1] = gv
				vol.rgba[i*4+2] = bv
				vol.rgba[i*4+3] = 1.0  // full alpha, same as computeFieldView
			}
			vol.loaded = true
		} else {
			// Standard mode: shared complex/gauge-aware dispatch
			computeStandardView(fl.sfa, fl.volN, &fl.workVol, mode, fl.fixedMax)
		}

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
	if fl.sfa.totalFrames == 0 {
		return // no frames available
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
	if fl.sfa.totalFrames == 0 {
		return
	}
	if params.curFrame >= int(fl.sfa.totalFrames) {
		params.curFrame = int(fl.sfa.totalFrames) - 1
	}
	if params.curFrame < 0 {
		params.curFrame = 0
	}
	// Sync fixed scale from viewParams
	if params.fixedScale {
		fl.fixedMax = params.fixedScaleVal
	} else {
		fl.fixedMax = 0
	}
	fl.startLoad(params.curFrame)
}

// ============================================================
// VecStream frame loader (async, similar to SFA frameLoader)
// ============================================================

type vsFrameLoader struct {
	vs            *vecstreamFile
	displayFrames []vsDisplayFrame
	volN          int
	fieldIdx      int        // which field to reconstruct
	mu            sync.Mutex
	ready         *volumeData
	readyT        float64
	readyF        int
	readyFType    uint8 // frame type of the loaded frame
	busy          atomic.Bool
	workVol       volumeData
}

func (vfl *vsFrameLoader) isLoading() bool {
	return vfl.busy.Load()
}

func (vfl *vsFrameLoader) startLoad(frameIdx int) {
	if vfl.busy.Load() {
		return
	}
	vfl.busy.Store(true)
	fi := vfl.fieldIdx
	fmt.Printf("Loading vecstream frame %d/%d (field %d)...\n", frameIdx+1, len(vfl.displayFrames), fi)
	go func() {
		defer vfl.busy.Store(false)

		if frameIdx < 0 || frameIdx >= len(vfl.displayFrames) {
			fmt.Printf("Frame %d out of range\n", frameIdx)
			return
		}
		df := vfl.displayFrames[frameIdx]

		// Find the index entry for the requested field at this display frame
		targetIdx := df.indexIdx // default: primary field
		for _, ei := range df.entryIdx {
			if int(vfl.vs.index[ei].fieldIdx) == fi {
				targetIdx = ei
				break
			}
		}

		t0 := time.Now()
		voxels, err := vfl.vs.reconstruct(targetIdx, fi)
		if err != nil {
			fmt.Printf("Reconstruct error: %v\n", err)
			return
		}
		dt := time.Since(t0)

		// Convert single-field voxels to RGBA volume for display
		// R = value, G = |value|, B = 0 (single field mode)
		n := vfl.volN
		nn := n * n * n
		vol := &vfl.workVol
		vol.n = n
		vol.loaded = true
		if len(vol.rgba) != nn*4 {
			vol.rgba = make([]float32, nn*4)
		}

		// Find max for normalization
		var maxVal float32
		for i := 0; i < nn; i++ {
			av := f32abs(voxels[i])
			if av > maxVal {
				maxVal = av
			}
		}
		invMax := float32(1)
		if maxVal*0.3 > 1e-20 {
			invMax = 1.0 / (maxVal * 0.3)
		}

		// Map to RGBA: R=positive, G=|value|, B=negative
		for i := 0; i < nn; i++ {
			v := voxels[i]
			av := f32abs(v) * invMax

			var rv, bv float32
			if v > 0 {
				rv = av
			} else {
				bv = av
			}
			gv := av * 0.3 // dim |value| for context

			j := i * 4
			vol.rgba[j+0] = rv
			vol.rgba[j+1] = gv
			vol.rgba[j+2] = bv
			vol.rgba[j+3] = 1.0
		}

		fType := vfl.vs.index[targetIdx].frameType
		fmt.Printf("  Reconstructed in %v (max=%.4g, type=%s)\n", dt, maxVal, vsFrameTypeName(fType))

		vfl.mu.Lock()
		vfl.ready = vol
		vfl.readyT = df.time
		vfl.readyF = frameIdx
		vfl.readyFType = fType
		vfl.mu.Unlock()
	}()
}

func (vfl *vsFrameLoader) consume() (*volumeData, float64, int, uint8, bool) {
	vfl.mu.Lock()
	defer vfl.mu.Unlock()
	if vfl.ready == nil {
		return nil, 0, 0, 0, false
	}
	v := vfl.ready
	t := vfl.readyT
	f := vfl.readyF
	ft := vfl.readyFType
	vfl.ready = nil
	return v, t, f, ft, true
}

func navigateVSFrame(vfl *vsFrameLoader, params *viewParams, delta int) {
	if vfl.isLoading() {
		return
	}
	nFrames := len(vfl.displayFrames)
	if nFrames == 0 {
		return
	}
	next := params.curFrame + delta
	if next < 0 {
		next = 0
	}
	if next >= nFrames {
		next = nFrames - 1
	}
	if next == params.curFrame {
		return
	}
	params.curFrame = next
	vfl.startLoad(next)
}

func loadVSFrame(vfl *vsFrameLoader, params *viewParams) {
	nFrames := len(vfl.displayFrames)
	if nFrames == 0 {
		return
	}
	if params.curFrame >= nFrames {
		params.curFrame = nFrames - 1
	}
	if params.curFrame < 0 {
		params.curFrame = 0
	}
	vfl.fieldIdx = params.vsFieldIdx
	vfl.startLoad(params.curFrame)
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

		// Extract columns and compute view (complex/gauge-aware)
		computeStandardView(sfa, int(sfa.Nx), vol, params.viewMode, 0)

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
	fmt.Println("  4              Field view  (R=|P|, G=phi^2, B=theta^2)")
	fmt.Println("  5              Velocity view (R=|v|, G=|vx|, B=|vy|)")
	fmt.Println("  6              Accel view  (R=cbrt|P|, G=theta^2, B=|v|)")
	fmt.Println("  8              U(1) gauge view (R=|E|, G=|Phi|^2, B=|A| links) [30-col files]")
	fmt.Println("  9              Charge view (R=+rhoQ, B=-rhoQ, G=|Phi|^2) [24/30-col files]")
	fmt.Println("  (complex files: field/velocity views use |Phi_a| moduli — phase-invariant)")
	fmt.Println("  M              Toggle PGA mode: spectrum / rotation (8-component v56 files)")
	fmt.Println("  N              Advance PGA rotation phase (when in rotation mode)")
	fmt.Println("  V              Toggle vec3d wireframe overlay")
	fmt.Println("  + / -          Brightness up/down")
	fmt.Println("  O / P          Opacity down/up")
	fmt.Println("  B              Toggle white/dark background")
	fmt.Println("  S              Save screenshot (WebP)")
	fmt.Println("  Shift+A        Export animated WebP (2x2 composite)")
	fmt.Println("  7              Cycle vecstream field (when .vecstream loaded)")
	fmt.Println("  R              Reset view")
	fmt.Println("  Escape / Q     Quit")
	fmt.Println()
}

func main() {
	fmt.Println("volview3: GPU volume ray marcher (OpenGL 4.3)")

	// Parse flags (before positional args)
	cpuprofile := flag.String("cpuprofile", "", "write CPU profile to file")
	snapshot := flag.String("snapshot", "", "Save snapshot(s) and exit. Modes: 'all' (all frames to dir), 'N' (single frame), 'N,M,...' (specific frames), 'animation' (animated webp)")
	outdir := flag.String("outdir", ".", "Output directory for snapshots")
	outfile := flag.String("out", "", "Output file path (for single frame or animation)")
	azimuthFlag := flag.Float64("azimuth", 0.8, "Camera azimuth in radians")
	elevationFlag := flag.Float64("elevation", 1.2, "Camera elevation in radians")
	distanceFlag := flag.Float64("distance", 2.8, "Camera distance")
	brightnessFlag := flag.Float64("brightness", 1.0, "Render brightness")
	opacityFlag := flag.Float64("opacity", 3.0, "Render opacity")
	widthFlag := flag.Int("width", 900, "Render width in pixels")
	heightFlag := flag.Int("height", 700, "Render height in pixels")
	bgWhiteFlag := flag.Bool("bg-white", false, "Use white background")
	compositeFlag := flag.Bool("composite", false, "Render 2x2 composite (4 views)")
	vecFlag := flag.String("vec", "", "Vec3D wireframe overlay file (.vec3d)")
	pgaModeFlag := flag.Int("pga-mode", 0, "PGA view mode: 0=spectrum, 1=rotation (8-comp v56 files)")
	pgaPhaseFlag := flag.Int("pga-phase", 0, "PGA rotation phase (which 3-of-8 components → R,G,B)")
	viewFlag := flag.Int("view", 0, "View mode: 0=field, 1=velocity, 2=accel, 3=U(1) gauge, 4=charge")
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

	// Parse positional args (SFA, VEC3D, or VecStream path) from remaining args after flags
	args := flag.Args()
	var sfaPath string
	var vecPath string
	var vecstreamPath string
	if *vecFlag != "" {
		vecPath = *vecFlag
	}
	for _, arg := range args {
		if strings.HasSuffix(arg, ".vec3d") {
			if vecPath == "" {
				vecPath = arg
			}
		} else if strings.HasSuffix(arg, ".vecstream") {
			if vecstreamPath == "" {
				vecstreamPath = arg
			}
		} else if sfaPath == "" {
			sfaPath = arg
		}
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

	// Load VecStream if provided
	var vstream *vecstreamFile
	var vsDisplayFrames []vsDisplayFrame
	if vecstreamPath != "" {
		var err error
		vstream, err = vecstreamOpen(vecstreamPath)
		if err != nil {
			log.Fatalf("Failed to open VecStream: %v", err)
		}
		defer vstream.close()
		vsDisplayFrames = vstream.buildDisplayFrames()
		fmt.Printf("VecStream: %d display frames, %d fields\n", len(vsDisplayFrames), vstream.nFields)
	}

	// Load Vec3D if provided
	var vec *vec3dFile
	if vecPath != "" {
		var err error
		vec, err = vec3dOpen(vecPath)
		if err != nil {
			log.Fatalf("Failed to open VEC3D: %v", err)
		}
	}

	if sfa == nil && vec == nil && vstream == nil {
		fmt.Fprintf(os.Stderr, "Usage: volview3 [-cpuprofile file] [-vec file.vec3d] <file.sfa|file.vec3d|file.vecstream>\n")
		os.Exit(1)
	}

	var volN int
	if sfa != nil {
		volN = int(sfa.Nx)
	} else if vstream != nil {
		volN = int(vstream.Nx)
	} else if vec != nil {
		volN = int(vec.Nx)
	}
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
	if *snapshot != "" {
		glfw.WindowHint(glfw.Visible, glfw.False)
	}

	winW, winH := defaultWinW, defaultWinH
	if *snapshot != "" {
		winW, winH = *widthFlag, *heightFlag
	}
	window, err := glfw.CreateWindow(winW, winH, "volview3 -- GPU Volume Viewer", nil, nil)
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

	// Initialize line GPU resources for vec3d wireframe
	var lineGPU *lineGPUResources
	if vec != nil {
		var err error
		lineGPU, err = initLineGPUResources()
		if err != nil {
			log.Fatalf("Init line GPU: %v", err)
		}
		wireData, nVerts := vec.buildWireframeData()
		lineGPU.uploadWireframe(wireData, nVerts)
		fmt.Printf("Wireframe: %d vertices uploaded\n", nVerts)
	}

	fmt.Printf("Volume: %dx%dx%d loaded\n", volN, volN, volN)
	printHelp()

	params := newViewParams()

	// Apply CLI camera overrides
	params.azimuth = float32(*azimuthFlag)
	params.elevation = float32(*elevationFlag)
	params.distance = float32(*distanceFlag)
	params.brightness = float32(*brightnessFlag)
	params.pgaMode = *pgaModeFlag
	params.pgaPhase = *pgaPhaseFlag
	params.viewMode = *viewFlag
	params.opacity = float32(*opacityFlag)
	params.bgWhite = *bgWhiteFlag
	params.updateCamera()

	var fl *frameLoader // declared early for callbacks; initialized after snapshot block

	// === Headless snapshot mode (runs before async loader) ===
	if *snapshot != "" && sfa != nil {
		renderW := *widthFlag
		renderH := *heightFlag

		os.MkdirAll(*outdir, 0755)

		// Parse frame list
		var frameList []int
		if *snapshot == "all" {
			for i := 0; i < int(sfa.totalFrames); i++ {
				frameList = append(frameList, i)
			}
		} else if *snapshot == "animation" {
			for i := 0; i < int(sfa.totalFrames); i++ {
				frameList = append(frameList, i)
			}
		} else {
			// Parse comma-separated frame indices
			for _, s := range strings.Split(*snapshot, ",") {
				s = strings.TrimSpace(s)
				if s == "" {
					continue
				}
				var n int
				fmt.Sscanf(s, "%d", &n)
				if n >= 0 && n < int(sfa.totalFrames) {
					frameList = append(frameList, n)
				}
			}
		}

		if len(frameList) == 0 {
			fmt.Fprintf(os.Stderr, "No valid frames to render\n")
			os.Exit(1)
		}

		// Create offscreen FBO
		var fbo, fboTex uint32
		fboW, fboH := renderW, renderH
		if *compositeFlag {
			fboW, fboH = renderW*2, renderH*2
		}
		gl.GenFramebuffers(1, &fbo)
		gl.BindFramebuffer(gl.FRAMEBUFFER, fbo)
		gl.GenTextures(1, &fboTex)
		gl.BindTexture(gl.TEXTURE_2D, fboTex)
		gl.TexImage2D(gl.TEXTURE_2D, 0, gl.RGBA8, int32(fboW), int32(fboH), 0, gl.RGBA, gl.UNSIGNED_BYTE, nil)
		gl.TexParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR)
		gl.TexParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR)
		gl.FramebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, fboTex, 0)

		if gl.CheckFramebufferStatus(gl.FRAMEBUFFER) != gl.FRAMEBUFFER_COMPLETE {
			log.Fatalf("FBO not complete")
		}

		gl.ActiveTexture(gl.TEXTURE0)
		gl.BindTexture(gl.TEXTURE_3D, gpu.volTexture)

		var animImages []image.Image
		var animDurations []uint
		var animDisposals []uint

		for _, fi := range frameList {
			frameTime, err := sfa.readFrame(uint32(fi))
			if err != nil {
				fmt.Fprintf(os.Stderr, "Frame %d error: %v\n", fi, err)
				continue
			}

			// Auto-detect v56 PGA file (8 columns, first named "M0_…").
			snapColName0 := strings.TrimRight(string(sfa.columns[0].name[:]), "\x00")
			isPGASnap := sfa.nColumns == 8 && strings.HasPrefix(snapColName0, "M0")
			if isPGASnap {
				var pgaCols [8][]float32
				for k := 0; k < 8; k++ {
					pgaCols[k] = sfa.extractColumnF32(k)
				}
				computePGAView(int(sfa.Nx), pgaCols, vol, params.pgaMode, params.pgaPhase)
			} else {
				computeStandardView(sfa, int(sfa.Nx), vol, params.viewMode, 0)
			}
			gpu.uploadVolume(vol)

			gl.BindFramebuffer(gl.FRAMEBUFFER, fbo)
			gl.Clear(gl.COLOR_BUFFER_BIT)

			if *compositeFlag {
				halfW, halfH := renderW, renderH
				renderQuadrant(gpu, params, 0, halfH, halfW, halfH, 1, 1, 1, false)
				renderQuadrant(gpu, params, halfW, halfH, halfW, halfH, 1, 0, 1, false)
				renderQuadrant(gpu, params, 0, 0, halfW, halfH, 1, 1, 1, true)
				renderQuadrant(gpu, params, halfW, 0, halfW, halfH, 0, 1, 0, false)
			} else {
				mr, mg, mb := float32(1), float32(1), float32(1)
				renderQuadrant(gpu, params, 0, 0, fboW, fboH, mr, mg, mb, params.bgWhite)
			}

			img := readFramebuffer(0, 0, fboW, fboH)

			if *snapshot == "animation" {
				animImages = append(animImages, img)
				animDurations = append(animDurations, 200)
				animDisposals = append(animDisposals, 1)
				fmt.Printf("  Frame %d/%d t=%.3f\n", fi+1, int(sfa.totalFrames), frameTime)
			} else {
				// Save individual image
				outPath := *outfile
				if outPath == "" {
					outPath = fmt.Sprintf("%s/frame_%04d_t%.3f.webp", *outdir, fi, frameTime)
				}
				f, err := os.Create(outPath)
				if err != nil {
					fmt.Fprintf(os.Stderr, "Create %s: %v\n", outPath, err)
					continue
				}
				if err := nativewebp.Encode(f, img, nil); err != nil {
					fmt.Fprintf(os.Stderr, "Encode %s: %v\n", outPath, err)
				}
				f.Close()
				fmt.Printf("Saved %s (frame %d, t=%.3f, %dx%d)\n", outPath, fi, frameTime, fboW, fboH)
			}
		}

		// Save animation if in animation mode
		if *snapshot == "animation" && len(animImages) > 0 {
			outPath := *outfile
			if outPath == "" {
				outPath = fmt.Sprintf("%s/animation.webp", *outdir)
			}
			f, err := os.Create(outPath)
			if err != nil {
				log.Fatalf("Create animation: %v", err)
			}
			ani := &nativewebp.Animation{
				Images:          animImages,
				Durations:       animDurations,
				Disposals:       animDisposals,
				LoopCount:       0,
				BackgroundColor: 0xFF0D0D12,
			}
			if err := nativewebp.EncodeAll(f, ani, nil); err != nil {
				log.Fatalf("Encode animation: %v", err)
			}
			f.Close()
			fmt.Printf("Animation saved: %s (%dx%d, %d frames)\n", outPath, fboW, fboH, len(animImages))
		}

		gl.BindFramebuffer(gl.FRAMEBUFFER, 0)
		gl.DeleteTextures(1, &fboTex)
		gl.DeleteFramebuffers(1, &fbo)
		fmt.Println("Snapshot mode complete")
		os.Exit(0)
	}

	// Start async frame loader (only for interactive mode — snapshot exits above)
	var vfl *vsFrameLoader // vecstream frame loader
	if sfa != nil && fl == nil {
		total := volN * volN * volN
		fl = &frameLoader{sfa: sfa, volN: volN, viewMode: *viewFlag}
		fl.workVol.absPBuf = make([]float32, total)
		fl.workVol.phi2Buf = make([]float32, total)
		fl.workVol.theta2Buf = make([]float32, total)
		fl.workVol.rgba = make([]float32, total*4)
		loadSFAFrame(fl, params)
	} else if vstream != nil {
		total := volN * volN * volN
		vfl = &vsFrameLoader{
			vs:            vstream,
			displayFrames: vsDisplayFrames,
			volN:          volN,
		}
		vfl.workVol.rgba = make([]float32, total*4)
		loadVSFrame(vfl, params)
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
		case glfw.KeyV:
			params.showVec = !params.showVec
			if params.showVec {
				fmt.Println("Vec3D overlay: ON")
			} else {
				fmt.Println("Vec3D overlay: OFF")
			}
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
		case glfw.Key4:
			if fl != nil && params.viewMode != 0 {
				params.viewMode = 0
				fl.viewMode = 0
				fmt.Println("View: FIELD (R=|P|, G=phi^2, B=theta^2)")
				loadSFAFrame(fl, params)
			}
		case glfw.Key5:
			if fl != nil && params.viewMode != 1 {
				params.viewMode = 1
				fl.viewMode = 1
				fmt.Println("View: VELOCITY (R=|v|, G=|vx|, B=|vy|)")
				loadSFAFrame(fl, params)
			}
		case glfw.Key6:
			if fl != nil && params.viewMode != 2 {
				params.viewMode = 2
				fl.viewMode = 2
				fmt.Println("View: ACCEL (R=cbrt|P|, G=theta^2, B=|v|)")
				loadSFAFrame(fl, params)
			}
		case glfw.Key8:
			if fl != nil && params.viewMode != 3 {
				params.viewMode = 3
				fl.viewMode = 3
				fmt.Println("View: U(1) GAUGE (R=|E|, G=|Phi|^2, B=|A| links)")
				loadSFAFrame(fl, params)
			}
		case glfw.Key9:
			if fl != nil && params.viewMode != 4 {
				params.viewMode = 4
				fl.viewMode = 4
				fmt.Println("View: CHARGE (R=+rhoQ, B=-rhoQ, G=|Phi|^2 context)")
				loadSFAFrame(fl, params)
			}
		case glfw.KeyM:
			// Toggle PGA mode (only meaningful for 8-component PGA files).
			if fl != nil {
				params.pgaMode = (params.pgaMode + 1) % 2
				fl.pgaMode = params.pgaMode
				if params.pgaMode == 0 {
					fmt.Println("PGA view: SPECTRUM (all 8 components → RGB via hue)")
				} else {
					fmt.Printf("PGA view: ROTATION (phase=%d → M[%d,%d,%d] → R,G,B)\n",
						params.pgaPhase, params.pgaPhase%8,
						(params.pgaPhase+1)%8, (params.pgaPhase+2)%8)
				}
				loadSFAFrame(fl, params)
			}
		case glfw.KeyN:
			// Advance PGA rotation phase.
			if fl != nil && params.pgaMode == 1 {
				params.pgaPhase = (params.pgaPhase + 1) % 8
				fl.pgaPhase = params.pgaPhase
				fmt.Printf("PGA rotation: phase=%d → M[%d,%d,%d] → R,G,B\n",
					params.pgaPhase, params.pgaPhase%8,
					(params.pgaPhase+1)%8, (params.pgaPhase+2)%8)
				loadSFAFrame(fl, params)
			}
		case glfw.KeyF:
			params.fixedScale = !params.fixedScale
			if params.fixedScale {
				fmt.Printf("Scale: FIXED (max=%.4f) — use [ ] to adjust\n", params.fixedScaleVal)
			} else {
				fmt.Println("Scale: AUTO")
			}
			if fl != nil {
				fl.fixedMax = 0
				if params.fixedScale {
					fl.fixedMax = params.fixedScaleVal
				}
				loadSFAFrame(fl, params)
			}
		case glfw.KeyLeftBracket: // [ key — decrease fixed scale (zoom in)
			params.fixedScaleVal *= 0.5
			if params.fixedScaleVal < 1e-8 {
				params.fixedScaleVal = 1e-8
			}
			if params.fixedScale {
				fmt.Printf("Scale: fixed max=%.6f\n", params.fixedScaleVal)
				if fl != nil {
					fl.fixedMax = params.fixedScaleVal
					loadSFAFrame(fl, params)
				}
			}
		case glfw.KeyRightBracket: // ] key — increase fixed scale (zoom out)
			params.fixedScaleVal *= 2.0
			if params.fixedScale {
				fmt.Printf("Scale: fixed max=%.6f\n", params.fixedScaleVal)
				if fl != nil {
					fl.fixedMax = params.fixedScaleVal
					loadSFAFrame(fl, params)
				}
			}
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
			} else if vfl != nil {
				navigateVSFrame(vfl, params, -1)
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
				} else if vfl != nil {
					navigateVSFrame(vfl, params, -1)
				}
			}
		case glfw.KeyRight, glfw.KeyD:
			if fl != nil {
				navigateFrame(fl, params, 1)
			} else if vfl != nil {
				navigateVSFrame(vfl, params, 1)
			}
		case glfw.KeyHome:
			if fl != nil {
				params.curFrame = 0
				loadSFAFrame(fl, params)
			} else if vfl != nil {
				params.curFrame = 0
				loadVSFrame(vfl, params)
			}
		case glfw.KeyEnd:
			if fl != nil {
				params.curFrame = int(sfa.totalFrames) - 1
				loadSFAFrame(fl, params)
			} else if vfl != nil {
				params.curFrame = len(vfl.displayFrames) - 1
				loadVSFrame(vfl, params)
			}
		// Vecstream field cycling
		case glfw.Key7:
			if vfl != nil {
				params.vsFieldIdx = (params.vsFieldIdx + 1) % int(vstream.nFields)
				fmt.Printf("VecStream field: %d/%d\n", params.vsFieldIdx, vstream.nFields)
				// Reset reconstruction cache for new field
				loadVSFrame(vfl, params)
			}
		}
	})

	// -- Main loop --
	for !window.ShouldClose() {
		glfw.PollEvents()

		// Check async frame loader (SFA)
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

		// Check async frame loader (VecStream)
		if vfl != nil {
			if newVol, frameTime, frameIdx, fType, ok := vfl.consume(); ok {
				if newVol.n != gpu.volSize {
					gpu.createVolumeTexture(newVol.n)
				}
				gpu.uploadVolume(newVol)
				fmt.Printf("Frame %d/%d t=%.3f [%s] field=%d (ready)\n",
					frameIdx+1, len(vfl.displayFrames), frameTime,
					vsFrameTypeName(fType), vfl.fieldIdx)
			}
		}

		// Animation (SFA)
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

		// Animation (VecStream)
		if params.playing && vfl != nil && !vfl.isLoading() {
			now := time.Now()
			if now.Sub(lastAnimTime) > animInterval {
				nextFrame := params.curFrame + 1
				if nextFrame >= len(vfl.displayFrames) {
					nextFrame = 0
				}
				params.curFrame = nextFrame
				loadVSFrame(vfl, params)
				lastAnimTime = now
			}
		}

		// Render
		w, h := window.GetFramebufferSize()
		gl.Viewport(0, 0, int32(w), int32(h))
		gl.Clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)

		if sfa != nil || vstream != nil {
			gpu.updateUniforms(params, float32(w), float32(h))
			gpu.render()
		} else if !params.bgWhite {
			// Vec3d-only mode: clear to dark background
			gl.ClearColor(0.05, 0.05, 0.07, 1.0)
			gl.Clear(gl.COLOR_BUFFER_BIT)
		} else {
			gl.ClearColor(1.0, 1.0, 1.0, 1.0)
			gl.Clear(gl.COLOR_BUFFER_BIT)
		}

		// Wireframe overlay
		if lineGPU != nil && params.showVec {
			lineGPU.render(params, float32(w), float32(h))
		}

		window.SwapBuffers()
	}

	fmt.Println("Exit")
}
