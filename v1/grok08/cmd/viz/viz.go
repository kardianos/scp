package main

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"sort"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/ebitenutil"
	"github.com/hajimehoshi/ebiten/v2/inpututil"
)

type Dimensions struct {
	Height int
	Width  int
	Length int
}

type Game struct {
	fields     map[int][]float64 // timestep -> flattened field
	dims       Dimensions
	maxT       int
	sliceDim   int   // 0=X, 1=Y, 2=Z, 3=T
	slicePos   []int // Positions for [X,Y,Z,T]
	holdLeft   bool
	holdRight  bool
	shift      bool
	sliceImage *ebiten.Image
}

func NewGame(filePath string) *Game {
	dims, fields := loadFile(filePath)
	keys := make([]int, 0, len(fields))
	for k := range fields {
		keys = append(keys, k)
	}
	sort.Ints(keys)
	maxT := keys[len(keys)-1]

	g := &Game{
		fields:   fields,
		dims:     dims,
		maxT:     maxT,
		sliceDim: 3, // Start with T as variable
		slicePos: []int{dims.Width / 2, dims.Height / 2, dims.Length / 2, maxT / 2},
	}
	return g
}

func loadFile(filePath string) (Dimensions, map[int][]float64) {
	data, err := os.ReadFile(filePath)
	if err != nil {
		log.Fatal(err)
	}
	reader := bytes.NewReader(data)

	magic := make([]byte, 8)
	reader.Read(magic)
	if string(magic) != "SCPV\x00\x00\x00\x00" {
		log.Fatal("Invalid magic")
	}

	dims := Dimensions{}
	fields := make(map[int][]float64)

	for {
		var typeID int64
		err = binary.Read(reader, binary.BigEndian, &typeID)
		if err == io.EOF {
			break
		} else if err != nil {
			log.Fatal(err)
		}

		var nameLen int16
		binary.Read(reader, binary.BigEndian, &nameLen)
		nameBytes := make([]byte, nameLen)
		reader.Read(nameBytes)
		name := string(nameBytes)

		var valLen int64
		binary.Read(reader, binary.BigEndian, &valLen)
		valBytes := make([]byte, valLen)
		reader.Read(valBytes)

		switch typeID {
		case 0: // height
			dims.Height = int(binary.BigEndian.Uint32(valBytes))
		case 1: // width
			dims.Width = int(binary.BigEndian.Uint32(valBytes))
		case 2: // length
			dims.Length = int(binary.BigEndian.Uint32(valBytes))
		case 4: // field data
			if name[0] == 't' {
				t := 0
				fmt.Sscanf(name, "t%d", &t)
				field := make([]float64, len(valBytes)/8)
				for i := 0; i < len(field); i++ {
					field[i] = math.Float64frombits(binary.LittleEndian.Uint64(valBytes[i*8 : (i+1)*8]))
				}
				fields[t] = field
			}
		}
	}
	return dims, fields
}

func (g *Game) Update() error {
	// Handle key presses
	if inpututil.IsKeyJustPressed(ebiten.Key1) {
		g.sliceDim = 0
	}
	if inpututil.IsKeyJustPressed(ebiten.Key2) {
		g.sliceDim = 1
	}
	if inpututil.IsKeyJustPressed(ebiten.Key3) {
		g.sliceDim = 2
	}
	if inpututil.IsKeyJustPressed(ebiten.Key4) {
		g.sliceDim = 3
	}
	g.shift = ebiten.IsKeyPressed(ebiten.KeyShiftLeft) || ebiten.IsKeyPressed(ebiten.KeyShiftRight)

	if inpututil.IsKeyJustPressed(ebiten.KeyLeft) || g.holdLeft {
		g.holdLeft = true
		g.slicePos[g.sliceDim]--
		if !g.shift {
			g.holdLeft = false
		}
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyRight) || g.holdRight {
		g.holdRight = true
		g.slicePos[g.sliceDim]++
		if !g.shift {
			g.holdRight = false
		}
	}
	if inpututil.IsKeyJustReleased(ebiten.KeyLeft) {
		g.holdLeft = false
	}
	if inpututil.IsKeyJustReleased(ebiten.KeyRight) {
		g.holdRight = false
	}

	// Clamp positions
	maxPos := []int{g.dims.Width - 1, g.dims.Height - 1, g.dims.Length - 1, g.maxT}
	g.slicePos[g.sliceDim] = max(0, min(g.slicePos[g.sliceDim], maxPos[g.sliceDim]))

	return nil
}

func (g *Game) Draw(screen *ebiten.Image) {
	// Render slice
	viewW, viewH := getSliceDimensions(g)
	if g.sliceImage == nil || g.sliceImage.Bounds().Dx() != viewW || g.sliceImage.Bounds().Dy() != viewH {
		g.sliceImage = ebiten.NewImage(viewW, viewH)
	}

	field := g.getSlice()

	pixels := make([]byte, viewW*viewH*4)
	minV, maxV := field[0], field[0]
	for _, v := range field {
		if v < minV {
			minV = v
		}
		if v > maxV {
			maxV = v
		}
	}
	rangeV := maxV - minV
	if rangeV == 0 {
		rangeV = 1
	}
	for j := 0; j < viewH; j++ {
		for i := 0; i < viewW; i++ {
			val := (field[j*viewW+i] - minV) / rangeV
			idx := (j*viewW + i) * 4
			pixels[idx] = uint8(255 * val)         // R
			pixels[idx+1] = 0                      // G
			pixels[idx+2] = uint8(255 * (1 - val)) // B
			pixels[idx+3] = 255                    // A
		}
	}
	g.sliceImage.WritePixels(pixels)

	op := &ebiten.DrawImageOptions{}
	scale := math.Min(float64(screen.Bounds().Dx())/float64(viewW), float64(screen.Bounds().Dy())/float64(viewH))
	op.GeoM.Scale(scale, scale)
	screen.DrawImage(g.sliceImage, op)

	// Overlay text
	dimNames := []string{"X", "Y", "Z", "T"}
	screenXY := getScreenXY(g.sliceDim)
	fixedAxes := getFixedAxes(g.sliceDim, g.slicePos, dimNames)
	varAxis := fmt.Sprintf("VAR %s=%d", dimNames[g.sliceDim], g.slicePos[g.sliceDim])
	msg := fmt.Sprintf("Screen XY: %s | %s | %s", screenXY, fixedAxes, varAxis)
	ebitenutil.DebugPrint(screen, msg)
}

func (g *Game) Layout(outsideWidth, outsideHeight int) (int, int) {
	return 800, 800
}

func getSliceDimensions(g *Game) (w, h int) {
	switch g.sliceDim {
	case 0:
		return g.dims.Height, g.dims.Length // YZ
	case 1:
		return g.dims.Width, g.dims.Length // XZ
	case 2:
		return g.dims.Width, g.dims.Height // XY
	case 3:
		return g.dims.Width, g.dims.Height // XY at T
	}
	return g.dims.Width, g.dims.Height
}

func (g *Game) getSlice() []float64 {
	field := g.fields[g.slicePos[3]] // Fixed T for dim<3, but var for dim=3
	w, h := getSliceDimensions(g)
	slice := make([]float64, w*h)
	for j := 0; j < h; j++ {
		for i := 0; i < w; i++ {
			var x, y, z int
			switch g.sliceDim {
			case 0: // YZ slice at fixed X
				x, y, z = g.slicePos[0], i, j
			case 1: // XZ slice at fixed Y
				x, y, z = i, g.slicePos[1], j
			case 2: // XY slice at fixed Z
				x, y, z = i, j, g.slicePos[2]
			case 3: // XY slice at fixed T (but T is var, wait no: for T var, it's XY at Z fixed, but user intent XY varying T? Wait, adjust
				x, y, z = i, j, g.slicePos[2] // Assume fixed Z for T slice
			}
			idx := x + y*g.dims.Width + z*g.dims.Width*g.dims.Height
			slice[j*w+i] = field[idx]
		}
	}
	return slice
}

func getScreenXY(dim int) string {
	// dimNames := []string{"X", "Y", "Z", "T"}
	switch dim {
	case 0:
		return "YZ"
	case 1:
		return "XZ"
	case 2:
		return "XY"
	case 3:
		return "XY" // For T var
	}
	return "XY"
}

func getFixedAxes(dim int, pos []int, dimNames []string) string {
	fixed := ""
	for i := 0; i < 4; i++ {
		if i != dim && i != 3 { // Exclude T if not fixed
			fixed += fmt.Sprintf("@%s=%d, ", dimNames[i], pos[i])
		}
	}
	return fixed[:len(fixed)-2] // Trim comma
}

func main() {
	game := NewGame("field.scpv")
	ebiten.SetWindowSize(800, 800)
	ebiten.SetWindowTitle("Field Viz - Slice XYZT")
	if err := ebiten.RunGame(game); err != nil {
		log.Fatal(err)
	}
}
