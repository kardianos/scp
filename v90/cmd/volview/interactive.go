// interactive.go — the interactive viewer (-i): rotate the simulation,
// scrub time, switch field aspects. Basis: sfa/volview (GL + GLFW,
// emission compositing); here as point-sprite splats over .fcs frames.
//
// Controls (printed on start and with H):
//   drag             orbit          wheel        zoom
//   shift+drag       pan            ,  /  .      step frame back/fwd
//   space            play/pause     [  /  ]      play speed
//   1..8             channel: es em ee x r phase fa1 thd
//   L                links on/off   T            tag tint on/off
//   A                time-average   N            per-frame / global norm
//   + / -            point scale    G / shift+G  gamma
//   S                screenshot     H            help    ESC/Q  quit
package main

import (
	"fmt"
	"image"
	"image/png"
	"math"
	"os"
	"runtime"
	"sort"
	"strings"
	"time"
	"unsafe"

	"github.com/go-gl/gl/v4.1-core/gl"
	"github.com/go-gl/glfw/v3.3/glfw"
)

func init() { runtime.LockOSThread() }

const vertSrc = `#version 410 core
layout (location = 0) in vec3 pos;
layout (location = 1) in float val;   // channel value in [0,1] (or -1 = dark)
layout (location = 2) in float tag;
uniform mat4 mvp;
uniform float ptScale;   // pixels per world unit at distance 1
layout (location = 3) in float radius;
out float vVal;
out float vTag;
void main() {
    gl_Position = mvp * vec4(pos, 1.0);
    float w = max(gl_Position.w, 1e-3);
    gl_PointSize = clamp(ptScale * radius / w, 1.0, 220.0);
    vVal = val;
    vTag = tag;
}`

const fragSrc = `#version 410 core
in float vVal;
in float vTag;
uniform float gamma;
uniform int cyclic;
uniform int tagTint;
out vec4 frag;
vec3 ramp(float v) { // heat: black -> orange -> white
    return vec3(min(1.0, 2.2*v),
                clamp(1.6*v - 0.25, 0.0, 1.0),
                clamp(2.5*v - 1.5, 0.0, 1.0));
}
vec3 wheel(float v) { // cyclic phase wheel (hue)
    float h = fract(v) * 6.0;
    vec3 c = clamp(abs(mod(vec3(h) + vec3(0.0, 4.0, 2.0), 6.0) - 3.0) - 1.0, 0.0, 1.0);
    return c * 0.9 + 0.05;
}
void main() {
    vec2 d = gl_PointCoord * 2.0 - 1.0;
    float r2 = dot(d, d);
    if (r2 > 1.0) discard;
    if (vVal < 0.0) discard;              // sentinel: dark cell
    float fall = 1.0 - r2;                 // soft sphere profile
    float v = pow(clamp(vVal, 0.0, 1.0), gamma);
    vec3 col = cyclic == 1 ? wheel(vVal) * (0.25 + 0.75*v) : ramp(v);
    if (tagTint == 1 && vTag > 0.5)
        col = vec3(0.2*col.r, 0.5 + 0.5*col.g, 0.5 + 0.5*col.b);
    float a = (0.15 + 0.7*v) * fall;
    frag = vec4(col * a, a);               // premultiplied
}`

const lineVertSrc = `#version 410 core
layout (location = 0) in vec3 pos;
layout (location = 1) in float w;
uniform mat4 mvp;
out float vW;
void main() { gl_Position = mvp * vec4(pos, 1.0); vW = w; }`

const lineFragSrc = `#version 410 core
in float vW;
out vec4 frag;
void main() {
    float a = clamp(vW, 0.0, 1.0) * 0.55;
    frag = vec4(vec3(0.45, 0.75, 1.0) * a, a);
}`

func compile(kind uint32, src string) uint32 {
	sh := gl.CreateShader(kind)
	cs, free := gl.Strs(src + "\x00")
	gl.ShaderSource(sh, 1, cs, nil)
	free()
	gl.CompileShader(sh)
	var ok int32
	gl.GetShaderiv(sh, gl.COMPILE_STATUS, &ok)
	if ok == gl.FALSE {
		var n int32
		gl.GetShaderiv(sh, gl.INFO_LOG_LENGTH, &n)
		log := strings.Repeat("\x00", int(n)+1)
		gl.GetShaderInfoLog(sh, n, nil, gl.Str(log))
		fmt.Fprintln(os.Stderr, "shader:", log)
		os.Exit(1)
	}
	return sh
}

func program(vs, fs string) uint32 {
	p := gl.CreateProgram()
	v := compile(gl.VERTEX_SHADER, vs)
	f := compile(gl.FRAGMENT_SHADER, fs)
	gl.AttachShader(p, v)
	gl.AttachShader(p, f)
	gl.LinkProgram(p)
	var ok int32
	gl.GetProgramiv(p, gl.LINK_STATUS, &ok)
	if ok == gl.FALSE {
		var n int32
		gl.GetProgramiv(p, gl.INFO_LOG_LENGTH, &n)
		log := strings.Repeat("\x00", int(n)+1)
		gl.GetProgramInfoLog(p, n, nil, gl.Str(log))
		fmt.Fprintln(os.Stderr, "link:", log)
		os.Exit(1)
	}
	gl.DeleteShader(v)
	gl.DeleteShader(f)
	return p
}

// mat4 helpers (column-major)
type mat4 [16]float32

func perspective(fovy, aspect, near, far float32) mat4 {
	f := float32(1.0 / math.Tan(float64(fovy)/2))
	var m mat4
	m[0] = f / aspect
	m[5] = f
	m[10] = (far + near) / (near - far)
	m[11] = -1
	m[14] = 2 * far * near / (near - far)
	return m
}

func mul(a, b mat4) mat4 {
	var m mat4
	for c := 0; c < 4; c++ {
		for r := 0; r < 4; r++ {
			var s float32
			for k := 0; k < 4; k++ {
				s += a[k*4+r] * b[c*4+k]
			}
			m[c*4+r] = s
		}
	}
	return m
}

func lookAt(eye, ctr [3]float32) mat4 {
	fw := [3]float32{ctr[0] - eye[0], ctr[1] - eye[1], ctr[2] - eye[2]}
	fl := float32(math.Sqrt(float64(fw[0]*fw[0] + fw[1]*fw[1] + fw[2]*fw[2])))
	fw = [3]float32{fw[0] / fl, fw[1] / fl, fw[2] / fl}
	up := [3]float32{0, 0, 1}
	s := [3]float32{fw[1]*up[2] - fw[2]*up[1], fw[2]*up[0] - fw[0]*up[2], fw[0]*up[1] - fw[1]*up[0]}
	sl := float32(math.Sqrt(float64(s[0]*s[0] + s[1]*s[1] + s[2]*s[2])))
	if sl < 1e-6 {
		s = [3]float32{1, 0, 0}
		sl = 1
	}
	s = [3]float32{s[0] / sl, s[1] / sl, s[2] / sl}
	u := [3]float32{s[1]*fw[2] - s[2]*fw[1], s[2]*fw[0] - s[0]*fw[2], s[0]*fw[1] - s[1]*fw[0]}
	var m mat4
	m[0], m[4], m[8] = s[0], s[1], s[2]
	m[1], m[5], m[9] = u[0], u[1], u[2]
	m[2], m[6], m[10] = -fw[0], -fw[1], -fw[2]
	m[15] = 1
	m[12] = -(s[0]*eye[0] + s[1]*eye[1] + s[2]*eye[2])
	m[13] = -(u[0]*eye[0] + u[1]*eye[1] + u[2]*eye[2])
	m[14] = fw[0]*eye[0] + fw[1]*eye[1] + fw[2]*eye[2]
	return m
}

type viewer struct {
	frames   []*frame
	avg      *frame
	channels []string
	chanIdx  int
	frameIdx int
	playing bool
	speed   float64 // frames per second
	useAvg  bool
	globalN bool
	links   bool
	tagTint bool
	gamma   float32
	ptScale float32
	// camera
	yaw, pitch float32
	dist       float32
	panX, panY float32 // in world units, camera plane
	// interaction
	lastX, lastY float64
	dragging     bool
	panning      bool
	// norm cache: [channel][frame] -> vmax ; global: [channel]
	vmaxFrame map[int][]float64
	vmaxGlob  map[int]float64
	dirty     bool
	shot      int
	selfcheck bool
}

func p99max(vals []float64) float64 {
	s := append([]float64(nil), vals...)
	sort.Float64s(s)
	v := 1e-12
	if len(s) > 0 {
		v = math.Max(v, s[len(s)*99/100])
	}
	return v
}

func (vw *viewer) frameFor() *frame {
	if vw.useAvg && vw.avg != nil {
		return vw.avg
	}
	return vw.frames[vw.frameIdx]
}

func (vw *viewer) channelVals(f *frame) ([]float32, float64) {
	name := vw.channels[vw.chanIdx]
	raw := make([]float64, len(f.cells))
	for i := range f.cells {
		raw[i] = f.chanVal(i, name)
	}
	var vmax float64 = 1
	if channelCyclic(name) {
		// cyclic channels are already in [0,1) (or -1 sentinel)
	} else if name == "fa1" {
		// diverging: map [-m, m] -> [0,1]
		m := 1e-12
		for _, v := range raw {
			if a := math.Abs(v); a > m {
				m = a
			}
		}
		for i, v := range raw {
			raw[i] = 0.5 + 0.5*v/m
		}
	} else {
		if vw.globalN {
			if v, ok := vw.vmaxGlob[vw.chanIdx]; ok {
				vmax = v
			} else {
				var all []float64
				for _, fr := range vw.frames {
					for i := range fr.cells {
						all = append(all, fr.chanVal(i, name))
					}
				}
				vmax = p99max(all)
				vw.vmaxGlob[vw.chanIdx] = vmax
			}
		} else {
			vmax = p99max(raw)
		}
		for i, v := range raw {
			raw[i] = v / vmax
		}
	}
	out := make([]float32, len(raw))
	for i, v := range raw {
		out[i] = float32(v)
	}
	return out, vmax
}

func runInteractive(frames []*frame, selfcheck bool) {
	vw := &viewer{
		frames: frames, avg: buildAvgFrame(frames),
		channels: runChannels(frames[0]),
		speed: 8, gamma: 0.7, ptScale: 900,
		yaw: 0.6, pitch: 0.5, dist: float32(frames[0].L) * 1.6,
		tagTint: true, links: len(frames[0].links) > 0,
		vmaxFrame: map[int][]float64{}, vmaxGlob: map[int]float64{},
		dirty: true, selfcheck: selfcheck,
	}

	if err := glfw.Init(); err != nil {
		fmt.Fprintln(os.Stderr, "glfw:", err)
		os.Exit(1)
	}
	defer glfw.Terminate()
	glfw.WindowHint(glfw.ContextVersionMajor, 4)
	glfw.WindowHint(glfw.ContextVersionMinor, 1)
	glfw.WindowHint(glfw.OpenGLProfile, glfw.OpenGLCoreProfile)
	glfw.WindowHint(glfw.OpenGLForwardCompatible, glfw.True)
	glfw.WindowHint(glfw.Samples, 4)
	win, err := glfw.CreateWindow(1100, 900, "fcs volview", nil, nil)
	if err != nil {
		fmt.Fprintln(os.Stderr, "window:", err)
		os.Exit(1)
	}
	win.MakeContextCurrent()
	if err := gl.Init(); err != nil {
		fmt.Fprintln(os.Stderr, "gl:", err)
		os.Exit(1)
	}
	glfw.SwapInterval(1)

	prog := program(vertSrc, fragSrc)
	lprog := program(lineVertSrc, lineFragSrc)

	var vao, vboPos, vboVal, vboTag, vboRad uint32
	gl.GenVertexArrays(1, &vao)
	gl.BindVertexArray(vao)
	gl.GenBuffers(1, &vboPos)
	gl.GenBuffers(1, &vboVal)
	gl.GenBuffers(1, &vboTag)
	gl.GenBuffers(1, &vboRad)

	var lvao, lvboPos, lvboW uint32
	gl.GenVertexArrays(1, &lvao)
	gl.BindVertexArray(lvao)
	gl.GenBuffers(1, &lvboPos)
	gl.GenBuffers(1, &lvboW)

	upload := func() {
		f := vw.frameFor()
		n := len(f.cells)
		pos := make([]float32, 3*n)
		rad := make([]float32, n)
		tag := make([]float32, n)
		for i, c := range f.cells {
			pos[3*i], pos[3*i+1], pos[3*i+2] = c.x, c.y, c.z
			rad[i] = c.r
			tag[i] = c.tag
		}
		vals, _ := vw.channelVals(f)
		gl.BindVertexArray(vao)
		gl.BindBuffer(gl.ARRAY_BUFFER, vboPos)
		gl.BufferData(gl.ARRAY_BUFFER, 4*len(pos), gl.Ptr(pos), gl.DYNAMIC_DRAW)
		gl.EnableVertexAttribArray(0)
		gl.VertexAttribPointer(0, 3, gl.FLOAT, false, 0, nil)
		gl.BindBuffer(gl.ARRAY_BUFFER, vboVal)
		gl.BufferData(gl.ARRAY_BUFFER, 4*len(vals), gl.Ptr(vals), gl.DYNAMIC_DRAW)
		gl.EnableVertexAttribArray(1)
		gl.VertexAttribPointer(1, 1, gl.FLOAT, false, 0, nil)
		gl.BindBuffer(gl.ARRAY_BUFFER, vboTag)
		gl.BufferData(gl.ARRAY_BUFFER, 4*len(tag), gl.Ptr(tag), gl.DYNAMIC_DRAW)
		gl.EnableVertexAttribArray(2)
		gl.VertexAttribPointer(2, 1, gl.FLOAT, false, 0, nil)
		gl.BindBuffer(gl.ARRAY_BUFFER, vboRad)
		gl.BufferData(gl.ARRAY_BUFFER, 4*len(rad), gl.Ptr(rad), gl.DYNAMIC_DRAW)
		gl.EnableVertexAttribArray(3)
		gl.VertexAttribPointer(3, 1, gl.FLOAT, false, 0, nil)
		// links
		if len(f.links) > 0 {
			lp := make([]float32, 0, 6*len(f.links))
			lw := make([]float32, 0, 2*len(f.links))
			half := float32(f.L / 2)
			for _, lk := range f.links {
				a, b := f.cells[lk.i], f.cells[lk.j]
				// skip links that wrap the periodic box (drawn across it)
				if abs32(a.x-b.x) > half || abs32(a.y-b.y) > half || abs32(a.z-b.z) > half {
					continue
				}
				w := lk.gg
				if lk.a <= 0 {
					w = 0.15 // dying channel: faint
				}
				lp = append(lp, a.x, a.y, a.z, b.x, b.y, b.z)
				lw = append(lw, w, w)
			}
			gl.BindVertexArray(lvao)
			gl.BindBuffer(gl.ARRAY_BUFFER, lvboPos)
			gl.BufferData(gl.ARRAY_BUFFER, 4*len(lp), gl.Ptr(lp), gl.DYNAMIC_DRAW)
			gl.EnableVertexAttribArray(0)
			gl.VertexAttribPointer(0, 3, gl.FLOAT, false, 0, nil)
			gl.BindBuffer(gl.ARRAY_BUFFER, lvboW)
			gl.BufferData(gl.ARRAY_BUFFER, 4*len(lw), gl.Ptr(lw), gl.DYNAMIC_DRAW)
			gl.EnableVertexAttribArray(1)
			gl.VertexAttribPointer(1, 1, gl.FLOAT, false, 0, nil)
		}
		vw.dirty = false
	}

	help := func() {
		fmt.Println("volview -i: drag orbit | shift+drag pan | wheel zoom | , . step | space play | [ ] speed")
		fmt.Println("  1..9/C channel (" + strings.Join(vw.channels, " ") + ") | L links | T tag tint | A avg | N norm | +/- size | G gamma | S shot | Q quit")
	}
	help()

	win.SetCursorPosCallback(func(_ *glfw.Window, x, y float64) {
		if vw.dragging {
			dx, dy := float32(x-vw.lastX), float32(y-vw.lastY)
			if vw.panning {
				s := vw.dist * 0.0012
				vw.panX -= dx * s
				vw.panY += dy * s
			} else {
				vw.yaw += dx * 0.008
				vw.pitch += dy * 0.008
				if vw.pitch > 1.55 {
					vw.pitch = 1.55
				}
				if vw.pitch < -1.55 {
					vw.pitch = -1.55
				}
			}
		}
		vw.lastX, vw.lastY = x, y
	})
	win.SetMouseButtonCallback(func(w *glfw.Window, b glfw.MouseButton, a glfw.Action, mods glfw.ModifierKey) {
		if b == glfw.MouseButtonLeft || b == glfw.MouseButtonRight || b == glfw.MouseButtonMiddle {
			vw.dragging = a == glfw.Press
			vw.panning = b != glfw.MouseButtonLeft || mods&glfw.ModShift != 0
		}
	})
	win.SetScrollCallback(func(_ *glfw.Window, _, dy float64) {
		vw.dist *= float32(math.Pow(0.9, dy))
		if vw.dist < 2 {
			vw.dist = 2
		}
	})
	win.SetKeyCallback(func(w *glfw.Window, k glfw.Key, _ int, a glfw.Action, mods glfw.ModifierKey) {
		if a != glfw.Press && a != glfw.Repeat {
			return
		}
		switch k {
		case glfw.KeyEscape, glfw.KeyQ:
			w.SetShouldClose(true)
		case glfw.KeySpace:
			vw.playing = !vw.playing
		case glfw.KeyComma, glfw.KeyLeft:
			vw.frameIdx = (vw.frameIdx + len(vw.frames) - 1) % len(vw.frames)
			vw.dirty = true
		case glfw.KeyPeriod, glfw.KeyRight:
			vw.frameIdx = (vw.frameIdx + 1) % len(vw.frames)
			vw.dirty = true
		case glfw.KeyLeftBracket:
			vw.speed = math.Max(1, vw.speed/1.5)
		case glfw.KeyRightBracket:
			vw.speed = math.Min(120, vw.speed*1.5)
		case glfw.KeyL:
			vw.links = !vw.links
		case glfw.KeyT:
			vw.tagTint = !vw.tagTint
		case glfw.KeyA:
			vw.useAvg = !vw.useAvg
			vw.dirty = true
		case glfw.KeyN:
			vw.globalN = !vw.globalN
			vw.dirty = true
		case glfw.KeyEqual, glfw.KeyKPAdd:
			vw.ptScale *= 1.15
		case glfw.KeyMinus, glfw.KeyKPSubtract:
			vw.ptScale /= 1.15
		case glfw.KeyG:
			if mods&glfw.ModShift != 0 {
				vw.gamma *= 1.15
			} else {
				vw.gamma /= 1.15
			}
		case glfw.KeyS:
			vw.shot++
		case glfw.KeyH:
			help()
		case glfw.KeyC:
			vw.chanIdx = (vw.chanIdx + 1) % len(vw.channels)
			vw.dirty = true
		default:
			if k >= glfw.Key1 && k <= glfw.Key9 {
				idx := int(k - glfw.Key1)
				if idx < len(vw.channels) {
					vw.chanIdx = idx
					vw.dirty = true
				}
			}
		}
	})

	gl.Enable(gl.PROGRAM_POINT_SIZE)
	gl.Enable(gl.BLEND)
	gl.BlendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA) // premultiplied
	gl.Disable(gl.DEPTH_TEST)

	lastStep := time.Now()
	nframes := 0
	for !win.ShouldClose() {
		if vw.playing && time.Since(lastStep).Seconds() > 1.0/vw.speed {
			vw.frameIdx = (vw.frameIdx + 1) % len(vw.frames)
			vw.dirty = true
			lastStep = time.Now()
		}
		if vw.dirty {
			upload()
			f := vw.frameFor()
			mode := "frame"
			if vw.useAvg {
				mode = "AVG"
			}
			win.SetTitle(fmt.Sprintf("fcs volview — t=%.2f  [%d/%d %s]  ch=%s  norm=%s",
				f.t, vw.frameIdx+1, len(vw.frames), mode, vw.channels[vw.chanIdx],
				map[bool]string{true: "global", false: "frame"}[vw.globalN]))
		}
		ww, wh := win.GetFramebufferSize()
		gl.Viewport(0, 0, int32(ww), int32(wh))
		gl.ClearColor(0.02, 0.015, 0.02, 1)
		gl.Clear(gl.COLOR_BUFFER_BIT)

		f := vw.frameFor()
		c := float32(f.L / 2)
		cp, sp := float32(math.Cos(float64(vw.pitch))), float32(math.Sin(float64(vw.pitch)))
		cy, sy := float32(math.Cos(float64(vw.yaw))), float32(math.Sin(float64(vw.yaw)))
		eye := [3]float32{c + vw.dist*cp*cy + vw.panX, c + vw.dist*cp*sy + vw.panY, c + vw.dist*sp}
		ctr := [3]float32{c + vw.panX, c + vw.panY, c}
		proj := perspective(0.9, float32(ww)/float32(wh), 0.1, 10*vw.dist)
		mvp := mul(proj, lookAt(eye, ctr))

		if vw.links && len(f.links) > 0 {
			gl.UseProgram(lprog)
			gl.UniformMatrix4fv(gl.GetUniformLocation(lprog, gl.Str("mvp\x00")), 1, false, &mvp[0])
			gl.BindVertexArray(lvao)
			nl := 0
			half := float32(f.L / 2)
			for _, lk := range f.links {
				a, b := f.cells[lk.i], f.cells[lk.j]
				if abs32(a.x-b.x) > half || abs32(a.y-b.y) > half || abs32(a.z-b.z) > half {
					continue
				}
				nl++
			}
			gl.DrawArrays(gl.LINES, 0, int32(2*nl))
		}

		gl.UseProgram(prog)
		gl.UniformMatrix4fv(gl.GetUniformLocation(prog, gl.Str("mvp\x00")), 1, false, &mvp[0])
		gl.Uniform1f(gl.GetUniformLocation(prog, gl.Str("ptScale\x00")), vw.ptScale)
		gl.Uniform1f(gl.GetUniformLocation(prog, gl.Str("gamma\x00")), vw.gamma)
			cyc := int32(0)
		if channelCyclic(vw.channels[vw.chanIdx]) {
			cyc = 1
		}
		gl.Uniform1i(gl.GetUniformLocation(prog, gl.Str("cyclic\x00")), cyc)
		tt := int32(0)
		if vw.tagTint {
			tt = 1
		}
		gl.Uniform1i(gl.GetUniformLocation(prog, gl.Str("tagTint\x00")), tt)
		gl.BindVertexArray(vao)
		gl.DrawArrays(gl.POINTS, 0, int32(len(f.cells)))

		if vw.shot > 0 {
			vw.shot--
			img := image.NewRGBA(image.Rect(0, 0, ww, wh))
			gl.ReadPixels(0, 0, int32(ww), int32(wh), gl.RGBA, gl.UNSIGNED_BYTE, unsafe.Pointer(&img.Pix[0]))
			// flip vertically
			for y := 0; y < wh/2; y++ {
				a := img.Pix[y*img.Stride : y*img.Stride+img.Stride]
				b := img.Pix[(wh-1-y)*img.Stride : (wh-1-y)*img.Stride+img.Stride]
				for i := range a {
					a[i], b[i] = b[i], a[i]
				}
			}
			name := fmt.Sprintf("volview_shot_%d.png", time.Now().Unix())
			if fo, err := os.Create(name); err == nil {
				png.Encode(fo, img)
				fo.Close()
				fmt.Println("# screenshot ->", name)
			}
		}

		win.SwapBuffers()
		glfw.PollEvents()
		nframes++
		if vw.selfcheck && nframes == 8 {
			vw.shot = 1
		}
		if vw.selfcheck && nframes >= 10 {
			fmt.Println("# selfcheck ok: rendered", nframes, "frames,",
				len(f.cells), "cells,", len(f.links), "links")
			break
		}
	}
}

func abs32(x float32) float32 {
	if x < 0 {
		return -x
	}
	return x
}
