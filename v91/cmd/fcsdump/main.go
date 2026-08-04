// fcsdump — FCS v3 stream reader for the pad campaign (pure analysis tool).
// Modes: -mode stats (per-frame summaries), cells (full per-cell rows),
// links (full per-link rows). -hot sets the hotspot Em threshold.
package main

import (
	"bufio"
	"encoding/binary"
	"flag"
	"fmt"
	"io"
	"math"
	"os"

	"github.com/klauspost/compress/zstd"
)

var (
	mode    = flag.String("mode", "stats", "stats|cells|links")
	hot     = flag.Float64("hot", 0.5, "hotspot Em threshold for stats")
	tagonly = flag.Bool("tagonly", false, "cells mode: tagged cells only")
)

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
	n := len(body) / (4 * stride)
	o := out[hdr:]
	for c := 0; c < stride; c++ {
		for i := 0; i < n; i++ {
			copy(o[(i*stride+c)*4:(i*stride+c)*4+4], body[(c*n+i)*4:])
		}
	}
	return out
}

func main() {
	flag.Parse()
	if flag.NArg() < 1 {
		fmt.Fprintln(os.Stderr, "usage: fcsdump [-mode stats|cells|links] file.fcs")
		os.Exit(2)
	}
	f, err := os.Open(flag.Arg(0))
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	r := bufio.NewReaderSize(f, 1<<20)
	zdec, _ := zstd.NewReader(nil)
	magic := make([]byte, 8)
	if _, err := io.ReadFull(r, magic); err != nil || string(magic[:4]) != "FCS1" {
		fmt.Fprintln(os.Stderr, "not an FCS v3 stream")
		os.Exit(1)
	}
	f32 := func(b []byte) float64 { return float64(math.Float32frombits(binary.LittleEndian.Uint32(b))) }
	f64 := func(b []byte) float64 { return math.Float64frombits(binary.LittleEndian.Uint64(b)) }
	u32 := func(b []byte) uint32 { return binary.LittleEndian.Uint32(b) }
	if *mode == "cells" {
		fmt.Println("t i x y z r es em ee xload tag")
	} else if *mode == "links" {
		fmt.Println("t i j d A lem gg")
	} else {
		fmt.Println("kind t N sumEm meanEm varEm maxEm nHot sumEe sumEs meanXload tagN tagEm | links: nl meanGG nLive meanD")
	}
	for {
		tag := make([]byte, 4)
		if _, err := io.ReadFull(r, tag); err != nil {
			break
		}
		lb := make([]byte, 8)
		if _, err := io.ReadFull(r, lb); err != nil {
			break
		}
		n := binary.LittleEndian.Uint64(lb)
		p := make([]byte, n)
		if _, err := io.ReadFull(r, p); err != nil {
			break
		}
		ct := string(tag)
		if ct == "CMPD" {
			inner := string(p[:4])
			rawlen := binary.LittleEndian.Uint64(p[4:12])
			codec := p[12]
			if codec != 1 {
				continue
			}
			raw, err := zdec.DecodeAll(p[13:], nil)
			if err != nil || uint64(len(raw)) != rawlen {
				fmt.Fprintln(os.Stderr, "decode error", err)
				continue
			}
			hdr, stride := 20, 12
			if inner == "LINK" {
				hdr, stride = 12, 6
			}
			p = unColT(unshuffle4(raw), hdr, stride)
			ct = inner
		}
		switch ct {
		case "CELL":
			t := f64(p[0:8])
			N := int(u32(p[16:20]))
			rec := p[20:]
			var sEm, sEe, sEs, sXl, vEm, mEm, tEm float64
			var nHot, tagN int
			for i := 0; i < N; i++ {
				b := rec[i*48:]
				es, em, ee, xl, tg := f32(b[16:]), f32(b[20:]), f32(b[24:]), f32(b[28:]), f32(b[32:])
				if *mode == "cells" {
					if !*tagonly || tg > 0.5 {
						fmt.Printf("%.2f %d %.4f %.4f %.4f %.4f %.5f %.5f %.5f %.5f %.0f\n",
							t, i, f32(b[0:]), f32(b[4:]), f32(b[8:]), f32(b[12:]), es, em, ee, xl, tg)
					}
					continue
				}
				sEm += em
				sEe += ee
				sEs += es
				sXl += xl
				if em > mEm {
					mEm = em
				}
				if em > *hot {
					nHot++
				}
				if tg > 0.5 {
					tagN++
					tEm += em
				}
			}
			if *mode == "stats" {
				mean := sEm / float64(N)
				for i := 0; i < N; i++ {
					em := f32(rec[i*48+20:])
					vEm += (em - mean) * (em - mean)
				}
				vEm /= float64(N)
				fmt.Printf("CELL %.2f %d %.4f %.6f %.6f %.4f %d %.4f %.4f %.6f %d %.4f\n",
					t, N, sEm, mean, vEm, mEm, nHot, sEe, sEs, sXl/float64(N), tagN, tEm)
			}
		case "LINK":
			t := f64(p[0:8])
			nl := int(u32(p[8:12]))
			rec := p[12:]
			var sGG, sD float64
			var live int
			for i := 0; i < nl; i++ {
				b := rec[i*24:]
				ii, jj := u32(b[0:4]), u32(b[4:8])
				d, A, lem, gg := f32(b[8:]), f32(b[12:]), f32(b[16:]), f32(b[20:])
				if *mode == "links" {
					fmt.Printf("%.2f %d %d %.4f %.4f %.5f %.4f\n", t, ii, jj, d, A, lem, gg)
					continue
				}
				sGG += gg
				sD += d
				if gg > 0.5 {
					live++
				}
			}
			if *mode == "stats" && nl > 0 {
				fmt.Printf("LINK %.2f %d %.4f %d %.4f\n", t, nl, sGG/float64(nl), live, sD/float64(nl))
			}
		}
	}
}
