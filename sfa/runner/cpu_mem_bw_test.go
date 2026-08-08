package main

import (
	"math"
	"testing"
)

func almostEq(a, b, tol float64) bool {
	return math.Abs(a-b) <= tol
}

func TestTheoreticalPeak(t *testing.T) {
	// Dual-channel DDR4-3200: 2 * 3200 * 8 / 1000 = 51.2
	if g := theoreticalPeakGBs(2, 3200); !almostEq(g, 51.2, 0.01) {
		t.Fatalf("dual DDR4-3200: got %v want 51.2", g)
	}
	// EPYC Rome 8ch DDR4-3200 = 204.8
	if g := theoreticalPeakGBs(8, 3200); !almostEq(g, 204.8, 0.01) {
		t.Fatalf("8ch DDR4-3200: got %v want 204.8", g)
	}
	// Genoa 12ch DDR5-4800 = 460.8
	if g := theoreticalPeakGBs(12, 4800); !almostEq(g, 460.8, 0.01) {
		t.Fatalf("12ch DDR5-4800: got %v want 460.8", g)
	}
	// Turin 12ch DDR5-6400 = 614.4
	if g := theoreticalPeakGBs(12, 6400); !almostEq(g, 614.4, 0.01) {
		t.Fatalf("12ch DDR5-6400: got %v want 614.4", g)
	}
}

func TestLookupCommonVastNames(t *testing.T) {
	cases := []struct {
		name      string
		hostCores float64
		familySub string
		peakMin   float64
		peakMax   float64
		sockets   int
	}{
		{"Xeon® E5-2680 v4", 56, "E5-2600 v4", 76, 77, 2}, // dual-socket from cores
		{"Xeon® E5-2673 v4", 28, "E5-2600 v4", 76, 77, 2},
		{"AMD EPYC 7742 64-Core Processor", 64, "Rome", 204, 205, 1},
		{"AMD EPYC 7742 64-Core Processor", 128, "Rome", 204, 205, 2}, // dual-socket cap
		{"AMD EPYC 9354 32-Core Processor", 128, "Genoa", 460, 461, 2}, // not 4S
		{"AMD EPYC 7K62 48-Core Processor", 48, "7K62", 204, 205, 1},
		{"AMD EPYC 7B13 64-Core Processor", 64, "7B13", 204, 205, 1},
		{"AMD EPYC 7302P 16-Core Processor", 16, "Rome", 204, 205, 1},
		{"AMD EPYC 9354 32-Core Processor", 32, "Genoa", 460, 461, 1},
		{"AMD Ryzen Threadripper PRO 3955WX 16-Cores", 16, "PRO 3000", 204, 205, 1},
		{"AMD Ryzen Threadripper 2970WX 24-Core Processor", 24, "Threadripper 2000", 93, 94, 1},
		{"AMD Ryzen 5 5600X 6-Core Processor", 6, "Ryzen 5000", 51, 52, 1},
		{"AMD Ryzen 7 7800X3D 8-Core Processor", 8, "Ryzen 7000", 83, 84, 1},
		{"Xeon® Gold 6138", 40, "Skylake", 127, 129, 1},
		{"Xeon® Gold 6430", 64, "Sapphire", 307, 308, 1},
		{"Xeon® W-2133", 12, "W-2100", 85, 86, 1},
		{"13th Gen Core™ i7-13700", 16, "13th/14th", 89, 90, 1},
		{"11th Gen Core™ i7-11700KF", 8, "11th", 51, 52, 1},
		{"Core™ i7-5930K", 6, "X99", 76, 77, 1},
	}
	for _, tc := range cases {
		bw := LookupCPUMemBW(tc.name, tc.hostCores)
		if bw == nil {
			t.Errorf("%q: no match", tc.name)
			continue
		}
		if tc.familySub != "" && !containsFold(bw.Family, tc.familySub) {
			t.Errorf("%q: family %q does not contain %q", tc.name, bw.Family, tc.familySub)
		}
		if bw.PeakGBs < tc.peakMin || bw.PeakGBs > tc.peakMax {
			t.Errorf("%q: peak %.2f not in [%.1f,%.1f] (family=%s)",
				tc.name, bw.PeakGBs, tc.peakMin, tc.peakMax, bw.Family)
		}
		if tc.sockets > 0 && bw.SocketsEst != tc.sockets {
			t.Errorf("%q: sockets=%d want %d (hostCores=%.0f modelCores=%d)",
				tc.name, bw.SocketsEst, tc.sockets, tc.hostCores, bw.ModelCores)
		}
		if bw.HostPeakGBs < bw.PeakGBs-0.01 {
			t.Errorf("%q: host peak %.2f < socket peak %.2f", tc.name, bw.HostPeakGBs, bw.PeakGBs)
		}
	}
}

func containsFold(s, sub string) bool {
	return len(sub) == 0 ||
		(len(s) >= len(sub) &&
			(func() bool {
				for i := 0; i+len(sub) <= len(s); i++ {
					if equalFoldASCII(s[i:i+len(sub)], sub) {
						return true
					}
				}
				return false
			})())
}

func equalFoldASCII(a, b string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := 0; i < len(a); i++ {
		ca, cb := a[i], b[i]
		if ca >= 'A' && ca <= 'Z' {
			ca += 'a' - 'A'
		}
		if cb >= 'A' && cb <= 'Z' {
			cb += 'a' - 'A'
		}
		if ca != cb {
			return false
		}
	}
	return true
}

func TestEmptyCPUName(t *testing.T) {
	if LookupCPUMemBW("", 0) != nil {
		t.Fatal("empty name should return nil")
	}
}
