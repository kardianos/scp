package main

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"
	"sync"
)

// CPUMemBW is the resolved DRAM bandwidth estimate for a host CPU.
// Units: GB/s (decimal, as in datasheet "theoretical peak").
//
// Vast does not expose DRAM STREAM. We derive peak from the published
// channel count × max official DRAM rate for the matched SKU/family:
//
//	peak_GBs = channels × (MT/s) × 8 / 1000
//
// STREAM triad typically lands at ~70–85% of peak when channels are fully
// populated; we surface both peak and an estimated STREAM floor.
type CPUMemBW struct {
	Match       string  `json:"match"`        // rule that hit
	Family      string  `json:"family"`       // e.g. "EPYC Rome (7002)"
	Channels    int     `json:"channels"`     // per socket
	MemMTs      int     `json:"mem_mts"`      // official max MT/s
	MemTech     string  `json:"mem_tech"`     // DDR4/DDR5/…
	PeakGBs     float64 `json:"peak_gbs"`     // theoretical peak per socket
	StreamEstGBs float64 `json:"stream_est_gbs"` // ~0.78 × peak (typical full-pop)
	SocketsEst  int     `json:"sockets_est"`  // 1 or 2 from host core count
	HostPeakGBs float64 `json:"host_peak_gbs"` // peak × sockets
	ModelCores  int     `json:"model_cores,omitempty"`
	Source      string  `json:"source"` // "builtin" | "override" | "formula"
	Confidence  string  `json:"confidence"` // high|medium|low
	Note        string  `json:"note,omitempty"`
}

// Theory: channels * MT/s * 8 bytes / 1000 → GB/s.
func theoreticalPeakGBs(channels, mts int) float64 {
	if channels <= 0 || mts <= 0 {
		return 0
	}
	return float64(channels) * float64(mts) * 8.0 / 1000.0
}

// streamFactor is a conservative STREAM-triad fraction of theoretical peak
// when all channels are populated with rated DIMMs.
const streamFactor = 0.78

// cpuFamily describes a generation's memory subsystem. Matched by regex on
// the Vast cpu_name string (after unicode cleanup).
type cpuFamily struct {
	// Match is a case-insensitive Go regexp. First match wins — order matters.
	Match string
	// Compiled at init.
	re *regexp.Regexp

	Family     string
	Channels   int
	MemMTs     int
	MemTech    string
	// ModelCores is typical full-SKU core count for dual-socket detection;
	// 0 means "unknown / do not estimate sockets from cores".
	ModelCores int
	// ModelCoresFromName: if true, parse "N-Core" / "N-Cores" from the name.
	ModelCoresFromName bool
	Confidence         string
	Note               string
}

// builtinCPUFamilies is ordered most-specific first.
// Numbers are official max per-socket configs (datasheet peak).
// Sources: AMD EPYC datasheets, Intel ARK, ServeTheHome channel×rate tables.
var builtinCPUFamilies = []cpuFamily{
	// --- AMD EPYC 9005 Turin ---
	{Match: `EPYC\s*9[0-9]{2}5`, Family: "EPYC Turin (9005)", Channels: 12, MemMTs: 6400, MemTech: "DDR5",
		ModelCoresFromName: true, Confidence: "high",
		Note: "12ch DDR5-6400 → 614.4 GB/s/socket theoretical"},
	// --- AMD EPYC 9004 Genoa / Bergamo / Siena ---
	{Match: `EPYC\s*9[0-9]{2}4`, Family: "EPYC Genoa/Bergamo (9004)", Channels: 12, MemMTs: 4800, MemTech: "DDR5",
		ModelCoresFromName: true, Confidence: "high",
		Note: "12ch DDR5-4800 → 460.8 GB/s/socket theoretical"},
	// Cloud / OEM Milan-Rome codenames (must beat generic 7xxx)
	{Match: `EPYC\s*7B13`, Family: "EPYC Milan cloud (7B13)", Channels: 8, MemMTs: 3200, MemTech: "DDR4",
		ModelCores: 64, Confidence: "high", Note: "GCP Milan-class 8ch DDR4-3200"},
	{Match: `EPYC\s*7B12`, Family: "EPYC Rome cloud (7B12)", Channels: 8, MemMTs: 3200, MemTech: "DDR4",
		ModelCores: 64, Confidence: "high", Note: "GCP Rome-class 8ch DDR4-3200"},
	{Match: `EPYC\s*7C13`, Family: "EPYC Milan cloud (7C13)", Channels: 8, MemMTs: 3200, MemTech: "DDR4",
		ModelCores: 64, Confidence: "high", Note: "Cloud Milan-class 8ch DDR4-3200"},
	{Match: `EPYC\s*7R32`, Family: "EPYC Rome cloud (7R32)", Channels: 8, MemMTs: 3200, MemTech: "DDR4",
		ModelCores: 48, Confidence: "high", Note: "AWS Rome-class 8ch DDR4-3200"},
	{Match: `EPYC\s*7V12`, Family: "EPYC Rome cloud (7V12)", Channels: 8, MemMTs: 3200, MemTech: "DDR4",
		ModelCores: 64, Confidence: "high", Note: "Azure Rome-class 8ch DDR4-3200"},
	{Match: `EPYC\s*7K62`, Family: "EPYC Rome cloud (7K62)", Channels: 8, MemMTs: 3200, MemTech: "DDR4",
		ModelCores: 48, Confidence: "medium", Note: "Alibaba Rome-class 8ch DDR4-3200"},
	// --- AMD EPYC 7003 Milan ---
	{Match: `EPYC\s*7[0-9]{2}3P?\b`, Family: "EPYC Milan (7003)", Channels: 8, MemMTs: 3200, MemTech: "DDR4",
		ModelCoresFromName: true, Confidence: "high",
		Note: "8ch DDR4-3200 → 204.8 GB/s/socket theoretical"},
	// --- AMD EPYC 7002 Rome ---
	{Match: `EPYC\s*7[0-9]{2}2P?\b`, Family: "EPYC Rome (7002)", Channels: 8, MemMTs: 3200, MemTech: "DDR4",
		ModelCoresFromName: true, Confidence: "high",
		Note: "8ch DDR4-3200 → 204.8 GB/s/socket theoretical"},
	// --- AMD EPYC 7001 Naples ---
	{Match: `EPYC\s*7[0-9]{2}1P?\b`, Family: "EPYC Naples (7001)", Channels: 8, MemMTs: 2666, MemTech: "DDR4",
		ModelCoresFromName: true, Confidence: "high",
		Note: "8ch DDR4-2666 → 170.6 GB/s/socket theoretical"},
	// Generic EPYC fallback (assume Rome/Milan 8ch — most marketplace hosts)
	{Match: `EPYC`, Family: "EPYC (unknown gen, assume 8ch DDR4-3200)", Channels: 8, MemMTs: 3200, MemTech: "DDR4",
		ModelCoresFromName: true, Confidence: "low",
		Note: "generation not parsed; assumed Rome/Milan-class"},

	// --- Threadripper PRO 7000 (WRX90, 8ch DDR5) ---
	{Match: `Threadripper\s+PRO\s+7[0-9]{3}`, Family: "Threadripper PRO 7000", Channels: 8, MemMTs: 5200, MemTech: "DDR5",
		ModelCoresFromName: true, Confidence: "high",
		Note: "8ch DDR5-5200 → 332.8 GB/s theoretical"},
	// --- Threadripper PRO 5000 ---
	{Match: `Threadripper\s+PRO\s+5[0-9]{3}`, Family: "Threadripper PRO 5000", Channels: 8, MemMTs: 3200, MemTech: "DDR4",
		ModelCoresFromName: true, Confidence: "high",
		Note: "8ch DDR4-3200 → 204.8 GB/s theoretical"},
	// --- Threadripper PRO 3000 ---
	{Match: `Threadripper\s+PRO\s+3[0-9]{3}`, Family: "Threadripper PRO 3000", Channels: 8, MemMTs: 3200, MemTech: "DDR4",
		ModelCoresFromName: true, Confidence: "high",
		Note: "8ch DDR4-3200 → 204.8 GB/s theoretical"},
	// --- Threadripper 9000 (Shimada Peak, 4ch DDR5) ---
	{Match: `Threadripper\s+9[0-9]{3}`, Family: "Threadripper 9000", Channels: 4, MemMTs: 6400, MemTech: "DDR5",
		ModelCoresFromName: true, Confidence: "high",
		Note: "4ch DDR5-6400 → 204.8 GB/s theoretical"},
	// --- Threadripper 7000 (non-PRO, 4ch DDR5) ---
	{Match: `Threadripper\s+7[0-9]{3}`, Family: "Threadripper 7000", Channels: 4, MemMTs: 5200, MemTech: "DDR5",
		ModelCoresFromName: true, Confidence: "high",
		Note: "4ch DDR5-5200 → 166.4 GB/s theoretical"},
	// --- Threadripper 3000 ---
	{Match: `Threadripper\s+3[0-9]{3}`, Family: "Threadripper 3000", Channels: 4, MemMTs: 3200, MemTech: "DDR4",
		ModelCoresFromName: true, Confidence: "high",
		Note: "4ch DDR4-3200 → 102.4 GB/s theoretical"},
	// --- Threadripper 2000 ---
	{Match: `Threadripper\s+2[0-9]{3}`, Family: "Threadripper 2000", Channels: 4, MemMTs: 2933, MemTech: "DDR4",
		ModelCoresFromName: true, Confidence: "high",
		Note: "4ch DDR4-2933 → 93.9 GB/s theoretical"},
	{Match: `Threadripper`, Family: "Threadripper (unknown)", Channels: 4, MemMTs: 3200, MemTech: "DDR4",
		ModelCoresFromName: true, Confidence: "low", Note: "assumed 4ch DDR4-3200"},

	// --- Ryzen 9000 (Zen5, dual-channel DDR5) ---
	{Match: `Ryzen\s+(?:9|7|5|3)\s*9[0-9]{3}`, Family: "Ryzen 9000", Channels: 2, MemMTs: 5600, MemTech: "DDR5",
		ModelCoresFromName: true, Confidence: "high",
		Note: "2ch DDR5-5600 → 89.6 GB/s theoretical (desktop)"},
	// --- Ryzen 7000 ---
	{Match: `Ryzen\s+(?:9|7|5|3)\s*7[0-9]{3}`, Family: "Ryzen 7000", Channels: 2, MemMTs: 5200, MemTech: "DDR5",
		ModelCoresFromName: true, Confidence: "high",
		Note: "2ch DDR5-5200 → 83.2 GB/s theoretical (desktop)"},
	// --- Ryzen 5000 ---
	{Match: `Ryzen\s+(?:9|7|5|3)\s*5[0-9]{3}`, Family: "Ryzen 5000", Channels: 2, MemMTs: 3200, MemTech: "DDR4",
		ModelCoresFromName: true, Confidence: "high",
		Note: "2ch DDR4-3200 → 51.2 GB/s theoretical (desktop)"},
	// --- Ryzen 3000 ---
	{Match: `Ryzen\s+(?:9|7|5|3)\s*3[0-9]{3}`, Family: "Ryzen 3000", Channels: 2, MemMTs: 3200, MemTech: "DDR4",
		ModelCoresFromName: true, Confidence: "high",
		Note: "2ch DDR4-3200 → 51.2 GB/s theoretical (desktop)"},
	{Match: `Ryzen`, Family: "Ryzen (unknown, assume dual-channel DDR4)", Channels: 2, MemMTs: 3200, MemTech: "DDR4",
		ModelCoresFromName: true, Confidence: "low"},

	// Intel Xeon Scalable model numbers: [tier digit][gen digit][xx]
	//   gen digit: 1=Skylake, 2=Cascade, 3=Ice Lake, 4=SPR, 5=EMR
	//   tier digit: 8=Platinum, 6=Gold high, 5=Gold/Silver, 4=Silver/Bronze
	// Order: classic gens first. "Xeon 6" (Granite) uses letter-suffixed SKUs (6700P…).
	// --- 5th gen Emerald Rapids (85xx / 65xx / 55xx) ---
	{Match: `Xeon.*(?:Platinum|Gold|Silver|Bronze)\s*[568]5[0-9]{2}\b`, Family: "Xeon Emerald Rapids (5th)", Channels: 8, MemMTs: 5600, MemTech: "DDR5",
		Confidence: "high", Note: "8ch DDR5-5600 → 358.4 GB/s/socket"},
	// --- 4th gen Sapphire Rapids (84xx / 64xx / 54xx / 44xx) ---
	{Match: `Xeon.*(?:Platinum|Gold|Silver|Bronze)\s*[4568]4[0-9]{2}\b`, Family: "Xeon Sapphire Rapids (4th)", Channels: 8, MemMTs: 4800, MemTech: "DDR5",
		Confidence: "high", Note: "8ch DDR5-4800 → 307.2 GB/s/socket"},
	// --- 3rd gen Ice Lake (83xx / 63xx / 53xx) ---
	{Match: `Xeon.*(?:Platinum|Gold|Silver|Bronze)\s*[568]3[0-9]{2}\b`, Family: "Xeon Ice Lake (3rd)", Channels: 8, MemMTs: 3200, MemTech: "DDR4",
		Confidence: "high", Note: "8ch DDR4-3200 → 204.8 GB/s/socket"},
	// --- 2nd gen Cascade Lake (82xx / 62xx / 52xx) ---
	{Match: `Xeon.*(?:Platinum|Gold|Silver|Bronze)\s*[568]2[0-9]{2}\b`, Family: "Xeon Cascade Lake (2nd)", Channels: 6, MemMTs: 2933, MemTech: "DDR4",
		Confidence: "high", Note: "6ch DDR4-2933 → 140.8 GB/s/socket"},
	// --- 1st gen Skylake-SP (81xx / 61xx / 51xx / 41xx) ---
	{Match: `Xeon.*(?:Platinum|Gold|Silver|Bronze)\s*[4568]1[0-9]{2}\b`, Family: "Xeon Skylake-SP (1st)", Channels: 6, MemMTs: 2666, MemTech: "DDR4",
		Confidence: "high", Note: "6ch DDR4-2666 → 128.0 GB/s/socket"},
	// --- Xeon 6 Granite Rapids (letter-suffixed SKUs: 6700P, 6780E, …) ---
	{Match: `Xeon.*\b6[0-9]{3}[A-Z]\b`, Family: "Xeon 6 (Granite Rapids)", Channels: 8, MemMTs: 6400, MemTech: "DDR5",
		Confidence: "medium", Note: "8ch DDR5-6400 → 409.6 GB/s/socket (mainstream SP)"},
	// Generic Gold/Platinum without parseable number
	{Match: `Xeon.*(?:Platinum|Gold)`, Family: "Xeon Scalable (unknown gen)", Channels: 6, MemMTs: 2666, MemTech: "DDR4",
		Confidence: "low", Note: "assumed Skylake-class 6ch"},

	// --- Xeon W workstation ---
	{Match: `Xeon®?\s*W-33[0-9]{2}`, Family: "Xeon W-3300", Channels: 8, MemMTs: 3200, MemTech: "DDR4",
		Confidence: "high", Note: "8ch DDR4-3200 → 204.8 GB/s"},
	{Match: `Xeon®?\s*W-32[0-9]{2}`, Family: "Xeon W-3200", Channels: 6, MemMTs: 2666, MemTech: "DDR4",
		Confidence: "high", Note: "6ch DDR4-2666 → 128.0 GB/s"},
	{Match: `Xeon®?\s*W-22[0-9]{2}`, Family: "Xeon W-2200", Channels: 4, MemMTs: 2933, MemTech: "DDR4",
		Confidence: "high", Note: "4ch DDR4-2933 → 93.9 GB/s"},
	{Match: `Xeon®?\s*W-21[0-9]{2}`, Family: "Xeon W-2100", Channels: 4, MemMTs: 2666, MemTech: "DDR4",
		Confidence: "high", Note: "4ch DDR4-2666 → 85.3 GB/s"},

	// --- Xeon E5 v4 / v3 / v2 (Haswell/Broadwell/Ivy Bridge EP) ---
	{Match: `E5-26[0-9]{2}\s*v4`, Family: "Xeon E5-2600 v4", Channels: 4, MemMTs: 2400, MemTech: "DDR4",
		Confidence: "high", Note: "4ch DDR4-2400 → 76.8 GB/s/socket (dual-socket common on Vast)"},
	{Match: `E5-16[0-9]{2}\s*v4`, Family: "Xeon E5-1600 v4", Channels: 4, MemMTs: 2400, MemTech: "DDR4",
		Confidence: "high", Note: "4ch DDR4-2400 → 76.8 GB/s"},
	{Match: `E5-26[0-9]{2}\s*v3`, Family: "Xeon E5-2600 v3", Channels: 4, MemMTs: 2133, MemTech: "DDR4",
		Confidence: "high", Note: "4ch DDR4-2133 → 68.3 GB/s/socket"},
	{Match: `E5-16[0-9]{2}\s*v3`, Family: "Xeon E5-1600 v3", Channels: 4, MemMTs: 2133, MemTech: "DDR4",
		Confidence: "high", Note: "4ch DDR4-2133 → 68.3 GB/s"},
	{Match: `E5-26[0-9]{2}\s*v2`, Family: "Xeon E5-2600 v2", Channels: 4, MemMTs: 1866, MemTech: "DDR3",
		Confidence: "high", Note: "4ch DDR3-1866 → 59.7 GB/s/socket"},
	{Match: `E5-26[0-9]{2}\b`, Family: "Xeon E5-2600 (v1/Sandy Bridge)", Channels: 4, MemMTs: 1600, MemTech: "DDR3",
		Confidence: "medium", Note: "4ch DDR3-1600 → 51.2 GB/s/socket"},
	{Match: `E5-16[0-9]{2}`, Family: "Xeon E5-1600", Channels: 4, MemMTs: 1866, MemTech: "DDR3",
		Confidence: "medium"},

	// --- Intel Core Ultra 200 (Arrow Lake) ---
	{Match: `Core™?\s*Ultra\s*[579]\s*2[0-9]{2}`, Family: "Core Ultra 200", Channels: 2, MemMTs: 6400, MemTech: "DDR5",
		Confidence: "high", Note: "2ch DDR5-6400 → 102.4 GB/s theoretical"},
	// --- 14th / 13th / 12th gen Core ---
	{Match: `(?:13th|14th)\s+Gen\s+Core|Core™?\s*i[3579]-1[34][0-9]{3}`, Family: "Core 13th/14th gen", Channels: 2, MemMTs: 5600, MemTech: "DDR5",
		Confidence: "medium", Note: "2ch DDR5-5600 → 89.6 GB/s (DDR4 configs lower)"},
	{Match: `12th\s+Gen\s+Core|Core™?\s*i[3579]-12[0-9]{3}`, Family: "Core 12th gen", Channels: 2, MemMTs: 4800, MemTech: "DDR5",
		Confidence: "medium", Note: "2ch DDR5-4800 → 76.8 GB/s (DDR4 configs lower)"},
	{Match: `11th\s+Gen\s+Core|Core™?\s*i[3579]-11[0-9]{3}`, Family: "Core 11th gen", Channels: 2, MemMTs: 3200, MemTech: "DDR4",
		Confidence: "high", Note: "2ch DDR4-3200 → 51.2 GB/s"},
	{Match: `10th\s+Gen\s+Core|Core™?\s*i[3579]-10[0-9]{3}`, Family: "Core 10th gen", Channels: 2, MemMTs: 2933, MemTech: "DDR4",
		Confidence: "high", Note: "2ch DDR4-2933 → 46.9 GB/s"},
	{Match: `Core™?\s*i[3579]-9[0-9]{3}`, Family: "Core 9th gen", Channels: 2, MemMTs: 2666, MemTech: "DDR4",
		Confidence: "high", Note: "2ch DDR4-2666 → 42.7 GB/s"},
	{Match: `Core™?\s*i[3579]-8[0-9]{3}`, Family: "Core 8th gen", Channels: 2, MemMTs: 2666, MemTech: "DDR4",
		Confidence: "high"},
	{Match: `Core™?\s*i[3579]-7[0-9]{3}`, Family: "Core 7th gen", Channels: 2, MemMTs: 2400, MemTech: "DDR4",
		Confidence: "high"},
	{Match: `Core™?\s*i[3579]-6[0-9]{3}`, Family: "Core 6th gen", Channels: 2, MemMTs: 2133, MemTech: "DDR4",
		Confidence: "high"},
	// X99 HEDT (i7-58xx/68xx/69xx) — 4 channel
	{Match: `Core™?\s*i7-[56][89][0-9]{2}`, Family: "Core i7 X99 HEDT", Channels: 4, MemMTs: 2400, MemTech: "DDR4",
		Confidence: "high", Note: "4ch DDR4-2400 → 76.8 GB/s (X99)"},
	{Match: `Core™?\s*i[3579]-`, Family: "Intel Core (unknown gen)", Channels: 2, MemMTs: 2400, MemTech: "DDR4",
		Confidence: "low"},
	{Match: `Pentium`, Family: "Pentium", Channels: 2, MemMTs: 2400, MemTech: "DDR4",
		Confidence: "low"},
}

var (
	cpuBWOnce     sync.Once
	cpuBWFamilies []cpuFamily
	// Exact/substring overrides loaded from ~/.scp-runner/cpu_mem_bw.json
	cpuBWOverrides []cpuBWOverride
)

// cpuBWOverride is a user-supplied exact or regex override.
// Place at ~/.scp-runner/cpu_mem_bw.json:
//
//	{"overrides":[
//	  {"match":"EPYC 7742","channels":8,"mem_mts":3200,"mem_tech":"DDR4","family":"EPYC Rome","model_cores":64}
//	]}
type cpuBWOverride struct {
	Match      string `json:"match"` // substring or /regex/
	Family     string `json:"family"`
	Channels   int    `json:"channels"`
	MemMTs     int    `json:"mem_mts"`
	MemTech    string `json:"mem_tech"`
	ModelCores int    `json:"model_cores"`
	Note       string `json:"note"`
	IsRegex    bool   `json:"-"`
	re         *regexp.Regexp
}

func initCPUMemBW() {
	cpuBWOnce.Do(func() {
		cpuBWFamilies = make([]cpuFamily, len(builtinCPUFamilies))
		copy(cpuBWFamilies, builtinCPUFamilies)
		for i := range cpuBWFamilies {
			cpuBWFamilies[i].re = regexp.MustCompile(`(?i)` + cpuBWFamilies[i].Match)
		}
		loadCPUMemBWOverrides()
	})
}

func loadCPUMemBWOverrides() {
	home, err := os.UserHomeDir()
	if err != nil {
		return
	}
	paths := []string{
		filepath.Join(home, ".scp-runner", "cpu_mem_bw.json"),
		filepath.Join(home, ".config", "scp-runner", "cpu_mem_bw.json"),
	}
	for _, p := range paths {
		data, err := os.ReadFile(p)
		if err != nil {
			continue
		}
		var doc struct {
			Overrides []cpuBWOverride `json:"overrides"`
		}
		if err := json.Unmarshal(data, &doc); err != nil {
			fmt.Fprintf(os.Stderr, "scp-runner: cpu_mem_bw override %s: %v\n", p, err)
			continue
		}
		for i := range doc.Overrides {
			o := &doc.Overrides[i]
			m := o.Match
			if strings.HasPrefix(m, "/") && strings.HasSuffix(m, "/") && len(m) > 2 {
				o.IsRegex = true
				o.re = regexp.MustCompile(`(?i)` + m[1:len(m)-1])
			}
			cpuBWOverrides = append(cpuBWOverrides, *o)
		}
		fmt.Fprintf(os.Stderr, "scp-runner: loaded %d CPU mem-bw overrides from %s\n", len(doc.Overrides), p)
		return
	}
}

var (
	reCoreCount = regexp.MustCompile(`(?i)(\d+)\s*-?\s*Cores?\b`)
	reCleanup   = regexp.MustCompile(`\s+`)
)

// normalizeCPUName strips trademark glyphs and collapses whitespace.
func normalizeCPUName(name string) string {
	r := strings.NewReplacer(
		"®", "", "™", "", "©", "",
		"（", "(", "）", ")",
	)
	s := r.Replace(name)
	s = reCleanup.ReplaceAllString(strings.TrimSpace(s), " ")
	return s
}

// LookupCPUMemBW maps a Vast cpu_name (+ optional host total core count) to
// theoretical DRAM bandwidth. hostCores is Vast's cpu_cores field (whole machine).
func LookupCPUMemBW(cpuName string, hostCores float64) *CPUMemBW {
	initCPUMemBW()
	name := normalizeCPUName(cpuName)
	if name == "" {
		return nil
	}

	// 1) User overrides first.
	for _, o := range cpuBWOverrides {
		hit := false
		if o.IsRegex && o.re != nil {
			hit = o.re.MatchString(name)
		} else {
			hit = strings.Contains(strings.ToUpper(name), strings.ToUpper(o.Match))
		}
		if !hit {
			continue
		}
		return finishCPUMemBW(name, o.Match, o.Family, o.Channels, o.MemMTs, o.MemTech,
			o.ModelCores, hostCores, "override", "high", o.Note)
	}

	// 2) Builtin family table.
	for _, f := range cpuBWFamilies {
		if f.re == nil || !f.re.MatchString(name) {
			continue
		}
		cores := f.ModelCores
		if f.ModelCoresFromName {
			if m := reCoreCount.FindStringSubmatch(name); len(m) == 2 {
				if n, err := strconv.Atoi(m[1]); err == nil {
					cores = n
				}
			}
		}
		return finishCPUMemBW(name, f.Match, f.Family, f.Channels, f.MemMTs, f.MemTech,
			cores, hostCores, "builtin", f.Confidence, f.Note)
	}
	return nil
}

func finishCPUMemBW(name, match, family string, ch, mts int, tech string,
	modelCores int, hostCores float64, source, conf, note string) *CPUMemBW {

	peak := theoreticalPeakGBs(ch, mts)
	sockets := 1
	// EPYC "P" and many cloud SKUs are single-socket only.
	singleSocketOnly := regexp.MustCompile(`(?i)EPYC\s*7[0-9A-Z]{2,4}P\b`).MatchString(name) ||
		regexp.MustCompile(`(?i)EPYC\s*9[0-9]{3}P\b`).MatchString(name)

	// Dual-socket estimate: host core count ≈ 2× model core count.
	// Cap at 2 — Vast GPU hosts are almost never 4S, and multi-GPU single-socket
	// machines often advertise inflated cpu_cores that would falsely multiply DRAM BW.
	if !singleSocketOnly && modelCores > 0 && hostCores >= float64(modelCores)*1.8 {
		sockets = 2
		if note != "" {
			note += "; "
		}
		note += fmt.Sprintf("assumed dual-socket (host cores %.0f ≥ 1.8× model %d)", hostCores, modelCores)
	}
	// E5-2600 class on Vast is almost always dual-socket when cpu_cores ≥ 28.
	if sockets == 1 && modelCores == 0 && hostCores >= 28 &&
		strings.Contains(strings.ToUpper(name), "E5-26") {
		sockets = 2
		if note != "" {
			note += "; "
		}
		note += "assumed dual-socket from host core count ≥ 28"
	}
	if singleSocketOnly {
		sockets = 1
	}

	out := &CPUMemBW{
		Match:        match,
		Family:       family,
		Channels:     ch,
		MemMTs:       mts,
		MemTech:      tech,
		PeakGBs:      peak,
		StreamEstGBs: peak * streamFactor,
		SocketsEst:   sockets,
		HostPeakGBs:  peak * float64(sockets),
		ModelCores:   modelCores,
		Source:       source,
		Confidence:   conf,
		Note:         note,
	}
	return out
}

// FormatCPUMemBW is a short human summary.
func FormatCPUMemBW(b *CPUMemBW) string {
	if b == nil {
		return "cpu_mem_bw=?"
	}
	return fmt.Sprintf("dram~%.0fGB/s peak (×%dsock=%.0f host) stream~%.0f [%s %dch %s-%d conf=%s]",
		b.PeakGBs, b.SocketsEst, b.HostPeakGBs, b.StreamEstGBs*float64(b.SocketsEst),
		b.Family, b.Channels, b.MemTech, b.MemMTs, b.Confidence)
}
