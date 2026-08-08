package main

import (
	"context"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"sync"
	"time"
)

const vastAPIBase = "https://console.vast.ai/api/v0"

// VastClient talks to the Vast.ai REST API.
type VastClient struct {
	apiKey     string
	httpClient *http.Client
}

// VastOffer is a GPU rental offer. Fields map to Vast /bundles/ response keys.
// Bandwidth units (from Vast docs):
//   - disk_bw:     disk read bandwidth in MB/s
//   - gpu_mem_bw:  GPU memory bandwidth in GB/s
//   - pcie_bw:     PCIe bandwidth (CPU↔GPU)
//   - inet_*:      network in Mb/s (CLI) / MB/s (REST docs vary)
// There is NO system-DRAM STREAM bandwidth field in the Vast API.
type VastOffer struct {
	ID                 int     `json:"id"`
	MachineID          int     `json:"machine_id"`
	GPUName            string  `json:"gpu_name"`
	NumGPUs            int     `json:"num_gpus"`
	DPHTot             float64 `json:"dph_total"`
	DiskSpace          float64 `json:"disk_space"`
	GPUMemMB           int64   `json:"gpu_ram"`
	Reliability        float64 `json:"reliability2"`
	InetDown           float64 `json:"inet_down"`
	InetUp             float64 `json:"inet_up"`
	Geolocation        string  `json:"geolocation"`
	CPUName            string  `json:"cpu_name"`
	CPUArch            string  `json:"cpu_arch"`
	CPUCores           float64 `json:"cpu_cores"`
	CPUCoresEffective  float64 `json:"cpu_cores_effective"`
	CPUGHz             float64 `json:"cpu_ghz"`
	CPURamMB           float64 `json:"cpu_ram"` // API returns MB
	HasAVX             int     `json:"has_avx"`
	DiskBW             float64 `json:"disk_bw"`    // MB/s
	GPUMemBW           float64 `json:"gpu_mem_bw"` // GB/s
	PCIeBW             float64 `json:"pcie_bw"`
	DiskName           string  `json:"disk_name"`
	MoboName           string  `json:"mobo_name"`
	PCIGen             float64 `json:"pci_gen"`
	// Filled client-side from cpu_mem_bw lookup (not from Vast API).
	CPUClass       string  `json:"cpu_class,omitempty"`
	CPUNote        string  `json:"cpu_note,omitempty"`
	CPUMemPeakGBs  float64 `json:"cpu_mem_peak_gbs,omitempty"`  // theoretical peak / socket
	CPUMemHostGBs  float64 `json:"cpu_mem_host_gbs,omitempty"`  // peak × estimated sockets
	CPUMemStreamGBs float64 `json:"cpu_mem_stream_gbs,omitempty"` // ~0.78 × host peak
	CPUMemChannels int     `json:"cpu_mem_channels,omitempty"`
	CPUMemMTs      int     `json:"cpu_mem_mts,omitempty"`
	CPUMemTech     string  `json:"cpu_mem_tech,omitempty"`
	CPUMemConf     string  `json:"cpu_mem_confidence,omitempty"`
}

// machineBlacklist is an in-memory, session-scoped set of Vast machine IDs
// that failed provisioning (ready-timeout, SSH auth failure, onstart failure).
// Blacklisted machines are excluded from subsequent offer selection.
var machineBlacklist = struct {
	mu  sync.Mutex
	ids map[int]bool
}{ids: make(map[int]bool)}

// blacklistMachine records a failed machine ID for this session.
func blacklistMachine(machineID int) {
	if machineID == 0 {
		return
	}
	machineBlacklist.mu.Lock()
	machineBlacklist.ids[machineID] = true
	machineBlacklist.mu.Unlock()
}

// isMachineBlacklisted reports whether a machine ID failed earlier this session.
func isMachineBlacklisted(machineID int) bool {
	machineBlacklist.mu.Lock()
	defer machineBlacklist.mu.Unlock()
	return machineBlacklist.ids[machineID]
}

// allowedGPUs defines the hardcoded set of acceptable GPU models with their
// minimum VRAM in MB and CUDA compute capability.
var allowedGPUs = []struct {
	Name   string // Vast.ai gpu_name (exact match)
	MinRAM int64  // Minimum VRAM in MB to accept
	Arch   string // CUDA arch (sm_XX) for reference
}{
	{"Tesla V100", 16 * 1024, "sm_70"},
	{"A100 SXM4", 40 * 1024, "sm_80"},
	{"A100 PCIE", 40 * 1024, "sm_80"},
	{"L40S", 45 * 1024, "sm_89"},
	{"H100 SXM", 80 * 1024, "sm_90"},
	{"H100 PCIE", 80 * 1024, "sm_90"},
	{"RTX PRO 4500", 32 * 1024, "sm_120"},
	{"RTX PRO 4000", 24 * 1024, "sm_120"},
	{"B200", 179 * 1024, "sm_100"},
}

// allowedRegions restricts provisioning to these country codes.
var allowedRegions = []string{"US", "CA"}

// searchFilterKeys is the set of Vast /bundles/ filter fields we forward to the API.
// cpu_name is intentionally absent: Vast returns it but rejects it as a filter.
var searchFilterKeys = map[string]bool{
	"verified": true, "external": true, "rented": true, "rentable": true,
	"gpu_name": true, "num_gpus": true, "gpu_ram": true, "gpu_total_ram": true,
	"gpu_mem_bw": true, "gpu_arch": true, "compute_cap": true, "cuda_max_good": true,
	"dph_total": true, "min_bid": true, "reliability": true, "reliability2": true,
	"disk_space": true, "disk_bw": true, "duration": true,
	"inet_down": true, "inet_up": true, "inet_down_cost": true, "inet_up_cost": true,
	"cpu_arch": true, "cpu_cores": true, "cpu_cores_effective": true,
	"cpu_ghz": true, "cpu_ram": true, "has_avx": true,
	"geolocation": true, "datacenter": true, "machine_id": true, "host_id": true,
	"pcie_bw": true, "pci_gen": true, "bw_nvlink": true, "total_flops": true,
	"dlperf": true, "static_ip": true, "direct_port_count": true,
	"driver_version": true, "ubuntu_version": true, "os_version": true,
	"vms_enabled": true, "verification": true, "mobo_name": true,
}

// localOnlyKeys are filter tokens handled client-side, not sent to Vast.
var localOnlyKeys = map[string]bool{
	"min_ram": true, "any_gpu": true, "whitelist": true,
	"region": true, "cpu_name": true, "cpu_contains": true,
	"max_dph": true, "max_price": true,
	"min_disk_bw": true, "min_cpu_cores": true, "min_cpu_ghz": true,
	"min_cpu_ram": true, "min_pcie_bw": true, "min_gpu_mem_bw": true,
	"min_cpu_mem_bw": true, "min_dram_bw": true, "min_mem_bw": true,
	"min_inet_down": true, "limit": true, "order": true,
}

// annotateCPU fills DRAM bandwidth fields from the local CPU lookup table.
func annotateCPU(o *VastOffer) {
	bw := LookupCPUMemBW(o.CPUName, o.CPUCores)
	if bw == nil {
		o.CPUClass = "unknown"
		o.CPUNote = "no CPU mem-bw match; add an override in ~/.scp-runner/cpu_mem_bw.json"
		return
	}
	o.CPUClass = bw.Family
	o.CPUNote = bw.Note
	o.CPUMemPeakGBs = bw.PeakGBs
	o.CPUMemHostGBs = bw.HostPeakGBs
	o.CPUMemStreamGBs = bw.StreamEstGBs * float64(bw.SocketsEst)
	o.CPUMemChannels = bw.Channels
	o.CPUMemMTs = bw.MemMTs
	o.CPUMemTech = bw.MemTech
	o.CPUMemConf = bw.Confidence
}

// VastInstance is a running instance.
type VastInstance struct {
	ID         int     `json:"id"`
	MachineID  int     `json:"machine_id"`
	Status     string  `json:"actual_status"`
	SSHHost    string  `json:"ssh_host"`
	SSHPort    int     `json:"ssh_port"`
	GPUName    string  `json:"gpu_name"`
	CurState   string  `json:"cur_state"`
	ContractID int     `json:"contract_id,omitempty"`
	DPHTotal   float64 `json:"dph_total,omitempty"`
}

// NewVastClient creates a client, reading the API key from ~/.vast_api_key.
func NewVastClient() (*VastClient, error) {
	home, err := os.UserHomeDir()
	if err != nil {
		return nil, fmt.Errorf("home dir: %w", err)
	}
	keyPaths := []string{
		filepath.Join(home, ".vast_api_key"),
		filepath.Join(home, ".config", "vastai", "vast_api_key"),
	}
	var data []byte
	for _, kp := range keyPaths {
		d, err := os.ReadFile(kp)
		if err == nil {
			data = d
			break
		}
	}
	if data == nil {
		key := os.Getenv("VAST_API_KEY")
		if key == "" {
			return nil, fmt.Errorf("no API key found in %v (and VAST_API_KEY not set)", keyPaths)
		}
		data = []byte(key)
	}
	return &VastClient{
		apiKey:     strings.TrimSpace(string(data)),
		httpClient: &http.Client{Timeout: 30 * time.Second},
	}, nil
}

func (v *VastClient) doRequest(ctx context.Context, method, path string, body io.Reader) ([]byte, error) {
	u := vastAPIBase + path
	req, err := http.NewRequestWithContext(ctx, method, u, body)
	if err != nil {
		return nil, fmt.Errorf("new request: %w", err)
	}
	req.Header.Set("Authorization", "Bearer "+v.apiKey)
	if body != nil {
		req.Header.Set("Content-Type", "application/json")
	}

	resp, err := v.httpClient.Do(req)
	if err != nil {
		return nil, fmt.Errorf("http do: %w", err)
	}
	defer resp.Body.Close()

	data, err := io.ReadAll(resp.Body)
	if err != nil {
		return nil, fmt.Errorf("read body: %w", err)
	}
	if resp.StatusCode >= 400 {
		return nil, fmt.Errorf("vast api %s %s: %d: %s", method, path, resp.StatusCode, string(data))
	}
	return data, nil
}

// gpuFilterMapToString converts a map[string]string GPU filter to the legacy
// whitespace-separated string format used by resolveSearchSpec.
func gpuFilterMapToString(m map[string]string) string {
	var parts []string
	for k, v := range m {
		// Preserve operator if the value already starts with one.
		if strings.HasPrefix(v, ">=") || strings.HasPrefix(v, "<=") ||
			strings.HasPrefix(v, "!=") || strings.HasPrefix(v, ">") ||
			strings.HasPrefix(v, "<") || strings.HasPrefix(v, "=") {
			parts = append(parts, k+v)
		} else {
			parts = append(parts, k+"="+v)
		}
	}
	if len(parts) == 0 {
		return "gpu_name=Tesla_V100 num_gpus=1 disk_space>=20"
	}
	return strings.Join(parts, " ")
}

// searchSpec is the resolved Vast query plus client-side post-filters.
type searchSpec struct {
	query       map[string]any
	minGPURamMB int64
	// Post-filters (client-side).
	enforceGPUWhitelist bool
	enforceRegion       bool
	regions             []string
	cpuNameContains     string
	minDiskBW           float64 // MB/s; 0 = no filter
	minCPUCores         float64
	minCPUGHz           float64
	minCPURamGB         float64
	minPCIeBW           float64
	minGPUMemBW         float64 // GB/s
	minCPUMemBW         float64 // GB/s theoretical host DRAM (from local table)
	minInetDown         float64
	maxDPH              float64 // 0 = no filter
	limit               int
}

// resolveSearchSpec parses a filter string into a Vast API query and post-filters.
//
// Supported tokens (field op value), ops: >= <= != = > <
//
// GPU / price / disk / network (forwarded to Vast when recognized):
//
//	gpu_name=Tesla_V100  num_gpus=1  gpu_ram>=16000 (MB in API)
//	dph_total<0.5  disk_bw>400  cpu_cores_effective>=8  cpu_ghz>3
//	cpu_ram>=32000 (MB)  has_avx=1  cpu_arch=amd64  pcie_bw>10
//	gpu_mem_bw>400  inet_down>500  geolocation=US
//
// Convenience aliases (local):
//
//	min_ram=32           → min GPU VRAM GB (post-filter)
//	max_dph=0.5          → dph_total <= 0.5
//	max_price=0.5        → same
//	min_disk_bw=400      → disk_bw >= 400 (MB/s)  ← closest to "memory bandwidth MB/s"
//	min_cpu_cores=8      → cpu_cores_effective >= 8
//	min_cpu_ghz=3        → cpu_ghz >= 3
//	min_cpu_ram=32       → cpu_ram >= 32 GB (converted to MB for API)
//	min_pcie_bw=10       → pcie_bw >= 10
//	min_gpu_mem_bw=400   → gpu_mem_bw >= 400 (GB/s GPU VRAM)
//	min_cpu_mem_bw=100   → host DRAM theoretical peak >= 100 GB/s (local CPU table)
//	min_dram_bw=100      → alias for min_cpu_mem_bw
//	min_mem_bw=100       → alias for min_cpu_mem_bw
//	min_inet_down=500    → inet_down >= 500
//	cpu_name=EPYC        → post-filter substring (API cannot filter cpu_name)
//	cpu_contains=EPYC    → same
//	any_gpu=1            → skip GPU allowlist (needed for cheap consumer GPUs / CPU work)
//	whitelist=0          → same as any_gpu=1
//	region=any           → skip US/CA region lock
//	region=US,CA,GB      → restrict to these country codes
//	limit=20             → cap returned offers after filter
//
// Defaults for provision path: verified/rentable, order by dph_total asc, type on-demand.
// GPU whitelist + NA region are ON by default (set any_gpu=1 / region=any to disable).
func resolveSearchSpec(filter string) searchSpec {
	spec := searchSpec{
		query: map[string]any{
			"verified": map[string]any{"eq": true},
			"external": map[string]any{"eq": false},
			"rented":   map[string]any{"eq": false},
			// rentable=true is REQUIRED: without it /bundles/ returns catalog machines
			// that are not currently rentable, and every create then fails with
			// "404/3603 no_such_ask ... not available". Matches vast-python's default.
			"rentable": map[string]any{"eq": true},
			"order":    [][]string{{"dph_total", "asc"}},
			"type":     "on-demand",
		},
		enforceGPUWhitelist: true,
		enforceRegion:       true,
		regions:             append([]string(nil), allowedRegions...),
		limit:               50,
	}

	var requestedGPU string

	for _, token := range strings.Fields(filter) {
		// Aliases without operators first.
		if strings.HasPrefix(token, "min_ram=") {
			val := token[len("min_ram="):]
			var gb int64
			fmt.Sscanf(val, "%d", &gb)
			spec.minGPURamMB = gb * 1024
			continue
		}
		if strings.HasPrefix(token, "any_gpu=") {
			v := strings.ToLower(token[len("any_gpu="):])
			if v == "1" || v == "true" || v == "yes" {
				spec.enforceGPUWhitelist = false
			}
			continue
		}
		if strings.HasPrefix(token, "whitelist=") {
			v := strings.ToLower(token[len("whitelist="):])
			if v == "0" || v == "false" || v == "no" || v == "off" {
				spec.enforceGPUWhitelist = false
			}
			continue
		}
		if strings.HasPrefix(token, "region=") {
			val := token[len("region="):]
			if strings.EqualFold(val, "any") || val == "*" {
				spec.enforceRegion = false
			} else {
				spec.enforceRegion = true
				spec.regions = nil
				for _, r := range strings.Split(val, ",") {
					r = strings.TrimSpace(r)
					if r != "" {
						spec.regions = append(spec.regions, strings.ToUpper(r))
					}
				}
			}
			continue
		}
		if strings.HasPrefix(token, "cpu_name=") || strings.HasPrefix(token, "cpu_contains=") {
			if i := strings.Index(token, "="); i > 0 {
				spec.cpuNameContains = strings.ReplaceAll(token[i+1:], "_", " ")
			}
			continue
		}
		if strings.HasPrefix(token, "max_dph=") || strings.HasPrefix(token, "max_price=") {
			i := strings.Index(token, "=")
			if f, err := strconv.ParseFloat(token[i+1:], 64); err == nil {
				spec.maxDPH = f
				spec.query["dph_total"] = map[string]any{"lte": f}
			}
			continue
		}
		if strings.HasPrefix(token, "min_disk_bw=") {
			if f, err := strconv.ParseFloat(token[len("min_disk_bw="):], 64); err == nil {
				spec.minDiskBW = f
				spec.query["disk_bw"] = map[string]any{"gte": f}
			}
			continue
		}
		if strings.HasPrefix(token, "min_cpu_cores=") {
			if f, err := strconv.ParseFloat(token[len("min_cpu_cores="):], 64); err == nil {
				spec.minCPUCores = f
				spec.query["cpu_cores_effective"] = map[string]any{"gte": f}
			}
			continue
		}
		if strings.HasPrefix(token, "min_cpu_ghz=") {
			if f, err := strconv.ParseFloat(token[len("min_cpu_ghz="):], 64); err == nil {
				spec.minCPUGHz = f
				spec.query["cpu_ghz"] = map[string]any{"gte": f}
			}
			continue
		}
		if strings.HasPrefix(token, "min_cpu_ram=") {
			// Alias takes GB; Vast API cpu_ram is MB.
			if f, err := strconv.ParseFloat(token[len("min_cpu_ram="):], 64); err == nil {
				spec.minCPURamGB = f
				spec.query["cpu_ram"] = map[string]any{"gte": f * 1024}
			}
			continue
		}
		if strings.HasPrefix(token, "min_pcie_bw=") {
			if f, err := strconv.ParseFloat(token[len("min_pcie_bw="):], 64); err == nil {
				spec.minPCIeBW = f
				spec.query["pcie_bw"] = map[string]any{"gte": f}
			}
			continue
		}
		if strings.HasPrefix(token, "min_gpu_mem_bw=") {
			if f, err := strconv.ParseFloat(token[len("min_gpu_mem_bw="):], 64); err == nil {
				spec.minGPUMemBW = f
				spec.query["gpu_mem_bw"] = map[string]any{"gte": f}
			}
			continue
		}
		// System DRAM bandwidth from local CPU model table (GB/s host peak).
		if strings.HasPrefix(token, "min_cpu_mem_bw=") ||
			strings.HasPrefix(token, "min_dram_bw=") ||
			strings.HasPrefix(token, "min_mem_bw=") {
			i := strings.Index(token, "=")
			if f, err := strconv.ParseFloat(token[i+1:], 64); err == nil {
				// Accept accidental MB/s inputs: values > 2000 are treated as MB/s → GB/s.
				if f > 2000 {
					f = f / 1000.0
				}
				spec.minCPUMemBW = f
			}
			continue
		}
		if strings.HasPrefix(token, "min_inet_down=") {
			if f, err := strconv.ParseFloat(token[len("min_inet_down="):], 64); err == nil {
				spec.minInetDown = f
				spec.query["inet_down"] = map[string]any{"gte": f}
			}
			continue
		}
		if strings.HasPrefix(token, "limit=") {
			if n, err := strconv.Atoi(token[len("limit="):]); err == nil && n > 0 {
				spec.limit = n
			}
			continue
		}

		// Generic field op value.
		matched := false
		for _, op := range []struct {
			sym string
			api string
		}{
			{">=", "gte"}, {"<=", "lte"}, {"!=", "neq"}, {"=", "eq"}, {">", "gt"}, {"<", "lt"},
		} {
			idx := strings.Index(token, op.sym)
			if idx <= 0 {
				continue
			}
			key := token[:idx]
			valStr := token[idx+len(op.sym):]
			if key == "gpu_name" {
				valStr = strings.ReplaceAll(valStr, "_", " ")
				requestedGPU = valStr
			}
			// cpu_name is not a Vast filter — treat as contains post-filter.
			if key == "cpu_name" || key == "cpu_contains" {
				spec.cpuNameContains = strings.ReplaceAll(valStr, "_", " ")
				matched = true
				break
			}
			if localOnlyKeys[key] {
				matched = true
				break
			}
			if !searchFilterKeys[key] && key != "order" && key != "type" && key != "limit" {
				// Still forward unknown keys — Vast may accept new ones.
			}
			spec.query[key] = map[string]any{op.api: coerceValue(valStr)}
			matched = true
			break
		}
		if !matched {
			fmt.Fprintf(os.Stderr, "scp-runner: ignoring unparsed filter token %q\n", token)
		}
	}

	// If a specific GPU was requested, look up its min RAM from the allowlist.
	if requestedGPU != "" && spec.minGPURamMB == 0 {
		for _, g := range allowedGPUs {
			if g.Name == requestedGPU {
				spec.minGPURamMB = g.MinRAM
				break
			}
		}
	}

	return spec
}

// resolveGPUSpec is the provision-path wrapper (keeps older call sites working).
func resolveGPUSpec(filter string) (query map[string]any, minRAM int64) {
	spec := resolveSearchSpec(filter)
	return spec.query, spec.minGPURamMB
}

// coerceValue converts a token value to bool/float/int when possible.
func coerceValue(s string) any {
	switch strings.ToLower(s) {
	case "true":
		return true
	case "false":
		return false
	}
	if i, err := strconv.ParseInt(s, 10, 64); err == nil {
		return i
	}
	if f, err := strconv.ParseFloat(s, 64); err == nil {
		return f
	}
	// Underscores → spaces for string enums (gpu names etc. already handled).
	return strings.ReplaceAll(s, "_", " ")
}

// isAllowedGPU checks if a GPU name is in the allowedGPUs whitelist.
func isAllowedGPU(name string) bool {
	for _, g := range allowedGPUs {
		if g.Name == name {
			return true
		}
	}
	return false
}

// isAllowedRegion checks if a geolocation string ends with an allowed country code.
func isAllowedRegion(geo string, regions []string) bool {
	if len(regions) == 0 {
		return true
	}
	for _, r := range regions {
		if strings.HasSuffix(geo, ", "+r) || strings.HasSuffix(geo, ","+r) ||
			strings.HasSuffix(geo, "_"+r) || strings.EqualFold(geo, r) {
			return true
		}
		// "California, US" / "Wisconsin, US" / bare ", US"
		if strings.HasSuffix(strings.ToUpper(geo), ", "+strings.ToUpper(r)) {
			return true
		}
	}
	// Empty geolocation: some hosts omit country; treat as fail when enforcing.
	return false
}

// filterOffers applies whitelist / region / VRAM / local thresholds.
func filterOffers(offers []VastOffer, spec searchSpec) []VastOffer {
	var filtered []VastOffer
	for _, o := range offers {
		if spec.enforceGPUWhitelist && !isAllowedGPU(o.GPUName) {
			continue
		}
		if spec.minGPURamMB > 0 && o.GPUMemMB < spec.minGPURamMB {
			continue
		}
		if spec.enforceRegion && !isAllowedRegion(o.Geolocation, spec.regions) {
			continue
		}
		if spec.maxDPH > 0 && o.DPHTot > spec.maxDPH {
			continue
		}
		if spec.minDiskBW > 0 && o.DiskBW < spec.minDiskBW {
			continue
		}
		if spec.minCPUCores > 0 && o.CPUCoresEffective < spec.minCPUCores {
			continue
		}
		if spec.minCPUGHz > 0 && o.CPUGHz < spec.minCPUGHz {
			continue
		}
		if spec.minCPURamGB > 0 {
			ramGB := o.CPURamMB
			if ramGB > 512 { // MB
				ramGB = ramGB / 1024
			}
			if ramGB < spec.minCPURamGB {
				continue
			}
		}
		if spec.minPCIeBW > 0 && o.PCIeBW < spec.minPCIeBW {
			continue
		}
		if spec.minGPUMemBW > 0 && o.GPUMemBW < spec.minGPUMemBW {
			continue
		}
		if spec.minInetDown > 0 && o.InetDown < spec.minInetDown {
			continue
		}
		if spec.cpuNameContains != "" &&
			!strings.Contains(strings.ToUpper(o.CPUName), strings.ToUpper(spec.cpuNameContains)) {
			continue
		}
		annotateCPU(&o)
		if spec.minCPUMemBW > 0 {
			// Prefer host-level peak (accounts for dual-socket estimate).
			hostBW := o.CPUMemHostGBs
			if hostBW <= 0 {
				hostBW = o.CPUMemPeakGBs
			}
			if hostBW < spec.minCPUMemBW {
				continue
			}
		}
		filtered = append(filtered, o)
	}
	if spec.limit > 0 && len(filtered) > spec.limit {
		filtered = filtered[:spec.limit]
	}
	return filtered
}

// SearchOffers finds available GPU offers matching a filter query.
// Results are post-filtered against the allowedGPUs whitelist, minimum VRAM,
// and North America region restriction (unless any_gpu=1 / region=any).
func (v *VastClient) SearchOffers(ctx context.Context, filter string) ([]VastOffer, error) {
	spec := resolveSearchSpec(filter)
	body, err := json.Marshal(spec.query)
	if err != nil {
		return nil, fmt.Errorf("marshal query: %w", err)
	}

	data, err := v.doRequest(ctx, "POST", "/bundles/", strings.NewReader(string(body)))
	if err != nil {
		return nil, fmt.Errorf("search offers: %w", err)
	}

	var result struct {
		Offers []VastOffer `json:"offers"`
	}
	if err := json.Unmarshal(data, &result); err != nil {
		return nil, fmt.Errorf("parse offers: %w", err)
	}

	// Annotate even pre-filter for debug; filterOffers re-annotates kept ones.
	filtered := filterOffers(result.Offers, spec)
	if len(filtered) == 0 && len(result.Offers) > 0 {
		fmt.Fprintf(os.Stderr, "scp-runner: %d offers found but none passed filters (whitelist=%v min_vram=%d MB region=%v max_dph=%g min_disk_bw=%g)\n",
			len(result.Offers), spec.enforceGPUWhitelist, spec.minGPURamMB, spec.regions, spec.maxDPH, spec.minDiskBW)
		for i, o := range result.Offers {
			if i >= 5 {
				fmt.Fprintf(os.Stderr, "  ... and %d more\n", len(result.Offers)-5)
				break
			}
			fmt.Fprintf(os.Stderr, "  rejected: %s cpu=%q disk_bw=%.0f $%.3f/hr %s\n",
				o.GPUName, strings.TrimSpace(o.CPUName), o.DiskBW, o.DPHTot, o.Geolocation)
		}
	}
	return filtered, nil
}

// CreateInstanceResult holds the response from creating a Vast.ai instance.
type CreateInstanceResult struct {
	Success     bool   `json:"success"`
	NewContract int    `json:"new_contract"`
	Error       string `json:"error,omitempty"`
	Msg         string `json:"msg,omitempty"`
}

// CreateInstance rents a GPU instance. Returns the instance/contract ID.
func (v *VastClient) CreateInstance(ctx context.Context, offerID int, image string, diskGB int, onstart string) (int, error) {
	// Match vast-python's create payload: target_state/runtype/client_id ensure the
	// rented instance actually boots an SSH server and transitions to running.
	payload := map[string]any{
		"client_id":    "me",
		"image":        image,
		"disk":         float64(diskGB),
		"onstart":      onstart,
		"runtype":      "ssh",
		"target_state": "running",
	}
	body, err := json.Marshal(payload)
	if err != nil {
		return 0, fmt.Errorf("marshal: %w", err)
	}

	data, err := v.doRequest(ctx, "PUT", fmt.Sprintf("/asks/%d/", offerID), strings.NewReader(string(body)))
	if err != nil {
		return 0, fmt.Errorf("create instance: %w", err)
	}

	var result CreateInstanceResult
	if err := json.Unmarshal(data, &result); err != nil {
		return 0, fmt.Errorf("parse response: %w: %s", err, string(data))
	}
	if !result.Success || result.NewContract == 0 {
		return 0, fmt.Errorf("create failed: %s %s", result.Error, result.Msg)
	}
	return result.NewContract, nil
}

// ShowInstances lists all running instances.
func (v *VastClient) ShowInstances(ctx context.Context) ([]VastInstance, error) {
	data, err := v.doRequest(ctx, "GET", "/instances/", nil)
	if err != nil {
		return nil, fmt.Errorf("show instances: %w", err)
	}

	var result struct {
		Instances []VastInstance `json:"instances"`
	}
	if err := json.Unmarshal(data, &result); err != nil {
		// Try as direct array.
		var instances []VastInstance
		if err2 := json.Unmarshal(data, &instances); err2 != nil {
			return nil, fmt.Errorf("parse instances: %w", err)
		}
		return instances, nil
	}
	return result.Instances, nil
}

// DestroyInstance terminates a running instance.
func (v *VastClient) DestroyInstance(ctx context.Context, instanceID int) error {
	_, err := v.doRequest(ctx, "DELETE", fmt.Sprintf("/instances/%d/", instanceID), nil)
	if err != nil {
		return fmt.Errorf("destroy instance %d: %w", instanceID, err)
	}
	return nil
}

// FormatOfferLine is a one-line human summary for CLI/logs.
func FormatOfferLine(o VastOffer) string {
	ramGB := o.CPURamMB
	if ramGB > 512 {
		ramGB = ramGB / 1024
	}
	dram := "?"
	if o.CPUMemHostGBs > 0 {
		dram = fmt.Sprintf("%.0fGB/s", o.CPUMemHostGBs)
	}
	return fmt.Sprintf("$%.4f/hr  dram=%s  disk_bw=%.0fMB/s  vCPU=%.1f@%.1fGHz  RAM=%.0fGB  pcie=%.1f  gpu_mem=%.0fGB/s  cpu=%q [%s]  gpu=%s×%d  %s  id=%d mach=%d",
		o.DPHTot, dram, o.DiskBW, o.CPUCoresEffective, o.CPUGHz, ramGB, o.PCIeBW, o.GPUMemBW,
		strings.TrimSpace(o.CPUName), o.CPUClass, o.GPUName, o.NumGPUs, o.Geolocation, o.ID, o.MachineID)
}
