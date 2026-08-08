package main

import (
	"encoding/json"
	"fmt"
	"reflect"
	"strings"
)

// --- Tool parameter and result structs ---

type SimSetupParams struct {
	Name      string            `json:"name" desc:"Instance name/ID (e.g. 'gpu1', 'local_test')" required:"true"`
	Executor  string            `json:"executor" desc:"Execution type: local or remote" required:"true"`
	Host      string            `json:"host" desc:"SSH host for existing remote instance"`
	Port      int               `json:"port" desc:"SSH port for existing remote instance"`
	GPUFilter map[string]string `json:"gpu_filter" desc:"Offer filter key=value (ops: >=,<=,!=,=,>,<). GPU: gpu_name, min_ram(GB), num_gpus. CPU: min_cpu_cores, min_cpu_ghz, min_cpu_ram(GB), cpu_name, has_avx. DRAM: min_cpu_mem_bw(GB/s theoretical host peak from local CPU table). Disk/GPU BW: min_disk_bw(MB/s), min_pcie_bw, min_gpu_mem_bw(GB/s). Price: max_dph. Escape: any_gpu=1, region=any|US,CA. Default GPU allowlist: Tesla_V100,A100_*,L40S,H100_*,B200,RTX_PRO_*"`
	DiskGB    int               `json:"disk_gb" desc:"Disk space in GB to provision (required for remote)" required:"true"`
	WorkDir   string            `json:"work_dir" desc:"Working directory"`
}

// SimSearchOffersParams searches Vast.ai without provisioning.
// Prefer structured fields; filter_str accepts the same token language as gpu_filter / CLI.
type SimSearchOffersParams struct {
	// FilterStr is a whitespace-separated query, e.g.
	// "max_dph=0.5 min_disk_bw=400 min_cpu_cores=8 any_gpu=1 region=US,CA"
	FilterStr string `json:"filter" desc:"Whitespace-separated filter tokens (same language as gpu_filter / CLI search). Example: max_dph=0.5 min_disk_bw=400 min_cpu_cores=8 any_gpu=1"`
	// Structured knobs (merged into filter; override filter_str on conflict for the same key).
	MaxDPH       float64 `json:"max_dph" desc:"Max $/hour (dph_total)"`
	MinDiskBW    float64 `json:"min_disk_bw" desc:"Min disk read bandwidth MB/s (Vast disk_bw; closest API field to 'memory bandwidth in MB/s')"`
	MinCPUCores  float64 `json:"min_cpu_cores" desc:"Min effective vCPUs"`
	MinCPUGHz    float64 `json:"min_cpu_ghz" desc:"Min CPU clock GHz"`
	MinCPURamGB  float64 `json:"min_cpu_ram_gb" desc:"Min system RAM in GB"`
	MinPCIeBW    float64 `json:"min_pcie_bw" desc:"Min PCIe bandwidth (CPU↔GPU)"`
	MinGPUMemBW  float64 `json:"min_gpu_mem_bw" desc:"Min GPU memory bandwidth GB/s"`
	MinCPUMemBW  float64 `json:"min_cpu_mem_bw" desc:"Min host DRAM theoretical peak GB/s from local CPU model table (channels×MT/s×8/1000 × sockets). Aliases: min_dram_bw, min_mem_bw. Overrides in ~/.scp-runner/cpu_mem_bw.json"`
	CPUContains  string  `json:"cpu_contains" desc:"Substring match on cpu_name (client-side; Vast API cannot filter cpu_name)"`
	GPUName      string  `json:"gpu_name" desc:"GPU model (underscores ok: RTX_3090)"`
	MinGPURamGB  float64 `json:"min_gpu_ram_gb" desc:"Min GPU VRAM in GB"`
	AnyGPU       bool    `json:"any_gpu" desc:"If true, skip the production GPU allowlist (needed for cheap consumer GPUs / CPU-oriented hosts)"`
	Region       string  `json:"region" desc:"Country codes comma-separated (default US,CA) or 'any'"`
	Limit        int     `json:"limit" desc:"Max offers to return (default 20)"`
}

type SimSearchOffersResult struct {
	Count  int              `json:"count"`
	Filter string           `json:"filter_used"`
	Note   string           `json:"note,omitempty"`
	Offers []VastOfferBrief `json:"offers"`
}

// VastOfferBrief is a compact offer row for MCP/CLI output.
type VastOfferBrief struct {
	ID                int     `json:"id"`
	MachineID         int     `json:"machine_id"`
	DPHTotal          float64 `json:"dph_total"`
	GPUName           string  `json:"gpu_name"`
	NumGPUs           int     `json:"num_gpus"`
	GPUMemMB          int64   `json:"gpu_ram_mb"`
	GPUMemBW          float64 `json:"gpu_mem_bw_gbs"`
	CPUName           string  `json:"cpu_name"`
	CPUClass          string  `json:"cpu_class,omitempty"`
	CPUNote           string  `json:"cpu_note,omitempty"`
	CPUMemPeakGBs     float64 `json:"cpu_mem_peak_gbs,omitempty"`  // theoretical / socket
	CPUMemHostGBs     float64 `json:"cpu_mem_host_gbs,omitempty"`  // × estimated sockets
	CPUMemStreamGBs   float64 `json:"cpu_mem_stream_gbs,omitempty"` // ~0.78 × host peak
	CPUMemChannels    int     `json:"cpu_mem_channels,omitempty"`
	CPUMemMTs         int     `json:"cpu_mem_mts,omitempty"`
	CPUMemTech        string  `json:"cpu_mem_tech,omitempty"`
	CPUMemConf        string  `json:"cpu_mem_confidence,omitempty"`
	CPUCoresEffective float64 `json:"cpu_cores_effective"`
	CPUGHz            float64 `json:"cpu_ghz"`
	CPURamGB          float64 `json:"cpu_ram_gb"`
	DiskBW            float64 `json:"disk_bw_mbs"`
	PCIeBW            float64 `json:"pcie_bw"`
	InetDown          float64 `json:"inet_down"`
	Geolocation       string  `json:"geolocation"`
	Reliability       float64 `json:"reliability"`
	DiskSpace         float64 `json:"disk_space_gb"`
}

type SimSetupResult struct {
	Status     string  `json:"status"`
	Name       string  `json:"name"`
	Type       string  `json:"type"`
	WorkDir    string  `json:"work_dir,omitempty"`
	Host       string  `json:"host,omitempty"`
	Port       int     `json:"port,omitempty"`
	GPUName    string  `json:"gpu_name,omitempty"`
	InstanceID int     `json:"instance_id,omitempty"`
	MachineID  int     `json:"machine_id,omitempty"`
	DPHTotal   float64 `json:"dph_total,omitempty"`
}

type SimStatusParams struct {
	Name string `json:"name" desc:"Instance name (omit to list all instances)"`
}

type SimTeardownParams struct {
	Name  string `json:"name" desc:"Instance name to teardown" required:"true"`
	Force bool   `json:"force" desc:"Force teardown even if output files have not been downloaded"`
}
type SimTeardownResult struct {
	Status string `json:"status"`
}

type SimBuildParams struct {
	Name    string   `json:"name" desc:"Instance name" required:"true"`
	Sources []string `json:"sources" desc:"Source file paths" required:"true"`
	Cmd     string   `json:"cmd" desc:"Build command (optional, auto-detected if omitted)"`
}

type SimRunParams struct {
	Name           string  `json:"name" desc:"Instance name" required:"true"`
	Config         string  `json:"config" desc:"Simulation config content (inline)" required:"true"`
	ID             string  `json:"id" desc:"Run identifier" required:"true"`
	NotifyInterval float64 `json:"notify_interval" desc:"Send progress notifications every N seconds (0=disabled)"`
	Wait           bool    `json:"wait" desc:"Block until run completes and return final status (default: false)"`
	AutoDownload   string  `json:"auto_download" desc:"Local directory for automatic incremental downloads during the run"`
	OnComplete     string  `json:"on_complete" desc:"Local shell command run (bash -c) when the run reaches a terminal state; env: SCP_RUN_ID, SCP_STATE, SCP_INSTANCE"`
	NoQueue        bool    `json:"no_queue" desc:"Local runs only: bypass the nproc-aware concurrency queue and start immediately"`
}

type SimRunResult struct {
	RunID    string `json:"run_id"`
	Status   string `json:"status"`
	Instance string `json:"instance"`
	Executor string `json:"executor"`
	Init     string `json:"init,omitempty"`
	InitSFA  string `json:"init_sfa,omitempty"`
}

type SimListTemplatesParams struct{}

type SimListTemplatesResult struct {
	Templates []TemplateInfo `json:"templates"`
}

type TemplateInfo struct {
	Path    string `json:"path"`
	Name    string `json:"name"`
	Size    int64  `json:"size"`
	ModTime string `json:"mod_time"`
}

type SimRunStatusParams struct {
	Name string `json:"name" desc:"Instance name" required:"true"`
	ID   string `json:"id" desc:"Run identifier" required:"true"`
}

type SimRunCancelParams struct {
	Name string `json:"name" desc:"Instance name" required:"true"`
	ID   string `json:"id" desc:"Run identifier" required:"true"`
}

type SimRunCancelResult struct {
	Status string `json:"status"`
}

type SimUploadParams struct {
	Name       string `json:"name" desc:"Instance name" required:"true"`
	LocalPath  string `json:"local_path" desc:"Local file path" required:"true"`
	RemotePath string `json:"remote_path" desc:"Remote destination path" required:"true"`
}

type SimUploadResult struct {
	Status string `json:"status"`
}

type SimDownloadParams struct {
	Name       string `json:"name" desc:"Instance name" required:"true"`
	RemotePath string `json:"remote_path" desc:"Remote file path" required:"true"`
	LocalPath  string `json:"local_path" desc:"Local destination path" required:"true"`
	Wait       bool   `json:"wait" desc:"Block until the transfer completes, verify local size == remote size, and return {remote_bytes, local_bytes, verified} (default: false, async)"`
}

type SimDownloadResult struct {
	DownloadID  string `json:"download_id"`
	Status      string `json:"status"`
	RemoteBytes int64  `json:"remote_bytes,omitempty"`
	LocalBytes  int64  `json:"local_bytes,omitempty"`
	Verified    bool   `json:"verified,omitempty"`
	Error       string `json:"error,omitempty"`
}

type SimDownloadStatusParams struct {
	Name string `json:"name" desc:"Instance name" required:"true"`
	ID   string `json:"id" desc:"Download identifier" required:"true"`
}

type SimListFilesParams struct {
	Name    string `json:"name" desc:"Instance name" required:"true"`
	Pattern string `json:"pattern" desc:"File glob pattern" required:"true"`
}

type SimListFilesResult struct {
	Files []FileInfo `json:"files"`
}

type SimExecParams struct {
	Name      string `json:"name" desc:"Instance name" required:"true"`
	Cmd       string `json:"cmd" desc:"Shell command to execute" required:"true"`
	TimeoutMs int    `json:"timeout_ms" desc:"Timeout in milliseconds (default 30000)"`
}

type SimExecResult struct {
	Output string `json:"output"`
	Status string `json:"status"`
	Error  string `json:"error,omitempty"`
}

// --- Reflection-based JSON Schema generation ---

func structToJSONSchema(t reflect.Type) map[string]any {
	if t.Kind() == reflect.Pointer {
		t = t.Elem()
	}

	props := map[string]any{}
	var required []string

	for i := 0; i < t.NumField(); i++ {
		field := t.Field(i)
		jsonTag := field.Tag.Get("json")
		name := strings.Split(jsonTag, ",")[0]
		if name == "" || name == "-" {
			continue
		}

		prop := map[string]any{}

		switch field.Type.Kind() {
		case reflect.String:
			prop["type"] = "string"
		case reflect.Int, reflect.Int64, reflect.Float64:
			prop["type"] = "number"
		case reflect.Bool:
			prop["type"] = "boolean"
		case reflect.Slice:
			prop["type"] = "array"
			if field.Type.Elem().Kind() == reflect.String {
				prop["items"] = map[string]any{"type": "string"}
			}
		default:
			prop["type"] = "object"
		}

		if desc := field.Tag.Get("desc"); desc != "" {
			prop["description"] = desc
		}
		if field.Tag.Get("required") == "true" {
			required = append(required, name)
		}

		props[name] = prop
	}

	schema := map[string]any{
		"type":       "object",
		"properties": props,
	}
	if len(required) > 0 {
		schema["required"] = required
	}
	return schema
}

// unmarshalParams converts a map[string]any into a typed struct via JSON round-trip.
// It pre-processes string values that look like JSON arrays or objects, which can
// happen when MCP clients serialize complex types as strings.
func unmarshalParams(params map[string]any, target any) error {
	// Pre-process: if a value is a string that looks like a JSON array, object,
	// boolean, or number, decode it so the round-trip produces the correct Go types.
	// This handles MCP clients that serialize complex types as strings.
	for k, v := range params {
		if s, ok := v.(string); ok {
			s = strings.TrimSpace(s)
			if (strings.HasPrefix(s, "[") && strings.HasSuffix(s, "]")) ||
				(strings.HasPrefix(s, "{") && strings.HasSuffix(s, "}")) {
				var decoded any
				if err := json.Unmarshal([]byte(s), &decoded); err == nil {
					params[k] = decoded
				}
			} else if s == "true" {
				params[k] = true
			} else if s == "false" {
				params[k] = false
			}
		}
	}

	data, err := json.Marshal(params)
	if err != nil {
		return fmt.Errorf("marshal params: %w", err)
	}
	return json.Unmarshal(data, target)
}
