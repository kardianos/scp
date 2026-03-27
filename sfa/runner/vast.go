package main

import (
	"context"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"strings"
	"time"
)

const vastAPIBase = "https://console.vast.ai/api/v0"

// VastClient talks to the Vast.ai REST API.
type VastClient struct {
	apiKey     string
	httpClient *http.Client
}

// VastOffer is a GPU rental offer.
type VastOffer struct {
	ID          int     `json:"id"`
	GPUName     string  `json:"gpu_name"`
	NumGPUs     int     `json:"num_gpus"`
	DPHTot      float64 `json:"dph_total"`
	DiskSpace   float64 `json:"disk_space"`
	GPUMemMB    int64   `json:"gpu_ram"`
	Reliability float64 `json:"reliability2"`
	InetDown    float64 `json:"inet_down"`
	Geolocation string  `json:"geolocation"`
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

// VastInstance is a running instance.
type VastInstance struct {
	ID         int    `json:"id"`
	Status     string `json:"actual_status"`
	SSHHost    string `json:"ssh_host"`
	SSHPort    int    `json:"ssh_port"`
	GPUName    string `json:"gpu_name"`
	CurState   string `json:"cur_state"`
	ContractID int    `json:"contract_id,omitempty"`
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

// resolveGPUSpec parses a gpu_filter string and returns the API query, the
// minimum VRAM requirement (in MB), and whether to enforce region filtering.
//
// The gpu_filter can be:
//   - A gpu_name filter: "gpu_name=Tesla_V100" — looked up in allowedGPUs
//   - A min_ram filter:  "min_ram=32" — picks cheapest allowed GPU with >= 32 GB
//   - Raw Vast.ai tokens for backwards compat (but still post-filtered)
//
// All queries enforce region=North America and the allowedGPUs whitelist.
func resolveGPUSpec(filter string) (query map[string]any, minRAM int64) {
	query = map[string]any{
		"verified": map[string]any{"eq": true},
		"external": map[string]any{"eq": false},
		"rented":   map[string]any{"eq": false},
		"order":    [][]string{{"dph_total", "asc"}},
		"type":     "on-demand",
	}

	var requestedGPU string

	for _, token := range strings.Fields(filter) {
		// Handle min_ram=XX (in GB) as a VRAM floor.
		if strings.HasPrefix(token, "min_ram=") {
			val := token[len("min_ram="):]
			var gb int64
			fmt.Sscanf(val, "%d", &gb)
			minRAM = gb * 1024
			continue
		}

		// Try operators in order: >=, <=, !=, =, >, <
		for _, op := range []struct {
			sym string
			api string
		}{
			{">=", "gte"}, {"<=", "lte"}, {"!=", "neq"}, {"=", "eq"}, {">", "gt"}, {"<", "lt"},
		} {
			if idx := strings.Index(token, op.sym); idx > 0 {
				key := token[:idx]
				val := token[idx+len(op.sym):]
				if key == "gpu_name" {
					val = strings.ReplaceAll(val, "_", " ")
					requestedGPU = val
				}
				query[key] = map[string]any{op.api: val}
				break
			}
		}
	}

	// If a specific GPU was requested, look up its min RAM from the allowlist.
	if requestedGPU != "" && minRAM == 0 {
		for _, g := range allowedGPUs {
			if g.Name == requestedGPU {
				minRAM = g.MinRAM
				break
			}
		}
	}

	return query, minRAM
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
func isAllowedRegion(geo string) bool {
	for _, r := range allowedRegions {
		if strings.HasSuffix(geo, ", "+r) || strings.HasSuffix(geo, ","+r) {
			return true
		}
	}
	return false
}

// filterOffers applies the allowedGPUs whitelist, minimum VRAM, and region
// restrictions to a list of offers.
func filterOffers(offers []VastOffer, minRAM int64) []VastOffer {
	var filtered []VastOffer
	for _, o := range offers {
		if !isAllowedGPU(o.GPUName) {
			continue
		}
		if minRAM > 0 && o.GPUMemMB < minRAM {
			continue
		}
		if !isAllowedRegion(o.Geolocation) {
			continue
		}
		filtered = append(filtered, o)
	}
	return filtered
}

// SearchOffers finds available GPU offers matching a filter query.
// Results are post-filtered against the allowedGPUs whitelist, minimum VRAM,
// and North America region restriction.
func (v *VastClient) SearchOffers(ctx context.Context, filter string) ([]VastOffer, error) {
	query, minRAM := resolveGPUSpec(filter)
	body, err := json.Marshal(query)
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

	filtered := filterOffers(result.Offers, minRAM)
	if len(filtered) == 0 && len(result.Offers) > 0 {
		// Log what we rejected to help debug.
		fmt.Fprintf(os.Stderr, "scp-runner: %d offers found but none passed filters (gpu whitelist, min_ram=%d MB, region=NA)\n",
			len(result.Offers), minRAM)
		for i, o := range result.Offers {
			if i >= 5 {
				fmt.Fprintf(os.Stderr, "  ... and %d more\n", len(result.Offers)-5)
				break
			}
			fmt.Fprintf(os.Stderr, "  rejected: %s (%d MB) %s $%.3f/hr\n",
				o.GPUName, o.GPUMemMB, o.Geolocation, o.DPHTot)
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
	payload := map[string]any{
		"client_id": "me",
		"image":     image,
		"disk":      diskGB,
		"onstart":   onstart,
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
