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
}

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

// parseFilter converts a CLI-style filter string like "gpu_name=Tesla_V100 num_gpus=1 disk_space>=20"
// into the Vast.ai API JSON query format: {"gpu_name": {"eq": "Tesla V100"}, "num_gpus": {"eq": "1"}, ...}
func parseFilter(filter string) map[string]any {
	query := map[string]any{
		"verified": map[string]any{"eq": true},
		"external": map[string]any{"eq": false},
		"rented":   map[string]any{"eq": false},
		"order":    [][]string{{"dph_total", "asc"}},
		"type":     "on-demand",
	}

	for _, token := range strings.Fields(filter) {
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
				// Replace underscores with spaces in GPU names.
				if key == "gpu_name" {
					val = strings.ReplaceAll(val, "_", " ")
				}
				query[key] = map[string]any{op.api: val}
				break
			}
		}
	}
	return query
}

// SearchOffers finds available GPU offers matching a filter query.
func (v *VastClient) SearchOffers(ctx context.Context, filter string) ([]VastOffer, error) {
	query := parseFilter(filter)
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
	return result.Offers, nil
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
