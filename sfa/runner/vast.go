package main

import (
	"context"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"net/url"
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
	keyPath := filepath.Join(home, ".vast_api_key")
	data, err := os.ReadFile(keyPath)
	if err != nil {
		// Fall back to env var.
		key := os.Getenv("VAST_API_KEY")
		if key == "" {
			return nil, fmt.Errorf("read %s: %w (and VAST_API_KEY not set)", keyPath, err)
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

// SearchOffers finds available GPU offers matching a filter query.
func (v *VastClient) SearchOffers(ctx context.Context, filter string) ([]VastOffer, error) {
	q := url.Values{}
	q.Set("q", filter)
	q.Set("order", "dph_total")
	q.Set("type", "on-demand")

	data, err := v.doRequest(ctx, "GET", "/bundles/?"+q.Encode(), nil)
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

// CreateInstance rents a GPU instance.
func (v *VastClient) CreateInstance(ctx context.Context, offerID int, image string, diskGB int, onstart string) (*VastInstance, error) {
	payload := map[string]any{
		"client_id": "me",
		"image":     image,
		"disk":      diskGB,
		"onstart":   onstart,
	}
	body, err := json.Marshal(payload)
	if err != nil {
		return nil, fmt.Errorf("marshal: %w", err)
	}

	data, err := v.doRequest(ctx, "PUT", fmt.Sprintf("/asks/%d/", offerID), strings.NewReader(string(body)))
	if err != nil {
		return nil, fmt.Errorf("create instance: %w", err)
	}

	var inst VastInstance
	if err := json.Unmarshal(data, &inst); err != nil {
		return nil, fmt.Errorf("parse instance: %w", err)
	}
	return &inst, nil
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
