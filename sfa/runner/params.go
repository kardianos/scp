package main

import (
	"encoding/json"
	"fmt"
	"reflect"
	"strings"
)

// --- Tool parameter and result structs ---

type SimSetupParams struct {
	Executor  string `json:"executor" desc:"Execution type: local or remote" required:"true"`
	Host      string `json:"host" desc:"SSH host for existing remote instance"`
	Port      int    `json:"port" desc:"SSH port for existing remote instance"`
	GPUFilter string `json:"gpu_filter" desc:"Vast.ai GPU filter (default: V100 16GB). Examples: 'gpu_name=Tesla_V100', 'gpu_name=RTX_4090'"`
	WorkDir   string `json:"work_dir" desc:"Working directory"`
}

type SimSetupResult struct {
	Status     string `json:"status"`
	Type       string `json:"type"`
	WorkDir    string `json:"work_dir,omitempty"`
	Host       string `json:"host,omitempty"`
	Port       int    `json:"port,omitempty"`
	GPUName    string `json:"gpu_name,omitempty"`
	InstanceID int    `json:"instance_id,omitempty"`
}

type SimStatusParams struct{}

type SimTeardownParams struct{}
type SimTeardownResult struct {
	Status string `json:"status"`
}

type SimBuildParams struct {
	Sources []string `json:"sources" desc:"Source file paths" required:"true"`
	Cmd     string   `json:"cmd" desc:"Build command (optional, auto-detected if omitted)"`
}

type SimRunParams struct {
	Config         string  `json:"config" desc:"Simulation config content (inline)" required:"true"`
	ID             string  `json:"id" desc:"Run identifier" required:"true"`
	NotifyInterval float64 `json:"notify_interval" desc:"Send progress notifications every N seconds (0=disabled)"`
	Wait           bool    `json:"wait" desc:"Block until run completes and return final status (default: false)"`
}

type SimRunResult struct {
	RunID    string `json:"run_id"`
	Status   string `json:"status"`
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
	ID string `json:"id" desc:"Run identifier" required:"true"`
}

type SimRunCancelParams struct {
	ID string `json:"id" desc:"Run identifier" required:"true"`
}

type SimRunCancelResult struct {
	Status string `json:"status"`
}

type SimUploadParams struct {
	LocalPath  string `json:"local_path" desc:"Local file path" required:"true"`
	RemotePath string `json:"remote_path" desc:"Remote destination path" required:"true"`
}

type SimUploadResult struct {
	Status string `json:"status"`
}

type SimDownloadParams struct {
	RemotePath string `json:"remote_path" desc:"Remote file path" required:"true"`
	LocalPath  string `json:"local_path" desc:"Local destination path" required:"true"`
}

type SimDownloadResult struct {
	DownloadID string `json:"download_id"`
	Status     string `json:"status"`
}

type SimDownloadStatusParams struct {
	ID string `json:"id" desc:"Download identifier" required:"true"`
}

type SimListFilesParams struct {
	Pattern string `json:"pattern" desc:"File glob pattern" required:"true"`
}

type SimListFilesResult struct {
	Files []FileInfo `json:"files"`
}

type SimExecParams struct {
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
func unmarshalParams(params map[string]any, target any) error {
	data, err := json.Marshal(params)
	if err != nil {
		return fmt.Errorf("marshal params: %w", err)
	}
	return json.Unmarshal(data, target)
}
