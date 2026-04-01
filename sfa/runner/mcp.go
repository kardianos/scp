package main

import (
	"bufio"
	"context"
	"encoding/json"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"sync"
	"time"
)

// --- JSON-RPC 2.0 types ---

type jsonRPCRequest struct {
	JSONRPC string         `json:"jsonrpc"`
	ID      any            `json:"id,omitempty"`
	Method  string         `json:"method"`
	Params  map[string]any `json:"params,omitempty"`
}

type jsonRPCResponse struct {
	JSONRPC string          `json:"jsonrpc"`
	ID      any             `json:"id,omitempty"`
	Result  json.RawMessage `json:"result,omitempty"`
	Error   *jsonRPCError   `json:"error,omitempty"`
}

type jsonRPCError struct {
	Code    int    `json:"code"`
	Message string `json:"message"`
}

// --- MCP tool definition ---

type ToolDef struct {
	Name        string
	Description string
	InputType   reflect.Type
	Handler     func(ctx context.Context, params any) (any, error)
}

// --- MCP server ---

type MCPServer struct {
	tools    []ToolDef
	executor Executor
	mu       sync.RWMutex
	writeMu  sync.Mutex
	writer   *json.Encoder
	monitor  *Monitor
	state    *RunnerState
	protoLog *os.File // protocol log file (nil = disabled)
}

// prefixWriter prepends a prefix to each Write call for protocol logging.
type prefixWriter struct {
	prefix string
	w      io.Writer
}

func (pw *prefixWriter) Write(p []byte) (int, error) {
	pw.w.Write([]byte(pw.prefix))
	return pw.w.Write(p)
}

func NewMCPServer(executor Executor) *MCPServer {
	s := &MCPServer{
		executor: executor,
	}
	s.tools = s.buildToolTable()
	return s
}

func (s *MCPServer) Serve(ctx context.Context) error {
	// Set up protocol logging to ~/.scp-runner/protocol.log
	var protoLog *os.File
	if dir, err := stateDir(); err == nil {
		logPath := filepath.Join(dir, "protocol.log")
		if f, err := os.OpenFile(logPath, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0644); err == nil {
			protoLog = f
			defer f.Close()
			fmt.Fprintf(f, "=== scp-runner protocol log started %s ===\n", time.Now().Format(time.RFC3339))
		}
	}
	s.protoLog = protoLog

	// Tee stdout through the log.
	var out io.Writer = os.Stdout
	if protoLog != nil {
		out = io.MultiWriter(os.Stdout, &prefixWriter{prefix: ">>> ", w: protoLog})
	}
	s.writer = json.NewEncoder(out)

	reader := bufio.NewReader(os.Stdin)

	for {
		select {
		case <-ctx.Done():
			return ctx.Err()
		default:
		}

		line, err := reader.ReadBytes('\n')
		if err != nil {
			if err == io.EOF {
				return nil
			}
			return fmt.Errorf("read stdin: %w", err)
		}

		// Log incoming message.
		if protoLog != nil {
			fmt.Fprintf(protoLog, "<<< %s", line)
		}

		var req jsonRPCRequest
		if err := json.Unmarshal(line, &req); err != nil {
			s.sendError(nil, -32700, fmt.Sprintf("parse error: %v", err))
			continue
		}

		// Notifications (no ID) can be async; requests need synchronous response.
		if req.ID == nil {
			go s.handleRequest(ctx, &req)
		} else {
			s.handleRequest(ctx, &req)
		}
	}
}

func (s *MCPServer) handleRequest(ctx context.Context, req *jsonRPCRequest) {
	switch req.Method {
	case "initialize":
		s.handleInitialize(req)
	case "tools/list":
		s.handleToolsList(req)
	case "tools/call":
		s.handleToolsCall(ctx, req)
	case "notifications/initialized":
		// Client notification, no response needed.
	default:
		s.sendError(req.ID, -32601, fmt.Sprintf("method not found: %s", req.Method))
	}
}

func (s *MCPServer) handleInitialize(req *jsonRPCRequest) {
	result := map[string]any{
		"protocolVersion": "2025-11-25",
		"capabilities": map[string]any{
			"tools":   map[string]any{},
			"logging": map[string]any{},
			"experimental": map[string]any{
				"claude/channel": map[string]any{},
			},
		},
		"serverInfo": map[string]any{
			"name":    "scp-runner",
			"version": "0.2.0",
		},
		"instructions": "Events from the scp-runner channel arrive as <channel source=\"scp-runner\" event=\"...\">. " +
			"Events include: run_complete, run_failed, download_complete, download_failed, teardown_blocked. " +
			"When you receive a run_complete event, download the output files. " +
			"When you receive a run_failed event, diagnose the error and decide whether to retry. " +
			"When you receive a teardown_blocked event, download the listed files before retrying teardown. " +
			"When you receive a download_complete event, note the file and continue with analysis or the next step.",
	}
	s.sendResult(req.ID, result)
}

func (s *MCPServer) handleToolsList(req *jsonRPCRequest) {
	var toolList []map[string]any
	for _, t := range s.tools {
		toolList = append(toolList, map[string]any{
			"name":        t.Name,
			"description": t.Description,
			"inputSchema": structToJSONSchema(t.InputType),
		})
	}
	s.sendResult(req.ID, map[string]any{"tools": toolList})
}

func (s *MCPServer) handleToolsCall(ctx context.Context, req *jsonRPCRequest) {
	params := req.Params
	toolName, _ := params["name"].(string)
	toolArgs, _ := params["arguments"].(map[string]any)
	if toolArgs == nil {
		toolArgs = map[string]any{}
	}

	// Data-driven dispatch: find tool by name.
	for _, t := range s.tools {
		if t.Name != toolName {
			continue
		}

		// Create a new instance of the input type and unmarshal params into it.
		input := reflect.New(t.InputType).Interface()
		if err := unmarshalParams(toolArgs, input); err != nil {
			s.sendResult(req.ID, map[string]any{
				"content": []map[string]any{
					{"type": "text", "text": fmt.Sprintf("error: invalid params: %v", err)},
				},
				"isError": true,
			})
			return
		}

		result, err := t.Handler(ctx, input)
		if err != nil {
			s.sendResult(req.ID, map[string]any{
				"content": []map[string]any{
					{"type": "text", "text": fmt.Sprintf("error: %v", err)},
				},
				"isError": true,
			})
			return
		}

		text, err := json.MarshalIndent(result, "", "  ")
		if err != nil {
			s.sendError(req.ID, -32603, fmt.Sprintf("marshal result: %v", err))
			return
		}
		s.sendResult(req.ID, map[string]any{
			"content": []map[string]any{
				{"type": "text", "text": string(text)},
			},
		})
		return
	}

	s.sendError(req.ID, -32602, fmt.Sprintf("tool not found: %s", toolName))
}

func (s *MCPServer) sendResult(id any, result any) {
	data, err := json.Marshal(result)
	if err != nil {
		log.Printf("marshal result: %v", err)
		return
	}
	resp := jsonRPCResponse{
		JSONRPC: "2.0",
		ID:      id,
		Result:  json.RawMessage(data),
	}
	s.writeMu.Lock()
	defer s.writeMu.Unlock()
	s.writer.Encode(resp)
}

func (s *MCPServer) sendError(id any, code int, msg string) {
	resp := jsonRPCResponse{
		JSONRPC: "2.0",
		ID:      id,
		Error:   &jsonRPCError{Code: code, Message: msg},
	}
	s.writeMu.Lock()
	defer s.writeMu.Unlock()
	s.writer.Encode(resp)
}

// sendNotification sends a server-initiated JSON-RPC notification (no id, no response expected).
func (s *MCPServer) sendNotification(method string, params any) {
	data, err := json.Marshal(params)
	if err != nil {
		log.Printf("marshal notification: %v", err)
		return
	}
	msg := map[string]any{
		"jsonrpc": "2.0",
		"method":  method,
		"params":  json.RawMessage(data),
	}
	s.writeMu.Lock()
	defer s.writeMu.Unlock()
	s.writer.Encode(msg)
}

// sendLogMessage sends an MCP logging notification (notifications/message).
// This is the standard MCP method that clients like Claude Code surface to users.
func (s *MCPServer) sendLogMessage(level string, message string) {
	s.sendNotification("notifications/message", map[string]any{
		"level":  level,
		"logger": "scp-runner",
		"data":   message,
	})
}

// handleRunDone is called directly from the executor goroutine when a run reaches
// terminal state. It sends a channel event immediately, without waiting for the
// 5-second monitor tick.
func (s *MCPServer) handleRunDone(id string, info RunInfo) {
	switch info.Status {
	case RunComplete:
		msg := fmt.Sprintf("Run %s complete — sim_time=%.1f/%.1f wall=%.1fs", id, info.SimTime, info.TotalTime, info.WallSecs)
		s.sendChannelEvent("run_complete", map[string]string{
			"run_id":    id,
			"sim_time":  fmt.Sprintf("%.1f", info.SimTime),
			"wall_secs": fmt.Sprintf("%.1f", info.WallSecs),
		}, msg)
	case RunFailed:
		msg := fmt.Sprintf("Run %s FAILED (%.1fs): %s", id, info.WallSecs, info.Error)
		s.sendChannelEvent("run_failed", map[string]string{
			"run_id": id,
			"error":  info.Error,
		}, msg)
	case RunCancelled:
		msg := fmt.Sprintf("Run %s cancelled at sim_time=%.1f/%.1f", id, info.SimTime, info.TotalTime)
		s.sendChannelEvent("run_cancelled", map[string]string{
			"run_id": id,
		}, msg)
	}
}

// sendChannelEvent sends a channel notification that arrives in the Claude session
// as a <channel source="scp-runner" event="..." ...> tag. This interrupts the session
// and lets the agent react to events like run completion, errors, or download status.
func (s *MCPServer) sendChannelEvent(event string, meta map[string]string, content string) {
	m := map[string]string{"event": event}
	for k, v := range meta {
		m[k] = v
	}
	s.sendNotification("notifications/claude/channel", map[string]any{
		"content": content,
		"meta":    m,
	})
}

// --- Tool table ---

func (s *MCPServer) buildToolTable() []ToolDef {
	return []ToolDef{
		{
			Name:        "sim_setup",
			Description: "Create or connect to an execution environment (local or remote GPU)",
			InputType:   reflect.TypeOf(SimSetupParams{}),
			Handler:     s.handleSetup,
		},
		{
			Name:        "sim_status",
			Description: "Get execution environment status (GPU, CPU, disk, active runs)",
			InputType:   reflect.TypeOf(SimStatusParams{}),
			Handler:     s.handleStatus,
		},
		{
			Name:        "sim_teardown",
			Description: "Destroy execution environment (terminates remote instances)",
			InputType:   reflect.TypeOf(SimTeardownParams{}),
			Handler:     s.handleTeardown,
		},
		{
			Name:        "sim_build",
			Description: "Compile simulation kernel from source files",
			InputType:   reflect.TypeOf(SimBuildParams{}),
			Handler:     s.handleBuild,
		},
		{
			Name:        "sim_run",
			Description: "Start a simulation run",
			InputType:   reflect.TypeOf(SimRunParams{}),
			Handler:     s.handleRun,
		},
		{
			Name:        "sim_run_status",
			Description: "Check simulation run progress (sim_time, wall_time, status, last diagnostic)",
			InputType:   reflect.TypeOf(SimRunStatusParams{}),
			Handler:     s.handleRunStatus,
		},
		{
			Name:        "sim_run_cancel",
			Description: "Cancel a running simulation",
			InputType:   reflect.TypeOf(SimRunCancelParams{}),
			Handler:     s.handleRunCancel,
		},
		{
			Name:        "sim_upload",
			Description: "Upload a file to the execution environment",
			InputType:   reflect.TypeOf(SimUploadParams{}),
			Handler:     s.handleUpload,
		},
		{
			Name:        "sim_download",
			Description: "Download result files from the execution environment (rsync with append-verify)",
			InputType:   reflect.TypeOf(SimDownloadParams{}),
			Handler:     s.handleDownload,
		},
		{
			Name:        "sim_download_status",
			Description: "Check download progress",
			InputType:   reflect.TypeOf(SimDownloadStatusParams{}),
			Handler:     s.handleDownloadStatus,
		},
		{
			Name:        "sim_list_files",
			Description: "List files in the execution environment",
			InputType:   reflect.TypeOf(SimListFilesParams{}),
			Handler:     s.handleListFiles,
		},
		{
			Name:        "sim_exec",
			Description: "Execute an arbitrary command in the execution environment",
			InputType:   reflect.TypeOf(SimExecParams{}),
			Handler:     s.handleExec,
		},
		{
			Name:        "sim_list_templates",
			Description: "List available template SFA files for use with init=template",
			InputType:   reflect.TypeOf(SimListTemplatesParams{}),
			Handler:     s.handleListTemplates,
		},
	}
}

// --- Tool handlers ---

func (s *MCPServer) handleSetup(ctx context.Context, raw any) (any, error) {
	p := raw.(*SimSetupParams)

	workDir := p.WorkDir
	if workDir == "" {
		workDir = "/tmp/scp-runner"
	}

	switch ExecType(p.Executor) {
	case ExecLocal:
		local := NewLocalExecutor(workDir)
		if s.monitor != nil {
			local.OnRunDone = s.monitor.Wakeup
		}
		local.OnRunComplete = s.handleRunDone
		if err := local.Setup(ctx); err != nil {
			return nil, fmt.Errorf("local setup: %w", err)
		}
		s.mu.Lock()
		s.executor = local
		s.mu.Unlock()
		// Clear remote instance from persisted state for local executor.
		if s.state != nil {
			s.state.ClearInstance()
		}
		return &SimSetupResult{Status: "ok", Type: "local", WorkDir: workDir}, nil

	case ExecRemote:
		remote := NewRemoteExecutor(RemoteConfig{WorkDir: workDir})
		if s.monitor != nil {
			remote.OnRunDone = s.monitor.Wakeup
		}
		remote.OnRunComplete = s.handleRunDone
		remote.OnDownloadDone = func(id, remotePath, localPath string, err error) {
			if err != nil {
				s.sendChannelEvent("download_failed", map[string]string{
					"download_id": id,
					"remote_path": remotePath,
					"local_path":  localPath,
					"error":       err.Error(),
				}, fmt.Sprintf("Download FAILED: %s -> %s: %v", remotePath, localPath, err))
			} else {
				s.sendChannelEvent("download_complete", map[string]string{
					"download_id": id,
					"remote_path": remotePath,
					"local_path":  localPath,
				}, fmt.Sprintf("Download complete: %s -> %s", remotePath, localPath))
			}
		}
		if err := remote.Setup(ctx); err != nil {
			return nil, fmt.Errorf("remote setup: %w", err)
		}
		if p.Host != "" && p.Port > 0 {
			// Connect to existing instance.
			if err := remote.Connect(ctx, p.Host, p.Port); err != nil {
				return nil, fmt.Errorf("connect: %w", err)
			}
		} else {
			// Auto-provision: find/create a Vast.ai instance.
			if err := remote.Provision(ctx, p.GPUFilter, p.DiskGB); err != nil {
				return nil, fmt.Errorf("provision: %w", err)
			}
		}
		s.mu.Lock()
		s.executor = remote
		s.mu.Unlock()
		// Persist instance state.
		if s.state != nil {
			s.state.SetInstance(&InstanceState{
				ID:      remote.instanceID,
				Host:    remote.sshHost,
				Port:    remote.sshPort,
				GPUName: remote.gpuName,
			})
		}
		return &SimSetupResult{
			Status:     "ok",
			Type:       "remote",
			Host:       remote.sshHost,
			Port:       remote.sshPort,
			GPUName:    remote.gpuName,
			InstanceID: remote.instanceID,
		}, nil

	default:
		return nil, fmt.Errorf("unknown type: %s (use 'local' or 'remote')", p.Executor)
	}
}

func (s *MCPServer) getExecutor() (Executor, error) {
	s.mu.RLock()
	defer s.mu.RUnlock()
	if s.executor == nil {
		return nil, fmt.Errorf("no executor configured — call sim_setup first")
	}
	return s.executor, nil
}

func (s *MCPServer) handleStatus(ctx context.Context, _ any) (any, error) {
	ex, err := s.getExecutor()
	if err != nil {
		return nil, err
	}
	return ex.Status(ctx)
}

func (s *MCPServer) handleTeardown(ctx context.Context, raw any) (any, error) {
	p := raw.(*SimTeardownParams)
	ex, err := s.getExecutor()
	if err != nil {
		return nil, err
	}

	// Safety check: refuse teardown if there are undownloaded output files.
	if !p.Force {
		if pending := s.pendingDownloads(); len(pending) > 0 {
			msg := fmt.Sprintf("TEARDOWN BLOCKED: %d output file(s) not yet downloaded:\n", len(pending))
			for _, f := range pending {
				msg += fmt.Sprintf("  - %s\n", f)
			}
			msg += "Use force=true to override, or download files first."
			s.sendChannelEvent("teardown_blocked", map[string]string{
				"count": fmt.Sprintf("%d", len(pending)),
			}, msg)
			return nil, fmt.Errorf("teardown blocked: %d output file(s) not downloaded — use force=true to override", len(pending))
		}
	}

	if err := ex.Teardown(ctx); err != nil {
		return nil, fmt.Errorf("teardown: %w", err)
	}
	s.mu.Lock()
	s.executor = nil
	s.mu.Unlock()
	// Clear persisted state.
	if s.state != nil {
		s.state.ClearInstance()
	}
	return &SimTeardownResult{Status: "ok"}, nil
}

// pendingDownloads returns a list of remote output files that have not been
// downloaded locally. It checks the persisted run state for output files and
// verifies whether corresponding local files exist with non-zero size.
func (s *MCPServer) pendingDownloads() []string {
	if s.state == nil {
		return nil
	}
	s.state.mu.Lock()
	runs := make(map[string]*PersistRun, len(s.state.Runs))
	for k, v := range s.state.Runs {
		runs[k] = v
	}
	s.state.mu.Unlock()

	var pending []string
	for _, run := range runs {
		if run.AutoDownload == "" {
			continue
		}
		for _, remote := range run.OutputFiles {
			// Compute expected local path.
			base := filepath.Base(remote)
			local := filepath.Join(run.AutoDownload, base)
			info, err := os.Stat(local)
			if err != nil || info.Size() == 0 {
				pending = append(pending, remote+" -> "+local)
			}
		}
	}

	// Also check active downloads that haven't completed.
	ex, err := s.getExecutor()
	if err != nil {
		return pending
	}
	downloads := make(map[string]*DownloadInfo)
	collectDownloads(ex, downloads)
	for _, dl := range downloads {
		if dl.Status == "running" || dl.Status == "started" {
			pending = append(pending, dl.RemotePath+" -> "+dl.LocalPath+" (in progress)")
		}
	}

	return pending
}

func (s *MCPServer) handleBuild(ctx context.Context, raw any) (any, error) {
	p := raw.(*SimBuildParams)

	ex, err := s.getExecutor()
	if err != nil {
		return nil, err
	}
	if len(p.Sources) == 0 {
		return nil, fmt.Errorf("sources required")
	}
	result, err := ex.Build(ctx, p.Sources, p.Cmd)
	if err != nil {
		return nil, err
	}
	// Persist binary path on successful build.
	if s.state != nil && result != nil && result.Status != "failed" {
		s.state.SetBinary(result.Binary)
	}
	return result, nil
}

// validInitModes lists the init modes supported by the simulation kernel.
var validInitModes = map[string]string{
	"oscillon": "Analytical oscillon initialization",
	"braid":    "Braid (topological winding) initialization",
	"sfa":      "Load full state from an SFA file (init_sfa=path, init_frame=N)",
	"exec":     "Run an external program to generate initial conditions (init_exec=cmd)",
	"template": "Stamp a small template SFA into an analytical background (init_sfa=path)",
}

// replaceConfigValue replaces the value of a key in key=value config content.
func replaceConfigValue(config, key, newVal string) string {
	var lines []string
	for _, line := range strings.Split(config, "\n") {
		trimmed := strings.TrimSpace(line)
		parts := strings.SplitN(trimmed, "=", 2)
		if len(parts) == 2 && strings.TrimSpace(parts[0]) == key {
			lines = append(lines, key+"="+newVal)
		} else {
			lines = append(lines, line)
		}
	}
	return strings.Join(lines, "\n")
}

// parseConfigValue extracts a key's value from key=value config content.
func parseConfigValue(config, key string) string {
	for _, line := range strings.Split(config, "\n") {
		line = strings.TrimSpace(line)
		if line == "" || line[0] == '#' {
			continue
		}
		parts := strings.SplitN(line, "=", 2)
		if len(parts) == 2 && strings.TrimSpace(parts[0]) == key {
			return strings.TrimSpace(parts[1])
		}
	}
	return ""
}

func (s *MCPServer) handleRun(ctx context.Context, raw any) (any, error) {
	p := raw.(*SimRunParams)

	ex, err := s.getExecutor()
	if err != nil {
		return nil, err
	}
	if p.Config == "" || p.ID == "" {
		return nil, fmt.Errorf("config and id required")
	}

	// Validate init mode if config is key=value content.
	var initMode, initSFA string
	if isConfigContent(p.Config) {
		initMode = parseConfigValue(p.Config, "init")
		initSFA = parseConfigValue(p.Config, "init_sfa")

		if initMode == "" {
			return nil, fmt.Errorf("config missing 'init' key — valid modes: oscillon, braid, sfa, exec, template")
		}
		if _, ok := validInitModes[initMode]; !ok {
			modes := make([]string, 0, len(validInitModes))
			for m := range validInitModes {
				modes = append(modes, m)
			}
			return nil, fmt.Errorf("unknown init mode %q — valid modes: %s", initMode, strings.Join(modes, ", "))
		}

		// For template and sfa modes, validate that init_sfa is provided and the file exists.
		if initMode == "template" || initMode == "sfa" {
			if initSFA == "" {
				return nil, fmt.Errorf("init=%s requires init_sfa=<path> to be set in config — use sim_list_templates to see available templates", initMode)
			}
			// Resolve bare names (e.g. "proton") to sfa/templates/proton.sfa.
			initSFA = resolveTemplatePath(initSFA)
			// Rewrite the config with the resolved absolute path.
			p.Config = replaceConfigValue(p.Config, "init_sfa", initSFA)

			// For local executor, check if the file exists.
			if ex.Type() == ExecLocal {
				if _, err := os.Stat(initSFA); err != nil {
					return nil, fmt.Errorf("init_sfa file not found: %s — use sim_list_templates to see available templates", initSFA)
				}
			}
		}
	}

	notifyInterval := time.Duration(p.NotifyInterval * float64(time.Second))

	// Use auto-download variant for remote executor when auto_download is set.
	if p.AutoDownload != "" {
		if remote, ok := ex.(*RemoteExecutor); ok {
			if err := remote.RunWithAutoDownload(ctx, p.Config, p.ID, notifyInterval, p.AutoDownload); err != nil {
				return nil, err
			}
		} else {
			// Local executor: ignore auto_download, just run normally.
			if err := ex.Run(ctx, p.Config, p.ID, notifyInterval); err != nil {
				return nil, err
			}
		}
	} else {
		if err := ex.Run(ctx, p.Config, p.ID, notifyInterval); err != nil {
			return nil, err
		}
	}

	// Persist run state.
	if s.state != nil {
		var outputFiles []string
		if isConfigContent(p.Config) {
			if out := parseConfigValue(p.Config, "output"); out != "" {
				outputFiles = append(outputFiles, out)
			}
			if diag := parseConfigValue(p.Config, "diag_file"); diag != "" {
				outputFiles = append(outputFiles, diag)
			}
		}
		s.state.SetRun(p.ID, &PersistRun{
			Status:       "running",
			AutoDownload: p.AutoDownload,
			OutputFiles:  outputFiles,
			RemoteDir:    "/root",
		})
	}

	if !p.Wait {
		return &SimRunResult{
			RunID:    p.ID,
			Status:   "started",
			Executor: string(ex.Type()),
			Init:     initMode,
			InitSFA:  initSFA,
		}, nil
	}

	// Block until run completes, polling every 2 seconds.
	ticker := time.NewTicker(2 * time.Second)
	defer ticker.Stop()
	for {
		select {
		case <-ctx.Done():
			return nil, ctx.Err()
		case <-ticker.C:
			info := ex.RunStatus(p.ID)
			if info == nil {
				return nil, fmt.Errorf("run %s disappeared", p.ID)
			}
			switch info.Status {
			case RunComplete, RunFailed, RunCancelled:
				return info, nil
			}
		}
	}
}

func (s *MCPServer) handleRunStatus(_ context.Context, raw any) (any, error) {
	p := raw.(*SimRunStatusParams)

	ex, err := s.getExecutor()
	if err != nil {
		return nil, err
	}
	if p.ID == "" {
		return nil, fmt.Errorf("id required")
	}
	info := ex.RunStatus(p.ID)
	if info == nil {
		return nil, fmt.Errorf("run %s not found", p.ID)
	}
	return info, nil
}

func (s *MCPServer) handleRunCancel(_ context.Context, raw any) (any, error) {
	p := raw.(*SimRunCancelParams)

	ex, err := s.getExecutor()
	if err != nil {
		return nil, err
	}
	if p.ID == "" {
		return nil, fmt.Errorf("id required")
	}
	if err := ex.RunCancel(p.ID); err != nil {
		return nil, err
	}
	return &SimRunCancelResult{Status: "cancelled"}, nil
}

func (s *MCPServer) handleUpload(ctx context.Context, raw any) (any, error) {
	p := raw.(*SimUploadParams)

	ex, err := s.getExecutor()
	if err != nil {
		return nil, err
	}
	if p.LocalPath == "" || p.RemotePath == "" {
		return nil, fmt.Errorf("local_path and remote_path required")
	}
	if err := ex.Upload(ctx, p.LocalPath, p.RemotePath); err != nil {
		return nil, err
	}
	return &SimUploadResult{Status: "ok"}, nil
}

func (s *MCPServer) handleDownload(ctx context.Context, raw any) (any, error) {
	p := raw.(*SimDownloadParams)

	ex, err := s.getExecutor()
	if err != nil {
		return nil, err
	}
	if p.RemotePath == "" || p.LocalPath == "" {
		return nil, fmt.Errorf("remote_path and local_path required")
	}
	id, err := ex.Download(ctx, p.RemotePath, p.LocalPath)
	if err != nil {
		return nil, err
	}
	return &SimDownloadResult{DownloadID: id, Status: "started"}, nil
}

func (s *MCPServer) handleDownloadStatus(_ context.Context, raw any) (any, error) {
	p := raw.(*SimDownloadStatusParams)

	ex, err := s.getExecutor()
	if err != nil {
		return nil, err
	}
	if p.ID == "" {
		return nil, fmt.Errorf("id required")
	}
	info := ex.DownloadStatus(p.ID)
	if info == nil {
		return nil, fmt.Errorf("download %s not found", p.ID)
	}
	return info, nil
}

func (s *MCPServer) handleListFiles(ctx context.Context, raw any) (any, error) {
	p := raw.(*SimListFilesParams)

	ex, err := s.getExecutor()
	if err != nil {
		return nil, err
	}
	if p.Pattern == "" {
		return nil, fmt.Errorf("pattern required")
	}
	files, err := ex.ListFiles(ctx, p.Pattern)
	if err != nil {
		return nil, err
	}
	return &SimListFilesResult{Files: files}, nil
}

// findTemplatesDir locates the sfa/templates/ directory by walking up from cwd.
func findTemplatesDir() string {
	cwd, _ := os.Getwd()
	for d := cwd; d != "/"; d = filepath.Dir(d) {
		td := filepath.Join(d, "sfa", "templates")
		if info, err := os.Stat(td); err == nil && info.IsDir() {
			return td
		}
	}
	return ""
}

// resolveTemplatePath resolves an init_sfa value to an absolute path.
// If the value is already absolute and exists, it is returned as-is.
// Otherwise, it checks the sfa/templates/ directory for a match
// (with or without .sfa extension).
func resolveTemplatePath(initSFA string) string {
	// Already absolute and exists — use directly.
	if filepath.IsAbs(initSFA) {
		if _, err := os.Stat(initSFA); err == nil {
			return initSFA
		}
	}

	td := findTemplatesDir()
	if td == "" {
		return initSFA
	}

	// Try exact name in templates dir.
	candidate := filepath.Join(td, initSFA)
	if _, err := os.Stat(candidate); err == nil {
		return candidate
	}

	// Try with .sfa extension.
	if !strings.HasSuffix(initSFA, ".sfa") {
		candidate = filepath.Join(td, initSFA+".sfa")
		if _, err := os.Stat(candidate); err == nil {
			return candidate
		}
	}

	return initSFA
}

func (s *MCPServer) handleListTemplates(_ context.Context, _ any) (any, error) {
	td := findTemplatesDir()
	if td == "" {
		return nil, fmt.Errorf("sfa/templates/ directory not found — expected at <project>/sfa/templates/")
	}

	entries, err := os.ReadDir(td)
	if err != nil {
		return nil, fmt.Errorf("read templates dir: %w", err)
	}

	var templates []TemplateInfo
	for _, e := range entries {
		if e.IsDir() || !strings.HasSuffix(e.Name(), ".sfa") {
			continue
		}
		info, err := e.Info()
		if err != nil {
			continue
		}
		templates = append(templates, TemplateInfo{
			Path:    filepath.Join(td, e.Name()),
			Name:    strings.TrimSuffix(e.Name(), ".sfa"),
			Size:    info.Size(),
			ModTime: info.ModTime().Format(time.RFC3339),
		})
	}

	return &SimListTemplatesResult{Templates: templates}, nil
}

func (s *MCPServer) handleExec(ctx context.Context, raw any) (any, error) {
	p := raw.(*SimExecParams)

	ex, err := s.getExecutor()
	if err != nil {
		return nil, err
	}
	if p.Cmd == "" {
		return nil, fmt.Errorf("cmd required")
	}
	timeoutMs := p.TimeoutMs
	if timeoutMs == 0 {
		timeoutMs = 30000
	}
	timeout := time.Duration(timeoutMs) * time.Millisecond
	out, err := ex.Exec(ctx, p.Cmd, timeout)
	if err != nil {
		return &SimExecResult{Output: out, Error: err.Error()}, nil
	}
	return &SimExecResult{Output: out, Status: "ok"}, nil
}
