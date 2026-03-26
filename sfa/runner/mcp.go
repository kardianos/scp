package main

import (
	"bufio"
	"context"
	"encoding/json"
	"fmt"
	"io"
	"log"
	"os"
	"reflect"
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
}

func NewMCPServer(executor Executor) *MCPServer {
	s := &MCPServer{
		executor: executor,
	}
	s.tools = s.buildToolTable()
	return s
}

func (s *MCPServer) Serve(ctx context.Context) error {
	s.writer = json.NewEncoder(os.Stdout)
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
		"protocolVersion": "2024-11-05",
		"capabilities": map[string]any{
			"tools": map[string]any{},
		},
		"serverInfo": map[string]any{
			"name":    "scp-runner",
			"version": "0.1.0",
		},
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
		if err := local.Setup(ctx); err != nil {
			return nil, fmt.Errorf("local setup: %w", err)
		}
		s.mu.Lock()
		s.executor = local
		s.mu.Unlock()
		return &SimSetupResult{Status: "ok", Type: "local", WorkDir: workDir}, nil

	case ExecRemote:
		remote := NewRemoteExecutor(RemoteConfig{WorkDir: workDir})
		if err := remote.Setup(ctx); err != nil {
			return nil, fmt.Errorf("remote setup: %w", err)
		}
		if p.Host != "" && p.Port > 0 {
			if err := remote.Connect(ctx, p.Host, p.Port); err != nil {
				return nil, fmt.Errorf("connect: %w", err)
			}
		}
		s.mu.Lock()
		s.executor = remote
		s.mu.Unlock()
		return &SimSetupResult{Status: "ok", Type: "remote", Host: p.Host}, nil

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

func (s *MCPServer) handleTeardown(ctx context.Context, _ any) (any, error) {
	ex, err := s.getExecutor()
	if err != nil {
		return nil, err
	}
	if err := ex.Teardown(ctx); err != nil {
		return nil, fmt.Errorf("teardown: %w", err)
	}
	s.mu.Lock()
	s.executor = nil
	s.mu.Unlock()
	return &SimTeardownResult{Status: "ok"}, nil
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
	return ex.Build(ctx, p.Sources, p.Cmd)
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
	if err := ex.Run(ctx, p.Config, p.ID); err != nil {
		return nil, err
	}
	return &SimRunResult{RunID: p.ID, Status: "started", Executor: string(ex.Type())}, nil
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
