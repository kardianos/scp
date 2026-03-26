package main

import (
	"context"
	"log"
	"os"
	"os/signal"
	"syscall"
)

func main() {
	log.SetFlags(log.Ltime | log.Lshortfile)
	log.SetOutput(os.Stderr)

	ctx, cancel := signal.NotifyContext(context.Background(), os.Interrupt, syscall.SIGTERM)
	defer cancel()

	// Create MCP server with no executor initially (sim_setup will configure one).
	server := NewMCPServer(nil)

	// Start background monitor and link it to the server.
	monitor := NewMonitor(server, "/tmp/scp-runner")
	server.monitor = monitor
	go monitor.Run(ctx)

	log.Println("scp-runner: MCP server starting on stdio")
	if err := server.Serve(ctx); err != nil {
		if ctx.Err() != nil {
			log.Println("scp-runner: shutdown")
			return
		}
		log.Fatalf("scp-runner: %v", err)
	}
}
