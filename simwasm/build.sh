#!/bin/fish

set -x GOOS js
set -x GOARCH wasm
go build -o simulation.wasm simulation.go