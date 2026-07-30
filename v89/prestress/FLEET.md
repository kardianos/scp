# PRESTRESS fleet (Go)

```bash
cd v89/prestress
go build -o fleet .
./fleet -dry -wave 2-4       # plan remaining program
./fleet -wave 1|2|3|4|all|2-4 -max 2 -threads 4
./fleet -force               # re-run even if logs present
```

Skips jobs whose `runs/<name>.log` already has `# RESULT conservation` or
`done <name>`. Tracks only children it starts (no process greps).

Waves defined in `fleet.go`: W1 frozen, W2 plast+topo, W3 discriminators,
W4 harden / m-family (single-foam P19-lite).
