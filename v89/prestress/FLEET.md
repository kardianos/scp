# PRESTRESS fleet (Go)

```bash
cd v89/prestress
go build -o fleet .
./fleet -dry                 # plan
./fleet -wave 1 -max 2 -threads 4
./fleet -force               # re-run even if logs present
```

Skips jobs whose `runs/<name>.log` already has `# RESULT conservation` or
`done <name>`. Tracks only children it starts (no process greps).

Wave 1 list is in `fleet.go` (`wave1()`).
