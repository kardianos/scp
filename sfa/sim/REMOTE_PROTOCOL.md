# Remote GPU Execution Protocol

## CRITICAL: Instance Lifecycle Rules

### Rule 1: NEVER destroy in the same command as download
Downloads and destroy MUST be separate steps with verification between them.

### Rule 2: Verify downloads BEFORE destroying
After every download, verify:
```bash
# On remote: get expected size
ssh remote "ls -l file.sfa | awk '{print \$5}'"
# Locally: check actual size
ls -l local_copy.sfa | awk '{print $5}'
# Compare: must match exactly
```

### Rule 3: Use a dedicated monitor agent for long runs
For any run > 30 minutes:
1. Launch the simulation via nohup
2. Launch a MONITOR AGENT that:
   - Checks simulation progress every 5 minutes
   - When sim completes: runs analysis + conversion on remote
   - Verifies analysis/conversion completed
   - Downloads all files with size verification
   - Only AFTER all downloads verified → destroys instance
3. The monitor agent is the ONLY entity that destroys the instance

### Rule 4: Analysis and conversion happen BEFORE download
On the remote, BEFORE downloading anything:
1. Run analyze_sfa → verify JSON exists and is non-empty
2. Run sfa_convert → verify f16 file exists and size is reasonable
3. Get file sizes: `ls -l *.sfa *.json *.tsv > manifest.txt`
4. Download manifest.txt FIRST
5. Download each file and verify size against manifest
6. Only after ALL files verified → destroy

### Download Verification Script
```bash
verify_download() {
    local remote_path=$1 local_path=$2 ssh_cmd=$3
    local remote_size=$($ssh_cmd "stat -c%s $remote_path" 2>/dev/null)
    local local_size=$(stat -c%s "$local_path" 2>/dev/null)
    if [ "$remote_size" = "$local_size" ] && [ -n "$remote_size" ]; then
        echo "VERIFIED: $local_path ($local_size bytes)"
        return 0
    else
        echo "MISMATCH: remote=$remote_size local=$local_size"
        return 1
    fi
}
```

### Template: Safe Remote Execution
```bash
# 1. Run simulation
ssh remote "nohup ./scp_sim_cuda config.cfg > log.txt 2>&1 &"

# 2. Monitor until complete (separate agent or polling loop)
while ! ssh remote "grep -q DONE log.txt"; do sleep 60; done

# 3. Analyze and convert ON REMOTE
ssh remote "
    ./analyze_sfa output.sfa --json analysis.json
    ./sfa_convert output.sfa output_f16.sfa --f16
    ls -l *.sfa *.json *.tsv > manifest.txt
"

# 4. Download with verification
scp remote:manifest.txt ./
for file in $(awk '{print $NF}' manifest.txt); do
    scp remote:$file ./
    verify_download "$file" "./$file" "ssh remote"
done

# 5. ONLY after all verified
if all_verified; then
    vastai destroy instance $ID
else
    echo "DOWNLOAD INCOMPLETE — instance kept alive"
fi
```
