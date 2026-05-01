#!/usr/bin/env fish
# run_bench.fish — Benchmark harness for foam_sim.
#
# Runs a fixed config (bench.cfg) on a fixed mesh (default foam500k.bin)
# and records timing + correctness signals to results.tsv. Each run is
# tagged with a version label so we can compare optimizations.
#
# Usage:
#   run_bench.fish <version_label> [<mesh.bin>] [<notes...>]
#   ./run_bench.fish baseline
#   ./run_bench.fish opt1_face_aos
#   ./run_bench.fish opt2_morton
#   ./run_bench.fish opt3_combined
#
# Outputs:
#   results.tsv           — append-only TSV of all runs
#   logs/log_<ts>.txt     — full stdout from each run
#   diag/diag_<ts>.tsv    — per-run diagnostics

set bench_dir (realpath (dirname (status --current-filename)))
set foam_dir (dirname $bench_dir)
set archive /space/scp/v55/bench

set tag $argv[1]
if test -z "$tag"
    echo "Usage: run_bench.fish <tag> [<mesh.bin>] [<notes...>]"
    exit 1
end

# 2nd arg is mesh ONLY if it's a .bin file that exists; else it's notes.
set mesh_arg $argv[2]
if test -z "$mesh_arg"
    set mesh $foam_dir/foam500k.bin
    set notes (string join " " $argv[3..-1])
else if string match -qr '\.bin$' $mesh_arg
    if string match -q '/*' $mesh_arg
        set mesh $mesh_arg
    else
        set mesh $foam_dir/$mesh_arg
    end
    set notes (string join " " $argv[3..-1])
else
    set mesh $foam_dir/foam500k.bin
    set notes (string join " " $argv[2..-1])
end
set mesh_name (basename $mesh .bin)

# Verify mesh and binary exist
if not test -e $mesh
    echo "ERROR: mesh not found: $mesh"
    exit 1
end
if not test -e $foam_dir/foam_sim
    echo "ERROR: foam_sim binary missing — run `make` first"
    exit 1
end

mkdir -p $bench_dir/logs $bench_dir/diag $archive
set ts (date -u +%Y%m%dT%H%M%S)
set log $bench_dir/logs/log_$ts.txt
set diag $bench_dir/diag/diag_$ts.tsv

# Get binary mtime for hash-of-binary comparison
set bin_mtime (stat -c '%Y' $foam_dir/foam_sim)
set bin_size (stat -c '%s' $foam_dir/foam_sim)

# Get N_cells from mesh header (offset 16: uint32 after magic+ver+L=4+4+8)
set n_cells (od -An -tu4 -j16 -N4 $mesh | string trim)

echo "==================================================================="
echo "Bench: $tag"
echo "  mesh:    $mesh_name ($n_cells cells)"
echo "  binary:  foam_sim ($bin_size bytes, mtime $bin_mtime)"
echo "  notes:   $notes"
echo "==================================================================="

# Run, capturing stdout
cd $foam_dir
set start_ns (date +%s%N)
env OMP_NUM_THREADS=8 ./foam_sim $mesh $bench_dir/bench.cfg > $log 2>&1
set rc $status
set end_ns (date +%s%N)
set wall_sec (math "($end_ns - $start_ns) / 1000000000")

if test $rc -ne 0
    echo "FAILED (exit $rc); see $log"
    exit $rc
end

# Move diag file out of foam_dir into bench/diag/ for archival
if test -e $foam_dir/bench_diag.tsv
    mv $foam_dir/bench_diag.tsv $diag
end

# Parse metrics from log
set ms_step (grep -oE '[0-9]+\.?[0-9]*ms/step' $log | tail -1 | string replace 'ms/step' '')
set final_drift (grep 'drift' $log | tail -1 | grep -oE '[+-][0-9]+\.[0-9]+%' | string replace '%' '')
set final_theta_rms (grep 'theta' $log | tail -1 | grep -oE 'θ_rms=[0-9]+\.[0-9]+' | string replace 'θ_rms=' '')
set complete_wall (grep 'Complete:' $log | grep -oE '[0-9]+s')
set n_steps (grep 'n_steps=' $log | head -1 | grep -oE 'n_steps=[0-9]+' | string replace 'n_steps=' '')

# Derive ms/step if missing
if test -z "$ms_step" -a -n "$n_steps"
    set ms_step (math "$wall_sec * 1000 / $n_steps")
end

set results $bench_dir/results.tsv
if not test -e $results
    echo -e "timestamp\tversion\tmesh\tn_cells\tT\tthreads\twall_sec\tms_per_step\tn_steps\tfinal_drift_pct\tfinal_theta_rms\tnotes" > $results
end
echo -e "$ts\t$tag\t$mesh_name\t$n_cells\t10\t8\t$wall_sec\t$ms_step\t$n_steps\t$final_drift\t$final_theta_rms\t$notes" >> $results

# Mirror to archive
cp $results $archive/
cp $log $archive/

echo ""
echo "Result:"
echo "  wall:     $wall_sec sec"
echo "  ms/step:  $ms_step"
echo "  n_steps:  $n_steps"
echo "  drift:    $final_drift%"
echo "  θ_rms:    $final_theta_rms"
echo "  log:      $log"
echo ""
echo "All results: $results"
