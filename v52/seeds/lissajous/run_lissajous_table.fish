#!/usr/bin/env fish
# run_lissajous_table.fish — 5×5×5 Lissajous phase sweep
#
# Superimposes 3 rotated blob copies with phase shifts Δ₀,Δ₁,Δ₂ ∈ {0, π/3, 2π/3, π, 4π/3}
# Runs each as N=64 periodic BC simulation for T=100
# Outputs to /space/scp/v52/lissajous/

set BASE_SEED (pwd)/v52/chirality_test/base_seed.sfa
set GEN (pwd)/v52/chirality_test/gen_lissajous_seed
set SIM (pwd)/bin/scp_sim
set OUTDIR /space/scp/v52/lissajous
set SEEDDIR $OUTDIR/seeds
set RESULTS $OUTDIR/results.tsv

mkdir -p $OUTDIR $SEEDDIR

# Phase values: 0, π/3, 2π/3, π, 4π/3
set PHASES 0.0000 1.0472 2.0944 3.1416 4.1888
set LABELS 0 p3 2p3 p 4p3

# Header for results
echo -e "d0\td1\td2\tE_init\tE_final\tdrift\tphi_max_i\tphi_max_f\tP_int_i\tP_int_f\ttheta_rms_f\tstatus" > $RESULTS

set total (math (count $PHASES) \* (count $PHASES) \* (count $PHASES))
set run 0

for i in (seq (count $PHASES))
    for j in (seq (count $PHASES))
        for k in (seq (count $PHASES))
            set d0 $PHASES[$i]; set d1 $PHASES[$j]; set d2 $PHASES[$k]
            set l0 $LABELS[$i]; set l1 $LABELS[$j]; set l2 $LABELS[$k]
            set name "d{$l0}_{$l1}_{$l2}"
            set run (math $run + 1)

            set seed_file $SEEDDIR/$name.sfa
            set out_file $OUTDIR/$name.sfa
            set diag_file $OUTDIR/$name.tsv
            set cfg_file $OUTDIR/$name.cfg

            # Skip if output already exists
            if test -f $diag_file
                echo "[$run/$total] $name — exists, skipping"
                continue
            end

            echo -n "[$run/$total] $name (Δ=$d0,$d1,$d2) ... "

            # Generate seed
            $GEN $BASE_SEED $seed_file $d0 $d1 $d2 > /dev/null 2>&1

            # Write config
            echo "\
N=64
L=8.3
T=100
dt_factor=0.025
m=1.5
m_theta=0
eta=0.5
mu=-41.345
kappa=50
bc_type=2
damp_width=0
damp_rate=0
precision=f16
init=sfa
init_sfa=$seed_file
output=$out_file
diag_file=$diag_file
snap_dt=25
diag_dt=2" > $cfg_file

            # Run simulation
            set result (OMP_NUM_THREADS=8 $SIM $cfg_file 2>&1)

            # Extract metrics
            set init_line (echo $result | grep -oP 'INIT: E_total=\K[^ ]+')
            set init_phi (echo $result | grep -oP 'INIT:.*phi_max=\K[^ ]+')
            set init_pint (echo $result | grep -oP 'INIT:.*P_int=\K[^ ]+')

            set final_e (echo $result | grep "=== COMPLETE" -A3 | grep -oP 'E_total=\K[^ ]+')
            set drift (echo $result | grep "=== COMPLETE" -A3 | grep -oP 'drift \K[^)]+')
            set final_phi (echo $result | grep "=== COMPLETE" -A5 | grep -oP 'phi_max=\K[^ ]+')
            set final_pint (echo $result | grep "=== COMPLETE" -A5 | grep -oP 'P_int=\K[^ ]+')
            set final_trms (echo $result | grep "=== COMPLETE" -A5 | grep -oP 'theta_rms=\K[^ ]+')

            # Determine status
            set status "dissolved"
            if test -n "$final_phi"
                set phi_val (echo $final_phi | tr -d '%')
                # Check if phi_max > 0.3 (soliton survives)
                if math "$phi_val > 0.3" > /dev/null 2>&1
                    set status "alive"
                end
            end

            echo -e "$d0\t$d1\t$d2\t$init_line\t$final_e\t$drift\t$init_phi\t$final_phi\t$init_pint\t$final_pint\t$final_trms\t$status" >> $RESULTS

            echo "$status (φ=$final_phi P=$final_pint)"

            # Clean up seed and config to save space
            rm -f $seed_file $cfg_file
        end
    end
end

echo ""
echo "=== COMPLETE ==="
echo "Results: $RESULTS"
echo "SFA files: $OUTDIR/"

# Summary
set alive (grep -c "alive" $RESULTS)
set dead (grep -c "dissolved" $RESULTS)
echo "Survivors: $alive / $total"
echo "Dissolved: $dead / $total"
