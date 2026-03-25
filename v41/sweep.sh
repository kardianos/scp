#!/bin/bash
# V41: Sweep 72 seed configurations (36 UDD + 36 UUD), rank by predicted stability
set -e
cd /home/d/code/scp/v41

mkdir -p seeds results

echo "chirality	A	R_tube	theta_init	v_profile	P_peak	P_opt	theta_ratio	v_ratio	rho_ratio	E_pot	P_conc	force_bal	S_pred	path" > sweep_results.tsv

N=192; L=25

idx=0
for chiral in UDD UUD; do
for A in 0.3 0.5 0.7; do
for R in 2.5 3.5 4.5; do
for vp in zero breathing contracting mixed; do
    idx=$((idx + 1))
    id=$(printf "%s_%03d" "$chiral" "$idx")
    seed="seeds/${id}.sfa"

    ./construct_seed -N $N -L $L -A $A -R $R \
        -theta_init 0.5 -v_profile $vp \
        -chirality $chiral -o "$seed" -precision f32 \
        >> sweep_results.tsv 2>/dev/null

    printf "\r%d/72 %s" $idx "$id"
done
done
done
done

echo ""
echo ""
echo "=== TOP 3 UDD (by S_pred) ==="
grep "^UDD" sweep_results.tsv | sort -t$'\t' -k14 -rn | head -3
echo ""
echo "=== TOP 3 UUD (by S_pred) ==="
grep "^UUD" sweep_results.tsv | sort -t$'\t' -k14 -rn | head -3
echo ""

# Extract the top 3 seed paths for each chirality
echo "=== Winners ==="
mkdir -p winners
for chiral in UDD UUD; do
    rank=0
    grep "^${chiral}" sweep_results.tsv | sort -t$'\t' -k14 -rn | head -3 | while read line; do
        rank=$((rank + 1))
        path=$(echo "$line" | cut -f15)
        A=$(echo "$line" | cut -f2)
        R=$(echo "$line" | cut -f3)
        vp=$(echo "$line" | cut -f5)
        score=$(echo "$line" | cut -f14)
        dest="winners/${chiral}_rank${rank}.sfa"
        cp "$path" "$dest"
        echo "${chiral} #${rank}: S_pred=${score} A=${A} R=${R} vp=${vp} → ${dest}"
    done
done

echo ""
echo "Done. Results in sweep_results.tsv, winners in winners/"
