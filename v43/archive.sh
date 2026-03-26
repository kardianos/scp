#!/bin/bash
# v43 archival script — upload to rclone remote and free disk space
# Run from /home/d/code/scp/
set -e

REMOTE="scpsfa:scpsfa"
LOCAL_ARCHIVE="/space/scp"

echo "=== Disk before ==="
df -h /home/d/code/scp/ | tail -1

# --- Step 1: Upload v43 SFA files to remote ---
echo ""
echo "=== Uploading v43 to remote ==="

# Large formation SFA — copy to remote AND local archive
echo "Uploading proton_formation.sfa (25 GB)..."
rclone copy v43/proton_formation/proton_formation.sfa $REMOTE/v43/proton_formation/ --progress
mkdir -p $LOCAL_ARCHIVE/v43/proton_formation
mv v43/proton_formation/proton_formation.sfa $LOCAL_ARCHIVE/v43/proton_formation/
echo "  Moved to $LOCAL_ARCHIVE (still on remote)"

# Analytical seeds — upload, keep local (small)
for f in proton_analytical_l2.sfa proton_analytical_l3.sfa neutron_analytical_l2.sfa uud_seed.sfa; do
    if [ -f "v43/proton_formation/$f" ]; then
        echo "Uploading $f..."
        rclone copy "v43/proton_formation/$f" $REMOTE/v43/proton_formation/
    fi
done

# Phase templates — upload, keep local (small)
for f in proton_averaged.sfa proton_peak.sfa proton_trough.sfa proton_rising.sfa proton_falling.sfa proton_template.sfa; do
    if [ -f "v43/proton_formation/$f" ]; then
        echo "Uploading $f..."
        rclone copy "v43/proton_formation/$f" $REMOTE/v43/proton_formation/
    fi
done

# --- Step 2: Upload webp files ---
echo ""
echo "=== Uploading webp files ==="
for f in v43/proton_formation/*.webp; do
    echo "Uploading $(basename $f)..."
    rclone copy "$f" $REMOTE/v43/proton_formation/
done

# --- Step 3: Upload tsv and md files ---
echo ""
echo "=== Uploading tsv/md analysis files ==="
find v43 -name "*.tsv" -o -name "*.md" | while read f; do
    dir=$(dirname "$f")
    echo "  $f -> $REMOTE/$dir/"
    rclone copy "$f" "$REMOTE/$dir/"
done

# --- Step 4: Archive and remove v34 SFA files ---
echo ""
echo "=== Archiving v34 SFA files ==="

# braid_hires.sfa — already on remote, just delete local
if [ -f "v34/torsion_coupling/data/braid_hires.sfa" ]; then
    echo "Verifying braid_hires.sfa on remote..."
    remote_size=$(rclone size $REMOTE/v34/braid_hires.sfa --json 2>/dev/null | grep -o '"bytes":[0-9]*' | grep -o '[0-9]*')
    local_size=$(stat -c%s v34/torsion_coupling/data/braid_hires.sfa)
    if [ "$remote_size" = "$local_size" ]; then
        echo "  Verified (remote=$remote_size, local=$local_size). Deleting local."
        rm v34/torsion_coupling/data/braid_hires.sfa
    else
        echo "  Size mismatch! Remote=$remote_size Local=$local_size. Keeping local."
    fi
fi

# Other v34 SFA files — upload then delete
for f in v34/torsion_coupling/data/cosserat.sfa \
         v34/torsion_coupling/data/compress_test.sfa \
         v34/torsion_coupling/theta_characterize/data/force_6field_opp.sfa \
         v34/torsion_coupling/theta_characterize/data/force_6field_same.sfa \
         v34/torsion_coupling/theta_characterize/data/winding_neg.sfa \
         v34/torsion_coupling/theta_characterize/data/winding_pos.sfa \
         v34/torsion_coupling/theta_characterize/data/force_3field.sfa; do
    if [ -f "$f" ]; then
        base=$(basename "$f")
        subdir=$(dirname "$f" | sed 's|v34/torsion_coupling/||')
        echo "Uploading $base..."
        rclone copy "$f" "$REMOTE/v34/$subdir/"
        echo "  Verifying..."
        remote_size=$(rclone size "$REMOTE/v34/$subdir/$base" --json 2>/dev/null | grep -o '"bytes":[0-9]*' | grep -o '[0-9]*')
        local_size=$(stat -c%s "$f")
        if [ "$remote_size" = "$local_size" ]; then
            echo "  Verified. Deleting local."
            rm "$f"
        else
            echo "  Size mismatch! Keeping local."
        fi
    fi
done

echo ""
echo "=== Disk after ==="
df -h /home/d/code/scp/ | tail -1

echo ""
echo "=== Done ==="
