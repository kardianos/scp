#!/bin/bash
set -euo pipefail
cd /home/d/code/scp/v93/theory/qgraph
export PATH="/home/d/.local/bin:$PATH"
# default model (fable/max not available on this install)
claude --permission-mode bypassPermissions \
  --add-dir /home/d/code/scp/v93/theory \
  --add-dir /home/d/code/scp/v93/laws \
  --add-dir /space/scp/v93 \
  -p "$(cat /home/d/code/scp/v93/theory/qgraph/BRIEF_claude.md)" \
  > /home/d/code/scp/v93/theory/qgraph/REVIEW_claude.md \
  2> /home/d/code/scp/v93/theory/qgraph/REVIEW_claude.err
ec=$?
echo EXIT:$ec >> /home/d/code/scp/v93/theory/qgraph/REVIEW_claude.err
touch /home/d/code/scp/v93/theory/qgraph/REVIEW_claude.done
exit $ec
