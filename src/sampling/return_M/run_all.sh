#!/usr/bin/env bash
set -euo pipefail

# 使い方:
#   ./run_all.sh 8   # ワーカー数
NPROCS="${1:-8}"

ROOT="$(cd "$(dirname "$0")" && pwd)"

run_dir () {
  local d="$1"
  echo "============================================================"
  echo "[DIR] $d"
  echo "[RUN] julia -p${NPROCS} run.jl"
  echo "------------------------------------------------------------"
  pushd "$ROOT/$d" >/dev/null
  julia -p"${NPROCS}" run.jl | tee "run_$(date +%Y%m%d_%H%M%S).log"
  popd >/dev/null
}

# 図A: decision comparison
run_dir "compare-decision"

# 図B: size sweep
run_dir "sweep-size"

# 図C: annealing schedule comparison
run_dir "compare-anneal"

echo "============================================================"
echo "ALL DONE."
echo "Data saved under:"
echo "  $ROOT/compare-decision/data"
echo "  $ROOT/sweep-size/data"
echo "  $ROOT/compare-anneal/data"
