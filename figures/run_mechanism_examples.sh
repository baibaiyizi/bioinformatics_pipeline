#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
for kind in drivers regulatory cpg reactions; do
  args=(--kind "$kind" --input "$ROOT/examples/mechanism_$kind.tsv" --output "$ROOT/previews/mechanism_$kind" --data-kind real_PFOS_source)
  if [[ "$kind" == drivers || "$kind" == reactions ]]; then
    args+=(--secondary "$ROOT/examples/mechanism_${kind}_secondary.tsv")
  fi
  "$PYTHON_BIN" "$ROOT/scripts/mechanism_recipes.py" "${args[@]}"
done
