#!/usr/bin/env bash
set -euo pipefail
TEMPLATE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="$TEMPLATE_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
exec "${SCRNA_PYTHON:-python}" -m scrna_template "$@"
