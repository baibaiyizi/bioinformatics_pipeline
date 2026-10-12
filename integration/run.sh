#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "${INTEGRATION_PYTHON:-python3}" "$ROOT_DIR/scripts/integrate.py" --root "$ROOT_DIR" "$@"
