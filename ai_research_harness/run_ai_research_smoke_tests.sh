#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ ! -x "${ROOT_DIR}/.venv/bin/python" ]]; then
  echo "缺少隔离环境：先运行 uv sync --extra test" >&2
  exit 4
fi

"${ROOT_DIR}/.venv/bin/python" -m pytest "${ROOT_DIR}/tests" -q
"${ROOT_DIR}/.venv/bin/python" -m compileall -q "${ROOT_DIR}/src"
bash -n "${ROOT_DIR}/run_ai_research.sh"
