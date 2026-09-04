#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ ! -x "${ROOT_DIR}/.venv/bin/python" ]]; then
  echo "缺少隔离环境：先在 ${ROOT_DIR} 运行 uv sync --extra test" >&2
  exit 4
fi

exec "${ROOT_DIR}/.venv/bin/python" -m ai_research_harness.cli "$@"
