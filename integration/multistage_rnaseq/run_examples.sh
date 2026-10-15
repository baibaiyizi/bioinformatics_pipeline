#!/usr/bin/env bash
set -euo pipefail
root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bash "$root_dir/examples/perturbation/run_example.sh"
bash "$root_dir/examples/splicing/run_example.sh"
