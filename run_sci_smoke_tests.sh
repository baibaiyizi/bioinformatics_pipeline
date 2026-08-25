#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python -m unittest discover -s "${ROOT_DIR}/tests" -v

for pipeline in ml singlecell_multiome spatial_metabolomics network_pharmacology molecular_modeling; do
  python -m unittest discover -s "${ROOT_DIR}/${pipeline}/tests" -v
done
