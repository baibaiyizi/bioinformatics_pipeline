#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$ROOT/.cache/matplotlib}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
"$PYTHON_BIN" "$ROOT/scripts/render.py" --kind ora --input "$ROOT/examples/marker_ora.tsv" --output "$ROOT/previews/python_ora" --value-label 'Marker ORA; annotation only'
"$PYTHON_BIN" "$ROOT/scripts/render.py" --kind marker --input "$ROOT/examples/marker_expression.tsv" --output "$ROOT/previews/python_marker" --value-label 'Within-gene cluster z-score; ≤−2.5 / ≥2.5 capped'
"$PYTHON_BIN" "$ROOT/scripts/render.py" --kind evidence_matrix --input "$ROOT/examples/program_matrix.tsv" --output "$ROOT/previews/python_program_matrix" --value-label 'Module-wise visual z-score; descriptive n=1/library'
"$PYTHON_BIN" "$ROOT/scripts/render.py" --kind volcano --input "$ROOT/examples/metabolomics_volcano.tsv" --output "$ROOT/previews/python_volcano" --value-label 'log2FC (PFOS − WT); frozen limma effect'
"$PYTHON_BIN" "$ROOT/scripts/render.py" --kind composition --input "$ROOT/examples/composition.tsv" --output "$ROOT/previews/python_composition" --value-label 'Captured cell proportion; n=1/library'
Rscript "$ROOT/scripts/render_examples.R"
"$PYTHON_BIN" "$ROOT/scripts/dense_recipes.py" --kind program_nes --input "$ROOT/examples/program_nes_display.tsv" --output "$ROOT/previews/python_program_nes"
"$PYTHON_BIN" "$ROOT/scripts/dense_recipes.py" --kind calibration --input "$ROOT/examples/permutation_null.tsv" --tests "$ROOT/examples/permutation_tests.tsv" --output "$ROOT/previews/python_calibration"
"$PYTHON_BIN" "$ROOT/scripts/dense_recipes.py" --kind locus --input "$ROOT/examples/shared_locus.tsv" --output "$ROOT/previews/python_locus"
bash "$ROOT/run_mechanism_examples.sh"
"$PYTHON_BIN" "$ROOT/scripts/update_previews.py"
