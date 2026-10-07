#!/usr/bin/env bash
# 仅在无既有结果的独立副本运行；fixture 是合成 methylKit 文本，不是比对产物。
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
for folder in 04methykit result metadata; do [[ ! -e "$ROOT_DIR/$folder" ]] || { echo "示例要求全新副本: $folder 已存在" >&2; exit 2; }; done
export WGBS_INPUT_MODE=cached
export WGBS_SAMPLES="$ROOT_DIR/examples/synthetic_core/samples.tsv"
export WGBS_GTF="$ROOT_DIR/examples/synthetic_core/synthetic.gtf"
source "$ROOT_DIR/scripts/load_config.sh"
# 显式登记随后加载的合成 staged 输入；没有声明 FASTQ/比对已执行。
python3 "$ROOT_DIR/scripts/input_audit.py" --root "$ROOT_DIR" --prefix WGBS --record upstream
mkdir -p "$ROOT_DIR/04methykit"
cp "$ROOT_DIR"/examples/synthetic_core/*_CpG.methylKit "$ROOT_DIR/04methykit/"
printf 'input_layer\torigin\texecution\nCpG_methylKit\tsynthetic_300_loci_4_samples\tfixture_loaded\nFASTQ_alignment_extraction\tnone\tnot_run\n' > "$ROOT_DIR/metadata/synthetic_input_origin.tsv"
bash "$ROOT_DIR/run_upstream.sh" --steps 12,13
bash "$ROOT_DIR/run_wgbs.sh" --steps 1,3
