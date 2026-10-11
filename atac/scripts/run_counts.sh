#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT_DIR/scripts/load_config.sh"
CHECK=0
REPORT_ONLY=0
while (( $# )); do
  case "$1" in
    --help|-h) echo 'counts 模式：原始整数 consensus fragments → DESeq2 → 差异表/诊断图/报告。--check/--dry-run 只读；--steps 5、--all 或无参数执行；--report-only 刷新报告。'; exit 0 ;;
    --list) echo '05 consensus_counts_DESeq2（原始区间）+ HTML report；BAM/DiffBind/注释/footprint 不在 counts 模式运行'; exit 0 ;;
    --check|--dry-run) CHECK=1; shift ;;
    --all) shift ;;
    --report-only) REPORT_ONLY=1; shift ;;
    --steps) [[ $# -ge 2 && "$2" =~ ^(5|05|diff|counts)$ ]] || { echo 'counts 模式只接受 --steps 5'; exit 2; }; shift 2 ;;
    *) echo "counts 模式不支持: $1" >&2; exit 2 ;;
  esac
done
python3 "$ROOT_DIR/scripts/input_audit.py" --root "$ROOT_DIR" --prefix ATAC
Rscript --vanilla "$ROOT_DIR/scripts/run_counts.R" --check
if [[ "$CHECK" == 1 ]]; then exit 0; fi
if [[ "$REPORT_ONLY" == 1 ]]; then
 python3 "$ROOT_DIR/scripts/input_audit.py" --root "$ROOT_DIR" --prefix ATAC --verify downstream
 python3 - "$ROOT_DIR/result/execution.tsv" <<'PYCHECK'
import csv, sys
with open(sys.argv[1], newline='') as handle:
    rows = [r for r in csv.DictReader(handle, delimiter='\t') if r['module'] == '05_consensus_counts_DESeq2']
if not rows or rows[-1]['status'] != 'COMPLETED':
    sys.exit('UNVERIFIED: counts 模块最新状态不是 COMPLETED，拒绝将旧表作为当前报告输入')
PYCHECK
else
 python3 "$ROOT_DIR/scripts/input_audit.py" --root "$ROOT_DIR" --prefix ATAC --record downstream
 mkdir -p "$ROOT_DIR/result"
 [[ -f "$ROOT_DIR/result/execution.tsv" ]] || printf 'time\tmodule\tstatus\n' > "$ROOT_DIR/result/execution.tsv"
 printf '%s\t05_consensus_counts_DESeq2\tSTARTED\n' "$(date -Is)" >> "$ROOT_DIR/result/execution.tsv"
 CURRENT_STEP=05_consensus_counts_DESeq2
 trap 'code=$?; if [[ $code != 0 ]]; then printf "%s\t%s\tFAILED\n" "$(date -Is)" "$CURRENT_STEP" >> "$ROOT_DIR/result/execution.tsv"; fi' EXIT
 Rscript --vanilla "$ROOT_DIR/scripts/run_counts.R"
 printf '%s\t05_consensus_counts_DESeq2\tCOMPLETED\n' "$(date -Is)" >> "$ROOT_DIR/result/execution.tsv"
fi
CURRENT_STEP=report
printf '%s\treport\tSTARTED\n' "$(date -Is)" >> "$ROOT_DIR/result/execution.tsv"
trap 'code=$?; if [[ $code != 0 ]]; then printf "%s\treport\tFAILED\n" "$(date -Is)" >> "$ROOT_DIR/result/execution.tsv"; fi' EXIT
RESULT_REPORT_ROOT="$ROOT_DIR" Rscript --vanilla "$ROOT_DIR/rmd/generate_result_report.R"
printf '%s\treport\tCOMPLETED\n' "$(date -Is)" >> "$ROOT_DIR/result/execution.tsv"
