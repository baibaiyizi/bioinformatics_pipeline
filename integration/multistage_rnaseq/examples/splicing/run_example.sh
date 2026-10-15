#!/usr/bin/env bash
set -euo pipefail

example_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
template_root="$(cd "$example_dir/../.." && pwd)"

Rscript "$template_root/scripts/splicing_trajectories.R" \
  --config-dir "$example_dir" \
  --rnaseq-root "$example_dir/input" \
  --output "$example_dir/result"

Rscript - "$example_dir/result" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
root <- args[[1L]]
stats <- read.delim(file.path(root, "tables", "target_event_stage_statistics.tsv"), check.names = FALSE)
psi <- read.delim(file.path(root, "tables", "target_event_psi_long.tsv"), check.names = FALSE, na.strings = "NA")
validation <- read.delim(file.path(root, "tables", "validation.tsv"), check.names = FALSE)
sashimi <- read.delim(file.path(root, "tables", "sashimi_status.tsv"), check.names = FALSE)
stopifnot(abs(stats$delta_psi_comparison_minus_reference[stats$stage == "S1" & stats$count_mode == "JC"] - 0.55) < 1e-9)
stopifnot(abs(stats$delta_psi_comparison_minus_reference[stats$stage == "S2" & stats$count_mode == "JC"] + 0.2) < 1e-9)
stopifnot(sum(is.na(psi$psi)) == 1L)
stopifnot(any(grepl("partial keys map", validation$detail)))
stopifnot(identical(sashimi$status[[1L]], "missing_bam"))
cat("example_validation=PASS\n")
RS
