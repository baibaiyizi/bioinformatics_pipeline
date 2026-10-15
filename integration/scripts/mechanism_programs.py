#!/usr/bin/env python3
"""Compare leading-edge membership across independently analysed RNA contexts.

This template consumes completed DE/GSEA outputs.  It does not run differential
expression or GSEA.  Overlap tests use pathway members with a finite DE statistic
in both contexts as their universe.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


DE_COLUMNS = {"gene_id", "SYMBOL", "log2FoldChange", "stat", "pvalue", "padj"}
GSEA_COLUMNS = {
    "pathway", "pval", "padj_context_collection", "NES", "size",
    "leadingEdge", "context", "collection",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_context(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("--context must be NAME=DE_TABLE")
    name, path = value.split("=", 1)
    name = name.strip()
    if not name:
        raise argparse.ArgumentTypeError("context name cannot be empty")
    return name, Path(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fgsea", required=True, type=Path)
    parser.add_argument("--gene-sets-rds", required=True, type=Path)
    parser.add_argument("--context", required=True, action="append", type=parse_context,
                        help="repeat NAME=DE_TABLE at least twice")
    parser.add_argument("--collection", required=True, help="list name inside bundle$libraries")
    parser.add_argument("--comparison-label", required=True, help="effect direction, e.g. treated - control")
    parser.add_argument("--assembly", required=True, help="reference assembly shared by all contexts")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-pathways", type=int, default=0,
                        help="optional exact pathway-count gate; 0 disables it")
    return parser.parse_args()


def export_members(rds_path: Path, collection: str) -> pd.DataFrame:
    code = r'''
args <- commandArgs(trailingOnly=TRUE)
x <- readRDS(args[[1]])
if (is.null(x$libraries[[args[[2]]]])) stop("collection absent from bundle: ", args[[2]])
sets <- x$libraries[[args[[2]]]]
rows <- do.call(rbind, lapply(names(sets), function(nm) {
  data.frame(pathway=nm, gene=unique(as.character(sets[[nm]])), stringsAsFactors=FALSE)
}))
write.table(rows, args[[3]], sep="\t", quote=FALSE, row.names=FALSE, na="")
'''
    with tempfile.TemporaryDirectory(prefix="program-turnover-") as temp_dir:
        target = Path(temp_dir) / "members.tsv"
        subprocess.run(
            ["Rscript", "-e", code, str(rds_path), collection, str(target)],
            check=True, capture_output=True, text=True,
        )
        return pd.read_csv(target, sep="\t", dtype=str).drop_duplicates()


def load_de(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, low_memory=False)
    missing = DE_COLUMNS - set(raw.columns)
    if missing:
        raise ValueError(f"{path} missing DE columns: {sorted(missing)}")
    for column in ("stat", "log2FoldChange", "pvalue", "padj"):
        raw[column] = pd.to_numeric(raw[column], errors="coerce")
    raw["SYMBOL"] = raw["SYMBOL"].astype("string").str.strip()
    raw = raw.loc[raw["SYMBOL"].notna() & raw["SYMBOL"].ne("") & np.isfinite(raw["stat"])].copy()
    raw["abs_stat"] = raw["stat"].abs()
    raw = raw.sort_values(["SYMBOL", "abs_stat", "gene_id"], ascending=[True, False, True])
    return raw.drop_duplicates("SYMBOL", keep="first").set_index("SYMBOL", drop=False)


def sign(value: object) -> int | None:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(value) or value == 0:
        return None
    return 1 if value > 0 else -1


def adjust_bh(values: pd.Series) -> np.ndarray:
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(float)
    out = np.full(len(numeric), np.nan)
    keep = np.isfinite(numeric)
    if keep.any():
        out[keep] = multipletests(numeric[keep], method="fdr_bh")[1]
    return out


def main() -> None:
    args = parse_args()
    fgsea_path = args.fgsea.resolve(strict=True)
    rds_path = args.gene_sets_rds.resolve(strict=True)
    contexts = [(name, path.resolve(strict=True)) for name, path in args.context]
    context_names = [name for name, _ in contexts]
    if len(context_names) < 2 or len(set(context_names)) != len(context_names):
        raise ValueError("provide at least two uniquely named contexts")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    fgsea = pd.read_csv(fgsea_path, sep="\t", low_memory=False)
    missing = GSEA_COLUMNS - set(fgsea.columns)
    if missing:
        raise ValueError(f"{fgsea_path} missing GSEA columns: {sorted(missing)}")
    effects = fgsea.loc[
        fgsea["collection"].eq(args.collection) & fgsea["context"].isin(context_names)
    ].copy()
    numeric_gsea = ["pval", "padj_context_collection", "NES", "size"]
    for column in numeric_gsea:
        effects[column] = pd.to_numeric(effects[column], errors="coerce")
    if not np.isfinite(effects[numeric_gsea].to_numpy(float)).all():
        raise ValueError("selected GSEA rows require finite pval, padj_context_collection, NES, and size")
    if not effects["pval"].between(0, 1).all() or not effects["padj_context_collection"].between(0, 1).all():
        raise ValueError("selected GSEA pval and padj_context_collection must be between 0 and 1")
    if not effects["size"].gt(0).all():
        raise ValueError("selected GSEA size must be positive")
    duplicates = effects.duplicated(["context", "pathway"], keep=False)
    if duplicates.any():
        raise ValueError("GSEA table has duplicate context-pathway rows")

    memberships = export_members(rds_path, args.collection)
    pathways = sorted(set(memberships["pathway"]))
    if args.expected_pathways and len(pathways) != args.expected_pathways:
        raise ValueError(f"expected {args.expected_pathways} pathways, observed {len(pathways)}")
    expected_rows = len(pathways) * len(context_names)
    if len(effects) != expected_rows or set(effects["pathway"]) != set(pathways):
        raise ValueError(f"GSEA table must contain every pathway in every context; expected {expected_rows}, observed {len(effects)}")
    effects = effects.sort_values(["pathway", "context"]).reset_index(drop=True)
    effects["raw_p_supported"] = pd.to_numeric(effects["pval"], errors="coerce").lt(0.05)
    effects["bh_supported"] = pd.to_numeric(effects["padj_context_collection"], errors="coerce").lt(0.05)
    effects["comparison_label"] = args.comparison_label
    effects["assembly"] = args.assembly
    effects["inference_scope"] = "independent_context_descriptive_not_longitudinal"

    de_map = {name: load_de(path) for name, path in contexts}
    set_map = memberships.groupby("pathway")["gene"].apply(lambda x: set(x.dropna())).to_dict()
    redundancy_rows: list[dict[str, object]] = []
    for set_a, set_b in itertools.combinations(pathways, 2):
        members_a, members_b = set_map[set_a], set_map[set_b]
        overlap = members_a & members_b
        union = members_a | members_b
        redundancy_rows.append({
            "setA": set_a,
            "setB": set_b,
            "nA": len(members_a),
            "nB": len(members_b),
            "overlap": len(overlap),
            "union": len(union),
            "Jaccard": len(overlap) / len(union) if union else np.nan,
        })
    gene_set_redundancy = pd.DataFrame(redundancy_rows)
    expected_set_pairs = len(pathways) * (len(pathways) - 1) // 2
    if len(gene_set_redundancy) != expected_set_pairs:
        raise RuntimeError(f"expected {expected_set_pairs} frozen gene-set pairs, observed {len(gene_set_redundancy)}")
    effect_lookup = {(row.context, row.pathway): row for row in effects.itertuples(index=False)}
    leading_map: dict[tuple[str, str], set[str]] = {}
    member_rows: list[dict[str, object]] = []
    for context in context_names:
        de = de_map[context]
        for pathway in pathways:
            leading = {x for x in str(effect_lookup[(context, pathway)].leadingEdge).split(";") if x and x != "nan"}
            leading_map[(context, pathway)] = leading
            for gene in sorted(set_map[pathway]):
                row = de.loc[gene] if gene in de.index else None
                member_rows.append({
                    "pathway": pathway, "context": context, "gene": gene,
                    "testable_in_context": row is not None, "leading_edge": gene in leading,
                    "log2FoldChange": np.nan if row is None else row.log2FoldChange,
                    "stat": np.nan if row is None else row.stat,
                    "pvalue": np.nan if row is None else row.pvalue,
                    "padj": np.nan if row is None else row.padj,
                    "gene_id": "" if row is None else row.gene_id,
                })

    pair_rows: list[dict[str, object]] = []
    pair_member_rows: list[dict[str, object]] = []
    for context_a, context_b in itertools.combinations(context_names, 2):
        de_a, de_b = de_map[context_a], de_map[context_b]
        for pathway in pathways:
            members = set_map[pathway]
            tested_a, tested_b = members & set(de_a.index), members & set(de_b.index)
            common = tested_a & tested_b
            le_a_all, le_b_all = leading_map[(context_a, pathway)], leading_map[(context_b, pathway)]
            le_a, le_b = le_a_all & common, le_b_all & common
            shared, union = le_a & le_b, le_a | le_b
            overlap_p = float(hypergeom.sf(len(shared) - 1, len(common), len(le_a), len(le_b))) if common else np.nan
            common_tested = concordant = 0
            for gene in common:
                a, b = sign(de_a.at[gene, "log2FoldChange"]), sign(de_b.at[gene, "log2FoldChange"])
                if a is not None and b is not None:
                    common_tested += 1
                    concordant += int(a == b)
            pair_rows.append({
                "pathway": pathway, "context_a": context_a, "context_b": context_b,
                "pair": f"{context_a}__{context_b}", "pathway_members_n": len(members),
                "testable_a_n": len(tested_a), "testable_b_n": len(tested_b),
                "common_testable_pathway_n": len(common),
                "leading_edge_a_total_n": len(le_a_all), "leading_edge_b_total_n": len(le_b_all),
                "leading_edge_a_common_universe_n": len(le_a), "leading_edge_b_common_universe_n": len(le_b),
                "shared_leading_edge_n": len(shared), "a_only_leading_edge_n": len(le_a - le_b),
                "b_only_leading_edge_n": len(le_b - le_a), "union_leading_edge_n": len(union),
                "jaccard": len(shared) / len(union) if union else np.nan,
                "common_pathway_effect_direction_tested_n": common_tested,
                "common_pathway_effect_direction_concordant_n": concordant,
                "common_pathway_effect_direction_concordance_rate": concordant / common_tested if common_tested else np.nan,
                "NES_a": effect_lookup[(context_a, pathway)].NES,
                "NES_b": effect_lookup[(context_b, pathway)].NES,
                "NES_direction_concordant": sign(effect_lookup[(context_a, pathway)].NES) == sign(effect_lookup[(context_b, pathway)].NES),
                "overlap_pvalue": overlap_p,
                "overlap_null": "hypergeometric_within_pathway_members_testable_in_both_contexts",
                "inference_scope": "independent_context_descriptive_not_longitudinal",
            })
            for gene in sorted(union):
                pair_member_rows.append({
                    "pathway": pathway, "context_a": context_a, "context_b": context_b,
                    "pair": f"{context_a}__{context_b}", "gene": gene,
                    "leading_edge_status": "shared" if gene in shared else f"{context_a}_only" if gene in le_a else f"{context_b}_only",
                    "common_testable_pathway_member": gene in common,
                    "log2FoldChange_a": de_a.at[gene, "log2FoldChange"] if gene in de_a.index else np.nan,
                    "stat_a": de_a.at[gene, "stat"] if gene in de_a.index else np.nan,
                    "pvalue_a": de_a.at[gene, "pvalue"] if gene in de_a.index else np.nan,
                    "padj_a": de_a.at[gene, "padj"] if gene in de_a.index else np.nan,
                    "log2FoldChange_b": de_b.at[gene, "log2FoldChange"] if gene in de_b.index else np.nan,
                    "stat_b": de_b.at[gene, "stat"] if gene in de_b.index else np.nan,
                    "pvalue_b": de_b.at[gene, "pvalue"] if gene in de_b.index else np.nan,
                    "padj_b": de_b.at[gene, "padj"] if gene in de_b.index else np.nan,
                })

    pair_summary = pd.DataFrame(pair_rows)
    pair_summary["overlap_bh_prespecified_family"] = adjust_bh(pair_summary["overlap_pvalue"])
    pair_summary["overlap_raw_p_supported"] = pair_summary["overlap_pvalue"].lt(0.05)
    pair_summary["overlap_bh_supported"] = pair_summary["overlap_bh_prespecified_family"].lt(0.05)
    outputs = {
        "program_context_effects.tsv": effects,
        "program_driver_members.tsv": pd.DataFrame(member_rows),
        "program_pairwise_turnover.tsv": pair_summary,
        "program_pair_driver_members.tsv": pd.DataFrame(pair_member_rows),
        "program_geneset_redundancy.tsv": gene_set_redundancy,
    }
    for name, frame in outputs.items():
        frame.to_csv(output / name, sep="\t", index=False, na_rep="NA")

    receipt = {
        "status": "completed", "generated_at": datetime.now(timezone.utc).isoformat(),
        "comparison_label": args.comparison_label, "assembly": args.assembly,
        "collection": args.collection, "contexts": context_names,
        "input_contract": {"DE_columns": sorted(DE_COLUMNS), "GSEA_columns": sorted(GSEA_COLUMNS)},
        "multiplicity": f"BH over {len(pair_summary)} prespecified pathway-context-pair overlap tests",
        "overlap_universe": "pathway members with finite DE statistic in both contexts",
        "geneset_redundancy": {
            "rows": len(gene_set_redundancy),
            "membership_scope": "all frozen collection members; no DE-testability or leading-edge filter",
            "significance_tested": False,
            "independent_vote_interpretation_allowed": False,
        },
        "input_sha256": {
            "fgsea": {"path": str(fgsea_path), "sha256": sha256(fgsea_path)},
            "gene_sets_rds": {"path": str(rds_path), "sha256": sha256(rds_path)},
            "contexts": {
                name: {"path": str(path), "sha256": sha256(path)} for name, path in contexts
            },
        },
        "code_sha256": sha256(Path(__file__).resolve()),
        "output_sha256": {name: sha256(output / name) for name in outputs},
        "rows": {name: len(frame) for name, frame in outputs.items()},
        "claim_boundary": "Cross-context sharing is descriptive and is not a longitudinal trajectory or causal regulatory chain. Frozen gene-set Jaccard records membership redundancy only; overlapping pathways are not independent votes.",
    }
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
