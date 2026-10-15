#!/usr/bin/env python3
"""Download or freeze OBOX and PABPN1 perturbation references into tidy tables."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import shutil
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


OBOX_URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs41586-023-06428-3/MediaObjects/"
    "41586_2023_6428_MOESM3_ESM.xlsx"
)
PABPN1_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174032/suppl/"
    "GSE174032_Pabpn1.FPKM.normbyERCC.txt.gz"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--obox-xlsx", type=Path)
    parser.add_argument("--pabpn-gz", type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def freeze_source(source: Path | None, url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source is None:
        with urllib.request.urlopen(url, timeout=120) as response, destination.open("wb") as output:
            shutil.copyfileobj(response, output)
    else:
        source = source.resolve()
        if not source.is_file():
            raise FileNotFoundError(source)
        if source != destination.resolve():
            shutil.copy2(source, destination)


def clean_symbols(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip().replace({"": pd.NA})


def prepare_obox(path: Path, table_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_excel(path, sheet_name="Table S6", header=None, engine="openpyxl")
    columns = [
        "gene_symbol", "source_base_mean", "source_log2fc", "source_pvalue",
        "source_padj", "source_logp", "author_group", "author_class",
    ]
    pieces: list[pd.DataFrame] = []
    exclusions: list[pd.DataFrame] = []
    for stage, start in (("E2C", 0), ("L2C", 9)):
        frame = raw.iloc[3:, start : start + 8].copy()
        frame.columns = columns
        frame["gene_symbol"] = clean_symbols(frame["gene_symbol"])
        frame = frame[frame["gene_symbol"].notna()].copy()
        for name in ("source_base_mean", "source_log2fc", "source_pvalue", "source_padj", "source_logp"):
            frame[name] = pd.to_numeric(frame[name], errors="coerce")
        frame["author_group"] = frame["author_group"].astype("string").str.strip()
        frame = frame[frame["author_group"].isin(["up-regulated", "down-regulated"])].copy()
        frame["source_stage"] = stage
        frame["direction"] = frame["author_group"].map(
            {"up-regulated": "up_in_mzKO", "down-regulated": "down_in_mzKO"}
        )
        frame["signature_id"] = "OBOX_mzKO_" + stage + "_" + frame["direction"].str.replace("_in_mzKO", "", regex=False)
        frame = frame.sort_values(
            ["gene_symbol", "source_padj", "source_log2fc"],
            ascending=[True, True, False],
            key=lambda x: x.abs() if x.name == "source_log2fc" else x,
            kind="mergesort",
        )
        duplicated = frame.duplicated("gene_symbol", keep="first")
        if duplicated.any():
            duplicate_rows = frame.loc[duplicated, ["signature_id", "source_stage", "direction", "gene_symbol"]].copy()
            duplicate_rows["reason"] = "duplicate_gene_symbol_lower_priority_source_row"
            exclusions.append(duplicate_rows)
        pieces.append(frame.loc[~duplicated])
    result = pd.concat(pieces, ignore_index=True)
    result = result[
        ["signature_id", "source_stage", "direction", "gene_symbol", "source_log2fc",
         "source_padj", "source_base_mean", "author_group", "author_class"]
    ].sort_values(["source_stage", "direction", "gene_symbol"], kind="mergesort")
    excluded = (
        pd.concat(exclusions, ignore_index=True)
        if exclusions
        else pd.DataFrame(columns=["signature_id", "source_stage", "direction", "gene_symbol", "reason"])
    )
    result.to_csv(table_dir / "obox_signatures.tsv", sep="\t", index=False)
    excluded.to_csv(table_dir / "obox_signature_exclusions.tsv", sep="\t", index=False)
    return result, excluded


def prepare_pabpn1(path: Path, table_dir: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t", compression="gzip")
    required = {
        "tracking_id", "PaKD.2C.1", "PaKD.2C.2", "PaWT.2C.1", "PaWT.2C.2",
        "PaKD.8C.1", "PaKD.8C.2", "PaWT.8C.1", "PaWT.8C.2",
    }
    missing = sorted(required.difference(raw.columns))
    if missing:
        raise ValueError(f"PABPN1 matrix is missing columns: {missing}")
    raw["tracking_id"] = clean_symbols(raw["tracking_id"])
    raw = raw[raw["tracking_id"].notna()].copy()
    records: list[pd.DataFrame] = []
    for source_stage, project_stage in (("2C", "L2C"), ("8C", "8C")):
        kd = [f"PaKD.{source_stage}.1", f"PaKD.{source_stage}.2"]
        wt = [f"PaWT.{source_stage}.1", f"PaWT.{source_stage}.2"]
        frame = pd.DataFrame({
            "stage": project_stage,
            "source_stage": source_stage,
            "gene_symbol": raw["tracking_id"],
            "mean_KD": raw[kd].mean(axis=1),
            "mean_WT": raw[wt].mean(axis=1),
        })
        frame["pass_expression"] = frame[["mean_KD", "mean_WT"]].max(axis=1) >= 1.0
        frame["effect_eps1"] = np.log2((frame["mean_KD"] + 1.0) / (frame["mean_WT"] + 1.0))
        frame["effect_eps0_1"] = np.log2((frame["mean_KD"] + 0.1) / (frame["mean_WT"] + 0.1))
        records.append(frame)
    result = pd.concat(records, ignore_index=True).sort_values(["stage", "gene_symbol"], kind="mergesort")
    result.to_csv(table_dir / "pabpn1_effects.tsv", sep="\t", index=False)
    return result


def main() -> None:
    args = parse_args()
    output = args.output.resolve()
    raw_dir = output / "raw"
    table_dir = output / "tables"
    raw_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    obox_path = raw_dir / "OBOX_Nature2023_MOESM3.xlsx"
    pabpn_path = raw_dir / "PABPN1_GSE174032_FPKM_normbyERCC.txt.gz"
    freeze_source(args.obox_xlsx, OBOX_URL, obox_path)
    freeze_source(args.pabpn_gz, PABPN1_URL, pabpn_path)
    obox, obox_exclusions = prepare_obox(obox_path, table_dir)
    pabpn1 = prepare_pabpn1(pabpn_path, table_dir)
    retrieved = datetime.now(timezone.utc).isoformat()
    manifest = pd.DataFrame([
        {
            "reference_id": "OBOX_mzKO_Nature_2023", "accession": "GSE215813",
            "doi": "10.1038/s41586-023-06428-3", "source_url": OBOX_URL,
            "raw_file": str(obox_path.relative_to(output)), "sha256": sha256(obox_path),
            "derived_table": "tables/obox_signatures.tsv", "derived_sha256": sha256(table_dir / "obox_signatures.tsv"),
            "retrieved_utc": retrieved,
            "role": "Author-classified mzKO up/down signatures from Table S6; exact gene-symbol mapping.",
        },
        {
            "reference_id": "PABPN1_depletion_NAR_2022", "accession": "GSE174032",
            "doi": "10.1093/nar/gkab1213", "source_url": PABPN1_URL,
            "raw_file": str(pabpn_path.relative_to(output)), "sha256": sha256(pabpn_path),
            "derived_table": "tables/pabpn1_effects.tsv", "derived_sha256": sha256(table_dir / "pabpn1_effects.tsv"),
            "retrieved_utc": retrieved,
            "role": "ERCC-normalized FPKM effect-size reference; source 2C is mapped to project L2C; no external P/FDR inferred.",
        },
    ])
    manifest.to_csv(table_dir / "reference_manifest.tsv", sep="\t", index=False)
    summary = pd.DataFrame([
        {"metric": "obox_signature_rows", "value": len(obox)},
        {"metric": "obox_signature_count", "value": obox["signature_id"].nunique()},
        {"metric": "obox_duplicate_rows_excluded", "value": len(obox_exclusions)},
        {"metric": "pabpn1_stage_gene_rows", "value": len(pabpn1)},
        {"metric": "pabpn1_expression_eligible_rows", "value": int(pabpn1["pass_expression"].sum())},
    ])
    summary.to_csv(table_dir / "preparation_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
