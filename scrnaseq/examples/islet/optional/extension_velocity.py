"""BAM-derived, beta-only splice counting and separately gated state dynamics.

All generated BAMs are below the extension run. The parent run is read only.
"""
from __future__ import annotations

import gzip
import json
import os
import shutil
import subprocess
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from .runtime import atomic_json, atomic_text, atomic_tsv, run_command, sha256


REPEAT_ALIAS_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chromAlias.txt"


def _contigs(path: Path) -> set[str]:
    values: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t", 1)
            if len(fields) != 2 or not fields[0]:
                raise ValueError(f"Malformed tab-delimited reference row: {path}")
            values.add(fields[0])
    if not values:
        raise ValueError(f"Reference contains no contigs: {path}")
    return values


def normalize_repeat_mask(source: Path, alias_path: Path, fasta_fai: Path, gtf: Path, target: Path) -> dict:
    """Translate the UCSC mm39 mask to the exact Ensembl/Cell Ranger contigs."""
    mapping: dict[str, str] = {}
    reverse: dict[str, str] = {}
    with alias_path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or not fields[0] or not fields[2]:
                raise ValueError("Malformed official mm39 chromosome alias row")
            ucsc, ensembl = fields[0], fields[2]
            if ucsc in mapping or ensembl in reverse:
                raise ValueError("Official chromosome aliases must be one-to-one")
            mapping[ucsc] = ensembl
            reverse[ensembl] = ucsc
    if not mapping:
        raise ValueError("Official chromosome alias table is empty")

    reference_contigs = _contigs(fasta_fai)
    gtf_contigs = _contigs(gtf)
    if set(reverse) != reference_contigs:
        raise ValueError("Official Ensembl aliases do not exactly cover the Cell Ranger FASTA")

    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + ".partial")
    temporary.unlink(missing_ok=True)
    source_contigs: set[str] = set()
    normalized_contigs: set[str] = set()
    try:
        with source.open(encoding="utf-8") as src, temporary.open("w", encoding="utf-8") as dest:
            for line in src:
                if line.startswith("#") or not line.strip():
                    dest.write(line)
                    continue
                fields = line.split("\t", 1)
                if len(fields) != 2 or fields[0] not in mapping:
                    raise ValueError(f"Repeat-mask contig lacks an official mm39 alias: {fields[0]}")
                source_contigs.add(fields[0])
                normalized_contigs.add(mapping[fields[0]])
                dest.write(mapping[fields[0]] + "\t" + fields[1])
        if source_contigs != set(mapping):
            raise ValueError("UCSC repeat mask and official mm39 alias contigs differ")
        if normalized_contigs != reference_contigs or not gtf_contigs <= normalized_contigs:
            raise ValueError("Normalized repeat mask does not cover the Cell Ranger FASTA/GTF contigs")
        temporary.replace(target)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return {"status": "PASS", "official_alias_url": REPEAT_ALIAS_URL,
            "source_mask": str(source.resolve()), "source_mask_sha256": sha256(source),
            "alias": str(alias_path.resolve()), "alias_sha256": sha256(alias_path),
            "fasta_fai": str(fasta_fai.resolve()), "fasta_fai_sha256": sha256(fasta_fai),
            "gtf": str(gtf.resolve()), "gtf_sha256": sha256(gtf),
            "normalized_mask": str(target.resolve()), "normalized_mask_sha256": sha256(target),
            "source_contigs": len(source_contigs), "normalized_contigs": len(normalized_contigs),
            "gtf_contigs": len(gtf_contigs)}


def barcode_from_global(value: str, library: str) -> str:
    prefix = library + "|"
    if not value.startswith(prefix):
        raise ValueError(f"Barcode/library mismatch: {value}, {library}")
    barcode = value[len(prefix):]
    if not barcode.endswith("-1") or not set(barcode[:-2]) <= set("ACGT"):
        raise ValueError(f"Unexpected 10x barcode: {barcode}")
    return barcode


def loom_barcode(value: str, library: str) -> str:
    barcode = value.split(":", 1)[-1]
    if barcode.endswith("x"):
        barcode = barcode[:-1] + "-1"
    if not barcode.endswith("-1"):
        raise ValueError(f"Unexpected loom barcode: {value}")
    return library + "|" + barcode


def _fingerprints(paths: list[Path]) -> dict:
    return {"files": [{"path": str(p.resolve()), "sha256": sha256(p)} for p in paths]}


def _verify_fingerprints(payload: dict) -> None:
    for item in payload["files"]:
        if sha256(Path(item["path"])) != item["sha256"]:
            raise RuntimeError(f"Execution input changed: {item['path']}")


def count_velocity(input_h5ad: Path, source_run: Path, output_dir: Path, reference: Path) -> dict:
    """Use independently locked velocyto, uint32 counts, and exact Cell Ranger GTF."""
    root = Path(__file__).resolve().parents[2]
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    env_root = root / "environments/velocity"
    counter = env_root / ".venv/bin/velocyto"
    samtools = shutil.which("samtools")
    if samtools is None:
        raise FileNotFoundError("samtools must be explicitly available on PATH")
    if not counter.is_file():
        run_command(["uv", "sync", "--frozen", "--project", env_root])
    source_gtf = reference / "genes/genes.gtf.gz"
    source_mask = reference.parent / "ucsc_repeatmasker/mm39.repeat_mask.gtf"
    alias = root / "reference/velocity/mm39.chromAlias.tsv"
    fasta_fai = reference / "fasta/genome.fa.fai"
    frozen = _fingerprints([input_h5ad, source_gtf, source_mask, alias, fasta_fai, env_root / "uv.lock"])
    frozen["started_at"] = datetime.now().isoformat()
    manifest = output_dir / "counting_inputs.json"
    if manifest.exists():
        old = json.loads(manifest.read_text())
        if old["files"] != frozen["files"]:
            raise RuntimeError("Velocity counting inputs changed; use a new run")
    else:
        atomic_json(manifest, frozen)
    source_hashes = pd.read_csv(source_run / "00_audit/output_hashes.tsv", sep="\t").set_index("path")
    backed = ad.read_h5ad(input_h5ad, backed="r")
    obs = backed.obs.loc[backed.obs.cell_type.eq("beta")].copy()
    backed.file.close()
    if obs.empty:
        raise ValueError("No beta cells")
    gtf = output_dir / "cellranger.genes.gtf"
    if not gtf.exists():
        temporary = gtf.with_suffix(".gtf.tmp")
        with gzip.open(source_gtf, "rb") as src, temporary.open("wb") as dest:
            shutil.copyfileobj(src, dest, 16 * 1024 * 1024)
        temporary.replace(gtf)
    mask = output_dir / "cellranger.repeat_mask.gtf"
    mask_receipt = normalize_repeat_mask(source_mask, alias, fasta_fai, gtf, mask)
    atomic_json(output_dir / "repeat_mask_preparation.json", mask_receipt)
    env = os.environ.copy()
    env.update({"PATH": str(Path(samtools).parent) + ":" + env.get("PATH", ""),
                "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "4"})
    rows = []
    for library, cells in obs.groupby("library_id", observed=True):
        library = str(library)
        folder = output_dir / library
        folder.mkdir(exist_ok=True)
        barcodes = folder / "beta_barcodes.tsv"
        atomic_text(barcodes, "\n".join(barcode_from_global(x, library) for x in cells.index) + "\n")
        relative_bam = f"01_cellranger/{library}/outs/possorted_genome_bam.bam"
        bam = source_run / relative_bam
        expected = source_hashes.loc[relative_bam]
        if bam.stat().st_size != int(expected["bytes"]):
            raise ValueError(f"Parent BAM size changed: {bam}")
        # Actual BAM digest is obtained once here, not once per symlink/output alias.
        receipt = folder / "bam_input.json"
        stat = bam.stat()
        identity = {"path": str(bam.resolve()), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns,
                    "sha256": str(expected["sha256"])}
        if receipt.exists():
            if json.loads(receipt.read_text()) != identity:
                raise RuntimeError("BAM identity changed")
        else:
            if sha256(bam) != identity["sha256"]:
                raise ValueError("BAM checksum differs from parent audit")
            atomic_json(receipt, identity)
        filtered = folder / "beta.coordinate.bam"
        filter_receipt = folder / "filter_completed.json"
        if not filter_receipt.exists():
            if shutil.disk_usage(output_dir).free < 100 * 1024**3:
                raise RuntimeError("Velocity counting requires >=100 GiB free workspace")
            temp_bam = folder / "beta.coordinate.partial.bam"
            run_command([samtools, "view", "-@", "4", "-b", "-D", f"CB:{barcodes}",
                         "-o", temp_bam, bam], env=env, log_file=folder / "filter.log")
            run_command([samtools, "quickcheck", "-v", temp_bam])
            temp_bam.replace(filtered)
            atomic_json(filter_receipt, {"barcodes_sha256": sha256(barcodes),
                                         "filtered_bytes": filtered.stat().st_size})
        else:
            state = json.loads(filter_receipt.read_text())
            if state["barcodes_sha256"] != sha256(barcodes) or filtered.stat().st_size != state["filtered_bytes"]:
                raise RuntimeError("Filtered BAM resume signature mismatch")
            run_command([samtools, "quickcheck", "-v", filtered])
        loom = folder / f"{library}.loom"
        completion = folder / "count_completed.json"
        if completion.exists():
            state = json.loads(completion.read_text())
            if not loom.exists() or sha256(loom) != state["loom_sha256"]:
                raise RuntimeError("Completed loom has changed")
        else:
            # A partial cell-sorted BAM must never be reused after interruption.
            sorted_bam = folder / "cellsorted_beta.coordinate.bam"
            if sorted_bam.exists():
                run_command([samtools, "quickcheck", "-v", sorted_bam])
            run_command([counter, "run", "-b", barcodes, "-o", folder, "-e", library,
                         "-m", mask, "-@", "4", "--samtools-memory", "4096", "-t", "uint32",
                         filtered, gtf], env=env, log_file=folder / "velocyto.log")
            if not loom.exists():
                raise RuntimeError(f"velocyto returned without {loom}")
            import loompy
            with loompy.connect(str(loom), mode="r") as ds:
                if not {"spliced", "unspliced", "ambiguous"} <= set(ds.layers.keys()):
                    raise ValueError("Incomplete loom layers")
                if ds.shape[1] != len(cells):
                    raise ValueError("Whitelist/loom cell count mismatch")
            atomic_json(completion, {"loom_sha256": sha256(loom), "cells": len(cells),
                                      "finished_at": datetime.now().isoformat()})
        rows.append({"library_id": library, "beta_cells": len(cells), "loom": str(loom),
                     "status": "completed", "counts_dtype": "uint32"})
    _verify_fingerprints(frozen)
    atomic_tsv(output_dir / "counting_status.tsv", pd.DataFrame(rows))
    return {"status": "completed", "libraries": rows}
