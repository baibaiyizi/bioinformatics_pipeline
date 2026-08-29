from __future__ import annotations

import math
from typing import Any

import networkx as nx
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from .contract import ContractError


LANES = ["disease_association", "ppi_proximity", "pathway_support", "structure_availability", "direct_compound_target"]


def build_multidigraph(nodes: pd.DataFrame, edges: pd.DataFrame) -> nx.MultiDiGraph:
    if nodes["node_id"].duplicated().any():
        raise ContractError("标准实体 node_id 重复")
    evidence_key = ["source", "target", "relation", "database", "database_version", "evidence_id", "evidence_layer"]
    if edges[evidence_key].astype(str).duplicated().any():
        raise ContractError("同一版本的证据边重复，禁止简单累加")
    known = set(nodes["node_id"])
    unknown = (set(edges["source"]) | set(edges["target"])) - known
    if unknown:
        raise ContractError(f"证据边引用未知实体: {sorted(unknown)[:5]}")
    graph = nx.MultiDiGraph(project="SCI network pharmacology", claim_boundary="prioritization_only")
    for row in nodes.fillna("").to_dict("records"):
        node_id = str(row.pop("node_id"))
        graph.add_node(node_id, **{k: str(v) for k, v in row.items()})
    for index, row in enumerate(edges.fillna("").to_dict("records")):
        source, target = str(row.pop("source")), str(row.pop("target"))
        key = f"e{index:07d}"
        graph.add_edge(source, target, key=key, **{k: (float(v) if k == "score" and v != "" else str(v)) for k, v in row.items()})
    return graph


def rank_candidates(targets: pd.DataFrame, lane_scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if targets.empty:
        raise ContractError("没有可排序 target")
    if lane_scores.empty:
        lane_scores = pd.DataFrame(columns=["target_id", "lane", "raw_score"])
    unknown = set(lane_scores.get("lane", [])) - set(LANES)
    if unknown:
        raise ContractError(f"未知 evidence lane: {sorted(unknown)}")
    if lane_scores.duplicated(["target_id", "lane"]).any():
        raise ContractError("同一 target/lane 只能有一个预先计算的分数")
    wide = targets[["target_id", "symbol", "evidence_layer"]].drop_duplicates("target_id").set_index("target_id")
    raw = lane_scores.pivot(index="target_id", columns="lane", values="raw_score") if not lane_scores.empty else pd.DataFrame(index=wide.index)
    raw = raw.reindex(index=wide.index, columns=LANES)
    percentiles = pd.DataFrame(index=raw.index)
    for lane in LANES:
        observed = raw[lane].dropna()
        percentiles[f"{lane}_percentile"] = raw[lane].rank(method="average", pct=True) if not observed.empty else math.nan
    available = percentiles.notna()
    wide["evidence_lane_count"] = available.sum(axis=1)
    wide["eligible"] = wide["evidence_lane_count"] >= 2
    wide["priority_score"] = percentiles.mean(axis=1, skipna=True).where(wide["eligible"])
    result = pd.concat([wide, raw.add_suffix("_raw"), percentiles], axis=1).reset_index()
    result = result.sort_values(["eligible", "priority_score", "target_id"], ascending=[False, False, True], na_position="last")

    loo_rows: list[dict[str, Any]] = []
    for target_id in percentiles.index:
        for omitted in LANES:
            cols = [f"{lane}_percentile" for lane in LANES if lane != omitted]
            values = percentiles.loc[target_id, cols].dropna()
            loo_rows.append({"target_id": target_id, "omitted_lane": omitted, "remaining_lane_count": len(values), "leave_one_lane_score": float(values.mean()) if len(values) >= 2 else math.nan})
    return result, pd.DataFrame(loo_rows)


def reactome_ora(selected: set[str], background_pathways: dict[str, list[dict[str, Any]]], species_taxid: int) -> pd.DataFrame:
    columns = ["pathway_id", "pathway_name", "species_taxid", "selected_hits", "selected_total", "background_hits", "background_total", "p_value", "fdr_bh"]
    if not background_pathways:
        raise ContractError("没有显式 Reactome 背景宇宙，禁止富集")
    background = set(background_pathways)
    if not selected <= background:
        raise ContractError("Reactome selected IDs 必须全部属于显式背景宇宙")
    prefix = "R-HSA-" if species_taxid == 9606 else "R-MMU-"
    membership: dict[str, set[str]] = {}
    names: dict[str, str] = {}
    for accession, pathways in background_pathways.items():
        for pathway in pathways:
            pathway_id = str(pathway.get("stId") or "")
            if not pathway_id.startswith(prefix):
                continue
            membership.setdefault(pathway_id, set()).add(accession)
            names[pathway_id] = str(pathway.get("displayName") or pathway.get("name") or pathway_id)
    rows = []
    for pathway_id, members in membership.items():
        hits = selected & members
        if not hits:
            continue
        p_value = hypergeom.sf(len(hits) - 1, len(background), len(members), len(selected))
        rows.append({"pathway_id": pathway_id, "pathway_name": names[pathway_id], "species_taxid": species_taxid, "selected_hits": len(hits), "selected_total": len(selected), "background_hits": len(members), "background_total": len(background), "p_value": float(p_value)})
    if not rows:
        return pd.DataFrame(columns=columns)
    fdr = multipletests([row["p_value"] for row in rows], method="fdr_bh")[1]
    for row, adjusted in zip(rows, fdr, strict=True):
        row["fdr_bh"] = float(adjusted)
    return pd.DataFrame(rows)[columns].sort_values(["fdr_bh", "p_value", "pathway_id"])


def string_threshold_sensitivity(edges: pd.DataFrame, targets: pd.DataFrame, thresholds: list[int]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    symbols = set(targets["symbol"].dropna().astype(str))
    for network_type in ("physical", "functional"):
        subset = edges[(edges["database"] == "STRING") & (edges["network_type"] == network_type)]
        for threshold in thresholds:
            kept = subset[subset["score"] * 1000 >= threshold]
            degree = {symbol: 0 for symbol in symbols}
            for row in kept.to_dict("records"):
                source_symbol = str(row.get("source_symbol", ""))
                target_symbol = str(row.get("target_symbol", ""))
                if source_symbol in degree:
                    degree[source_symbol] += 1
                if target_symbol in degree:
                    degree[target_symbol] += 1
            for symbol, value in degree.items():
                rows.append({"network_type": network_type, "required_score": threshold, "symbol": symbol, "degree": value})
    return pd.DataFrame(rows)
