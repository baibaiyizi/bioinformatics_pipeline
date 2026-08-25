#!/usr/bin/env python3
"""Build a versioned evidence graph and candidate prioritization artifacts."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import pandas as pd
import yaml


class ContractError(ValueError):
    pass


TIER_WEIGHT = {"experimental": 1.0, "curated": 0.9, "genetic": 0.9, "clinical": 1.0, "computational": 0.5, "text_mined": 0.3}


def resolve(base: Path, raw: str) -> Path:
    path = Path(raw).expanduser()
    return path if path.is_absolute() else (base / path).resolve()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    required = {"nodes", "edges", "output_dir"}
    missing = sorted(required - set(config))
    if missing:
        raise ContractError(f"配置缺少字段: {', '.join(missing)}")
    dropout = float(config.get("edge_dropout", 0.1))
    if not 0 <= dropout < 1:
        raise ContractError("edge_dropout 必须位于 [0,1)")
    return config


def load_inputs(nodes_path: Path, edges_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    nodes = pd.read_csv(nodes_path, dtype=str).fillna("")
    edges = pd.read_csv(edges_path)
    node_required = {"node_id", "node_type", "canonical_id", "species"}
    edge_required = {"source", "target", "relation", "evidence_source", "evidence_id", "source_version", "evidence_tier", "confidence"}
    missing_nodes = sorted(node_required - set(nodes))
    missing_edges = sorted(edge_required - set(edges))
    if missing_nodes or missing_edges:
        raise ContractError(f"字段缺失 nodes={missing_nodes}, edges={missing_edges}")
    if nodes["node_id"].eq("").any() or nodes["node_id"].duplicated().any():
        raise ContractError("node_id 必须非空且唯一")
    molecular = nodes["node_type"].isin(["gene", "protein", "target"])
    if nodes.loc[molecular, "species"].eq("").any():
        raise ContractError("gene/protein/target 节点必须显式记录 species/taxon")
    if nodes["canonical_id"].eq("").any():
        raise ContractError("所有节点必须有 canonical_id")
    if edges[list(edge_required - {"confidence"})].astype(str).eq("").any().any():
        raise ContractError("证据边的关系、来源、版本、证据 ID 和层级不得为空")
    edge_nodes = set(edges["source"].astype(str)) | set(edges["target"].astype(str))
    unknown = sorted(edge_nodes - set(nodes["node_id"]))
    if unknown:
        raise ContractError(f"边表引用未知节点: {', '.join(unknown[:10])}")
    edges["confidence"] = pd.to_numeric(edges["confidence"], errors="raise")
    if not edges["confidence"].between(0, 1).all():
        raise ContractError("confidence 必须位于 [0,1]")
    unknown_tiers = sorted(set(edges["evidence_tier"]) - set(TIER_WEIGHT))
    if unknown_tiers:
        raise ContractError(f"未知 evidence_tier: {', '.join(unknown_tiers)}")
    evidence_key = ["source", "target", "relation", "evidence_source", "source_version", "evidence_id"]
    if edges[evidence_key].astype(str).duplicated().any():
        raise ContractError("同一关系中的版本化 evidence 记录重复，禁止重复计权")
    edges["evidence_weight"] = edges["confidence"] * edges["evidence_tier"].map(TIER_WEIGHT)
    return nodes, edges


def build_graph(nodes: pd.DataFrame, edges: pd.DataFrame, directed: bool) -> nx.Graph:
    graph: nx.Graph = nx.DiGraph() if directed else nx.Graph()
    for row in nodes.to_dict("records"):
        graph.add_node(str(row.pop("node_id")), **{k: str(v) for k, v in row.items()})
    collapsed = edges.groupby(["source", "target", "relation"], as_index=False).agg(
        weight=("evidence_weight", "sum"),
        evidence_count=("evidence_id", "nunique"),
        evidence_sources=("evidence_source", lambda x: "|".join(sorted(set(map(str, x))))),
        source_versions=("source_version", lambda x: "|".join(sorted(set(map(str, x))))),
    )
    for row in collapsed.to_dict("records"):
        source, target = str(row.pop("source")), str(row.pop("target"))
        if graph.has_edge(source, target):
            graph[source][target]["weight"] += float(row["weight"])
            graph[source][target]["evidence_count"] += int(row["evidence_count"])
            graph[source][target]["relations"] = graph[source][target]["relations"] + "|" + str(row["relation"])
        else:
            graph.add_edge(source, target, weight=float(row["weight"]), evidence_count=int(row["evidence_count"]), relations=str(row["relation"]), evidence_sources=str(row["evidence_sources"]), source_versions=str(row["source_versions"]))
    return graph


def graph_metrics(graph: nx.Graph) -> pd.DataFrame:
    if graph.number_of_nodes() == 0:
        raise ContractError("证据图为空")
    undirected = graph.to_undirected()
    degree = nx.degree_centrality(undirected)
    between = nx.betweenness_centrality(undirected, weight=None, normalized=True)
    pagerank = nx.pagerank(graph, weight="weight")
    weighted_degree = dict(undirected.degree(weight="weight"))
    return pd.DataFrame({
        "node_id": list(graph.nodes),
        "degree_centrality": [degree[n] for n in graph.nodes],
        "betweenness_centrality": [between[n] for n in graph.nodes],
        "pagerank": [pagerank[n] for n in graph.nodes],
        "weighted_degree": [weighted_degree[n] for n in graph.nodes],
    })


def robustness(graph: nx.Graph, repeats: int, dropout: float, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    nodes = list(graph.nodes)
    edges = list(graph.edges)
    records = []
    for repeat in range(repeats):
        keep = rng.random(len(edges)) >= dropout
        sampled = graph.__class__()
        sampled.add_nodes_from(graph.nodes(data=True))
        sampled.add_edges_from([(*edges[i], graph.edges[edges[i]]) for i in range(len(edges)) if keep[i]])
        scores = dict(sampled.to_undirected().degree(weight="weight"))
        ranking = pd.Series(scores).rank(ascending=False, method="average")
        records.extend({"repeat": repeat, "node_id": node, "weighted_degree": scores[node], "rank": ranking[node]} for node in nodes)
    frame = pd.DataFrame(records)
    return frame.groupby("node_id", as_index=False).agg(rank_median=("rank", "median"), rank_iqr=("rank", lambda x: x.quantile(0.75) - x.quantile(0.25)), weighted_degree_mean=("weighted_degree", "mean"))


def candidate_table(graph: nx.Graph, nodes: pd.DataFrame, metrics: pd.DataFrame, candidate_types: set[str], disease_types: set[str]) -> pd.DataFrame:
    node_types = nodes.set_index("node_id")["node_type"].to_dict()
    diseases = [n for n in graph if node_types.get(n) in disease_types]
    candidates = [n for n in graph if node_types.get(n) in candidate_types]
    undirected = graph.to_undirected()
    rows = []
    for candidate in candidates:
        distances = []
        for disease in diseases:
            try:
                distances.append(nx.shortest_path_length(undirected, candidate, disease))
            except nx.NetworkXNoPath:
                pass
        rows.append({"node_id": candidate, "nearest_disease_steps": min(distances) if distances else np.nan, "reachable_diseases": len(distances)})
    table = pd.DataFrame(rows)
    if table.empty:
        return pd.DataFrame(columns=["node_id", "nearest_disease_steps", "reachable_diseases", "priority_score"])
    table = table.merge(metrics, on="node_id", how="left")
    distance_score = 1 / table["nearest_disease_steps"].replace(0, 1)
    def rescale(series: pd.Series) -> pd.Series:
        span = series.max() - series.min()
        return pd.Series(0.0, index=series.index) if span == 0 else (series - series.min()) / span
    table["priority_score"] = 0.5 * rescale(table["weighted_degree"]) + 0.3 * rescale(distance_score.fillna(0)) + 0.2 * rescale(table["reachable_diseases"].astype(float))
    return table.sort_values("priority_score", ascending=False)


def write_cards(path: Path, candidates: pd.DataFrame, nodes: pd.DataFrame, edges: pd.DataFrame) -> None:
    node_lookup = nodes.set_index("node_id").to_dict("index")
    parts = ["# 候选证据卡", "", "> 网络排名仅用于候选优先级，不等于疗效、因果机制或直接结合证据。", ""]
    for row in candidates.head(30).to_dict("records"):
        node = row["node_id"]
        related = edges[(edges["source"].astype(str) == node) | (edges["target"].astype(str) == node)]
        parts.extend([
            f"## {node}", "",
            f"- 类型：{node_lookup[node]['node_type']}",
            f"- canonical ID：{node_lookup[node]['canonical_id']}",
            f"- 优先级分数：{row['priority_score']:.4f}",
            f"- 最近 disease/phenotype 路径长度：{row['nearest_disease_steps']}",
            f"- 独立 evidence ID 数：{related['evidence_id'].nunique()}",
            f"- 来源版本：{' | '.join(sorted(set(related['source_version'].astype(str))))}",
            "- 解释边界：需要回看逐条证据并用独立实验验证。", "",
        ])
    path.write_text("\n".join(parts), encoding="utf-8")


def run(config_path: Path) -> Path:
    config_path = config_path.resolve()
    config = load_config(config_path)
    nodes_path = resolve(config_path.parent, config["nodes"])
    edges_path = resolve(config_path.parent, config["edges"])
    out = resolve(config_path.parent, config["output_dir"])
    if not nodes_path.is_file() or not edges_path.is_file():
        raise ContractError("nodes 或 edges 输入不存在")
    nodes, edges = load_inputs(nodes_path, edges_path)
    graph = build_graph(nodes, edges, bool(config.get("directed", False)))
    metrics = graph_metrics(graph)
    robust = robustness(graph, int(config.get("robustness_repeats", 50)), float(config.get("edge_dropout", 0.1)), int(config.get("seed", 20260824)))
    candidates = candidate_table(graph, nodes, metrics, set(config.get("candidate_types", ["compound", "drug"])), set(config.get("disease_types", ["disease", "phenotype"])))
    out.mkdir(parents=True, exist_ok=True)
    nodes.to_csv(out / "nodes_validated.csv", index=False)
    edges.to_csv(out / "evidence_edges_validated.csv", index=False)
    metrics.to_csv(out / "node_metrics.csv", index=False)
    robust.to_csv(out / "robustness.csv", index=False)
    candidates.to_csv(out / "candidate_priority.csv", index=False)
    nx.write_graphml(graph, out / "evidence_graph.graphml")
    write_cards(out / "candidate_evidence_cards.md", candidates, nodes, edges)
    coverage = edges.groupby(["evidence_source", "source_version", "evidence_tier"], as_index=False).agg(n_evidence=("evidence_id", "nunique"), n_relations=("relation", "nunique"))
    coverage.to_csv(out / "evidence_coverage.csv", index=False)
    audit = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "nodes": str(nodes_path), "nodes_sha256": sha256(nodes_path),
        "edges": str(edges_path), "edges_sha256": sha256(edges_path),
        "n_nodes": graph.number_of_nodes(), "n_edges_collapsed": graph.number_of_edges(), "n_evidence_rows": len(edges),
        "seed": int(config.get("seed", 20260824)), "edge_dropout": float(config.get("edge_dropout", 0.1)),
        "claim_boundary": "centrality and network proximity are prioritization signals, not efficacy or causal mechanism",
    }
    (out / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    (out / "effective_config.yml").write_text(yaml.safe_dump(config, allow_unicode=True, sort_keys=False), encoding="utf-8")
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    try:
        print(run(args.config))
    except (ContractError, OSError, yaml.YAMLError, pd.errors.ParserError) as exc:
        print(f"[CONTRACT_ERROR] {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
