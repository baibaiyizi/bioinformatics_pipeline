from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import networkx as nx
import pandas as pd
import yaml

from .adapters import AlphaFold, ChEMBL, Ensembl, OpenTargets, PubChem, RCSB, Reactome, STRING, UniProt
from .analysis import build_multidigraph, rank_candidates, reactome_ora, string_threshold_sensitivity
from .contract import ContractError, load_config, resolve
from .http import AuditClient, Transport


class EvidenceStore:
    def __init__(self):
        self.nodes: dict[str, dict[str, Any]] = {}
        self.edges: list[dict[str, Any]] = []

    def node(self, node_id: str, node_type: str, canonical_id: str, species_taxid: int | str = "", label: str = "") -> str:
        record = {"node_id": node_id, "node_type": node_type, "canonical_id": canonical_id, "species_taxid": species_taxid, "label": label}
        previous = self.nodes.get(node_id)
        if previous is not None and previous != record:
            raise ContractError(f"实体 {node_id} 的规范化身份冲突")
        self.nodes[node_id] = record
        return node_id

    def edge(self, source: str, target: str, relation: str, database: str, database_version: str, evidence_id: str, evidence_layer: str, species_taxid: int | str, score: float | None, *, lane: str = "", network_type: str = "", source_symbol: str = "", target_symbol: str = "") -> None:
        self.edges.append({"source": source, "target": target, "relation": relation, "database": database, "database_version": database_version, "evidence_id": evidence_id, "evidence_layer": evidence_layer, "species_taxid": species_taxid, "score": score if score is not None else math.nan, "lane": lane, "network_type": network_type, "source_symbol": source_symbol, "target_symbol": target_symbol})


def _one_mouse_ortholog(ensembl: Ensembl, human_id: str) -> tuple[str, str]:
    homologs = [row for row in ensembl.mouse_orthologs(human_id) if str(row.get("type", "")) == "ortholog_one2one"]
    if len(homologs) != 1:
        raise ContractError(f"human→mouse ortholog 映射不是一对一: {human_id}, n={len(homologs)}")
    target_id = str(homologs[0]["target"]["id"])
    lookup = ensembl.lookup(target_id)
    return target_id, str(lookup.get("display_name") or target_id)


def _uniprot_from_xrefs(rows: list[dict[str, Any]], target_id: str) -> str | None:
    reviewed = sorted({str(row.get("primary_id")) for row in rows if row.get("primary_id") and "SWISSPROT" in str(row.get("dbname", "")).upper()})
    if len(reviewed) > 1:
        raise ContractError(f"Ensembl→reviewed UniProt 映射一对多: {target_id}: {reviewed}")
    return reviewed[0] if reviewed else None


def _source_version(http: AuditClient, source: str, api_surface: str) -> str:
    versions = {str(row["database_version"]) for row in http.records if row["source"] == source and row["database_version"] != "not_reported"}
    return "|".join(sorted(versions)) if versions else f"{api_surface};release=not_reported"


def run(config_path: Path, *, transport: Transport | None = None) -> Path:
    config_path = config_path.resolve()
    config = load_config(config_path)
    out = resolve(config_path.parent, str(config["output_dir"]))
    out.mkdir(parents=True, exist_ok=True)
    http = AuditClient(out, float(config["retrieval"]["timeout_seconds"]), transport)
    ot, chembl, pubchem = OpenTargets(http), ChEMBL(http), PubChem(http)
    string, ensembl = STRING(http), Ensembl(http)
    reactome, uniprot, rcsb, alphafold = Reactome(http), UniProt(http), RCSB(http), AlphaFold(http)
    taxid = int(config["species_taxid"])
    max_targets = int(config["retrieval"]["max_targets"])
    http.set_source_version("opentargets", "api-v4;platform-release=not_reported")
    http.set_source_version("pubchem", "PUG-REST;data-release=not_reported")
    http.set_source_version("rcsb", "search-api-v2;archive-release=not_reported")
    http.set_source_version("alphafold", "prediction-api;record-versioned")
    ensembl.version()
    string_version = string.version()
    chembl_version = chembl.version() if any(compound["id_type"] == "chembl" for compound in config["compounds"]) else "not_used"
    reactome_version = reactome.version() if config["analysis"]["reactome_background_uniprot"] else "not_used"
    store = EvidenceStore()
    project_node = store.node(f"project:{config['project_id']}", "project", config["project_id"], taxid, config["project_id"])

    target_records: dict[str, dict[str, Any]] = {}
    disease_scores: dict[str, float] = {}
    direct_scores: dict[str, float] = {}
    direct_layers: dict[str, set[str]] = {}
    compound_records: list[dict[str, Any]] = []
    chembl_activities: list[dict[str, Any]] = []
    chembl_mechanisms: list[dict[str, Any]] = []
    id_mappings: list[dict[str, Any]] = []
    chembl_parent_by_input: dict[str, str] = {}
    compound_node_by_input: dict[str, str] = {}

    def register_target(target_id: str, symbol: str, layer: str) -> None:
        current = target_records.get(target_id)
        if current and current["symbol"] not in {symbol, target_id} and symbol not in {"", target_id}:
            raise ContractError(f"target symbol 映射冲突: {target_id}")
        target_records[target_id] = {"target_id": target_id, "symbol": symbol or target_id, "species_taxid": taxid, "evidence_layer": layer if not current else (current["evidence_layer"] if current["evidence_layer"] == layer else "mixed_explicit_layers")}
        store.node(f"target:{target_id}", "target", target_id, taxid, symbol or target_id)

    for seed in config["target_seeds"]:
        target_id = seed["ensembl_gene_id"]
        lookup = ensembl.lookup(target_id)
        register_target(target_id, str(lookup.get("display_name") or target_id), "direct_species")
        store.edge(project_node, f"target:{target_id}", "target_seed", seed["evidence_source"], "user_provided", seed["evidence_id"], "direct_species", taxid, 1.0)

    for disease_id in config["disease_ids"]:
        disease_node = store.node(f"disease:{disease_id}", "disease", disease_id, "", disease_id)
        rows = ot.disease_targets(disease_id, int(config["retrieval"]["page_size"]), max_targets)
        for row in rows:
            human_id, symbol, layer = row["target_id"], row["symbol"], "direct_species"
            target_id = human_id
            if taxid == 10090:
                target_id, symbol = _one_mouse_ortholog(ensembl, human_id)
                layer = "human_ortholog"
                id_mappings.append({"source_id": human_id, "target_id": target_id, "mapping": "Ensembl ortholog", "layer": layer})
            register_target(target_id, symbol, layer)
            disease_scores[target_id] = max(disease_scores.get(target_id, 0), float(row["score"]))
            store.edge(disease_node, f"target:{target_id}", "associated_with", "Open Targets", "api-v4", f"{disease_id}:{human_id}", layer, taxid, float(row["score"]), lane="disease_association")

    for compound in config["compounds"]:
        compound_id, id_type = compound["id"], compound["id_type"]
        identity: dict[str, Any]
        if id_type == "chembl":
            molecule = chembl.molecule(compound_id)
            hierarchy = molecule.get("molecule_hierarchy") or {}
            parent_id = str(hierarchy.get("parent_chembl_id") or compound_id)
            chembl_parent_by_input[compound_id] = parent_id
            structures = molecule.get("molecule_structures") or {}
            inchikey = str(structures.get("standard_inchi_key") or "")
            if not inchikey:
                raise ContractError(f"ChEMBL {compound_id} 缺少标准 InChIKey，不能核验 PubChem 精确身份")
            identity = pubchem.identity("inchikey", inchikey)
            chembl_mechanisms.extend(chembl.paged("mechanism", "mechanisms", {"parent_molecule_chembl_id": parent_id}))
            drug_targets = ot.drug_targets(compound_id)
        else:
            identity = pubchem.identity(id_type, compound_id)
            drug_targets = []
        pubchem_id = str(identity["CID"])
        node_id = store.node(f"compound:CID{pubchem_id}", "compound", f"CID:{pubchem_id}", "", compound_id)
        compound_node_by_input[compound_id] = node_id
        compound_records.append({"input_id": compound_id, "input_id_type": id_type, "chembl_parent_id": chembl_parent_by_input.get(compound_id, ""), "node_id": node_id, "pubchem_cid": pubchem_id, "inchikey": identity.get("InChIKey", ""), "canonical_smiles": identity.get("ConnectivitySMILES") or identity.get("CanonicalSMILES") or "", "isomeric_smiles": identity.get("SMILES") or identity.get("IsomericSMILES") or "", "formal_charge": identity.get("Charge", "")})
        for row in drug_targets:
            human_id, symbol, layer = row["target_id"], row["symbol"], "direct_species"
            target_id = human_id
            if taxid == 10090:
                target_id, symbol = _one_mouse_ortholog(ensembl, human_id)
                layer = "human_ortholog"
                id_mappings.append({"source_id": human_id, "target_id": target_id, "mapping": "Ensembl ortholog", "layer": layer})
            register_target(target_id, symbol, layer)
            direct_scores[target_id] = 1.0
            direct_layers.setdefault(target_id, set()).add(layer)
            evidence_id = f"{compound_id}:{human_id}"
            store.edge(node_id, f"target:{target_id}", "targets", "Open Targets", "api-v4;mechanismsOfAction", evidence_id, layer, taxid, 1.0, lane="direct_compound_target")

    if len(target_records) > max_targets:
        ordered = sorted(target_records, key=lambda target: (target not in direct_scores, -disease_scores.get(target, 0), target))[:max_targets]
        keep = set(ordered)
        target_records = {key: value for key, value in target_records.items() if key in keep}
        disease_scores = {key: value for key, value in disease_scores.items() if key in keep}
        direct_scores = {key: value for key, value in direct_scores.items() if key in keep}

    targets = pd.DataFrame(target_records.values())
    if targets.empty:
        raise ContractError("公开检索与 target_seeds 均未形成可分析 target")
    symbols = sorted({symbol for symbol in targets["symbol"] if symbol})
    required_score = min(int(config["retrieval"]["string_required_score"]), min(config["analysis"]["string_thresholds"]))
    string_rows: list[dict[str, Any]] = []
    symbol_to_target = targets.set_index("symbol")["target_id"].to_dict()
    for network_type in ("physical", "functional"):
        for row in string.network(symbols, taxid, network_type, required_score):
            source_symbol, target_symbol = str(row.get("preferredName_A") or ""), str(row.get("preferredName_B") or "")
            if source_symbol not in symbol_to_target or target_symbol not in symbol_to_target:
                continue
            source_id, target_id = symbol_to_target[source_symbol], symbol_to_target[target_symbol]
            score = float(row.get("score") or 0)
            edge = {"source": f"target:{source_id}", "target": f"target:{target_id}", "relation": "physical_interaction" if network_type == "physical" else "functional_association", "database": "STRING", "database_version": string_version, "evidence_id": f"{network_type}:{row.get('stringId_A')}:{row.get('stringId_B')}", "evidence_layer": "direct_species", "species_taxid": taxid, "score": score, "lane": "ppi_proximity" if network_type == "physical" else "", "network_type": network_type, "source_symbol": source_symbol, "target_symbol": target_symbol}
            store.edges.append(edge)
            string_rows.append(edge)

    ppi_scores = {target_id: 0.0 for target_id in target_records}
    for row in string_rows:
        if row["network_type"] != "physical":
            continue
        for symbol in (row["source_symbol"], row["target_symbol"]):
            target_id = symbol_to_target[symbol]
            ppi_scores[target_id] = max(ppi_scores[target_id], float(row["score"]))

    uniprot_by_target: dict[str, str] = {}
    for row in targets.sort_values("target_id").to_dict("records"):
        target_id = row["target_id"]
        accession = _uniprot_from_xrefs(ensembl.xrefs(target_id), target_id)
        id_mappings.append({"source_id": target_id, "target_id": accession or "none", "mapping": "Ensembl xref reviewed UniProt", "layer": row["evidence_layer"]})
        if not accession:
            continue
        protein_record = uniprot.entry(accession)
        if str(protein_record.get("primaryAccession") or "") != accession:
            raise ContractError(f"UniProt accession 身份不匹配: {accession}")
        uniprot_taxid = int((protein_record.get("organism") or {}).get("taxonId") or 0)
        if uniprot_taxid != taxid:
            raise ContractError(f"UniProt 物种不匹配: {accession}, expected={taxid}, observed={uniprot_taxid}")
        uniprot_by_target[target_id] = accession

    for target_id, accession in uniprot_by_target.items():
        chembl_target = chembl.single_protein_target(accession, taxid) if chembl_parent_by_input else None
        if not chembl_target:
            continue
        target_chembl_id = str(chembl_target["target_chembl_id"])
        id_mappings.append({"source_id": accession, "target_id": target_chembl_id, "mapping": "ChEMBL single protein target", "layer": target_records[target_id]["evidence_layer"]})
        for compound_id, parent_id in chembl_parent_by_input.items():
            node_id = compound_node_by_input[compound_id]
            mechanisms = [row for row in chembl_mechanisms if str(row.get("parent_molecule_chembl_id") or "") == parent_id and str(row.get("target_chembl_id") or "") == target_chembl_id]
            for index, mechanism in enumerate(mechanisms):
                evidence_id = f"{parent_id}:{target_chembl_id}:mechanism:{index}:{mechanism.get('action_type') or 'unknown'}"
                store.edge(node_id, f"target:{target_id}", "curated_mechanism_target", "ChEMBL", chembl_version, evidence_id, target_records[target_id]["evidence_layer"], taxid, 1.0, lane="direct_compound_target")
                direct_scores[target_id] = 1.0
                direct_layers.setdefault(target_id, set()).add(target_records[target_id]["evidence_layer"])
            activities = chembl.paged("activity", "activities", {"parent_molecule_chembl_id": parent_id, "target_chembl_id": target_chembl_id, "pchembl_value__isnull": "false"}, max_records=10000)
            for activity in activities:
                relation = str(activity.get("standard_relation") or "")
                valid = bool(activity.get("pchembl_value") not in {None, ""} and relation == "=" and not activity.get("data_validity_comment") and int(activity.get("potential_duplicate") or 0) == 0 and int(activity.get("target_tax_id") or 0) == taxid)
                normalized = {**activity, "sci_target_ensembl_id": target_id, "sci_target_uniprot": accession, "sci_target_chembl_id": target_chembl_id, "sci_usable_direct_evidence": valid, "sci_exclusion_reason": "none" if valid else "non_exact_or_flagged_activity"}
                chembl_activities.append(normalized)
                if not valid:
                    continue
                activity_id = str(activity.get("activity_id") or f"{activity.get('assay_chembl_id', 'assay')}:{activity.get('pchembl_value')}")
                store.edge(node_id, f"target:{target_id}", "measured_activity", "ChEMBL", chembl_version, f"activity:{activity_id}", target_records[target_id]["evidence_layer"], taxid, float(activity["pchembl_value"]), lane="direct_compound_target")
                direct_scores[target_id] = 1.0
                direct_layers.setdefault(target_id, set()).add(target_records[target_id]["evidence_layer"])

    structure_rows: list[dict[str, Any]] = []
    structure_scores: dict[str, float] = {}
    structure_target_ids = sorted(target_records, key=lambda target: (target not in direct_scores, -disease_scores.get(target, 0), target))[: int(config["retrieval"]["max_structure_checks"])]
    for target_id in structure_target_ids:
        row = target_records[target_id]
        accession = uniprot_by_target.get(target_id)
        if not accession:
            continue
        pdb_ids = rcsb.by_uniprot(accession)
        if pdb_ids:
            structure_scores[target_id] = 1.0
            structure_rows.append({"target_id": target_id, "uniprot": accession, "structure_source": "RCSB", "structure_id": pdb_ids[0], "structure_query_scope": "top_1_availability", "modeling_role": "candidate_requires_holo_validation"})
            structure_node = store.node(f"structure:PDB:{pdb_ids[0]}", "structure", pdb_ids[0], taxid, pdb_ids[0])
            store.edge(f"target:{target_id}", structure_node, "has_experimental_structure", "RCSB PDB", "search-api-v2", f"{accession}:{pdb_ids[0]}", row["evidence_layer"], taxid, 1.0, lane="structure_availability")
        else:
            predictions = alphafold.prediction(accession)
            if predictions:
                structure_scores[target_id] = 0.5
                model_id = str(predictions[0].get("modelEntityId") or f"AF-{accession}-F1")
                structure_rows.append({"target_id": target_id, "uniprot": accession, "structure_source": "AlphaFold", "structure_id": model_id, "structure_query_scope": "used_only_after_no_RCSB_top_1_hit", "modeling_role": "background_only_no_redocking"})
                structure_node = store.node(f"structure:AF:{model_id}", "predicted_structure", model_id, taxid, model_id)
                store.edge(f"target:{target_id}", structure_node, "has_predicted_structure", "AlphaFold DB", f"api-v{predictions[0].get('latestVersion', 'unknown')}", f"{accession}:{model_id}", row["evidence_layer"], taxid, 0.5, lane="structure_availability")

    pathway_results = pd.DataFrame(columns=["pathway_id", "pathway_name", "species_taxid", "selected_hits", "selected_total", "background_hits", "background_total", "p_value", "fdr_bh"])
    pathway_scores: dict[str, float] = {}
    background = list(config["analysis"]["reactome_background_uniprot"])
    selected_accessions = set(uniprot_by_target.values())
    if background:
        universe = list(dict.fromkeys(background))
        if not selected_accessions <= set(universe):
            raise ContractError("Reactome 显式背景必须包含所有已映射候选 UniProt")
        background_map = {accession: reactome.pathways(accession) for accession in universe}
        pathway_results = reactome_ora(selected_accessions, background_map, taxid)
        significant = set(pathway_results.loc[pathway_results["fdr_bh"] <= 0.05, "pathway_id"]) if not pathway_results.empty else set()
        for target_id, accession in uniprot_by_target.items():
            memberships = {str(row.get("stId") or "") for row in background_map[accession]}
            pathway_scores[target_id] = float(len(memberships & significant))
            for pathway_id in memberships & significant:
                pathway_node = store.node(f"pathway:{pathway_id}", "pathway", pathway_id, taxid, pathway_id)
                store.edge(f"target:{target_id}", pathway_node, "member_of", "Reactome", reactome_version, f"{accession}:{pathway_id}", target_records[target_id]["evidence_layer"], taxid, 1.0, lane="pathway_support")

    lane_rows: list[dict[str, Any]] = []
    def add_lane(values: dict[str, float], lane: str, available: bool = True) -> None:
        if available:
            for target_id, value in values.items():
                score = float(value)
                if target_id in target_records and math.isfinite(score) and score > 0:
                    lane_rows.append({"target_id": target_id, "lane": lane, "raw_score": score})
    add_lane(disease_scores, "disease_association", bool(config["disease_ids"]))
    add_lane(ppi_scores, "ppi_proximity", any(row["network_type"] == "physical" for row in string_rows))
    add_lane(pathway_scores, "pathway_support", bool(background))
    add_lane(structure_scores, "structure_availability", bool(structure_rows))
    add_lane(direct_scores, "direct_compound_target", bool(direct_scores))
    lane_scores = pd.DataFrame(lane_rows, columns=["target_id", "lane", "raw_score"])
    ranked, leave_one = rank_candidates(targets, lane_scores)

    nodes = pd.DataFrame(store.nodes.values())
    edges = pd.DataFrame(store.edges)
    graph = build_multidigraph(nodes, edges)
    out.mkdir(parents=True, exist_ok=True)
    nodes.to_csv(out / "standard_entities.tsv", sep="\t", index=False)
    edges.to_csv(out / "evidence_edges.tsv", sep="\t", index=False)
    (pd.DataFrame(id_mappings) if id_mappings else pd.DataFrame(columns=["source_id", "target_id", "mapping", "layer"])).to_csv(out / "id_mapping.tsv", sep="\t", index=False)
    (pd.DataFrame(compound_records) if compound_records else pd.DataFrame(columns=["input_id", "input_id_type", "chembl_parent_id", "node_id", "pubchem_cid", "inchikey", "canonical_smiles", "isomeric_smiles", "formal_charge"])).to_csv(out / "compound_identity.tsv", sep="\t", index=False)
    (pd.DataFrame(chembl_activities) if chembl_activities else pd.DataFrame(columns=["activity_id", "assay_chembl_id", "pchembl_value", "standard_relation", "sci_target_ensembl_id", "sci_target_uniprot", "sci_target_chembl_id", "sci_usable_direct_evidence", "sci_exclusion_reason"])).to_csv(out / "chembl_activity.tsv", sep="\t", index=False)
    (pd.DataFrame(chembl_mechanisms) if chembl_mechanisms else pd.DataFrame(columns=["parent_molecule_chembl_id", "target_chembl_id", "action_type", "mechanism_of_action"])).to_csv(out / "chembl_mechanism.tsv", sep="\t", index=False)
    pathway_results.to_csv(out / "reactome_ora.tsv", sep="\t", index=False)
    ranked.to_csv(out / "candidate_priority.tsv", sep="\t", index=False)
    lane_scores.to_csv(out / "evidence_lane_scores.tsv", sep="\t", index=False)
    leave_one.to_csv(out / "leave_one_lane_out.tsv", sep="\t", index=False)
    sensitivity = string_threshold_sensitivity(edges, targets, config["analysis"]["string_thresholds"])
    sensitivity.to_csv(out / "string_threshold_sensitivity.tsv", sep="\t", index=False)
    nx.write_graphml(graph, out / "evidence_multidigraph.graphml")
    coverage = edges.groupby(["database", "database_version", "evidence_layer", "lane"], dropna=False, as_index=False).agg(n_evidence=("evidence_id", "nunique"), n_edges=("evidence_id", "size"))
    coverage.to_csv(out / "source_coverage.tsv", sep="\t", index=False)

    structure_lookup = {row["target_id"]: row for row in structure_rows}
    handoff_rows = []
    for compound in compound_records:
        for target in ranked[ranked["eligible"]].to_dict("records"):
            structure = structure_lookup.get(target["target_id"], {})
            handoff_rows.append({"compound_input_id": compound["input_id"], "pubchem_cid": compound["pubchem_cid"], "target_ensembl_id": target["target_id"], "target_symbol": target["symbol"], "species_taxid": taxid, "evidence_layer": target["evidence_layer"], "priority_score": target["priority_score"], "evidence_lane_count": target["evidence_lane_count"], "direct_evidence_layer": "|".join(sorted(direct_layers.get(target["target_id"], set()))) or "none", "structure_source": structure.get("structure_source", "none"), "structure_id": structure.get("structure_id", "none"), "handoff_status": "candidate_only_requires_modeling_preflight"})
    (pd.DataFrame(handoff_rows) if handoff_rows else pd.DataFrame(columns=["compound_input_id", "pubchem_cid", "target_ensembl_id", "target_symbol", "species_taxid", "evidence_layer", "priority_score", "evidence_lane_count", "direct_evidence_layer", "structure_source", "structure_id", "handoff_status"])).to_csv(out / "modeling_handoff.tsv", sep="\t", index=False)

    http.write_manifest()
    audit = {"project_id": config["project_id"], "created_at": datetime.now(timezone.utc).isoformat(), "species_taxid": taxid, "n_targets": len(targets), "n_entities": graph.number_of_nodes(), "n_evidence_edges": graph.number_of_edges(), "eligible_targets": int(ranked["eligible"].sum()), "database_versions": {source: _source_version(http, source, surface) for source, surface in {"opentargets": "api-v4", "chembl": "api-data", "pubchem": "PUG-REST", "string": "REST", "ensembl": "REST", "uniprot": "REST", "reactome": "ContentService", "rcsb": "search-v2", "alphafold": "prediction-api"}.items() if any(row["source"] == source for row in http.records)}, "pagination_calls": len(http.records), "claim_boundary": "Network rank is candidate priority only, not efficacy, causality, binding affinity, stability, or mechanism.", "species_boundary": "Human disease evidence in mouse projects is stored only in the human_ortholog layer."}
    (out / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    (out / "effective_config.yml").write_text(yaml.safe_dump(config, allow_unicode=True, sort_keys=False), encoding="utf-8")
    return out
