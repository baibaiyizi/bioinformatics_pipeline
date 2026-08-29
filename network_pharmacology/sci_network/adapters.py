from __future__ import annotations

import html
import re
from typing import Any

from .http import AuditClient, RetrievalError


OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"


class OpenTargets:
    def __init__(self, http: AuditClient):
        self.http = http

    def disease_targets(self, disease_id: str, page_size: int, max_targets: int) -> list[dict[str, Any]]:
        query = """query DiseaseTargets($id:String!,$page:Pagination!){disease(efoId:$id){id name associatedTargets(page:$page){count rows{target{id approvedSymbol} score datasourceScores{id score}}}}}"""
        rows: list[dict[str, Any]] = []
        expected: int | None = None
        page_index = 0
        while len(rows) < max_targets:
            variables = {"id": disease_id, "page": {"index": page_index, "size": page_size}}
            payload, _ = self.http.request_json("opentargets", f"disease_{disease_id}_{page_index}", "POST", OT_URL, json_body={"query": query, "variables": variables}, page=variables["page"])
            disease = payload.get("data", {}).get("disease")
            if not disease:
                raise RetrievalError(f"Open Targets 未识别 disease_id={disease_id}")
            block = disease["associatedTargets"]
            count = int(block["count"])
            if expected is None:
                expected = count
            elif expected != count:
                raise RetrievalError("Open Targets 分页期间 count 发生变化")
            page_rows = block.get("rows") or []
            rows.extend(page_rows)
            if not page_rows or len(rows) >= expected:
                break
            page_index += 1
        if expected is not None and expected <= max_targets and len(rows) != expected:
            raise RetrievalError(f"Open Targets 分页不完整: expected={expected}, retrieved={len(rows)}")
        return [{"disease_id": disease_id, "target_id": r["target"]["id"], "symbol": r["target"].get("approvedSymbol", ""), "score": float(r.get("score") or 0), "datasource_scores": r.get("datasourceScores") or []} for r in rows[:max_targets]]

    def drug_targets(self, chembl_id: str) -> list[dict[str, Any]]:
        # 2026-08 的 Drug schema 已移除 linkedTargets；直接使用无分页参数、
        # 完整返回的 mechanismsOfAction.rows，避免请求不存在的旧字段。
        query = """query DrugMechanisms($id:String!){drug(chemblId:$id){id name mechanismsOfAction{rows{mechanismOfAction actionType targets{id approvedSymbol}}}}}"""
        payload, _ = self.http.request_json("opentargets", f"drug_mechanisms_{chembl_id}", "POST", OT_URL, json_body={"query": query, "variables": {"id": chembl_id}}, page={"selection": "all_mechanismsOfAction_rows"})
        drug = payload.get("data", {}).get("drug")
        if not drug:
            raise RetrievalError(f"Open Targets 未识别 drug={chembl_id}")
        merged: dict[str, dict[str, Any]] = {}
        for mechanism in (drug.get("mechanismsOfAction") or {}).get("rows") or []:
            for target in mechanism.get("targets") or []:
                record = merged.setdefault(target["id"], {"compound_id": chembl_id, "target_id": target["id"], "symbol": target.get("approvedSymbol", ""), "mechanisms": [], "action_types": []})
                record["mechanisms"].append(mechanism.get("mechanismOfAction") or "")
                record["action_types"].append(mechanism.get("actionType") or "")
        return [{**row, "mechanisms": sorted(set(row["mechanisms"])), "action_types": sorted(set(row["action_types"]))} for row in merged.values()]


class ChEMBL:
    BASE = "https://www.ebi.ac.uk/chembl/api/data"

    def __init__(self, http: AuditClient):
        self.http = http

    def molecule(self, chembl_id: str) -> dict[str, Any]:
        payload, _ = self.http.request_json("chembl", f"molecule_{chembl_id}", "GET", f"{self.BASE}/molecule/{chembl_id}.json")
        if not payload or payload.get("molecule_chembl_id") != chembl_id:
            raise RetrievalError(f"ChEMBL molecule 身份不匹配: {chembl_id}")
        return payload

    def version(self) -> str:
        payload, _ = self.http.request_json("chembl", "status", "GET", f"{self.BASE}/status.json")
        version = str(payload.get("chembl_db_version") or "not_reported")
        self.http.set_source_version("chembl", version)
        return version

    def paged(self, endpoint: str, key: str, params: dict[str, Any], max_records: int = 1000) -> list[dict[str, Any]]:
        offset, limit, rows, expected = 0, 200, [], None
        while len(rows) < max_records:
            request_params = {**params, "format": "json", "limit": limit, "offset": offset}
            payload, _ = self.http.request_json("chembl", f"{endpoint}_{offset}", "GET", f"{self.BASE}/{endpoint}", params=request_params, page={"offset": offset, "limit": limit})
            meta = payload.get("page_meta") or {}
            count = int(meta.get("total_count", 0))
            if expected is None:
                expected = count
                if expected > max_records:
                    raise RetrievalError(f"ChEMBL {endpoint} total_count={expected} 超过显式检索上限 {max_records}；必须收窄查询或提高上限，禁止静默截断")
            elif expected != count:
                raise RetrievalError(f"ChEMBL {endpoint} 分页期间 total_count 变化")
            page_rows = payload.get(key) or []
            rows.extend(page_rows)
            if not page_rows or len(rows) >= count:
                break
            offset += limit
        if expected is not None and expected <= max_records and len(rows) != expected:
            raise RetrievalError(f"ChEMBL {endpoint} 分页不完整: expected={expected}, retrieved={len(rows)}")
        return rows

    def single_protein_target(self, accession: str, taxid: int) -> dict[str, Any] | None:
        rows = self.paged("target", "targets", {"target_components__accession": accession, "target_type": "SINGLE PROTEIN"}, max_records=100)
        matches = []
        for row in rows:
            components = row.get("target_components") or []
            accessions = {str(component.get("accession") or "") for component in components}
            if row.get("target_type") == "SINGLE PROTEIN" and int(row.get("tax_id") or 0) == taxid and accessions == {accession}:
                matches.append(row)
        if len(matches) > 1:
            raise RetrievalError(f"ChEMBL 单蛋白 target 映射一对多: {accession}, n={len(matches)}")
        return matches[0] if matches else None


class PubChem:
    BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

    def __init__(self, http: AuditClient):
        self.http = http

    def identity(self, id_type: str, value: str) -> dict[str, Any]:
        namespace = {"pubchem_cid": "cid", "inchikey": "inchikey"}.get(id_type)
        if not namespace:
            raise RetrievalError("PubChem 精确身份查询需要 CID 或 InChIKey")
        # PUG-REST 会在 property 响应中自动附带 CID；CID 本身不是可请求的
        # property 名，把它放入字段列表会得到 HTTP 400。
        fields = "MolecularFormula,MolecularWeight,ConnectivitySMILES,SMILES,InChIKey,Charge"
        payload, _ = self.http.request_json("pubchem", f"identity_{value}", "GET", f"{self.BASE}/compound/{namespace}/{value}/property/{fields}/JSON")
        rows = (payload.get("PropertyTable") or {}).get("Properties") or []
        if len(rows) != 1:
            raise RetrievalError(f"PubChem 身份解析非一对一: {id_type}={value}, n={len(rows)}")
        return rows[0]


class STRING:
    BASE = "https://string-db.org/api/json/network"

    def __init__(self, http: AuditClient):
        self.http = http

    def network(self, symbols: list[str], taxid: int, network_type: str, required_score: int) -> list[dict[str, Any]]:
        if network_type not in {"physical", "functional"}:
            raise ValueError("network_type must be physical or functional")
        if not symbols:
            return []
        params = {"identifiers": "\r".join(symbols), "species": taxid, "network_type": network_type, "required_score": required_score, "add_nodes": 0}
        payload, _ = self.http.request_json("string", f"{network_type}_{required_score}", "GET", self.BASE, params=params)
        if not isinstance(payload, list):
            raise RetrievalError("STRING network 响应不是列表")
        return payload

    def version(self) -> str:
        payload, _ = self.http.request_json("string", "version", "GET", "https://string-db.org/api/json/version")
        rows = payload if isinstance(payload, list) else []
        version = str(rows[0].get("string_version") if rows else "not_reported")
        self.http.set_source_version("string", version)
        return version


class Ensembl:
    BASE = "https://asia.ensembl.org/biomart/martservice"

    def __init__(self, http: AuditClient):
        self.http = http

    @staticmethod
    def _dataset(ensembl_id: str) -> str:
        return "mmusculus_gene_ensembl" if ensembl_id.startswith("ENSMUSG") else "hsapiens_gene_ensembl"

    def _query(self, request_id: str, dataset: str, filter_name: str, filter_value: str, attributes: list[str]) -> list[list[str]]:
        attribute_xml = "".join(f'<Attribute name="{html.escape(attribute, quote=True)}"/>' for attribute in attributes)
        query = f'<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6"><Dataset name="{html.escape(dataset, quote=True)}" interface="default"><Filter name="{html.escape(filter_name, quote=True)}" value="{html.escape(filter_value, quote=True)}"/>{attribute_xml}</Dataset></Query>'
        payload, _ = self.http.request_text("ensembl", request_id, "GET", self.BASE, params={"query": query}, page={"selection": "exact_identifier", "pagination": "not_applicable"})
        if payload.startswith("Query ERROR"):
            raise RetrievalError(f"Ensembl BioMart query 失败: {request_id}")
        rows = [line.rstrip("\r").split("\t") for line in payload.splitlines() if line.strip()]
        if any(len(row) != len(attributes) for row in rows):
            raise RetrievalError(f"Ensembl BioMart 列数不匹配: {request_id}")
        return rows

    def mouse_orthologs(self, human_id: str) -> list[dict[str, Any]]:
        rows = self._query(f"ortholog_{human_id}", "hsapiens_gene_ensembl", "ensembl_gene_id", human_id, ["mmusculus_homolog_ensembl_gene", "mmusculus_homolog_orthology_type"])
        return [{"type": row[1], "target": {"id": row[0]}} for row in rows if row[0]]

    def xrefs(self, ensembl_id: str) -> list[dict[str, Any]]:
        rows = self._query(f"xrefs_{ensembl_id}", self._dataset(ensembl_id), "ensembl_gene_id", ensembl_id, ["ensembl_gene_id", "uniprotswissprot"])
        return [{"dbname": "Uniprot/SWISSPROT", "primary_id": row[1]} for row in rows if row[0] == ensembl_id and row[1]]

    def lookup(self, ensembl_id: str) -> dict[str, Any]:
        rows = self._query(f"lookup_{ensembl_id}", self._dataset(ensembl_id), "ensembl_gene_id", ensembl_id, ["ensembl_gene_id", "external_gene_name"])
        if len(rows) != 1 or rows[0][0] != ensembl_id:
            raise RetrievalError(f"Ensembl lookup 身份不匹配: {ensembl_id}")
        return {"id": rows[0][0], "display_name": rows[0][1] or ensembl_id}

    def version(self) -> str:
        payload, _ = self.http.request_text("ensembl", "version", "GET", self.BASE, params={"type": "registry"}, page={"selection": "mart_registry", "pagination": "not_applicable"})
        match = re.search(r'database="ensembl_mart_(\d+)"', payload)
        if not match:
            raise RetrievalError("Ensembl BioMart registry 未报告可核验 release")
        version = f"release-{match.group(1)}"
        self.http.set_source_version("ensembl", version)
        return version


class Reactome:
    BASE = "https://reactome.org/ContentService/data/mapping/UniProt"

    def __init__(self, http: AuditClient):
        self.http = http

    def pathways(self, uniprot: str) -> list[dict[str, Any]]:
        payload, _ = self.http.request_json("reactome", f"pathways_{uniprot}", "GET", f"{self.BASE}/{uniprot}/pathways")
        if not isinstance(payload, list):
            raise RetrievalError("Reactome mapping 响应不是列表")
        return payload

    def version(self) -> str:
        payload, _ = self.http.request_json("reactome", "version", "GET", "https://reactome.org/ContentService/data/database/version")
        version = str(payload)
        self.http.set_source_version("reactome", version)
        return version


class UniProt:
    BASE = "https://rest.uniprot.org/uniprotkb"

    def __init__(self, http: AuditClient):
        self.http = http

    def entry(self, accession: str) -> dict[str, Any]:
        payload, _ = self.http.request_json("uniprot", f"entry_{accession}", "GET", f"{self.BASE}/{accession}", params={"format": "json"})
        return payload


class RCSB:
    URL = "https://search.rcsb.org/rcsbsearch/v2/query"

    def __init__(self, http: AuditClient):
        self.http = http

    def by_uniprot(self, accession: str) -> list[str]:
        # 这里回答的是“是否至少有一个实验结构”，不是构建完整 PDB 数据集；
        # 因而显式做 top-1 查询，并保留 total_count 供来源覆盖审计。
        body = {"query": {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession", "operator": "exact_match", "value": accession}}, "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": 1}, "sort": [{"sort_by": "score", "direction": "desc"}]}}
        payload, _ = self.http.request_json("rcsb", f"uniprot_{accession}", "POST", self.URL, json_body=body, page={"start": 0, "rows": 1, "selection": "top_1_structure_availability"})
        rows = payload.get("result_set") or []
        return [str(row["identifier"]) for row in rows]


class AlphaFold:
    BASE = "https://alphafold.ebi.ac.uk/api/prediction"

    def __init__(self, http: AuditClient):
        self.http = http

    def prediction(self, accession: str) -> list[dict[str, Any]]:
        payload, _ = self.http.request_json("alphafold", f"prediction_{accession}", "GET", f"{self.BASE}/{accession}")
        return payload if isinstance(payload, list) else []
