from __future__ import annotations

import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import httpx


class RetrievalError(RuntimeError):
    """A public-source request failed or was incomplete."""


Transport = Callable[..., tuple[int, dict[str, str], Any]]
RETRYABLE_STATUS = {429, 500, 502, 503, 504}


def _canonical_bytes(payload: Any) -> bytes:
    return json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


class AuditClient:
    def __init__(self, output_dir: Path, timeout: float = 30, transport: Transport | None = None):
        self.output_dir = output_dir
        self.timeout = timeout
        self.transport = transport
        self.records: list[dict[str, Any]] = []
        self.source_versions: dict[str, str] = {}
        self.raw_dir = output_dir / "provenance" / "raw"

    def set_source_version(self, source: str, version: str) -> None:
        value = str(version).strip()
        if value:
            self.source_versions[source] = value

    def request_json(self, source: str, request_id: str, method: str, url: str, *, params: dict[str, Any] | None = None, json_body: dict[str, Any] | None = None, page: dict[str, Any] | None = None) -> tuple[Any, dict[str, str]]:
        method = method.upper()
        status = 0
        headers: dict[str, str] = {}
        payload: Any = None
        max_attempts = 4
        for attempt in range(max_attempts):
            try:
                if self.transport is not None:
                    status, headers, payload = self.transport(method=method, url=url, params=params, json=json_body)
                    headers = {str(k).lower(): str(v) for k, v in headers.items()}
                else:
                    with httpx.Client(timeout=self.timeout, follow_redirects=True, headers={"Accept": "application/json", "User-Agent": "sci-network-pharmacology/0.2"}) as client:
                        response = client.request(method, url, params=params, json=json_body)
                    status = response.status_code
                    headers = {k.lower(): v for k, v in response.headers.items()}
                    if status == 204:
                        payload = {}
                    else:
                        try:
                            payload = response.json()
                        except ValueError as exc:
                            if status in RETRYABLE_STATUS and attempt < max_attempts - 1:
                                time.sleep(min(float(headers.get("retry-after", "1") or 1), 5.0))
                                continue
                            raise RetrievalError(f"{source} 返回非 JSON 响应: HTTP {status}") from exc
            except httpx.TransportError as exc:
                if attempt < max_attempts - 1:
                    time.sleep(1.0)
                    continue
                raise RetrievalError(f"{source} 网络请求失败，重试 {max_attempts} 次后仍不可用: {type(exc).__name__}") from exc
            if status not in RETRYABLE_STATUS or attempt == max_attempts - 1:
                break
            time.sleep(min(float(headers.get("retry-after", "1") or 1), 5.0))
        if status < 200 or status >= 300:
            raise RetrievalError(f"{source} 请求失败: HTTP {status}, request_id={request_id}")
        if isinstance(payload, dict) and payload.get("errors"):
            raise RetrievalError(f"{source} GraphQL 返回 errors, request_id={request_id}")

        digest = hashlib.sha256(_canonical_bytes(payload)).hexdigest()
        record = {"source": source, "request_id": request_id, "accessed_at": datetime.now(timezone.utc).isoformat(), "method": method, "url": url, "params": params or {}, "json": json_body or {}, "status": status, "response_sha256": digest, "database_version": headers.get("x-uniprot-release") or headers.get("x-api-version") or headers.get("x-release") or self.source_versions.get(source, "not_reported"), "page": page or {}}
        self.records.append(record)
        self.raw_dir.mkdir(parents=True, exist_ok=True)
        (self.raw_dir / f"{len(self.records):04d}_{source}_{request_id}.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
        return payload, headers

    def request_text(self, source: str, request_id: str, method: str, url: str, *, params: dict[str, Any] | None = None, page: dict[str, Any] | None = None) -> tuple[str, dict[str, str]]:
        method = method.upper()
        status = 0
        headers: dict[str, str] = {}
        payload = ""
        max_attempts = 4
        for attempt in range(max_attempts):
            try:
                if self.transport is not None:
                    status, headers, raw_payload = self.transport(method=method, url=url, params=params, json=None)
                    headers = {str(k).lower(): str(v) for k, v in headers.items()}
                    payload = str(raw_payload)
                else:
                    with httpx.Client(timeout=self.timeout, follow_redirects=True, headers={"Accept": "text/plain, application/xml", "User-Agent": "sci-network-pharmacology/0.2"}) as client:
                        response = client.request(method, url, params=params)
                    status = response.status_code
                    headers = {k.lower(): v for k, v in response.headers.items()}
                    payload = response.text
            except httpx.TransportError as exc:
                if attempt < max_attempts - 1:
                    time.sleep(1.0)
                    continue
                raise RetrievalError(f"{source} 网络请求失败，重试 {max_attempts} 次后仍不可用: {type(exc).__name__}") from exc
            if status not in RETRYABLE_STATUS or attempt == max_attempts - 1:
                break
            time.sleep(min(float(headers.get("retry-after", "1") or 1), 5.0))
        if status < 200 or status >= 300:
            raise RetrievalError(f"{source} 请求失败: HTTP {status}, request_id={request_id}")
        digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
        record = {"source": source, "request_id": request_id, "accessed_at": datetime.now(timezone.utc).isoformat(), "method": method, "url": url, "params": params or {}, "json": {}, "status": status, "response_sha256": digest, "database_version": headers.get("x-api-version") or headers.get("x-release") or self.source_versions.get(source, "not_reported"), "page": page or {}}
        self.records.append(record)
        self.raw_dir.mkdir(parents=True, exist_ok=True)
        (self.raw_dir / f"{len(self.records):04d}_{source}_{request_id}.txt").write_text(payload, encoding="utf-8")
        return payload, headers

    def write_manifest(self) -> Path:
        path = self.output_dir / "provenance" / "requests.jsonl"
        path.parent.mkdir(parents=True, exist_ok=True)
        for record in self.records:
            if record["database_version"] == "not_reported" and record["source"] in self.source_versions:
                record["database_version"] = self.source_versions[record["source"]]
        path.write_text("".join(json.dumps(row, ensure_ascii=False) + "\n" for row in self.records), encoding="utf-8")
        return path
