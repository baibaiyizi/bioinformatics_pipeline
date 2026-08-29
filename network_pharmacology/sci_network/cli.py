from __future__ import annotations

import argparse
import sys
from pathlib import Path

import httpx
import yaml

from .contract import ContractError
from .http import RetrievalError
from .pipeline import run


def main() -> int:
    parser = argparse.ArgumentParser(description="Seed-driven public evidence network pharmacology")
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    try:
        print(run(args.config))
    except (ContractError, RetrievalError, OSError, yaml.YAMLError, httpx.HTTPError) as exc:
        print(f"[PIPELINE_ERROR] {exc}", file=sys.stderr)
        return 2
    return 0
