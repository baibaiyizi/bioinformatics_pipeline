from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import yaml

from .contract import ContractError
from .runtime import CommandError
from .workflow import analyze, dock, md_prepare, md_production, md_smoke, preflight


STAGES = {"preflight": preflight, "dock": dock, "md-prepare": md_prepare, "md-smoke": md_smoke, "md-production": md_production, "analyze": analyze}


def main() -> int:
    parser = argparse.ArgumentParser(description="Vina to GROMACS staged workflow")
    parser.add_argument("stage", choices=STAGES)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    try:
        print(STAGES[args.stage](args.config))
    except (ContractError, CommandError, OSError, yaml.YAMLError, pd.errors.ParserError, ImportError) as exc:
        print(f"[STAGE_ERROR] {exc}", file=sys.stderr)
        return 2
    return 0
