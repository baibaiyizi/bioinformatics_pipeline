from pathlib import Path

from ai_research_harness.cli import build_parser


def test_cli_surface() -> None:
    parser = build_parser()
    assert parser.parse_args(["doctor"]).command == "doctor"
    args = parser.parse_args(["run", "--project", "/tmp/project", "--prompt", "audit"])
    assert args.prompt == "audit"
    args = parser.parse_args(["evolve", "--project", "/tmp/project", "--run-id", "r"])
    assert args.scope == "project" and args.concurrency == 3
    args = parser.parse_args(["gate", "--candidate", "recuris-v1"])
    assert args.scope == "global"
