#!/usr/bin/env python3
"""复制一个分析模板并记录实际源码版本；拒绝覆盖现有项目。"""
import argparse
import csv
import hashlib
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parent


def copy_template(module: str, destination: Path):
    source = ROOT / module
    if Path(module).name != module or not (source / "README.md").is_file():
        raise ValueError("请提供模板库中的单个模块名称")
    if destination.exists():
        raise FileExistsError(f"目标已存在，拒绝覆盖: {destination}")
    destination = destination.resolve()
    if destination.is_relative_to(ROOT):
        raise ValueError("项目副本应位于模板库之外")
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    dirty = bool(subprocess.check_output(["git", "status", "--porcelain", "--", module], cwd=ROOT))
    shutil.copytree(source, destination, symlinks=True,
                    ignore=shutil.ignore_patterns(".git", ".venv*", "__pycache__", ".cache", ".pytest_cache", "result", "runs", "history"))
    with (destination / "template_snapshot.tsv").open("w") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["module", "source_commit", "source_dirty", "relative_path", "sha256"])
        for path in sorted(destination.rglob("*")):
            if path.is_file() and not path.is_symlink() and path.name != "template_snapshot.tsv":
                writer.writerow([module, commit, str(dirty).lower(), str(path.relative_to(destination)),
                                 hashlib.sha256(path.read_bytes()).hexdigest()])
    print(f"已复制 {module} 到 {destination}；修改配置及样本表后先运行模块预检。")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("module")
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    copy_template(args.module, args.destination)


if __name__ == "__main__":
    main()
