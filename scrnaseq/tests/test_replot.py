import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
from PIL import Image
import pytest

from scrna_template.replot import replot
from scrna_template.runtime import sha256


def _fixture(tmp_path):
    root = tmp_path / "project with spaces"
    source = root / "runs/frozen"
    for name in ("00_audit", "07_source_data", "08_figures"):
        (source / name).mkdir(parents=True)
    (source / "00_audit/run_result.json").write_text(json.dumps({"status": "PASS"}))
    frame = pd.DataFrame({"UMAP1": [0., 1., 2., 3.], "UMAP2": [1., 0., 3., 2.], "label": ["WT", "WT", "PFOS", "PFOS"]})
    path = source / "07_source_data/umap_condition.tsv"
    frame.to_csv(path, sep="\t", index=False)
    Image.new("RGB", (10, 10), "white").save(source / "08_figures/umap_condition.png")
    pd.DataFrame([{"path": "07_source_data/umap_condition.tsv", "sha256": sha256(path)}]).to_csv(source / "00_audit/output_sha256.tsv", sep="\t", index=False)
    config = SimpleNamespace(root=root, resolve=lambda p: root / p,
        data={'statistics': {'control': 'WT', 'case': 'PFOS'}},
        samples=lambda: pd.DataFrame({'library_id': ['arbitrary1', 'arbitrary2'], 'condition': ['WT', 'PFOS']}))
    return config, path


def test_pure_replot_preserves_numerical_source_and_fixed_canvas(tmp_path):
    config, source = _fixture(tmp_path)
    before = sha256(source)
    assert replot(config, "runs/frozen", check=True)["status"] == "PASS"
    assert not (config.root / "result").exists()
    receipt = replot(config, "runs/frozen")
    assert receipt["statistics_recomputed"] is False
    assert receipt["figure_pairs"] == 1
    assert sha256(source) == before
    with Image.open(config.root / "result/02_atlas/plots/umap_condition.png") as im:
        assert im.mode == "RGB" and abs(im.width - 2161) <= 1


def test_changed_source_does_not_replace_existing_result(tmp_path):
    config, source = _fixture(tmp_path)
    result = config.root / "result"
    result.mkdir()
    marker = result / "keep.txt"
    marker.write_text("existing valid result")
    source.write_text(source.read_text().replace("WT", "changed"))
    with pytest.raises(RuntimeError, match="签名不一致"):
        replot(config, "runs/frozen")
    assert marker.read_text() == "existing valid result"


def test_configured_condition_names_and_library_metadata_have_stable_colors():
    from scrna_template.plotting import sample_palettes
    config = SimpleNamespace(data={'statistics': {'control': 'vehicle', 'case': 'exposed'}},
        samples=lambda: pd.DataFrame({'library_id': ['library_without_group_prefix'], 'condition': ['vehicle']}))
    conditions, libraries = sample_palettes(config)
    assert conditions == {'vehicle': '#5F8FB3', 'exposed': '#CE7C8E'}
    assert libraries['library_without_group_prefix'] == conditions['vehicle']
