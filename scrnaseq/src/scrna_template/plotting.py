from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from .figure_style import apply_style, save_figure as export_figure, COLORS, WIDTH_MM, PROFILE

apply_style()


CELL_TYPE_COLORS = PROFILE["palettes"]["cell_type"]
CLUSTER_COLORS = PROFILE["palettes"]["clusters"]


def sample_palettes(config):
    stats = config.data["statistics"]
    conditions = {str(stats["control"]): COLORS["WT"], str(stats["case"]): COLORS["PFOS"]}
    libraries = {str(row.library_id): conditions[str(row.condition)] for row in config.samples().itertuples()}
    return conditions, libraries


def category_color(category: str) -> str:
    if category in COLORS:
        return COLORS[category]
    condition = category.split("_", 1)[0]
    if condition in COLORS:
        return COLORS[condition]
    if category in CELL_TYPE_COLORS:
        return CELL_TYPE_COLORS[category]
    if category in CLUSTER_COLORS:
        return CLUSTER_COLORS[category]
    raise ValueError(f"未在项目 profile 注册类别颜色: {category}")


def ordered_categories(values):
    observed = set(map(str, values))
    named = list(COLORS) + list(CELL_TYPE_COLORS) + list(CLUSTER_COLORS)
    return [key for key in named if key in observed] + sorted(observed - set(named))




def save_figure(fig: plt.Figure, target_stem: Path, dpi: int = 300) -> list[Path]:
    # Keep panel composition and readable height; recompute margins at final width.
    import matplotlib.text
    apply_style()
    matplotlib.rcParams["savefig.bbox"] = None
    height_mm = max(90, fig.get_figheight() * 25.4)
    fig.set_size_inches(WIDTH_MM / 25.4, height_mm / 25.4, forward=True)
    for text in fig.findobj(match=matplotlib.text.Text):
        text.set_fontfamily("Liberation Sans")
        text.set_fontsize(7 if text.get_fontsize() <= 7 else 8)
    for ax in fig.axes:
        if ax.get_title():
            ax.set_title("\n".join(__import__("textwrap").wrap(ax.get_title(), width=64)), fontsize=8)
    if fig._suptitle is not None:
        fig._suptitle.set_text("\n".join(__import__("textwrap").wrap(fig._suptitle.get_text(), width=100)))
        fig._suptitle.set_fontsize(8)
    try:
        fig.tight_layout(pad=1.3)
    except ValueError:
        pass
    target_stem.parent.mkdir(parents=True, exist_ok=True)
    outputs = [target_stem.with_suffix(".png"), target_stem.with_suffix(".pdf")]
    export_figure(fig, outputs[1], height_mm=height_mm)
    plt.close(fig)
    return outputs


def categorical_embedding(embedding: np.ndarray, values: pd.Series, title: str, target_stem: Path, color_map: dict[str, str], dpi: int = 300) -> list[Path]:
    series = values.astype(str)
    unknown = set(series) - set(color_map)
    if unknown:
        raise ValueError(f"类别缺固定颜色: {sorted(unknown)}")
    categories = [name for name in color_map if name in set(series)]
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    for index, category in enumerate(categories):
        mask = series.to_numpy() == category
        color = color_map[category]
        ax.scatter(embedding[mask, 0], embedding[mask, 1], s=3, alpha=0.75, linewidths=0, color=color, label=category, rasterized=True)
    ax.set_title(title)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(markerscale=3, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    return save_figure(fig, target_stem, dpi)


def stacked_composition(frame: pd.DataFrame, target_stem: Path, dpi: int = 300) -> list[Path]:
    pivot = frame.pivot(index="library_id", columns="cell_type", values="proportion").fillna(0)
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    pivot = pivot.reindex(columns=ordered_categories(pivot.columns))
    colors = [category_color(str(category)) for category in pivot.columns]
    pivot.plot(kind="bar", stacked=True, ax=ax, color=colors, width=0.72)
    ax.set_ylabel("Captured cell proportion")
    ax.set_xlabel("Library")
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    return save_figure(fig, target_stem, dpi)
