from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


CELL_TYPE_COLORS = {
    "alpha": "#0072B2",
    "beta": "#E69F00",
    "delta": "#CC79A7",
    "endothelial": "#56B4E9",
    "immune": "#D55E00",
    "mesenchymal": "#009E73",
    "unknown": "#999999",
}


def save_figure(fig: plt.Figure, target_stem: Path, dpi: int = 300) -> list[Path]:
    target_stem.parent.mkdir(parents=True, exist_ok=True)
    outputs = [target_stem.with_suffix(".png"), target_stem.with_suffix(".pdf")]
    fig.savefig(outputs[0], dpi=dpi, bbox_inches="tight")
    fig.savefig(outputs[1], bbox_inches="tight")
    plt.close(fig)
    return outputs


def categorical_embedding(embedding: np.ndarray, values: pd.Series, title: str, target_stem: Path, dpi: int = 300) -> list[Path]:
    series = values.astype(str)
    categories = sorted(series.unique())
    cmap = plt.get_cmap("tab20")
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    for index, category in enumerate(categories):
        mask = series.to_numpy() == category
        color = CELL_TYPE_COLORS.get(category, cmap(index % 20))
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
    cmap = plt.get_cmap("tab20")
    colors = [CELL_TYPE_COLORS.get(category, cmap(index % 20)) for index, category in enumerate(pivot.columns.astype(str))]
    pivot.plot(kind="bar", stacked=True, ax=ax, color=colors, width=0.72)
    ax.set_ylabel("Captured cell proportion")
    ax.set_xlabel("Library")
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    return save_figure(fig, target_stem, dpi)
