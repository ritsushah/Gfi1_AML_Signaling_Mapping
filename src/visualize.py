"""Publication-quality figures for the pipeline."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

from .config import FIGURES

sns.set_theme(style="whitegrid", context="paper", font_scale=1.15)
plt.rcParams["figure.dpi"] = 150
plt.rcParams["savefig.bbox"] = "tight"
plt.rcParams["savefig.facecolor"] = "white"


def plot_gfi1_distribution(meta: pd.DataFrame, outfile: str = "gfi1_status_distribution.png") -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # Histogram of GFI1 expression
    for status, color in [("Low", "#d62728"), ("Intermediate", "#7f7f7f"), ("High", "#2ca02c")]:
        subset = meta.loc[meta["GFI1_status"] == status, "GFI1_expr"]
        axes[0].hist(subset, bins=30, alpha=0.7, label=status, color=color)
    axes[0].set_xlabel("GFI1 expression")
    axes[0].set_ylabel("Number of cells")
    axes[0].set_title("GFI1 expression by status")
    axes[0].legend()

    # Cell-type composition of GFI1-Low cells
    low = meta[meta["GFI1_status"] == "Low"]
    ct_counts = low["cell_type"].value_counts()
    axes[1].barh(ct_counts.index, ct_counts.values, color="#1f77b4")
    axes[1].set_xlabel("Cell count")
    axes[1].set_title("Cell-type composition of GFI1-Low cells")
    axes[1].invert_yaxis()

    path = FIGURES / outfile
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_top_secreted(secreted: pd.DataFrame, top_n: int = 12, outfile: str = "top_secreted_candidates.png") -> Path:
    df = secreted.head(top_n).copy()
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ["#d62728" if x > 1 else "#ff7f0e" for x in df["log2FC"]]
    ax.barh(df["gene"], df["log2FC"], color=colors)
    ax.set_xlabel("log2 fold-change (GFI1-Low / GFI1-High)")
    ax.set_title(f"Top {top_n} up-regulated secreted candidates")
    ax.axvline(0, color="black", linewidth=0.8)
    ax.invert_yaxis()
    path = FIGURES / outfile
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_feature_importance(importance: pd.Series, top_n: int = 15, outfile: str = "ml_feature_importance.png") -> Path:
    top = importance.head(top_n)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(top.index, top.values, color="#2ca02c")
    ax.set_xlabel("XGBoost feature importance")
    ax.set_title(f"Top {top_n} genes ranked by predictive importance for high-risk status")
    ax.invert_yaxis()
    path = FIGURES / outfile
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_lr_heatmap(interactions: pd.DataFrame, outfile: str = "lr_interaction_heatmap.png") -> Path | None:
    if interactions.empty:
        return None
    # Pivot to ligand×receptor mean score (averaged over cell-type pairs)
    pivot = interactions.groupby(["ligand", "receptor"])["interaction_score"].mean().unstack(fill_value=0)
    fig, ax = plt.subplots(figsize=(9, 6))
    sns.heatmap(pivot, annot=True, fmt=".2f", cmap="YlOrRd", ax=ax, linewidths=0.5)
    ax.set_title("Ligand–Receptor interaction scores\n(GFI1-Low leukemic → niche cells)")
    path = FIGURES / outfile
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path
