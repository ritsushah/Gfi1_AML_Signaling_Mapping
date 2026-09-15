"""Stage 1 – Differential expression between GFI1-low and GFI1-high leukemic cells."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

from .config import (
    SECRETED_CANDIDATES, DEG_PADJ_THRESHOLD, DEG_LOGFC_THRESHOLD,
    TABLES, SENDER_TYPES,
)


def run_deg(
    expr: pd.DataFrame,
    meta: pd.DataFrame,
    group_col: str = "GFI1_status",
    low_label: str = "Low",
    high_label: str = "High",
    restrict_to_senders: bool = True,
) -> pd.DataFrame:
    """
    Wilcoxon rank-sum test per gene: GFI1-Low vs GFI1-High.

    Returns a DataFrame with log2FC, p-value, adjusted p-value, and flags.
    """
    if restrict_to_senders and "cell_type" in meta.columns:
        mask = meta["cell_type"].isin(SENDER_TYPES)
        expr = expr.loc[mask]
        meta = meta.loc[mask]

    low_idx = meta.index[meta[group_col] == low_label]
    high_idx = meta.index[meta[group_col] == high_label]

    if len(low_idx) < 10 or len(high_idx) < 10:
        raise ValueError("Too few cells in Low or High GFI1 groups for reliable DEG.")

    records = []
    for gene in expr.columns:
        a = expr.loc[low_idx, gene].values
        b = expr.loc[high_idx, gene].values
        # log2 fold-change of means
        mean_a = np.mean(a) + 1e-9
        mean_b = np.mean(b) + 1e-9
        log2fc = np.log2(mean_a / mean_b)
        try:
            stat, pval = mannwhitneyu(a, b, alternative="two-sided")
        except ValueError:
            pval = 1.0
        records.append({
            "gene": gene,
            "mean_GFI1_low": mean_a,
            "mean_GFI1_high": mean_b,
            "log2FC": log2fc,
            "pvalue": pval,
        })

    deg = pd.DataFrame(records)
    _, padj, _, _ = multipletests(deg["pvalue"], method="fdr_bh")
    deg["padj"] = padj
    deg["significant"] = (deg["padj"] < DEG_PADJ_THRESHOLD) & (deg["log2FC"] > DEG_LOGFC_THRESHOLD)
    deg["is_secreted"] = deg["gene"].isin(SECRETED_CANDIDATES)
    deg = deg.sort_values(["significant", "log2FC"], ascending=[False, False])
    return deg.reset_index(drop=True)


def prioritize_secreted(deg: pd.DataFrame) -> pd.DataFrame:
    """Return significantly up-regulated secreted candidates."""
    return deg[deg["significant"] & deg["is_secreted"]].copy()


def save_deg_results(deg: pd.DataFrame, secreted: pd.DataFrame) -> None:
    deg.to_csv(TABLES / "deg_gfi1_low_vs_high.csv", index=False)
    secreted.to_csv(TABLES / "secreted_candidates.csv", index=False)
    print(f"[DEG] Full table  → {TABLES / 'deg_gfi1_low_vs_high.csv'}")
    print(f"[DEG] Secreted    → {TABLES / 'secreted_candidates.csv'}  ({len(secreted)} genes)")
