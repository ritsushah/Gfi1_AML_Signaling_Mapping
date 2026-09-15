"""
Synthetic multi-cell-type expression data that mimics Gfi1-low AML niche remodeling.

This module generates a self-contained dataset so the full pipeline can run
without external downloads. Real-data loaders can replace these functions later.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path

from .config import RANDOM_SEED, SECRETED_CANDIDATES, DATA_PROCESSED


def generate_sc_expression(
    n_cells: int = 4000,
    n_background_genes: int = 800,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Create a synthetic single-cell expression matrix + cell metadata.

    Returns
    -------
    expr : DataFrame
        cells × genes (log-normalized-like continuous values)
    meta : DataFrame
        cell_id, cell_type, sample_id, condition, GFI1_expr, GFI1_status
    """
    rng = np.random.default_rng(RANDOM_SEED)

    # Gene universe (include receptors so LR scoring works)
    receptors = ["CXCR4", "CXCR1", "CXCR2", "CCR2", "IL6R", "FLT1", "KDR",
                 "TGFBR1", "TGFBR2", "CSF1R", "KIT", "PDGFRA", "TEK", "CD44"]
    core_genes = ["GFI1"] + SECRETED_CANDIDATES + receptors
    background = [f"BG_{i:04d}" for i in range(n_background_genes)]
    genes = core_genes + background
    n_genes = len(genes)

    # Cell-type composition
    cell_types = rng.choice(
        ["Leukemic_Blast", "HSC_MPP", "GMP", "Monocyte", "Macrophage",
         "T_cell", "MSC", "Endothelial", "Erythroid"],
        size=n_cells,
        p=[0.32, 0.07, 0.09, 0.10, 0.06, 0.08, 0.10, 0.08, 0.10],
    )

    samples = [f"AML_{i:02d}" for i in range(1, 9)] + ["Healthy_01", "Healthy_02"]
    sample_id = rng.choice(samples, size=n_cells)
    condition = np.where(pd.Series(sample_id).str.startswith("AML"), "AML", "Healthy")

    # Base expression
    X = rng.exponential(scale=0.6, size=(n_cells, n_genes))

    # Gene index helpers
    gfi1_idx = genes.index("GFI1")
    secreted_idx = {g: genes.index(g) for g in SECRETED_CANDIDATES if g in genes}

    # Make GFI1 lower in a subset of leukemic cells
    leukemic_mask = np.isin(cell_types, ["Leukemic_Blast", "HSC_MPP", "GMP"])
    low_gfi1_mask = leukemic_mask & (rng.random(n_cells) < 0.45)
    X[low_gfi1_mask, gfi1_idx] *= 0.25          # strong down-regulation
    X[~low_gfi1_mask & leukemic_mask, gfi1_idx] *= 1.4

    # When GFI1 is low, up-regulate a subset of secreted factors (the biological signal)
    upregulated = ["CXCL12", "CXCL8", "CCL2", "IL6", "VEGFA", "TGFB1",
                   "CSF1", "SPP1", "MMP9", "ANGPT2", "S100A8", "S100A9"]
    for g in upregulated:
        if g in secreted_idx:
            X[low_gfi1_mask, secreted_idx[g]] *= rng.uniform(2.5, 4.5, size=low_gfi1_mask.sum())

    # Add modest cell-type specific baselines + receptor expression on niche cells
    for ct, boost_genes in [
        ("MSC", ["CXCL12", "KITLG", "ANGPT1", "CXCR4", "TGFBR1", "CD44"]),
        ("Endothelial", ["VEGFA", "ANGPT2", "KDR", "FLT1", "TEK"]),
        ("Monocyte", ["IL1B", "TNF", "CCL3", "CSF1R", "CCR2", "IL6R"]),
        ("Macrophage", ["CSF1R", "CD44", "TGFBR2"]),
    ]:
        mask = cell_types == ct
        for g in boost_genes:
            if g in genes:
                idx = genes.index(g)
                X[mask, idx] *= 2.2

    # Log-like transform + small noise
    X = np.log1p(X) + rng.normal(0, 0.05, size=X.shape)
    X = np.clip(X, 0, None)

    expr = pd.DataFrame(X, columns=genes,
                        index=[f"Cell_{i:05d}" for i in range(n_cells)])

    # GFI1 status
    gfi1_vals = expr["GFI1"].values
    low_t = np.quantile(gfi1_vals, 0.30)
    high_t = np.quantile(gfi1_vals, 0.70)
    status = np.full(n_cells, "Intermediate", dtype=object)
    status[gfi1_vals <= low_t] = "Low"
    status[gfi1_vals >= high_t] = "High"

    meta = pd.DataFrame({
        "cell_id": expr.index,
        "cell_type": cell_types,
        "sample_id": sample_id,
        "condition": condition,
        "GFI1_expr": gfi1_vals,
        "GFI1_status": status,
    }).set_index("cell_id")

    return expr, meta


def generate_bulk_clinical(n_patients: int = 250) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Synthetic bulk expression + clinical outcome table for ML prioritization.
    """
    rng = np.random.default_rng(RANDOM_SEED + 1)

    # Use the same secreted genes as features
    genes = SECRETED_CANDIDATES.copy()
    n_genes = len(genes)

    X = rng.normal(loc=0, scale=1, size=(n_patients, n_genes))

    # Ground-truth drivers that influence risk
    true_drivers = ["CXCL12", "IL6", "VEGFA", "TGFB1", "CSF1", "SPP1", "MMP9"]
    driver_idx = [genes.index(g) for g in true_drivers if g in genes]

    risk_score = X[:, driver_idx].sum(axis=1) + rng.normal(0, 0.4, n_patients)
    high_risk = (risk_score > np.median(risk_score)).astype(int)
    # Higher risk → shorter survival
    survival_time = np.exp(3.2 - 0.45 * risk_score) + rng.exponential(0.3, n_patients)
    event = rng.binomial(1, 0.65, n_patients)

    expr = pd.DataFrame(X, columns=genes,
                        index=[f"Patient_{i:03d}" for i in range(n_patients)])
    clinical = pd.DataFrame({
        "patient_id": expr.index,
        "high_risk": high_risk,
        "survival_months": survival_time,
        "event": event,
        "age": rng.integers(22, 78, n_patients),
        "sex": rng.choice(["M", "F"], n_patients),
    }).set_index("patient_id")

    return expr, clinical, true_drivers


def save_demo_data() -> dict[str, Path]:
    """Generate and persist demo datasets. Returns paths."""
    expr_sc, meta_sc = generate_sc_expression()
    expr_bulk, clin, drivers = generate_bulk_clinical()

    paths = {}
    paths["sc_expr"] = DATA_PROCESSED / "demo_sc_expression.csv"
    paths["sc_meta"] = DATA_PROCESSED / "demo_sc_metadata.csv"
    paths["bulk_expr"] = DATA_PROCESSED / "demo_bulk_expression.csv"
    paths["bulk_clin"] = DATA_PROCESSED / "demo_bulk_clinical.csv"
    paths["drivers"] = DATA_PROCESSED / "demo_true_drivers.txt"

    expr_sc.to_csv(paths["sc_expr"])
    meta_sc.to_csv(paths["sc_meta"])
    expr_bulk.to_csv(paths["bulk_expr"])
    clin.to_csv(paths["bulk_clin"])
    paths["drivers"].write_text("\n".join(drivers))

    return paths
