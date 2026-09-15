#!/usr/bin/env python3
"""
End-to-end pipeline for Gfi1-regulated intercellular signaling mapping in AML.

Usage
-----
    python -m src.run_pipeline
    # or
    python src/run_pipeline.py
"""

from __future__ import annotations

import sys
from pathlib import Path

# Allow running as script
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.data_generation import save_demo_data, generate_sc_expression, generate_bulk_clinical
from src.differential_expression import run_deg, prioritize_secreted, save_deg_results
from src.ml_prioritization import train_risk_model, save_ml_results
from src.cell_communication import score_interactions, save_communication_results
from src.visualize import (
    plot_gfi1_distribution, plot_top_secreted,
    plot_feature_importance, plot_lr_heatmap,
)
from src.config import TABLES, FIGURES, RESULTS
from src.report import write_summary_report

import pandas as pd


def main():
    print("=" * 60)
    print("  Gfi1-AML Intercellular Signaling Mapping Pipeline")
    print("=" * 60)

    # ------------------------------------------------------------------
    # 0. Data
    # ------------------------------------------------------------------
    print("\n[0] Generating self-contained demo datasets...")
    paths = save_demo_data()
    expr_sc = pd.read_csv(paths["sc_expr"], index_col=0)
    meta_sc = pd.read_csv(paths["sc_meta"], index_col=0)
    expr_bulk = pd.read_csv(paths["bulk_expr"], index_col=0)
    clin = pd.read_csv(paths["bulk_clin"], index_col=0)
    true_drivers = paths["drivers"].read_text().strip().splitlines()
    print(f"    scRNA-like matrix : {expr_sc.shape[0]} cells × {expr_sc.shape[1]} genes")
    print(f"    Bulk clinical set : {expr_bulk.shape[0]} patients × {expr_bulk.shape[1]} genes")

    # ------------------------------------------------------------------
    # 1. Differential expression
    # ------------------------------------------------------------------
    print("\n[1] Differential expression (GFI1-Low vs High leukemic cells)...")
    deg = run_deg(expr_sc, meta_sc)
    secreted = prioritize_secreted(deg)
    save_deg_results(deg, secreted)
    print(f"    Significant secreted candidates: {len(secreted)}")
    if len(secreted):
        print("    Top 5:", ", ".join(secreted["gene"].head(5).tolist()))

    # ------------------------------------------------------------------
    # 2. Machine-learning prioritization
    # ------------------------------------------------------------------
    print("\n[2] Machine-learning prioritization of clinical drivers...")
    pipe, importance, metrics = train_risk_model(expr_bulk, clin)
    save_ml_results(pipe, importance, metrics)

    # Recovery check (demo only)
    recovered = set(importance.head(10).index) & set(true_drivers)
    print(f"    Recovered true drivers in top-10: {sorted(recovered)}")

    # ------------------------------------------------------------------
    # 3. Cell–cell communication
    # ------------------------------------------------------------------
    print("\n[3] Scoring ligand–receptor interactions (GFI1-Low → niche)...")
    interactions = score_interactions(expr_sc, meta_sc)
    save_communication_results(interactions)
    if len(interactions):
        print("    Top interaction:", 
              f"{interactions.iloc[0]['ligand']} → {interactions.iloc[0]['receptor']} "
              f"({interactions.iloc[0]['sender']} → {interactions.iloc[0]['receiver']})")

    # ------------------------------------------------------------------
    # 4. Figures
    # ------------------------------------------------------------------
    print("\n[4] Generating figures...")
    plot_gfi1_distribution(meta_sc)
    if len(secreted):
        plot_top_secreted(secreted)
    plot_feature_importance(importance)
    plot_lr_heatmap(interactions)
    print(f"    Figures written to {FIGURES}")

    # ------------------------------------------------------------------
    # 5. Compile polished report
    # ------------------------------------------------------------------
    print("\n[5] Compiling summary report...")
    report_path = write_summary_report(
        deg=deg,
        secreted=secreted,
        importance=importance,
        metrics=metrics,
        interactions=interactions,
        true_drivers=true_drivers,
    )
    print(f"    Report → {report_path}")

    print("\n" + "=" * 60)
    print("  Pipeline finished successfully.")
    print(f"  Results directory: {RESULTS}")
    print("=" * 60)


if __name__ == "__main__":
    main()
