"""Compile a polished Markdown summary report from pipeline outputs."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pandas as pd
from tabulate import tabulate

from .config import RESULTS, TABLES, FIGURES


def write_summary_report(
    deg: pd.DataFrame,
    secreted: pd.DataFrame,
    importance: pd.Series,
    metrics: dict,
    interactions: pd.DataFrame,
    true_drivers: list[str] | None = None,
) -> Path:
    """Generate RESULTS/SUMMARY_REPORT.md."""

    lines = []
    lines.append("# Gfi1-Regulated Intercellular Signaling in AML")
    lines.append("## Computational Pipeline Summary Report")
    lines.append("")
    lines.append(f"**Generated:** {datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')}")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("### 1. Hypothesis")
    lines.append("")
    lines.append(
        "Loss or down-regulation of the transcriptional repressor **Gfi1** in acute myeloid "
        "leukemia (AML) cells leads to de-repression of secreted ligands (cytokines, chemokines, "
        "extracellular-matrix factors). These ligands engage receptors on bone-marrow stromal "
        "and endothelial cells, remodeling the niche into a protective microenvironment."
    )
    lines.append("")
    lines.append("This report summarizes an end-to-end *in-silico* mapping of that axis.")
    lines.append("")

    # ----- DEG -----
    lines.append("### 2. Differential Expression (GFI1-Low vs GFI1-High)")
    lines.append("")
    n_sig = deg["significant"].sum()
    n_sec = len(secreted)
    lines.append(f"- Total genes tested: **{len(deg)}**")
    lines.append(f"- Significantly up-regulated (padj < 0.05 & log2FC > 0.5): **{n_sig}**")
    lines.append(f"- Of which annotated as secreted / extracellular: **{n_sec}**")
    lines.append("")
    if n_sec > 0:
        lines.append("**Top secreted candidates**")
        lines.append("")
        top = secreted.head(10)[["gene", "log2FC", "padj", "mean_GFI1_low", "mean_GFI1_high"]]
        lines.append(tabulate(top, headers="keys", tablefmt="github", floatfmt=".3f", showindex=False))
        lines.append("")
        lines.append("![Top secreted](figures/top_secreted_candidates.png)")
        lines.append("")

    # ----- ML -----
    lines.append("### 3. Machine-Learning Prioritization")
    lines.append("")
    lines.append(
        "An XGBoost classifier was trained to predict high-risk clinical status from the "
        "expression of candidate secreted genes."
    )
    lines.append("")
    lines.append("| Metric | Value |")
    lines.append("|--------|-------|")
    lines.append(f"| Test AUC | {metrics['test_auc']:.3f} |")
    lines.append(f"| Test Average Precision | {metrics['test_ap']:.3f} |")
    lines.append(f"| 5-fold CV AUC | {metrics['cv_auc_mean']:.3f} ± {metrics['cv_auc_std']:.3f} |")
    lines.append("")
    lines.append("**Top genes by feature importance**")
    lines.append("")
    top_imp = importance.head(12).reset_index()
    top_imp.columns = ["gene", "importance"]
    lines.append(tabulate(top_imp, headers="keys", tablefmt="github", floatfmt=".4f", showindex=False))
    lines.append("")
    lines.append("![Feature importance](figures/ml_feature_importance.png)")
    lines.append("")

    if true_drivers:
        recovered = sorted(set(importance.head(10).index) & set(true_drivers))
        lines.append(f"*Demo recovery check:* {len(recovered)}/{len(true_drivers)} planted drivers recovered in top-10 ({', '.join(recovered)}).")
        lines.append("")

    # ----- Communication -----
    lines.append("### 4. Intercellular Ligand–Receptor Map")
    lines.append("")
    lines.append(
        "Ligand–receptor scores were computed between GFI1-Low leukemic/progenitor cells "
        "(senders) and niche populations (MSC, endothelial, monocyte/macrophage)."
    )
    lines.append("")
    if len(interactions) > 0:
        lines.append(f"- Scored interactions: **{len(interactions)}**")
        lines.append("")
        lines.append("**Highest-scoring pairs**")
        lines.append("")
        top_lr = interactions.head(12)[
            ["ligand", "receptor", "sender", "receiver", "interaction_score"]
        ]
        lines.append(tabulate(top_lr, headers="keys", tablefmt="github", floatfmt=".3f", showindex=False))
        lines.append("")
        lines.append("![LR heatmap](figures/lr_interaction_heatmap.png)")
        lines.append("")
    else:
        lines.append("*No interactions passed the expression threshold.*")
        lines.append("")

    # ----- Structural note -----
    lines.append("### 5. Structural Modeling & Virtual Screening (Next Step)")
    lines.append("")
    lines.append(
        "The highest-priority ligand–receptor pairs identified above are candidates for "
        "AlphaFold3 / ColabFold complex prediction followed by AutoDock Vina virtual screening "
        "of FDA-approved compound libraries. See `docs/structural_followup.md` for a concrete protocol."
    )
    lines.append("")

    # ----- Files -----
    lines.append("### 6. Output Files")
    lines.append("")
    lines.append("| File | Description |")
    lines.append("|------|-------------|")
    lines.append("| `results/tables/deg_gfi1_low_vs_high.csv` | Full differential-expression table |")
    lines.append("| `results/tables/secreted_candidates.csv` | Prioritized secreted ligands |")
    lines.append("| `results/tables/ml_feature_importance.csv` | ML ranking of clinical drivers |")
    lines.append("| `results/tables/ml_performance.csv` | Model performance metrics |")
    lines.append("| `results/tables/lr_interaction_scores.csv` | Ligand–receptor interaction scores |")
    lines.append("| `results/figures/*.png` | All publication-ready figures |")
    lines.append("| `results/models/xgboost_risk_model.joblib` | Trained risk model |")
    lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("*This report was auto-generated by the Gfi1-AML Signaling Mapping pipeline.*")
    lines.append("")

    report_path = RESULTS / "SUMMARY_REPORT.md"
    report_path.write_text("\n".join(lines))
    return report_path
