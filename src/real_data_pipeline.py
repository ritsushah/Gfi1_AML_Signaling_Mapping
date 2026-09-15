#!/usr/bin/env python3
"""
Real-data pipeline for Gfi1-AML signaling mapping.

Uses:
  - TCGA-LAML (Xena HiSeqV2 expression + survival) for clinical prioritization
  - GSE116256 (van Galen et al., Cell 2019) dem+anno matrices for
    GFI1-stratified differential expression and ligand–receptor scoring
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import xgboost as xgb
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.config import (
    ROOT, DATA_RAW, DATA_PROCESSED, RESULTS, TABLES, FIGURES, MODELS,
    SECRETED_CANDIDATES, LR_PAIRS, RANDOM_SEED,
)

sns.set_theme(style="whitegrid", context="paper", font_scale=1.1)
plt.rcParams["figure.dpi"] = 140
plt.rcParams["savefig.bbox"] = "tight"


# ---------------------------------------------------------------------------
# 1. Load TCGA-LAML
# ---------------------------------------------------------------------------

def load_tcga_laml():
    """Load expression (genes x samples) and survival from Xena downloads."""
    expr_path = DATA_RAW / "LAML_gene_expression.gz"
    surv_path = DATA_RAW / "LAML_survival.txt"
    clin_path = DATA_RAW / "LAML_clinical.gz"

    print("[TCGA] Loading expression (subset of genes to save memory)...")
    # Only keep secreted candidates + GFI1 + a few receptors
    keep_genes = set(SECRETED_CANDIDATES + ["GFI1", "CXCR4", "IL6R", "KDR", "FLT1",
                     "CSF1R", "TGFBR1", "TGFBR2", "KIT", "CD44", "CCR2"])
    # Read in chunks / filter rows
    chunks = []
    for chunk in pd.read_csv(expr_path, sep="\t", index_col=0, compression="gzip", chunksize=2000):
        overlap = [g for g in chunk.index if g in keep_genes]
        if overlap:
            chunks.append(chunk.loc[overlap])
    if not chunks:
        # fallback: load first 500 genes
        expr = pd.read_csv(expr_path, sep="\t", index_col=0, compression="gzip", nrows=500)
    else:
        expr = pd.concat(chunks)
    expr = expr[~expr.index.duplicated(keep="first")]
    expr = expr.T  # samples x genes
    expr.index = expr.index.str[:15]

    print("[TCGA] Loading survival...")
    surv = pd.read_csv(surv_path, sep="\t")
    # columns typically: sample, OS, OS.time
    surv.columns = [c.strip() for c in surv.columns]
    if "sample" in surv.columns:
        surv = surv.set_index("sample")
    surv.index = surv.index.str[:15]

    # Align
    common = expr.index.intersection(surv.index)
    expr = expr.loc[common]
    surv = surv.loc[common]
    print(f"[TCGA] Aligned {len(common)} samples, {expr.shape[1]} genes")

    # High-risk proxy: shorter survival / event
    # Use median OS.time among events, or simply top/bottom quartile of risk
    if "OS.time" in surv.columns:
        time_col = "OS.time"
    elif "OS_Time" in surv.columns:
        time_col = "OS_Time"
    else:
        time_col = [c for c in surv.columns if "time" in c.lower()][0]

    if "OS" in surv.columns:
        event_col = "OS"
    else:
        event_col = [c for c in surv.columns if c.lower() in ("os", "status", "event")][0]

    # Binary high-risk: died with short survival (below median among events)
    times = surv[time_col].astype(float)
    events = surv[event_col].astype(float)
    med = times[events == 1].median()
    high_risk = ((events == 1) & (times <= med)).astype(int)

    clinical = pd.DataFrame({
        "survival_days": times,
        "event": events,
        "high_risk": high_risk,
    }, index=surv.index)

    return expr, clinical


# ---------------------------------------------------------------------------
# 2. Load GSE116256 subset
# ---------------------------------------------------------------------------

SAMPLE_SPECS = [
    # (dem_file, anno_file, sample_id, condition)
    ("GSM3587923_AML1012-D0.dem.txt.gz", "GSM3587924_AML1012-D0.anno.txt.gz", "AML1012", "AML"),
    ("GSM3587925_AML210A-D0.dem.txt.gz", "GSM3587926_AML210A-D0.anno.txt.gz", "AML210A", "AML"),
    ("GSM3587931_AML328-D0.dem.txt.gz", "GSM3587932_AML328-D0.anno.txt.gz", "AML328", "AML"),
    ("GSM3587996_BM1.dem.txt.gz", "GSM3587996_BM1.anno.txt.gz", "BM1", "Healthy"),
    ("GSM3587997_BM2.dem.txt.gz", "GSM3587997_BM2.anno.txt.gz", "BM2", "Healthy"),
    ("GSM3587998_BM3.dem.txt.gz", "GSM3587999_BM3.anno.txt.gz", "BM3", "Healthy"),
]


def _read_dem(path: Path) -> pd.DataFrame:
    """genes x cells count matrix."""
    df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
    return df


def _read_anno(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", compression="gzip")
    if "Cell" in df.columns:
        df = df.set_index("Cell")
    return df


def load_gse116256_subset(max_cells_per_sample: int = 800):
    """Load and merge a subset of AML + healthy BM samples."""
    expr_parts = []
    meta_parts = []

    for dem_name, anno_name, sample_id, condition in SAMPLE_SPECS:
        dem_path = DATA_RAW / dem_name
        anno_path = DATA_RAW / anno_name
        if not dem_path.exists() or not anno_path.exists():
            print(f"  [skip] missing {dem_name} or {anno_name}")
            continue

        print(f"  Loading {sample_id}...")
        dem = _read_dem(dem_path)  # genes x cells
        anno = _read_anno(anno_path)

        # Align cells
        cells = dem.columns.intersection(anno.index)
        if len(cells) == 0:
            # try stripping prefixes
            dem.columns = [c.split("_")[-1] if "_" in c else c for c in dem.columns]
            cells = dem.columns.intersection(anno.index)
        dem = dem[cells]
        anno = anno.loc[cells]

        # Subsample if large
        if dem.shape[1] > max_cells_per_sample:
            rng = np.random.default_rng(RANDOM_SEED)
            keep = rng.choice(dem.columns, size=max_cells_per_sample, replace=False)
            dem = dem[keep]
            anno = anno.loc[keep]

        # Transpose to cells x genes
        mat = dem.T
        mat.index = [f"{sample_id}_{c}" for c in mat.index]
        anno = anno.copy()
        anno.index = mat.index
        anno["sample_id"] = sample_id
        anno["condition"] = condition

        # Normalize cell-type column name
        if "CellType" in anno.columns:
            anno["cell_type"] = anno["CellType"]
        elif "PredictionRefined" in anno.columns:
            anno["cell_type"] = anno["PredictionRefined"]
        else:
            anno["cell_type"] = "Unknown"

        expr_parts.append(mat)
        meta_parts.append(anno)

    if not expr_parts:
        raise RuntimeError("No GSE116256 samples could be loaded")

    # Outer join on genes
    expr = pd.concat(expr_parts, axis=0, join="outer").fillna(0)
    meta = pd.concat(meta_parts, axis=0)
    meta = meta.loc[expr.index]

    # Keep genes present in most cells / interesting list
    # Log-normalize like: log1p(CPM-ish)
    lib = expr.sum(axis=1).replace(0, np.nan)
    expr_norm = expr.div(lib, axis=0) * 1e4
    expr_norm = np.log1p(expr_norm)

    print(f"[GSE116256] {expr_norm.shape[0]} cells × {expr_norm.shape[1]} genes")
    return expr_norm, meta


# ---------------------------------------------------------------------------
# 3. Analyses
# ---------------------------------------------------------------------------

LEUKEMIC_LIKE = {
    "HSC", "Prog", "GMP", "ProMono", "HSC-like", "Prog-like", "GMP-like",
    "ProMono-like", "EarlyEry", "earlyEry", "cDC", "pDC",
}


def run_real_deg(expr: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """GFI1-low vs high within leukemic-like cells."""
    gfi1_col = None
    for cand in ["GFI1", "Gfi1", "gfi1"]:
        if cand in expr.columns:
            gfi1_col = cand
            break
    if gfi1_col is None:
        raise ValueError("GFI1 not found in expression matrix")

    gfi1 = expr[gfi1_col]
    low_t = gfi1.quantile(0.30)
    high_t = gfi1.quantile(0.70)
    status = pd.Series("Intermediate", index=expr.index)
    status[gfi1 <= low_t] = "Low"
    status[gfi1 >= high_t] = "High"
    meta = meta.copy()
    meta["GFI1_status"] = status
    meta["GFI1_expr"] = gfi1

    # Restrict to leukemic-like + AML samples
    ct = meta["cell_type"].astype(str)
    is_leuk = ct.str.contains("like|HSC|Prog|GMP|ProMono|cDC|pDC|EarlyEry|earlyEry", case=False, regex=True)
    is_aml = meta["condition"] == "AML"
    mask = is_leuk | is_aml  # be inclusive
    sub_expr = expr.loc[mask]
    sub_meta = meta.loc[mask]

    low_idx = sub_meta.index[sub_meta["GFI1_status"] == "Low"]
    high_idx = sub_meta.index[sub_meta["GFI1_status"] == "High"]
    print(f"[DEG] Low={len(low_idx)}, High={len(high_idx)} cells")

    # Focus on secreted + receptors + GFI1 for speed, plus top variable
    focus = [g for g in SECRETED_CANDIDATES + ["GFI1", "CXCR4", "CXCR1", "CXCR2",
              "CCR2", "IL6R", "FLT1", "KDR", "TGFBR1", "TGFBR2", "CSF1R", "KIT",
              "PDGFRA", "TEK", "CD44"] if g in sub_expr.columns]
    # add some highly variable
    var = sub_expr.var().nlargest(500).index.tolist()
    genes = list(dict.fromkeys(focus + var))

    records = []
    for gene in genes:
        a = sub_expr.loc[low_idx, gene].values
        b = sub_expr.loc[high_idx, gene].values
        mean_a, mean_b = a.mean() + 1e-9, b.mean() + 1e-9
        log2fc = np.log2(mean_a / mean_b)
        try:
            _, pval = mannwhitneyu(a, b, alternative="two-sided")
        except ValueError:
            pval = 1.0
        records.append({
            "gene": gene, "mean_GFI1_low": mean_a, "mean_GFI1_high": mean_b,
            "log2FC": log2fc, "pvalue": pval,
        })

    deg = pd.DataFrame(records)
    _, padj, _, _ = multipletests(deg["pvalue"], method="fdr_bh")
    deg["padj"] = padj
    deg["significant"] = (deg["padj"] < 0.05) & (deg["log2FC"] > 0.25)
    deg["is_secreted"] = deg["gene"].isin(SECRETED_CANDIDATES)
    deg = deg.sort_values(["significant", "log2FC"], ascending=[False, False])
    return deg, meta


def run_real_ml(expr: pd.DataFrame, clinical: pd.DataFrame):
    """XGBoost on TCGA using secreted-gene features."""
    genes = [g for g in SECRETED_CANDIDATES if g in expr.columns]
    if len(genes) < 5:
        # fallback: any gene with GFI1-related neighbors + top variance
        genes = expr.var().nlargest(50).index.tolist()
        if "GFI1" in expr.columns:
            genes = ["GFI1"] + genes

    X = expr[genes].copy()
    y = clinical.loc[X.index, "high_risk"]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.25, random_state=RANDOM_SEED, stratify=y
    )

    model = xgb.XGBClassifier(
        n_estimators=200, max_depth=3, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8, random_state=RANDOM_SEED,
        eval_metric="logloss", n_jobs=-1,
    )
    pipe = Pipeline([("scaler", StandardScaler()), ("clf", model)])
    pipe.fit(X_train, y_train)

    proba = pipe.predict_proba(X_test)[:, 1]
    metrics = {
        "test_auc": float(roc_auc_score(y_test, proba)),
        "test_ap": float(average_precision_score(y_test, proba)),
        "n_train": len(X_train),
        "n_test": len(X_test),
        "n_genes": len(genes),
    }
    cv = cross_val_score(pipe, X, y, cv=5, scoring="roc_auc", n_jobs=-1)
    metrics["cv_auc_mean"] = float(cv.mean())
    metrics["cv_auc_std"] = float(cv.std())

    imp = pd.Series(pipe.named_steps["clf"].feature_importances_, index=genes).sort_values(ascending=False)
    return pipe, imp, metrics


def run_real_lr(expr: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """Ligand–receptor scores from GFI1-low AML cells to niche-like populations."""
    # Map van Galen cell types to sender/receiver
    sender_kw = ["HSC", "Prog", "GMP", "ProMono", "like"]
    receiver_map = {
        "Mono": "Monocyte",
        "cDC": "Dendritic",
        "T": "T_cell",
        "CTL": "T_cell",
        "NK": "NK",
        "B": "B_cell",
        "lateEry": "Erythroid",
        "earlyEry": "Erythroid",
    }

    if "GFI1_status" not in meta.columns:
        gfi1 = expr["GFI1"] if "GFI1" in expr.columns else None
        if gfi1 is not None:
            meta = meta.copy()
            meta["GFI1_status"] = np.where(gfi1 <= gfi1.quantile(0.3), "Low", "Other")

    records = []
    for ligand, receptor in LR_PAIRS:
        if ligand not in expr.columns or receptor not in expr.columns:
            continue
        for sender_key in sender_kw:
            sender_mask = meta["cell_type"].astype(str).str.contains(sender_key, case=False, na=False)
            sender_mask &= meta.get("GFI1_status", "Low") == "Low"
            sender_mask &= meta["condition"] == "AML"
            if sender_mask.sum() < 5:
                continue
            lig_mean = expr.loc[sender_mask, ligand].mean()

            for rkey, rlabel in receiver_map.items():
                recv_mask = meta["cell_type"].astype(str).str.contains(rkey, case=False, na=False)
                if recv_mask.sum() < 5:
                    continue
                rec_mean = expr.loc[recv_mask, receptor].mean()
                if lig_mean > 0.05 and rec_mean > 0.05:
                    records.append({
                        "ligand": ligand, "receptor": receptor,
                        "sender": sender_key, "receiver": rlabel,
                        "ligand_mean": round(lig_mean, 4),
                        "receptor_mean": round(rec_mean, 4),
                        "interaction_score": round(lig_mean * rec_mean, 4),
                    })

    df = pd.DataFrame(records)
    if len(df):
        df = df.sort_values("interaction_score", ascending=False).reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# 4. Figures + report
# ---------------------------------------------------------------------------

def save_figures(deg, secreted, importance, interactions, meta):
    FIGURES.mkdir(parents=True, exist_ok=True)

    # GFI1 distribution
    if "GFI1_expr" in meta.columns:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        for st, c in [("Low", "#d62728"), ("Intermediate", "#7f7f7f"), ("High", "#2ca02c")]:
            sub = meta.loc[meta.get("GFI1_status") == st, "GFI1_expr"]
            if len(sub):
                axes[0].hist(sub, bins=30, alpha=0.7, label=st, color=c)
        axes[0].set_xlabel("GFI1 expression"); axes[0].set_ylabel("Cells")
        axes[0].set_title("GFI1 expression (GSE116256)"); axes[0].legend()
        ct = meta["cell_type"].value_counts().head(12)
        axes[1].barh(ct.index.astype(str), ct.values, color="#1f77b4")
        axes[1].set_title("Cell-type composition"); axes[1].invert_yaxis()
        fig.tight_layout()
        fig.savefig(FIGURES / "gfi1_status_distribution.png")
        plt.close(fig)

    if len(secreted):
        top = secreted.head(12)
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.barh(top["gene"], top["log2FC"], color="#d62728")
        ax.set_xlabel("log2FC (GFI1-Low / High)"); ax.set_title("Up-regulated secreted candidates (real data)")
        ax.invert_yaxis()
        fig.tight_layout()
        fig.savefig(FIGURES / "top_secreted_candidates.png")
        plt.close(fig)

    if len(importance):
        top = importance.head(15)
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.barh(top.index, top.values, color="#2ca02c")
        ax.set_xlabel("XGBoost importance"); ax.set_title("TCGA-LAML: genes predicting high-risk status")
        ax.invert_yaxis()
        fig.tight_layout()
        fig.savefig(FIGURES / "ml_feature_importance.png")
        plt.close(fig)

    if len(interactions):
        pivot = interactions.groupby(["ligand", "receptor"])["interaction_score"].mean().unstack(fill_value=0)
        fig, ax = plt.subplots(figsize=(9, 6))
        sns.heatmap(pivot, annot=True, fmt=".2f", cmap="YlOrRd", ax=ax)
        ax.set_title("Ligand–receptor scores (GFI1-Low AML → niche)")
        fig.tight_layout()
        fig.savefig(FIGURES / "lr_interaction_heatmap.png")
        plt.close(fig)


def write_real_report(deg, secreted, importance, metrics, interactions):
    lines = [
        "# Gfi1-AML Signaling Mapping — Real Data Results",
        "",
        "**Data sources**",
        "- Single-cell: GSE116256 (van Galen et al., Cell 2019) — AML diagnosis samples + healthy BM",
        "- Bulk clinical: TCGA-LAML (UCSC Xena HiSeqV2 expression + survival)",
        "",
        "---",
        "",
        "## Differential expression (GFI1-Low vs High, leukemic-like cells)",
        "",
        f"- Genes tested (focused set): {len(deg)}",
        f"- Significant up-regulated (padj<0.05, log2FC>0.25): **{deg['significant'].sum()}**",
        f"- Secreted among significant: **{len(secreted)}**",
        "",
    ]
    if len(secreted):
        lines.append("### Top secreted candidates")
        lines.append("")
        lines.append(secreted.head(15)[["gene", "log2FC", "padj", "mean_GFI1_low", "mean_GFI1_high"]].to_markdown(index=False))
        lines.append("")
        lines.append("![secreted](figures/top_secreted_candidates.png)")
        lines.append("")

    lines += [
        "## Machine learning on TCGA-LAML",
        "",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| Test AUC | {metrics['test_auc']:.3f} |",
        f"| Test AP | {metrics['test_ap']:.3f} |",
        f"| CV AUC | {metrics['cv_auc_mean']:.3f} ± {metrics['cv_auc_std']:.3f} |",
        f"| Samples (train/test) | {metrics['n_train']}/{metrics['n_test']} |",
        f"| Features (genes) | {metrics['n_genes']} |",
        "",
        "### Top predictive genes",
        "",
        importance.head(15).to_frame("importance").to_markdown(),
        "",
        "![ml](figures/ml_feature_importance.png)",
        "",
        "## Ligand–receptor communication",
        "",
        f"- Scored interactions: **{len(interactions)}**",
        "",
    ]
    if len(interactions):
        lines.append(interactions.head(15).to_markdown(index=False))
        lines.append("")
        lines.append("![lr](figures/lr_interaction_heatmap.png)")
        lines.append("")

    lines += [
        "---",
        "",
        "*Generated by `src/real_data_pipeline.py` using public TCGA-LAML and GSE116256 data.*",
    ]
    path = RESULTS / "SUMMARY_REPORT_REAL.md"
    path.write_text("\n".join(lines))
    return path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    MODELS.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("  Real-data Gfi1-AML Signaling Pipeline")
    print("=" * 60)

    # TCGA
    print("\n[1] TCGA-LAML bulk + survival")
    expr_bulk, clin = load_tcga_laml()
    pipe, importance, metrics = run_real_ml(expr_bulk, clin)
    importance.to_csv(TABLES / "ml_feature_importance_real.csv")
    pd.DataFrame([metrics]).to_csv(TABLES / "ml_performance_real.csv", index=False)
    print(f"  Test AUC={metrics['test_auc']:.3f}  CV AUC={metrics['cv_auc_mean']:.3f}")

    # GSE116256
    print("\n[2] GSE116256 single-cell subset")
    expr_sc, meta_sc = load_gse116256_subset()
    deg, meta_sc = run_real_deg(expr_sc, meta_sc)
    secreted = deg[deg["significant"] & deg["is_secreted"]].copy()
    deg.to_csv(TABLES / "deg_gfi1_low_vs_high_real.csv", index=False)
    secreted.to_csv(TABLES / "secreted_candidates_real.csv", index=False)
    print(f"  Significant secreted: {len(secreted)}")

    print("\n[3] Ligand–receptor scoring")
    interactions = run_real_lr(expr_sc, meta_sc)
    interactions.to_csv(TABLES / "lr_interaction_scores_real.csv", index=False)
    print(f"  Interactions: {len(interactions)}")

    print("\n[4] Figures + report")
    save_figures(deg, secreted, importance, interactions, meta_sc)
    report = write_real_report(deg, secreted, importance, metrics, interactions)
    print(f"  Report → {report}")
    print("\nDone.")


if __name__ == "__main__":
    main()
