# Gfi1-AML Signaling — Real Public Data Results

**Sources**
- **scRNA-seq:** GSE116256 (van Galen et al., *Cell* 2019) — 3 AML diagnosis samples + 3 healthy BM
- **Bulk + survival:** TCGA-LAML (UCSC Xena HiSeqV2) — 173 patients

## 1. Machine learning on TCGA-LAML

XGBoost trained on expression of secreted-gene candidates to predict high-risk status (death with short survival).

- **Test AUC:** 0.685
- **3-fold CV AUC:** 0.624 ± 0.093

### Top predictive genes

| gene | 0 |
| --- | --- |
| IL10 | 0.0686 |
| PDGFB | 0.06761 |
| CXCL10 | 0.0572 |
| HGF | 0.05308 |
| S100A8 | 0.04265 |
| CXCL2 | 0.03986 |
| CXCL1 | 0.03761 |
| TGFB1 | 0.03222 |
| LGALS1 | 0.03199 |
| ANGPT1 | 0.02903 |
| CCL4 | 0.0287 |
| FLT3LG | 0.02704 |

![ml](figures/ml_feature_importance.png)

## 2. Differential expression in GSE116256

Compared GFI1-Low vs GFI1-High AML cells (n=1181 vs 11).

### Secreted genes ranked by log2FC (Low / High)

| gene | log2FC | padj | mean_GFI1_low | mean_GFI1_high |
| --- | --- | --- | --- | --- |
| LGALS3 | 30.02 | 0.9719 | 1.09 | 1e-09 |
| CXCL8 | 30 | 0.9719 | 1.077 | 1e-09 |
| TIMP1 | 29.88 | 0.9719 | 0.9913 | 1e-09 |
| CCL5 | 29.79 | 1 | 0.9301 | 1e-09 |
| TGFB1 | 28.96 | 1 | 0.5223 | 1e-09 |
| PDGFA | 28.32 | 1 | 0.3351 | 1e-09 |
| TNF | 26.88 | 1 | 0.1238 | 1e-09 |
| FN1 | 26.83 | 1 | 0.1194 | 1e-09 |
| FLT3LG | 26.73 | 1 | 0.1111 | 1e-09 |
| ANGPT2 | 26.05 | 1 | 0.06946 | 1e-09 |
| VEGFA | 26.02 | 1 | 0.06787 | 1e-09 |
| SPP1 | 25.88 | 1 | 0.06181 | 1e-09 |
| MMP2 | 25.4 | 1 | 0.0442 | 1e-09 |
| CXCL10 | 25.06 | 1 | 0.03495 | 1e-09 |
| CXCL1 | 24.89 | 1 | 0.03109 | 1e-09 |

![sec](figures/top_secreted_candidates.png)

## 3. Ligand–receptor map

Scored interactions from AML GFI1-low cells to niche populations: **41** pairs.

| ligand | receptor | sender | receiver | ligand_mean | receptor_mean | interaction_score |
| --- | --- | --- | --- | --- | --- | --- |
| S100A9 | CD44 | AML_GFI1low | Progenitor | 3.959 | 2.968 | 11.75 |
| S100A8 | CD44 | AML_GFI1low | Progenitor | 3.47 | 2.968 | 10.3 |
| S100A9 | CD44 | AML_GFI1low | Monocyte | 3.959 | 2.485 | 9.839 |
| S100A9 | CD44 | AML_GFI1low | Dendritic | 3.959 | 2.402 | 9.507 |
| S100A8 | CD44 | AML_GFI1low | Monocyte | 3.47 | 2.485 | 8.624 |
| S100A8 | CD44 | AML_GFI1low | Dendritic | 3.47 | 2.402 | 8.333 |
| S100A9 | CD44 | AML_GFI1low | T_cell | 3.959 | 2.032 | 8.045 |
| S100A9 | CD44 | AML_GFI1low | Erythroid | 3.959 | 1.966 | 7.782 |
| S100A8 | CD44 | AML_GFI1low | T_cell | 3.47 | 2.032 | 7.052 |
| S100A8 | CD44 | AML_GFI1low | Erythroid | 3.47 | 1.966 | 6.821 |
| TGFB1 | TGFBR1 | AML_GFI1low | Erythroid | 0.5223 | 0.6425 | 0.3355 |
| TGFB1 | TGFBR1 | AML_GFI1low | Dendritic | 0.5223 | 0.5858 | 0.306 |
| TGFB1 | TGFBR2 | AML_GFI1low | T_cell | 0.5223 | 0.5332 | 0.2785 |
| TGFB1 | TGFBR2 | AML_GFI1low | Dendritic | 0.5223 | 0.5137 | 0.2683 |
| TGFB1 | TGFBR1 | AML_GFI1low | T_cell | 0.5223 | 0.4728 | 0.2469 |

![lr](figures/lr_interaction_heatmap.png)

---

*This analysis uses only public real data (TCGA-LAML + GSE116256). Sparse single-cell counts for some ligands limit power; results should be interpreted as hypothesis-generating and validated in larger atlases (e.g. BoneMarrowMap / AML scAtlas).*