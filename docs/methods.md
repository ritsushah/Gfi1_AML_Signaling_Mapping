# Methods

## Synthetic data generation (demo mode)

A multi-cell-type expression matrix is generated that embeds the core biological hypothesis:

- A subset of leukemic / progenitor cells are assigned low GFI1 expression.
- In those cells a defined set of secreted factors (CXCL12, IL6, VEGFA, TGFB1, CSF1, SPP1, MMP9, …) is up-regulated.
- Niche cell types (MSC, endothelial, monocytes) express the corresponding receptors.
- A parallel bulk expression + clinical outcome table is created in which a subset of the same secreted genes drive high-risk status and shorter survival.

This design allows the full pipeline to be executed and validated without external downloads while still recovering the planted signal.

## Differential expression

- Cells are stratified into GFI1-Low / Intermediate / High by expression quantiles.
- Analysis is restricted to leukemic and early-progenitor populations.
- Per-gene Wilcoxon rank-sum tests are performed; p-values are adjusted by Benjamini–Hochberg FDR.
- Genes with padj < 0.05 and log2FC > 0.5 are called significantly up-regulated.
- A curated secreted / extracellular gene list is used to prioritize candidates (expandable with GO:0005576 or UniProt keyword searches).

## Machine-learning prioritization

- Expression of candidate secreted genes is used as features.
- An XGBoost classifier predicts binary high-risk status.
- Performance is reported as test-set AUC / average precision and 5-fold cross-validated AUC.
- Feature importance rankings identify the strongest clinical drivers.

## Ligand–receptor communication

- A lightweight mean-expression product score is computed for each (ligand, receptor, sender, receiver) combination.
- Senders are restricted to GFI1-Low leukemic / progenitor cells; receivers are MSC, endothelial and monocyte/macrophage populations.
- For production analyses this step should be replaced by CellPhoneDB, LIANA or NicheNet statistical frameworks.

## Reporting

All tables, figures and a narrative Markdown summary are written automatically to `results/`.
