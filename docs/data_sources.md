# Public Data Sources

## Single-cell RNA-seq (AML + bone-marrow niche)

| Resource | Accession / Link | Notes |
|----------|------------------|-------|
| van Galen et al. 2019 | [GSE116256](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116256) | Classic Seq-Well dataset (16 AML + 5 healthy) |
| BoneMarrowMap / Zeng 2025 | [Zenodo](https://zenodo.org/records/18528442) | ~1.8 M cells from 21 studies, projected onto a common atlas |
| AML scAtlas (Whittle et al.) | [CELLxGENE](https://cellxgene.cziscience.com/collections/071b706a-7ea7-47a4-bddf-6457725839fc) / Figshare | ~750 k cells, 159 AML + 51 healthy |
| Additional cohorts | Search GEO for “AML scRNA-seq” | Many smaller studies can be integrated |

## Bulk expression + clinical outcome

| Resource | Link | Notes |
|----------|------|-------|
| TCGA-LAML | [GDC](https://portal.gdc.cancer.gov/) / [Xena](https://xenabrowser.net/) | RNA-seq + survival |
| BeatAML2 | [BeatAML](https://biodev.github.io/BeatAML2/) | Large functional + clinical cohort |
| GEO microarray / RNA-seq | GSE12417, GSE37642, GSE71014, … | Useful for external validation |

## Gfi1-focused experimental models

Search GEO / ArrayExpress for mouse *Gfi1* knockout or knockdown hematopoietic datasets. These are valuable for validating direct transcriptional targets.

## Practical download tips

```bash
# Example: GEO series matrix (may still need supplementary count files)
pip install GEOparse
python -c "
import GEOparse
gse = GEOparse.get_GEO(geo='GSE116256', destdir='data/raw')
print(gse.phenotype_data.head())
"

# For large atlases prefer the processed RDS / h5ad objects provided by the authors
# rather than re-processing raw FASTQs.
```

Once you have count matrices and metadata, convert them to the simple `pandas` format used by this pipeline (`cells × genes` expression + metadata DataFrame with at least `cell_type` and a GFI1 expression column) and the existing analysis modules will run unchanged.
