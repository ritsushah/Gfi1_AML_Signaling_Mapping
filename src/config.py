"""Central configuration for the Gfi1-AML signaling pipeline."""

from pathlib import Path

# Paths
ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = ROOT / "data" / "raw"
DATA_PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"
TABLES = RESULTS / "tables"
MODELS = RESULTS / "models"

# Ensure directories exist
for d in [DATA_RAW, DATA_PROCESSED, FIGURES, TABLES, MODELS]:
    d.mkdir(parents=True, exist_ok=True)

# Analysis parameters
RANDOM_SEED = 42
GFI1_LOW_QUANTILE = 0.30
GFI1_HIGH_QUANTILE = 0.70
DEG_PADJ_THRESHOLD = 0.05
DEG_LOGFC_THRESHOLD = 0.5

# Starter list of secreted / extracellular factors relevant to niche signaling
# Expand with GO:0005576 or UniProt KW-0964 for production runs
SECRETED_CANDIDATES = [
    "CXCL12", "CXCL8", "CXCL10", "CCL2", "CCL5", "CCL3", "CCL4",
    "IL6", "IL1B", "IL10", "TNF", "IFNG",
    "VEGFA", "VEGFC", "PDGFA", "PDGFB", "FGF2", "HGF",
    "TGFB1", "TGFB2", "BMP4",
    "CSF1", "CSF2", "CSF3", "KITLG", "FLT3LG",
    "ANXA1", "S100A8", "S100A9", "LGALS1", "LGALS3",
    "THBS1", "SPP1", "FN1", "COL1A1", "MMP2", "MMP9", "TIMP1",
    "ANGPT1", "ANGPT2", "CXCL1", "CXCL2",
]

# Demo ligand-receptor pairs (extend with CellPhoneDB / OmniPath)
LR_PAIRS = [
    ("CXCL12", "CXCR4"),
    ("CXCL8", "CXCR1"),
    ("CXCL8", "CXCR2"),
    ("CCL2", "CCR2"),
    ("IL6", "IL6R"),
    ("VEGFA", "FLT1"),
    ("VEGFA", "KDR"),
    ("TGFB1", "TGFBR1"),
    ("TGFB1", "TGFBR2"),
    ("CSF1", "CSF1R"),
    ("KITLG", "KIT"),
    ("PDGFA", "PDGFRA"),
    ("ANGPT1", "TEK"),
    ("SPP1", "CD44"),
]

SENDER_TYPES = ["Leukemic_Blast", "HSC_MPP", "GMP"]
RECEIVER_TYPES = ["MSC", "Endothelial", "Monocyte", "Macrophage"]
