# Central path definitions for the photobleaching proteomics analysis.
#
# Every script imports its inputs and outputs from here, so nothing in the
# repository depends on where it was cloned to or on the current working
# directory:
#
#     python src/oe_analysis.py        # works from any directory
#
# One database search per tissue. `data/oe_proteomics/` is the current search
# of the 17 OE runs, and it is the search behind every committed OE figure,
# table and notebook output. An earlier search of the same 17 raw runs exists
# outside this repository; it quantified them near-identically but carried
# obsolete gene symbols, and nothing here reads it.

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

# ── Inputs ───────────────────────────────────────────────────────────────────
DATA_DIR      = PROJECT_ROOT / "data"
OE_DATA       = DATA_DIR / "oe_proteomics"             # olfactory epithelium (17 TSV)
SC_DATA       = DATA_DIR / "spinal_cord_proteomics"    # spinal cord          (4 TSV)

REFERENCE_DIR = DATA_DIR / "reference"
OLF_XLSX      = REFERENCE_DIR / "olfactory_genes.xlsx"        # ZFIN + Daniocell sheets
OLF_CSV       = REFERENCE_DIR / "olfactory_genes.csv"         # flat list used by the notebook
SC_GENES_XLSX = REFERENCE_DIR / "spinal_cord_gene_lists.xlsx"

# ── Outputs ──────────────────────────────────────────────────────────────────
RESULTS_DIR = PROJECT_ROOT / "results"


class ResultDir:
    """A results/<group>/ folder that splits figures from tables.

    Path arithmetic routes by file extension, so call sites stay unchanged:

        OUTDIR / "UMAP_OE.pdf"      -> results/oe/figures/UMAP_OE.pdf
        OUTDIR / "clusters.csv"     -> results/oe/tables/clusters.csv
    """

    TABLE_SUFFIXES = {".csv", ".tsv", ".xlsx", ".txt"}

    def __init__(self, root):
        self.root = Path(root)

    def __truediv__(self, name):
        sub = "tables" if Path(str(name)).suffix.lower() in self.TABLE_SUFFIXES else "figures"
        target = self.root / sub
        target.mkdir(parents=True, exist_ok=True)
        return target / str(name)

    def mkdir(self, *args, **kwargs):
        for sub in ("figures", "tables"):
            (self.root / sub).mkdir(parents=True, exist_ok=True)

    def __fspath__(self):
        return str(self.root)

    def __str__(self):
        return str(self.root)

    def __repr__(self):
        return f"ResultDir({str(self.root)!r})"


OE_RESULTS  = ResultDir(RESULTS_DIR / "oe")
SC_RESULTS  = ResultDir(RESULTS_DIR / "spinal_cord")
DEG_RESULTS = ResultDir(RESULTS_DIR / "deg")
