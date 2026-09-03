# Results

Each analysis writes `figures/` (PDF + PNG, 300 dpi, Arial) and `tables/` (CSV). Everything
here reproduces the panel of record from the data in `data/` — see "Pinned outputs" in the
main README for the two steps that needed to be stored rather than recomputed.

## `oe/` — `src/oe_analysis.py`, `src/umap_tuning*.py` (reads `data/oe_proteomics/`)

| Figure | What it shows |
|---|---|
| `UMAP_OE` | published sample embedding, Leiden clusters (tuned parameters — see the main README) |
| `UMAP_OE_large` | published embedding with oversized fonts and markers for slides |
| `Pearson_Correlation_OE` | 17 x 17 sample correlation, clustered |
| `Heatmap_All_Proteins_All_Samples` | all proteins, 50 k-means gene clusters x 17 samples |
| `Heatmap_All_Proteins_PB_Only` | the same, 11 PB samples only |
| `Heatmap_OE_Proteins_All_Batches` | 51 detected olfactory markers x 17 samples |
| `Heatmap_OE_Proteins_Batch1` / `Batch2` | markers within one batch (48 / 32 genes) |
| `*_PB_Only` | the PB-sample subset of each heatmap (46 / 44 / 28 genes) |
| `GO_Enrichment_OE` | top terms of the 207 significant ones, for the 258-gene PB-neg set |

Tables: `Pearson_Correlation_OE.csv`, `Heatmap_All_Proteins_cluster_assignments.csv`
(gene -> k-means cluster, 1,343 genes), `GO_Enrichment_OE_table.csv` (207 terms, pinned).

## `spinal_cord/` — `src/spine_analysis.py`

The same panel set for the 4 SC runs: `UMAP_SC`, `Pearson_Correlation_SC`,
`Heatmap_All_Proteins_*` (20 k-means clusters), `Heatmap_SC_Proteins_*` (7 detected spinal
cord markers), `Heatmap_SC_Proteins_3to6dpf_*` (2 markers), `GO_Enrichment_SC` (80 terms).

Tables also include `UMAP_SC_embedding.csv` — the pinned 4-point embedding the panel is drawn
from, because this UMAP call is not reproducible run to run (main README explains).

With one run per condition, the SC correlation and UMAP panels describe 4 points — they are
descriptive, not statistical.
