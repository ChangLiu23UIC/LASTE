# Manuscript figures

Curated copies of panels from `results/`, kept together so the paper's figure set is one
folder. They are copies, not the originals — regenerate in `results/`, then re-copy.

| File | Source |
|---|---|
| `UMAP_OE.pdf` | `results/oe/figures/` (tuned parameters — see the main README) |
| `Pearson_Correlation_OE.pdf` | `results/oe/figures/` |
| `Heatmap_OE_Proteins_All_Batches.pdf` | `results/oe/figures/` — 51 olfactory markers |
| `OE_Heatmap_All_Proteins_All_Samples.pdf` | `results/oe/figures/Heatmap_All_Proteins_All_Samples.pdf` |
| `GO_Enrichment_OE.pdf` | `results/oe/figures/` |
| `UMAP_SC.pdf`, `Pearson_Correlation_SC.pdf`, `Heatmap_SC_Proteins_All_Samples.pdf`, `GO_Enrichment_SC.pdf` | `results/spinal_cord/figures/` |

The OE all-protein heatmap is renamed on copy (`OE_` prefix) so it stays distinguishable from
its spinal cord counterpart, which has the same name in `results/`.

Every panel here is produced by a script in `src/`, so the whole set regenerates with two
commands (see the main README). Nothing in this folder depends on a notebook.

Two things that used to sit here have moved to `archive/figures/`:
`GO_Enrichment_STRINGdb_bubble.png` (from `archive/notebooks/deg_analysis.ipynb` — still
redrawable by hand from `data/oe_proteomics/`, no network needed) and the
`*_AllGenes_KMeans.pdf` DEG panels (from an archived notebook run against a search that is not
in this repository, so not reproducible at all). Pull either back in if the paper needs it.
