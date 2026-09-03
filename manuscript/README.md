# Manuscript figures

Curated copies of panels from `results/`, kept together so the paper's figure set is one
folder. They are copies, not the originals — regenerate in `results/`, then re-copy.

| File | Source |
|---|---|
| `UMAP_OE.pdf` | `results/oe/figures/` — tuned parameters, see the main README |
| `Pearson_Correlation_OE.pdf` | `results/oe/figures/` |
| `Heatmap_OE_Proteins_All_Batches.pdf` | `results/oe/figures/` — 51 olfactory markers |
| `OE_Heatmap_All_Proteins_All_Samples.pdf` | `results/oe/figures/Heatmap_All_Proteins_All_Samples.pdf` |
| `GO_Enrichment_OE.pdf` | `results/oe/figures/` |
| `Pearson_Correlation_SC.pdf` | `results/spinal_cord/figures/` |
| `Heatmap_SC_Proteins_All_Samples.pdf` | `results/spinal_cord/figures/` |
| `GO_Enrichment_SC.pdf` | `results/spinal_cord/figures/` |

The OE all-protein heatmap is renamed on copy (`OE_` prefix) so it stays distinguishable from
its spinal cord counterpart, which has the same name in `results/`.

There is no SC UMAP panel here. `src/spine_analysis.py` writes one to
`results/spinal_cord/figures/UMAP_SC.*`, but with one run per condition it describes 4 points
and carries no cluster structure, so it was not carried into the paper's figure set.

Every panel here is produced by a script in `src/`, so the whole set regenerates with the two
commands in the main README. Nothing in this folder depends on a notebook.
