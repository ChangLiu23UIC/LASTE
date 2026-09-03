# Photobleaching Proteomics — Zebrafish Olfactory Epithelium and Spinal Cord

Label-free proteomics of photobleached (PB) versus non-photobleached (NPB) cell fractions
in *Danio rerio*, across two tissues:

- **OE** — olfactory epithelium, two batches (`OE1`, `OE2`), 17 runs
- **SC** — spinal cord, single batch, 4 runs

Quantification is NSAF (normalised spectral abundance factor), computed from spectral counts
and protein length and restricted to *Danio rerio* entries. Downstream: sample correlation,
UMAP embedding with Leiden grouping, hierarchically clustered heatmaps, differential
abundance testing, and GO enrichment.

Every figure and table under `results/oe/` and `results/spinal_cord/` reproduces the
corresponding panel of record. Rerunning the scripts regenerates them pixel-for-pixel; the
two places where that needed help are described under "Pinned outputs" below.

## Repository layout

```
data/
  oe_proteomics/             17 TSV — olfactory epithelium
  spinal_cord_proteomics/     4 TSV — spinal cord
  reference/                 olfactory and spinal cord gene reference lists
src/
  paths.py                   all input/output locations (single source of truth)
  oe_analysis.py             OE pipeline  -> results/oe/
  spine_analysis.py          SC pipeline  -> results/spinal_cord/
  umap_tuning.py             OE UMAP with tunable parameters (CLI)
  umap_tuning_large.py       same, oversized fonts/markers for slide-sized panels
results/
  oe/, spinal_cord/          figures/ (PDF + PNG) and tables/ (CSV) per analysis
manuscript/
  figures/                   curated panels selected for the paper
archive/                     notebooks and their figures, patch scripts, screenshots
```

Every script resolves its paths through `src/paths.py`, so all of them can be run from any
working directory. Each folder has its own README with the detail for that layer.

## One OE search

`data/oe_proteomics/` is the only OE dataset here, and every script reads it. A later search of the same 17 raw runs exists outside this
repository. It quantifies the runs near-identically (Spearman rho 0.93–0.98 per run on NSAF
for shared symbols) but annotates them with modern ZFIN symbols where this one still has
`zgc:` / `si:` / legacy names, so marker-list overlap differs sharply: 51 olfactory markers
are detected here against 113 there. That search is deliberately not in the repository —
its filenames are identical to these and its contents are not, and every figure of record
came from the search kept here.

## Sample naming

`{tissue}{batch}_{treatment}_{fraction}_R{replicate}` — e.g. `OE1_PB_neg_R2`.

| Token | Meaning |
|---|---|
| `OE1` / `OE2` | olfactory epithelium batch 1 (`1-27-*` runs) / batch 2 (`5-22-*` runs) |
| `SC` | spinal cord (single batch) |
| `PB` | photobleached |
| `NPB` | non-photobleached control (raw files use the legacy prefix `BY`) |
| `pos` / `neg` | sorted fractions: non-bleached cells / cells the photobleached signal shifts into |
| `R1`–`R3` | replicate |

The raw-filename to sample-label mapping lives in `FILE_MAPPING` at the top of each script.

## Quickstart

```bash
python -m venv .venv && .venv/Scripts/activate      # Windows; use bin/activate on Unix
pip install -r requirements.txt

python src/oe_analysis.py         # -> results/oe/{figures,tables}/
python src/spine_analysis.py      # -> results/spinal_cord/{figures,tables}/
```

The differential-abundance work lives in `archive/notebooks/` and is not part of the pipeline —
see `archive/README.md`.

`src/oe_analysis.py` does not draw a UMAP. The published `UMAP_OE.*` panel comes from
`src/umap_tuning.py`, so a pipeline rerun leaves it intact. It was generated with:

```bash
python src/umap_tuning.py --n_neighbors 3 --min_dist 0.1 --metric correlation \
    --leiden_resolution 4.0 --leiden_neighbors 3 --random_state 12 \
    --repulsion_strength 3.05 --spread 2.1
```

`src/umap_tuning_large.py` takes the same arguments and wrote `UMAP_OE_large.*`.

## Pinned outputs

Two steps in the pipeline are not reproducible on their own, so the published answer is
stored and reused. Delete the file named below to recompute that step from scratch.

- **`results/spinal_cord/tables/UMAP_SC_embedding.csv`** — UMAP on 4 points with
  `n_neighbors=2` is degenerate. Across fresh processes it settles on one of about three
  different layouts, and `random_state`, `np.random.seed`, `NUMBA_NUM_THREADS=1` and
  `PYTHONHASHSEED` all fail to fix it. The coordinates of record are stored here and reused.
  (The 17-sample OE embedding is bit-for-bit reproducible; this affects only SC.)
- **`results/*/tables/GO_Enrichment_*_table.csv`** — g:Profiler is a live service and its
  annotation database moves: the same query returned 207 OE terms in May 2026 and 203 today.
  The tables of record are stored and the figures are drawn from them, so no network is
  needed for a rerun. Re-querying is one deletion away.

With one run per condition, the SC correlation and UMAP panels describe 4 points — they are
descriptive, not statistical, and the SC UMAP in particular carries no cluster structure.

## Analysis parameters of record

- **NSAF matrix** — 1,343 gene symbols x 17 OE samples; 214 x 4 for SC. Log10 transformed,
  z-scored per gene, `Danio rerio` entries only.
- **PB-neg enriched set** — genes with NSAF > 0 in at least 3 of the 6 OE `PB_neg` samples:
  **258 genes**. This is the GO enrichment query, and an experimental definition rather than
  the `olfactory_genes.xlsx` marker list.
- **Olfactory marker overlap** — of the 1,177 unique symbols in `olfactory_genes.xlsx`,
  **51 are detected** in the OE proteome. Spinal cord: 7 of the 1,462 `All Genes_List`
  symbols, 2 of the 256 in the `3-6 dpf_List` subset.
- **GO enrichment** — g:Profiler, organism `drerio`, gene symbols lower-cased, no custom
  background, FDR < 0.05. 207 significant terms for OE (258-gene query), 80 for SC (81-gene
  query). Mixed-case symbols return zero significant terms, hence the lower-casing.
- **UMAP** — proteins detected in at least 3 samples; script defaults `n_neighbors=4`,
  `min_dist=0.1`, `metric=correlation`, `random_state=42`; the published panel uses the tuned
  arguments above. Leiden runs on the 2-D embedding, not on the protein-space distances, so
  the clusters describe the picture rather than being independent evidence for it.
- **Heatmaps** — top 60 most variable proteins, row z-scores, hierarchical clustering on both
  axes; k-means gene clustering (k=50 for OE, k=20 for SC) for the all-protein panels.
- **Differential abundance** (notebook) — Welch t-test on log-transformed NSAF with
  Benjamini-Hochberg correction, k-means grouping of the resulting gene sets.

## Notes

- All figures are Arial at 300 dpi, written as both PDF (vector) and PNG.
- The per-sample `_LABEL_*` dictionaries in the UMAP scripts are hand overrides tuned to this
  embedding. They are correct for the data in `data/oe_proteomics/` and meaningless for any
  other — clear them if the data or the UMAP parameters change.
- The differential-abundance notebooks and every figure they produced are in `archive/`.
  Nothing in `results/` or `manuscript/` depends on them: every committed panel comes from a
  script in `src/`.
- No licence file is included yet — add one before making the repository public.
