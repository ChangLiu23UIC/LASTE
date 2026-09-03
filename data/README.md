# Data

## `oe_proteomics/` — olfactory epithelium, 17 runs

FragPipe protein-level TSVs, two batches:

| Prefix | Batch | Runs |
|---|---|---|
| `1-27-*` | `OE1` | 9 (2 NPB-neg, 2 NPB-pos, 3 PB-neg, 2 PB-pos) |
| `5-22-*` | `OE2` | 8 (1 NPB-neg, 1 NPB-pos, 3 PB-neg, 3 PB-pos) |

`BY` in a filename is the legacy name for the non-photobleached control (`NPB`).
`1-27-PB-pos-1` and `5-22-BY-neg-2` do not exist — those runs were not acquired.

Everything in the repository reads this folder: `src/oe_analysis.py`, both UMAP scripts, and
`notebooks/deg_analysis.ipynb`. It yields 1,343 gene symbols across the 17 runs.

### The other search, deliberately not included

The same 17 raw runs were searched a second time, later. That search is **not** in this
repository. The comparison, in case it comes up again:

| | this folder | the later search |
|---|---|---|
| Gene symbols across all runs | 1,343 | 1,334 |
| Symbol overlap (case-insensitive) | 643 shared / 699 only here / 691 only there | |
| Legacy placeholders (`zgc:` / `si:` / `wu:`) | 401 | 68 |
| Obsolete names | `actc`, `aldoa`, `anx`, `alpha2(i)` | current ZFIN: `acta1a`, `actn3a`, `abat` |
| Olfactory markers detected (of 1,177) | 51 | 113 |
| Quantification agreement | Spearman rho 0.93–0.98 per run on NSAF for shared symbols | |

The two measure the runs almost identically and *annotate* them very differently. Because the
analysis matches proteins to marker lists by gene symbol, swapping one for the other changes
the marker heatmaps, the GO query and the UMAP layout. The filenames are identical in both, so
if the other search is ever brought back, keep it in its own folder and join on `Protein ID`,
never on `Gene`.

## `spinal_cord_proteomics/` — spinal cord, 4 runs

`BY_neg`, `BY_pos`, `PB_neg`, `PB_pos` — one run per condition, single batch, no replicates.
Read as `SC_NPB_neg`, `SC_NPB_pos`, `SC_PB_neg`, `SC_PB_pos`. 214 gene symbols.

## `reference/`

| File | Used by | Contents |
|---|---|---|
| `olfactory_genes.xlsx` | `src/oe_analysis.py` | sheet `daniocell+zfin`, 1,177 unique olfactory gene symbols |
| `olfactory_genes.csv` | `notebooks/deg_analysis.ipynb` | flat olfactory gene list |
| `spinal_cord_gene_lists.xlsx` | `src/spine_analysis.py` | sheets `All Genes_List` (1,462 symbols) and `3-6 dpf_List` (256), column `ZF_Spinal Cord` |

## TSV columns used

`Gene`, `Length`, `Total Spectral Count`, `Organism`. NSAF per run is `SAF / sum(SAF)`, where
`SAF = Total Spectral Count / Length`, computed over `Organism == "Danio rerio"` rows only and
summed per gene symbol.
