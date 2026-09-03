# OE Proteomics Analysis - LASTE methodology
#
# Run from anywhere:  python src/oe_analysis.py
# Inputs  : data/oe_proteomics/ + data/reference/olfactory_genes.xlsx
# Outputs : results/oe/{figures,tables}/

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import umap.umap_ as umap
from scipy.stats import f_oneway
from statsmodels.stats.multitest import multipletests
from adjustText import adjust_text
from gprofiler import GProfiler

warnings.filterwarnings('ignore')

# Style exactly as in umap_cluster.py / heatmap_quantify.py
sns.set(style="whitegrid", font_scale=1.1)

from paths import OE_RESULTS, OE_DATA, OLF_XLSX

OUTDIR = OE_RESULTS
OUTDIR.mkdir(exist_ok=True)

DATA_DIR = OE_DATA

# File → new sample label mapping
FILE_MAPPING = {
    '1-27-BY-neg-1.tsv': 'OE1_NPB_neg_R1',
    '1-27-BY-neg-2.tsv': 'OE1_NPB_neg_R2',
    '1-27-BY-pos-1.tsv': 'OE1_NPB_pos_R1',
    '1-27-BY-pos-2.tsv': 'OE1_NPB_pos_R2',
    '1-27-PB-neg-1.tsv': 'OE1_PB_neg_R1',
    '1-27-PB-neg-2.tsv': 'OE1_PB_neg_R2',
    '1-27-PB-neg-3.tsv': 'OE1_PB_neg_R3',
    '1-27-PB-pos-2.tsv': 'OE1_PB_pos_R2',
    '1-27-PB-pos-3.tsv': 'OE1_PB_pos_R3',
    '5-22-BY-neg-1.tsv': 'OE2_NPB_neg_R1',
    '5-22-BY-pos-1.tsv': 'OE2_NPB_pos_R1',
    '5-22-PB-neg-1.tsv': 'OE2_PB_neg_R1',
    '5-22-PB-neg-2.tsv': 'OE2_PB_neg_R2',
    '5-22-PB-neg-3.tsv': 'OE2_PB_neg_R3',
    '5-22-PB-pos-1.tsv': 'OE2_PB_pos_R1',
    '5-22-PB-pos-2.tsv': 'OE2_PB_pos_R2',
    '5-22-PB-pos-3.tsv': 'OE2_PB_pos_R3',
}


# =============================================================
# NSAF Computation  (from umap_cluster.py / heatmap_quantify.py)
# =============================================================
def compute_nsaf(filepath, organism="Danio rerio", sep="\t"):
    df = pd.read_csv(filepath, sep=sep, low_memory=False)
    df = df[df["Organism"] == organism]
    if df.empty:
        return pd.Series(dtype=float)
    df["SAF"]  = df["Total Spectral Count"] / df["Length"].replace(0, np.nan)
    df["NSAF"] = df["SAF"] / df["SAF"].sum()
    result = df[["Gene", "NSAF"]].groupby("Gene").sum()
    return result["NSAF"]


def build_nsaf_matrix(file_mapping, data_dir):
    """Build genes x samples NSAF matrix with renamed columns."""
    result_dict = {}
    for fname, label in file_mapping.items():
        fpath = data_dir / fname
        series = compute_nsaf(fpath)
        if not series.empty:
            result_dict[label] = series
    # genes x samples (rows=genes, cols=samples)
    return pd.DataFrame(result_dict).fillna(0)


# =============================================================
# Group inference from new sample names  (from heatmap_quantify.py)
# =============================================================
def infer_group(name: str) -> str:
    n = name.lower()
    # check NPB before PB to avoid substring collision
    if "npb" in n and "pos" in n: return "NPB-pos"
    if "npb" in n and "neg" in n: return "NPB-neg"
    if "_pb_" in n and "pos" in n: return "PB-pos"
    if "_pb_" in n and "neg" in n: return "PB-neg"
    return "Other"


GROUP_PALETTE = {
    "PB-pos":  "red",
    "PB-neg":  "pink",
    "NPB-pos": "royalblue",
    "NPB-neg": "lightblue",
    "Other":   "gray",
}
GROUP_ORDER = ["PB-pos", "PB-neg", "NPB-pos", "NPB-neg", "Other"]


# =============================================================
# LOAD & NORMALIZE  (from everything_together.py)
# =============================================================
print("Loading NSAF matrix...")
nsaf_matrix = build_nsaf_matrix(FILE_MAPPING, DATA_DIR)
# nsaf_matrix: genes x samples
print(f"  {nsaf_matrix.shape[0]} genes x {nsaf_matrix.shape[1]} samples")

# Log10 transform + z-score per gene (row-wise)  [from heatmap_quantify.py]
nsaf_log = np.log10(nsaf_matrix + 1e-6)   # genes x samples

# Drop constant genes before z-scoring
row_var = nsaf_log.var(axis=1)
nsaf_log = nsaf_log.loc[row_var > 0]
print(f"  After removing constant genes: {nsaf_log.shape[0]}")

# Z-score per gene (StandardScaler row-wise)  [from everything_together.py]
nsaf_z = pd.DataFrame(
    StandardScaler().fit_transform(nsaf_log.T).T,   # scale per gene across samples
    index=nsaf_log.index,
    columns=nsaf_log.columns
)

# Transpose for UMAP: samples x genes
X_samples = nsaf_z.T   # samples x genes
samples    = list(X_samples.index)

# L2-normalize rows before UMAP  [from everything_together.py]
X_norm = normalize(X_samples.values, norm="l2", axis=1)


# =============================================================
# UMAP
# =============================================================
# The published embedding is produced by src/umap_tuning.py with the tuned
# parameters recorded in the README, so this pipeline does not draw one.



def _push_col_dend(g, gap_frac=0.01):
    """After tick_top(), move the col-dendrogram up so it clears the x-tick labels."""
    try:
        g.fig.canvas.draw()
        renderer = g.fig.canvas.get_renderer()
        t_inv = g.fig.transFigure.inverted()
        tops = [t_inv.transform(lbl.get_window_extent(renderer))[1, 1]
                for lbl in g.ax_heatmap.get_xticklabels() if lbl.get_text()]
        if not tops:
            return 1.02
        p_d = g.ax_col_dendrogram.get_position()
        new_bot = max(tops) + gap_frac
        g.ax_col_dendrogram.set_position(
            [p_d.x0, new_bot, p_d.width, p_d.height])
        return new_bot + p_d.height + 0.05
    except Exception:
        return 1.02


# =============================================================
# PEARSON CORRELATION  (from everything_together.py + pearson_correlation.py)
# =============================================================
print("Computing Pearson correlation...")

# nsaf_z: genes x samples — .corr() correlates COLUMNS = samples → 17x17 sample-sample matrix
sample_corr = nsaf_z.corr(method="pearson")
sample_corr.to_csv(OUTDIR / "Pearson_Correlation_OE.csv")
print(f"  Correlation matrix: {sample_corr.shape} (should be 17x17 samples)")

g_corr = sns.clustermap(
    sample_corr,
    cmap="coolwarm", center=0,
    annot=True, fmt=".2f",
    annot_kws={"size": 10, "fontweight": "bold", "fontfamily": "Arial"},
    row_cluster=True, col_cluster=True,
    figsize=(15, 13),
    linewidths=0.3, linecolor="#cccccc",
    cbar_pos=(1.03, 0.25, 0.025, 0.50),
    cbar_kws={"label": ""},
    dendrogram_ratio=(0.10, 0.05),
)
g_corr.ax_cbar.set_ylabel("Pearson r", fontsize=16, fontweight="bold", fontfamily="Arial", labelpad=12)
g_corr.ax_cbar.yaxis.set_label_position("right")
g_corr.ax_cbar.set_title("")
g_corr.ax_cbar.tick_params(labelsize=16, right=True, left=False, labelright=True, labelleft=False)
for t in g_corr.ax_cbar.get_yticklabels():
    t.set_fontweight("bold"); t.set_fontfamily("Arial")

g_corr.ax_heatmap.xaxis.tick_top()
g_corr.ax_heatmap.xaxis.set_label_position("top")
plt.setp(g_corr.ax_heatmap.get_xticklabels(), rotation=45, ha="left",
         fontsize=16, fontweight="bold", fontfamily="Arial")
_corr_title_y = _push_col_dend(g_corr)
plt.setp(g_corr.ax_heatmap.get_yticklabels(), rotation=0,
         fontsize=16, fontweight="bold", fontfamily="Arial")
g_corr.ax_heatmap.set_xlabel("Sample", fontsize=16, fontweight="bold", fontfamily="Arial")
g_corr.ax_heatmap.set_ylabel("Sample", fontsize=16, fontweight="bold", fontfamily="Arial")
g_corr.fig.suptitle("Sample-Sample Pearson Correlation", y=_corr_title_y, fontsize=16, fontweight="bold", fontfamily="Arial")

g_corr.savefig(OUTDIR / "Pearson_Correlation_OE.pdf", dpi=300, bbox_inches="tight")
g_corr.savefig(OUTDIR / "Pearson_Correlation_OE.png", dpi=300, bbox_inches="tight")
plt.close("all")
print("  Saved Pearson_Correlation_OE")


# =============================================================
# HEATMAP HELPER  (from heatmap_quantify.py plot_heatmap_grouped)
# =============================================================
def plot_heatmap(nsaf_df_genes_x_samples, title, savename,
                 selected_genes=None, keep_common_only=False, min_value=0.0):
    """
    nsaf_df_genes_x_samples: genes x samples raw NSAF matrix.
    Applies log10, z-score, groups by treatment, saves clustermap.
    """
    df = nsaf_df_genes_x_samples.copy()

    # Filter to selected genes first
    if selected_genes is not None:
        overlap = df.index.intersection(selected_genes)
        print(f"  [{savename}] OE gene overlap: {len(overlap)}/{len(selected_genes)}")
        df = df.loc[overlap]
        if df.empty:
            print(f"  [{savename}] No genes found — skipping."); return

    # Keep common genes if requested
    if keep_common_only:
        mask = (df > min_value).all(axis=1)
        df = df.loc[mask]
        if df.empty:
            print(f"  [{savename}] No common genes — skipping."); return

    # Log10 + drop constant rows  [from heatmap_quantify.py]
    df_log = np.log10(df + 1e-6)
    rv = df_log.var(axis=1)
    df_log = df_log.loc[rv > 0]
    if df_log.empty:
        print(f"  [{savename}] All constant after log — skipping."); return

    # Z-score per gene (row-wise)
    df_z = pd.DataFrame(
        StandardScaler().fit_transform(df_log.T).T,
        index=df_log.index, columns=df_log.columns
    )

    n_g  = len(df_z)
    fw   = max(18, len(df_z.columns) * 0.9 + 5)
    fh   = max(14, n_g * 0.26 + 4)

    g = sns.clustermap(
        df_z,
        cmap="magma",
        vmin=-3, vmax=3,
        row_cluster=True,
        col_cluster=True,
        figsize=(fw, fh),
        cbar_pos=(1.03, 0.25, 0.025, 0.50),
        cbar_kws={"label": ""},
        linewidths=0.0,
        dendrogram_ratio=(0.15, 0.05),
    )

    # Colorbar — label on RIGHT side so it doesn't overlap gene names
    g.ax_cbar.set_ylabel("Z-scored", fontsize=16,
                          fontweight="bold", fontfamily="Arial", labelpad=12)
    g.ax_cbar.yaxis.set_label_position("right")
    g.ax_cbar.set_title("")
    g.ax_cbar.set_xlabel("")
    g.ax_cbar.tick_params(labelsize=16, right=True, left=False,
                           labelright=True, labelleft=False)
    for t in g.ax_cbar.get_yticklabels():
        t.set_fontweight("bold"); t.set_fontfamily("Arial")

    tick_fs = max(12, min(16, 400 // n_g))
    g.ax_heatmap.xaxis.tick_top()
    g.ax_heatmap.xaxis.set_label_position("top")
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="left", fontsize=16,
             fontweight="bold", fontfamily="Arial")
    _title_y = _push_col_dend(g)
    plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=tick_fs,
             fontweight="bold", fontfamily="Arial")
    g.ax_heatmap.set_xlabel("", fontsize=18)
    g.ax_heatmap.set_ylabel("Genes",   fontsize=18, fontweight="bold", fontfamily="Arial")

    g.fig.suptitle(title, y=_title_y, fontsize=16, fontweight="bold", fontfamily="Arial")

    g.savefig(OUTDIR / f"{savename}.pdf", dpi=300, bbox_inches="tight")
    g.savefig(OUTDIR / f"{savename}.png", dpi=300, bbox_inches="tight")
    plt.close("all")
    print(f"  Saved {savename}  ({n_g} genes x {len(df_z.columns)} samples)")


# =============================================================
# LOAD OLFACTORY GENES  (from heatmap.py / heatmap_quantify.py)
# =============================================================
print("\nLoading olfactory genes from olf.xlsx...")
olf_df = pd.read_excel(OLF_XLSX, sheet_name="daniocell+zfin")
olf_vals = olf_df.values.ravel()
olf_genes = {str(x).strip().lower() for x in olf_vals if pd.notna(x) and str(x).strip()}
print(f"  Unique olfactory genes: {len(olf_genes)}")

# Match to NSAF matrix (case-insensitive)
nsaf_genes_lower = {g.lower(): g for g in nsaf_matrix.index}
olf_matched = [nsaf_genes_lower[g] for g in olf_genes if g in nsaf_genes_lower]
print(f"  Matched in proteomics: {len(olf_matched)}")


# =============================================================
# HEATMAPS
# =============================================================
print("\nGenerating heatmaps...")

# 1. All proteins, all samples — KMeans gene clustering (50 cluster representatives)
# Each of the 1343 genes is assigned to one of 50 clusters; rows show cluster mean profiles.
print("  Heatmap 1: KMeans gene clustering (k=50 cluster means)...")
from sklearn.cluster import KMeans as _KM

# Build z-scored matrix for clustering (genes x samples)
_nsaf_log = np.log10(nsaf_matrix + 1e-6)
_rv       = _nsaf_log.var(axis=1)
_nsaf_log = _nsaf_log.loc[_rv > 0]
_nsaf_z   = pd.DataFrame(
    StandardScaler().fit_transform(_nsaf_log.T).T,
    index=_nsaf_log.index, columns=_nsaf_log.columns
)

N_GENE_CLUSTERS = 50
_km      = _KM(n_clusters=N_GENE_CLUSTERS, random_state=42, n_init=20)
_gc      = _km.fit_predict(_nsaf_z.values)  # cluster per gene

# Compute mean z-score profile per cluster.
# Row label = top 3 most variable genes in the cluster + gene count.
_rows = []
for cl in range(N_GENE_CLUSTERS):
    mask = _gc == cl
    if mask.sum() == 0: continue
    cluster_df = _nsaf_z.loc[mask]
    mean_row   = cluster_df.mean(axis=0)
    # Pick up to 3 genes with highest variance across samples as representatives
    top_genes  = cluster_df.var(axis=1).nlargest(min(2, mask.sum())).index.tolist()
    rep_name   = " / ".join(top_genes)
    mean_row.name = f"{rep_name}  (n={mask.sum()})"
    _rows.append(mean_row)
_cluster_means = pd.DataFrame(_rows)   # clusters x samples (already z-scored means)

_fw = max(18, len(_cluster_means.columns) * 0.9 + 5)
_fh = max(14, len(_cluster_means) * 0.26 + 4)

_g = sns.clustermap(
    _cluster_means, cmap="magma", vmin=-3, vmax=3,
    row_cluster=True, col_cluster=True,
    figsize=(_fw, _fh),
    cbar_pos=(1.08, 0.25, 0.025, 0.50),
    cbar_kws={"label": ""},
    linewidths=0.0, dendrogram_ratio=(0.15, 0.05),
)
_g.ax_cbar.set_ylabel("Z-scored", fontsize=16, fontweight="bold", fontfamily="Arial", labelpad=12)
_g.ax_cbar.yaxis.set_label_position("right")
_g.ax_cbar.set_title("")
_g.ax_cbar.tick_params(labelsize=16, right=True, left=False, labelright=True, labelleft=False)
for _t in _g.ax_cbar.get_yticklabels(): _t.set_fontweight("bold"); _t.set_fontfamily("Arial")
_g.ax_heatmap.xaxis.tick_top()
_g.ax_heatmap.xaxis.set_label_position("top")
plt.setp(_g.ax_heatmap.get_xticklabels(), rotation=45, ha="left", fontsize=16, fontweight="bold", fontfamily="Arial")
_title_y_g = _push_col_dend(_g)
plt.setp(_g.ax_heatmap.get_yticklabels(), fontsize=16, fontweight="bold", fontfamily="Arial")
_g.ax_heatmap.set_xlabel("", fontsize=18)
_g.ax_heatmap.set_ylabel("Gene Cluster", fontsize=18, fontweight="bold", fontfamily="Arial")
_g.fig.suptitle(f"All Proteins - KMeans Gene Clusters (k={N_GENE_CLUSTERS})", y=_title_y_g, fontsize=16, fontweight="bold", fontfamily="Arial")
_g.savefig(OUTDIR / "Heatmap_All_Proteins_All_Samples.pdf", dpi=300, bbox_inches="tight")
_g.savefig(OUTDIR / "Heatmap_All_Proteins_All_Samples.png", dpi=300, bbox_inches="tight")
plt.close("all")
# Save cluster membership table
pd.DataFrame({"Gene": _nsaf_z.index, "Cluster": _gc + 1}).to_csv(
    OUTDIR / "Heatmap_All_Proteins_cluster_assignments.csv", index=False)
print(f"  Saved Heatmap_All_Proteins_All_Samples  ({N_GENE_CLUSTERS} gene clusters x {len(_cluster_means.columns)} samples)")

# 2. OE-specific, all batches
plot_heatmap(nsaf_matrix, "OE-Specific Proteins - Batch 1 & 2",
             "Heatmap_OE_Proteins_All_Batches",
             selected_genes=olf_matched)

# 3. OE-specific, batch 1 only
oe1_cols = [c for c in nsaf_matrix.columns if c.startswith("OE1")]
plot_heatmap(nsaf_matrix[oe1_cols],
             "OE-Specific Proteins - Batch 1 (OE1)",
             "Heatmap_OE_Proteins_Batch1",
             selected_genes=olf_matched)

# 4. OE-specific, batch 2 only
oe2_cols = [c for c in nsaf_matrix.columns if c.startswith("OE2")]
plot_heatmap(nsaf_matrix[oe2_cols],
             "OE-Specific Proteins - Batch 2 (OE2)",
             "Heatmap_OE_Proteins_Batch2",
             selected_genes=olf_matched)

# =============================================================
# PB-ONLY HEATMAPS (same 4 panels, restricted to PB samples)
# =============================================================
print("\nGenerating PB-only heatmaps...")
pb_cols   = [c for c in nsaf_matrix.columns if "_PB_" in c]
pb_oe1    = [c for c in pb_cols if c.startswith("OE1")]
pb_oe2    = [c for c in pb_cols if c.startswith("OE2")]
print(f"  PB samples: {len(pb_cols)} total  |  OE1: {len(pb_oe1)}  |  OE2: {len(pb_oe2)}")

# 5. All proteins, PB only (KMeans gene clusters)
print("  Heatmap 5: All proteins, PB only (KMeans k=50)...")
_nsaf_log_pb = np.log10(nsaf_matrix[pb_cols] + 1e-6)
_nsaf_log_pb = _nsaf_log_pb.loc[_nsaf_log_pb.var(axis=1) > 0]
_nsaf_z_pb   = pd.DataFrame(
    StandardScaler().fit_transform(_nsaf_log_pb.T).T,
    index=_nsaf_log_pb.index, columns=_nsaf_log_pb.columns
)
from sklearn.cluster import KMeans as _KM2
_gc_pb   = _KM2(n_clusters=N_GENE_CLUSTERS, random_state=42, n_init=20).fit_predict(_nsaf_z_pb.values)
_rows_pb = []
for cl in range(N_GENE_CLUSTERS):
    mask = _gc_pb == cl
    if mask.sum() == 0: continue
    cluster_df = _nsaf_z_pb.loc[mask]
    mean_row   = cluster_df.mean(axis=0)
    top_genes  = cluster_df.var(axis=1).nlargest(min(2, mask.sum())).index.tolist()
    mean_row.name = " / ".join(top_genes) + f"  (n={mask.sum()})"
    _rows_pb.append(mean_row)
_cm_pb = pd.DataFrame(_rows_pb)
_fw_pb  = max(14, len(_cm_pb.columns) * 0.9 + 5)
_fh_pb  = max(14, len(_cm_pb) * 0.26 + 4)
_g5 = sns.clustermap(
    _cm_pb, cmap="magma", vmin=-3, vmax=3,
    row_cluster=True, col_cluster=True,
    figsize=(_fw_pb, _fh_pb),
    cbar_pos=(1.08, 0.25, 0.025, 0.50), cbar_kws={"label": ""},
    linewidths=0.0, dendrogram_ratio=(0.15, 0.05),
)
_g5.ax_cbar.set_ylabel("Z-scored", fontsize=16, fontweight="bold", fontfamily="Arial", labelpad=12)
_g5.ax_cbar.yaxis.set_label_position("right")
_g5.ax_cbar.set_title("")
_g5.ax_cbar.tick_params(labelsize=16, right=True, left=False, labelright=True, labelleft=False)
for _t in _g5.ax_cbar.get_yticklabels(): _t.set_fontweight("bold"); _t.set_fontfamily("Arial")
_g5.ax_heatmap.xaxis.tick_top()
_g5.ax_heatmap.xaxis.set_label_position("top")
plt.setp(_g5.ax_heatmap.get_xticklabels(), rotation=45, ha="left", fontsize=16, fontweight="bold", fontfamily="Arial")
_title_y_g5 = _push_col_dend(_g5)
plt.setp(_g5.ax_heatmap.get_yticklabels(), fontsize=16, fontweight="bold", fontfamily="Arial")
_g5.ax_heatmap.set_xlabel("", fontsize=18)
_g5.ax_heatmap.set_ylabel("Gene Cluster", fontsize=18, fontweight="bold", fontfamily="Arial")
_g5.fig.suptitle(f"All Proteins - PB Only - KMeans Gene Clusters (k={N_GENE_CLUSTERS})", y=_title_y_g5, fontsize=16, fontweight="bold", fontfamily="Arial")
_g5.savefig(OUTDIR / "Heatmap_All_Proteins_PB_Only.pdf", dpi=300, bbox_inches="tight")
_g5.savefig(OUTDIR / "Heatmap_All_Proteins_PB_Only.png", dpi=300, bbox_inches="tight")
plt.close("all")
print(f"  Saved Heatmap_All_Proteins_PB_Only  ({len(_cm_pb)} clusters x {len(_cm_pb.columns)} samples)")

# 6. OE-specific, all batches, PB only
plot_heatmap(nsaf_matrix[pb_cols],
             "OE-Specific Proteins (PB-neg) - Batch 1 & 2 - PB Only",
             "Heatmap_OE_Proteins_All_Batches_PB_Only",
             selected_genes=olf_matched)

# 7. OE-specific, batch 1, PB only
plot_heatmap(nsaf_matrix[pb_oe1],
             "OE-Specific Proteins - Batch 1 (OE1) - PB Only",
             "Heatmap_OE_Proteins_Batch1_PB_Only",
             selected_genes=olf_matched)

# 8. OE-specific, batch 2, PB only
plot_heatmap(nsaf_matrix[pb_oe2],
             "OE-Specific Proteins - Batch 2 (OE2) - PB Only",
             "Heatmap_OE_Proteins_Batch2_PB_Only",
             selected_genes=olf_matched)


# =============================================================
# GO ENRICHMENT  (from enrichment_analysis.py pattern + gprofiler)
# =============================================================
print("\nRunning GO enrichment...")

# PB-neg enriched genes: detected in >= 3 PB_neg samples
pb_neg_cols = [c for c in nsaf_matrix.columns if "_PB_neg" in c]
detected    = (nsaf_matrix[pb_neg_cols] > 0).sum(axis=1)
enriched_ids = detected[detected >= 3].index.tolist()
go_genes     = [g.lower() for g in enriched_ids if g.lower() not in ["nan", ""]]
go_genes     = list(set(go_genes))
print(f"  Query: {len(go_genes)} PB-neg enriched genes")

# g:Profiler is a live service and its annotation database moves: the same query
# returns a different term list months later. The published panel is therefore pinned
# to the result of record — if the table already exists it is reused, otherwise the API
# is queried and the answer saved. Delete GO_Enrichment_OE_table.csv to re-query.
_GO_CSV = OUTDIR / "GO_Enrichment_OE_table.csv"
if Path(_GO_CSV).exists():
    go_res = pd.read_csv(_GO_CSV)
    print(f"  reusing pinned enrichment from GO_Enrichment_OE_table.csv ({len(go_res)} terms)")
else:
    gp = GProfiler(return_dataframe=True)
    try:
        go_res = gp.profile(
            organism="drerio", query=go_genes,
            sources=["GO:BP", "GO:MF", "GO:CC"],
            significance_threshold_method="fdr",
            user_threshold=0.05,
        )
        print(f"  Significant GO terms: {len(go_res)}")
    except Exception as e:
        print(f"  g:Profiler error: {e}"); go_res = pd.DataFrame()

if not go_res.empty:
    go_res.to_csv(OUTDIR / "GO_Enrichment_OE_table.csv", index=False)
    top = go_res.sort_values("p_value").head(20).copy()
    top["-log10p"]    = -np.log10(top["p_value"].clip(lower=1e-300))
    top["term_short"] = top["name"].apply(lambda x: x if len(x) <= 50 else x[:47] + "...")

    src_col = {"GO:BP": "#377EB8", "GO:MF": "#FF7F00", "GO:CC": "#4DAF4A"}
    bar_col = [src_col.get(s, "gray") for s in top["source"]]

    fig, ax = plt.subplots(figsize=(11, 7))
    y = np.arange(len(top))
    ax.barh(y, top["-log10p"], color=bar_col, edgecolor="black", linewidth=0.7, height=0.65)
    ax.set_yticks(y)
    ax.set_yticklabels(top["term_short"].values, fontsize=16,
                       fontweight="bold", fontfamily="Arial")
    ax.invert_yaxis()
    ax.set_xlabel("-log10(adjusted p-value)", fontsize=16,
                  fontweight="bold", fontfamily="Arial")
    ax.set_title("GO Enrichment - OE-Enriched Proteins (PB-neg)", fontsize=16,
                 fontweight="bold", fontfamily="Arial")
    ax.tick_params(axis="x", labelsize=16)
    for t in ax.get_xticklabels():
        t.set_fontweight("bold"); t.set_fontfamily("Arial")
    ax.grid(True, axis="x", alpha=0.3, linestyle="--")
    ax.set_axisbelow(True)
    leg_h = [plt.matplotlib.patches.Patch(color=v, label=k) for k, v in src_col.items()]
    ax.legend(handles=leg_h, fontsize=16, frameon=True,
              loc="center left", bbox_to_anchor=(1.02, 0.5))
    plt.tight_layout()
    fig.savefig(OUTDIR / "GO_Enrichment_OE.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "GO_Enrichment_OE.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Saved GO_Enrichment_OE")
else:
    print("  No significant GO terms.")

print(f"\nDone. All results -> {OUTDIR}/")
