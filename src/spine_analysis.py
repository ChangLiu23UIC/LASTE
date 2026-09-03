# Spinal Cord Proteomics Analysis
# Same pipeline as oe_analysis.py, adapted for the spinal cord dataset
# 4 samples: SC_NPB_neg, SC_NPB_pos, SC_PB_neg, SC_PB_pos
#
# Run from anywhere:  python src/spine_analysis.py
# Inputs  : data/spinal_cord_proteomics/ + data/reference/*.xlsx
# Outputs : results/spinal_cord/{figures,tables}/

import os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.cluster import KMeans
import umap.umap_ as umap
from adjustText import adjust_text
from gprofiler import GProfiler

warnings.filterwarnings("ignore")
sns.set(style="whitegrid", font_scale=1.1)

from paths import SC_RESULTS, SC_DATA, OLF_XLSX, SC_GENES_XLSX

OUTDIR   = SC_RESULTS
OUTDIR.mkdir(exist_ok=True)

DATA_DIR = SC_DATA

# File -> new sample label  (SC = spinal cord, single batch)
FILE_MAPPING = {
    "BY_neg.tsv": "SC_NPB_neg",
    "BY_pos.tsv": "SC_NPB_pos",
    "PB_neg.tsv": "SC_PB_neg",
    "PB_pos.tsv": "SC_PB_pos",
}
SAMPLE_ORDER = list(FILE_MAPPING.values())

GROUP_PALETTE = {
    "PB-pos":  "red",
    "PB-neg":  "pink",
    "NPB-pos": "royalblue",
    "NPB-neg": "lightblue",
    "Other":   "gray",
}
GROUP_ORDER_LIST = ["PB-pos", "PB-neg", "NPB-pos", "NPB-neg", "Other"]

SHAPE = {"PB": "o", "NPB": "s"}

N_GENE_CLUSTERS = 20   # fewer clusters — only 4 samples


def infer_group(name):
    n = name.lower()
    if "npb" in n and "pos" in n: return "NPB-pos"
    if "npb" in n and "neg" in n: return "NPB-neg"
    if "_pb_" in n and "pos" in n: return "PB-pos"
    if "_pb_" in n and "neg" in n: return "PB-neg"
    # also handle bare PB_neg / PB_pos naming without leading underscore
    if "pb_neg" in n: return "PB-neg"
    if "pb_pos" in n: return "PB-pos"
    return "Other"


def compute_nsaf(filepath, organism="Danio rerio", sep="\t"):
    df = pd.read_csv(filepath, sep=sep, low_memory=False)
    df = df[df["Organism"] == organism].copy()
    if df.empty: return pd.Series(dtype=float)
    df["SAF"]  = df["Total Spectral Count"] / df["Length"].replace(0, np.nan)
    df["NSAF"] = df["SAF"] / df["SAF"].sum()
    return df[["Gene", "NSAF"]].groupby("Gene").sum()["NSAF"]


# =============================================================================
# Load & normalise
# =============================================================================
print("Loading Spine-Data...")
raw = {}
for fname, label in FILE_MAPPING.items():
    s = compute_nsaf(DATA_DIR / fname)
    if not s.empty:
        raw[label] = s

nsaf_matrix = pd.DataFrame(raw).fillna(0)
print(f"  {nsaf_matrix.shape[0]} genes x {nsaf_matrix.shape[1]} samples")
print(f"  Samples: {list(nsaf_matrix.columns)}")

nsaf_log = np.log10(nsaf_matrix + 1e-6)
nsaf_log = nsaf_log.loc[nsaf_log.var(axis=1) > 0]
nsaf_z   = pd.DataFrame(
    StandardScaler().fit_transform(nsaf_log.T).T,
    index=nsaf_log.index, columns=nsaf_log.columns
)

samples  = list(nsaf_z.columns)
X_norm   = normalize(nsaf_z.T.values, norm="l2")
X_norm_index = list(nsaf_z.T.index)   # row order of X_norm / embedding


# =============================================================================
# UMAP  (n_neighbors=2 — minimum for 4 samples)
# =============================================================================
print("\nRunning UMAP...")
# UMAP on 4 points with n_neighbors=2 is degenerate: this call is NOT reproducible.
# Across fresh processes it settles on one of ~3 different layouts, and neither
# random_state, np.random.seed, NUMBA_NUM_THREADS=1 nor PYTHONHASHSEED fixes it. The
# published panel is therefore pinned: the coordinates are written once and reused, so
# the figure is stable for anyone who reruns this script. Delete the CSV to redraw.
_EMB_CSV = OUTDIR / "UMAP_SC_embedding.csv"
if Path(_EMB_CSV).exists():
    _emb_df   = pd.read_csv(_EMB_CSV, index_col=0).loc[list(X_norm_index)]
    embedding = _emb_df[["UMAP1", "UMAP2"]].to_numpy()
    print("  reusing pinned embedding from UMAP_SC_embedding.csv")
else:
    reducer   = umap.UMAP(n_neighbors=2, min_dist=0.1, metric="cosine",
                          n_components=2, random_state=42)
    embedding = reducer.fit_transform(X_norm)
    pd.DataFrame(embedding, index=list(X_norm_index),
                 columns=["UMAP1", "UMAP2"]).to_csv(_EMB_CSV)
    print("  computed a fresh embedding and pinned it to UMAP_SC_embedding.csv")

pal      = sns.color_palette("tab10", len(samples))
pt_color = {s: pal[i] for i, s in enumerate(samples)}

fig, ax = plt.subplots(figsize=(9, 7), constrained_layout=True)
ax.set_facecolor("white"); fig.patch.set_facecolor("white")
ax.grid(False)
ax.spines["top"].set_visible(False);  ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(2);   ax.spines["bottom"].set_linewidth(2)

for i, s in enumerate(samples):
    polarity  = "pos" if "pos" in s.lower() else "neg"
    treatment = "NPB" if "npb" in s.lower() else "PB"
    marker    = SHAPE[treatment]
    color     = pt_color[s]
    if polarity == "pos":
        ax.scatter(embedding[i,0], embedding[i,1],
                   facecolors=color, edgecolors="black",
                   marker=marker, s=200, linewidth=0.8, alpha=0.92, zorder=3)
    else:
        ax.scatter(embedding[i,0], embedding[i,1],
                   facecolors="none", edgecolors=color,
                   marker=marker, s=200, linewidth=2.5, alpha=0.92, zorder=3)

# Label each point directly (only 4 points — no need for adjustText)
for i, s in enumerate(samples):
    ax.text(embedding[i,0], embedding[i,1] + 0.15, s,
            fontsize=16, ha="center", va="bottom",
            fontfamily="Arial", fontweight="bold", color=pt_color[s])

# Legend
sep = plt.Line2D([], [], linestyle="none", label="")
pt_handles = [plt.Line2D([0],[0], marker="o", color="w",
               markerfacecolor=pt_color[s], markeredgecolor="black",
               markersize=11, label=s) for s in samples]
fill_handles = [
    plt.Line2D([0],[0], marker="o", color="w", markerfacecolor="#555",
               markeredgecolor="black", markersize=11, label="pos (solid)"),
    plt.Line2D([0],[0], marker="o", color="w", markerfacecolor="none",
               markeredgecolor="#555", markeredgewidth=2.5, markersize=11, label="neg (hollow)"),
]
shape_handles = [
    plt.Line2D([0],[0], marker="o", color="#555", markersize=11,
               linestyle="None", label="PB (circle)"),
    plt.Line2D([0],[0], marker="s", color="#555", markersize=11,
               linestyle="None", label="NPB (square)"),
]
leg = ax.legend(handles=pt_handles+[sep]+fill_handles+[sep]+shape_handles,
                bbox_to_anchor=(1.02,1), loc="upper left",
                frameon=True, framealpha=1.0, edgecolor="black")
for t in leg.get_texts():
    t.set_fontsize(13); t.set_fontweight("bold"); t.set_fontfamily("Arial")

ax.set_xlabel("UMAP-1", fontsize=16, fontweight="bold", fontfamily="Arial")
ax.set_ylabel("UMAP-2", fontsize=16, fontweight="bold", fontfamily="Arial")
ax.set_title("Spinal Cord Proteomics - UMAP", fontsize=16, fontweight="bold", fontfamily="Arial")
ax.tick_params(axis="both", labelsize=16, width=2, length=6)
for tick in ax.get_xticklabels() + ax.get_yticklabels():
    tick.set_fontweight("bold"); tick.set_fontfamily("Arial")

fig.savefig(OUTDIR / "UMAP_SC.pdf", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(OUTDIR / "UMAP_SC.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("  Saved UMAP_SC")


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


# =============================================================================
# Pearson Correlation
# =============================================================================
print("Computing Pearson correlation...")
sample_corr = nsaf_z.corr(method="pearson")
sample_corr.to_csv(OUTDIR / "Pearson_Correlation_SC.csv")

g_corr = sns.clustermap(
    sample_corr,
    cmap="coolwarm", center=0,
    annot=True, fmt=".2f",
    annot_kws={"size": 14, "fontweight": "bold", "fontfamily": "Arial"},
    row_cluster=True, col_cluster=True,
    figsize=(8, 7),
    linewidths=0.5, linecolor="#cccccc",
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

g_corr.savefig(OUTDIR / "Pearson_Correlation_SC.pdf", dpi=300, bbox_inches="tight")
g_corr.savefig(OUTDIR / "Pearson_Correlation_SC.png", dpi=300, bbox_inches="tight")
plt.close("all")
print("  Saved Pearson_Correlation_SC")


# =============================================================================
# Heatmap helper
# =============================================================================
def plot_heatmap(nsaf_df, title, savename, selected_genes=None, max_genes=60):
    df = nsaf_df.copy()
    if selected_genes is not None:
        overlap = df.index.intersection(selected_genes)
        print(f"  [{savename}] SC gene overlap: {len(overlap)}/{len(selected_genes)}")
        df = df.loc[overlap]
        if df.empty: print(f"  [{savename}] No genes — skipping."); return

    df_log = np.log10(df + 1e-6)
    df_log = df_log.loc[df_log.var(axis=1) > 0]
    if df_log.empty: print(f"  [{savename}] All constant — skipping."); return

    df_z = pd.DataFrame(
        StandardScaler().fit_transform(df_log.T).T,
        index=df_log.index, columns=df_log.columns
    )

    # Most variable if too many
    if len(df_z) > max_genes:
        df_z = df_z.loc[df_z.var(axis=1).nlargest(max_genes).index]

    n_g = len(df_z)
    fw = max(10, len(df_z.columns) * 1.2 + 5)
    fh = max(10, n_g * 0.26 + 4)

    g = sns.clustermap(
        df_z, cmap="magma", vmin=-3, vmax=3,
        row_cluster=True, col_cluster=True,
        figsize=(fw, fh),
        cbar_pos=(1.03, 0.25, 0.025, 0.50), cbar_kws={"label": ""},
        linewidths=0.0, dendrogram_ratio=(0.15, 0.05),
    )
    g.ax_cbar.set_ylabel("Z-scored", fontsize=16, fontweight="bold", fontfamily="Arial", labelpad=12)
    g.ax_cbar.yaxis.set_label_position("right")
    g.ax_cbar.set_title("")
    g.ax_cbar.tick_params(labelsize=16, right=True, left=False, labelright=True, labelleft=False)
    for t in g.ax_cbar.get_yticklabels(): t.set_fontweight("bold"); t.set_fontfamily("Arial")

    tick_fs = max(12, min(16, 400 // n_g))
    g.ax_heatmap.xaxis.tick_top()
    g.ax_heatmap.xaxis.set_label_position("top")
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="left", fontsize=16, fontweight="bold", fontfamily="Arial")
    _title_y = _push_col_dend(g)
    plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=tick_fs, fontweight="bold", fontfamily="Arial")
    g.ax_heatmap.set_xlabel("", fontsize=18)
    g.ax_heatmap.set_ylabel("Genes",   fontsize=18, fontweight="bold", fontfamily="Arial")

    g.fig.suptitle(title, y=_title_y, fontsize=16, fontweight="bold", fontfamily="Arial")
    g.savefig(OUTDIR / f"{savename}.pdf", dpi=300, bbox_inches="tight")
    g.savefig(OUTDIR / f"{savename}.png", dpi=300, bbox_inches="tight")
    plt.close("all")
    print(f"  Saved {savename}  ({n_g} genes x {len(df_z.columns)} samples)")


# =============================================================================
# Spinal cord genes from data/reference/spinal_cord_gene_lists.xlsx
# =============================================================================
print(f"\nLoading spinal cord genes from {SC_GENES_XLSX.name}...")

df_sc_all = pd.read_excel(SC_GENES_XLSX, sheet_name="All Genes_List")
sc_all_raw = set(df_sc_all["ZF_Spinal Cord"].dropna().str.strip().str.lower())

df_sc_dpf = pd.read_excel(SC_GENES_XLSX, sheet_name="3-6 dpf_List")
_dpf_col  = df_sc_dpf.columns[0]   # header row is itself a gene
sc_dpf_raw = set([_dpf_col.strip().lower()]
                 + df_sc_dpf[_dpf_col].dropna().str.strip().str.lower().tolist())

print(f"  All SC genes: {len(sc_all_raw)},  3-6 dpf subset: {len(sc_dpf_raw)}")

nsaf_lower     = {g.lower(): g for g in nsaf_matrix.index}
sc_matched     = [nsaf_lower[g] for g in sc_all_raw if g in nsaf_lower]
sc_matched_dpf = [nsaf_lower[g] for g in sc_dpf_raw  if g in nsaf_lower]
print(f"  Matched in proteomics - All: {len(sc_matched)},  3-6 dpf: {len(sc_matched_dpf)}")


# =============================================================================
# Heatmaps — all samples
# =============================================================================
print("\nGenerating heatmaps...")

# 1. All proteins, all samples (KMeans gene clusters)
print("  Heatmap 1: All proteins, KMeans gene clustering...")
_nsaf_log = np.log10(nsaf_matrix + 1e-6)
_nsaf_log = _nsaf_log.loc[_nsaf_log.var(axis=1) > 0]
_nsaf_z   = pd.DataFrame(
    StandardScaler().fit_transform(_nsaf_log.T).T,
    index=_nsaf_log.index, columns=_nsaf_log.columns
)
_gc   = KMeans(n_clusters=N_GENE_CLUSTERS, random_state=42, n_init=20).fit_predict(_nsaf_z.values)
_rows = []
for cl in range(N_GENE_CLUSTERS):
    mask = _gc == cl
    if mask.sum() == 0: continue
    cdf  = _nsaf_z.loc[mask]
    mr   = cdf.mean(axis=0)
    top  = cdf.var(axis=1).nlargest(min(2, mask.sum())).index.tolist()
    mr.name = " / ".join(top) + f"  (n={mask.sum()})"
    _rows.append(mr)
_cm  = pd.DataFrame(_rows)
_fw  = max(10, len(_cm.columns) * 1.2 + 5)
_fh  = max(10, len(_cm) * 0.26 + 4)
_g1  = sns.clustermap(
    _cm, cmap="magma", vmin=-3, vmax=3,
    row_cluster=True, col_cluster=True,
    figsize=(_fw, _fh), cbar_pos=(1.08, 0.25, 0.025, 0.50), cbar_kws={"label": ""},
    linewidths=0.0, dendrogram_ratio=(0.15, 0.05),
)
_g1.ax_cbar.set_ylabel("Z-scored", fontsize=16, fontweight="bold", fontfamily="Arial", labelpad=12)
_g1.ax_cbar.yaxis.set_label_position("right")
_g1.ax_cbar.set_title("")
_g1.ax_cbar.tick_params(labelsize=16, right=True, left=False, labelright=True, labelleft=False)
for _t in _g1.ax_cbar.get_yticklabels(): _t.set_fontweight("bold"); _t.set_fontfamily("Arial")
_g1.ax_heatmap.xaxis.tick_top()
_g1.ax_heatmap.xaxis.set_label_position("top")
plt.setp(_g1.ax_heatmap.get_xticklabels(), rotation=45, ha="left", fontsize=16, fontweight="bold", fontfamily="Arial")
_title_y_g1 = _push_col_dend(_g1)
plt.setp(_g1.ax_heatmap.get_yticklabels(), fontsize=16, fontweight="bold", fontfamily="Arial")
_g1.ax_heatmap.set_xlabel("", fontsize=18)
_g1.ax_heatmap.set_ylabel("Gene Cluster", fontsize=18, fontweight="bold", fontfamily="Arial")
_g1.fig.suptitle(f"All Proteins - SC - KMeans Gene Clusters (k={N_GENE_CLUSTERS})",
                 y=_title_y_g1, fontsize=16, fontweight="bold", fontfamily="Arial")
_g1.savefig(OUTDIR / "Heatmap_All_Proteins_All_Samples.pdf", dpi=300, bbox_inches="tight")
_g1.savefig(OUTDIR / "Heatmap_All_Proteins_All_Samples.png", dpi=300, bbox_inches="tight")
plt.close("all")
pd.DataFrame({"Gene": _nsaf_z.index, "Cluster": _gc+1}).to_csv(
    OUTDIR / "Heatmap_All_Proteins_cluster_assignments.csv", index=False)
print(f"  Saved Heatmap_All_Proteins_All_Samples  ({len(_cm)} clusters x {len(_cm.columns)} samples)")

# 2. SC-specific (all ZFIN genes), all samples
plot_heatmap(nsaf_matrix, "SC Genes (All ZFIN) - All Samples",
             "Heatmap_SC_Proteins_All_Samples", selected_genes=sc_matched)

# 3. SC-specific (3-6 dpf subset), all samples
if sc_matched_dpf:
    plot_heatmap(nsaf_matrix, "SC Genes (3-6 dpf subset) - All Samples",
                 "Heatmap_SC_Proteins_3to6dpf_All_Samples", selected_genes=sc_matched_dpf)

# PB-only heatmaps
print("\nGenerating PB-only heatmaps...")
pb_cols = [c for c in nsaf_matrix.columns if "PB" in c and "NPB" not in c]
print(f"  PB samples: {pb_cols}")

plot_heatmap(nsaf_matrix[pb_cols], "All Proteins - SC - PB Only",
             "Heatmap_All_Proteins_PB_Only")

plot_heatmap(nsaf_matrix[pb_cols], "SC Genes (All ZFIN) - PB Only",
             "Heatmap_SC_Proteins_PB_Only", selected_genes=sc_matched)

if sc_matched_dpf:
    plot_heatmap(nsaf_matrix[pb_cols], "SC Genes (3-6 dpf subset) - PB Only",
                 "Heatmap_SC_Proteins_3to6dpf_PB_Only", selected_genes=sc_matched_dpf)


# =============================================================================
# GO Enrichment
# =============================================================================
print("\nRunning GO enrichment...")
pb_neg_cols = [c for c in nsaf_matrix.columns if "PB" in c and "neg" in c.lower() and "NPB" not in c]
detected    = (nsaf_matrix[pb_neg_cols] > 0).sum(axis=1)
go_ids      = detected[detected >= 1].index.tolist()
go_genes    = list(set(g.lower() for g in go_ids if g.lower() not in ["nan", ""]))
print(f"  Query: {len(go_genes)} PB-neg detected genes")

# g:Profiler is a live service and its annotation database moves: the same query
# returns a different term list months later. The published panel is therefore pinned
# to the result of record — if the table already exists it is reused, otherwise the API
# is queried and the answer saved. Delete GO_Enrichment_SC_table.csv to re-query.
_GO_CSV = OUTDIR / "GO_Enrichment_SC_table.csv"
if Path(_GO_CSV).exists():
    go_res = pd.read_csv(_GO_CSV)
    print(f"  reusing pinned enrichment from GO_Enrichment_SC_table.csv ({len(go_res)} terms)")
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
    go_res.to_csv(OUTDIR / "GO_Enrichment_SC_table.csv", index=False)
    top = go_res.sort_values("p_value").head(20).copy()
    top["-log10p"]    = -np.log10(top["p_value"].clip(lower=1e-300))
    top["term_short"] = top["name"].apply(lambda x: x if len(x)<=48 else x[:45]+"...")
    src_col = {"GO:BP": "#377EB8", "GO:MF": "#FF7F00", "GO:CC": "#4DAF4A"}
    bar_col = [src_col.get(s, "#777") for s in top["source"]]
    fig, ax = plt.subplots(figsize=(12, 7), constrained_layout=True)
    y = np.arange(len(top))
    ax.barh(y, top["-log10p"], color=bar_col, edgecolor="black", linewidth=0.7, height=0.65)
    ax.set_yticks(y)
    ax.set_yticklabels(top["term_short"].values, fontsize=16, fontweight="bold", fontfamily="Arial")
    ax.invert_yaxis()
    ax.set_xlabel("-log10(adjusted p-value)", fontsize=16, fontweight="bold", fontfamily="Arial")
    ax.set_title("GO Enrichment - SC PB-neg Enriched Proteins", fontsize=16, fontweight="bold", fontfamily="Arial")
    ax.tick_params(axis="x", labelsize=16)
    for t in ax.get_xticklabels(): t.set_fontweight("bold"); t.set_fontfamily("Arial")
    ax.grid(True, axis="x", alpha=0.3, linestyle="--"); ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(2);  ax.spines["bottom"].set_linewidth(2)
    leg_h = [plt.matplotlib.patches.Patch(color=v, label=k) for k,v in src_col.items()]
    ax.legend(handles=leg_h, fontsize=16, frameon=True,
              loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.savefig(OUTDIR / "GO_Enrichment_SC.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "GO_Enrichment_SC.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Saved GO_Enrichment_SC")
else:
    print("  No significant GO terms found.")

print(f"\nDone. All results -> {OUTDIR}/")
