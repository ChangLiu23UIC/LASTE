# =============================================================================
#  OE UMAP — parameters can be set here OR passed via CLI (CLI takes priority)
#
#  Usage examples:
#    python src/umap_tuning.py
#    python src/umap_tuning.py --n_neighbors 5 --min_dist 0.1
#    python src/umap_tuning.py --metric cosine --leiden_resolution 0.5
#    python src/umap_tuning.py --n_neighbors 3 --min_dist 0.1 --repulsion_strength 2.0 --spread 1.5
#    python src/umap_tuning.py --n_neighbors 3 --min_dist 0.05 --out results/oe/figures/UMAP_test
#
#  Published panel (results/oe/figures/UMAP_OE.pdf) was produced with:
#    python src/umap_tuning.py --n_neighbors 3 --min_dist 0.1 --metric correlation \
#        --leiden_resolution 4.0 --leiden_neighbors 3 --random_state 12 \
#        --repulsion_strength 3.05 --spread 2.1
#
#  Key repulsion parameters:
#    --repulsion_strength  push points apart (default 1.0, try 1.5-3.0)
#    --spread              overall scale of embedding (default 1.0, try 1.5-2.0)
#    --min_dist            min spacing between points (default 0.1, try 0.2-0.5)
# =============================================================================

import argparse

from paths import RESULTS_DIR

_DEFAULT_OUT = str(RESULTS_DIR / "oe" / "figures" / "UMAP_OE")

_p = argparse.ArgumentParser(description="OE UMAP with tunable parameters")
_p.add_argument("--n_neighbors",        type=int,   default=4,             help="UMAP neighbours (default 4)")
_p.add_argument("--min_dist",           type=float, default=0.1,           help="UMAP min_dist (default 0.1)")
_p.add_argument("--metric",             type=str,   default="correlation",  help="UMAP metric: cosine|euclidean|correlation (default correlation)")
_p.add_argument("--spread",             type=float, default=1.0,           help="UMAP spread — overall embedding scale (default 1.0, try 1.5-2.0)")
_p.add_argument("--repulsion_strength", type=float, default=1.0,           help="UMAP repulsion between points (default 1.0, try 1.5-3.0)")
_p.add_argument("--random_state",       type=int,   default=42,            help="Random seed (default 42)")
_p.add_argument("--leiden_resolution",  type=float, default=1.0,           help="Leiden resolution (default 1.0)")
_p.add_argument("--leiden_neighbors",   type=int,   default=3,             help="Leiden KNN neighbours (default 3)")
_p.add_argument("--label_split_x",      type=float, default=None,          help="X threshold for label side (auto if omitted)")
_p.add_argument("--label_offset",       type=float, default=0.2,           help="Label nudge in data units (default 0.2)")
_p.add_argument("--out",                type=str,   default=_DEFAULT_OUT,  help="Output path without extension (default results/oe/figures/UMAP_OE)")
_args = _p.parse_args()

N_NEIGHBORS        = _args.n_neighbors
MIN_DIST           = _args.min_dist
METRIC             = _args.metric
SPREAD             = _args.spread
REPULSION_STRENGTH = _args.repulsion_strength
RANDOM_STATE       = _args.random_state
LEIDEN_RESOLUTION  = _args.leiden_resolution
LEIDEN_N_NEIGHBORS = _args.leiden_neighbors
LABEL_H_OFFSET     = _args.label_offset
OUT_FILE           = _args.out
# LABEL_SPLIT_X resolved after embedding (auto = midpoint of x range)
_LABEL_SPLIT_X_ARG = _args.label_split_x

import warnings, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.neighbors import NearestNeighbors
from adjustText import adjust_text
import umap.umap_ as umap
import igraph as ig
import leidenalg

# sklearn >= 1.6 removed 'force_all_finite'; patch umap's local binding directly
import inspect as _inspect, sklearn.utils.validation as _skval
if "force_all_finite" not in _inspect.signature(_skval.check_array).parameters:
    _orig_ca = _skval.check_array
    def _ca_compat(*args, **kw):
        v = kw.pop("force_all_finite", None)
        if v is not None:
            kw.setdefault("ensure_all_finite", v)
        return _orig_ca(*args, **kw)
    umap.check_array = _ca_compat  # umap is the umap.umap_ module object here

warnings.filterwarnings("ignore")
sns.set(style="whitegrid", font_scale=1.1)
from paths import OE_DATA

DATA_DIR = OE_DATA
FILE_MAPPING = {
    "1-27-BY-neg-1.tsv": "OE1_NPB_neg_R1",
    "1-27-BY-neg-2.tsv": "OE1_NPB_neg_R2",
    "1-27-BY-pos-1.tsv": "OE1_NPB_pos_R1",
    "1-27-BY-pos-2.tsv": "OE1_NPB_pos_R2",
    "1-27-PB-neg-1.tsv": "OE1_PB_neg_R1",
    "1-27-PB-neg-2.tsv": "OE1_PB_neg_R2",
    "1-27-PB-neg-3.tsv": "OE1_PB_neg_R3",
    "1-27-PB-pos-2.tsv": "OE1_PB_pos_R2",
    "1-27-PB-pos-3.tsv": "OE1_PB_pos_R3",
    "5-22-BY-neg-1.tsv": "OE2_NPB_neg_R1",
    "5-22-BY-pos-1.tsv": "OE2_NPB_pos_R1",
    "5-22-PB-neg-1.tsv": "OE2_PB_neg_R1",
    "5-22-PB-neg-2.tsv": "OE2_PB_neg_R2",
    "5-22-PB-neg-3.tsv": "OE2_PB_neg_R3",
    "5-22-PB-pos-1.tsv": "OE2_PB_pos_R1",
    "5-22-PB-pos-2.tsv": "OE2_PB_pos_R2",
    "5-22-PB-pos-3.tsv": "OE2_PB_pos_R3",
}
SAMPLE_ORDER = list(FILE_MAPPING.values())

SHAPE = {"PB": "o", "NPB": "s"}

# ── Load & normalise ──────────────────────────────────────────────────────────
def compute_nsaf(filepath, organism="Danio rerio", sep="\t"):
    df = pd.read_csv(filepath, sep=sep, low_memory=False)
    df = df[df["Organism"] == organism].copy()
    if df.empty: return pd.Series(dtype=float)
    df["SAF"]  = df["Total Spectral Count"] / df["Length"].replace(0, np.nan)
    df["NSAF"] = df["SAF"] / df["SAF"].sum()
    return df[["Gene","NSAF"]].groupby("Gene").sum()["NSAF"]

print("Loading data...")
raw = {}
for fname, label in FILE_MAPPING.items():
    s = compute_nsaf(DATA_DIR / fname)
    if not s.empty:
        raw[label] = s

nsaf_matrix = pd.DataFrame(raw).fillna(0)          # genes x samples
nsaf_log    = np.log10(nsaf_matrix + 1e-6)
nsaf_log    = nsaf_log.loc[nsaf_log.var(axis=1) > 0]
nsaf_z      = pd.DataFrame(
    StandardScaler().fit_transform(nsaf_log.T).T,
    index=nsaf_log.index, columns=nsaf_log.columns
)
X_samples = nsaf_z.T                                # samples x genes
X_norm    = normalize(X_samples.values, norm="l2")
samples   = list(X_samples.index)
print(f"  {nsaf_z.shape[0]} genes x {len(samples)} samples")

# ── UMAP ──────────────────────────────────────────────────────────────────────
print(f"\nRunning UMAP  (n_neighbors={N_NEIGHBORS}, min_dist={MIN_DIST}, metric={METRIC})...")
reducer   = umap.UMAP(n_neighbors=N_NEIGHBORS, min_dist=MIN_DIST,
                      metric=METRIC, n_components=2, random_state=RANDOM_STATE,
                      spread=SPREAD, repulsion_strength=REPULSION_STRENGTH)
embedding = reducer.fit_transform(X_norm)

# ── Leiden clustering ─────────────────────────────────────────────────────────
print(f"Running Leiden  (n_neighbors={LEIDEN_N_NEIGHBORS}, resolution={LEIDEN_RESOLUTION})...")
n_nbrs = max(2, min(LEIDEN_N_NEIGHBORS, len(samples) - 1))
X_l    = normalize(embedding, norm="l2")
nbrs   = NearestNeighbors(n_neighbors=n_nbrs, metric="cosine").fit(X_l)
knn    = nbrs.kneighbors_graph(mode="connectivity")
knn    = knn.minimum(knn.T)
src, tgt = knn.nonzero()
G      = ig.Graph(edges=list(zip(src.tolist(), tgt.tolist())), directed=False)
part   = leidenalg.find_partition(G, leidenalg.RBConfigurationVertexPartition,
                                  resolution_parameter=LEIDEN_RESOLUTION, seed=RANDOM_STATE)
labels = np.array(part.membership)
n_cl   = len(np.unique(labels))
print(f"  {n_cl} clusters found")

cluster_names = {}
for cl in np.unique(labels):
    idx     = np.where(labels == cl)[0]
    samps   = [samples[i] for i in idx]
    batches = [s.split("_")[0] for s in samps]
    treats  = [s.split("_")[1] for s in samps]
    pols    = [s.split("_")[2] for s in samps]
    name = max(set(batches), key=batches.count)
    if len(set(treats)) == 1:
        name += "_" + treats[0]
    if len(set(pols)) == 1:
        name += "_" + pols[0]
    cluster_names[cl] = name
print(f"  Names: {cluster_names}")

pal      = sns.color_palette("tab10", n_cl)
cl_color = {cl: pal[i] for i, cl in enumerate(np.unique(labels))}

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(11, 9), constrained_layout=True)
ax.set_facecolor("white"); fig.patch.set_facecolor("white")
ax.grid(False)
ax.spines["top"].set_visible(False);  ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(2);   ax.spines["bottom"].set_linewidth(2)

for i, s in enumerate(samples):
    cl        = labels[i]
    color     = cl_color[cl]
    polarity  = s.split("_")[2]
    treatment = s.split("_")[1]
    marker    = SHAPE[treatment]
    if polarity == "pos":
        ax.scatter(embedding[i,0], embedding[i,1],
                   facecolors=color, edgecolors="black",
                   marker=marker, s=160, linewidth=0.8, alpha=0.92, zorder=3)
    else:
        ax.scatter(embedding[i,0], embedding[i,1],
                   facecolors="none", edgecolors=color,
                   marker=marker, s=160, linewidth=2.2, alpha=0.92, zorder=3)

# Auto-compute split if not provided via CLI
LABEL_SPLIT_X = _LABEL_SPLIT_X_ARG if _LABEL_SPLIT_X_ARG is not None \
                else (embedding[:,0].min() + embedding[:,0].max()) / 2

# Per-sample side overrides ("right" = text to the right, "left" = text to the left)
_LABEL_SIDE_OVERRIDES = {
    "OE2_PB_neg_R3": "right",
    "OE1_PB_pos_R3": "right",
    "OE1_PB_pos_R2": "left",
    "OE1_PB_neg_R1": "right",
}

# Extra x nudge on top of the base offset (positive = further right, negative = further left)
_LABEL_X_NUDGE = {
    "OE2_PB_neg_R3":  0.8,
    "OE1_PB_pos_R2": -1.5,
}

# Extra y nudge (positive = further up, negative = further down)
_LABEL_Y_NUDGE = {}

# Samples that should be placed at the same y level as the point (not below)
_LABEL_SAME_Y = {"OE1_PB_pos_R2"}

# Labels — below each point; x<-10 → bottom-right, x>=-10 → bottom-left
texts, tx, ty = [], [], []
for i, s in enumerate(samples):
    x_pt, y_pt = embedding[i,0], embedding[i,1]
    label_color = cl_color[labels[i]]   # match the point's Leiden color
    if s in _LABEL_SIDE_OVERRIDES:
        if _LABEL_SIDE_OVERRIDES[s] == "right":
            x0, ha = x_pt + LABEL_H_OFFSET, "left"
        else:
            x0, ha = x_pt - LABEL_H_OFFSET, "right"
    elif x_pt < -10:
        x0, ha = x_pt + LABEL_H_OFFSET, "left"   # bottom right of point
    else:
        x0, ha = x_pt - LABEL_H_OFFSET, "right"  # bottom left of point
    x0 += _LABEL_X_NUDGE.get(s, 0.0)
    if s in _LABEL_SAME_Y:
        y0, va = y_pt + _LABEL_Y_NUDGE.get(s, 0.0), "center"
    else:
        y0, va = y_pt - LABEL_H_OFFSET + _LABEL_Y_NUDGE.get(s, 0.0), "top"
    texts.append(ax.text(x0, y0, s, fontsize=7.5, ha=ha, va=va,
                         fontfamily="Arial", color=label_color, fontweight="bold"))
    tx.append(x_pt); ty.append(y_pt)

adjust_text(texts, ax=ax,
            x=np.array(tx), y=np.array(ty),
            target_x=np.array(tx), target_y=np.array(ty),
            only_move={"text":"xy","static":"xy","explode":"xy","pull":"xy"},
            force_static=(0.5,0.5), force_explode=(0.2,0.4),
            force_pull=(0.02,0.02), force_text=(0.1,0.2),
            expand=(1.05,1.2), min_arrow_len=4,
            arrowprops=dict(arrowstyle="-", color="#cccccc", lw=0.5))

# Legend — section headers as blank bold entries, then items underneath
def header(txt):
    return plt.Line2D([], [], linestyle="none", label=txt)

sep = plt.Line2D([], [], linestyle="none", label="")

cl_handles = [plt.Line2D([0],[0], marker="o", color="w",
               markerfacecolor=cl_color[cl], markeredgecolor="black",
               markersize=11, label=cluster_names[cl])
              for cl in np.unique(labels)]

fill_handles = [
    plt.Line2D([0],[0], marker="o", color="w", markerfacecolor="#555555",
               markeredgecolor="black", markersize=11, label="pos  (solid)"),
    plt.Line2D([0],[0], marker="o", color="w", markerfacecolor="none",
               markeredgecolor="#555555", markeredgewidth=2.2,
               markersize=11, label="neg  (hollow)"),
]
shape_handles = [
    plt.Line2D([0],[0], marker="o", color="#555555", markersize=11,
               linestyle="None", markeredgecolor="#555555", label="PB  (circle)"),
    plt.Line2D([0],[0], marker="s", color="#555555", markersize=11,
               linestyle="None", markeredgecolor="#555555", label="NPB  (square)"),
]

all_handles = (
    [header("— Leiden Clusters —")] + cl_handles +
    [sep, header("— Polarity —")]   + fill_handles +
    [sep, header("— Treatment —")]  + shape_handles
)
leg = ax.legend(handles=all_handles,
                bbox_to_anchor=(1.02,1), loc="upper left",
                frameon=True, framealpha=1.0, edgecolor="black",
                handlelength=1.5)
for t in leg.get_texts():
    txt = t.get_text()
    if txt.startswith("—"):          # section header
        t.set_fontsize(16); t.set_fontstyle("italic")
        t.set_fontweight("bold"); t.set_color("#444444")
    else:
        t.set_fontsize(16); t.set_fontweight("bold")
    t.set_fontfamily("Arial")

ax.set_xlabel("UMAP-1", fontsize=18, fontweight="bold", fontfamily="Arial")
ax.set_ylabel("UMAP-2", fontsize=18, fontweight="bold", fontfamily="Arial")
ax.set_title("OE Proteomics - UMAP", fontsize=20, fontweight="bold", fontfamily="Arial")
ax.tick_params(axis="both", labelsize=14, width=2, length=6)
for tick in ax.get_xticklabels() + ax.get_yticklabels():
    tick.set_fontweight("bold"); tick.set_fontfamily("Arial")

Path(OUT_FILE).parent.mkdir(parents=True, exist_ok=True)
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT_FILE}.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"\nSaved  {OUT_FILE}.pdf / .png")
print(f"Parameters used:  n_neighbors={N_NEIGHBORS}  min_dist={MIN_DIST}  metric={METRIC}  spread={SPREAD}  repulsion={REPULSION_STRENGTH}  resolution={LEIDEN_RESOLUTION}")
