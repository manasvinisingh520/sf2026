"""
End-to-end pipeline to:
1) read an .h5ad
2) locate TF activity scores
3) compute cluster × TF summary (mean/median, delta vs rest, effect size)
4) rank top TFs per cluster
5) assign astrocyte subtype labels from TF signatures (rule-based, editable)
6) output CSVs + diagnostic plots

Requirements:
  pip install scanpy anndata pandas numpy scipy seaborn matplotlib
"""

import re
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import ranksums
import matplotlib.pyplot as plt

# -------------------------
# USER SETTINGS
# -------------------------
H5AD_PATH = "data/GRN/EC_cells_tf_scores.h5ad"

# Column name in adata.obs containing cluster labels (try these in order)
CLUSTER_CANDIDATES = ["cluster", "leiden", "gremlin_cluster", "seurat_clusters"]

# Column name in adata.obs containing pathology stage (optional; for plots)
STAGE_CANDIDATES = ["stage", "pathology_stage", "braak", "braak_stage"]

# TF activity location candidates:
# - in obs columns (e.g. "TF_STAT3" or "STAT3_activity")
# - in obsm matrices (e.g. "X_tf", "X_tf_activity", "tf_activity")
OBSM_CANDIDATES = ["X_tf_activity", "tf_activity", "X_TF", "X_regulonAUC", "X_aucell"]
OBS_REGEX_CANDIDATES = [
    r"^TF[_:]",          # TF_STAT3, TF:STAT3
    r"_tf$",
    r"_activity$",
    r"^regulon_",
    r"^AUC_",
]

# Neighborhood module scoring is NOT needed here since you already have TF activity scores.
# This script assumes TF activity is already computed and stored.

# Label rules: customize as needed for astrocytes
LABEL_RULES = {
    "Reactive_inflammatory": {"STAT3", "RELA", "NFKB1", "NFKB2", "CEBPB", "CEBPD", "IRF1"},
    "Stress_IER": {"FOS", "JUN", "JUNB", "ATF3", "FOSL2", "EGR1"},
    "Interferon": {"IRF7", "IRF9", "STAT1", "STAT2"},
    "Homeostatic_identity": {"SOX9", "NFIA", "NFIB", "NFIX", "HES1", "HEY1"},
    "Proteostasis_UPR": {"ATF4", "XBP1", "DDIT3", "CREB3L2"},
}

TOP_TFS_PER_CLUSTER = 15
MIN_CELLS_PER_CLUSTER = 30  # warn if below

# -------------------------
# HELPERS
# -------------------------
def _find_first_present(candidates, available):
    for c in candidates:
        if c in available:
            return c
    return None

def infer_cluster_col(adata: sc.AnnData) -> str:
    col = _find_first_present(CLUSTER_CANDIDATES, adata.obs.columns)
    if col is None:
        raise ValueError(
            f"Could not find cluster column in adata.obs. "
            f"Tried: {CLUSTER_CANDIDATES}. Available: {list(adata.obs.columns)[:50]}..."
        )
    return col

def infer_stage_col(adata: sc.AnnData) -> str | None:
    return _find_first_present(STAGE_CANDIDATES, adata.obs.columns)

def extract_tf_activity(adata: sc.AnnData) -> pd.DataFrame:
    """
    Returns a DataFrame: index=cells, columns=TFs, values=activity.
    Tries obsm candidates first, then obs column regex matches.
    """
    # 1) Try obsm
    for key in OBSM_CANDIDATES:
        if key in adata.obsm.keys():
            mat = adata.obsm[key]
            # infer TF names
            tf_names = None
            # common places TF names live
            for uns_key in ["tf_names", "TF_names", f"{key}_names", "regulon_names"]:
                if uns_key in adata.uns:
                    tf_names = list(adata.uns[uns_key])
                    break
            if tf_names is None:
                # fallback: if var has TF names (rare for obsm), otherwise make generic names
                tf_names = [f"TF_{i}" for i in range(mat.shape[1])]
            if mat.shape[1] != len(tf_names):
                tf_names = [f"{key}_{i}" for i in range(mat.shape[1])]
            return pd.DataFrame(mat, index=adata.obs_names, columns=tf_names)

    # 2) Try obs regex
    obs_cols = list(adata.obs.columns)
    matched = set()
    for pattern in OBS_REGEX_CANDIDATES:
        rx = re.compile(pattern)
        for c in obs_cols:
            if rx.search(c):
                matched.add(c)

    if len(matched) == 0:
        raise ValueError(
            "Could not locate TF activity in adata.obsm or adata.obs. "
            f"Tried obsm keys: {OBSM_CANDIDATES} and obs regex: {OBS_REGEX_CANDIDATES}."
        )

    # Use matched columns; attempt to clean TF names
    tf_df = adata.obs.loc[:, sorted(matched)].copy()
    # clean names: TF_STAT3 -> STAT3, regulon_STAT3 -> STAT3, etc.
    clean = []
    for c in tf_df.columns:
        cc = re.sub(r"^(TF[_:]|regulon_|AUC_)", "", c)
        cc = re.sub(r"(_activity|_tf)$", "", cc)
        clean.append(cc)
    tf_df.columns = clean
    tf_df.index = adata.obs_names
    return tf_df

def cluster_tf_summary(tf_activity: pd.DataFrame, clusters: pd.Series) -> pd.DataFrame:
    """
    Compute cluster × TF:
      - mean_activity
      - median_activity
      - delta_mean = mean(cluster) - mean(rest)
      - effect_size (Cohen's d) using means and stds
      - p_value via Wilcoxon rank-sum (cluster vs rest)
      - BH-adjusted p
    """
    # Ensure alignment
    clusters = clusters.loc[tf_activity.index]
    clust_vals = clusters.astype(str)

    out_rows = []
    for cl in sorted(clust_vals.unique()):
        idx_in = clust_vals == cl
        idx_out = ~idx_in
        n_in = int(idx_in.sum())
        n_out = int(idx_out.sum())
        if n_in < MIN_CELLS_PER_CLUSTER:
            print(f"[warn] cluster {cl} has only {n_in} cells (threshold {MIN_CELLS_PER_CLUSTER}).")

        X_in = tf_activity.loc[idx_in]
        X_out = tf_activity.loc[idx_out]

        mean_in = X_in.mean(axis=0)
        mean_out = X_out.mean(axis=0)
        med_in = X_in.median(axis=0)
        med_out = X_out.median(axis=0)

        # Cohen's d
        sd_in = X_in.std(axis=0, ddof=1).replace(0, np.nan)
        sd_out = X_out.std(axis=0, ddof=1).replace(0, np.nan)
        pooled = np.sqrt(((n_in - 1) * (sd_in**2) + (n_out - 1) * (sd_out**2)) / (n_in + n_out - 2))
        cohen_d = (mean_in - mean_out) / pooled

        # Wilcoxon rank-sum per TF (vectorized loop; OK for moderate TF counts)
        pvals = []
        for tf in tf_activity.columns:
            try:
                stat, p = ranksums(X_in[tf].values, X_out[tf].values)
            except Exception:
                p = np.nan
            pvals.append(p)
        pvals = np.array(pvals, dtype=float)
        # BH correction
        order = np.argsort(pvals)
        ranked = np.empty_like(order)
        ranked[order] = np.arange(1, len(pvals) + 1)
        padj = pvals * len(pvals) / ranked
        padj = np.minimum.accumulate(padj[np.argsort(order)])[np.argsort(np.argsort(order))]
        padj = np.clip(padj, 0, 1)

        df_cl = pd.DataFrame({
            "cluster": cl,
            "TF": tf_activity.columns,
            "n_in": n_in,
            "n_out": n_out,
            "mean_cluster": mean_in.values,
            "mean_rest": mean_out.values,
            "median_cluster": med_in.values,
            "median_rest": med_out.values,
            "delta_mean": (mean_in - mean_out).values,
            "cohen_d": cohen_d.values,
            "p_value": pvals,
            "p_adj": padj,
        })
        out_rows.append(df_cl)

    return pd.concat(out_rows, ignore_index=True)

def top_tfs_per_cluster(summary: pd.DataFrame, n=TOP_TFS_PER_CLUSTER) -> pd.DataFrame:
    """
    Rank TFs per cluster by delta_mean (or cohen_d).
    """
    return (summary
            .sort_values(["cluster", "delta_mean"], ascending=[True, False])
            .groupby("cluster")
            .head(n)
            .reset_index(drop=True))

def assign_labels_from_top_tfs(top_tfs: pd.DataFrame, rules=LABEL_RULES) -> pd.DataFrame:
    """
    Assign a program label per cluster based on overlap between top TFs and rule sets.
    Returns cluster -> label + score details.
    """
    rows = []
    for cl, g in top_tfs.groupby("cluster"):
        tfset = set(map(str.upper, g["TF"].tolist()))
        scores = {}
        for label, tfs in rules.items():
            scores[label] = len(tfset.intersection(set(map(str.upper, tfs))))
        best = max(scores, key=scores.get)
        rows.append({
            "cluster": cl,
            "assigned_label": best if scores[best] > 0 else "Uncertain",
            "rule_overlap": scores[best],
            "all_rule_overlaps": scores,
            "top_TFs": ", ".join(g["TF"].tolist()[:10]),
        })
    return pd.DataFrame(rows)

def plot_heatmap(summary: pd.DataFrame, top_tfs: pd.DataFrame, out_png="cluster_tf_heatmap.png"):
    """
    Heatmap of delta_mean for top TFs (union across clusters).
    """
    tf_union = sorted(set(top_tfs["TF"]))
    mat = (summary[summary["TF"].isin(tf_union)]
           .pivot(index="cluster", columns="TF", values="delta_mean")
           .fillna(0.0))
    plt.figure(figsize=(min(18, 0.35 * mat.shape[1] + 3), min(10, 0.6 * mat.shape[0] + 2)))
    plt.imshow(mat.values, aspect="auto")
    plt.yticks(range(mat.shape[0]), mat.index)
    plt.xticks(range(mat.shape[1]), mat.columns, rotation=90)
    plt.colorbar(label="delta_mean (cluster - rest)")
    plt.title("TF program specificity per cluster (delta mean activity)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

# -------------------------
# MAIN
# -------------------------
def main():
    adata = sc.read_h5ad(H5AD_PATH)

    cluster_col = infer_cluster_col(adata)
    stage_col = infer_stage_col(adata)

    clusters = adata.obs[cluster_col].copy()
    tf_activity = extract_tf_activity(adata)

    # Optional: restrict to a subset of cells if needed (e.g., EC astrocytes already filtered)
    # tf_activity = tf_activity.loc[adata.obs_names]

    print(f"Using cluster column: {cluster_col}")
    if stage_col:
        print(f"Found stage column: {stage_col}")
    print(f"TF activity shape: {tf_activity.shape} (cells x TFs)")

    summary = cluster_tf_summary(tf_activity, clusters)
    top = top_tfs_per_cluster(summary, n=TOP_TFS_PER_CLUSTER)
    labels = assign_labels_from_top_tfs(top)

    # Save outputs
    """summary.to_csv("cluster_tf_summary.csv", index=False)
    top.to_csv("top_TFs_per_cluster.csv", index=False)
    labels.to_csv("cluster_subtype_assignments.csv", index=False)"""

    # Heatmap diagnostic
    plot_heatmap(summary, top, out_png="cluster_TF_delta_heatmap.png")

    # If you have a cell embedding, you can plot TF activity on UMAP as well:
    # Example: if adata.obsm has "X_umap" or "X_GREmLN_UMAP"
    # (Uncomment and adapt keys)
    #
    # if "X_umap" in adata.obsm:
    #     sc.pl.embedding(adata, basis="umap", color=[cluster_col, stage_col] if stage_col else [cluster_col])
    #
    # Add top TF activities to obs for plotting:
    # for tf in top.groupby("cluster").head(3)["TF"].unique():
    #     if tf in tf_activity.columns:
    #         adata.obs[f"TFact_{tf}"] = tf_activity[tf].values
    # sc.pl.umap(adata, color=[f"TFact_{tf}" for tf in top["TF"].unique()[:8]])

    print("Done. Wrote:")
    print(" - cluster_tf_summary.csv")
    print(" - top_TFs_per_cluster.csv")
    print(" - cluster_subtype_assignments.csv")
    print(" - cluster_TF_delta_heatmap.png")

if __name__ == "__main__":
    main()