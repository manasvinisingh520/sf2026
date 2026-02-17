"""Scratch: patient counts, or extract best model per (attn, loss, dim) from expression experiments CSV."""

import argparse
import os
import pandas as pd
import anndata as ad
from config import REGION_TO_TAB, DATA_DIR
from utils import read_excel_columns
from read_tensorboard_files import check_overfitting


def extract_best_per_attn_loss_dim(csv_path: str, embedding: bool = None):
    """Read experiments CSV and return the best model (by F1) per group.
    If embedding=True or CSV has 'gene_encoder', group by (gene_encoder, attention, loss_ptau, pca_components).
    Otherwise (expression) group by (attention, loss_ptau, pca_components).
    """
    df = pd.read_csv(csv_path)
    if "avg_f1" not in df.columns:
        raise ValueError(f"CSV must have column 'avg_f1'. Found: {list(df.columns)}")
    df = df.dropna(subset=["avg_f1"])
    use_embedding = embedding if embedding is not None else ("gene_encoder" in df.columns)
    if use_embedding:
        group_cols = ["gene_encoder", "attention", "loss_ptau", "pca_components"]
    else:
        group_cols = ["attention", "loss_ptau", "pca_components"]
    for c in group_cols:
        if c not in df.columns:
            raise ValueError(f"CSV must have column '{c}'. Found: {list(df.columns)}")
    df = df.copy()
    # Normalize no-PCA for expression (embedding uses 2,4,8,16 only)
    if not use_embedding:
        no_pca = (df["pca_components"].isna()) | (df["pca_components"].astype(str).str.lower().isin(("none", "nan", "")))
        df = df.astype({"pca_components": object})
        df.loc[no_pca, "pca_components"] = "None"
    # Best row per group by F1
    idx = df.groupby(group_cols)["avg_f1"].idxmax()
    # INSERT_YOUR_CODE
    overfits = []
    for i in idx:
        row = df.loc[i]
        logdir = tensorboard_dir_for_row(row)
        try:
            overfit = check_overfitting(logdir)  # single run path already includes _testfoldN
        except Exception as e:
            print(f"Warning: Could not check overfitting for {logdir}: {e}")
            overfit = None
        overfits.append(overfit)
    df.loc[idx, "overfit"] = overfits
    best = df.loc[idx].sort_values(group_cols).reset_index(drop=True)
    return best


def tensorboard_dir_for_row(row, log_base=None):
    """Build TensorBoard directory path for a best-model row.
    Example: I:\\sf2026\\Astrocytes_V2_run2\\V2_None_mse_geneenc_drop0.5_wd1.0_geTrue_attnFalse_mse_pca16_all_strong_testfold8
    """
    region = row["region"]
    if log_base is None:
        log_base = os.path.join("I:\\", "sf2026", f"Astrocytes_{region}_run2")
    pooling = row.get("pooling", "mean")
    loss_ptau = row["loss_ptau"]
    gene_encoder = row.get("gene_encoder", False)
    attention = row.get("attention", False)
    dropout = row["dropout"]
    wd = row["weight_decay"]
    exp_id = row.get("exp_id", row.get("run", ""))
    test_fold_idx = row.get("test_fold_idx")
    mname = "_geneenc" if gene_encoder else ("_attn" if attention else "")
    name = f"{region}_{pooling}_{loss_ptau}{mname}_drop{dropout}_wd{wd}_{exp_id}"
    if test_fold_idx is not None:
        name += f"_testfold{test_fold_idx}"
    return os.path.join(log_base, name)


def print_patient_counts():
    """Print the number of unique patients for each Astrocytes region."""
    regions = ["EC", "ITG", "PFC", "V2", "V1"]
    
    print("="*60)
    print("Number of patients per Astrocytes region")
    print("="*60)
    
    for region in regions:
        # Try to load from perCell AnnData file first
        anndata_path = os.path.join('data', 'model_data', f'{region}_AnnData_perCell.h5ad')
        
        if os.path.exists(anndata_path):
            try:
                adata = ad.read_h5ad(anndata_path)
                n_patients = adata.obs['SampleName'].nunique()
                print(f"{region}: {n_patients} patients")
            except Exception as e:
                print(f"{region}: Error loading AnnData - {e}")
        else:
            # Fallback to metadata file
            metadata_path = os.path.join(DATA_DIR, "2025-11-16_Astrocytes_metadata.xlsx")
            tab_index = REGION_TO_TAB.get(region)
            
            if tab_index is not None and os.path.exists(metadata_path):
                try:
                    metadata = read_excel_columns(
                        metadata_path,
                        columns=['cell_annotation', 'SampleName'],
                        sheet_name=tab_index
                    )
                    n_patients = metadata['SampleName'].nunique()
                    print(f"{region}: {n_patients} patients (from metadata)")
                except Exception as e:
                    print(f"{region}: Error loading metadata - {e}")
            else:
                print(f"{region}: No data file found")
    
    print("="*60)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scratch: patient counts or best model per (attn, loss, dim)")
    parser.add_argument("--csv", type=str, default=None, help="Experiments CSV; if set, extract best per combo: Ex. ")
    parser.add_argument("--embedding", action="store_true", help="CSV is from run_overfitting_experiments (embedding); group by ge, attn, loss, pca")
    parser.add_argument("--patients", action="store_true", help="Print patient counts per region (default if no --csv)")
    args = parser.parse_args()

    if args.csv:
        if not os.path.isfile(args.csv):
            print(f"File not found: {args.csv}")
            exit(1)
        best = extract_best_per_attn_loss_dim(args.csv, embedding=args.embedding)
        if args.embedding or "gene_encoder" in best.columns:
            print("Best model per (gene_encoder, attention, loss_ptau, pca_components):")
            print("=" * 80)
            for _, row in best.iterrows():
                print(f"  ge={row['gene_encoder']}, attn={row['attention']}, loss={row['loss_ptau']}, pca={row['pca_components']}  ->  "
                      f"run={row['run']}  F1={row['avg_f1']:.4f}, MSE={row['avg_mse']:.4f}, MAE={row['avg_mae']:.4f}, overfits={row['overfit']}")
        else:
            print("Best model per (attention, loss_ptau, pca_components):")
            print("=" * 80)
            for _, row in best.iterrows():
                print(f"  attn={row['attention']}, loss={row['loss_ptau']}, pca={row['pca_components']}  ->  "
                      f"run={row['run']}  F1={row['avg_f1']:.4f}, MSE={row['avg_mse']:.4f}, MAE={row['avg_mae']:.4f}, overfits={row['overfit']}")
        print("=" * 80)
    else:
        print_patient_counts()
