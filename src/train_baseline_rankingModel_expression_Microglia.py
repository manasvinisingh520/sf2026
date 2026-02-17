"""Train Microglia baseline on raw expression folds (CERAD + Braak classification).

Loads data from create_raw_expression_data.py with --cell_type Microglia:
    {region}_Microglia_raw_expression_data.pkl

Returns (avg_f1, 0, 0, fold_results) for compatibility with run_overfitting_experiments_microglia_expression.py.
"""

import argparse
import os
import pickle

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from torch.utils.tensorboard import SummaryWriter
from torch.utils.data import DataLoader

from baseline_model_Microglia import BaselineRankingModelMicroglia, CERAD_CLASS_ORDER, BRAAK_TO_IDX
from train_model_Microglia import (
    EmbeddingDatasetMicroglia,
    train_epoch_microglia,
    evaluate_microglia,
)


def main(args):
    data_path = os.path.join("data", "model_data", f"{args.region}_Microglia_raw_expression_data.pkl")
    print(f"Loading Microglia raw-expression data from {data_path}")

    with open(data_path, "rb") as f:
        expression_folds, patient_info_folds, ptau_folds, thal_folds = pickle.load(f)

    n_folds = len(expression_folds)
    print(f"Loaded {n_folds} folds")

    device = torch.device(args.device)
    criterion_cerad = nn.CrossEntropyLoss()
    criterion_braak = nn.CrossEntropyLoss()
    cerad_to_idx = {s: i for i, s in enumerate(CERAD_CLASS_ORDER)}
    n_cerad_classes = len(CERAD_CLASS_ORDER)
    n_braak_classes = len(BRAAK_TO_IDX)

    all_fold_patient_data = []
    for fold_idx in range(n_folds):
        fold_X = np.array(expression_folds[fold_idx], dtype=np.float32)
        fold_ptau = ptau_folds[fold_idx]
        fold_thal = thal_folds[fold_idx]
        fold_cerad = np.array([cerad_to_idx.get(str(t).strip().lower(), 0) for t in fold_thal])
        fold_braak = np.array([BRAAK_TO_IDX.get(int(round(p)), 0) for p in fold_ptau])
        all_fold_patient_data.append({"embeddings": fold_X, "cerad": fold_cerad, "braak": fold_braak})
        print(f"Fold {fold_idx + 1}: X={fold_X.shape}, patients={len(fold_ptau)}")

    fold_results = []
    for test_fold_idx in range(n_folds):
        print(f"\n{'='*60}\nTest Fold {test_fold_idx + 1}/{n_folds}\n{'='*60}")
        train_idx = [i for i in range(n_folds) if i != test_fold_idx]
        train_embeddings = np.concatenate([all_fold_patient_data[i]["embeddings"] for i in train_idx], axis=0)
        train_cerad = np.concatenate([all_fold_patient_data[i]["cerad"] for i in train_idx], axis=0)
        train_braak = np.concatenate([all_fold_patient_data[i]["braak"] for i in train_idx], axis=0)
        test_embeddings = all_fold_patient_data[test_fold_idx]["embeddings"]
        test_cerad = all_fold_patient_data[test_fold_idx]["cerad"]
        test_braak = all_fold_patient_data[test_fold_idx]["braak"]

        scaler = StandardScaler()
        train_embeddings = scaler.fit_transform(train_embeddings)
        test_embeddings = scaler.transform(test_embeddings)
        pca_components = getattr(args, "pca_components", None)
        if pca_components is not None:
            pca = PCA(n_components=pca_components)
            train_embeddings = pca.fit_transform(train_embeddings)
            test_embeddings = pca.transform(test_embeddings)

        train_dataset = EmbeddingDatasetMicroglia(train_embeddings, train_cerad, train_braak)
        test_dataset = EmbeddingDatasetMicroglia(test_embeddings, test_cerad, test_braak)
        train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True)
        test_loader = DataLoader(test_dataset, batch_size=args.batch_size, shuffle=False)

        n_genes = train_embeddings.shape[1]
        model = BaselineRankingModelMicroglia(
            input_dim=n_genes,
            input_shape=None,
            n_cerad_classes=n_cerad_classes,
            n_braak_classes=n_braak_classes,
            use_attention=False,
            num_heads=1,
            dropout=args.dropout,
        ).to(device)
        optimizer = optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)

        log_base = getattr(args, "log_base", "runs")
        dropout_str = f"_drop{args.dropout}"
        wd_str = f"_wd{args.weight_decay}"
        exp_id_str = f"_{args.exp_id}" if getattr(args, "exp_id", None) else ""
        log_dir = os.path.join(log_base, f"{args.region}_Microglia_expr{dropout_str}{wd_str}{exp_id_str}_testfold{test_fold_idx+1}")
        if os.path.exists(log_dir):
            try:
                import shutil
                shutil.rmtree(log_dir)
            except (PermissionError, OSError):
                pass
        writer = SummaryWriter(log_dir)

        best_test_f1 = float("-inf")
        best_epoch = 0
        for epoch in range(args.epochs):
            train_loss = train_epoch_microglia(
                model, train_loader, criterion_cerad, criterion_braak, optimizer, device
            )
            test_loss, test_f1_cerad, test_f1_braak = evaluate_microglia(
                model, test_loader, criterion_cerad, criterion_braak, device
            )
            test_f1 = (test_f1_cerad + test_f1_braak) / 2.0
            writer.add_scalar("Loss/Train", train_loss, epoch)
            writer.add_scalar("Loss/Test", test_loss, epoch)
            writer.add_scalar("Metrics/Test_F1_cerad", test_f1_cerad, epoch)
            writer.add_scalar("Metrics/Test_F1_braak", test_f1_braak, epoch)
            writer.add_scalar("Metrics/Test_F1", test_f1, epoch)
            if test_f1 > best_test_f1:
                best_test_f1 = test_f1
                best_epoch = epoch
            if (epoch + 1) % 10 == 0:
                print(f"Epoch {epoch+1}/{args.epochs} - Train Loss: {train_loss:.4f}, Test F1: {test_f1:.4f} (CERAD: {test_f1_cerad:.4f}, Braak: {test_f1_braak:.4f})")
        writer.close()
        print(f"\nTest Fold {test_fold_idx + 1} Best F1: {best_test_f1:.4f} at epoch {best_epoch + 1}")
        fold_results.append({"test_fold": test_fold_idx + 1, "best_f1": best_test_f1, "best_epoch": best_epoch + 1, "best_test_mse": 0.0, "best_test_mae": 0.0})

    avg_f1 = np.mean([r["best_f1"] for r in fold_results])
    print(f"\n{'='*60}\nLeave-One-Fold-Out Summary (Microglia Expression)\n{'='*60}")
    for r in fold_results:
        print(f"\033[92mTest Fold {r['test_fold']}: F1={r['best_f1']:.4f}\033[0m")
    print(f"\nAverage F1: {avg_f1:.4f}")
    return avg_f1, 0.0, 0.0, fold_results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train Microglia baseline on raw expression (CERAD + Braak)")
    parser.add_argument("--region", type=str, required=True)
    parser.add_argument("--batch_size", type=int, default=32)
    parser.add_argument("--epochs", type=int, default=1000)
    parser.add_argument("--lr", type=float, default=0.001)
    parser.add_argument("--device", type=str, default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--dropout", type=float, default=0.2)
    parser.add_argument("--weight_decay", type=float, default=1e-2)
    parser.add_argument("--pca_components", type=int, default=None)
    parser.add_argument("--log_base", type=str, default="runs")
    parser.add_argument("--exp_id", type=str, default=None)
    args = parser.parse_args()
    main(args)
