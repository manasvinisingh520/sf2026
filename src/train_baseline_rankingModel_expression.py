"""Train baseline ranking model on patient-level raw expression folds.

This script reuses the model + training utilities from `src/train_model.py` and `src/baseline_model.py`,
but loads raw expression data produced by `src/create_raw_expression_data.py`.

Expected pickle format (per region):
    (expression_folds, patient_info_folds, ptau_folds, thal_folds)

Where:
    expression_folds: list[np.ndarray], each (n_patients_in_fold, n_genes)
    ptau_folds: list[list[float]], each length n_patients_in_fold
    thal_folds: list[list[int/str]], each length n_patients_in_fold
"""

import argparse
import os
import pickle

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.preprocessing import StandardScaler
from torch.utils.tensorboard import SummaryWriter
from torch.utils.data import DataLoader

from baseline_model import BaselineRankingModel
from train_model import EmbeddingDataset, train_epoch, evaluate


def main(args):
    data_path = os.path.join(
        "data",
        "model_data",
        f"{args.region}_raw_expression_data.pkl",
    )
    print(f"Loading raw-expression data from {data_path}")

    with open(data_path, "rb") as f:
        expression_folds, patient_info_folds, ptau_folds, thal_folds = pickle.load(f)

    n_folds = len(expression_folds)
    print(f"Loaded {n_folds} folds")

    if args.attention:
        print("Note: --attention is ignored for raw expression data (1D input); attention only applies to 2D (n_genes, n_dims) inputs.")

    device = torch.device(args.device)
    print(f"Using device: {device}")

    # Loss functions
    if args.loss_ptau == "mse":
        criterion_ptau = nn.MSELoss()
    elif args.loss_ptau == "mae":
        criterion_ptau = nn.L1Loss()
    elif args.loss_ptau == "huber":
        criterion_ptau = nn.HuberLoss()
    else:
        raise ValueError(f"Unknown loss_ptau: {args.loss_ptau}")

    criterion_thal = nn.CrossEntropyLoss()

    # Thal mapping across all folds (manual mapping, consistent with baseline script)
    all_thal_labels = []
    for thal_fold in thal_folds:
        all_thal_labels.extend(thal_fold)
    all_thal_labels_int = [int(x) for x in all_thal_labels]
    thal_classes = sorted(set(all_thal_labels_int))
    thal_to_idx = {c: i for i, c in enumerate(thal_classes)}
    idx_to_thal = {i: c for c, i in thal_to_idx.items()}
    n_thal_classes = len(thal_classes)

    print(f"Thal classes (sorted): {thal_classes}")
    print(f"Thal to index mapping: {thal_to_idx}")

    # Prepare fold data (already patient-level)
    all_fold_patient_data = []
    for fold_idx in range(n_folds):
        fold_X = np.array(expression_folds[fold_idx], dtype=np.float32)  # (n_patients, n_genes)
        fold_ptau = np.array(ptau_folds[fold_idx], dtype=np.float32)
        fold_thal_int = [int(t) for t in thal_folds[fold_idx]]
        fold_thal = np.array([thal_to_idx[t] for t in fold_thal_int], dtype=np.int64)

        # log1p ptau target
        fold_ptau = np.log1p(fold_ptau)

        # keep patient_info_folds for debugging/order checks
        _fold_patient_info = patient_info_folds[fold_idx]
        if len(_fold_patient_info) != len(fold_X):
            print(
                f"Warning: fold {fold_idx}: patient_info has {len(_fold_patient_info)} patients, "
                f"but X has {len(fold_X)} rows. Ensure create_raw_expression_data.py preserves order."
            )

        all_fold_patient_data.append({"embeddings": fold_X, "ptau": fold_ptau, "thal": fold_thal})

        print(f"Fold {fold_idx + 1}: X={fold_X.shape}, patients={len(fold_ptau)}")

    # Leave-one-fold-out CV
    fold_results = []
    for test_fold_idx in range(n_folds):
        print(f"\n{'='*60}")
        print(f"Test Fold {test_fold_idx + 1}/{n_folds}")
        print(f"Training on folds: {[i+1 for i in range(n_folds) if i != test_fold_idx]}")
        print(f"{'='*60}")

        train_fold_indices = [i for i in range(n_folds) if i != test_fold_idx]
        train_embeddings = np.concatenate([all_fold_patient_data[i]["embeddings"] for i in train_fold_indices], axis=0)
        train_ptau = np.concatenate([all_fold_patient_data[i]["ptau"] for i in train_fold_indices], axis=0)
        train_thal = np.concatenate([all_fold_patient_data[i]["thal"] for i in train_fold_indices], axis=0)

        test_embeddings = all_fold_patient_data[test_fold_idx]["embeddings"]
        test_ptau = all_fold_patient_data[test_fold_idx]["ptau"]
        test_thal = all_fold_patient_data[test_fold_idx]["thal"]

        # StandardScaler normalization (fit on train, apply to train + test)
        scaler = StandardScaler()
        train_embeddings = scaler.fit_transform(train_embeddings)
        test_embeddings = scaler.transform(test_embeddings)
        print("Applied StandardScaler normalization (fit on train, applied to train and test)")

        print(f"Train: {len(train_embeddings)} patients, Test: {len(test_embeddings)} patients")

        train_dataset = EmbeddingDataset(train_embeddings, train_ptau, train_thal)
        test_dataset = EmbeddingDataset(test_embeddings, test_ptau, test_thal)

        train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True)
        test_loader = DataLoader(test_dataset, batch_size=args.batch_size, shuffle=False)

        model = BaselineRankingModel(
            input_dim=train_embeddings.shape[1],
            input_shape=None,
            n_thal_classes=n_thal_classes,
            use_attention=args.attention,
            num_heads=args.attention_heads if args.attention else 1,
            dropout=args.dropout
        ).to(device)

        optimizer = optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)

        attention_suffix = "_attn" if args.attention else ""
        dropout_str = f"_drop{args.dropout}"
        wd_str = f"_wd{args.weight_decay}"
        log_dir = os.path.join(
            "runs",
            f"{args.region}_{args.cell_type}_expr_{args.loss_ptau}{attention_suffix}{dropout_str}{wd_str}_testfold{test_fold_idx+1}",
        )
        writer = SummaryWriter(log_dir)

        best_test_f1 = float("-inf")
        best_epoch = 0
        best_test_mse = None
        best_test_mae = None

        for epoch in range(args.epochs):
            train_loss, train_ptau_loss, train_thal_loss = train_epoch(
                model, train_loader, criterion_ptau, criterion_thal, optimizer, device
            )
            test_loss, test_mse, test_mae, test_f1 = evaluate(
                model, test_loader, criterion_ptau, criterion_thal, device, idx_to_thal
            )

            writer.add_scalar("Loss/Train", train_loss, epoch)
            writer.add_scalar("Loss/Test", test_loss, epoch)
            writer.add_scalar("Loss/Train_Ptau", train_ptau_loss, epoch)
            writer.add_scalar("Loss/Train_Thal", train_thal_loss, epoch)
            writer.add_scalar("Metrics/Test_MSE", test_mse, epoch)
            writer.add_scalar("Metrics/Test_MAE", test_mae, epoch)
            writer.add_scalar("Metrics/Test_F1", test_f1, epoch)

            if test_f1 > best_test_f1:
                best_test_f1 = test_f1
                best_epoch = epoch
                best_test_mse = test_mse
                best_test_mae = test_mae

            if (epoch + 1) % 10 == 0:
                print(
                    f"Epoch {epoch+1}/{args.epochs} - Train Loss: {train_loss:.4f}, "
                    f"Test Loss: {test_loss:.4f}, Test F1: {test_f1:.4f}, Test MAE: {test_mae:.4f}"
                )

        writer.close()

        print(f"\nTest Fold {test_fold_idx + 1} Best F1: {best_test_f1:.4f} at epoch {best_epoch + 1}")
        print(f"  MSE at best F1: {best_test_mse:.4f}, MAE at best F1: {best_test_mae:.4f}")

        fold_results.append(
            {
                "test_fold": test_fold_idx + 1,
                "best_f1": best_test_f1,
                "best_epoch": best_epoch + 1,
                "best_test_mse": best_test_mse,
                "best_test_mae": best_test_mae,
            }
        )

    print(f"\n{'='*60}")
    print("Leave-One-Fold-Out Cross-Validation Summary (Expression)")
    print(f"{'='*60}")
    for result in fold_results:
        print(
            f"\033[92mTest Fold {result['test_fold']}: F1={result['best_f1']:.4f}, "
            f"MSE={result['best_test_mse']:.4f}, MAE={result['best_test_mae']:.4f}\033[0m"
        )

    avg_f1 = np.mean([r["best_f1"] for r in fold_results])
    avg_mse = np.mean([r["best_test_mse"] for r in fold_results])
    avg_mae = np.mean([r["best_test_mae"] for r in fold_results])

    print(f"\nAverage F1: {avg_f1:.4f}")
    print(f"Average MSE: {avg_mse:.4f}")
    print(f"Average MAE: {avg_mae:.4f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train baseline ranking model on raw expression folds")
    parser.add_argument("--region", type=str, required=True, help="Region name (e.g., EC, ITG, PFC)")
    parser.add_argument(
        "--cell_type",
        type=str,
        default='Astrocytes',
        choices=["Astrocytes", "Microglia"],
        help="Cell type (Astrocytes or Microglia, default: Astrocytes)",
    )
    parser.add_argument(
        "--loss_ptau",
        type=str,
        default="mse",
        choices=["mse", "mae", "huber"],
        help="Loss function for ptau prediction",
    )
    parser.add_argument("--batch_size", type=int, default=32, help="Batch size")
    parser.add_argument("--epochs", type=int, default=1000, help="Number of epochs")
    parser.add_argument("--lr", type=float, default=0.001, help="Learning rate")
    parser.add_argument(
        "--device",
        type=str,
        default="cuda" if torch.cuda.is_available() else "cpu",
        help="Device to use (cuda or cpu)",
    )
    parser.add_argument("--dropout", type=float, default=0.2, help="Dropout rate (default: 0.2)")
    parser.add_argument("--weight_decay", type=float, default=1e-2, help="Weight decay for Adam optimizer (default: 1e-2)")
    parser.add_argument("--attention", action="store_true", help="Use self-attention layer (ignored for expression data: 1D input)")
    parser.add_argument("--attention_heads", type=int, default=1, help="Number of attention heads when --attention is set (default: 1)")

    main(parser.parse_args())

