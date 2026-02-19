"""Train model with k-fold cross-validation.
Predicts p-tau (regression) and Thal category (classification).
Uses BaselineRankingModel or GeneEncoderModel.
"""

import argparse
import pickle
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from torch.utils.tensorboard import SummaryWriter
from sklearn.metrics import f1_score, mean_squared_error, mean_absolute_error
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import os
from tqdm import tqdm

from baseline_model import BaselineRankingModel
from gene_encoder_model import GeneEncoderModel, gene_entropy_loss


class EmbeddingDataset(Dataset):
    """Dataset for embeddings with ptau and thal targets."""
    def __init__(self, embeddings, ptau, thal):
        self.embeddings = torch.FloatTensor(embeddings)
        self.ptau = torch.FloatTensor(ptau)
        self.thal = torch.LongTensor(thal)

    def __len__(self):
        return len(self.embeddings)

    def __getitem__(self, idx):
        return self.embeddings[idx], self.ptau[idx], self.thal[idx]


def apply_pooling(embeddings, method='mean'):
    """Apply pooling across axis=1 (n_genes) to get (n_samples, n_dims) shape.

    Args:
        embeddings: numpy array of shape (n_samples, n_genes, n_dims)
        method: 'mean' or 'max'

    Returns:
        Pooled embeddings of shape (n_samples, n_dims)
    """
    if method == 'mean':
        return np.mean(embeddings, axis=1)
    elif method == 'max':
        return np.max(embeddings, axis=1)
    else:
        raise ValueError(f"Unknown method: {method}")



def train_epoch(model, dataloader, criterion_ptau, criterion_thal, optimizer, device, lambda_entropy=0.0):
    """Train for one epoch."""
    model.train()
    total_loss = 0
    ptau_losses = []
    thal_losses = []

    for embeddings, ptau, thal in dataloader:
        embeddings = embeddings.to(device)
        ptau = ptau.to(device)
        thal = thal.to(device)

        optimizer.zero_grad()
        out = model(embeddings)
        ptau_pred, thal_pred = out[0], out[1]
        w = out[2] if len(out) == 3 else None

        loss_thal = criterion_thal(thal_pred, thal)
        loss = loss_thal
        if criterion_ptau is not None:
            loss_ptau = criterion_ptau(ptau_pred, ptau)
            loss = loss + loss_ptau
            ptau_losses.append(loss_ptau.item())
        else:
            ptau_losses.append(0.0)
        if lambda_entropy > 0 and w is not None:
            loss = loss + lambda_entropy * gene_entropy_loss(w)

        loss.backward()
        optimizer.step()

        total_loss += loss.item()
        thal_losses.append(loss_thal.item())
    return total_loss / len(dataloader), np.mean(ptau_losses), np.mean(thal_losses)


def evaluate(model, dataloader, criterion_ptau, criterion_thal, device, idx_to_thal):
    """Evaluate model on validation set."""
    model.eval()
    total_loss = 0
    ptau_preds = []
    ptau_targets = []
    thal_preds = []
    thal_targets = []

    with torch.no_grad():
        for embeddings, ptau, thal in dataloader:
            embeddings = embeddings.to(device)
            ptau = ptau.to(device)
            thal = thal.to(device)

            out = model(embeddings)
            ptau_pred, thal_pred = out[0], out[1]
            thal_class_preds = torch.argmax(thal_pred, dim=1)

            loss_thal = criterion_thal(thal_pred, thal)
            loss = loss_thal
            if criterion_ptau is not None:
                loss_ptau = criterion_ptau(ptau_pred, ptau)
                loss = loss + loss_ptau
            total_loss += loss.item()

            ptau_preds.extend(ptau_pred.cpu().numpy())
            ptau_targets.extend(ptau.cpu().numpy())
            thal_preds.extend(thal_class_preds.cpu().numpy())
            thal_targets.extend(thal.cpu().numpy())

    # Convert back to original labels for F1 score using manual mapping
    thal_preds_labels = [idx_to_thal[i] for i in thal_preds]
    thal_targets_labels = [idx_to_thal[i] for i in thal_targets]

    # Metrics
    mse = mean_squared_error(ptau_targets, ptau_preds)
    mae = mean_absolute_error(ptau_targets, ptau_preds)
    f1 = f1_score(thal_targets_labels, thal_preds_labels, average='weighted')

    return total_loss / len(dataloader), mse, mae, f1


def main(args):
    data_path = os.path.join('data', 'model_data', f'{args.region}_data_for_model_training.pkl')
    print(f"Loading data from {data_path}")

    with open(data_path, 'rb') as f:
        embeddings_folds, patient_info_folds, ptau_folds, thal_folds = pickle.load(f)

    print(f"Loaded {len(embeddings_folds)} folds")

    if args.gene_encoder and args.pooling is not None:
        raise ValueError("--gene_encoder requires pooling=None (no gene pooling); omit --pooling")

    device = torch.device(args.device)
    print(f"Using device: {device}")

    if args.loss_ptau is None or (isinstance(args.loss_ptau, str) and args.loss_ptau.lower() == 'none'):
        criterion_ptau = None
    elif args.gene_encoder:
        criterion_ptau = nn.HuberLoss()
    elif args.loss_ptau == 'mse':
        criterion_ptau = nn.MSELoss()
    elif args.loss_ptau == 'mae':
        criterion_ptau = nn.L1Loss()
    elif args.loss_ptau == 'huber':
        criterion_ptau = nn.HuberLoss()

    all_thal_labels = []
    for thal_fold in thal_folds:
        if isinstance(thal_fold, (list, np.ndarray)):
            all_thal_labels.extend(thal_fold)
        else:
            all_thal_labels.extend(thal_fold.values if hasattr(thal_fold, 'values') else thal_fold)
    all_thal_labels_int = [int(x) for x in all_thal_labels]
    thal_classes = sorted(set(all_thal_labels_int))
    thal_to_idx = {c: i for i, c in enumerate(thal_classes)}
    idx_to_thal = {i: c for c, i in thal_to_idx.items()}
    n_thal_classes = len(thal_classes)
    criterion_thal = nn.CrossEntropyLoss()

    # Prepare patient-level data for all folds first
    n_folds = len(embeddings_folds)
    all_fold_patient_data = []  # List of dicts, one per fold

    # Pipeline: cells -> genes -> patient embedding -> StandardScaler -> PCA -> model
    print(f"\n{'='*60}")
    print("Pipeline: cells -> genes -> patient embedding (StandardScaler+PCA in CV)")
    print(f"{'='*60}")

    for fold_idx in range(n_folds):
        print(f"\nPreparing fold {fold_idx + 1}/{n_folds}...")

        # Step 1: cells -> genes (optional pooling across genes)
        fold_embeddings = embeddings_folds[fold_idx]  # (n_cells, n_genes, n_dims)
        fold_ptau = ptau_folds[fold_idx]
        fold_thal = thal_folds[fold_idx]
        fold_patient_info = patient_info_folds[fold_idx]

        if args.pooling is None:
            pooled_embeddings = fold_embeddings  # Keep (n_cells, n_genes, n_dims)
        else:
            pooled_embeddings = apply_pooling(fold_embeddings, method=args.pooling)  # (n_cells, n_dims)
        print(f"  Pooled embeddings shape: {pooled_embeddings.shape}")

        # Step 2: cells -> patient embedding (mean across cells per patient)
        patient_embeddings = []
        patient_ptau = []
        patient_thal = []

        cell_idx = 0
        for patient_idx, (patient_name, cell_names) in enumerate(fold_patient_info.items()):
            n_cells_patient = len(cell_names)
            patient_cell_embeddings = pooled_embeddings[cell_idx:cell_idx + n_cells_patient]

            # Aggregate patient cells (mean pooling across cells)
            patient_emb = np.mean(patient_cell_embeddings, axis=0)

            patient_embeddings.append(patient_emb)
            patient_ptau.append(fold_ptau[patient_idx])
            patient_thal.append(fold_thal[patient_idx])

            cell_idx += n_cells_patient

        patient_embeddings = np.array(patient_embeddings)
        patient_ptau = np.array(patient_ptau)
        patient_thal = np.array(patient_thal)

        print(f"  Patient embeddings shape: {patient_embeddings.shape}")

        patient_ptau_log = np.log1p(patient_ptau)
        patient_thal_int = [int(t) for t in patient_thal]
        patient_thal_encoded = np.array([thal_to_idx[t] for t in patient_thal_int])
        all_fold_patient_data.append({
            'embeddings': patient_embeddings,
            'ptau': patient_ptau_log,
            'thal': patient_thal_encoded
        })

        print(f"  Number of patients: {len(patient_ptau)}")

    # Leave-one-fold-out cross-validation
    fold_results = []

    for test_fold_idx in range(n_folds):
        print(f"\n{'='*60}")
        print(f"Test Fold {test_fold_idx + 1}/{n_folds}")
        print(f"Training on folds: {[i+1 for i in range(n_folds) if i != test_fold_idx]}")
        print(f"Testing on fold: {test_fold_idx + 1}")
        print(f"{'='*60}")

        train_fold_indices = [i for i in range(n_folds) if i != test_fold_idx]
        train_embeddings = np.concatenate([all_fold_patient_data[i]['embeddings'] for i in train_fold_indices], axis=0)
        train_ptau = np.concatenate([all_fold_patient_data[i]['ptau'] for i in train_fold_indices], axis=0)
        train_thal = np.concatenate([all_fold_patient_data[i]['thal'] for i in train_fold_indices], axis=0)
        test_embeddings = all_fold_patient_data[test_fold_idx]['embeddings']
        test_ptau = all_fold_patient_data[test_fold_idx]['ptau']
        test_thal = all_fold_patient_data[test_fold_idx]['thal']

        n_train, n_test = train_embeddings.shape[0], test_embeddings.shape[0]

        if len(train_embeddings.shape) == 3:
            n_genes, n_dims = train_embeddings.shape[1], train_embeddings.shape[2]

            train_2d = train_embeddings.reshape(n_train, -1)
            test_2d = test_embeddings.reshape(n_test, -1)

            scaler = StandardScaler()
            train_scaled = scaler.fit_transform(train_2d)
            test_scaled = scaler.transform(test_2d)

            train_scaled = train_scaled.reshape(n_train, n_genes, n_dims)
            test_scaled = test_scaled.reshape(n_test, n_genes, n_dims)

            if args.pca_components is not None:
                train_for_pca = train_scaled.reshape(-1, n_dims)
                test_for_pca = test_scaled.reshape(-1, n_dims)
                pca = PCA(n_components=args.pca_components)
                train_pca = pca.fit_transform(train_for_pca)
                test_pca = pca.transform(test_for_pca)
                train_embeddings = train_pca.reshape(n_train, n_genes, args.pca_components)
                test_embeddings = test_pca.reshape(n_test, n_genes, args.pca_components)
                print(f"Applied StandardScaler (fit on train) + PCA 16->{args.pca_components} (fit on scaled train), transform train+test")
            else:
                train_embeddings = train_scaled
                test_embeddings = test_scaled
                print(f"Applied StandardScaler (fit on train), transform train+test (no PCA, keeping {n_dims} dims)")
        else:
            scaler = StandardScaler()
            train_scaled = scaler.fit_transform(train_embeddings)
            test_scaled = scaler.transform(test_embeddings)

            if args.pca_components is not None:
                pca = PCA(n_components=args.pca_components)
                train_embeddings = pca.fit_transform(train_scaled)
                test_embeddings = pca.transform(test_scaled)
                print(f"Applied StandardScaler (fit on train) + PCA 16->{args.pca_components} (fit on scaled train), transform train+test")
            else:
                train_embeddings = train_scaled
                test_embeddings = test_scaled
                print(f"Applied StandardScaler (fit on train), transform train+test (no PCA, keeping 16 dims)")

        print(f"Train: {len(train_embeddings)} patients, Test: {len(test_embeddings)} patients")

        train_dataset = EmbeddingDataset(train_embeddings, train_ptau, train_thal)
        test_dataset = EmbeddingDataset(test_embeddings, test_ptau, test_thal)

        train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True)
        test_loader = DataLoader(test_dataset, batch_size=args.batch_size, shuffle=False)

        if args.gene_encoder:
            if len(train_embeddings.shape) != 3:
                raise ValueError("--gene_encoder requires pooling=None (3D input: n_patients, n_genes, n_dims)")
            n_genes, n_dims = train_embeddings.shape[1], train_embeddings.shape[2]
            model = GeneEncoderModel(
                n_dims=n_dims,
                n_thal_classes=n_thal_classes,
                dropout=args.dropout,
                n_genes=n_genes,
            ).to(device)
        else:
            if len(train_embeddings.shape) == 3:
                input_shape = (train_embeddings.shape[1], train_embeddings.shape[2])
                input_dim = None
            else:
                input_shape = None
                input_dim = train_embeddings.shape[1]

            model = BaselineRankingModel(
                input_dim=input_dim,
                input_shape=input_shape,
                n_thal_classes=n_thal_classes,
                use_attention=args.attention,
                num_heads=args.attention_heads if args.attention else 1,
                dropout=args.dropout
            ).to(device)

        optimizer = optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)

        model_suffix = '_geneenc' if args.gene_encoder else ''
        attention_suffix = '_attn' if args.attention else ''
        dropout_str = f'_drop{args.dropout}'
        wd_str = f'_wd{args.weight_decay}'
        exp_id_str = f"_{args.exp_id}" if getattr(args, 'exp_id', None) else ''
        log_base = getattr(args, 'log_base', 'runs')
        log_dir = os.path.join(log_base, f'{args.region}_{args.pooling}_{args.loss_ptau}{model_suffix}{attention_suffix}{dropout_str}{wd_str}{exp_id_str}_testfold{test_fold_idx+1}')
        if os.path.exists(log_dir):
            try:
                import shutil
                shutil.rmtree(log_dir)
            except (PermissionError, OSError) as e:
                print(f"  Warning: Could not delete {log_dir}: {e}. TensorBoard may have it open. Continuing anyway.")
        writer = SummaryWriter(log_dir)

        best_test_f1 = float('-inf')
        best_epoch = 0
        best_test_mse = None
        best_test_mae = None

        for epoch in range(args.epochs):
            train_loss, train_ptau_loss, train_thal_loss = train_epoch(
                model, train_loader, criterion_ptau, criterion_thal, optimizer, device,
                lambda_entropy=args.lambda_entropy if args.gene_encoder else 0.0
            )
            test_loss, test_mse, test_mae, test_f1 = evaluate(
                model, test_loader, criterion_ptau, criterion_thal, device, idx_to_thal
            )

            writer.add_scalar('Loss/Train', train_loss, epoch)
            writer.add_scalar('Loss/Test', test_loss, epoch)
            writer.add_scalar('Loss/Train_Ptau', train_ptau_loss, epoch)
            writer.add_scalar('Loss/Train_Thal', train_thal_loss, epoch)
            writer.add_scalar('Metrics/Test_MSE', test_mse, epoch)
            writer.add_scalar('Metrics/Test_MAE', test_mae, epoch)
            writer.add_scalar('Metrics/Test_F1', test_f1, epoch)

            if test_f1 > best_test_f1:
                best_test_f1 = test_f1
                best_epoch = epoch
                best_test_mse = test_mse
                best_test_mae = test_mae

            if (epoch + 1) % 10 == 0:
                print(f"Epoch {epoch+1}/{args.epochs} - Train Loss: {train_loss:.4f}, "
                      f"Test Loss: {test_loss:.4f}, Test F1: {test_f1:.4f}, Test MAE: {test_mae:.4f}")

        writer.close()

        print(f"\nTest Fold {test_fold_idx + 1} Best F1: {best_test_f1:.4f} at epoch {best_epoch + 1}")
        print(f"  MSE at best F1: {best_test_mse:.4f}, MAE at best F1: {best_test_mae:.4f}")
        fold_results.append({
            'test_fold': test_fold_idx + 1,
            'best_f1': best_test_f1,
            'best_epoch': best_epoch + 1,
            'best_test_mse': best_test_mse,
            'best_test_mae': best_test_mae
        })

    print(f"\n{'='*60}")
    print("Leave-One-Fold-Out Cross-Validation Summary")
    print(f"{'='*60}")
    for result in fold_results:
        print(f"\033[92mTest Fold {result['test_fold']}: F1={result['best_f1']:.4f}, "
              f"MSE={result['best_test_mse']:.4f}, MAE={result['best_test_mae']:.4f}\033[0m")

    avg_f1 = np.mean([r['best_f1'] for r in fold_results])
    avg_mse = np.mean([r['best_test_mse'] for r in fold_results])
    avg_mae = np.mean([r['best_test_mae'] for r in fold_results])

    print(f"\nAverage F1: {avg_f1:.4f}")
    print(f"Average MSE: {avg_mse:.4f}")
    print(f"Average MAE: {avg_mae:.4f}")

    return avg_f1, avg_mse, avg_mae, fold_results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Train model (baseline or gene encoder)')
    parser.add_argument('--region', type=str, required=True, help='Region name (e.g., EC, ITG, PFC)')
    parser.add_argument('--pooling', type=str, choices=['mean', 'max'], default=None,
                       help='Pooling method: mean or max (default: None)')
    parser.add_argument('--loss_ptau', type=str, default='mse',
                       choices=['mse', 'mae', 'huber', 'none'],
                       help='Loss function for ptau prediction (none to omit ptau loss)')
    parser.add_argument('--batch_size', type=int, default=8, help='Batch size')
    parser.add_argument('--epochs', type=int, default=100, help='Number of epochs')
    parser.add_argument('--lr', type=float, default=0.001, help='Learning rate')
    parser.add_argument('--device', type=str, default='cuda' if torch.cuda.is_available() else 'cpu',
                       help='Device to use (cuda or cpu)')
    parser.add_argument('--attention', action='store_true',
                       help='Use self-attention layer in the model')
    parser.add_argument('--attention_heads', type=int, default=1,
                       help='Number of attention heads (default: 1)')
    parser.add_argument('--dropout', type=float, default=0.2,
                       help='Dropout rate (default: 0.2)')
    parser.add_argument('--weight_decay', type=float, default=1e-1,
                       help='Weight decay for Adam optimizer (default: 1e-1)')
    parser.add_argument('--pca_components', type=int, default=None,
                       help='Number of PCA components; if not set, keep 16 dims (no PCA)')
    parser.add_argument('--gene_encoder', action='store_true',
                       help='Use GeneEncoderModel (gene encoder + attention pooling); requires pooling=None')
    parser.add_argument('--lambda_entropy', type=float, default=1e-3,
                       help='Gene entropy regularization weight (default: 1e-3, gene_encoder only)')
    args = parser.parse_args()

    main(args)
