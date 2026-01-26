"""Train baseline ranking model with k-fold cross-validation.
Predicts p-tau (regression) and Thal category (classification).
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
import os
from tqdm import tqdm


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


class BaselineRankingModel(nn.Module):
    """Baseline model with two heads: ptau regression and thal classification."""
    def __init__(self, input_dim=16, hidden_dim1=64, hidden_dim2=32, n_thal_classes=None):
        super(BaselineRankingModel, self).__init__()
        
        # Shared layers
        self.shared = nn.Sequential(
            nn.Linear(input_dim, hidden_dim1),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(hidden_dim1, hidden_dim2),
            nn.ReLU(),
            nn.Dropout(0.2)
        )
        
        # Head 1: P-tau regression (ReLU output for log1p(p-tau))
        self.ptau_head = nn.Sequential(
            nn.Linear(hidden_dim2, 1),
            nn.ReLU()  # Ensures non-negative output for log1p
        )
        
        # Head 2: Thal classification (softmax)
        if n_thal_classes is None:
            raise ValueError("n_thal_classes must be specified")
        self.thal_head = nn.Sequential(
            nn.Linear(hidden_dim2, n_thal_classes),
            nn.Softmax(dim=1)
        )
    
    def forward(self, x):
        shared_features = self.shared(x)
        ptau_pred = self.ptau_head(shared_features).squeeze()
        thal_pred = self.thal_head(shared_features)  # (batch_size, n_thal_classes) with softmax
        return ptau_pred, thal_pred


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
        raise ValueError(f"Unknown pooling method: {method}. Use 'mean' or 'max'.")


def train_epoch(model, dataloader, criterion_ptau, criterion_thal, optimizer, device):
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
        ptau_pred, thal_pred = model(embeddings)
        
        # Compute losses
        loss_ptau = criterion_ptau(ptau_pred, ptau)
        # Cross-entropy loss for classification (takes softmax probabilities and class labels)
        loss_thal = criterion_thal(thal_pred, thal)
        loss = loss_ptau + loss_thal
        
        loss.backward()
        optimizer.step()
        
        total_loss += loss.item()
        ptau_losses.append(loss_ptau.item())
        thal_losses.append(loss_thal.item())
        print(f"PTAU Loss: {loss_ptau.item()}, Thal Loss: {loss_thal.item()}")
    
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
            
            ptau_pred, thal_pred = model(embeddings)
            
            loss_ptau = criterion_ptau(ptau_pred, ptau)
            # Cross-entropy loss for classification
            loss_thal = criterion_thal(thal_pred, thal)
            loss = loss_ptau + loss_thal
            
            total_loss += loss.item()
            
            ptau_preds.extend(ptau_pred.cpu().numpy())
            ptau_targets.extend(ptau.cpu().numpy())
            # Get predicted class from softmax probabilities
            thal_class_preds = torch.argmax(thal_pred, dim=1)
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


def main():
    parser = argparse.ArgumentParser(description='Train baseline ranking model')
    parser.add_argument('--region', type=str, required=True, help='Region name (e.g., EC, ITG, PFC)')
    parser.add_argument('--pooling', type=str, choices=['mean', 'max'], default='mean',
                       help='Pooling method: mean or max')
    parser.add_argument('--loss_ptau', type=str, default='mse',
                       choices=['mse', 'mae', 'huber'],
                       help='Loss function for ptau prediction')
    parser.add_argument('--batch_size', type=int, default=32, help='Batch size')
    parser.add_argument('--epochs', type=int, default=100, help='Number of epochs')
    parser.add_argument('--lr', type=float, default=0.001, help='Learning rate')
    parser.add_argument('--hidden_dim1', type=int, default=64, help='First hidden layer size')
    parser.add_argument('--hidden_dim2', type=int, default=32, help='Second hidden layer size')
    parser.add_argument('--device', type=str, default='cuda' if torch.cuda.is_available() else 'cpu',
                       help='Device to use (cuda or cpu)')
    
    args = parser.parse_args()
    
    # Load data
    data_path = os.path.join('data', 'model_data', f'{args.region}_data_for_model_training.pkl')
    print(f"Loading data from {data_path}")
    
    with open(data_path, 'rb') as f:
        embeddings_folds, patient_info_folds, ptau_folds, thal_folds = pickle.load(f)
    
    print(f"Loaded {len(embeddings_folds)} folds")
    
    # Setup device
    device = torch.device(args.device)
    print(f"Using device: {device}")
    
    # Setup loss functions
    if args.loss_ptau == 'mse':
        criterion_ptau = nn.MSELoss()
    elif args.loss_ptau == 'mae':
        criterion_ptau = nn.L1Loss()
    elif args.loss_ptau == 'huber':
        criterion_ptau = nn.HuberLoss()
    
    # Cross-entropy loss for classification (will be created after determining number of classes)
    
    # Prepare Thal mapping (manual mapping instead of LabelEncoder)
    all_thal_labels = []
    for thal_fold in thal_folds:
        if isinstance(thal_fold, (list, np.ndarray)):
            all_thal_labels.extend(thal_fold)
        else:
            all_thal_labels.extend(thal_fold.values if hasattr(thal_fold, 'values') else thal_fold)
    
    # Convert to numeric if they're strings
    all_thal_labels = [int(float(str(x))) if str(x).isdigit() or ('.' in str(x) and str(x).replace('.', '').isdigit()) else x for x in all_thal_labels]
    
    # Manual mapping: sorted classes to consecutive indices
    thal_classes = sorted(set(all_thal_labels))
    thal_to_idx = {c: i for i, c in enumerate(thal_classes)}
    idx_to_thal = {i: c for c, i in thal_to_idx.items()}
    n_thal_classes = len(thal_classes)
    
    print(f"Thal classes (sorted): {thal_classes}")
    print(f"Thal to index mapping: {thal_to_idx}")
    print(f"Number of Thal classes: {n_thal_classes}")
    print(f"Note: Categories may not be consecutive (e.g., 0, 1, 3, 4, 5)")
    
    # Cross-entropy loss for classification
    criterion_thal = nn.CrossEntropyLoss()
    
    # Prepare patient-level data for all folds first
    n_folds = len(embeddings_folds)
    all_fold_patient_data = []  # List of dicts, one per fold
    
    print(f"\n{'='*60}")
    print("Preparing patient-level data for all folds")
    print(f"{'='*60}")
    
    for fold_idx in range(n_folds):
        print(f"\nPreparing fold {fold_idx + 1}/{n_folds}...")
        
        # Prepare fold data
        fold_embeddings = embeddings_folds[fold_idx]  # (n_cells, n_genes, n_dims)
        fold_ptau = ptau_folds[fold_idx]  # List of values (one per patient)
        fold_thal = thal_folds[fold_idx]  # List of values (one per patient)
        fold_patient_info = patient_info_folds[fold_idx]  # Dict: patient_name -> cell_names
        
        # Apply pooling to get (n_cells, n_dims)
        pooled_embeddings = apply_pooling(fold_embeddings, method=args.pooling)
        print(f"  Pooled embeddings shape: {pooled_embeddings.shape}")
        
        # Map cell-level embeddings to patient-level (since ptau/thal are per-patient)
        # Embeddings are concatenated in patient order matching fold_patient_info
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
            
            # Get ptau and thal for this patient
            patient_ptau.append(fold_ptau[patient_idx])
            patient_thal.append(fold_thal[patient_idx])
            
            cell_idx += n_cells_patient
        
        patient_embeddings = np.array(patient_embeddings)
        patient_ptau = np.array(patient_ptau)
        patient_thal = np.array(patient_thal)
        
        # Apply log1p to ptau
        patient_ptau_log = np.log1p(patient_ptau)
        
        # Encode thal labels using manual mapping
        patient_thal_encoded = np.array([thal_to_idx[t] for t in patient_thal])
        
        print(f"  Patient embeddings shape: {patient_embeddings.shape}")
        print(f"  Number of patients: {len(patient_ptau)}")
        
        # Store fold data
        all_fold_patient_data.append({
            'embeddings': patient_embeddings,
            'ptau': patient_ptau_log,
            'thal': patient_thal_encoded
        })
    
    # Leave-one-fold-out cross-validation
    # For each fold i: train on all folds except i, test on fold i
    fold_results = []
    
    for test_fold_idx in range(n_folds):
        print(f"\n{'='*60}")
        print(f"Test Fold {test_fold_idx + 1}/{n_folds}")
        print(f"Training on folds: {[i+1 for i in range(n_folds) if i != test_fold_idx]}")
        print(f"Testing on fold: {test_fold_idx + 1}")
        print(f"{'='*60}")
        
        # Combine training folds (all except test fold)
        train_fold_indices = [i for i in range(n_folds) if i != test_fold_idx]
        train_embeddings = np.concatenate([all_fold_patient_data[i]['embeddings'] for i in train_fold_indices], axis=0)
        train_ptau = np.concatenate([all_fold_patient_data[i]['ptau'] for i in train_fold_indices], axis=0)
        train_thal = np.concatenate([all_fold_patient_data[i]['thal'] for i in train_fold_indices], axis=0)
        
        # Test fold data
        test_embeddings = all_fold_patient_data[test_fold_idx]['embeddings']
        test_ptau = all_fold_patient_data[test_fold_idx]['ptau']
        test_thal = all_fold_patient_data[test_fold_idx]['thal']
        
        print(f"Train: {len(train_embeddings)} patients, Test: {len(test_embeddings)} patients")
        
        # Create datasets and dataloaders
        train_dataset = EmbeddingDataset(train_embeddings, train_ptau, train_thal)
        test_dataset = EmbeddingDataset(test_embeddings, test_ptau, test_thal)
        
        train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True)
        test_loader = DataLoader(test_dataset, batch_size=args.batch_size, shuffle=False)
        
        # Initialize model
        model = BaselineRankingModel(
            input_dim=train_embeddings.shape[1],
            hidden_dim1=args.hidden_dim1,
            hidden_dim2=args.hidden_dim2,
            n_thal_classes=n_thal_classes
        ).to(device)
        
        optimizer = optim.Adam(model.parameters(), lr=args.lr)
        
        # TensorBoard writer
        log_dir = os.path.join('runs', f'{args.region}_{args.pooling}_{args.loss_ptau}_testfold{test_fold_idx+1}')
        writer = SummaryWriter(log_dir)
        
        # Training loop
        best_test_f1 = 0
        best_epoch = 0
        best_test_mse = None
        best_test_mae = None
        
        for epoch in range(args.epochs):
            train_loss, train_ptau_loss, train_thal_loss = train_epoch(
                model, train_loader, criterion_ptau, criterion_thal, optimizer, device
            )
            
            # Evaluate on test fold
            test_loss, test_mse, test_mae, test_f1 = evaluate(
                model, test_loader, criterion_ptau, criterion_thal, device, idx_to_thal
            )
            
            # Log to TensorBoard
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
    
    # Print summary
    print(f"\n{'='*60}")
    print("Leave-One-Fold-Out Cross-Validation Summary")
    print(f"{'='*60}")
    for result in fold_results:
        print(f"Test Fold {result['test_fold']}: F1={result['best_f1']:.4f}, "
              f"MSE={result['best_test_mse']:.4f}, MAE={result['best_test_mae']:.4f}")
    
    avg_f1 = np.mean([r['best_f1'] for r in fold_results])
    avg_mse = np.mean([r['best_test_mse'] for r in fold_results])
    avg_mae = np.mean([r['best_test_mae'] for r in fold_results])
    
    print(f"\nAverage F1: {avg_f1:.4f}")
    print(f"Average MSE: {avg_mse:.4f}")
    print(f"Average MAE: {avg_mae:.4f}")


if __name__ == "__main__":
    main()
