
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GREmLN-style Fine-tuning Template (PyTorch)
-------------------------------------------
This script fine-tunes a GREmLN-like model on single-cell RNA data with a graph-aware
inductive bias. It supports three objectives:
  1) Cell-state classification (e.g., disease vs healthy)  -> Cross-entropy
  2) Expression reconstruction (optional)                  -> MSE
  3) Graph smoothness regularization on gene embeddings    -> Laplacian penalty

You can run this file directly to test on synthetic data:
    python gremln_finetune_template.py --synthetic

Or adapt the `load_real_data(...)` function to return your real matrices:
  - X_train, X_val:        (num_cells, num_genes) expression (float32)
  - y_train, y_val:        (num_cells,) integer labels in [0..K-1]
  - A:                     (num_genes, num_genes) adjacency (0/1 or weights)
  - (optional) mask_genes: boolean mask of genes to include (shape num_genes)

Model summary
-------------
We implement a compact "graph-aware transformer-ish encoder":
  - Learnable gene token embeddings (num_genes x d_model)
  - Cell-specific token mixing via gated attention biased by adjacency (A)
  - Cell embedding via mean pooling over per-gene hidden states
  - Heads:
      * Classification head on pooled cell embedding
      * Reconstruction head mapping hidden states back to expression

IMPORTANT: This is a faithful *template* capturing the fine-tuning objectives and training
plumbing; it is NOT the official GREmLN code. It is designed to be clear, hackable, and
strongly typed for your dataset.

Author: ChatGPT (MIT Course 6 vibe ☺)
"""

import argparse
import math
import os
import random
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset, random_split

# -----------------------------
# Utilities
# -----------------------------

def set_seed(seed: int = 1337):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

def to_device(*tensors, device):
    return [t.to(device) for t in tensors]

def train_val_split(X, y, val_frac=0.15, seed=1337):
    n = X.shape[0]
    idx = np.arange(n)
    rng = np.random.default_rng(seed)
    rng.shuffle(idx)
    cut = int(n * (1 - val_frac))
    tr_idx, va_idx = idx[:cut], idx[cut:]
    return X[tr_idx], y[tr_idx], X[va_idx], y[va_idx]

# -----------------------------
# Data loading placeholders
# -----------------------------

def load_synthetic_data(num_cells=4000, num_genes=2000, n_classes=2, p_edge=0.001, seed=42):
    """
    Generates a toy dataset:
      - Expression matrix X ~ ReLU(W z + noise), where z depends on class
      - Sparse random adjacency A with p_edge
    """
    rng = np.random.default_rng(seed)
    # Random low-rank latent
    z = rng.standard_normal((num_cells, 8)).astype(np.float32)
    # Class labels affect z mean
    y = rng.integers(0, n_classes, size=num_cells)
    for k in range(n_classes):
        z[y == k] += (k + 1) * 0.5

    W = rng.standard_normal((8, num_genes)).astype(np.float32) * 0.5
    X = z @ W + 0.1 * rng.standard_normal((num_cells, num_genes)).astype(np.float32)
    X = np.clip(X, 0, None)  # ReLU-like nonnegativity
    # scale per cell
    X = X / (X.sum(axis=1, keepdims=True) + 1e-6) * 1e4
    X = np.log1p(X).astype(np.float32)

    # Random symmetric adjacency (undirected)
    A = (rng.random((num_genes, num_genes)) < p_edge).astype(np.float32)
    A = np.triu(A, 1)
    A = A + A.T
    np.fill_diagonal(A, 0.0)

    # Train/val split
    X_train, y_train, X_val, y_val = train_val_split(X, y, val_frac=0.2, seed=seed)
    return X_train, y_train, X_val, y_val, A

def load_real_data():
    """
    Replace this with your real loading logic.
    Must return: X_train, y_train, X_val, y_val, A
    Shapes:
      - X_* : (num_cells, num_genes) float32 matrix
      - y_* : (num_cells,) int labels in [0..K-1]
      - A   : (num_genes, num_genes) adjacency (float32)
    """
    raise NotImplementedError("Please implement load_real_data() for your dataset.")

# -----------------------------
# Graph utilities
# -----------------------------

def normalized_laplacian(A: torch.Tensor, eps: float = 1e-6) -> torch.Tensor:
    """
    Compute symmetric normalized Laplacian: L = I - D^{-1/2} A D^{-1/2}
    A: (G, G) adjacency, nonnegative
    """
    deg = A.sum(dim=1)
    D_inv_sqrt = torch.diag(torch.pow(deg + eps, -0.5))
    I = torch.eye(A.size(0), device=A.device, dtype=A.dtype)
    L = I - D_inv_sqrt @ A @ D_inv_sqrt
    return L

# -----------------------------
# Model
# -----------------------------

class GraphAwareMixer(nn.Module):
    """
    A lightweight, graph-biased token mixer.
    Given per-gene hidden states H (B, G, D) and adjacency A (G, G),
    we compute a gated message passing that respects A.

    H' = H + sigmoid(Gate) * (A_norm @ H) W_mlp

    Where A_norm is degree-normalized adjacency.
    """
    def __init__(self, d_model: int, dropout: float = 0.1):
        super().__init__()
        self.proj = nn.Linear(d_model, d_model)
        self.gate = nn.Linear(d_model, d_model)
        self.ffn = nn.Sequential(
            nn.Linear(d_model, 4 * d_model),
            nn.GELU(),
            nn.Linear(4 * d_model, d_model)
        )
        self.dropout = nn.Dropout(dropout)
        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)

    def forward(self, H: torch.Tensor, A_norm: torch.Tensor) -> torch.Tensor:
        # H: (B, G, D), A_norm: (G, G)
        B, G, D = H.shape
        # Graph-biased mixing
        gate = torch.sigmoid(self.gate(H))     # (B, G, D)
        mixed = torch.matmul(A_norm, H)        # (B, G, D) since broadcast A_norm over batch
        mixed = self.proj(mixed)               # (B, G, D)
        H2 = H + self.dropout(gate * mixed)
        H2 = self.norm1(H2)
        # FFN
        H3 = H2 + self.dropout(self.ffn(H2))
        H3 = self.norm2(H3)
        return H3

class GremlnLite(nn.Module):
    """
    GREmLN-like encoder with:
      - Learnable gene embeddings E_g (G, D)
      - Input projection of expression to D
      - N graph-aware mixer blocks
      - Cell embedding via mean pooling
      - Heads: classification (cell-level) and reconstruction (gene-level)
    """
    def __init__(self, num_genes: int, n_classes: int, d_model: int = 256,
                 n_blocks: int = 3, dropout: float = 0.1, recon_head: bool = True):
        super().__init__()
        self.num_genes = num_genes
        self.n_classes = n_classes
        self.d_model = d_model
        self.recon_head = recon_head

        self.gene_token = nn.Parameter(torch.randn(num_genes, d_model) * 0.02)
        self.expr_proj = nn.Linear(1, d_model)   # project each gene's scalar expr to D

        self.blocks = nn.ModuleList([GraphAwareMixer(d_model, dropout) for _ in range(n_blocks)])
        self.pool = nn.AdaptiveAvgPool1d(1)  # pool over gene dimension after permuting

        self.cls_head = nn.Sequential(
            nn.LayerNorm(d_model),
            nn.Linear(d_model, n_classes)
        )
        if self.recon_head:
            self.recon = nn.Sequential(
                nn.LayerNorm(d_model),
                nn.Linear(d_model, 1)  # per-gene reconstruction
            )

    def forward(self, X: torch.Tensor, A_norm: torch.Tensor) -> Tuple[torch.Tensor, Optional[torch.Tensor], torch.Tensor]:
        """
        X:      (B, G) float32 expression
        A_norm: (G, G) normalized adjacency (shared across batch)
        Returns:
          logits: (B, K)
          X_hat:  (B, G) if recon_head else None
          H:      (B, G, D) final per-gene hidden states (for analyses/regularization)
        """
        B, G = X.shape
        assert G == self.num_genes, "X gene dimension mismatch"

        # Expand expression per gene to (B, G, 1) then project to D
        X_in = X.unsqueeze(-1)               # (B, G, 1)
        X_emb = self.expr_proj(X_in)         # (B, G, D)

        # Add learnable gene tokens
        H = X_emb + self.gene_token.unsqueeze(0)  # (B, G, D)

        # Graph-aware mixer blocks
        for blk in self.blocks:
            H = blk(H, A_norm)

        # Pool to cell embedding (mean over genes)
        # pool expects (B, D, G)
        H_perm = H.permute(0, 2, 1)              # (B, D, G)
        cell_emb = self.pool(H_perm).squeeze(-1)  # (B, D)

        logits = self.cls_head(cell_emb)          # (B, K)

        X_hat = None
        if self.recon_head:
            X_hat = self.recon(H).squeeze(-1)    # (B, G)

        return logits, X_hat, H

# -----------------------------
# Losses
# -----------------------------

@dataclass
class LossWeights:
    lambda_class: float = 1.0
    lambda_recon: float = 1.0
    lambda_graph: float = 1.0

def laplacian_smoothness(H: torch.Tensor, L: torch.Tensor) -> torch.Tensor:
    """
    Laplacian smoothness: sum over batch of Tr(H^T L H) where H is (B, G, D)
    Encourages connected genes to have similar embeddings.
    """
    B, G, D = H.shape
    # (B, G, D) -> (B, D, G) for Tr calculation
    H_t = H.transpose(1, 2)               # (B, D, G)
    LH = torch.matmul(L, H)               # (B, G, D)
    # For each batch: Tr(H^T L H) = sum over d of h_d^T L h_d
    # Implement as elementwise product and sum
    smooth = (H_t.transpose(1,2) * LH).sum(dim=(1,2))  # (B,)
    return smooth.mean()

def compute_losses(logits, y, X_hat, X, H, L, weights: LossWeights):
    losses = {}
    # Classification
    if weights.lambda_class > 0:
        losses["class"] = F.cross_entropy(logits, y)
    else:
        losses["class"] = torch.tensor(0.0, device=logits.device)

    # Reconstruction
    if X_hat is not None and weights.lambda_recon > 0:
        losses["recon"] = F.mse_loss(X_hat, X)
    else:
        losses["recon"] = torch.tensor(0.0, device=logits.device)

    # Graph Laplacian smoothness
    if weights.lambda_graph > 0:
        losses["graph"] = laplacian_smoothness(H, L)
    else:
        losses["graph"] = torch.tensor(0.0, device=logits.device)

    total = (weights.lambda_class * losses["class"]
             + weights.lambda_recon * losses["recon"]
             + weights.lambda_graph * losses["graph"])
    return total, losses

# -----------------------------
# Training
# -----------------------------

def make_loaders(X_train, y_train, X_val, y_val, batch_size=64):
    train_ds = TensorDataset(torch.from_numpy(X_train), torch.from_numpy(y_train).long())
    val_ds   = TensorDataset(torch.from_numpy(X_val),   torch.from_numpy(y_val).long())
    train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True, drop_last=False)
    val_loader   = DataLoader(val_ds,   batch_size=batch_size, shuffle=False, drop_last=False)
    return train_loader, val_loader

def accuracy(logits, y):
    preds = logits.argmax(dim=1)
    return (preds == y).float().mean().item()

def train_model(X_train, y_train, X_val, y_val, A,
                d_model=256, n_blocks=3, dropout=0.1,
                lr=2e-3, weight_decay=1e-4, batch_size=64,
                max_epochs=50, patience=6,
                lambda_class=1.0, lambda_recon=1.0, lambda_graph=0.1,
                device: Optional[str] = None):
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    G = A.shape[0]
    n_classes = int(np.max(y_train) + 1)

    # Dataloaders
    train_loader, val_loader = make_loaders(X_train, y_train, X_val, y_val, batch_size=batch_size)

    # Model & graph matrices
    model = GremlnLite(num_genes=G, n_classes=n_classes, d_model=d_model, n_blocks=n_blocks, dropout=dropout, recon_head=(lambda_recon>0)).to(device)
    A_t = torch.from_numpy(A).to(device)
    # normalized adjacency for mixer
    deg = A_t.sum(dim=1, keepdim=True) + 1e-6
    A_norm = A_t / deg
    # Laplacian for smoothness
    L = normalized_laplacian(A_t)

    # Optimizer & scheduler
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=max_epochs)

    weights = LossWeights(lambda_class=lambda_class, lambda_recon=lambda_recon, lambda_graph=lambda_graph)

    best_val = float("inf")
    best_state = None
    bad_epochs = 0

    for epoch in range(1, max_epochs + 1):
        model.train()
        tr_loss = 0.0
        tr_acc  = 0.0
        n_tr_batches = 0

        for Xb, yb in train_loader:
            Xb, yb = to_device(Xb, yb, device=device)
            opt.zero_grad(set_to_none=True)
            logits, X_hat, H = model(Xb, A_norm)
            total, parts = compute_losses(logits, yb, X_hat, Xb, H, L, weights)
            total.backward()
            nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            opt.step()

            tr_loss += total.item()
            tr_acc  += accuracy(logits, yb)
            n_tr_batches += 1

        tr_loss /= max(1, n_tr_batches)
        tr_acc  /= max(1, n_tr_batches)

        # Validation
        model.eval()
        va_loss = 0.0
        va_acc  = 0.0
        n_va_batches = 0
        with torch.no_grad():
            for Xb, yb in val_loader:
                Xb, yb = to_device(Xb, yb, device=device)
                logits, X_hat, H = model(Xb, A_norm)
                total, parts = compute_losses(logits, yb, X_hat, Xb, H, L, weights)
                va_loss += total.item()
                va_acc  += accuracy(logits, yb)
                n_va_batches += 1
        va_loss /= max(1, n_va_batches)
        va_acc  /= max(1, n_va_batches)

        sched.step()

        print(f"Epoch {epoch:03d} | train loss {tr_loss:.4f} acc {tr_acc:.3f} | val loss {va_loss:.4f} acc {va_acc:.3f} | lr {sched.get_last_lr()[0]:.2e}")

        # Early stopping on val loss
        if va_loss + 1e-6 < best_val:
            best_val = va_loss
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            bad_epochs = 0
        else:
            bad_epochs += 1
            if bad_epochs >= patience:
                print(f"Early stopping at epoch {epoch}. Best val loss: {best_val:.4f}")
                break

    if best_state is not None:
        model.load_state_dict(best_state)

    return model

# -----------------------------
# Main
# -----------------------------

def main():
    parser = argparse.ArgumentParser(description="GREmLN-style fine-tuning")
    parser.add_argument("--synthetic", action="store_true", help="Use synthetic data")
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch", type=int, default=64)
    parser.add_argument("--d_model", type=int, default=256)
    parser.add_argument("--blocks", type=int, default=3)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--lr", type=float, default=2e-3)
    parser.add_argument("--wd", type=float, default=1e-4)
    parser.add_argument("--lambda_class", type=float, default=1.0)
    parser.add_argument("--lambda_recon", type=float, default=1.0)
    parser.add_argument("--lambda_graph", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=1337)
    args = parser.parse_args()

    set_seed(args.seed)

    if args.synthetic:
        X_train, y_train, X_val, y_val, A = load_synthetic_data(
            num_cells=4000, num_genes=2000, n_classes=2, p_edge=0.002, seed=42
        )
    else:
        X_train, y_train, X_val, y_val, A = load_real_data()

    model = train_model(
        X_train, y_train, X_val, y_val, A,
        d_model=args.d_model, n_blocks=args.blocks, dropout=args.dropout,
        lr=args.lr, weight_decay=args.wd, batch_size=args.batch,
        max_epochs=args.epochs, patience=6,
        lambda_class=args.lambda_class, lambda_recon=args.lambda_recon, lambda_graph=args.lambda_graph
    )

    # Save the fine-tuned weights
    os.makedirs("checkpoints", exist_ok=True)
    ckpt_path = f"checkpoints/gremln_lite_finetuned.pt"
    torch.save(model.state_dict(), ckpt_path)
    print(f"\n✅ Saved fine-tuned weights to: {ckpt_path}")

    # Also save the config for reproducibility
    with open("checkpoints/config.txt", "w") as f:
        f.write(str(vars(args)))
    print("📝 Saved training config to checkpoints/config.txt")

if __name__ == "__main__":
    main()
