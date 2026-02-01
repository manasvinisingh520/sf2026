"""Baseline ranking model: ptau regression + thal classification."""

import torch
import torch.nn as nn


class BaselineRankingModel(nn.Module):
    """Baseline model with two heads: ptau regression and thal classification."""

    def __init__(self, input_dim=16, input_shape=None, n_thal_classes=None, use_attention=False, num_heads=1, dropout=0.2):
        super(BaselineRankingModel, self).__init__()

        self.use_attention = use_attention

        # Hardcoded hidden layer dimensions
        hidden_dim1 = 8
        hidden_dim2 = 4
        hidden_dim3 = 4

        # Handle 2D input (n_genes, n_dims) or 1D input (flattened)
        if input_shape is not None and len(input_shape) == 2:
            # 2D input: (n_genes, n_dims) - can apply attention on genes
            n_genes, n_dims = input_shape
            flattened_dim = n_genes * n_dims
            self.flatten_input = nn.Flatten()

            if use_attention:
                # Self-attention on gene dimension before flattening
                # Input: (batch, n_genes, n_dims), Output: (batch, n_genes, n_dims)
                self.attention = nn.MultiheadAttention(
                    embed_dim=n_dims,
                    num_heads=num_heads,
                    batch_first=True,
                    dropout=0.1
                )
            else:
                self.attention = None
            # After attention (if used), flatten to (batch, n_genes * n_dims)
            first_layer_input = flattened_dim
        else:
            # 1D input: already flattened - attention not applicable
            self.flatten_input = nn.Identity()
            self.attention = None
            first_layer_input = input_dim

        # Shared layers
        self.shared = nn.Sequential(
            nn.Linear(first_layer_input, hidden_dim1),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim1, hidden_dim2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim2, hidden_dim3),
            nn.ReLU(),
            nn.Dropout(dropout)
        )

        # Head 1: P-tau regression (ReLU output for log1p(p-tau)), already have log1p applied in the data
        self.ptau_head = nn.Linear(hidden_dim3, 1)

        # Head 2: Thal classification (softmax)
        if n_thal_classes is None:
            raise ValueError("n_thal_classes must be specified")
        # do NOT apply softmax here, it will be applied in the loss function
        self.thal_head = nn.Linear(hidden_dim3, n_thal_classes)

    def forward(self, x):
        # Apply attention only for 2D input (batch, n_genes, n_dims)
        if self.use_attention and self.attention is not None and len(x.shape) == 3:
            # 2D input: (batch, n_genes, n_dims) - apply attention on gene dimension
            x_attn, _ = self.attention(x, x, x)  # Self-attention
            x = self.flatten_input(x_attn)
        else:
            # Flatten if input is 2D, otherwise pass through
            x = self.flatten_input(x)

        shared_features = self.shared(x)
        ptau_pred = self.ptau_head(shared_features).squeeze()
        thal_pred = self.thal_head(shared_features)  # (batch_size, n_thal_classes) with softmax
        return ptau_pred, thal_pred
