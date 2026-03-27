#!/usr/bin/env python3
"""
05q_cnn_nmd.py

CNN model predicting NMD status from full transcript sequence + structure.

9 input channels × 8192 positions:
  0-3: A, C, G, T one-hot
  4:   Exon-exon junction indicator
  5:   ATG position indicator (1 at the A of each ATG)
  6:   Stop codon in frame 0 (position % 3 == 0)
  7:   Stop codon in frame 1 (position % 3 == 1)
  8:   Stop codon in frame 2 (position % 3 == 2)

The frame channels give the CNN the building blocks for ORF reasoning:
ATG positions, stop codons in each reading frame, and junction positions.
The CNN must learn to associate these to predict NMD.

Input:  data_mashr/analysis_cache/cnn_data.tsv
Output: data_mashr/analysis_cache/cnn_nmd_results.tsv
        data_mashr/analysis_cache/cnn_nmd_model.pt
"""

import os
import time
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import roc_auc_score

os.chdir("/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/isopair_wrapper")

SEQ_LEN = 4096
N_CHANNELS = 9
HOLDOUT_CHRS = {"chr1", "chr3", "chr5", "chr7"}
BATCH_SIZE = 64   # smaller batch for larger input
EPOCHS = 50
LR = 1e-3
PATIENCE = 7
DEVICE = torch.device("mps" if torch.backends.mps.is_available() else "cpu")

print(f"=== CNN NMD Model (Full Transcript, 9 Channels) ===")
print(f"Device: {DEVICE}")
print(f"Sequence length: {SEQ_LEN} (covers ~70% of transcripts without truncation)")
print(f"Channels: {N_CHANNELS} (ACGT + junction + ATG + stop×3)")
print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

# ==============================================================================
# 1. Load and encode data
# ==============================================================================
print("Loading data...")
df = pd.read_csv("data_mashr/analysis_cache/cnn_data.tsv", sep="\t",
                  dtype={'junctions': str})
df['junctions'] = df['junctions'].fillna('')
print(f"  Isoforms: {len(df)}")
print(f"  NMD: {df['is_nmd'].sum()}  Non-NMD: {(~df['is_nmd'].astype(bool)).sum()}")

NUC_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
STOP_CODONS = {'TAA', 'TAG', 'TGA'}

def encode_transcript(seq_str, junc_str):
    """Encode transcript as 9 × SEQ_LEN matrix."""
    encoded = np.zeros((N_CHANNELS, SEQ_LEN), dtype=np.float32)
    seq_len = min(len(seq_str), SEQ_LEN)

    # Channels 0-3: nucleotide one-hot
    for i in range(seq_len):
        idx = NUC_MAP.get(seq_str[i], -1)
        if idx >= 0:
            encoded[idx, i] = 1.0

    # Channel 4: junction positions
    if junc_str:
        for pos_str in junc_str.split(','):
            pos = int(pos_str) - 1  # 1-based → 0-based
            if 0 <= pos < SEQ_LEN:
                encoded[4, pos] = 1.0

    # Channel 5: ATG positions
    for i in range(seq_len - 2):
        if seq_str[i:i+3] == 'ATG':
            encoded[5, i] = 1.0

    # Channels 6-8: stop codons by frame
    for i in range(seq_len - 2):
        codon = seq_str[i:i+3]
        if codon in STOP_CODONS:
            frame = i % 3
            encoded[6 + frame, i] = 1.0

    return encoded

print("Encoding transcripts (9 channels)...")
t0 = time.time()
sequences = np.stack([encode_transcript(s, j)
                       for s, j in zip(df['seq'], df['junctions'])])
elapsed = time.time() - t0
print(f"  Encoded in {elapsed:.1f} seconds")
print(f"  Shape: {sequences.shape}")
print(f"  Memory: {sequences.nbytes / 1e9:.1f} GB")

# Channel statistics
for ch, name in enumerate(['A', 'C', 'G', 'T', 'Junction', 'ATG',
                            'Stop_f0', 'Stop_f1', 'Stop_f2']):
    mean_count = sequences[:, ch, :].sum(axis=1).mean()
    print(f"    Ch{ch} ({name}): mean {mean_count:.1f} per transcript")

labels = df['is_nmd'].values.astype(np.float32)

# Train/test split
train_mask = ~df['chr'].isin(HOLDOUT_CHRS).values
test_mask = df['chr'].isin(HOLDOUT_CHRS).values
print(f"\n  Train: {train_mask.sum()} (NMD: {int(labels[train_mask].sum())})")
print(f"  Test:  {test_mask.sum()} (NMD: {int(labels[test_mask].sum())})")

# ==============================================================================
# 2. Dataset
# ==============================================================================
class SeqDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.from_numpy(X)
        self.y = torch.from_numpy(y)
    def __len__(self):
        return len(self.y)
    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]

train_ds = SeqDataset(sequences[train_mask], labels[train_mask])
test_ds = SeqDataset(sequences[test_mask], labels[test_mask])

n_pos = int(labels[train_mask].sum())
n_neg = int(train_mask.sum()) - n_pos
pos_weight = torch.tensor([n_neg / n_pos], dtype=torch.float32).to(DEVICE)

train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True)
test_loader = DataLoader(test_ds, batch_size=BATCH_SIZE, shuffle=False)

# ==============================================================================
# 3. Model — deeper CNN for longer input
# ==============================================================================
class NMD_CNN(nn.Module):
    """
    1D CNN on 9-channel transcript encoding (8192 positions).

    The 9 channels provide:
    - Sequence context (ACGT): for learning motifs like Kozak
    - Junction positions: structural information for EJC locations
    - ATG positions: where translation can initiate
    - Stop codons by frame: where ORFs terminate in each reading frame

    The CNN must learn to associate these features spatially to predict NMD.
    """
    def __init__(self):
        super().__init__()
        self.features = nn.Sequential(
            # Block 1: detect local motifs (Kozak, codon patterns)
            nn.Conv1d(N_CHANNELS, 64, kernel_size=15, padding=7),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.MaxPool1d(4),  # 8192 → 2048

            # Block 2: medium-range patterns
            nn.Conv1d(64, 128, kernel_size=9, padding=4),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.MaxPool1d(4),  # 2048 → 512

            # Block 3: longer-range associations (ATG ↔ stop ↔ junction)
            nn.Conv1d(128, 128, kernel_size=7, padding=3),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.MaxPool1d(4),  # 512 → 128

            # Block 4: global patterns
            nn.Conv1d(128, 64, kernel_size=5, padding=2),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.AdaptiveMaxPool1d(1),  # 128 → 1
        )
        self.classifier = nn.Sequential(
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Dropout(0.4),
            nn.Linear(32, 1),
        )

    def forward(self, x):
        x = self.features(x)
        x = x.squeeze(-1)
        return self.classifier(x)

model = NMD_CNN().to(DEVICE)
n_params = sum(p.numel() for p in model.parameters())
print(f"\nModel parameters: {n_params:,}")

# ==============================================================================
# 4. Training
# ==============================================================================
criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
optimizer = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=1e-5)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=3, factor=0.5)

best_auc = 0.0
best_epoch = 0
patience_counter = 0

print(f"\n=== Training ===\n")

for epoch in range(EPOCHS):
    model.train()
    train_loss = 0.0
    n_batches = 0
    for X_batch, y_batch in train_loader:
        X_batch, y_batch = X_batch.to(DEVICE), y_batch.to(DEVICE)
        optimizer.zero_grad()
        logits = model(X_batch).squeeze(-1)
        loss = criterion(logits, y_batch)
        loss.backward()
        optimizer.step()
        train_loss += loss.item()
        n_batches += 1
    train_loss /= n_batches

    model.eval()
    test_preds, test_labels_list = [], []
    with torch.no_grad():
        for X_batch, y_batch in test_loader:
            X_batch = X_batch.to(DEVICE)
            logits = model(X_batch).squeeze(-1)
            probs = torch.sigmoid(logits).cpu().numpy()
            test_preds.extend(probs)
            test_labels_list.extend(y_batch.numpy())

    test_auc = roc_auc_score(test_labels_list, test_preds)
    scheduler.step(1 - test_auc)
    lr = optimizer.param_groups[0]['lr']

    print(f"  Epoch {epoch+1:3d} | Loss: {train_loss:.4f} | "
          f"Test AUC: {test_auc:.4f} | LR: {lr:.6f}")

    if test_auc > best_auc:
        best_auc = test_auc
        best_epoch = epoch + 1
        patience_counter = 0
        torch.save(model.state_dict(), "data_mashr/analysis_cache/cnn_nmd_model.pt")
    else:
        patience_counter += 1
        if patience_counter >= PATIENCE:
            print(f"\n  Early stopping at epoch {epoch+1} (best: {best_epoch})")
            break

print(f"\n  Best test AUC: {best_auc:.4f} at epoch {best_epoch}")

# ==============================================================================
# 5. Final predictions
# ==============================================================================
print("\n=== Final Evaluation ===\n")

model.load_state_dict(torch.load("data_mashr/analysis_cache/cnn_nmd_model.pt",
                                  weights_only=True))
model.eval()

all_ds = SeqDataset(sequences, labels)
all_loader = DataLoader(all_ds, batch_size=BATCH_SIZE, shuffle=False)

all_preds = []
with torch.no_grad():
    for X_batch, _ in all_loader:
        X_batch = X_batch.to(DEVICE)
        logits = model(X_batch).squeeze(-1)
        probs = torch.sigmoid(logits).cpu().numpy()
        all_preds.extend(probs)

all_preds = np.array(all_preds)

train_auc = roc_auc_score(labels[train_mask], all_preds[train_mask])
test_auc = roc_auc_score(labels[test_mask], all_preds[test_mask])
print(f"  Train AUC: {train_auc:.4f}")
print(f"  Test AUC:  {test_auc:.4f}")

print(f"\n=== Model Comparison ===")
print(f"  CNN full transcript 9-channel:    {test_auc:.4f}")
# NOTE: Reference benchmarks below are from unified_model.rds and
# model_comparison.rds. If those models are re-run, update these values
# or load them dynamically. See FEATURE_DICTIONARY.md for feature definitions.
# Last verified: 2026-03-27
print(f"  Unified model (TD2, 11 feat):     0.897")
print(f"  Combined model (TD2+Ref, 24 feat): 0.940")

# Save predictions
results_df = df[['isoform_id', 'is_nmd', 'chr', 'tx_length']].copy()
results_df['cnn_pred'] = all_preds
results_df['in_holdout'] = test_mask
results_df.to_csv("data_mashr/analysis_cache/cnn_nmd_results.tsv",
                   sep="\t", index=False)

print(f"\n  Saved to data_mashr/analysis_cache/cnn_nmd_results.tsv")
print(f"\nCompleted: {time.strftime('%Y-%m-%d %H:%M:%S')}")
