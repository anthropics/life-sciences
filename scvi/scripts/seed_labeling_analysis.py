#!/usr/bin/env python3
"""
Semi-Supervised Cell Annotation using scANVI Seed Labeling

This script provides a complete workflow for annotating cells using marker
gene signatures and scANVI semi-supervised learning.

Reference: https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/seed_labeling.html
"""

import anndata as ad
import scanpy as sc
import scvi
import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
import argparse
import json
import warnings

# Import modular utilities
from seed_labeling_core import (
    load_marker_genes,
    score_cells_by_markers,
    select_seed_cells,
    create_seed_labels,
    train_scvi_model,
    train_scanvi_from_seeds,
    predict_labels,
    get_prediction_confidence,
    save_model
)

warnings.filterwarnings('ignore')

print("=" * 80)
print("Semi-Supervised Cell Annotation with scANVI Seed Labeling")
print("=" * 80)

# Default parameters
DEFAULT_N_SEED_CELLS = 50
DEFAULT_N_LATENT = 30
DEFAULT_N_LAYERS = 2
DEFAULT_MAX_EPOCHS_SCVI = 100
DEFAULT_MAX_EPOCHS_SCANVI = 25
DEFAULT_RESOLUTION = 1.0
DEFAULT_MIN_DIST = 0.3

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description='Semi-Supervised Cell Annotation using scANVI Seed Labeling',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
Examples:
  python3 seed_labeling_analysis.py data.h5ad --markers-file markers.json
  python3 seed_labeling_analysis.py data.h5ad --markers-file markers.json --batch-key batch
  python3 seed_labeling_analysis.py data.h5ad --markers-file markers.json --n-seed-cells 100
    """
)

# Required arguments
parser.add_argument('input_file', help='Input .h5ad file with raw counts')
parser.add_argument('--markers-file', required=True,
                    help='JSON file with marker gene signatures')

# Optional: batch correction
parser.add_argument('--batch-key', type=str, default=None,
                    help='Column for batch correction')

# Output
parser.add_argument('--output-dir', type=str, help='Output directory')

# Seed selection
parser.add_argument('--n-seed-cells', type=int, default=DEFAULT_N_SEED_CELLS,
                    help=f'Number of seed cells per type (default: {DEFAULT_N_SEED_CELLS})')
parser.add_argument('--min-score-percentile', type=float, default=95,
                    help='Minimum percentile for seed selection (default: 95)')

# Model architecture
parser.add_argument('--n-latent', type=int, default=DEFAULT_N_LATENT,
                    help=f'Latent space dimensions (default: {DEFAULT_N_LATENT})')
parser.add_argument('--n-layers', type=int, default=DEFAULT_N_LAYERS,
                    help=f'Number of hidden layers (default: {DEFAULT_N_LAYERS})')

# Training
parser.add_argument('--max-epochs-scvi', type=int, default=DEFAULT_MAX_EPOCHS_SCVI,
                    help=f'Max epochs for scVI (default: {DEFAULT_MAX_EPOCHS_SCVI})')
parser.add_argument('--max-epochs-scanvi', type=int, default=DEFAULT_MAX_EPOCHS_SCANVI,
                    help=f'Max epochs for scANVI (default: {DEFAULT_MAX_EPOCHS_SCANVI})')

# Downstream
parser.add_argument('--resolution', type=float, default=DEFAULT_RESOLUTION,
                    help=f'Leiden clustering resolution (default: {DEFAULT_RESOLUTION})')
parser.add_argument('--min-dist', type=float, default=DEFAULT_MIN_DIST,
                    help=f'UMAP min_dist parameter (default: {DEFAULT_MIN_DIST})')

# Advanced
parser.add_argument('--seed', type=int, default=0,
                    help='Random seed (default: 0)')

args = parser.parse_args()

# Set random seeds
scvi.settings.seed = args.seed
np.random.seed(args.seed)

# Configure PyTorch
if torch.cuda.is_available():
    torch.set_float32_matmul_precision("high")
    print(f"\nGPU detected: {torch.cuda.get_device_name(0)}")
else:
    print("\nNo GPU detected, using CPU (training may be slow)")

# Verify input files
if not os.path.exists(args.input_file):
    print(f"\nError: File '{args.input_file}' not found!")
    sys.exit(1)

if not os.path.exists(args.markers_file):
    print(f"\nError: Markers file '{args.markers_file}' not found!")
    sys.exit(1)

input_file = args.input_file
base_name = os.path.splitext(os.path.basename(input_file))[0]

# Set up output directory
if args.output_dir:
    output_dir = args.output_dir
else:
    output_dir = f"{base_name}_seed_labeling_results"

os.makedirs(output_dir, exist_ok=True)
print(f"\nOutput directory: {output_dir}")

# Display parameters
print(f"\nParameters:")
print(f"  Markers file: {args.markers_file}")
print(f"  Batch key: {args.batch_key or 'None'}")
print(f"  Seed cells per type: {args.n_seed_cells}")
print(f"  Model: n_latent={args.n_latent}, n_layers={args.n_layers}")
print(f"  Training: scVI={args.max_epochs_scvi} epochs, scANVI={args.max_epochs_scanvi} epochs")

# =============================================================================
# Step 1: Load data and markers
# =============================================================================
print("\n" + "=" * 80)
print("[1/7] Loading data and marker genes...")
print("=" * 80)

adata = ad.read_h5ad(input_file)
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

# Load marker genes
markers = load_marker_genes(args.markers_file)
print(f"\nLoaded markers for {len(markers)} cell types:")
for ct, genes in markers.items():
    n_pos = len(genes.get('positive', []))
    n_neg = len(genes.get('negative', []))
    print(f"  {ct}: {n_pos} positive, {n_neg} negative markers")

# Verify markers exist in data
all_markers = set()
for ct, genes in markers.items():
    all_markers.update(genes.get('positive', []))
    all_markers.update(genes.get('negative', []))

missing = all_markers - set(adata.var_names)
if missing:
    print(f"\nWarning: {len(missing)} markers not found in data: {missing}")
    # Filter markers to only include present genes
    for ct in markers:
        markers[ct]['positive'] = [g for g in markers[ct].get('positive', []) if g in adata.var_names]
        markers[ct]['negative'] = [g for g in markers[ct].get('negative', []) if g in adata.var_names]

# =============================================================================
# Step 2: Score cells by markers
# =============================================================================
print("\n" + "=" * 80)
print("[2/7] Scoring cells by marker gene expression...")
print("=" * 80)

scores = score_cells_by_markers(adata, markers)

# Save scores
scores_df = pd.DataFrame(scores, index=adata.obs_names)
scores_path = os.path.join(output_dir, 'marker_scores.csv')
scores_df.to_csv(scores_path)
print(f"  Saved marker scores to: {scores_path}")

# =============================================================================
# Step 3: Select seed cells
# =============================================================================
print("\n" + "=" * 80)
print("[3/7] Selecting seed cells...")
print("=" * 80)

seed_masks = select_seed_cells(
    scores,
    n_cells=args.n_seed_cells,
    min_percentile=args.min_score_percentile
)

# Create seed labels
create_seed_labels(adata, seed_masks, list(markers.keys()))

# Count seeds
seed_counts = adata.obs['seed_labels'].value_counts()
print(f"\nSeed label distribution:")
for label, count in seed_counts.items():
    print(f"  {label}: {count}")

# =============================================================================
# Step 4: Train scVI model
# =============================================================================
print("\n" + "=" * 80)
print("[4/7] Training scVI base model...")
print("=" * 80)

# Setup AnnData
scvi.model.SCVI.setup_anndata(
    adata,
    batch_key=args.batch_key,
    labels_key='seed_labels'
)

scvi_model = train_scvi_model(
    adata,
    n_latent=args.n_latent,
    n_layers=args.n_layers,
    max_epochs=args.max_epochs_scvi
)

# Save scVI model
scvi_path = os.path.join(output_dir, 'scvi_model')
save_model(scvi_model, scvi_path)
print(f"  Saved scVI model to: {scvi_path}")

# =============================================================================
# Step 5: Train scANVI model
# =============================================================================
print("\n" + "=" * 80)
print("[5/7] Training scANVI classifier with seed labels...")
print("=" * 80)

scanvi_model = train_scanvi_from_seeds(
    scvi_model,
    adata,
    labels_key='seed_labels',
    unlabeled_category='Unknown',
    max_epochs=args.max_epochs_scanvi
)

# Save scANVI model
scanvi_path = os.path.join(output_dir, 'scanvi_model')
save_model(scanvi_model, scanvi_path)
print(f"  Saved scANVI model to: {scanvi_path}")

# =============================================================================
# Step 6: Predict labels and compute confidence
# =============================================================================
print("\n" + "=" * 80)
print("[6/7] Predicting cell type labels...")
print("=" * 80)

# Get predictions
predictions = predict_labels(scanvi_model, adata)
adata.obs['predicted_type'] = predictions

# Get confidence scores
confidence = get_prediction_confidence(scanvi_model, adata)
adata.obs['prediction_confidence'] = confidence.max(axis=1)

# Store confidence per class
for i, ct in enumerate(scanvi_model.adata.obs['seed_labels'].cat.categories):
    if ct != 'Unknown':
        adata.obs[f'prob_{ct}'] = confidence[:, i]

# Print prediction summary
pred_counts = adata.obs['predicted_type'].value_counts()
print(f"\nPrediction distribution:")
for label, count in pred_counts.items():
    pct = count / adata.n_obs * 100
    print(f"  {label}: {count} ({pct:.1f}%)")

# Confidence summary
print(f"\nPrediction confidence:")
print(f"  Mean: {adata.obs['prediction_confidence'].mean():.3f}")
print(f"  Median: {adata.obs['prediction_confidence'].median():.3f}")
print(f"  >0.9: {(adata.obs['prediction_confidence'] > 0.9).sum()} cells")
print(f"  >0.7: {(adata.obs['prediction_confidence'] > 0.7).sum()} cells")

# =============================================================================
# Step 7: Downstream analysis and visualization
# =============================================================================
print("\n" + "=" * 80)
print("[7/7] Computing UMAP and generating visualizations...")
print("=" * 80)

# Get latent representation
SCANVI_LATENT_KEY = "X_scANVI"
adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(adata)

# Compute neighbors and UMAP
sc.pp.neighbors(adata, use_rep=SCANVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=args.min_dist)
sc.tl.leiden(adata, resolution=args.resolution)

# Set up plotting
plt.style.use('default')
sc.settings.set_figure_params(dpi=100, frameon=False)

# Plot 1: Seed labels
print("  Creating seed labels visualization...")
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(adata, color='seed_labels', ax=ax, show=False,
           title='Seed Labels (marker-based selection)')
plt.tight_layout()
seed_path = os.path.join(output_dir, 'seed_labels_umap.png')
plt.savefig(seed_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: {seed_path}")

# Plot 2: Predictions
print("  Creating predictions visualization...")
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(adata, color='predicted_type', ax=ax, show=False,
           title='Predicted Cell Types (scANVI)')
plt.tight_layout()
pred_path = os.path.join(output_dir, 'predictions_umap.png')
plt.savefig(pred_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: {pred_path}")

# Plot 3: Comparison
print("  Creating comparison visualization...")
fig, axes = plt.subplots(1, 2, figsize=(16, 7))
sc.pl.umap(adata, color='seed_labels', ax=axes[0], show=False,
           title='Seed Labels')
sc.pl.umap(adata, color='predicted_type', ax=axes[1], show=False,
           title='Predictions')
plt.tight_layout()
comp_path = os.path.join(output_dir, 'comparison_umap.png')
plt.savefig(comp_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: {comp_path}")

# Plot 4: Confidence distribution
print("  Creating confidence visualization...")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Histogram
axes[0].hist(adata.obs['prediction_confidence'], bins=50, edgecolor='black')
axes[0].axvline(0.9, color='red', linestyle='--', label='0.9 threshold')
axes[0].axvline(0.7, color='orange', linestyle='--', label='0.7 threshold')
axes[0].set_xlabel('Prediction Confidence')
axes[0].set_ylabel('Number of Cells')
axes[0].set_title('Confidence Distribution')
axes[0].legend()

# UMAP colored by confidence
sc.pl.umap(adata, color='prediction_confidence', ax=axes[1], show=False,
           cmap='viridis', title='Prediction Confidence')

plt.tight_layout()
conf_path = os.path.join(output_dir, 'prediction_confidence.png')
plt.savefig(conf_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: {conf_path}")

# =============================================================================
# Save annotated data
# =============================================================================
print("\n" + "=" * 80)
print("Saving annotated data...")
print("=" * 80)

output_h5ad = os.path.join(output_dir, f'{base_name}_annotated.h5ad')
adata.write(output_h5ad)
print(f"  Saved: {output_h5ad}")

# =============================================================================
# Summary
# =============================================================================
print("\n" + "=" * 80)
print("Seed Labeling Analysis Complete!")
print("=" * 80)

print(f"\nResults saved to: {output_dir}/")
print("\nGenerated files:")
print(f"  1. seed_labels_umap.png - Seed label visualization")
print(f"  2. predictions_umap.png - Predicted cell types")
print(f"  3. comparison_umap.png - Side-by-side comparison")
print(f"  4. prediction_confidence.png - Confidence distribution")
print(f"  5. marker_scores.csv - Cell scores per marker set")
print(f"  6. {base_name}_annotated.h5ad - Annotated dataset")
print(f"  7. scvi_model/ - Saved scVI model")
print(f"  8. scanvi_model/ - Saved scANVI model")

print(f"\nAnnotations stored in adata.obs:")
print(f"  - 'seed_labels': Original seed labels")
print(f"  - 'predicted_type': scANVI predictions")
print(f"  - 'prediction_confidence': Max probability")

print("\nNext steps:")
print("  - Validate predictions with marker expression")
print("  - Check low-confidence cells for manual curation")
print("  - Use predictions for downstream analysis")
