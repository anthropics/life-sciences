#!/usr/bin/env python3
"""
Reference-to-Query Label Transfer using scVI/scANVI

This script provides a complete workflow for integrating datasets and
transferring cell type annotations from a labeled reference to an
unlabeled query dataset.

Reference: https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/tabula_muris.html
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
import warnings

# Import modular utilities
from label_transfer_core import (
    normalize_gene_length,
    concatenate_datasets,
    select_hvg_across_batches,
    setup_combined_anndata,
    train_scvi_integration,
    train_scanvi_transfer,
    predict_labels,
    get_prediction_probabilities,
    evaluate_transfer,
    save_model
)

warnings.filterwarnings('ignore')

print("=" * 80)
print("Reference-to-Query Label Transfer with scVI/scANVI")
print("=" * 80)

# Default parameters
DEFAULT_N_LATENT = 30
DEFAULT_N_LAYERS = 2
DEFAULT_N_TOP_GENES = 2000
DEFAULT_MAX_EPOCHS_SCVI = 400
DEFAULT_MAX_EPOCHS_SCANVI = 20
DEFAULT_MIN_DIST = 0.3

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description='Reference-to-Query Label Transfer using scVI/scANVI',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
Examples:
  python3 label_transfer_analysis.py --reference ref.h5ad --query query.h5ad --labels-key cell_type
  python3 label_transfer_analysis.py --reference ref.h5ad --query query.h5ad --labels-key cell_type --batch-key tech
  python3 label_transfer_analysis.py --reference ss2.h5ad --query 10x.h5ad --labels-key cell_type --normalize-gene-length --reference-tech SS2
    """
)

# Required arguments
parser.add_argument('--reference', required=True, help='Path to annotated reference .h5ad file')
parser.add_argument('--query', required=True, help='Path to unannotated query .h5ad file')
parser.add_argument('--labels-key', required=True, help='Column in reference.obs with cell type labels')

# Output
parser.add_argument('--output-dir', type=str, default='label_transfer_results',
                    help='Output directory (default: label_transfer_results)')

# Batch correction
parser.add_argument('--batch-key', type=str, default='dataset',
                    help='Column name for batch information (default: dataset)')
parser.add_argument('--reference-batch', type=str, default='reference',
                    help='Batch label for reference (default: reference)')
parser.add_argument('--query-batch', type=str, default='query',
                    help='Batch label for query (default: query)')

# Technology-specific
parser.add_argument('--normalize-gene-length', action='store_true',
                    help='Apply gene length normalization (for SmartSeq2)')
parser.add_argument('--reference-tech', type=str, choices=['10x', 'SS2', 'other'],
                    default='10x', help='Reference technology (default: 10x)')
parser.add_argument('--query-tech', type=str, choices=['10x', 'SS2', 'other'],
                    default='10x', help='Query technology (default: 10x)')
parser.add_argument('--gene-length-file', type=str, default=None,
                    help='Path to gene length file for normalization')

# Model architecture
parser.add_argument('--n-latent', type=int, default=DEFAULT_N_LATENT,
                    help=f'Latent space dimensions (default: {DEFAULT_N_LATENT})')
parser.add_argument('--n-layers', type=int, default=DEFAULT_N_LAYERS,
                    help=f'Number of hidden layers (default: {DEFAULT_N_LAYERS})')
parser.add_argument('--n-top-genes', type=int, default=DEFAULT_N_TOP_GENES,
                    help=f'Number of HVGs to select (default: {DEFAULT_N_TOP_GENES})')

# Training
parser.add_argument('--max-epochs-scvi', type=int, default=DEFAULT_MAX_EPOCHS_SCVI,
                    help=f'Max epochs for scVI (default: {DEFAULT_MAX_EPOCHS_SCVI})')
parser.add_argument('--max-epochs-scanvi', type=int, default=DEFAULT_MAX_EPOCHS_SCANVI,
                    help=f'Max epochs for scANVI (default: {DEFAULT_MAX_EPOCHS_SCANVI})')

# Downstream
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
for path, name in [(args.reference, 'Reference'), (args.query, 'Query')]:
    if not os.path.exists(path):
        print(f"\nError: {name} file '{path}' not found!")
        sys.exit(1)

# Set up output directory
os.makedirs(args.output_dir, exist_ok=True)
print(f"\nOutput directory: {args.output_dir}")

# Display parameters
print(f"\nParameters:")
print(f"  Reference: {args.reference}")
print(f"  Query: {args.query}")
print(f"  Labels key: {args.labels_key}")
print(f"  Batch key: {args.batch_key}")
print(f"  Technologies: reference={args.reference_tech}, query={args.query_tech}")
print(f"  Model: n_latent={args.n_latent}, n_layers={args.n_layers}")
print(f"  HVGs: {args.n_top_genes}")

# =============================================================================
# Step 1: Load data
# =============================================================================
print("\n" + "=" * 80)
print("[1/9] Loading data...")
print("=" * 80)

reference = ad.read_h5ad(args.reference)
query = ad.read_h5ad(args.query)

print(f"Reference: {reference.n_obs} cells × {reference.n_vars} genes")
print(f"Query: {query.n_obs} cells × {query.n_vars} genes")

# Verify labels key exists in reference
if args.labels_key not in reference.obs.columns:
    print(f"\nError: Labels key '{args.labels_key}' not found in reference.obs")
    print(f"Available columns: {list(reference.obs.columns)}")
    sys.exit(1)

n_cell_types = reference.obs[args.labels_key].nunique()
print(f"Cell types in reference: {n_cell_types}")

# =============================================================================
# Step 2: Gene length normalization (if needed)
# =============================================================================
print("\n" + "=" * 80)
print("[2/9] Preprocessing...")
print("=" * 80)

if args.normalize_gene_length:
    if args.reference_tech == 'SS2':
        print("  Applying gene length normalization to reference (SmartSeq2)...")
        reference = normalize_gene_length(reference, args.gene_length_file)
    if args.query_tech == 'SS2':
        print("  Applying gene length normalization to query (SmartSeq2)...")
        query = normalize_gene_length(query, args.gene_length_file)
else:
    print("  No gene length normalization required")

# =============================================================================
# Step 3: Concatenate datasets
# =============================================================================
print("\n" + "=" * 80)
print("[3/9] Concatenating datasets...")
print("=" * 80)

# Add batch labels
reference.obs[args.batch_key] = args.reference_batch
query.obs[args.batch_key] = args.query_batch

# Add technology labels
reference.obs['tech'] = args.reference_tech
query.obs['tech'] = args.query_tech

# Concatenate
adata = concatenate_datasets(
    reference, query,
    labels_key=args.labels_key,
    batch_key=args.batch_key
)

print(f"Combined: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"  Reference cells: {(adata.obs[args.batch_key] == args.reference_batch).sum()}")
print(f"  Query cells: {(adata.obs[args.batch_key] == args.query_batch).sum()}")

# =============================================================================
# Step 4: HVG selection
# =============================================================================
print("\n" + "=" * 80)
print("[4/9] Selecting highly variable genes...")
print("=" * 80)

# Store counts
adata.layers["counts"] = adata.X.copy()

# Select HVGs
select_hvg_across_batches(
    adata,
    batch_key=args.batch_key,
    n_top_genes=args.n_top_genes
)

print(f"  Selected {adata.n_vars} highly variable genes")

# =============================================================================
# Step 5: Train scVI
# =============================================================================
print("\n" + "=" * 80)
print("[5/9] Training scVI integration model...")
print("=" * 80)

# Create transfer labels column (query cells are "Unknown")
TRANSFER_KEY = "cell_type_transfer"
adata.obs[TRANSFER_KEY] = adata.obs[args.labels_key].copy()
query_mask = adata.obs[args.batch_key] == args.query_batch
adata.obs.loc[query_mask, TRANSFER_KEY] = "Unknown"

# Setup AnnData
setup_combined_anndata(
    adata,
    batch_key=args.batch_key,
    labels_key=TRANSFER_KEY
)

# Train scVI
scvi_model = train_scvi_integration(
    adata,
    n_latent=args.n_latent,
    n_layers=args.n_layers,
    max_epochs=args.max_epochs_scvi
)

# Save scVI model
scvi_path = os.path.join(args.output_dir, 'scvi_model')
save_model(scvi_model, scvi_path)
print(f"  Saved scVI model to: {scvi_path}")

# Get scVI latent representation
SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation()

# =============================================================================
# Step 6: Train scANVI
# =============================================================================
print("\n" + "=" * 80)
print("[6/9] Training scANVI transfer model...")
print("=" * 80)

scanvi_model = train_scanvi_transfer(
    scvi_model,
    adata,
    labels_key=TRANSFER_KEY,
    unlabeled_category="Unknown",
    max_epochs=args.max_epochs_scanvi
)

# Save scANVI model
scanvi_path = os.path.join(args.output_dir, 'scanvi_model')
save_model(scanvi_model, scanvi_path)
print(f"  Saved scANVI model to: {scanvi_path}")

# Get scANVI latent representation
SCANVI_LATENT_KEY = "X_scANVI"
adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation()

# =============================================================================
# Step 7: Predict labels
# =============================================================================
print("\n" + "=" * 80)
print("[7/9] Predicting labels for query cells...")
print("=" * 80)

# Get predictions
predictions = predict_labels(scanvi_model, adata)
adata.obs['predicted_type'] = predictions

# Get confidence
probabilities = get_prediction_probabilities(scanvi_model, adata)
adata.obs['prediction_confidence'] = probabilities.max(axis=1)

# Summary for query cells
query_predictions = adata.obs.loc[query_mask, 'predicted_type']
query_confidence = adata.obs.loc[query_mask, 'prediction_confidence']

print(f"\nQuery cell predictions:")
for ct, count in query_predictions.value_counts().items():
    pct = count / query_mask.sum() * 100
    print(f"  {ct}: {count} ({pct:.1f}%)")

print(f"\nConfidence statistics (query cells):")
print(f"  Mean: {query_confidence.mean():.3f}")
print(f"  Median: {query_confidence.median():.3f}")
print(f"  >0.9: {(query_confidence > 0.9).sum()} cells")
print(f"  >0.7: {(query_confidence > 0.7).sum()} cells")

# =============================================================================
# Step 8: Evaluation
# =============================================================================
print("\n" + "=" * 80)
print("[8/9] Evaluating transfer accuracy...")
print("=" * 80)

# For reference cells, compare predictions to true labels
ref_mask = adata.obs[args.batch_key] == args.reference_batch
ref_true = adata.obs.loc[ref_mask, args.labels_key]
ref_pred = adata.obs.loc[ref_mask, 'predicted_type']

accuracy = (ref_true == ref_pred).mean()
print(f"  Reference cell accuracy: {accuracy:.3f}")

# Compute confusion matrix
metrics = evaluate_transfer(ref_true, ref_pred)
print(f"  Weighted F1 score: {metrics['f1_weighted']:.3f}")

# =============================================================================
# Step 9: Visualization
# =============================================================================
print("\n" + "=" * 80)
print("[9/9] Generating visualizations...")
print("=" * 80)

# Compute UMAP
sc.pp.neighbors(adata, use_rep=SCANVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=args.min_dist)

plt.style.use('default')
sc.settings.set_figure_params(dpi=100, frameon=False)

# Plot 1: Integration UMAP
print("  Creating integration visualization...")
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
sc.pl.umap(adata, color=args.batch_key, ax=axes[0], show=False,
           title='Colored by Dataset')
sc.pl.umap(adata, color=args.labels_key, ax=axes[1], show=False,
           title='Colored by Cell Type (Reference)')
plt.tight_layout()
int_path = os.path.join(args.output_dir, 'integration_umap.png')
plt.savefig(int_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: {int_path}")

# Plot 2: Label transfer UMAP
print("  Creating label transfer visualization...")
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
sc.pl.umap(adata, color='predicted_type', ax=axes[0], show=False,
           title='Predicted Cell Types')
sc.pl.umap(adata, color='prediction_confidence', ax=axes[1], show=False,
           cmap='viridis', title='Prediction Confidence')
plt.tight_layout()
lt_path = os.path.join(args.output_dir, 'label_transfer_umap.png')
plt.savefig(lt_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: {lt_path}")

# Plot 3: Confusion matrix
print("  Creating confusion matrix...")
conf_df = pd.crosstab(ref_true, ref_pred, normalize='index')

fig, ax = plt.subplots(figsize=(12, 10))
sns.heatmap(conf_df, annot=True, fmt='.2f', cmap='Blues', ax=ax)
ax.set_xlabel('Predicted')
ax.set_ylabel('True')
ax.set_title('Confusion Matrix (Reference Cells)')
plt.tight_layout()
cm_path = os.path.join(args.output_dir, 'confusion_matrix.png')
plt.savefig(cm_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: {cm_path}")

# Plot 4: Confidence distribution
print("  Creating confidence visualization...")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

axes[0].hist(query_confidence, bins=50, edgecolor='black', alpha=0.7)
axes[0].axvline(0.9, color='red', linestyle='--', label='0.9')
axes[0].axvline(0.7, color='orange', linestyle='--', label='0.7')
axes[0].set_xlabel('Confidence')
axes[0].set_ylabel('Number of Cells')
axes[0].set_title('Query Cell Confidence Distribution')
axes[0].legend()

# Confidence by predicted type
conf_by_type = adata.obs.loc[query_mask].groupby('predicted_type')['prediction_confidence'].mean()
conf_by_type.plot(kind='bar', ax=axes[1], color='steelblue')
axes[1].set_xlabel('Predicted Type')
axes[1].set_ylabel('Mean Confidence')
axes[1].set_title('Mean Confidence by Predicted Type')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
conf_path = os.path.join(args.output_dir, 'prediction_confidence.png')
plt.savefig(conf_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: {conf_path}")

# =============================================================================
# Save results
# =============================================================================
print("\n" + "=" * 80)
print("Saving results...")
print("=" * 80)

# Save transferred labels
labels_df = adata.obs.loc[query_mask, ['predicted_type', 'prediction_confidence']].copy()
labels_path = os.path.join(args.output_dir, 'transferred_labels.csv')
labels_df.to_csv(labels_path)
print(f"  Saved: {labels_path}")

# Save integrated data
output_h5ad = os.path.join(args.output_dir, 'integrated_data.h5ad')
adata.write(output_h5ad)
print(f"  Saved: {output_h5ad}")

# =============================================================================
# Summary
# =============================================================================
print("\n" + "=" * 80)
print("Label Transfer Complete!")
print("=" * 80)

print(f"\nResults saved to: {args.output_dir}/")
print("\nGenerated files:")
print(f"  1. integration_umap.png - Integration visualization")
print(f"  2. label_transfer_umap.png - Transfer results")
print(f"  3. confusion_matrix.png - Prediction accuracy")
print(f"  4. prediction_confidence.png - Confidence distribution")
print(f"  5. transferred_labels.csv - Predicted labels for query")
print(f"  6. integrated_data.h5ad - Combined annotated dataset")
print(f"  7. scvi_model/ - Saved scVI model")
print(f"  8. scanvi_model/ - Saved scANVI model")

print(f"\nTransfer accuracy: {accuracy:.1%}")
print(f"Query cells annotated: {query_mask.sum()}")

print("\nNext steps:")
print("  - Validate predictions with marker expression")
print("  - Filter low-confidence predictions if needed")
print("  - Use integrated embedding for downstream analysis")
