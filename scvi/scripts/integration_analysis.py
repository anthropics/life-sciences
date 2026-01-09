#!/usr/bin/env python3
"""
Single-Cell Data Integration Analysis using scvi-tools

This script provides a complete integration workflow using scVI (unsupervised)
and scANVI (semi-supervised) models following scvi-tools best practices.

Reference: https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/harmonization.html
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
from integration_core import (
    setup_anndata_scvi,
    train_scvi_model,
    train_scanvi_model,
    get_latent_representation,
    compute_neighbors_and_umap,
    cluster_cells,
    save_model
)
from integration_metrics import calculate_integration_metrics, format_metrics_table

warnings.filterwarnings('ignore')

print("=" * 80)
print("Single-Cell Data Integration Analysis")
print("Using scVI and scANVI from scvi-tools")
print("=" * 80)

# Default parameters
DEFAULT_N_LAYERS = 2
DEFAULT_N_LATENT = 30
DEFAULT_GENE_LIKELIHOOD = 'nb'
DEFAULT_MAX_EPOCHS_SCANVI = 20
DEFAULT_RESOLUTION = 1.0
DEFAULT_N_NEIGHBORS = 15
DEFAULT_MIN_DIST = 0.3
DEFAULT_UNLABELED = 'Unknown'

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description='Single-Cell Data Integration using scVI/scANVI',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
Examples:
  python3 integration_analysis.py data.h5ad --batch-key batch
  python3 integration_analysis.py data.h5ad --batch-key sample --labels-key cell_type
  python3 integration_analysis.py data.h5ad --batch-key batch --n-latent 50 --n-layers 3
  python3 integration_analysis.py data.h5ad --batch-key batch --skip-scanvi
    """
)

# Required arguments
parser.add_argument('input_file', help='Input .h5ad file with raw counts')
parser.add_argument('--batch-key', required=True, help='Column in adata.obs containing batch identifiers')

# Optional: cell type annotations
parser.add_argument('--labels-key', type=str, default=None,
                    help='Column with cell type annotations (enables scANVI)')
parser.add_argument('--skip-scanvi', action='store_true',
                    help='Skip scANVI even if labels are available')
parser.add_argument('--unlabeled-category', type=str, default=DEFAULT_UNLABELED,
                    help=f'Label for unlabeled cells (default: {DEFAULT_UNLABELED})')

# Output
parser.add_argument('--output-dir', type=str, help='Output directory')

# Model architecture
parser.add_argument('--n-layers', type=int, default=DEFAULT_N_LAYERS,
                    help=f'Number of hidden layers (default: {DEFAULT_N_LAYERS})')
parser.add_argument('--n-latent', type=int, default=DEFAULT_N_LATENT,
                    help=f'Latent space dimensions (default: {DEFAULT_N_LATENT})')
parser.add_argument('--gene-likelihood', type=str, default=DEFAULT_GENE_LIKELIHOOD,
                    choices=['nb', 'zinb', 'poisson'],
                    help=f'Gene likelihood distribution (default: {DEFAULT_GENE_LIKELIHOOD})')

# Training
parser.add_argument('--max-epochs-scvi', type=int, default=None,
                    help='Max epochs for scVI (default: auto)')
parser.add_argument('--max-epochs-scanvi', type=int, default=DEFAULT_MAX_EPOCHS_SCANVI,
                    help=f'Max epochs for scANVI (default: {DEFAULT_MAX_EPOCHS_SCANVI})')
parser.add_argument('--no-early-stopping', action='store_true',
                    help='Disable early stopping')

# Downstream analysis
parser.add_argument('--resolution', type=float, default=DEFAULT_RESOLUTION,
                    help=f'Leiden clustering resolution (default: {DEFAULT_RESOLUTION})')
parser.add_argument('--n-neighbors', type=int, default=DEFAULT_N_NEIGHBORS,
                    help=f'Number of neighbors for kNN graph (default: {DEFAULT_N_NEIGHBORS})')
parser.add_argument('--min-dist', type=float, default=DEFAULT_MIN_DIST,
                    help=f'UMAP min_dist parameter (default: {DEFAULT_MIN_DIST})')

# Advanced
parser.add_argument('--layer', type=str, default=None,
                    help='Layer containing counts (default: auto-detect)')
parser.add_argument('--seed', type=int, default=0,
                    help='Random seed for reproducibility (default: 0)')

args = parser.parse_args()

# Set random seeds
scvi.settings.seed = args.seed
np.random.seed(args.seed)

# Configure PyTorch for better performance
if torch.cuda.is_available():
    torch.set_float32_matmul_precision("high")
    print(f"\nGPU detected: {torch.cuda.get_device_name(0)}")
else:
    print("\nNo GPU detected, using CPU (training may be slow)")

# Verify input file
if not os.path.exists(args.input_file):
    print(f"\nError: File '{args.input_file}' not found!")
    sys.exit(1)

input_file = args.input_file
base_name = os.path.splitext(os.path.basename(input_file))[0]

# Set up output directory
if args.output_dir:
    output_dir = args.output_dir
else:
    output_dir = f"{base_name}_integration_results"

os.makedirs(output_dir, exist_ok=True)
print(f"\nOutput directory: {output_dir}")

# Display parameters
print(f"\nParameters:")
print(f"  Batch key: {args.batch_key}")
print(f"  Labels key: {args.labels_key or 'None (scVI only)'}")
print(f"  Model: n_layers={args.n_layers}, n_latent={args.n_latent}, gene_likelihood={args.gene_likelihood}")
print(f"  Training: early_stopping={not args.no_early_stopping}")
print(f"  Downstream: resolution={args.resolution}, n_neighbors={args.n_neighbors}, min_dist={args.min_dist}")

# =============================================================================
# Step 1: Load and prepare data
# =============================================================================
print("\n" + "=" * 80)
print("[1/7] Loading data...")
print("=" * 80)

adata = ad.read_h5ad(input_file)
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

# Verify batch key exists
if args.batch_key not in adata.obs.columns:
    print(f"\nError: Batch key '{args.batch_key}' not found in adata.obs")
    print(f"Available columns: {list(adata.obs.columns)}")
    sys.exit(1)

n_batches = adata.obs[args.batch_key].nunique()
print(f"Batches: {n_batches} unique values in '{args.batch_key}'")
print(f"  Batch counts: {dict(adata.obs[args.batch_key].value_counts())}")

# Check for labels
run_scanvi = False
if args.labels_key and not args.skip_scanvi:
    if args.labels_key not in adata.obs.columns:
        print(f"\nWarning: Labels key '{args.labels_key}' not found, skipping scANVI")
    else:
        run_scanvi = True
        n_labels = adata.obs[args.labels_key].nunique()
        print(f"Cell types: {n_labels} unique values in '{args.labels_key}'")

# Auto-detect counts layer
layer = args.layer
if layer is None:
    if 'counts' in adata.layers:
        layer = 'counts'
        print(f"Using counts from adata.layers['counts']")
    elif adata.X is not None:
        # Check if X contains integers (likely counts)
        if hasattr(adata.X, 'toarray'):
            sample = adata.X[:100].toarray()
        else:
            sample = adata.X[:100]
        if np.allclose(sample, sample.astype(int)):
            layer = None
            print(f"Using counts from adata.X")
        else:
            print("\nWarning: adata.X may contain normalized data, not raw counts")
            print("scVI requires raw count data. Proceeding anyway...")
            layer = None
else:
    print(f"Using counts from adata.layers['{layer}']")

# =============================================================================
# Step 2: Setup AnnData for scvi-tools
# =============================================================================
print("\n" + "=" * 80)
print("[2/7] Setting up data for scvi-tools...")
print("=" * 80)

setup_anndata_scvi(
    adata,
    batch_key=args.batch_key,
    layer=layer,
    labels_key=args.labels_key if run_scanvi else None
)

# =============================================================================
# Step 3: Train scVI model
# =============================================================================
print("\n" + "=" * 80)
print("[3/7] Training scVI model...")
print("=" * 80)

scvi_model = train_scvi_model(
    adata,
    n_layers=args.n_layers,
    n_latent=args.n_latent,
    gene_likelihood=args.gene_likelihood,
    max_epochs=args.max_epochs_scvi,
    early_stopping=not args.no_early_stopping
)

# Extract scVI latent representation
SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = get_latent_representation(scvi_model, adata)
print(f"  Stored latent representation in adata.obsm['{SCVI_LATENT_KEY}']")

# Save scVI model
scvi_model_path = os.path.join(output_dir, 'scvi_model')
save_model(scvi_model, scvi_model_path)
print(f"  Saved model to: {scvi_model_path}")

# =============================================================================
# Step 4: Train scANVI model (if labels available)
# =============================================================================
SCANVI_LATENT_KEY = None
scanvi_model = None

if run_scanvi:
    print("\n" + "=" * 80)
    print("[4/7] Training scANVI model...")
    print("=" * 80)

    scanvi_model = train_scanvi_model(
        scvi_model,
        adata,
        labels_key=args.labels_key,
        unlabeled_category=args.unlabeled_category,
        max_epochs=args.max_epochs_scanvi
    )

    # Extract scANVI latent representation
    SCANVI_LATENT_KEY = "X_scANVI"
    adata.obsm[SCANVI_LATENT_KEY] = get_latent_representation(scanvi_model, adata)
    print(f"  Stored latent representation in adata.obsm['{SCANVI_LATENT_KEY}']")

    # Save scANVI model
    scanvi_model_path = os.path.join(output_dir, 'scanvi_model')
    save_model(scanvi_model, scanvi_model_path)
    print(f"  Saved model to: {scanvi_model_path}")
else:
    print("\n" + "=" * 80)
    print("[4/7] Skipping scANVI (no labels provided or --skip-scanvi)")
    print("=" * 80)

# =============================================================================
# Step 5: Downstream analysis (neighbors, UMAP, clustering)
# =============================================================================
print("\n" + "=" * 80)
print("[5/7] Computing neighbors, UMAP, and clustering...")
print("=" * 80)

# Process scVI embedding
print("\n  Processing scVI embedding...")
compute_neighbors_and_umap(
    adata,
    use_rep=SCVI_LATENT_KEY,
    n_neighbors=args.n_neighbors,
    min_dist=args.min_dist,
    key_added='scvi'
)
cluster_cells(adata, resolution=args.resolution, neighbors_key='scvi_neighbors', key_added='leiden_scvi')

# Process scANVI embedding if available
if SCANVI_LATENT_KEY:
    print("\n  Processing scANVI embedding...")
    compute_neighbors_and_umap(
        adata,
        use_rep=SCANVI_LATENT_KEY,
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist,
        key_added='scanvi'
    )
    cluster_cells(adata, resolution=args.resolution, neighbors_key='scanvi_neighbors', key_added='leiden_scanvi')

# =============================================================================
# Step 6: Calculate integration metrics
# =============================================================================
print("\n" + "=" * 80)
print("[6/7] Calculating integration metrics...")
print("=" * 80)

metrics_results = []

# scVI metrics
print("\n  Computing metrics for scVI...")
scvi_metrics = calculate_integration_metrics(
    adata,
    batch_key=args.batch_key,
    label_key=args.labels_key,
    embed_key=SCVI_LATENT_KEY,
    cluster_key='leiden_scvi'
)
scvi_metrics['method'] = 'scVI'
metrics_results.append(scvi_metrics)

# scANVI metrics
if SCANVI_LATENT_KEY:
    print("  Computing metrics for scANVI...")
    scanvi_metrics = calculate_integration_metrics(
        adata,
        batch_key=args.batch_key,
        label_key=args.labels_key,
        embed_key=SCANVI_LATENT_KEY,
        cluster_key='leiden_scanvi'
    )
    scanvi_metrics['method'] = 'scANVI'
    metrics_results.append(scanvi_metrics)

# Create metrics DataFrame
metrics_df = pd.DataFrame(metrics_results)
metrics_df = metrics_df.set_index('method')

# Save metrics
metrics_path = os.path.join(output_dir, 'integration_metrics.csv')
metrics_df.to_csv(metrics_path)
print(f"\n  Saved metrics to: {metrics_path}")

# Print metrics table
print("\n" + format_metrics_table(metrics_df))

# =============================================================================
# Step 7: Generate visualizations
# =============================================================================
print("\n" + "=" * 80)
print("[7/7] Generating visualizations...")
print("=" * 80)

# Set up plotting style
plt.style.use('default')
sc.settings.set_figure_params(dpi=100, frameon=False)

# Plot scVI UMAP
print("\n  Creating scVI visualization...")
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

sc.pl.embedding(adata, basis='X_umap_scvi', color=args.batch_key,
                title='scVI - Colored by Batch', ax=axes[0], show=False)
sc.pl.embedding(adata, basis='X_umap_scvi', color='leiden_scvi',
                title='scVI - Clusters', ax=axes[1], show=False)

plt.tight_layout()
scvi_plot_path = os.path.join(output_dir, 'scvi_latent_umap.png')
plt.savefig(scvi_plot_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: {scvi_plot_path}")

# Plot scANVI UMAP
if SCANVI_LATENT_KEY:
    print("  Creating scANVI visualization...")
    n_cols = 3 if args.labels_key else 2
    fig, axes = plt.subplots(1, n_cols, figsize=(7*n_cols, 6))

    sc.pl.embedding(adata, basis='X_umap_scanvi', color=args.batch_key,
                    title='scANVI - Colored by Batch', ax=axes[0], show=False)
    sc.pl.embedding(adata, basis='X_umap_scanvi', color='leiden_scanvi',
                    title='scANVI - Clusters', ax=axes[1], show=False)
    if args.labels_key:
        sc.pl.embedding(adata, basis='X_umap_scanvi', color=args.labels_key,
                        title='scANVI - Cell Types', ax=axes[2], show=False)

    plt.tight_layout()
    scanvi_plot_path = os.path.join(output_dir, 'scanvi_latent_umap.png')
    plt.savefig(scanvi_plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {scanvi_plot_path}")

# Create comparison figure
print("  Creating comparison visualization...")
if SCANVI_LATENT_KEY:
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    sc.pl.embedding(adata, basis='X_umap_scvi', color=args.batch_key,
                    title='scVI - Batch', ax=axes[0, 0], show=False)
    sc.pl.embedding(adata, basis='X_umap_scvi', color='leiden_scvi',
                    title='scVI - Clusters', ax=axes[0, 1], show=False)
    sc.pl.embedding(adata, basis='X_umap_scanvi', color=args.batch_key,
                    title='scANVI - Batch', ax=axes[1, 0], show=False)
    sc.pl.embedding(adata, basis='X_umap_scanvi', color='leiden_scanvi',
                    title='scANVI - Clusters', ax=axes[1, 1], show=False)
else:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    sc.pl.embedding(adata, basis='X_umap_scvi', color=args.batch_key,
                    title='scVI - Batch', ax=axes[0], show=False)
    sc.pl.embedding(adata, basis='X_umap_scvi', color='leiden_scvi',
                    title='scVI - Clusters', ax=axes[1], show=False)

plt.tight_layout()
comparison_path = os.path.join(output_dir, 'integration_comparison.png')
plt.savefig(comparison_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: {comparison_path}")

# =============================================================================
# Save integrated data
# =============================================================================
print("\n" + "=" * 80)
print("Saving integrated data...")
print("=" * 80)

output_h5ad = os.path.join(output_dir, f'{base_name}_integrated.h5ad')
adata.write(output_h5ad)
print(f"  Saved: {output_h5ad}")

# =============================================================================
# Summary
# =============================================================================
print("\n" + "=" * 80)
print("Integration Analysis Complete!")
print("=" * 80)

print(f"\nResults saved to: {output_dir}/")
print("\nGenerated files:")
print(f"  1. scvi_latent_umap.png - scVI UMAP visualization")
if SCANVI_LATENT_KEY:
    print(f"  2. scanvi_latent_umap.png - scANVI UMAP visualization")
print(f"  3. integration_comparison.png - Side-by-side comparison")
print(f"  4. integration_metrics.csv - Quantitative benchmarks")
print(f"  5. {base_name}_integrated.h5ad - Integrated dataset")
print(f"  6. scvi_model/ - Saved scVI model")
if SCANVI_LATENT_KEY:
    print(f"  7. scanvi_model/ - Saved scANVI model")

print("\nIntegrated embeddings stored in adata.obsm:")
print(f"  - '{SCVI_LATENT_KEY}' - scVI latent representation")
if SCANVI_LATENT_KEY:
    print(f"  - '{SCANVI_LATENT_KEY}' - scANVI latent representation")

print("\nNext steps:")
print("  - Examine visualizations to assess integration quality")
print("  - Use integrated embedding for downstream analysis")
print("  - Consider differential expression with scVI's built-in DE")
print("  - For cell annotation, use scANVI predictions or CellTypist")
