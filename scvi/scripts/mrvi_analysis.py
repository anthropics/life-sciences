#!/usr/bin/env python3
"""
Multi-Sample Single-Cell Analysis using MrVI

This script provides a complete workflow for analyzing multi-sample scRNA-seq
data using MrVI (Multi-resolution Variational Inference) from scvi-tools.

Reference: https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/MrVI_tutorial.html
"""

import anndata as ad
import scanpy as sc
import scvi
from scvi.external import MRVI
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
from mrvi_core import (
    setup_anndata_mrvi,
    train_mrvi_model,
    get_latent_representation,
    compute_sample_distances,
    run_differential_expression,
    run_differential_abundance,
    save_model
)

warnings.filterwarnings('ignore')

print("=" * 80)
print("Multi-Sample Single-Cell Analysis with MrVI")
print("=" * 80)

# Default parameters
DEFAULT_N_LATENT = 30
DEFAULT_MAX_EPOCHS = 400
DEFAULT_RESOLUTION = 1.0
DEFAULT_MIN_DIST = 0.3
DEFAULT_N_TOP_GENES = 10000

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description='Multi-Sample Single-Cell Analysis using MrVI',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
Examples:
  python3 mrvi_analysis.py data.h5ad --sample-key patient_id
  python3 mrvi_analysis.py data.h5ad --sample-key donor --batch-key site
  python3 mrvi_analysis.py data.h5ad --sample-key patient --sample-cov-keys condition treatment
  python3 mrvi_analysis.py data.h5ad --sample-key donor --cell-type-key cell_type
    """
)

# Required arguments
parser.add_argument('input_file', help='Input .h5ad file with raw counts')
parser.add_argument('--sample-key', required=True,
                    help='Column in adata.obs containing sample identifiers')

# Optional: batch and cell type
parser.add_argument('--batch-key', type=str, default=None,
                    help='Column for nuisance covariates (e.g., sequencing batch)')
parser.add_argument('--cell-type-key', type=str, default=None,
                    help='Column for cell type grouping in distance analysis')

# Output
parser.add_argument('--output-dir', type=str, help='Output directory')

# Model architecture
parser.add_argument('--n-latent', type=int, default=DEFAULT_N_LATENT,
                    help=f'Latent space dimensions (default: {DEFAULT_N_LATENT})')

# Training
parser.add_argument('--max-epochs', type=int, default=DEFAULT_MAX_EPOCHS,
                    help=f'Maximum training epochs (default: {DEFAULT_MAX_EPOCHS})')
parser.add_argument('--no-early-stopping', action='store_true',
                    help='Disable early stopping')

# Downstream analysis
parser.add_argument('--resolution', type=float, default=DEFAULT_RESOLUTION,
                    help=f'Leiden clustering resolution (default: {DEFAULT_RESOLUTION})')
parser.add_argument('--min-dist', type=float, default=DEFAULT_MIN_DIST,
                    help=f'UMAP min_dist parameter (default: {DEFAULT_MIN_DIST})')

# DE/DA analysis
parser.add_argument('--sample-cov-keys', type=str, nargs='+', default=None,
                    help='Sample-level covariates for DE/DA analysis')

# Feature selection
parser.add_argument('--n-top-genes', type=int, default=DEFAULT_N_TOP_GENES,
                    help=f'Number of HVGs if not pre-selected (default: {DEFAULT_N_TOP_GENES})')
parser.add_argument('--skip-hvg', action='store_true',
                    help='Skip HVG selection (assume already done)')

# Advanced
parser.add_argument('--seed', type=int, default=0,
                    help='Random seed for reproducibility (default: 0)')

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
    output_dir = f"{base_name}_mrvi_results"

os.makedirs(output_dir, exist_ok=True)
print(f"\nOutput directory: {output_dir}")

# Display parameters
print(f"\nParameters:")
print(f"  Sample key: {args.sample_key}")
print(f"  Batch key: {args.batch_key or 'None'}")
print(f"  Cell type key: {args.cell_type_key or 'None'}")
print(f"  Model: n_latent={args.n_latent}")
print(f"  Training: max_epochs={args.max_epochs}, early_stopping={not args.no_early_stopping}")
if args.sample_cov_keys:
    print(f"  DE/DA covariates: {args.sample_cov_keys}")

# =============================================================================
# Step 1: Load and prepare data
# =============================================================================
print("\n" + "=" * 80)
print("[1/7] Loading data...")
print("=" * 80)

adata = ad.read_h5ad(input_file)
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

# Verify sample key exists
if args.sample_key not in adata.obs.columns:
    print(f"\nError: Sample key '{args.sample_key}' not found in adata.obs")
    print(f"Available columns: {list(adata.obs.columns)}")
    sys.exit(1)

n_samples = adata.obs[args.sample_key].nunique()
print(f"Samples: {n_samples} unique values in '{args.sample_key}'")

# Verify batch key if provided
if args.batch_key and args.batch_key not in adata.obs.columns:
    print(f"\nWarning: Batch key '{args.batch_key}' not found, ignoring")
    args.batch_key = None

# Verify cell type key if provided
if args.cell_type_key and args.cell_type_key not in adata.obs.columns:
    print(f"\nWarning: Cell type key '{args.cell_type_key}' not found, ignoring")
    args.cell_type_key = None

# =============================================================================
# Step 2: Feature selection
# =============================================================================
print("\n" + "=" * 80)
print("[2/7] Feature selection...")
print("=" * 80)

if not args.skip_hvg:
    if 'highly_variable' not in adata.var.columns:
        print(f"  Selecting {args.n_top_genes} highly variable genes...")
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=args.n_top_genes,
            flavor='seurat_v3',
            subset=True
        )
        print(f"  Selected {adata.n_vars} genes")
    else:
        n_hvg = adata.var['highly_variable'].sum()
        print(f"  Using {n_hvg} pre-selected highly variable genes")
        if not adata.var['highly_variable'].all():
            adata = adata[:, adata.var['highly_variable']].copy()
else:
    print(f"  Skipping HVG selection, using all {adata.n_vars} genes")

# =============================================================================
# Step 3: Setup and train MrVI
# =============================================================================
print("\n" + "=" * 80)
print("[3/7] Setting up MrVI...")
print("=" * 80)

setup_anndata_mrvi(adata, sample_key=args.sample_key, batch_key=args.batch_key)

print("\n" + "=" * 80)
print("[4/7] Training MrVI model...")
print("=" * 80)

model = train_mrvi_model(
    adata,
    n_latent=args.n_latent,
    max_epochs=args.max_epochs,
    early_stopping=not args.no_early_stopping
)

# Extract latent representation
MRVI_LATENT_KEY = "X_mrvi_u"
adata.obsm[MRVI_LATENT_KEY] = get_latent_representation(model, adata)
print(f"  Stored latent representation in adata.obsm['{MRVI_LATENT_KEY}']")

# Save model
model_path = os.path.join(output_dir, 'mrvi_model')
save_model(model, model_path)
print(f"  Saved model to: {model_path}")

# =============================================================================
# Step 5: Downstream analysis
# =============================================================================
print("\n" + "=" * 80)
print("[5/7] Computing neighbors, UMAP, and clustering...")
print("=" * 80)

# Compute neighbors and UMAP
sc.pp.neighbors(adata, use_rep=MRVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=args.min_dist)
sc.tl.leiden(adata, resolution=args.resolution)

n_clusters = adata.obs['leiden'].nunique()
print(f"  Found {n_clusters} clusters")

# =============================================================================
# Step 6: Sample distance analysis
# =============================================================================
print("\n" + "=" * 80)
print("[6/7] Computing sample distances...")
print("=" * 80)

distance_dir = os.path.join(output_dir, 'sample_distances')
os.makedirs(distance_dir, exist_ok=True)

if args.cell_type_key:
    print(f"  Computing distances grouped by '{args.cell_type_key}'...")
    distances = compute_sample_distances(model, groupby=args.cell_type_key)

    # Save distance matrices
    cell_types = adata.obs[args.cell_type_key].unique()
    for ct in cell_types:
        ct_safe = ct.replace('/', '_').replace(' ', '_')
        try:
            ct_dist = distances.loc[{f"{args.cell_type_key}_name": ct}]
            dist_df = pd.DataFrame(
                ct_dist.values,
                index=ct_dist.coords['sample_x'].values,
                columns=ct_dist.coords['sample_y'].values
            )
            dist_df.to_csv(os.path.join(distance_dir, f'{ct_safe}_distances.csv'))
        except Exception as e:
            print(f"    Warning: Could not save distances for {ct}: {e}")

    print(f"  Saved distance matrices to: {distance_dir}")

    # Create heatmaps
    print("  Creating distance heatmaps...")
    n_types = min(len(cell_types), 6)  # Limit to 6 for visualization
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for i, ct in enumerate(cell_types[:n_types]):
        try:
            ct_dist = distances.loc[{f"{args.cell_type_key}_name": ct}]
            dist_matrix = ct_dist.values
            sns.heatmap(dist_matrix, ax=axes[i], cmap='viridis',
                       xticklabels=False, yticklabels=False)
            axes[i].set_title(ct[:20])
        except:
            axes[i].set_visible(False)

    for i in range(n_types, 6):
        axes[i].set_visible(False)

    plt.tight_layout()
    heatmap_path = os.path.join(output_dir, 'sample_distance_heatmaps.png')
    plt.savefig(heatmap_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {heatmap_path}")
else:
    print("  Computing overall sample distances...")
    distances = compute_sample_distances(model, groupby=None)

    # Save overall distance matrix
    dist_df = pd.DataFrame(
        distances.values,
        index=distances.coords['sample_x'].values,
        columns=distances.coords['sample_y'].values
    )
    dist_df.to_csv(os.path.join(distance_dir, 'overall_distances.csv'))
    print(f"  Saved distances to: {distance_dir}")

# =============================================================================
# Step 7: DE and DA analysis (if covariates provided)
# =============================================================================
de_results = None
da_results = None

if args.sample_cov_keys:
    print("\n" + "=" * 80)
    print("[7/7] Differential expression and abundance analysis...")
    print("=" * 80)

    # Verify covariates exist in sample info
    valid_covs = []
    for cov in args.sample_cov_keys:
        if cov in model.sample_info.columns:
            valid_covs.append(cov)
        else:
            print(f"  Warning: Covariate '{cov}' not found in sample info, skipping")

    if valid_covs:
        # Differential Expression
        print(f"\n  Running differential expression for: {valid_covs}")
        de_dir = os.path.join(output_dir, 'de_results')
        os.makedirs(de_dir, exist_ok=True)

        try:
            de_results = run_differential_expression(model, valid_covs, store_lfc=True)

            # Save effect sizes
            for cov in de_results.effect_size.coords['covariate'].values:
                effect = de_results.effect_size.sel(covariate=cov).values
                cov_safe = cov.replace('/', '_').replace(' ', '_')
                np.save(os.path.join(de_dir, f'{cov_safe}_effect_sizes.npy'), effect)

            # Add first covariate effect to adata for visualization
            first_cov = de_results.effect_size.coords['covariate'].values[0]
            adata.obs['DE_effect_size'] = de_results.effect_size.sel(covariate=first_cov).values

            print(f"  Saved DE results to: {de_dir}")
        except Exception as e:
            print(f"  Warning: DE analysis failed: {e}")

        # Differential Abundance
        print(f"\n  Running differential abundance for: {valid_covs}")
        da_dir = os.path.join(output_dir, 'da_results')
        os.makedirs(da_dir, exist_ok=True)

        try:
            da_results = run_differential_abundance(model, valid_covs)

            # Save log probabilities
            for cov in valid_covs:
                log_probs_key = f"{cov}_log_probs"
                if hasattr(da_results, log_probs_key):
                    log_probs = getattr(da_results, log_probs_key)
                    cov_safe = cov.replace('/', '_').replace(' ', '_')

                    # Compute log ratio for first two categories
                    categories = log_probs.coords[cov].values
                    if len(categories) >= 2:
                        lfc = log_probs.loc[{cov: categories[1]}] - log_probs.loc[{cov: categories[0]}]
                        adata.obs['DA_lfc'] = lfc.values
                        np.save(os.path.join(da_dir, f'{cov_safe}_log_ratio.npy'), lfc.values)

            print(f"  Saved DA results to: {da_dir}")
        except Exception as e:
            print(f"  Warning: DA analysis failed: {e}")
else:
    print("\n" + "=" * 80)
    print("[7/7] Skipping DE/DA analysis (no covariates provided)")
    print("=" * 80)

# =============================================================================
# Visualizations
# =============================================================================
print("\n" + "=" * 80)
print("Generating visualizations...")
print("=" * 80)

plt.style.use('default')
sc.settings.set_figure_params(dpi=100, frameon=False)

# Main UMAP
print("  Creating UMAP visualization...")
color_keys = [args.sample_key, 'leiden']
if args.cell_type_key:
    color_keys.append(args.cell_type_key)

n_cols = len(color_keys)
fig, axes = plt.subplots(1, n_cols, figsize=(6*n_cols, 5))
if n_cols == 1:
    axes = [axes]

for ax, color in zip(axes, color_keys):
    sc.pl.umap(adata, color=color, ax=ax, show=False, title=color)

plt.tight_layout()
umap_path = os.path.join(output_dir, 'mrvi_umap.png')
plt.savefig(umap_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: {umap_path}")

# DE effect size UMAP
if 'DE_effect_size' in adata.obs.columns:
    print("  Creating DE effect size UMAP...")
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color='DE_effect_size', cmap='viridis', ax=ax, show=False,
               title='Differential Expression Effect Size')
    plt.tight_layout()
    de_path = os.path.join(output_dir, 'de_effect_umap.png')
    plt.savefig(de_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {de_path}")

# DA log fold change UMAP
if 'DA_lfc' in adata.obs.columns:
    print("  Creating DA log fold change UMAP...")
    fig, ax = plt.subplots(figsize=(8, 6))
    vmax = np.percentile(np.abs(adata.obs['DA_lfc']), 95)
    sc.pl.umap(adata, color='DA_lfc', cmap='coolwarm', ax=ax, show=False,
               vmin=-vmax, vmax=vmax, title='Differential Abundance (log ratio)')
    plt.tight_layout()
    da_path = os.path.join(output_dir, 'da_lfc_umap.png')
    plt.savefig(da_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {da_path}")

# =============================================================================
# Save processed data
# =============================================================================
print("\n" + "=" * 80)
print("Saving processed data...")
print("=" * 80)

output_h5ad = os.path.join(output_dir, f'{base_name}_mrvi.h5ad')
adata.write(output_h5ad)
print(f"  Saved: {output_h5ad}")

# =============================================================================
# Summary
# =============================================================================
print("\n" + "=" * 80)
print("MrVI Analysis Complete!")
print("=" * 80)

print(f"\nResults saved to: {output_dir}/")
print("\nGenerated files:")
print(f"  1. mrvi_umap.png - UMAP visualization")
print(f"  2. sample_distances/ - Sample distance matrices")
if args.cell_type_key:
    print(f"  3. sample_distance_heatmaps.png - Distance heatmaps by cell type")
if de_results is not None:
    print(f"  4. de_results/ - Differential expression results")
    print(f"  5. de_effect_umap.png - DE effect sizes on UMAP")
if da_results is not None:
    print(f"  6. da_results/ - Differential abundance results")
    print(f"  7. da_lfc_umap.png - DA log ratios on UMAP")
print(f"  8. {base_name}_mrvi.h5ad - Processed dataset")
print(f"  9. mrvi_model/ - Saved MrVI model")

print(f"\nLatent embedding stored in adata.obsm['{MRVI_LATENT_KEY}']")

print("\nNext steps:")
print("  - Examine sample distances for population-specific relationships")
print("  - Perform hierarchical clustering of samples")
print("  - Run gene set enrichment on DE genes")
print("  - Integrate results with clinical metadata")
