#!/usr/bin/env python3
"""
TotalVI CITE-Seq Analysis Pipeline

Complete workflow for analyzing CITE-Seq data using TotalVI from scvi-tools.
Performs joint RNA + protein dimensionality reduction, denoising, clustering,
and differential expression analysis.

Usage:
    python totalvi_analysis.py input.h5mu --batch-key batch
    python totalvi_analysis.py input.h5ad --protein-obsm-key protein_expression
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import matplotlib.pyplot as plt
import torch

# Import core functions
from totalvi_core import (
    setup_mudata_totalvi,
    setup_anndata_totalvi,
    train_totalvi_model,
    get_latent_representation,
    get_denoised_expression,
    get_foreground_probability,
    run_differential_expression,
    validate_citeseq_data,
    save_model,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="TotalVI CITE-Seq analysis pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Input/output
    parser.add_argument("input", help="Path to input file (.h5mu or .h5ad)")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory (default: <input_basename>_totalvi_results)",
    )

    # Data format options
    parser.add_argument(
        "--protein-obsm-key",
        default="protein_expression",
        help="For AnnData: obsm key containing protein data",
    )
    parser.add_argument(
        "--rna-mod-key",
        default="rna",
        help="For MuData: modality name for RNA",
    )
    parser.add_argument(
        "--protein-mod-key",
        default="prot",
        help="For MuData: modality name for protein",
    )

    # Batch correction
    parser.add_argument(
        "--batch-key",
        default=None,
        help="Column in obs for batch correction",
    )

    # Model architecture
    parser.add_argument(
        "--n-latent",
        type=int,
        default=20,
        help="Latent space dimensions",
    )
    parser.add_argument(
        "--n-layers-encoder",
        type=int,
        default=2,
        help="Number of encoder layers",
    )
    parser.add_argument(
        "--n-layers-decoder",
        type=int,
        default=1,
        help="Number of decoder layers",
    )

    # Training
    parser.add_argument(
        "--max-epochs",
        type=int,
        default=400,
        help="Maximum training epochs",
    )
    parser.add_argument(
        "--early-stopping",
        action="store_true",
        default=True,
        help="Enable early stopping",
    )
    parser.add_argument(
        "--no-early-stopping",
        action="store_false",
        dest="early_stopping",
        help="Disable early stopping",
    )

    # Downstream analysis
    parser.add_argument(
        "--resolution",
        type=float,
        default=0.5,
        help="Leiden clustering resolution",
    )
    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=15,
        help="Number of neighbors for kNN graph",
    )
    parser.add_argument(
        "--min-dist",
        type=float,
        default=0.3,
        help="UMAP minimum distance parameter",
    )

    # Differential expression
    parser.add_argument(
        "--de-clusters",
        nargs=2,
        default=None,
        help="Two cluster IDs to compare for DE (e.g., --de-clusters 0 1)",
    )
    parser.add_argument(
        "--de-delta",
        type=float,
        default=0.5,
        help="Effect size threshold for DE",
    )

    # Other
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed",
    )
    parser.add_argument(
        "--skip-de",
        action="store_true",
        help="Skip differential expression analysis",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    # Set seeds
    scvi.settings.seed = args.seed
    torch.set_float32_matmul_precision("high")
    sc.set_figure_params(figsize=(6, 6), frameon=False)

    # Setup output directory
    if args.output_dir is None:
        input_basename = Path(args.input).stem
        args.output_dir = f"{input_basename}_totalvi_results"
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Output directory: {args.output_dir}")

    # Load data
    print(f"\n{'='*60}")
    print("Loading data...")
    print(f"{'='*60}")

    input_path = Path(args.input)
    is_mudata = input_path.suffix in [".h5mu"]

    if is_mudata:
        import muon
        mdata = muon.read_h5mu(args.input)
        print(f"Loaded MuData: {list(mdata.mod.keys())}")

        rna = mdata.mod[args.rna_mod_key]
        protein = mdata.mod[args.protein_mod_key]
        data = mdata
    else:
        adata = sc.read_h5ad(args.input)
        print(f"Loaded AnnData: {adata.shape}")

        if args.protein_obsm_key not in adata.obsm:
            raise ValueError(
                f"Protein data not found in obsm['{args.protein_obsm_key}']. "
                "Specify --protein-obsm-key or use MuData format."
            )

        rna = adata
        protein = None  # Protein in obsm
        data = adata

    print(f"RNA: {rna.n_obs} cells × {rna.n_vars} genes")
    if protein is not None:
        print(f"Protein: {protein.n_obs} cells × {protein.n_vars} proteins")

    # Validate data
    print(f"\n{'='*60}")
    print("Validating data...")
    print(f"{'='*60}")

    if is_mudata:
        is_valid = validate_citeseq_data(mdata, args.rna_mod_key, args.protein_mod_key)
    else:
        is_valid = True  # Basic validation for AnnData
        print("AnnData format - skipping detailed validation")

    # Prepare data
    print(f"\n{'='*60}")
    print("Preparing data for TotalVI...")
    print(f"{'='*60}")

    # Convert sparse to dense if needed
    if protein is not None and hasattr(protein.X, 'toarray'):
        print("Converting protein matrix to dense...")
        protein.X = protein.X.toarray()

    if hasattr(rna.X, 'toarray'):
        print("Converting RNA matrix to dense...")
        rna.X = rna.X.toarray()

    # Ensure counts layer
    if "counts" not in rna.layers:
        print("Creating counts layer...")
        rna.layers["counts"] = rna.X.copy()

    # Setup for TotalVI
    print(f"\n{'='*60}")
    print("Setting up TotalVI...")
    print(f"{'='*60}")

    if is_mudata:
        setup_mudata_totalvi(
            mdata,
            rna_layer="counts",
            batch_key=args.batch_key,
            rna_mod_key=args.rna_mod_key,
            protein_mod_key=args.protein_mod_key,
        )
    else:
        setup_anndata_totalvi(
            adata,
            layer="counts",
            batch_key=args.batch_key,
            protein_obsm_key=args.protein_obsm_key,
        )

    # Train model
    print(f"\n{'='*60}")
    print("Training TotalVI model...")
    print(f"{'='*60}")

    model = train_totalvi_model(
        data,
        n_latent=args.n_latent,
        n_layers_encoder=args.n_layers_encoder,
        n_layers_decoder=args.n_layers_decoder,
        max_epochs=args.max_epochs,
        early_stopping=args.early_stopping,
    )

    # Plot training history
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    axes[0].plot(model.history["elbo_train"].values, label="Train")
    if "elbo_validation" in model.history:
        axes[0].plot(model.history["elbo_validation"].values, label="Validation")
    axes[0].set_xlabel("Epoch")
    axes[0].set_ylabel("ELBO")
    axes[0].legend()
    axes[0].set_title("TotalVI Training")

    axes[1].plot(model.history["reconstruction_loss_train"].values)
    axes[1].set_xlabel("Epoch")
    axes[1].set_ylabel("Reconstruction Loss")
    axes[1].set_title("Reconstruction Loss")

    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "training_history.png"), dpi=150)
    plt.close()
    print(f"Saved: training_history.png")

    # Extract latent representation
    print(f"\n{'='*60}")
    print("Extracting latent representation...")
    print(f"{'='*60}")

    latent = get_latent_representation(model)
    rna.obsm["X_totalVI"] = latent
    print(f"Latent shape: {latent.shape}")

    # Get denoised expression
    print(f"\n{'='*60}")
    print("Computing denoised expression...")
    print(f"{'='*60}")

    rna_denoised, protein_denoised = get_denoised_expression(model)
    rna.layers["denoised_rna"] = rna_denoised

    if protein is not None:
        protein.layers["denoised_protein"] = protein_denoised
        protein_fg_prob = get_foreground_probability(model)
        protein.layers["foreground_prob"] = protein_fg_prob
        print(f"Foreground probability range: {protein_fg_prob.min():.2f} - {protein_fg_prob.max():.2f}")

    # Downstream analysis
    print(f"\n{'='*60}")
    print("Running downstream analysis...")
    print(f"{'='*60}")

    sc.pp.neighbors(rna, use_rep="X_totalVI", n_neighbors=args.n_neighbors)
    sc.tl.umap(rna, min_dist=args.min_dist)
    sc.tl.leiden(rna, key_added="leiden_totalVI", resolution=args.resolution)

    print(f"Found {rna.obs['leiden_totalVI'].nunique()} clusters")

    # Visualization
    print(f"\n{'='*60}")
    print("Generating visualizations...")
    print(f"{'='*60}")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    sc.pl.umap(rna, color="leiden_totalVI", ax=axes[0], show=False, title="Clusters")

    if args.batch_key and args.batch_key in rna.obs.columns:
        sc.pl.umap(rna, color=args.batch_key, ax=axes[1], show=False, title="Batch")
    else:
        axes[1].set_visible(False)

    # Show first protein marker
    if protein is not None and protein.n_vars > 0:
        first_protein = protein.var_names[0]
        idx = 0
        rna.obs["_protein_viz"] = np.log1p(protein_denoised[:, idx])
        sc.pl.umap(rna, color="_protein_viz", ax=axes[2], show=False,
                   title=f"{first_protein} (denoised)", vmax="p99")
        del rna.obs["_protein_viz"]
    else:
        axes[2].set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "totalvi_umap.png"), dpi=150)
    plt.close()
    print(f"Saved: totalvi_umap.png")

    # Differential expression
    if not args.skip_de:
        print(f"\n{'='*60}")
        print("Running differential expression...")
        print(f"{'='*60}")

        if args.de_clusters:
            group1, group2 = args.de_clusters
        else:
            # Compare largest two clusters
            cluster_counts = rna.obs["leiden_totalVI"].value_counts()
            group1 = str(cluster_counts.index[0])
            group2 = str(cluster_counts.index[1])
            print(f"Comparing clusters {group1} vs {group2}")

        # Need to specify groupby with modality prefix for MuData
        if is_mudata:
            groupby = f"{args.rna_mod_key}:leiden_totalVI"
        else:
            groupby = "leiden_totalVI"

        de_results = run_differential_expression(
            model,
            groupby=groupby,
            group1=group1,
            group2=group2,
            delta=args.de_delta,
        )

        de_results.to_csv(os.path.join(args.output_dir, "differential_expression.csv"))
        print(f"Saved: differential_expression.csv")

        # Summary
        sig_de = de_results[de_results["is_de_fdr"]].shape[0]
        print(f"Significant DE features (FDR < 0.05): {sig_de}")

    # Save results
    print(f"\n{'='*60}")
    print("Saving results...")
    print(f"{'='*60}")

    model_dir = os.path.join(args.output_dir, "totalvi_model")
    save_model(model, model_dir)
    print(f"Saved: totalvi_model/")

    if is_mudata:
        import muon
        output_path = os.path.join(args.output_dir, f"{Path(args.input).stem}_totalvi.h5mu")
        mdata.write_h5mu(output_path)
    else:
        output_path = os.path.join(args.output_dir, f"{Path(args.input).stem}_totalvi.h5ad")
        adata.write_h5ad(output_path)

    print(f"Saved: {os.path.basename(output_path)}")

    print(f"\n{'='*60}")
    print("TotalVI analysis complete!")
    print(f"{'='*60}")
    print(f"Results saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
