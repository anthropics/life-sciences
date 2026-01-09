#!/usr/bin/env python3
"""
DestVI Spatial Deconvolution Pipeline

Complete workflow for deconvolving spatial transcriptomics using DestVI.

Usage:
    python destvi_analysis.py --reference ref.h5ad --spatial spatial.h5ad --labels-key cell_type
"""

import argparse
import os
from pathlib import Path

import numpy as np
import scanpy as sc
import scvi
import matplotlib.pyplot as plt
import torch

from destvi_core import (
    setup_reference_data,
    train_reference_model,
    setup_spatial_data,
    train_spatial_model,
    get_proportions,
    get_gamma_values,
    filter_shared_genes,
    select_hvg_genes,
    add_proportions_to_obs,
    add_gamma_to_obsm,
    save_reference_model,
    save_spatial_model,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="DestVI spatial deconvolution pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--reference", required=True, help="scRNA-seq reference .h5ad")
    parser.add_argument("--spatial", required=True, help="Spatial .h5ad")
    parser.add_argument("--labels-key", required=True, help="Cell type column in reference")
    parser.add_argument("--output-dir", default=None)

    # Gene selection
    parser.add_argument("--n-hvg", type=int, default=3000)

    # Training
    parser.add_argument("--ref-epochs", type=int, default=300)
    parser.add_argument("--spatial-epochs", type=int, default=2500)
    parser.add_argument("--weight-obs", action="store_true")

    parser.add_argument("--seed", type=int, default=0)

    return parser.parse_args()


def main():
    args = parse_args()

    scvi.settings.seed = args.seed
    torch.set_float32_matmul_precision("high")
    sc.set_figure_params(figsize=(6, 6), frameon=False)

    if args.output_dir is None:
        args.output_dir = f"{Path(args.spatial).stem}_destvi_results"
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Output: {args.output_dir}")

    # Load reference
    print(f"\n{'='*60}\nLoading reference...\n{'='*60}")
    sc_adata = sc.read_h5ad(args.reference)
    print(f"Reference: {sc_adata.n_obs} cells × {sc_adata.n_vars} genes")
    print(f"Cell types: {sc_adata.obs[args.labels_key].value_counts().to_dict()}")

    if "counts" not in sc_adata.layers:
        sc_adata.layers["counts"] = sc_adata.X.copy()

    # Load spatial
    print(f"\n{'='*60}\nLoading spatial...\n{'='*60}")
    st_adata = sc.read_h5ad(args.spatial)
    print(f"Spatial: {st_adata.n_obs} spots × {st_adata.n_vars} genes")

    if "counts" not in st_adata.layers:
        st_adata.layers["counts"] = st_adata.X.copy()

    # Filter shared genes
    print(f"\n{'='*60}\nFiltering genes...\n{'='*60}")
    sc_adata, st_adata = filter_shared_genes(sc_adata, st_adata)

    # HVG
    hvg = select_hvg_genes(sc_adata, n_top_genes=args.n_hvg, layer="counts")
    sc_adata = sc_adata[:, hvg].copy()
    st_adata = st_adata[:, hvg].copy()

    # Train reference
    print(f"\n{'='*60}\nTraining reference model...\n{'='*60}")
    setup_reference_data(sc_adata, layer="counts", labels_key=args.labels_key)
    ref_model = train_reference_model(
        sc_adata, weight_obs=args.weight_obs, max_epochs=args.ref_epochs
    )

    plt.figure(figsize=(8, 4))
    plt.plot(ref_model.history["elbo_train"].values)
    plt.xlabel("Epoch")
    plt.ylabel("ELBO")
    plt.title("Reference Model Training")
    plt.savefig(os.path.join(args.output_dir, "reference_training.png"), dpi=150)
    plt.close()

    # Train spatial
    print(f"\n{'='*60}\nTraining spatial model...\n{'='*60}")
    setup_spatial_data(st_adata, layer="counts")
    spatial_model = train_spatial_model(
        st_adata, ref_model, max_epochs=args.spatial_epochs
    )

    plt.figure(figsize=(8, 4))
    plt.plot(spatial_model.history["elbo_train"].values)
    plt.xlabel("Epoch")
    plt.ylabel("ELBO")
    plt.title("Spatial Model Training")
    plt.savefig(os.path.join(args.output_dir, "spatial_training.png"), dpi=150)
    plt.close()

    # Get results
    print(f"\n{'='*60}\nExtracting results...\n{'='*60}")
    proportions = get_proportions(spatial_model)
    add_proportions_to_obs(st_adata, proportions)

    gamma = get_gamma_values(spatial_model)
    add_gamma_to_obsm(st_adata, gamma)

    # Visualize proportions
    print(f"\n{'='*60}\nVisualizing...\n{'='*60}")
    cell_types = proportions.columns.tolist()
    n_types = min(6, len(cell_types))

    if "spatial" in st_adata.uns:
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()
        for i, ct in enumerate(cell_types[:n_types]):
            sc.pl.spatial(st_adata, color=f"prop_{ct}", ax=axes[i], show=False,
                         title=ct, vmax=1.0)
        plt.tight_layout()
        plt.savefig(os.path.join(args.output_dir, "proportions_spatial.png"), dpi=150)
        plt.close()
    else:
        print("No spatial coordinates - skipping spatial plot")

    # Save
    print(f"\n{'='*60}\nSaving results...\n{'='*60}")
    save_reference_model(ref_model, os.path.join(args.output_dir, "reference_model"))
    save_spatial_model(spatial_model, os.path.join(args.output_dir, "spatial_model"))

    proportions.to_csv(os.path.join(args.output_dir, "cell_type_proportions.csv"))
    st_adata.write_h5ad(os.path.join(args.output_dir, f"{Path(args.spatial).stem}_deconvolved.h5ad"))

    print(f"\n{'='*60}\nDestVI analysis complete!\n{'='*60}")
    print(f"Results: {args.output_dir}")


if __name__ == "__main__":
    main()
