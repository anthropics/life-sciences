#!/usr/bin/env python3
"""
PeakVI scATAC-seq Analysis Pipeline

Complete workflow for analyzing single-cell ATAC-seq data using PeakVI.

Usage:
    python peakvi_analysis.py input.h5ad
    python peakvi_analysis.py input.h5ad --batch-key batch --resolution 0.3
"""

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import matplotlib.pyplot as plt
import torch

from peakvi_core import (
    setup_anndata_peakvi,
    train_peakvi_model,
    get_latent_representation,
    run_differential_accessibility,
    get_cluster_markers,
    filter_peaks_by_detection,
    calculate_qc_metrics,
    save_model,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="PeakVI scATAC-seq analysis pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("input", help="Path to input .h5ad file")
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--batch-key", default=None)

    # Filtering
    parser.add_argument("--min-detection", type=float, default=0.05)

    # Model
    parser.add_argument("--n-latent", type=int, default=13)
    parser.add_argument("--n-hidden", type=int, default=128)
    parser.add_argument("--max-epochs", type=int, default=500)

    # Downstream
    parser.add_argument("--resolution", type=float, default=0.2)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--min-dist", type=float, default=0.2)

    # DA
    parser.add_argument("--skip-da", action="store_true")
    parser.add_argument("--n-top-markers", type=int, default=50)

    parser.add_argument("--seed", type=int, default=0)

    return parser.parse_args()


def main():
    args = parse_args()

    scvi.settings.seed = args.seed
    torch.set_float32_matmul_precision("high")
    sc.set_figure_params(figsize=(6, 6), frameon=False)

    if args.output_dir is None:
        args.output_dir = f"{Path(args.input).stem}_peakvi_results"
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Output: {args.output_dir}")

    # Load
    print(f"\n{'='*60}\nLoading data...\n{'='*60}")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} peaks")

    # QC
    calculate_qc_metrics(adata)

    # Filter
    print(f"\n{'='*60}\nFiltering peaks...\n{'='*60}")
    adata = filter_peaks_by_detection(adata, min_detection=args.min_detection)

    # Setup
    print(f"\n{'='*60}\nSetting up PeakVI...\n{'='*60}")
    setup_anndata_peakvi(adata, batch_key=args.batch_key)

    # Train
    print(f"\n{'='*60}\nTraining PeakVI...\n{'='*60}")
    model = train_peakvi_model(
        adata,
        n_latent=args.n_latent,
        n_hidden=args.n_hidden,
        max_epochs=args.max_epochs,
    )

    # Plot training
    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.plot(model.history["elbo_train"].values, label="Train")
    if "elbo_validation" in model.history:
        plt.plot(model.history["elbo_validation"].values, label="Validation")
    plt.xlabel("Epoch")
    plt.ylabel("ELBO")
    plt.legend()
    plt.title("Training")

    plt.subplot(1, 2, 2)
    plt.plot(model.history["reconstruction_loss_train"].values)
    plt.xlabel("Epoch")
    plt.ylabel("Reconstruction Loss")
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "training_history.png"), dpi=150)
    plt.close()

    # Latent
    print(f"\n{'='*60}\nExtracting latent...\n{'='*60}")
    adata.obsm["X_peakvi"] = get_latent_representation(model)

    # Downstream
    print(f"\n{'='*60}\nClustering...\n{'='*60}")
    sc.pp.neighbors(adata, use_rep="X_peakvi", n_neighbors=args.n_neighbors)
    sc.tl.umap(adata, min_dist=args.min_dist)
    sc.tl.leiden(adata, key_added="clusters_peakvi", resolution=args.resolution)
    print(f"Found {adata.obs['clusters_peakvi'].nunique()} clusters")

    # Visualize
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sc.pl.umap(adata, color="clusters_peakvi", ax=axes[0], show=False, title="Clusters")
    if args.batch_key and args.batch_key in adata.obs.columns:
        sc.pl.umap(adata, color=args.batch_key, ax=axes[1], show=False, title="Batch")
    else:
        sc.pl.umap(adata, color="n_peaks", ax=axes[1], show=False, title="Peaks/Cell")
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "peakvi_umap.png"), dpi=150)
    plt.close()

    # DA
    if not args.skip_da:
        print(f"\n{'='*60}\nFinding marker peaks...\n{'='*60}")
        markers = get_cluster_markers(
            model, adata, "clusters_peakvi", n_top=args.n_top_markers
        )

        all_markers = pd.concat(
            [df.assign(cluster=k) for k, df in markers.items()], ignore_index=False
        )
        all_markers.to_csv(os.path.join(args.output_dir, "marker_peaks.csv"))
        print(f"Saved: marker_peaks.csv")

    # Save
    print(f"\n{'='*60}\nSaving results...\n{'='*60}")
    save_model(model, os.path.join(args.output_dir, "peakvi_model"))
    adata.write_h5ad(os.path.join(args.output_dir, f"{Path(args.input).stem}_peakvi.h5ad"))

    print(f"\n{'='*60}\nPeakVI analysis complete!\n{'='*60}")
    print(f"Results: {args.output_dir}")


if __name__ == "__main__":
    main()
