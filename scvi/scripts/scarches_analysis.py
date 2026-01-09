#!/usr/bin/env python3
"""
scArches Reference Mapping Pipeline

Complete workflow for mapping query data to reference atlases using scArches.

Usage:
    python scarches_analysis.py --reference ref.h5ad --query query.h5ad --labels-key cell_type
    python scarches_analysis.py --query query.h5ad --reference-model-path ./ref_model --labels-key cell_type
"""

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import anndata
import matplotlib.pyplot as plt
import torch
from sklearn.metrics import classification_report

from scarches_core import (
    setup_reference_data,
    train_reference_model,
    train_classifier,
    prepare_query_data,
    map_query_to_reference,
    transfer_labels,
    verify_reference_preservation,
    get_latent_representation,
    save_model,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="scArches reference mapping pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Input
    parser.add_argument("--reference", default=None, help="Reference .h5ad (if training new model)")
    parser.add_argument("--reference-model-path", default=None, help="Path to existing reference model")
    parser.add_argument("--query", required=True, help="Query .h5ad to map")
    parser.add_argument("--labels-key", required=True, help="Cell type column in reference")
    parser.add_argument("--output-dir", default=None)

    # Reference model
    parser.add_argument("--n-latent", type=int, default=30)
    parser.add_argument("--n-layers", type=int, default=2)
    parser.add_argument("--ref-epochs", type=int, default=400)
    parser.add_argument("--batch-key", default=None)

    # Query mapping
    parser.add_argument("--query-epochs", type=int, default=200)

    # HVG
    parser.add_argument("--n-hvg", type=int, default=2000)

    parser.add_argument("--seed", type=int, default=0)

    return parser.parse_args()


def main():
    args = parse_args()

    scvi.settings.seed = args.seed
    torch.set_float32_matmul_precision("high")
    sc.set_figure_params(figsize=(6, 6), frameon=False)

    if args.output_dir is None:
        args.output_dir = f"{Path(args.query).stem}_scarches_results"
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Output: {args.output_dir}")

    # Reference model
    if args.reference_model_path:
        print(f"\n{'='*60}\nUsing existing reference model...\n{'='*60}")
        ref_model_path = args.reference_model_path

        # Load reference for classifier training
        if args.reference:
            adata_ref = sc.read_h5ad(args.reference)
            if "counts" not in adata_ref.layers:
                adata_ref.layers["counts"] = adata_ref.X.copy()
            model_ref = scvi.model.SCVI.load(ref_model_path, adata=adata_ref)
        else:
            model_ref = None
            adata_ref = None
    else:
        if not args.reference:
            raise ValueError("Must provide --reference or --reference-model-path")

        print(f"\n{'='*60}\nLoading reference...\n{'='*60}")
        adata_ref = sc.read_h5ad(args.reference)
        print(f"Reference: {adata_ref.n_obs} cells × {adata_ref.n_vars} genes")

        if "counts" not in adata_ref.layers:
            adata_ref.layers["counts"] = adata_ref.X.copy()

        # HVG selection
        sc.pp.highly_variable_genes(
            adata_ref, n_top_genes=args.n_hvg, batch_key=args.batch_key,
            flavor="seurat_v3", subset=True, layer="counts"
        )
        print(f"HVGs: {adata_ref.n_vars}")

        # Train reference
        print(f"\n{'='*60}\nTraining reference model...\n{'='*60}")
        setup_reference_data(adata_ref, layer="counts", batch_key=args.batch_key)
        model_ref = train_reference_model(
            adata_ref, n_latent=args.n_latent, n_layers=args.n_layers,
            max_epochs=args.ref_epochs
        )

        adata_ref.obsm["X_scVI"] = get_latent_representation(model_ref)

        ref_model_path = os.path.join(args.output_dir, "reference_model")
        save_model(model_ref, ref_model_path)

    # Train classifier
    if model_ref and adata_ref is not None:
        print(f"\n{'='*60}\nTraining classifier...\n{'='*60}")
        clf, _ = train_classifier(model_ref, adata_ref, labels_key=args.labels_key)
    else:
        clf = None

    # Load query
    print(f"\n{'='*60}\nLoading query...\n{'='*60}")
    adata_query = sc.read_h5ad(args.query)
    print(f"Query: {adata_query.n_obs} cells × {adata_query.n_vars} genes")

    if "counts" not in adata_query.layers:
        adata_query.layers["counts"] = adata_query.X.copy()

    # Prepare query
    print(f"\n{'='*60}\nPreparing query...\n{'='*60}")
    prepare_query_data(adata_query, ref_model_path)

    # Map query
    print(f"\n{'='*60}\nMapping query to reference...\n{'='*60}")
    model_query = map_query_to_reference(
        adata_query, ref_model_path, max_epochs=args.query_epochs
    )

    adata_query.obsm["X_scVI"] = get_latent_representation(model_query)

    # Transfer labels
    if clf:
        print(f"\n{'='*60}\nTransferring labels...\n{'='*60}")
        predictions, confidence = transfer_labels(clf, model_query, adata_query)
        adata_query.obs["predicted_cell_type"] = predictions
        adata_query.obs["prediction_confidence"] = confidence

        print(f"\nPredicted cell types:")
        print(adata_query.obs["predicted_cell_type"].value_counts())

    # Visualization
    print(f"\n{'='*60}\nVisualizing...\n{'='*60}")

    if adata_ref is not None and clf:
        # Combined visualization
        adata_combined = anndata.concat(
            [adata_ref, adata_query], label="dataset", keys=["reference", "query"]
        )
        adata_combined.obsm["X_scVI"] = model_query.get_latent_representation(adata_combined)

        sc.pp.neighbors(adata_combined, use_rep="X_scVI")
        sc.tl.umap(adata_combined)

        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        sc.pl.umap(adata_combined, color="dataset", ax=axes[0], show=False, title="Dataset")
        sc.pl.umap(adata_combined, color=args.labels_key, ax=axes[1], show=False,
                   title="Cell Types", na_color="lightgray")

        adata_combined.obs["display_type"] = adata_combined.obs[args.labels_key].copy()
        query_mask = adata_combined.obs["dataset"] == "query"
        adata_combined.obs.loc[query_mask, "display_type"] = \
            adata_combined.obs.loc[query_mask, "predicted_cell_type"]
        sc.pl.umap(adata_combined, color="display_type", ax=axes[2], show=False,
                   title="With Predictions")
        plt.tight_layout()
        plt.savefig(os.path.join(args.output_dir, "combined_umap.png"), dpi=150)
        plt.close()

        # Save combined
        adata_combined.write_h5ad(os.path.join(args.output_dir, "combined.h5ad"))

    # Confidence distribution
    if clf:
        plt.figure(figsize=(8, 4))
        plt.hist(adata_query.obs["prediction_confidence"], bins=50)
        plt.xlabel("Prediction Confidence")
        plt.ylabel("Cells")
        plt.axvline(0.7, color='r', linestyle='--', label='Threshold')
        plt.legend()
        plt.title("Label Transfer Confidence")
        plt.savefig(os.path.join(args.output_dir, "confidence_distribution.png"), dpi=150)
        plt.close()

    # Save
    print(f"\n{'='*60}\nSaving results...\n{'='*60}")
    save_model(model_query, os.path.join(args.output_dir, "query_model"))
    adata_query.write_h5ad(os.path.join(args.output_dir, f"{Path(args.query).stem}_mapped.h5ad"))

    if clf:
        pd.DataFrame({
            "cell_id": adata_query.obs_names,
            "predicted_type": predictions,
            "confidence": confidence,
        }).to_csv(os.path.join(args.output_dir, "predictions.csv"), index=False)

    print(f"\n{'='*60}\nscArches mapping complete!\n{'='*60}")
    print(f"Results: {args.output_dir}")


if __name__ == "__main__":
    main()
