#!/usr/bin/env python3
"""
Cross-Method Comparison for Spatial Deconvolution

Compare results from Stereoscope, DestVI, and Cell2location on the same dataset.

Usage:
    python compare_deconvolution_methods.py \
        --spatial spatial.h5ad \
        --reference scrna.h5ad \
        --labels-key cell_type

    # Skip specific methods
    python compare_deconvolution_methods.py \
        --spatial spatial.h5ad \
        --reference scrna.h5ad \
        --labels-key cell_type \
        --skip-cell2location
"""

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import matplotlib.pyplot as plt
import seaborn as sns
import torch
from scipy.stats import spearmanr, pearsonr

# Import method-specific models
from scvi.external import RNAStereoscope, SpatialStereoscope
from scvi.model import CondSCVI

try:
    from cell2location.models import Cell2location, RegressionModel
    from cell2location.utils.filtering import filter_genes
    CELL2LOC_AVAILABLE = True
except ImportError:
    CELL2LOC_AVAILABLE = False
    print("Cell2location not available. Install with: pip install cell2location")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare spatial deconvolution methods",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--spatial", required=True, help="Spatial .h5ad file")
    parser.add_argument("--reference", required=True, help="scRNA-seq reference .h5ad")
    parser.add_argument("--labels-key", required=True, help="Cell type column")
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--n-hvg", type=int, default=3000, help="Number of HVGs")
    parser.add_argument("--skip-stereoscope", action="store_true")
    parser.add_argument("--skip-destvi", action="store_true")
    parser.add_argument("--skip-cell2location", action="store_true")
    parser.add_argument("--seed", type=int, default=0)

    return parser.parse_args()


def prepare_data(st_adata, sc_adata, labels_key, n_hvg):
    """Prepare data for deconvolution methods."""
    # Ensure counts layer
    if "counts" not in st_adata.layers:
        st_adata.layers["counts"] = st_adata.X.copy()
    if "counts" not in sc_adata.layers:
        sc_adata.layers["counts"] = sc_adata.X.copy()

    # Remove MT genes
    st_adata.var["mt"] = st_adata.var_names.str.startswith("MT-")
    sc_adata.var["mt"] = sc_adata.var_names.str.startswith("MT-")
    st_adata = st_adata[:, ~st_adata.var["mt"]].copy()
    sc_adata = sc_adata[:, ~sc_adata.var["mt"]].copy()

    # Find shared genes
    shared_genes = np.intersect1d(st_adata.var_names, sc_adata.var_names)
    st_adata = st_adata[:, shared_genes].copy()
    sc_adata = sc_adata[:, shared_genes].copy()

    # HVG selection
    sc.pp.highly_variable_genes(
        sc_adata,
        n_top_genes=min(n_hvg, len(shared_genes)),
        layer="counts",
        flavor="seurat_v3",
        subset=True
    )

    # Align to HVGs
    hvg = sc_adata.var_names
    st_adata = st_adata[:, hvg].copy()

    print(f"Final genes: {len(hvg)}")
    return st_adata, sc_adata


def run_stereoscope(st_adata, sc_adata, labels_key, output_dir):
    """Run Stereoscope deconvolution."""
    print("\n" + "="*60)
    print("Running Stereoscope...")
    print("="*60)

    # Setup and train reference
    RNAStereoscope.setup_anndata(sc_adata, layer="counts", labels_key=labels_key)
    sc_model = RNAStereoscope(sc_adata)
    sc_model.train(max_epochs=100)

    # Plot training
    plt.figure(figsize=(8, 4))
    plt.plot(sc_model.history["elbo_train"].values)
    plt.xlabel("Epoch")
    plt.ylabel("ELBO")
    plt.title("Stereoscope Reference Training")
    plt.savefig(os.path.join(output_dir, "stereoscope_ref_training.png"), dpi=150)
    plt.close()

    # Setup and train spatial
    SpatialStereoscope.setup_anndata(st_adata, layer="counts")
    st_model = SpatialStereoscope.from_rna_model(st_adata, sc_model)
    st_model.train(max_epochs=2000)

    # Plot training
    plt.figure(figsize=(8, 4))
    plt.plot(st_model.history["elbo_train"].values)
    plt.xlabel("Epoch")
    plt.ylabel("ELBO")
    plt.title("Stereoscope Spatial Training")
    plt.savefig(os.path.join(output_dir, "stereoscope_spatial_training.png"), dpi=150)
    plt.close()

    # Get proportions
    proportions = st_model.get_proportions()
    sc_model.save(os.path.join(output_dir, "stereoscope_ref_model"), overwrite=True)
    st_model.save(os.path.join(output_dir, "stereoscope_spatial_model"), overwrite=True)

    return proportions


def run_destvi(st_adata, sc_adata, labels_key, output_dir):
    """Run DestVI deconvolution."""
    print("\n" + "="*60)
    print("Running DestVI...")
    print("="*60)

    from scvi.model import DestVI

    # Setup and train reference (CondSCVI)
    CondSCVI.setup_anndata(sc_adata, layer="counts", labels_key=labels_key)
    sc_model = CondSCVI(sc_adata)
    sc_model.train(max_epochs=300)

    # Plot training
    plt.figure(figsize=(8, 4))
    plt.plot(sc_model.history["elbo_train"].values)
    plt.xlabel("Epoch")
    plt.ylabel("ELBO")
    plt.title("DestVI Reference Training")
    plt.savefig(os.path.join(output_dir, "destvi_ref_training.png"), dpi=150)
    plt.close()

    # Setup and train spatial
    DestVI.setup_anndata(st_adata, layer="counts")
    st_model = DestVI.from_rna_model(st_adata, sc_model)
    st_model.train(max_epochs=2500)

    # Plot training
    plt.figure(figsize=(8, 4))
    plt.plot(st_model.history["elbo_train"].values)
    plt.xlabel("Epoch")
    plt.ylabel("ELBO")
    plt.title("DestVI Spatial Training")
    plt.savefig(os.path.join(output_dir, "destvi_spatial_training.png"), dpi=150)
    plt.close()

    # Get proportions
    proportions = st_model.get_proportions()
    sc_model.save(os.path.join(output_dir, "destvi_ref_model"), overwrite=True)
    st_model.save(os.path.join(output_dir, "destvi_spatial_model"), overwrite=True)

    return proportions


def run_cell2location(st_adata, sc_adata, labels_key, output_dir):
    """Run Cell2location deconvolution."""
    if not CELL2LOC_AVAILABLE:
        print("Cell2location not available, skipping...")
        return None

    print("\n" + "="*60)
    print("Running Cell2location...")
    print("="*60)

    # Gene filtering
    selected_genes = filter_genes(
        sc_adata,
        cell_count_cutoff=5,
        cell_percentage_cutoff2=0.03,
        nonz_mean_cutoff=1.12
    )
    sc_adata_c2l = sc_adata[:, selected_genes].copy()
    st_adata_c2l = st_adata[:, np.intersect1d(st_adata.var_names, selected_genes)].copy()

    # Setup and train regression model
    RegressionModel.setup_anndata(sc_adata_c2l, labels_key=labels_key)
    ref_model = RegressionModel(sc_adata_c2l)
    ref_model.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002)

    # Export signatures
    sc_adata_c2l = ref_model.export_posterior(
        sc_adata_c2l,
        sample_kwargs={"num_samples": 1000, "batch_size": 2500}
    )

    # Extract signature matrix
    if "means_per_cluster_mu_fg" in sc_adata_c2l.varm.keys():
        inf_aver = sc_adata_c2l.varm["means_per_cluster_mu_fg"][
            [f"means_per_cluster_mu_fg_{i}" for i in sc_adata_c2l.uns["mod"]["factor_names"]]
        ].copy()
    else:
        inf_aver = sc_adata_c2l.var[
            [f"means_per_cluster_mu_fg_{i}" for i in sc_adata_c2l.uns["mod"]["factor_names"]]
        ].copy()
    inf_aver.columns = sc_adata_c2l.uns["mod"]["factor_names"]

    # Align genes
    shared_genes = np.intersect1d(st_adata_c2l.var_names, inf_aver.index)
    st_adata_c2l = st_adata_c2l[:, shared_genes].copy()
    inf_aver = inf_aver.loc[shared_genes, :].copy()

    # Add sample key if missing
    if "sample" not in st_adata_c2l.obs.columns:
        st_adata_c2l.obs["sample"] = "sample_1"

    # Setup and train Cell2location
    Cell2location.setup_anndata(st_adata_c2l, batch_key="sample")
    spatial_model = Cell2location(
        st_adata_c2l,
        cell_state_df=inf_aver,
        N_cells_per_location=30,
        detection_alpha=200
    )
    spatial_model.train(max_epochs=30000, batch_size=None, train_size=1)

    # Export results
    st_adata_c2l = spatial_model.export_posterior(
        st_adata_c2l,
        sample_kwargs={"num_samples": 1000, "batch_size": st_adata_c2l.n_obs}
    )

    # Get proportions (normalize abundances to proportions)
    abundances = st_adata_c2l.obsm["q05_cell_abundance_w_sf"]
    proportions = pd.DataFrame(
        abundances / abundances.sum(axis=1, keepdims=True),
        index=st_adata_c2l.obs_names,
        columns=inf_aver.columns
    )

    ref_model.save(os.path.join(output_dir, "cell2location_ref_model"), overwrite=True)
    spatial_model.save(os.path.join(output_dir, "cell2location_spatial_model"), overwrite=True)

    return proportions


def compare_methods(results, output_dir, st_adata):
    """Generate comparison visualizations."""
    print("\n" + "="*60)
    print("Generating comparisons...")
    print("="*60)

    # Get common cell types across all methods
    all_cell_types = set.intersection(*[set(df.columns) for df in results.values()])
    cell_types = sorted(list(all_cell_types))

    methods = list(results.keys())
    n_methods = len(methods)

    # 1. Correlation heatmap for each cell type
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for i, ct in enumerate(cell_types[:6]):
        if i >= len(axes):
            break

        # Build correlation matrix
        ct_data = pd.DataFrame({
            m: results[m][ct] for m in methods
        })
        corr = ct_data.corr(method="spearman")

        sns.heatmap(corr, annot=True, cmap="RdYlGn", vmin=0, vmax=1,
                   ax=axes[i], fmt=".2f")
        axes[i].set_title(f"{ct}")

    # Hide unused axes
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    plt.suptitle("Method Agreement per Cell Type (Spearman r)", fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "method_correlation_heatmaps.png"), dpi=150)
    plt.close()

    # 2. Pairwise scatter plots for top cell types
    for ct in cell_types[:3]:
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))

        pairs = [(0, 1), (0, 2), (1, 2)] if n_methods >= 3 else [(0, 1)]
        for idx, (i, j) in enumerate(pairs):
            if i >= n_methods or j >= n_methods:
                continue

            m1, m2 = methods[i], methods[j]
            x = results[m1][ct].values
            y = results[m2][ct].values

            axes[idx].scatter(x, y, alpha=0.3, s=10)
            r, _ = spearmanr(x, y)
            axes[idx].set_xlabel(m1)
            axes[idx].set_ylabel(m2)
            axes[idx].set_title(f"r = {r:.3f}")

            # Add diagonal
            lims = [min(x.min(), y.min()), max(x.max(), y.max())]
            axes[idx].plot(lims, lims, "r--", alpha=0.5)

        plt.suptitle(f"{ct} Proportions: Method Comparison", fontsize=14)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"scatter_{ct}.png"), dpi=150)
        plt.close()

    # 3. Spatial visualization comparison
    if "spatial" in st_adata.obsm:
        for ct in cell_types[:3]:
            fig, axes = plt.subplots(1, n_methods, figsize=(5 * n_methods, 4))
            if n_methods == 1:
                axes = [axes]

            for idx, method in enumerate(methods):
                st_adata.obs[f"_{method}_{ct}"] = results[method][ct].values
                sc.pl.spatial(
                    st_adata,
                    color=f"_{method}_{ct}",
                    ax=axes[idx],
                    show=False,
                    title=f"{method}",
                    cmap="viridis",
                    vmax="p95"
                )

            plt.suptitle(f"{ct} - Spatial Comparison", fontsize=14)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"spatial_{ct}.png"), dpi=150)
            plt.close()

    # 4. Summary statistics table
    summary_data = []
    for ct in cell_types:
        row = {"cell_type": ct}
        for method in methods:
            row[f"{method}_mean"] = results[method][ct].mean()
            row[f"{method}_std"] = results[method][ct].std()

        # Pairwise correlations
        if n_methods >= 2:
            for i in range(n_methods):
                for j in range(i + 1, n_methods):
                    m1, m2 = methods[i], methods[j]
                    r, _ = spearmanr(results[m1][ct], results[m2][ct])
                    row[f"corr_{m1}_{m2}"] = r

        summary_data.append(row)

    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(os.path.join(output_dir, "method_comparison_summary.csv"), index=False)

    # 5. Overall agreement bar chart
    avg_correlations = []
    for method in methods:
        other_methods = [m for m in methods if m != method]
        corrs = []
        for ct in cell_types:
            for other in other_methods:
                r, _ = spearmanr(results[method][ct], results[other][ct])
                corrs.append(r)
        avg_correlations.append(np.mean(corrs))

    plt.figure(figsize=(8, 5))
    plt.bar(methods, avg_correlations, color=["#1f77b4", "#ff7f0e", "#2ca02c"][:n_methods])
    plt.ylabel("Average Correlation with Other Methods")
    plt.title("Method Agreement (Higher = More Consistent)")
    plt.ylim(0, 1)
    for i, v in enumerate(avg_correlations):
        plt.text(i, v + 0.02, f"{v:.3f}", ha="center")
    plt.savefig(os.path.join(output_dir, "method_agreement_overall.png"), dpi=150)
    plt.close()

    print(f"\nSummary saved to: {output_dir}/method_comparison_summary.csv")
    return summary_df


def main():
    args = parse_args()

    scvi.settings.seed = args.seed
    torch.set_float32_matmul_precision("high")
    sc.set_figure_params(figsize=(6, 6), frameon=False)

    if args.output_dir is None:
        args.output_dir = f"{Path(args.spatial).stem}_comparison"
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Output directory: {args.output_dir}")

    # Load data
    print("\nLoading data...")
    st_adata = sc.read_h5ad(args.spatial)
    sc_adata = sc.read_h5ad(args.reference)

    print(f"Spatial: {st_adata.n_obs} spots × {st_adata.n_vars} genes")
    print(f"Reference: {sc_adata.n_obs} cells × {sc_adata.n_vars} genes")

    # Prepare data
    st_adata, sc_adata = prepare_data(st_adata, sc_adata, args.labels_key, args.n_hvg)

    # Run methods
    results = {}

    if not args.skip_stereoscope:
        results["Stereoscope"] = run_stereoscope(
            st_adata.copy(), sc_adata.copy(), args.labels_key, args.output_dir
        )

    if not args.skip_destvi:
        results["DestVI"] = run_destvi(
            st_adata.copy(), sc_adata.copy(), args.labels_key, args.output_dir
        )

    if not args.skip_cell2location and CELL2LOC_AVAILABLE:
        c2l_result = run_cell2location(
            st_adata.copy(), sc_adata.copy(), args.labels_key, args.output_dir
        )
        if c2l_result is not None:
            results["Cell2location"] = c2l_result

    if len(results) < 2:
        print("\nNeed at least 2 methods for comparison. Exiting.")
        return

    # Save individual results
    for method, props in results.items():
        props.to_csv(os.path.join(args.output_dir, f"{method.lower()}_proportions.csv"))

    # Compare methods
    summary = compare_methods(results, args.output_dir, st_adata)

    print("\n" + "="*60)
    print("Comparison complete!")
    print("="*60)
    print(f"\nResults saved to: {args.output_dir}")
    print("\nGenerated files:")
    print("  - method_comparison_summary.csv")
    print("  - method_correlation_heatmaps.png")
    print("  - method_agreement_overall.png")
    print("  - scatter_<celltype>.png")
    print("  - spatial_<celltype>.png")


if __name__ == "__main__":
    main()
