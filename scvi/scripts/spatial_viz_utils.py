#!/usr/bin/env python3
"""
Advanced Visualization Utilities for Spatial Transcriptomics

Provides publication-quality visualization functions for:
- Cell type proportion maps
- Multi-panel comparisons
- Spatial statistics overlays
- Interactive exploration

Usage:
    from spatial_viz_utils import (
        plot_proportions_grid,
        plot_celltype_comparison,
        plot_dominant_celltype,
        plot_niche_composition,
        plot_spatial_correlation,
    )
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
import seaborn as sns
from typing import List, Optional, Union, Dict


def plot_proportions_grid(
    adata,
    proportions_key: str = "proportions",
    cell_types: Optional[List[str]] = None,
    ncols: int = 4,
    spot_size: float = 1.2,
    cmap: str = "viridis",
    vmax: Union[float, str] = 1.0,
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    dpi: int = 150
):
    """
    Plot cell type proportions in a grid layout.

    Parameters
    ----------
    adata : AnnData
        Spatial data with proportions in obsm
    proportions_key : str
        Key in obsm containing proportions DataFrame
    cell_types : list, optional
        Cell types to plot (default: all)
    ncols : int
        Number of columns in grid
    spot_size : float
        Size of spatial spots
    cmap : str
        Colormap name
    vmax : float or str
        Maximum color value (use "p95" for 95th percentile)
    title : str, optional
        Main title
    save_path : str, optional
        Path to save figure
    dpi : int
        DPI for saved figure
    """
    props = adata.obsm[proportions_key]
    if cell_types is None:
        cell_types = props.columns.tolist()

    nrows = int(np.ceil(len(cell_types) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
    axes = axes.flatten() if nrows > 1 or ncols > 1 else [axes]

    for idx, ct in enumerate(cell_types):
        adata.obs[f"_plot_{ct}"] = props[ct].values
        sc.pl.spatial(
            adata,
            color=f"_plot_{ct}",
            ax=axes[idx],
            show=False,
            title=ct,
            cmap=cmap,
            vmax=vmax if isinstance(vmax, float) else None,
            size=spot_size
        )

    # Hide unused axes
    for j in range(len(cell_types), len(axes)):
        axes[j].axis("off")

    if title:
        fig.suptitle(title, fontsize=14, y=1.02)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")
        plt.close()
    else:
        return fig


def plot_dominant_celltype(
    adata,
    proportions_key: str = "proportions",
    threshold: float = 0.1,
    spot_size: float = 1.2,
    palette: Optional[Dict[str, str]] = None,
    show_legend: bool = True,
    title: str = "Dominant Cell Type per Spot",
    save_path: Optional[str] = None,
    dpi: int = 150
):
    """
    Plot the dominant cell type in each spatial spot.

    Parameters
    ----------
    adata : AnnData
        Spatial data with proportions
    proportions_key : str
        Key in obsm for proportions
    threshold : float
        Minimum proportion to be considered dominant
    spot_size : float
        Size of spots
    palette : dict, optional
        Cell type to color mapping
    show_legend : bool
        Whether to show legend
    title : str
        Plot title
    save_path : str, optional
        Path to save figure
    dpi : int
        DPI for saved figure
    """
    props = adata.obsm[proportions_key]

    # Find dominant cell type
    dominant = props.idxmax(axis=1)
    max_prop = props.max(axis=1)

    # Mark low-confidence spots
    dominant[max_prop < threshold] = "Mixed/Uncertain"

    adata.obs["dominant_celltype"] = pd.Categorical(dominant)

    # Create palette if not provided
    if palette is None:
        cell_types = [ct for ct in props.columns if ct in dominant.unique()]
        colors = plt.cm.tab20(np.linspace(0, 1, len(cell_types)))
        palette = {ct: colors[i] for i, ct in enumerate(cell_types)}
        palette["Mixed/Uncertain"] = "lightgray"

    fig, ax = plt.subplots(figsize=(8, 8))
    sc.pl.spatial(
        adata,
        color="dominant_celltype",
        ax=ax,
        show=False,
        title=title,
        size=spot_size,
        palette=palette
    )

    if show_legend:
        handles = [mpatches.Patch(color=palette.get(ct, "gray"), label=ct)
                  for ct in adata.obs["dominant_celltype"].cat.categories]
        ax.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc="upper left")

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")
        plt.close()
    else:
        return fig


def plot_celltype_comparison(
    adata,
    cell_type: str,
    method_results: Dict[str, pd.DataFrame],
    spot_size: float = 1.2,
    cmap: str = "viridis",
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    dpi: int = 150
):
    """
    Compare a single cell type across multiple deconvolution methods.

    Parameters
    ----------
    adata : AnnData
        Spatial data
    cell_type : str
        Cell type to compare
    method_results : dict
        Dictionary of method_name -> proportions DataFrame
    spot_size : float
        Size of spots
    cmap : str
        Colormap
    title : str, optional
        Main title
    save_path : str, optional
        Path to save
    dpi : int
        Resolution
    """
    methods = list(method_results.keys())
    n_methods = len(methods)

    fig, axes = plt.subplots(1, n_methods + 1, figsize=(4 * (n_methods + 1), 4))

    # Global min/max for consistent coloring
    all_values = np.concatenate([method_results[m][cell_type].values for m in methods])
    vmax = np.percentile(all_values, 95)

    for idx, method in enumerate(methods):
        adata.obs[f"_cmp_{method}"] = method_results[method][cell_type].values
        sc.pl.spatial(
            adata,
            color=f"_cmp_{method}",
            ax=axes[idx],
            show=False,
            title=method,
            cmap=cmap,
            vmax=vmax,
            size=spot_size
        )

    # Add correlation scatter in last panel
    if n_methods >= 2:
        from scipy.stats import spearmanr

        ax = axes[-1]
        m1, m2 = methods[0], methods[1]
        x = method_results[m1][cell_type].values
        y = method_results[m2][cell_type].values

        ax.scatter(x, y, alpha=0.3, s=10)
        r, _ = spearmanr(x, y)
        ax.set_xlabel(m1)
        ax.set_ylabel(m2)
        ax.set_title(f"Correlation: r={r:.3f}")

        # Diagonal line
        lims = [0, max(x.max(), y.max())]
        ax.plot(lims, lims, "r--", alpha=0.5)
    else:
        axes[-1].axis("off")

    if title:
        fig.suptitle(title or f"{cell_type} Comparison", fontsize=14, y=1.02)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")
        plt.close()
    else:
        return fig


def plot_niche_composition(
    adata,
    cluster_key: str,
    proportions_key: str = "proportions",
    figsize: tuple = (12, 6),
    show_percentages: bool = True,
    title: str = "Niche Composition",
    save_path: Optional[str] = None,
    dpi: int = 150
):
    """
    Stacked bar chart showing cell type composition per spatial cluster/niche.

    Parameters
    ----------
    adata : AnnData
        Spatial data
    cluster_key : str
        Column in obs with cluster assignments
    proportions_key : str
        Key in obsm for proportions
    figsize : tuple
        Figure size
    show_percentages : bool
        Whether to show percentage labels
    title : str
        Plot title
    save_path : str, optional
        Path to save
    dpi : int
        Resolution
    """
    props = adata.obsm[proportions_key]
    clusters = adata.obs[cluster_key]

    # Aggregate by cluster
    cluster_props = props.groupby(clusters).mean()

    # Sort cell types by overall abundance
    cell_type_order = cluster_props.mean().sort_values(ascending=False).index

    # Create stacked bar chart
    fig, ax = plt.subplots(figsize=figsize)

    cluster_props[cell_type_order].plot(
        kind="bar",
        stacked=True,
        ax=ax,
        colormap="tab20"
    )

    ax.set_xlabel("Spatial Cluster")
    ax.set_ylabel("Mean Proportion")
    ax.set_title(title)
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", title="Cell Type")

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")
        plt.close()
    else:
        return fig


def plot_spatial_correlation(
    adata,
    proportions_key: str = "proportions",
    method: str = "spearman",
    figsize: tuple = (10, 8),
    cluster: bool = True,
    title: str = "Cell Type Co-localization",
    save_path: Optional[str] = None,
    dpi: int = 150
):
    """
    Plot correlation matrix showing cell type co-localization patterns.

    Parameters
    ----------
    adata : AnnData
        Spatial data
    proportions_key : str
        Key in obsm for proportions
    method : str
        Correlation method ("spearman" or "pearson")
    figsize : tuple
        Figure size
    cluster : bool
        Whether to cluster the heatmap
    title : str
        Plot title
    save_path : str, optional
        Path to save
    dpi : int
        Resolution
    """
    props = adata.obsm[proportions_key]
    corr_matrix = props.corr(method=method)

    fig, ax = plt.subplots(figsize=figsize)

    if cluster:
        g = sns.clustermap(
            corr_matrix,
            cmap="RdBu_r",
            center=0,
            vmin=-1,
            vmax=1,
            annot=True,
            fmt=".2f",
            figsize=figsize
        )
        g.fig.suptitle(title, y=1.02)
        fig = g.fig
    else:
        sns.heatmap(
            corr_matrix,
            cmap="RdBu_r",
            center=0,
            vmin=-1,
            vmax=1,
            annot=True,
            fmt=".2f",
            ax=ax
        )
        ax.set_title(title)
        plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")
        plt.close()
    else:
        return fig


def plot_proportion_distributions(
    adata,
    proportions_key: str = "proportions",
    cell_types: Optional[List[str]] = None,
    ncols: int = 4,
    figsize_per_panel: tuple = (3, 2.5),
    title: str = "Cell Type Proportion Distributions",
    save_path: Optional[str] = None,
    dpi: int = 150
):
    """
    Plot distribution of proportions for each cell type.

    Parameters
    ----------
    adata : AnnData
        Spatial data
    proportions_key : str
        Key in obsm for proportions
    cell_types : list, optional
        Cell types to plot
    ncols : int
        Number of columns
    figsize_per_panel : tuple
        Size of each subplot
    title : str
        Main title
    save_path : str, optional
        Path to save
    dpi : int
        Resolution
    """
    props = adata.obsm[proportions_key]
    if cell_types is None:
        cell_types = props.columns.tolist()

    nrows = int(np.ceil(len(cell_types) / ncols))
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows)
    )
    axes = axes.flatten()

    for idx, ct in enumerate(cell_types):
        values = props[ct].values
        axes[idx].hist(values, bins=50, edgecolor="black", alpha=0.7)
        axes[idx].axvline(values.mean(), color="red", linestyle="--", label=f"Mean: {values.mean():.3f}")
        axes[idx].set_title(ct, fontsize=10)
        axes[idx].set_xlabel("Proportion")
        axes[idx].legend(fontsize=8)

    # Hide unused axes
    for j in range(len(cell_types), len(axes)):
        axes[j].axis("off")

    fig.suptitle(title, fontsize=14)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")
        plt.close()
    else:
        return fig


def plot_spatial_with_tissue(
    adata,
    color: str,
    img_alpha: float = 0.5,
    spot_alpha: float = 0.8,
    spot_size: float = 1.5,
    cmap: str = "viridis",
    vmax: Union[float, str] = "p95",
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    dpi: int = 150
):
    """
    Plot spatial data overlaid on tissue image with adjustable transparency.

    Parameters
    ----------
    adata : AnnData
        Spatial data
    color : str
        Variable to color by
    img_alpha : float
        Transparency of tissue image
    spot_alpha : float
        Transparency of spots
    spot_size : float
        Size of spots
    cmap : str
        Colormap
    vmax : float or str
        Maximum color value
    title : str, optional
        Plot title
    save_path : str, optional
        Path to save
    dpi : int
        Resolution
    """
    fig, ax = plt.subplots(figsize=(8, 8))

    sc.pl.spatial(
        adata,
        color=color,
        ax=ax,
        show=False,
        title=title or color,
        cmap=cmap,
        vmax=vmax,
        size=spot_size,
        alpha_img=img_alpha,
        alpha=spot_alpha
    )

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")
        plt.close()
    else:
        return fig


def create_summary_figure(
    adata,
    proportions_key: str = "proportions",
    cluster_key: Optional[str] = None,
    top_n: int = 6,
    title: str = "Spatial Deconvolution Summary",
    save_path: Optional[str] = None,
    dpi: int = 150
):
    """
    Create a comprehensive summary figure combining multiple visualizations.

    Parameters
    ----------
    adata : AnnData
        Spatial data
    proportions_key : str
        Key in obsm for proportions
    cluster_key : str, optional
        Column with cluster assignments
    top_n : int
        Number of top cell types to show
    title : str
        Main title
    save_path : str, optional
        Path to save
    dpi : int
        Resolution
    """
    props = adata.obsm[proportions_key]

    # Get top cell types by mean proportion
    top_cts = props.mean().sort_values(ascending=False).head(top_n).index.tolist()

    fig = plt.figure(figsize=(20, 16))

    # Row 1: Spatial maps for top cell types
    for i, ct in enumerate(top_cts[:4]):
        ax = fig.add_subplot(4, 4, i + 1)
        adata.obs[f"_summary_{ct}"] = props[ct].values
        sc.pl.spatial(adata, color=f"_summary_{ct}", ax=ax, show=False, title=ct, cmap="viridis")

    # Row 2: Spatial maps continued + dominant cell type
    for i, ct in enumerate(top_cts[4:6]):
        ax = fig.add_subplot(4, 4, 5 + i)
        adata.obs[f"_summary_{ct}"] = props[ct].values
        sc.pl.spatial(adata, color=f"_summary_{ct}", ax=ax, show=False, title=ct, cmap="viridis")

    # Dominant cell type
    ax = fig.add_subplot(4, 4, 7)
    dominant = props.idxmax(axis=1)
    adata.obs["_dominant"] = pd.Categorical(dominant)
    sc.pl.spatial(adata, color="_dominant", ax=ax, show=False, title="Dominant Cell Type")

    # Distribution of dominant types
    ax = fig.add_subplot(4, 4, 8)
    adata.obs["_dominant"].value_counts().plot(kind="barh", ax=ax)
    ax.set_xlabel("Number of Spots")
    ax.set_title("Dominant Cell Type Counts")

    # Row 3: Co-localization heatmap
    ax = fig.add_subplot(4, 2, 5)
    corr = props[top_cts].corr(method="spearman")
    sns.heatmap(corr, cmap="RdBu_r", center=0, annot=True, fmt=".2f", ax=ax)
    ax.set_title("Cell Type Co-localization (Spearman r)")

    # Proportion distributions
    ax = fig.add_subplot(4, 2, 6)
    props[top_cts].boxplot(ax=ax)
    ax.set_ylabel("Proportion")
    ax.set_title("Proportion Distributions")
    ax.tick_params(axis="x", rotation=45)

    # Row 4: Niche composition (if clusters available)
    if cluster_key and cluster_key in adata.obs.columns:
        ax = fig.add_subplot(4, 1, 4)
        cluster_props = props.groupby(adata.obs[cluster_key]).mean()
        cluster_props[top_cts].plot(kind="bar", stacked=True, ax=ax, colormap="tab20")
        ax.set_xlabel("Cluster")
        ax.set_ylabel("Mean Proportion")
        ax.set_title("Cell Type Composition per Cluster")
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")

    fig.suptitle(title, fontsize=16, y=1.01)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")
        plt.close()
    else:
        return fig
