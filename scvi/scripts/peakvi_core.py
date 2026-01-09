#!/usr/bin/env python3
"""
PeakVI Core Functions

Modular functions for scATAC-seq analysis with PeakVI.
"""

import numpy as np
import scvi
from typing import Optional, Dict, List


def setup_anndata_peakvi(adata, batch_key: Optional[str] = None):
    """
    Register AnnData for PeakVI.

    Parameters
    ----------
    adata : AnnData
        Peak-by-cell matrix
    batch_key : str, optional
        Column in obs for batch correction
    """
    scvi.model.PEAKVI.setup_anndata(adata, batch_key=batch_key)
    print(f"AnnData registered for PeakVI")
    if batch_key:
        print(f"  Batch key: {batch_key}")


def train_peakvi_model(
    adata,
    n_latent: int = 13,
    n_hidden: int = 128,
    n_layers_encoder: int = 2,
    n_layers_decoder: int = 2,
    dropout_rate: float = 0.1,
    max_epochs: int = 500,
    early_stopping: bool = True,
    batch_size: int = 128,
) -> scvi.model.PEAKVI:
    """
    Train PeakVI model.

    Parameters
    ----------
    adata : AnnData
        Registered AnnData object
    n_latent : int
        Latent space dimensions
    n_hidden : int
        Hidden layer size
    n_layers_encoder : int
        Number of encoder layers
    n_layers_decoder : int
        Number of decoder layers
    dropout_rate : float
        Dropout for regularization
    max_epochs : int
        Maximum training epochs
    early_stopping : bool
        Enable early stopping
    batch_size : int
        Training batch size

    Returns
    -------
    scvi.model.PEAKVI
        Trained model
    """
    model = scvi.model.PEAKVI(
        adata,
        n_latent=n_latent,
        n_hidden=n_hidden,
        n_layers_encoder=n_layers_encoder,
        n_layers_decoder=n_layers_decoder,
        dropout_rate=dropout_rate,
    )

    print(f"\nModel architecture:")
    print(f"  Latent dimensions: {n_latent}")
    print(f"  Hidden size: {n_hidden}")
    print(f"  Encoder layers: {n_layers_encoder}")
    print(f"  Decoder layers: {n_layers_decoder}")

    model.train(
        max_epochs=max_epochs,
        early_stopping=early_stopping,
        batch_size=batch_size,
    )

    final_epoch = len(model.history["elbo_train"])
    print(f"\nTraining complete after {final_epoch} epochs")

    return model


def get_latent_representation(model: scvi.model.PEAKVI, adata=None) -> np.ndarray:
    """Extract latent representation."""
    return model.get_latent_representation(adata)


def run_differential_accessibility(
    model: scvi.model.PEAKVI,
    groupby: str,
    group1: str,
    group2: Optional[str] = None,
    batch_correction: bool = False,
    fdr_target: float = 0.05,
) -> "pd.DataFrame":
    """
    Run differential accessibility analysis.

    Parameters
    ----------
    model : scvi.model.PEAKVI
        Trained model
    groupby : str
        Column for comparison
    group1 : str
        Target group
    group2 : str, optional
        Reference group (None = vs rest)
    batch_correction : bool
        Account for batch effects
    fdr_target : float
        FDR threshold

    Returns
    -------
    pd.DataFrame
        DA results with prob_da, is_da_fdr, bayes_factor, etc.
    """
    return model.differential_accessibility(
        groupby=groupby,
        group1=group1,
        group2=group2,
        batch_correction=batch_correction,
        fdr_target=fdr_target,
    )


def get_cluster_markers(
    model: scvi.model.PEAKVI,
    adata,
    cluster_key: str,
    n_top: int = 50,
    prob_da_threshold: float = 0.8,
) -> Dict[str, "pd.DataFrame"]:
    """
    Find top DA peaks for each cluster vs rest.

    Parameters
    ----------
    model : scvi.model.PEAKVI
        Trained model
    adata : AnnData
        Data with cluster assignments
    cluster_key : str
        Column containing cluster labels
    n_top : int
        Number of top markers per cluster
    prob_da_threshold : float
        Minimum prob_da for marker consideration

    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary mapping cluster -> marker peaks DataFrame
    """
    clusters = adata.obs[cluster_key].unique()
    all_markers = {}

    for cluster in clusters:
        da = model.differential_accessibility(
            groupby=cluster_key,
            group1=str(cluster),
            group2=None,
        )
        markers = da[da["prob_da"] > prob_da_threshold].nlargest(n_top, "bayes_factor")
        all_markers[str(cluster)] = markers
        print(f"Cluster {cluster}: {len(markers)} marker peaks")

    return all_markers


def filter_peaks_by_detection(
    adata, min_detection: float = 0.05
) -> "anndata.AnnData":
    """
    Filter peaks based on detection rate.

    Parameters
    ----------
    adata : AnnData
        Peak matrix
    min_detection : float
        Minimum fraction of cells with peak detected

    Returns
    -------
    AnnData
        Filtered data
    """
    import scanpy as sc

    adata.var["n_cells"] = np.array((adata.X > 0).sum(axis=0)).flatten()
    adata.var["detection_rate"] = adata.var["n_cells"] / adata.n_obs

    min_cells = int(adata.n_obs * min_detection)
    n_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=min_cells)
    n_after = adata.n_vars

    print(f"Filtered peaks: {n_before} -> {n_after} (>{min_detection*100:.1f}% detection)")
    return adata


def calculate_qc_metrics(adata):
    """Calculate cell-level QC metrics for scATAC data."""
    adata.obs["n_peaks"] = np.array((adata.X > 0).sum(axis=1)).flatten()
    adata.obs["total_counts"] = np.array(adata.X.sum(axis=1)).flatten()
    print(f"QC metrics added: n_peaks, total_counts")


def save_model(model: scvi.model.PEAKVI, path: str, overwrite: bool = True):
    """Save PeakVI model."""
    model.save(path, overwrite=overwrite)
    print(f"Model saved to: {path}")


def load_model(path: str, adata) -> scvi.model.PEAKVI:
    """Load PeakVI model."""
    return scvi.model.PEAKVI.load(path, adata=adata)
