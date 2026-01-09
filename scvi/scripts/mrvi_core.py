"""
Core functions for MrVI multi-sample analysis.

This module provides modular utility functions for single-cell multi-sample
analysis using MrVI from scvi-tools. Functions can be used independently
or as part of the complete mrvi_analysis.py pipeline.
"""

import scvi
from scvi.external import MRVI
import scanpy as sc
import numpy as np
import pandas as pd
from typing import Optional, List, Union
import xarray as xr


def setup_anndata_mrvi(
    adata,
    sample_key: str,
    batch_key: Optional[str] = None,
    layer: Optional[str] = None,
    categorical_covariate_keys: Optional[List[str]] = None,
    continuous_covariate_keys: Optional[List[str]] = None
) -> None:
    """
    Register AnnData object with MrVI.

    This function prepares the AnnData object for use with MrVI by registering
    the count data, sample information, and optional covariates.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with raw counts.
    sample_key : str
        Key in adata.obs for sample identifiers.
    batch_key : str, optional
        Key in adata.obs for nuisance covariates (e.g., batch).
    layer : str, optional
        Key in adata.layers containing counts. If None, uses adata.X.
    categorical_covariate_keys : list, optional
        Additional categorical covariates.
    continuous_covariate_keys : list, optional
        Continuous covariates to include.

    Returns
    -------
    None
        Modifies adata in place with MrVI registration.
    """
    print(f"  Registering AnnData with MrVI...")
    print(f"    Sample key: {sample_key}")
    if batch_key:
        print(f"    Batch key: {batch_key}")

    MRVI.setup_anndata(
        adata,
        sample_key=sample_key,
        batch_key=batch_key,
        layer=layer,
        categorical_covariate_keys=categorical_covariate_keys,
        continuous_covariate_keys=continuous_covariate_keys
    )
    print("  AnnData registered successfully")


def train_mrvi_model(
    adata,
    n_latent: int = 30,
    max_epochs: int = 400,
    early_stopping: bool = True,
    batch_size: int = 256,
    **model_kwargs
):
    """
    Train an MrVI model for multi-sample analysis.

    MrVI learns sample-aware representations that capture both cell-intrinsic
    states and sample-specific effects.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix (must be registered with setup_anndata_mrvi first).
    n_latent : int, default=30
        Dimensionality of the latent space.
    max_epochs : int, default=400
        Maximum number of training epochs.
    early_stopping : bool, default=True
        Whether to use early stopping.
    batch_size : int, default=256
        Minibatch size for training.
    **model_kwargs
        Additional kwargs passed to MRVI model constructor.

    Returns
    -------
    MRVI
        Trained MrVI model.
    """
    print(f"  Initializing MrVI model...")
    print(f"    Latent dimensions: {n_latent}")

    model = MRVI(adata, n_latent=n_latent, **model_kwargs)

    print(f"  Training MrVI model (max {max_epochs} epochs)...")
    model.train(
        max_epochs=max_epochs,
        early_stopping=early_stopping,
        batch_size=batch_size
    )

    print(f"  Training complete!")

    return model


def get_latent_representation(model, adata=None) -> np.ndarray:
    """
    Extract the sample-independent latent representation (u) from MrVI.

    Parameters
    ----------
    model : MRVI
        Trained MrVI model.
    adata : AnnData, optional
        Data to get representation for. If None, uses training data.

    Returns
    -------
    np.ndarray
        Latent representation matrix (n_cells × n_latent).
    """
    return model.get_latent_representation(adata=adata)


def compute_sample_distances(
    model,
    groupby: Optional[str] = None,
    keep_cell: bool = False,
    batch_size: int = 32
) -> xr.DataArray:
    """
    Compute local sample distances.

    MrVI computes distances between samples in the learned latent space,
    optionally grouped by cell type or other categorical variables.

    Parameters
    ----------
    model : MRVI
        Trained MrVI model.
    groupby : str, optional
        Column in adata.obs to group cells by (e.g., 'cell_type').
        If provided, returns distances per group.
    keep_cell : bool, default=False
        If True, return cell-level distances (memory intensive).
        If False, return aggregated distances per group.
    batch_size : int, default=32
        Batch size for computing distances.

    Returns
    -------
    xr.DataArray
        Sample distance matrix/matrices.
    """
    print(f"  Computing sample distances...")
    if groupby:
        print(f"    Grouped by: {groupby}")

    distances = model.get_local_sample_distances(
        keep_cell=keep_cell,
        groupby=groupby,
        batch_size=batch_size
    )

    return distances


def run_differential_expression(
    model,
    sample_cov_keys: List[str],
    store_lfc: bool = True
):
    """
    Run differential expression analysis linked to sample covariates.

    MrVI's DE analysis provides cell-level effect sizes that capture how
    sample-level covariates affect gene expression.

    Parameters
    ----------
    model : MRVI
        Trained MrVI model.
    sample_cov_keys : list of str
        Sample-level covariates to test (e.g., ['condition', 'treatment']).
    store_lfc : bool, default=True
        Whether to store log fold changes.

    Returns
    -------
    object
        DE results with effect_size and optionally lfc attributes.
    """
    print(f"  Running differential expression...")
    print(f"    Covariates: {sample_cov_keys}")

    de_results = model.differential_expression(
        sample_cov_keys=sample_cov_keys,
        store_lfc=store_lfc
    )

    return de_results


def run_differential_abundance(
    model,
    sample_cov_keys: List[str]
):
    """
    Run differential abundance analysis.

    Detects compositional changes in cell states between conditions.

    Parameters
    ----------
    model : MRVI
        Trained MrVI model.
    sample_cov_keys : list of str
        Sample-level covariates to test.

    Returns
    -------
    object
        DA results with log probability ratios.
    """
    print(f"  Running differential abundance...")
    print(f"    Covariates: {sample_cov_keys}")

    da_results = model.differential_abundance(sample_cov_keys=sample_cov_keys)

    return da_results


def compute_da_log_ratio(
    da_results,
    covariate: str,
    numerator: str,
    denominator: str
) -> np.ndarray:
    """
    Compute log ratio of abundances between two conditions.

    Parameters
    ----------
    da_results : object
        Results from run_differential_abundance.
    covariate : str
        Name of the covariate.
    numerator : str
        Category to use as numerator.
    denominator : str
        Category to use as denominator.

    Returns
    -------
    np.ndarray
        Log ratio values per cell.
    """
    log_probs_key = f"{covariate}_log_probs"
    log_probs = getattr(da_results, log_probs_key)

    num_probs = log_probs.loc[{covariate: numerator}]
    denom_probs = log_probs.loc[{covariate: denominator}]

    return (num_probs - denom_probs).values


def save_model(model, path: str) -> None:
    """
    Save a trained MrVI model.

    Parameters
    ----------
    model : MRVI
        Trained model to save.
    path : str
        Directory path to save model.
    """
    model.save(path, overwrite=True)


def load_model(path: str, adata):
    """
    Load a saved MrVI model.

    Parameters
    ----------
    path : str
        Directory path containing saved model.
    adata : AnnData
        AnnData object (must match the data used for training).

    Returns
    -------
    MRVI
        Loaded model.
    """
    return MRVI.load(path, adata=adata)


def get_sample_info(model) -> pd.DataFrame:
    """
    Get sample-level information from the model.

    Parameters
    ----------
    model : MRVI
        Trained MrVI model.

    Returns
    -------
    pd.DataFrame
        DataFrame with sample metadata.
    """
    return model.sample_info


def reorder_covariate_categories(
    model,
    covariate: str,
    order: List[str]
) -> None:
    """
    Reorder categories of a covariate for DE/DA analysis.

    This affects which category is used as the reference in comparisons.

    Parameters
    ----------
    model : MRVI
        Trained MrVI model.
    covariate : str
        Name of the covariate to reorder.
    order : list of str
        Desired order of categories (first will be reference).

    Returns
    -------
    None
        Modifies model.sample_info in place.
    """
    model.sample_info[covariate] = model.sample_info[covariate].cat.reorder_categories(order)
    print(f"  Reordered '{covariate}' categories: {order}")


def extract_de_genes(
    de_results,
    covariate: str,
    effect_threshold: float = 0.5,
    top_n: Optional[int] = None
) -> pd.DataFrame:
    """
    Extract significant DE genes for a covariate.

    Parameters
    ----------
    de_results : object
        Results from run_differential_expression.
    covariate : str
        Covariate name to extract genes for.
    effect_threshold : float, default=0.5
        Minimum absolute effect size threshold.
    top_n : int, optional
        Return only top N genes by effect size.

    Returns
    -------
    pd.DataFrame
        DataFrame with gene names and effect sizes.
    """
    if not hasattr(de_results, 'lfc'):
        raise ValueError("DE results don't contain LFC. Run with store_lfc=True.")

    # Get log fold changes
    lfc = de_results.lfc.sel(covariate=covariate)

    # Compute mean absolute LFC across cells
    mean_lfc = lfc.mean(dim='cell').values
    genes = lfc.coords['gene'].values

    df = pd.DataFrame({
        'gene': genes,
        'mean_lfc': mean_lfc,
        'abs_lfc': np.abs(mean_lfc)
    })

    # Filter by threshold
    df = df[df['abs_lfc'] >= effect_threshold]

    # Sort by absolute effect
    df = df.sort_values('abs_lfc', ascending=False)

    if top_n:
        df = df.head(top_n)

    return df.reset_index(drop=True)
