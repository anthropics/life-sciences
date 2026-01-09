#!/usr/bin/env python3
"""
DestVI Core Functions

Modular functions for spatial transcriptomics deconvolution with DestVI.
"""

import numpy as np
import scvi
from scvi.model import CondSCVI, DestVI
from typing import Optional, Dict


def setup_reference_data(
    sc_adata,
    layer: str = "counts",
    labels_key: str = "cell_type",
):
    """
    Register scRNA-seq reference data for CondSCVI.

    Parameters
    ----------
    sc_adata : AnnData
        scRNA-seq reference with cell type annotations
    layer : str
        Layer containing raw counts
    labels_key : str
        Column with cell type labels
    """
    CondSCVI.setup_anndata(sc_adata, layer=layer, labels_key=labels_key)
    print(f"Reference data registered")
    print(f"  Layer: {layer}")
    print(f"  Labels key: {labels_key}")
    print(f"  Cell types: {sc_adata.obs[labels_key].nunique()}")


def train_reference_model(
    sc_adata,
    weight_obs: bool = False,
    max_epochs: int = 300,
    early_stopping: bool = True,
) -> CondSCVI:
    """
    Train reference model (CondSCVI).

    Parameters
    ----------
    sc_adata : AnnData
        Registered scRNA-seq reference
    weight_obs : bool
        Weight cells by inverse frequency (for imbalanced data)
    max_epochs : int
        Maximum training epochs
    early_stopping : bool
        Enable early stopping

    Returns
    -------
    CondSCVI
        Trained reference model
    """
    model = CondSCVI(sc_adata, weight_obs=weight_obs)

    print(f"\nTraining reference model (CondSCVI):")
    print(f"  Weight observations: {weight_obs}")
    print(f"  Max epochs: {max_epochs}")

    model.train(max_epochs=max_epochs, early_stopping=early_stopping)

    final_epoch = len(model.history["elbo_train"])
    print(f"Training complete after {final_epoch} epochs")

    return model


def setup_spatial_data(st_adata, layer: str = "counts"):
    """
    Register spatial data for DestVI.

    Parameters
    ----------
    st_adata : AnnData
        Spatial transcriptomics data
    layer : str
        Layer containing raw counts
    """
    DestVI.setup_anndata(st_adata, layer=layer)
    print(f"Spatial data registered")
    print(f"  Layer: {layer}")
    print(f"  Spots: {st_adata.n_obs}")


def train_spatial_model(
    st_adata,
    reference_model: CondSCVI,
    max_epochs: int = 2500,
    early_stopping: bool = True,
) -> DestVI:
    """
    Train spatial deconvolution model (DestVI).

    Parameters
    ----------
    st_adata : AnnData
        Registered spatial data
    reference_model : CondSCVI
        Trained reference model
    max_epochs : int
        Maximum training epochs (min 1000 recommended)
    early_stopping : bool
        Enable early stopping

    Returns
    -------
    DestVI
        Trained spatial model
    """
    model = DestVI.from_rna_model(st_adata, reference_model)

    print(f"\nTraining spatial model (DestVI):")
    print(f"  Max epochs: {max_epochs}")

    model.train(max_epochs=max_epochs, early_stopping=early_stopping)

    final_epoch = len(model.history["elbo_train"])
    print(f"Training complete after {final_epoch} epochs")

    return model


def get_proportions(model: DestVI) -> "pd.DataFrame":
    """
    Get cell type proportions for each spot.

    Returns
    -------
    pd.DataFrame
        Proportions matrix (spots x cell types)
    """
    return model.get_proportions()


def get_gamma_values(model: DestVI) -> Dict[str, np.ndarray]:
    """
    Get intra-cell-type variation (gamma values).

    Returns
    -------
    Dict[str, np.ndarray]
        Dictionary mapping cell type -> gamma array
    """
    return model.get_gamma()


def get_celltype_expression(
    model: DestVI,
    cell_type: str,
    st_adata,
    proportion_threshold: float = 0.1,
) -> "pd.DataFrame":
    """
    Get cell-type-specific expression for spots.

    Parameters
    ----------
    model : DestVI
        Trained spatial model
    cell_type : str
        Cell type name
    st_adata : AnnData
        Spatial data
    proportion_threshold : float
        Minimum proportion to include spot

    Returns
    -------
    pd.DataFrame
        Expression matrix for qualifying spots
    """
    proportions = model.get_proportions()

    if f"prop_{cell_type}" in st_adata.obs.columns:
        mask = st_adata.obs[f"prop_{cell_type}"] > proportion_threshold
    else:
        mask = proportions[cell_type] > proportion_threshold

    indices = np.where(mask)[0]
    print(f"Getting expression for {len(indices)} spots with {cell_type} > {proportion_threshold}")

    return model.get_scale_for_ct(cell_type, indices=indices)


def filter_shared_genes(sc_adata, st_adata):
    """
    Filter both datasets to shared genes.

    Returns
    -------
    tuple
        (filtered_sc_adata, filtered_st_adata)
    """
    shared_genes = np.intersect1d(sc_adata.var_names, st_adata.var_names)
    print(f"Shared genes: {len(shared_genes)}")

    sc_filtered = sc_adata[:, shared_genes].copy()
    st_filtered = st_adata[:, shared_genes].copy()

    return sc_filtered, st_filtered


def select_hvg_genes(sc_adata, n_top_genes: int = 3000, layer: str = "counts"):
    """
    Select highly variable genes from reference.

    Parameters
    ----------
    sc_adata : AnnData
        Reference data
    n_top_genes : int
        Number of HVGs to select
    layer : str
        Layer with counts

    Returns
    -------
    list
        HVG names
    """
    import scanpy as sc

    sc.pp.highly_variable_genes(
        sc_adata,
        n_top_genes=n_top_genes,
        subset=False,
        layer=layer,
        flavor="seurat_v3",
    )

    hvg = sc_adata.var_names[sc_adata.var.highly_variable].tolist()
    print(f"Selected {len(hvg)} highly variable genes")
    return hvg


def add_proportions_to_obs(st_adata, proportions: "pd.DataFrame"):
    """
    Add proportion columns to st_adata.obs for easy plotting.
    """
    st_adata.obsm["proportions"] = proportions
    for ct in proportions.columns:
        st_adata.obs[f"prop_{ct}"] = proportions[ct].values
    print(f"Added proportions for {len(proportions.columns)} cell types")


def add_gamma_to_obsm(st_adata, gamma_dict: Dict[str, np.ndarray]):
    """
    Add gamma values to st_adata.obsm.
    """
    for ct, gamma in gamma_dict.items():
        st_adata.obsm[f"{ct}_gamma"] = gamma
        st_adata.obs[f"{ct}_gamma_PC1"] = gamma[:, 0]
    print(f"Added gamma values for {len(gamma_dict)} cell types")


def save_reference_model(model: CondSCVI, path: str, overwrite: bool = True):
    """Save reference model."""
    model.save(path, overwrite=overwrite)
    print(f"Reference model saved to: {path}")


def save_spatial_model(model: DestVI, path: str, overwrite: bool = True):
    """Save spatial model."""
    model.save(path, overwrite=overwrite)
    print(f"Spatial model saved to: {path}")


def load_reference_model(path: str, adata) -> CondSCVI:
    """Load reference model."""
    return CondSCVI.load(path, adata=adata)


def load_spatial_model(path: str, adata) -> DestVI:
    """Load spatial model."""
    return DestVI.load(path, adata=adata)
