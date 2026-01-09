"""
Core integration functions for scVI/scANVI workflows.

This module provides modular utility functions for single-cell data integration
using scvi-tools. Functions can be used independently or as part of the
complete integration_analysis.py pipeline.
"""

import scvi
import scanpy as sc
import numpy as np
from typing import Optional, Literal


def setup_anndata_scvi(
    adata,
    batch_key: str,
    layer: Optional[str] = None,
    labels_key: Optional[str] = None,
    categorical_covariate_keys: Optional[list] = None,
    continuous_covariate_keys: Optional[list] = None
) -> None:
    """
    Register AnnData object with scvi-tools.

    This function prepares the AnnData object for use with scVI/scANVI models
    by registering the count data, batch information, and optional covariates.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with raw counts.
    batch_key : str
        Key in adata.obs for batch information.
    layer : str, optional
        Key in adata.layers containing counts. If None, uses adata.X.
    labels_key : str, optional
        Key in adata.obs for cell type labels (required for scANVI).
    categorical_covariate_keys : list, optional
        Keys for additional categorical covariates to include.
    continuous_covariate_keys : list, optional
        Keys for continuous covariates to include.

    Returns
    -------
    None
        Modifies adata in place with scvi-tools registration.
    """
    print(f"  Registering AnnData with scvi-tools...")
    print(f"    Batch key: {batch_key}")
    print(f"    Layer: {layer or 'X'}")
    if labels_key:
        print(f"    Labels key: {labels_key}")

    scvi.model.SCVI.setup_anndata(
        adata,
        layer=layer,
        batch_key=batch_key,
        labels_key=labels_key,
        categorical_covariate_keys=categorical_covariate_keys,
        continuous_covariate_keys=continuous_covariate_keys
    )
    print("  AnnData registered successfully")


def train_scvi_model(
    adata,
    n_layers: int = 2,
    n_latent: int = 30,
    gene_likelihood: Literal['nb', 'zinb', 'poisson'] = 'nb',
    max_epochs: Optional[int] = None,
    early_stopping: bool = True,
    batch_size: int = 128,
    plan_kwargs: Optional[dict] = None,
    **model_kwargs
):
    """
    Train an scVI model for unsupervised integration.

    scVI uses a variational autoencoder to learn a batch-corrected latent
    representation of single-cell gene expression data.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix (must be registered with setup_anndata_scvi first).
    n_layers : int, default=2
        Number of hidden layers in encoder and decoder networks.
    n_latent : int, default=30
        Dimensionality of the latent space.
    gene_likelihood : {'nb', 'zinb', 'poisson'}, default='nb'
        Distribution to model gene expression:
        - 'nb': Negative binomial (recommended for most data)
        - 'zinb': Zero-inflated negative binomial (for highly sparse data)
        - 'poisson': Poisson distribution
    max_epochs : int, optional
        Maximum number of training epochs. If None, automatically determined.
    early_stopping : bool, default=True
        Whether to use early stopping based on validation loss.
    batch_size : int, default=128
        Minibatch size for training.
    plan_kwargs : dict, optional
        Additional kwargs for training plan (e.g., learning rate).
    **model_kwargs
        Additional kwargs passed to SCVI model constructor.

    Returns
    -------
    scvi.model.SCVI
        Trained scVI model.
    """
    print(f"  Initializing scVI model...")
    print(f"    Architecture: {n_layers} layers, {n_latent} latent dims")
    print(f"    Gene likelihood: {gene_likelihood}")

    model = scvi.model.SCVI(
        adata,
        n_layers=n_layers,
        n_latent=n_latent,
        gene_likelihood=gene_likelihood,
        **model_kwargs
    )

    print(f"  Training scVI model...")
    train_kwargs = {
        'early_stopping': early_stopping,
        'batch_size': batch_size,
    }
    if max_epochs is not None:
        train_kwargs['max_epochs'] = max_epochs
    if plan_kwargs is not None:
        train_kwargs['plan_kwargs'] = plan_kwargs

    model.train(**train_kwargs)

    print(f"  Training complete!")
    print(f"    Final training loss: {model.history['elbo_train'].iloc[-1]:.4f}")
    if 'elbo_validation' in model.history:
        print(f"    Final validation loss: {model.history['elbo_validation'].iloc[-1]:.4f}")

    return model


def train_scanvi_model(
    scvi_model,
    adata,
    labels_key: str,
    unlabeled_category: str = 'Unknown',
    max_epochs: int = 20,
    n_samples_per_label: Optional[int] = None,
    **train_kwargs
):
    """
    Train an scANVI model for semi-supervised integration.

    scANVI extends scVI by incorporating cell type annotations, which helps
    preserve biological signal during integration. It can also predict labels
    for unlabeled cells.

    Parameters
    ----------
    scvi_model : scvi.model.SCVI
        Pre-trained scVI model to initialize from.
    adata : AnnData
        Annotated data matrix (same as used for scVI).
    labels_key : str
        Key in adata.obs containing cell type labels.
    unlabeled_category : str, default='Unknown'
        Label used for cells without annotations (will be predicted).
    max_epochs : int, default=20
        Maximum training epochs (typically fewer needed than scVI).
    n_samples_per_label : int, optional
        Number of cells per label to sample during training.
    **train_kwargs
        Additional kwargs passed to train method.

    Returns
    -------
    scvi.model.SCANVI
        Trained scANVI model.
    """
    print(f"  Initializing scANVI from pre-trained scVI...")
    print(f"    Labels key: {labels_key}")
    print(f"    Unlabeled category: {unlabeled_category}")

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        adata=adata,
        labels_key=labels_key,
        unlabeled_category=unlabeled_category
    )

    # Count labeled vs unlabeled
    n_labeled = (adata.obs[labels_key] != unlabeled_category).sum()
    n_unlabeled = (adata.obs[labels_key] == unlabeled_category).sum()
    print(f"    Labeled cells: {n_labeled}")
    print(f"    Unlabeled cells: {n_unlabeled}")

    print(f"  Training scANVI model (max {max_epochs} epochs)...")
    train_params = {'max_epochs': max_epochs}
    if n_samples_per_label is not None:
        train_params['n_samples_per_label'] = n_samples_per_label
    train_params.update(train_kwargs)

    scanvi_model.train(**train_params)

    print(f"  Training complete!")

    return scanvi_model


def get_latent_representation(model, adata=None) -> np.ndarray:
    """
    Extract latent representation from a trained scVI/scANVI model.

    Parameters
    ----------
    model : scvi.model.SCVI or scvi.model.SCANVI
        Trained model.
    adata : AnnData, optional
        Data to get representation for. If None, uses training data.

    Returns
    -------
    np.ndarray
        Latent representation matrix (n_cells × n_latent).
    """
    return model.get_latent_representation(adata=adata)


def predict_labels(scanvi_model, adata=None, soft: bool = False):
    """
    Predict cell type labels using a trained scANVI model.

    Parameters
    ----------
    scanvi_model : scvi.model.SCANVI
        Trained scANVI model.
    adata : AnnData, optional
        Data to predict labels for. If None, uses training data.
    soft : bool, default=False
        If True, return probability distribution over labels.

    Returns
    -------
    np.ndarray or pd.DataFrame
        Predicted labels (hard) or probabilities (soft).
    """
    if soft:
        return scanvi_model.predict(adata=adata, soft=True)
    return scanvi_model.predict(adata=adata)


def compute_neighbors_and_umap(
    adata,
    use_rep: str,
    n_neighbors: int = 15,
    min_dist: float = 0.3,
    spread: float = 1.0,
    key_added: Optional[str] = None
) -> None:
    """
    Compute k-nearest neighbor graph and UMAP embedding.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with latent representation.
    use_rep : str
        Key in adata.obsm containing the representation to use.
    n_neighbors : int, default=15
        Number of neighbors for kNN graph.
    min_dist : float, default=0.3
        UMAP minimum distance parameter.
    spread : float, default=1.0
        UMAP spread parameter.
    key_added : str, optional
        Suffix for keys added to adata. If None, uses default keys.

    Returns
    -------
    None
        Modifies adata in place.
    """
    neighbors_key = f"{key_added}_neighbors" if key_added else 'neighbors'
    umap_key = f"X_umap_{key_added}" if key_added else 'X_umap'

    print(f"    Computing neighbors (n={n_neighbors}) using {use_rep}...")
    sc.pp.neighbors(
        adata,
        use_rep=use_rep,
        n_neighbors=n_neighbors,
        key_added=neighbors_key if key_added else None
    )

    print(f"    Computing UMAP (min_dist={min_dist})...")
    sc.tl.umap(
        adata,
        min_dist=min_dist,
        spread=spread,
        neighbors_key=neighbors_key if key_added else None
    )

    # Rename UMAP coordinates if key_added specified
    if key_added:
        adata.obsm[umap_key] = adata.obsm['X_umap'].copy()


def cluster_cells(
    adata,
    resolution: float = 1.0,
    neighbors_key: Optional[str] = None,
    key_added: str = 'leiden'
) -> None:
    """
    Cluster cells using Leiden algorithm.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with neighbor graph.
    resolution : float, default=1.0
        Resolution parameter for Leiden clustering.
    neighbors_key : str, optional
        Key for neighbor graph in adata.obsp.
    key_added : str, default='leiden'
        Key to add cluster labels to adata.obs.

    Returns
    -------
    None
        Modifies adata in place.
    """
    print(f"    Clustering cells (resolution={resolution})...")
    sc.tl.leiden(
        adata,
        resolution=resolution,
        neighbors_key=neighbors_key,
        key_added=key_added
    )
    n_clusters = adata.obs[key_added].nunique()
    print(f"    Found {n_clusters} clusters")


def save_model(model, path: str) -> None:
    """
    Save a trained scVI/scANVI model.

    Parameters
    ----------
    model : scvi.model.SCVI or scvi.model.SCANVI
        Trained model to save.
    path : str
        Directory path to save model.
    """
    model.save(path, overwrite=True)


def load_scvi_model(path: str, adata):
    """
    Load a saved scVI model.

    Parameters
    ----------
    path : str
        Directory path containing saved model.
    adata : AnnData
        AnnData object (must match the data used for training).

    Returns
    -------
    scvi.model.SCVI
        Loaded model.
    """
    return scvi.model.SCVI.load(path, adata=adata)


def load_scanvi_model(path: str, adata):
    """
    Load a saved scANVI model.

    Parameters
    ----------
    path : str
        Directory path containing saved model.
    adata : AnnData
        AnnData object (must match the data used for training).

    Returns
    -------
    scvi.model.SCANVI
        Loaded model.
    """
    return scvi.model.SCANVI.load(path, adata=adata)


def prepare_query_data(query_adata, reference_model) -> None:
    """
    Prepare query data for mapping to a reference model.

    This aligns the query data with the reference model's features and
    registers it appropriately.

    Parameters
    ----------
    query_adata : AnnData
        Query data to prepare.
    reference_model : scvi.model.SCVI or scvi.model.SCANVI
        Reference model to map to.

    Returns
    -------
    None
        Modifies query_adata in place.
    """
    scvi.model.SCVI.prepare_query_anndata(query_adata, reference_model)


def transfer_labels_to_query(
    reference_model,
    query_adata,
    max_epochs: int = 100,
    weight_decay: float = 0.0
):
    """
    Transfer cell type labels from reference to query data.

    Uses scANVI's query-to-reference mapping to predict labels for query cells.

    Parameters
    ----------
    reference_model : scvi.model.SCANVI
        Trained scANVI model on reference data.
    query_adata : AnnData
        Query data (must be prepared with prepare_query_data).
    max_epochs : int, default=100
        Training epochs for query model.
    weight_decay : float, default=0.0
        Weight decay for query training.

    Returns
    -------
    tuple
        (query_model, predictions, latent_representation)
    """
    print("  Loading query data into reference model...")
    query_model = reference_model.load_query_data(query_adata)

    print(f"  Training query model ({max_epochs} epochs)...")
    query_model.train(
        max_epochs=max_epochs,
        plan_kwargs={'weight_decay': weight_decay}
    )

    predictions = query_model.predict()
    latent = query_model.get_latent_representation()

    return query_model, predictions, latent
