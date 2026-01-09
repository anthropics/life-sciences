"""
Core functions for scVI/scANVI label transfer workflow.

This module provides modular utility functions for integrating datasets
and transferring cell type labels from reference to query.
"""

import scvi
import scanpy as sc
import anndata
import numpy as np
import pandas as pd
from typing import Optional, Dict, Tuple
from sklearn.metrics import accuracy_score, f1_score, classification_report


def normalize_gene_length(
    adata,
    gene_length_file: Optional[str] = None,
    gene_lengths: Optional[pd.DataFrame] = None
):
    """
    Normalize counts by gene length for SmartSeq2 data.

    SmartSeq2 read counts are proportional to gene length (unlike UMI-based
    methods). This normalization corrects for that bias.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with raw counts.
    gene_length_file : str, optional
        Path to file with gene lengths (gene name in first column, length in second).
    gene_lengths : pd.DataFrame, optional
        DataFrame with gene lengths indexed by gene name.

    Returns
    -------
    AnnData
        Copy of adata with normalized counts.
    """
    adata = adata.copy()

    if gene_lengths is None and gene_length_file is not None:
        gene_lengths = pd.read_csv(
            gene_length_file,
            delimiter=None,  # Auto-detect
            header=None,
            index_col=0
        )
    elif gene_lengths is None:
        # Try to fetch default gene lengths
        try:
            gene_lengths = pd.read_csv(
                "https://raw.githubusercontent.com/chenlingantelope/HarmonizationSCANVI/master/data/gene_len.txt",
                delimiter=" ",
                header=None,
                index_col=0
            )
            print("  Using default gene length file")
        except Exception as e:
            raise ValueError(f"Could not load gene lengths: {e}")

    # Filter to genes present in data
    common_genes = gene_lengths.index.intersection(adata.var_names)
    adata = adata[:, common_genes].copy()
    gene_lengths = gene_lengths.loc[common_genes]

    # Normalize
    median_length = np.median(gene_lengths.iloc[:, 0].values)
    adata.X = adata.X / gene_lengths.iloc[:, 0].values * median_length
    adata.X = np.rint(adata.X)  # Round to integers

    print(f"  Normalized {len(common_genes)} genes by length")

    return adata


def concatenate_datasets(
    reference,
    query,
    labels_key: str,
    batch_key: str = 'dataset',
    unlabeled: str = 'Unknown'
) -> anndata.AnnData:
    """
    Concatenate reference and query datasets.

    Parameters
    ----------
    reference : AnnData
        Annotated reference dataset with cell type labels.
    query : AnnData
        Unannotated query dataset.
    labels_key : str
        Column in reference.obs containing cell type labels.
    batch_key : str, default='dataset'
        Column name for batch information.
    unlabeled : str, default='Unknown'
        Label to assign to query cells.

    Returns
    -------
    AnnData
        Concatenated dataset.
    """
    # Ensure query has labels column (as Unknown)
    if labels_key not in query.obs.columns:
        query.obs[labels_key] = unlabeled

    # Concatenate
    adata = anndata.concat([reference, query])

    # Make var_names unique
    adata.var_names_make_unique()

    print(f"  Concatenated: {adata.n_obs} cells × {adata.n_vars} genes")

    return adata


def select_hvg_across_batches(
    adata,
    batch_key: str,
    n_top_genes: int = 2000,
    flavor: str = 'seurat_v3'
) -> None:
    """
    Select highly variable genes accounting for batch structure.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    batch_key : str
        Column for batch information.
    n_top_genes : int, default=2000
        Number of HVGs to select.
    flavor : str, default='seurat_v3'
        HVG selection method.

    Returns
    -------
    None
        Modifies adata in place.
    """
    # Normalize for HVG selection only
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)

    # Select HVGs with batch correction
    sc.pp.highly_variable_genes(
        adata_norm,
        flavor=flavor,
        n_top_genes=n_top_genes,
        batch_key=batch_key
    )

    # Transfer HVG info to original adata
    adata.var['highly_variable'] = adata_norm.var['highly_variable']

    # Subset to HVGs
    adata._inplace_subset_var(adata.var['highly_variable'])

    print(f"  Selected {adata.n_vars} highly variable genes")


def setup_combined_anndata(
    adata,
    batch_key: str,
    labels_key: str,
    layer: str = 'counts'
) -> None:
    """
    Setup AnnData for scVI/scANVI training.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    batch_key : str
        Column for batch information.
    labels_key : str
        Column for cell type labels.
    layer : str, default='counts'
        Layer containing count data.

    Returns
    -------
    None
        Registers adata with scvi-tools.
    """
    print(f"  Setting up AnnData for scVI...")
    print(f"    Batch key: {batch_key}")
    print(f"    Labels key: {labels_key}")

    scvi.model.SCVI.setup_anndata(
        adata,
        layer=layer,
        batch_key=batch_key,
        labels_key=labels_key
    )


def train_scvi_integration(
    adata,
    n_latent: int = 30,
    n_layers: int = 2,
    max_epochs: int = 400,
    early_stopping: bool = True,
    **kwargs
):
    """
    Train scVI model for dataset integration.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix (must be set up first).
    n_latent : int, default=30
        Latent space dimensions.
    n_layers : int, default=2
        Number of hidden layers.
    max_epochs : int, default=400
        Maximum training epochs.
    early_stopping : bool, default=True
        Whether to use early stopping.
    **kwargs
        Additional arguments to SCVI.

    Returns
    -------
    scvi.model.SCVI
        Trained scVI model.
    """
    print(f"  Initializing scVI...")
    print(f"    n_latent={n_latent}, n_layers={n_layers}")

    model = scvi.model.SCVI(
        adata,
        n_latent=n_latent,
        n_layers=n_layers,
        **kwargs
    )

    print(f"  Training scVI ({max_epochs} max epochs)...")
    model.train(max_epochs=max_epochs, early_stopping=early_stopping)

    print(f"  Training complete!")
    return model


def train_scanvi_transfer(
    scvi_model,
    adata,
    labels_key: str,
    unlabeled_category: str = 'Unknown',
    max_epochs: int = 20,
    n_samples_per_label: int = 100,
    **kwargs
):
    """
    Train scANVI model for label transfer.

    Parameters
    ----------
    scvi_model : scvi.model.SCVI
        Pre-trained scVI model.
    adata : AnnData
        Annotated data matrix.
    labels_key : str
        Column with cell type labels.
    unlabeled_category : str, default='Unknown'
        Label for query cells.
    max_epochs : int, default=20
        Maximum training epochs.
    n_samples_per_label : int, default=100
        Samples per label during training.
    **kwargs
        Additional arguments to train.

    Returns
    -------
    scvi.model.SCANVI
        Trained scANVI model.
    """
    print(f"  Initializing scANVI from scVI...")

    # Count labeled vs unlabeled
    n_labeled = (adata.obs[labels_key] != unlabeled_category).sum()
    n_unlabeled = (adata.obs[labels_key] == unlabeled_category).sum()
    print(f"    Labeled cells: {n_labeled}")
    print(f"    Unlabeled cells: {n_unlabeled}")

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        adata=adata,
        unlabeled_category=unlabeled_category,
        labels_key=labels_key
    )

    print(f"  Training scANVI ({max_epochs} epochs)...")
    scanvi_model.train(
        max_epochs=max_epochs,
        n_samples_per_label=n_samples_per_label,
        **kwargs
    )

    print(f"  Training complete!")
    return scanvi_model


def predict_labels(scanvi_model, adata=None) -> np.ndarray:
    """
    Predict cell type labels.

    Parameters
    ----------
    scanvi_model : scvi.model.SCANVI
        Trained scANVI model.
    adata : AnnData, optional
        Data to predict labels for.

    Returns
    -------
    np.ndarray
        Predicted labels.
    """
    return scanvi_model.predict(adata=adata)


def get_prediction_probabilities(scanvi_model, adata=None) -> np.ndarray:
    """
    Get prediction probability distribution.

    Parameters
    ----------
    scanvi_model : scvi.model.SCANVI
        Trained scANVI model.
    adata : AnnData, optional
        Data to get probabilities for.

    Returns
    -------
    np.ndarray
        Probability matrix (n_cells × n_classes).
    """
    return scanvi_model.predict(adata=adata, soft=True)


def evaluate_transfer(
    true_labels,
    predicted_labels
) -> Dict[str, float]:
    """
    Evaluate label transfer accuracy.

    Parameters
    ----------
    true_labels : array-like
        Ground truth labels.
    predicted_labels : array-like
        Predicted labels.

    Returns
    -------
    dict
        Dictionary with accuracy metrics.
    """
    # Convert to arrays
    true_labels = np.asarray(true_labels)
    predicted_labels = np.asarray(predicted_labels)

    # Filter out unknowns
    mask = true_labels != 'Unknown'
    true_labels = true_labels[mask]
    predicted_labels = predicted_labels[mask]

    metrics = {
        'accuracy': accuracy_score(true_labels, predicted_labels),
        'f1_weighted': f1_score(true_labels, predicted_labels, average='weighted'),
        'f1_macro': f1_score(true_labels, predicted_labels, average='macro')
    }

    return metrics


def get_confusion_matrix(
    true_labels,
    predicted_labels,
    normalize: str = 'true'
) -> pd.DataFrame:
    """
    Compute confusion matrix as DataFrame.

    Parameters
    ----------
    true_labels : array-like
        Ground truth labels.
    predicted_labels : array-like
        Predicted labels.
    normalize : str, default='true'
        Normalization: 'true', 'pred', 'all', or None.

    Returns
    -------
    pd.DataFrame
        Confusion matrix.
    """
    return pd.crosstab(
        true_labels,
        predicted_labels,
        normalize='index' if normalize == 'true' else normalize
    )


def save_model(model, path: str) -> None:
    """
    Save a trained model.

    Parameters
    ----------
    model : scvi model
        Trained scVI or scANVI model.
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
        Path to saved model.
    adata : AnnData
        AnnData object.

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
        Path to saved model.
    adata : AnnData
        AnnData object.

    Returns
    -------
    scvi.model.SCANVI
        Loaded model.
    """
    return scvi.model.SCANVI.load(path, adata=adata)


def filter_by_confidence(
    adata,
    predictions: np.ndarray,
    probabilities: np.ndarray,
    threshold: float = 0.7,
    low_conf_label: str = 'Low_confidence'
) -> np.ndarray:
    """
    Filter predictions by confidence threshold.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    predictions : np.ndarray
        Predicted labels.
    probabilities : np.ndarray
        Prediction probabilities.
    threshold : float, default=0.7
        Confidence threshold.
    low_conf_label : str, default='Low_confidence'
        Label for low-confidence predictions.

    Returns
    -------
    np.ndarray
        Filtered predictions.
    """
    confidence = probabilities.max(axis=1)
    filtered = predictions.copy()
    filtered[confidence < threshold] = low_conf_label

    n_low = (confidence < threshold).sum()
    print(f"  {n_low} cells below confidence threshold {threshold}")

    return filtered
