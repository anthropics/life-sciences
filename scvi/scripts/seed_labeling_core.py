"""
Core functions for scANVI seed labeling workflow.

This module provides modular utility functions for semi-supervised cell
annotation using marker gene signatures and scANVI.
"""

import scvi
import scanpy as sc
import numpy as np
import pandas as pd
import json
from typing import Dict, List, Optional, Union
from pathlib import Path


def load_marker_genes(path: str) -> Dict[str, Dict[str, List[str]]]:
    """
    Load marker gene signatures from a JSON or CSV file.

    Parameters
    ----------
    path : str
        Path to marker gene file.
        JSON format: {"cell_type": {"positive": [...], "negative": [...]}}
        CSV format: cell_type,gene,direction columns

    Returns
    -------
    dict
        Dictionary mapping cell types to positive/negative gene lists.
    """
    path = Path(path)

    if path.suffix == '.json':
        with open(path) as f:
            markers = json.load(f)
    elif path.suffix == '.csv':
        df = pd.read_csv(path)
        markers = {}
        for ct in df['cell_type'].unique():
            ct_df = df[df['cell_type'] == ct]
            markers[ct] = {
                'positive': ct_df[ct_df['direction'] == 'positive']['gene'].tolist(),
                'negative': ct_df[ct_df['direction'] == 'negative']['gene'].tolist()
            }
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")

    return markers


def score_cells_by_markers(
    adata,
    markers: Dict[str, Dict[str, List[str]]],
    layer: Optional[str] = None
) -> Dict[str, np.ndarray]:
    """
    Score cells based on marker gene expression.

    Computes a score for each cell type based on expression of positive
    markers minus expression of negative markers.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with raw counts.
    markers : dict
        Dictionary mapping cell types to positive/negative gene lists.
    layer : str, optional
        Layer containing counts. If None, uses adata.X.

    Returns
    -------
    dict
        Dictionary mapping cell types to score arrays.
    """
    print("  Normalizing data for scoring...")

    # Create normalized copy for scoring
    normalized = adata.copy()
    if layer:
        normalized.X = normalized.layers[layer].copy()

    sc.pp.normalize_total(normalized, target_sum=1e4)
    sc.pp.log1p(normalized)

    # Get all marker genes
    all_markers = set()
    for ct, genes in markers.items():
        all_markers.update(genes.get('positive', []))
        all_markers.update(genes.get('negative', []))

    # Filter to present genes
    present_markers = all_markers & set(normalized.var_names)
    normalized = normalized[:, list(present_markers)].copy()

    # Scale for scoring
    sc.pp.scale(normalized)

    # Compute scores
    scores = {}
    for ct, gene_sets in markers.items():
        print(f"  Scoring {ct}...")
        score = np.zeros(normalized.n_obs)

        # Add positive markers
        for gene in gene_sets.get('positive', []):
            if gene in normalized.var_names:
                expr = np.array(normalized[:, gene].X).flatten()
                score += expr

        # Subtract negative markers
        for gene in gene_sets.get('negative', []):
            if gene in normalized.var_names:
                expr = np.array(normalized[:, gene].X).flatten()
                score -= expr

        scores[ct] = score

    return scores


def select_seed_cells(
    scores: Dict[str, np.ndarray],
    n_cells: int = 50,
    min_percentile: float = 95,
    method: str = 'top_n'
) -> Dict[str, np.ndarray]:
    """
    Select seed cells based on marker scores.

    Parameters
    ----------
    scores : dict
        Dictionary mapping cell types to score arrays.
    n_cells : int, default=50
        Number of cells to select per type.
    min_percentile : float, default=95
        Minimum score percentile for selection.
    method : str, default='top_n'
        Selection method: 'top_n' or 'percentile'.

    Returns
    -------
    dict
        Dictionary mapping cell types to boolean mask arrays.
    """
    seed_masks = {}
    n_total = len(list(scores.values())[0])

    for ct, score in scores.items():
        mask = np.zeros(n_total, dtype=bool)

        if method == 'top_n':
            # Select top N cells by score
            top_idx = np.argsort(score)[-n_cells:]
            mask[top_idx] = True
        elif method == 'percentile':
            # Select cells above percentile threshold
            threshold = np.percentile(score, min_percentile)
            mask = score >= threshold

        # Apply minimum percentile filter for top_n method
        if method == 'top_n':
            threshold = np.percentile(score, min_percentile)
            mask = mask & (score >= threshold)

        seed_masks[ct] = mask
        print(f"  {ct}: {mask.sum()} seed cells selected")

    return seed_masks


def create_seed_labels(
    adata,
    seed_masks: Dict[str, np.ndarray],
    cell_types: List[str],
    key_added: str = 'seed_labels',
    unlabeled: str = 'Unknown'
) -> None:
    """
    Create seed label column in AnnData.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    seed_masks : dict
        Dictionary mapping cell types to boolean masks.
    cell_types : list
        List of cell type names.
    key_added : str, default='seed_labels'
        Key to add to adata.obs.
    unlabeled : str, default='Unknown'
        Label for unlabeled cells.

    Returns
    -------
    None
        Modifies adata in place.
    """
    # Initialize all as unlabeled
    labels = np.array([unlabeled] * adata.n_obs)

    # Assign seed labels
    for ct in cell_types:
        if ct in seed_masks:
            labels[seed_masks[ct]] = ct

    adata.obs[key_added] = pd.Categorical(labels)


def train_scvi_model(
    adata,
    n_latent: int = 30,
    n_layers: int = 2,
    max_epochs: int = 100,
    early_stopping: bool = True,
    **kwargs
):
    """
    Train scVI base model.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix (must be set up with scvi.model.SCVI.setup_anndata).
    n_latent : int, default=30
        Latent space dimensions.
    n_layers : int, default=2
        Number of hidden layers.
    max_epochs : int, default=100
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
    print(f"  Initializing scVI model...")
    print(f"    n_latent={n_latent}, n_layers={n_layers}")

    model = scvi.model.SCVI(
        adata,
        n_latent=n_latent,
        n_layers=n_layers,
        **kwargs
    )

    print(f"  Training scVI ({max_epochs} epochs)...")
    model.train(max_epochs=max_epochs, early_stopping=early_stopping)

    print(f"  Training complete!")
    return model


def train_scanvi_from_seeds(
    scvi_model,
    adata,
    labels_key: str = 'seed_labels',
    unlabeled_category: str = 'Unknown',
    max_epochs: int = 25,
    **kwargs
):
    """
    Train scANVI model from pre-trained scVI using seed labels.

    Parameters
    ----------
    scvi_model : scvi.model.SCVI
        Pre-trained scVI model.
    adata : AnnData
        Annotated data matrix.
    labels_key : str, default='seed_labels'
        Key in adata.obs containing seed labels.
    unlabeled_category : str, default='Unknown'
        Label for unlabeled cells.
    max_epochs : int, default=25
        Maximum training epochs.
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
        unlabeled_category=unlabeled_category,
        labels_key=labels_key
    )

    print(f"  Training scANVI ({max_epochs} epochs)...")
    scanvi_model.train(max_epochs=max_epochs, **kwargs)

    print(f"  Training complete!")
    return scanvi_model


def predict_labels(
    scanvi_model,
    adata=None
) -> np.ndarray:
    """
    Predict cell type labels using trained scANVI model.

    Parameters
    ----------
    scanvi_model : scvi.model.SCANVI
        Trained scANVI model.
    adata : AnnData, optional
        Data to predict labels for. If None, uses training data.

    Returns
    -------
    np.ndarray
        Predicted labels for each cell.
    """
    return scanvi_model.predict(adata=adata)


def get_prediction_confidence(
    scanvi_model,
    adata=None
) -> np.ndarray:
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


def get_latent_representation(
    model,
    adata=None
) -> np.ndarray:
    """
    Get latent representation from scVI or scANVI model.

    Parameters
    ----------
    model : scvi.model.SCVI or scvi.model.SCANVI
        Trained model.
    adata : AnnData, optional
        Data to get representation for.

    Returns
    -------
    np.ndarray
        Latent representation (n_cells × n_latent).
    """
    return model.get_latent_representation(adata=adata)


def save_model(model, path: str) -> None:
    """
    Save a trained scVI/scANVI model.

    Parameters
    ----------
    model : scvi.model.SCVI or scvi.model.SCANVI
        Trained model.
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


def validate_markers_in_data(
    adata,
    markers: Dict[str, Dict[str, List[str]]]
) -> Dict[str, Dict[str, List[str]]]:
    """
    Validate and filter markers to those present in the data.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    markers : dict
        Dictionary of marker gene signatures.

    Returns
    -------
    dict
        Filtered markers containing only genes present in data.
    """
    gene_set = set(adata.var_names)
    filtered = {}

    for ct, genes in markers.items():
        filtered[ct] = {
            'positive': [g for g in genes.get('positive', []) if g in gene_set],
            'negative': [g for g in genes.get('negative', []) if g in gene_set]
        }

        missing_pos = set(genes.get('positive', [])) - gene_set
        missing_neg = set(genes.get('negative', [])) - gene_set

        if missing_pos:
            print(f"  Warning: {ct} missing positive markers: {missing_pos}")
        if missing_neg:
            print(f"  Warning: {ct} missing negative markers: {missing_neg}")

    return filtered


def refine_seeds_by_confidence(
    adata,
    predictions: np.ndarray,
    confidence: np.ndarray,
    threshold: float = 0.9,
    key_added: str = 'refined_seeds'
) -> None:
    """
    Create refined seed labels from high-confidence predictions.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    predictions : np.ndarray
        Predicted labels.
    confidence : np.ndarray
        Confidence matrix (n_cells × n_classes).
    threshold : float, default=0.9
        Minimum confidence for inclusion.
    key_added : str, default='refined_seeds'
        Key to add to adata.obs.

    Returns
    -------
    None
        Modifies adata in place.
    """
    max_conf = confidence.max(axis=1)
    high_conf_mask = max_conf >= threshold

    labels = np.array(['Unknown'] * adata.n_obs)
    labels[high_conf_mask] = predictions[high_conf_mask]

    adata.obs[key_added] = pd.Categorical(labels)

    n_refined = high_conf_mask.sum()
    print(f"  Created refined seeds: {n_refined} cells (>{threshold} confidence)")
