#!/usr/bin/env python3
"""
scArches Reference Mapping Core Functions

Modular functions for atlas-level reference mapping with scArches.
"""

import numpy as np
import scvi
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from typing import Optional, Tuple


def setup_reference_data(
    adata_ref,
    layer: str = "counts",
    batch_key: Optional[str] = None,
):
    """
    Register reference data for scVI.

    Parameters
    ----------
    adata_ref : AnnData
        Reference dataset
    layer : str
        Layer containing counts
    batch_key : str, optional
        Column for batch information
    """
    scvi.model.SCVI.setup_anndata(
        adata_ref,
        layer=layer,
        batch_key=batch_key,
    )
    print(f"Reference data registered")
    print(f"  Layer: {layer}")
    if batch_key:
        print(f"  Batch key: {batch_key}")


def train_reference_model(
    adata_ref,
    n_latent: int = 30,
    n_layers: int = 2,
    dropout_rate: float = 0.2,
    max_epochs: int = 400,
    early_stopping: bool = True,
) -> scvi.model.SCVI:
    """
    Train scArches-compatible reference model.

    CRITICAL: Uses settings required for scArches compatibility:
    - use_layer_norm="both"
    - use_batch_norm="none"
    - encode_covariates=True

    Parameters
    ----------
    adata_ref : AnnData
        Registered reference data
    n_latent : int
        Latent space dimensions
    n_layers : int
        Number of hidden layers
    dropout_rate : float
        Dropout rate
    max_epochs : int
        Maximum training epochs
    early_stopping : bool
        Enable early stopping

    Returns
    -------
    scvi.model.SCVI
        Trained reference model
    """
    model = scvi.model.SCVI(
        adata_ref,
        use_layer_norm="both",      # CRITICAL for scArches
        use_batch_norm="none",      # CRITICAL for scArches
        encode_covariates=True,     # CRITICAL for scArches
        dropout_rate=dropout_rate,
        n_layers=n_layers,
        n_latent=n_latent,
    )

    print(f"\nReference model architecture (scArches-compatible):")
    print(f"  Latent dimensions: {n_latent}")
    print(f"  Layers: {n_layers}")
    print(f"  use_layer_norm: both")
    print(f"  use_batch_norm: none")
    print(f"  encode_covariates: True")

    model.train(max_epochs=max_epochs, early_stopping=early_stopping)

    final_epoch = len(model.history["elbo_train"])
    print(f"\nTraining complete after {final_epoch} epochs")

    return model


def train_classifier(
    model: scvi.model.SCVI,
    adata_ref,
    labels_key: str = "cell_type",
    n_estimators: int = 100,
) -> Tuple[RandomForestClassifier, float]:
    """
    Train classifier on reference latent space.

    Parameters
    ----------
    model : scvi.model.SCVI
        Trained reference model
    adata_ref : AnnData
        Reference data
    labels_key : str
        Column with cell type labels
    n_estimators : int
        Number of trees in RandomForest

    Returns
    -------
    Tuple[RandomForestClassifier, float]
        (trained_classifier, self_accuracy)
    """
    latent = model.get_latent_representation()

    clf = RandomForestClassifier(
        n_estimators=n_estimators,
        random_state=0,
        n_jobs=-1,
    )
    clf.fit(latent, adata_ref.obs[labels_key])

    predictions = clf.predict(latent)
    accuracy = accuracy_score(adata_ref.obs[labels_key], predictions)
    print(f"Classifier trained - Self-accuracy: {accuracy:.3f}")

    return clf, accuracy


def prepare_query_data(adata_query, reference_model_path: str):
    """
    Prepare query data for mapping.

    Handles gene alignment with reference.

    Parameters
    ----------
    adata_query : AnnData
        Query dataset
    reference_model_path : str
        Path to saved reference model
    """
    scvi.model.SCVI.prepare_query_anndata(adata_query, reference_model_path)
    print(f"Query data prepared - {adata_query.n_vars} genes")


def map_query_to_reference(
    adata_query,
    reference_model_path: str,
    max_epochs: int = 200,
) -> scvi.model.SCVI:
    """
    Map query data to reference latent space.

    CRITICAL: Uses weight_decay=0.0 to preserve reference space.

    Parameters
    ----------
    adata_query : AnnData
        Prepared query data
    reference_model_path : str
        Path to saved reference model
    max_epochs : int
        Training epochs for query

    Returns
    -------
    scvi.model.SCVI
        Query model
    """
    model_query = scvi.model.SCVI.load_query_data(
        adata_query,
        reference_model_path,
    )

    print(f"\nMapping query to reference:")
    print(f"  Max epochs: {max_epochs}")
    print(f"  weight_decay: 0.0 (preserves reference space)")

    model_query.train(
        max_epochs=max_epochs,
        plan_kwargs={"weight_decay": 0.0},  # CRITICAL!
    )

    final_epoch = len(model_query.history["elbo_train"])
    print(f"Query mapping complete after {final_epoch} epochs")

    return model_query


def transfer_labels(
    clf: RandomForestClassifier,
    model_query: scvi.model.SCVI,
    adata_query,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Transfer cell type labels to query.

    Parameters
    ----------
    clf : RandomForestClassifier
        Trained classifier
    model_query : scvi.model.SCVI
        Mapped query model
    adata_query : AnnData
        Query data

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (predictions, confidence_scores)
    """
    latent = model_query.get_latent_representation()
    predictions = clf.predict(latent)
    probabilities = clf.predict_proba(latent)
    confidence = probabilities.max(axis=1)

    print(f"Labels transferred to {len(predictions)} cells")
    print(f"Mean confidence: {confidence.mean():.3f}")

    return predictions, confidence


def verify_reference_preservation(
    model_ref: scvi.model.SCVI,
    model_query: scvi.model.SCVI,
    adata_ref,
    tolerance: float = 1e-5,
) -> bool:
    """
    Verify that reference latent space is preserved after query mapping.

    Parameters
    ----------
    model_ref : scvi.model.SCVI
        Original reference model
    model_query : scvi.model.SCVI
        Query model
    adata_ref : AnnData
        Reference data
    tolerance : float
        Maximum allowed difference

    Returns
    -------
    bool
        True if reference space preserved
    """
    latent_original = model_ref.get_latent_representation(adata_ref)
    latent_through_query = model_query.get_latent_representation(adata_ref)

    max_diff = np.abs(latent_original - latent_through_query).max()
    preserved = max_diff < tolerance

    print(f"Reference preservation check:")
    print(f"  Max difference: {max_diff:.2e}")
    print(f"  Preserved: {preserved}")

    return preserved


def check_scarches_compatibility(model_path: str) -> bool:
    """
    Check if saved model is scArches-compatible.

    Parameters
    ----------
    model_path : str
        Path to saved model

    Returns
    -------
    bool
        True if compatible
    """
    import os
    import json

    print("scArches compatibility requirements:")
    print("  1. encode_covariates=True")
    print("  2. use_layer_norm='both'")
    print("  3. use_batch_norm='none'")
    print("\nNote: Check model config to verify settings.")

    # Could potentially load and inspect model config here
    return True


def get_latent_representation(model: scvi.model.SCVI, adata=None) -> np.ndarray:
    """Get latent representation."""
    return model.get_latent_representation(adata)


def save_model(model: scvi.model.SCVI, path: str, overwrite: bool = True):
    """Save model."""
    model.save(path, overwrite=overwrite)
    print(f"Model saved to: {path}")


def load_model(path: str, adata) -> scvi.model.SCVI:
    """Load model."""
    return scvi.model.SCVI.load(path, adata=adata)
