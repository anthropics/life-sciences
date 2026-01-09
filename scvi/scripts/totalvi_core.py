#!/usr/bin/env python3
"""
TotalVI Core Functions

Modular functions for CITE-Seq analysis with TotalVI.
Use these for custom workflows or integration with other pipelines.
"""

import numpy as np
import scvi
from typing import Optional, Tuple, Union


def setup_mudata_totalvi(
    mdata,
    rna_layer: str = "counts",
    protein_layer: Optional[str] = None,
    batch_key: Optional[str] = None,
    rna_mod_key: str = "rna",
    protein_mod_key: str = "prot",
):
    """
    Register MuData for TotalVI.

    Parameters
    ----------
    mdata : mudata.MuData
        MuData object with RNA and protein modalities
    rna_layer : str
        Layer containing RNA counts
    protein_layer : str, optional
        Layer containing protein counts (None uses .X)
    batch_key : str, optional
        Column in obs for batch correction
    rna_mod_key : str
        Modality name for RNA
    protein_mod_key : str
        Modality name for protein
    """
    modalities = {
        "rna_layer": rna_mod_key,
        "protein_layer": protein_mod_key,
    }

    if batch_key:
        modalities["batch_key"] = rna_mod_key

    scvi.model.TOTALVI.setup_mudata(
        mdata,
        rna_layer=rna_layer,
        protein_layer=protein_layer,
        batch_key=batch_key,
        modalities=modalities,
    )

    print(f"MuData registered for TotalVI")
    print(f"  RNA modality: {rna_mod_key}")
    print(f"  Protein modality: {protein_mod_key}")
    if batch_key:
        print(f"  Batch key: {batch_key}")


def setup_anndata_totalvi(
    adata,
    layer: str = "counts",
    batch_key: Optional[str] = None,
    protein_obsm_key: str = "protein_expression",
):
    """
    Register AnnData for TotalVI.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData with protein in obsm
    layer : str
        Layer containing RNA counts
    batch_key : str, optional
        Column in obs for batch correction
    protein_obsm_key : str
        Key in obsm containing protein expression
    """
    scvi.model.TOTALVI.setup_anndata(
        adata,
        layer=layer,
        batch_key=batch_key,
        protein_expression_obsm_key=protein_obsm_key,
    )

    print(f"AnnData registered for TotalVI")
    print(f"  RNA layer: {layer}")
    print(f"  Protein obsm key: {protein_obsm_key}")
    if batch_key:
        print(f"  Batch key: {batch_key}")


def train_totalvi_model(
    data,
    n_latent: int = 20,
    n_layers_encoder: int = 2,
    n_layers_decoder: int = 1,
    max_epochs: int = 400,
    early_stopping: bool = True,
    batch_size: int = 128,
) -> scvi.model.TOTALVI:
    """
    Train TotalVI model.

    Parameters
    ----------
    data : MuData or AnnData
        Registered data object
    n_latent : int
        Latent space dimensions
    n_layers_encoder : int
        Number of encoder layers
    n_layers_decoder : int
        Number of decoder layers
    max_epochs : int
        Maximum training epochs
    early_stopping : bool
        Enable early stopping
    batch_size : int
        Training batch size

    Returns
    -------
    scvi.model.TOTALVI
        Trained model
    """
    model = scvi.model.TOTALVI(
        data,
        n_latent=n_latent,
        n_layers_encoder=n_layers_encoder,
        n_layers_decoder=n_layers_decoder,
    )

    print(f"\nModel architecture:")
    print(f"  Latent dimensions: {n_latent}")
    print(f"  Encoder layers: {n_layers_encoder}")
    print(f"  Decoder layers: {n_layers_decoder}")

    print(f"\nTraining with:")
    print(f"  Max epochs: {max_epochs}")
    print(f"  Early stopping: {early_stopping}")
    print(f"  Batch size: {batch_size}")

    model.train(
        max_epochs=max_epochs,
        early_stopping=early_stopping,
        batch_size=batch_size,
    )

    final_epoch = len(model.history["elbo_train"])
    print(f"\nTraining complete after {final_epoch} epochs")

    return model


def get_latent_representation(
    model: scvi.model.TOTALVI,
    adata=None,
) -> np.ndarray:
    """
    Extract joint latent representation.

    Parameters
    ----------
    model : scvi.model.TOTALVI
        Trained model
    adata : AnnData, optional
        Data to get representation for (default: training data)

    Returns
    -------
    np.ndarray
        Latent representation (n_cells x n_latent)
    """
    return model.get_latent_representation(adata)


def get_denoised_expression(
    model: scvi.model.TOTALVI,
    adata=None,
    n_samples: int = 25,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get denoised RNA and protein expression.

    Parameters
    ----------
    model : scvi.model.TOTALVI
        Trained model
    adata : AnnData, optional
        Data to denoise (default: training data)
    n_samples : int
        Number of Monte Carlo samples

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (denoised_rna, denoised_protein)
    """
    rna_denoised, protein_denoised = model.get_normalized_expression(
        adata=adata,
        n_samples=n_samples,
        return_mean=True,
    )

    return rna_denoised, protein_denoised


def get_foreground_probability(
    model: scvi.model.TOTALVI,
    adata=None,
    n_samples: int = 25,
) -> np.ndarray:
    """
    Get protein foreground probability.

    Probability that protein signal is real vs background noise.

    Parameters
    ----------
    model : scvi.model.TOTALVI
        Trained model
    adata : AnnData, optional
        Data to analyze (default: training data)
    n_samples : int
        Number of Monte Carlo samples

    Returns
    -------
    np.ndarray
        Foreground probability (n_cells x n_proteins)
        Values near 1 = real signal, near 0 = background
    """
    return model.get_protein_foreground_probability(
        adata=adata,
        n_samples=n_samples,
        return_mean=True,
    )


def run_differential_expression(
    model: scvi.model.TOTALVI,
    groupby: str,
    group1: str,
    group2: str,
    delta: float = 0.5,
    batch_correction: bool = True,
    fdr_target: float = 0.05,
) -> "pd.DataFrame":
    """
    Run differential expression analysis.

    Parameters
    ----------
    model : scvi.model.TOTALVI
        Trained model
    groupby : str
        Column for comparison (use 'modality:column' for MuData)
    group1 : str
        First group ID
    group2 : str
        Second group ID
    delta : float
        Effect size threshold
    batch_correction : bool
        Account for batch effects
    fdr_target : float
        FDR threshold

    Returns
    -------
    pd.DataFrame
        DE results with columns:
        - is_de_fdr: Significant after FDR correction
        - bayes_factor: Evidence strength
        - lfc_mean: Log fold change
        - raw_normalized_mean1/2: Expression in each group
    """
    de_results = model.differential_expression(
        groupby=groupby,
        group1=group1,
        group2=group2,
        delta=delta,
        batch_correction=batch_correction,
        fdr_target=fdr_target,
    )

    return de_results


def filter_de_results(
    de_results: "pd.DataFrame",
    protein_names: Optional[list] = None,
    min_bayes_factor_protein: float = 0.7,
    min_bayes_factor_rna: float = 3.0,
    min_lfc: float = 0.0,
    min_nonzero_proportion: float = 0.1,
) -> Tuple["pd.DataFrame", "pd.DataFrame"]:
    """
    Filter DE results into protein and RNA markers.

    Parameters
    ----------
    de_results : pd.DataFrame
        Results from run_differential_expression
    protein_names : list, optional
        List of protein feature names
    min_bayes_factor_protein : float
        Minimum Bayes factor for protein markers
    min_bayes_factor_rna : float
        Minimum Bayes factor for RNA markers
    min_lfc : float
        Minimum log fold change
    min_nonzero_proportion : float
        Minimum proportion of cells with expression

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        (protein_markers, rna_markers)
    """
    # Filter significant
    sig = de_results[de_results["is_de_fdr"]]

    if protein_names is not None:
        protein_markers = sig[
            (sig["bayes_factor"] > min_bayes_factor_protein)
            & (sig["lfc_mean"] > min_lfc)
            & (sig.index.isin(protein_names))
        ]

        rna_markers = sig[
            (sig["bayes_factor"] > min_bayes_factor_rna)
            & (sig["non_zeros_proportion1"] > min_nonzero_proportion)
            & (sig["lfc_mean"] > min_lfc)
            & (~sig.index.isin(protein_names))
        ]
    else:
        protein_markers = sig[sig["bayes_factor"] > min_bayes_factor_protein]
        rna_markers = sig[sig["bayes_factor"] > min_bayes_factor_rna]

    return protein_markers, rna_markers


def validate_citeseq_data(
    mdata,
    rna_mod_key: str = "rna",
    protein_mod_key: str = "prot",
) -> bool:
    """
    Validate CITE-Seq data for TotalVI.

    Parameters
    ----------
    mdata : mudata.MuData
        MuData object
    rna_mod_key : str
        RNA modality name
    protein_mod_key : str
        Protein modality name

    Returns
    -------
    bool
        True if data passes validation
    """
    issues = []
    recommendations = []

    rna = mdata.mod[rna_mod_key]
    protein = mdata.mod[protein_mod_key]

    # Check for dense matrices (TotalVI requirement)
    if hasattr(protein.X, "toarray"):
        issues.append("Protein matrix is sparse - will convert to dense")

    # Check cell counts
    if rna.n_obs < 1000:
        recommendations.append(f"Low cell count ({rna.n_obs}). Results may be noisy.")

    # Check protein count
    if protein.n_vars < 10:
        recommendations.append(
            f"Few proteins ({protein.n_vars}). Consider if TotalVI is needed."
        )
    if protein.n_vars > 200:
        recommendations.append(
            f"Many proteins ({protein.n_vars}). May need longer training."
        )

    # Check for zeros
    protein_X = protein.X.toarray() if hasattr(protein.X, "toarray") else protein.X
    zero_proteins = (protein_X.sum(axis=0) == 0).sum()
    if zero_proteins > 0:
        issues.append(f"{zero_proteins} proteins have zero counts in all cells")

    # Check batch info
    if "batch" in rna.obs.columns:
        batch_sizes = rna.obs["batch"].value_counts()
        if batch_sizes.min() < 100:
            recommendations.append(
                "Some batches have <100 cells. May affect batch correction."
            )

    # Check cell alignment
    if set(rna.obs_names) != set(protein.obs_names):
        issues.append("Cell barcodes don't match between RNA and protein")

    print("CITE-Seq Data Validation:")
    print("-" * 40)

    if issues:
        for issue in issues:
            print(f"  WARNING: {issue}")

    if recommendations:
        for rec in recommendations:
            print(f"  NOTE: {rec}")

    if not issues and not recommendations:
        print("  All checks passed!")

    return len([i for i in issues if "sparse" not in i]) == 0


def save_model(model: scvi.model.TOTALVI, path: str, overwrite: bool = True):
    """Save TotalVI model."""
    model.save(path, overwrite=overwrite)
    print(f"Model saved to: {path}")


def load_model(path: str, adata) -> scvi.model.TOTALVI:
    """Load TotalVI model."""
    return scvi.model.TOTALVI.load(path, adata=adata)
