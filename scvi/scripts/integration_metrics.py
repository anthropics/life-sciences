"""
Integration quality metrics for single-cell data integration.

This module provides functions to calculate benchmarking metrics for
evaluating integration quality, following scib-metrics methodology.

Metrics include:
- Batch mixing metrics (how well batches are integrated)
- Bio-conservation metrics (how well biological signal is preserved)
"""

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import silhouette_score, normalized_mutual_info_score, adjusted_rand_score
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import issparse
from typing import Optional, List, Dict


def calculate_integration_metrics(
    adata,
    batch_key: str,
    label_key: Optional[str] = None,
    embed_key: str = 'X_scVI',
    cluster_key: Optional[str] = None
) -> Dict[str, float]:
    """
    Calculate comprehensive integration quality metrics.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with integration results.
    batch_key : str
        Key in adata.obs for batch information.
    label_key : str, optional
        Key in adata.obs for cell type labels.
    embed_key : str, default='X_scVI'
        Key in adata.obsm for integrated embedding.
    cluster_key : str, optional
        Key in adata.obs for cluster assignments.

    Returns
    -------
    dict
        Dictionary of metric names and values.
    """
    metrics = {}
    embedding = adata.obsm[embed_key]

    # Batch mixing metrics
    print(f"    Computing batch mixing metrics...")

    # Silhouette score (batch) - lower is better for integration
    batch_labels = adata.obs[batch_key].values
    try:
        batch_silhouette = silhouette_score(embedding, batch_labels)
        # Invert so higher is better (well-mixed batches)
        metrics['batch_silhouette'] = -batch_silhouette
    except Exception as e:
        print(f"      Warning: Could not compute batch silhouette: {e}")
        metrics['batch_silhouette'] = np.nan

    # iLISI (integration LISI) - higher is better
    try:
        ilisi = compute_lisi(embedding, adata.obs, batch_key)
        metrics['iLISI'] = ilisi.mean()
    except Exception as e:
        print(f"      Warning: Could not compute iLISI: {e}")
        metrics['iLISI'] = np.nan

    # Graph connectivity (batch)
    try:
        graph_conn = graph_connectivity(adata, batch_key)
        metrics['graph_connectivity'] = graph_conn
    except Exception as e:
        print(f"      Warning: Could not compute graph connectivity: {e}")
        metrics['graph_connectivity'] = np.nan

    # Bio-conservation metrics (if labels available)
    if label_key is not None and label_key in adata.obs.columns:
        print(f"    Computing bio-conservation metrics...")
        cell_labels = adata.obs[label_key].values

        # Silhouette score (cell type) - higher is better
        try:
            # Filter out any unlabeled cells
            valid_mask = ~pd.isna(cell_labels)
            if valid_mask.sum() > 100:
                label_silhouette = silhouette_score(
                    embedding[valid_mask],
                    cell_labels[valid_mask]
                )
                metrics['label_silhouette'] = label_silhouette
            else:
                metrics['label_silhouette'] = np.nan
        except Exception as e:
            print(f"      Warning: Could not compute label silhouette: {e}")
            metrics['label_silhouette'] = np.nan

        # cLISI (cell type LISI) - lower is better (cell types should be separated)
        try:
            clisi = compute_lisi(embedding, adata.obs, label_key)
            # Invert so higher is better
            metrics['cLISI'] = 1.0 / clisi.mean() if clisi.mean() > 0 else np.nan
        except Exception as e:
            print(f"      Warning: Could not compute cLISI: {e}")
            metrics['cLISI'] = np.nan

    # Clustering metrics (if cluster key provided)
    if cluster_key is not None and cluster_key in adata.obs.columns:
        print(f"    Computing clustering metrics...")
        cluster_labels = adata.obs[cluster_key].values

        if label_key is not None and label_key in adata.obs.columns:
            cell_labels = adata.obs[label_key].values

            # NMI (Normalized Mutual Information)
            try:
                valid_mask = ~pd.isna(cell_labels)
                nmi = normalized_mutual_info_score(
                    cell_labels[valid_mask],
                    cluster_labels[valid_mask]
                )
                metrics['NMI'] = nmi
            except Exception as e:
                print(f"      Warning: Could not compute NMI: {e}")
                metrics['NMI'] = np.nan

            # ARI (Adjusted Rand Index)
            try:
                valid_mask = ~pd.isna(cell_labels)
                ari = adjusted_rand_score(
                    cell_labels[valid_mask],
                    cluster_labels[valid_mask]
                )
                metrics['ARI'] = ari
            except Exception as e:
                print(f"      Warning: Could not compute ARI: {e}")
                metrics['ARI'] = np.nan

    # Compute overall score (weighted average)
    batch_metrics = ['batch_silhouette', 'iLISI', 'graph_connectivity']
    bio_metrics = ['label_silhouette', 'cLISI', 'NMI', 'ARI']

    batch_scores = [metrics.get(m, np.nan) for m in batch_metrics if m in metrics]
    bio_scores = [metrics.get(m, np.nan) for m in bio_metrics if m in metrics]

    batch_scores = [s for s in batch_scores if not np.isnan(s)]
    bio_scores = [s for s in bio_scores if not np.isnan(s)]

    if batch_scores:
        metrics['batch_score'] = np.mean(batch_scores)
    if bio_scores:
        metrics['bio_score'] = np.mean(bio_scores)
    if batch_scores and bio_scores:
        metrics['overall_score'] = 0.4 * metrics['batch_score'] + 0.6 * metrics['bio_score']

    return metrics


def compute_lisi(
    embedding: np.ndarray,
    metadata: pd.DataFrame,
    label_key: str,
    perplexity: int = 30
) -> np.ndarray:
    """
    Compute Local Inverse Simpson's Index (LISI).

    LISI measures the effective number of categories in local neighborhoods.
    - For batch: higher LISI = better mixing
    - For cell type: lower LISI = better preservation

    Parameters
    ----------
    embedding : np.ndarray
        Low-dimensional embedding (n_cells × n_dims).
    metadata : pd.DataFrame
        Cell metadata containing label column.
    label_key : str
        Column name for labels.
    perplexity : int, default=30
        Perplexity for computing local neighborhoods.

    Returns
    -------
    np.ndarray
        LISI scores for each cell.
    """
    n_cells = embedding.shape[0]
    n_neighbors = min(3 * perplexity, n_cells - 1)

    # Get k-nearest neighbors
    nn = NearestNeighbors(n_neighbors=n_neighbors, algorithm='ball_tree')
    nn.fit(embedding)
    distances, indices = nn.kneighbors(embedding)

    # Get labels
    labels = metadata[label_key].values
    unique_labels = np.unique(labels)
    label_to_idx = {l: i for i, l in enumerate(unique_labels)}

    # Compute LISI for each cell
    lisi_scores = np.zeros(n_cells)

    for i in range(n_cells):
        # Get neighbor labels
        neighbor_labels = labels[indices[i]]

        # Count label frequencies
        label_counts = np.zeros(len(unique_labels))
        for label in neighbor_labels:
            label_counts[label_to_idx[label]] += 1

        # Compute Simpson's Index and invert
        proportions = label_counts / label_counts.sum()
        simpson = np.sum(proportions ** 2)
        lisi_scores[i] = 1.0 / simpson if simpson > 0 else 1.0

    return lisi_scores


def graph_connectivity(adata, batch_key: str) -> float:
    """
    Compute graph connectivity metric for batch integration.

    Measures the fraction of cells whose k-nearest neighbors
    include cells from other batches.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with computed neighbors.
    batch_key : str
        Key in adata.obs for batch information.

    Returns
    -------
    float
        Graph connectivity score (0-1, higher is better).
    """
    # Use the connectivities matrix from neighbors
    if 'connectivities' not in adata.obsp:
        return np.nan

    connectivities = adata.obsp['connectivities']
    if issparse(connectivities):
        connectivities = connectivities.toarray()

    batches = adata.obs[batch_key].values
    n_cells = len(batches)

    # For each cell, check if it has neighbors from different batches
    cross_batch_connections = 0

    for i in range(n_cells):
        cell_batch = batches[i]
        neighbors = np.where(connectivities[i] > 0)[0]

        if len(neighbors) > 0:
            neighbor_batches = batches[neighbors]
            if any(neighbor_batches != cell_batch):
                cross_batch_connections += 1

    return cross_batch_connections / n_cells


def compare_embeddings(
    adata,
    embed_keys: List[str],
    batch_key: str,
    label_key: Optional[str] = None,
    cluster_keys: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Compare multiple integration embeddings using benchmarking metrics.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with multiple embeddings.
    embed_keys : list of str
        Keys in adata.obsm for embeddings to compare.
    batch_key : str
        Key in adata.obs for batch information.
    label_key : str, optional
        Key in adata.obs for cell type labels.
    cluster_keys : list of str, optional
        Keys in adata.obs for cluster assignments (one per embedding).

    Returns
    -------
    pd.DataFrame
        DataFrame with metrics for each embedding.
    """
    results = []

    if cluster_keys is None:
        cluster_keys = [None] * len(embed_keys)

    for embed_key, cluster_key in zip(embed_keys, cluster_keys):
        if embed_key not in adata.obsm:
            print(f"  Warning: {embed_key} not found in adata.obsm, skipping")
            continue

        print(f"\n  Evaluating {embed_key}...")
        metrics = calculate_integration_metrics(
            adata,
            batch_key=batch_key,
            label_key=label_key,
            embed_key=embed_key,
            cluster_key=cluster_key
        )
        metrics['embedding'] = embed_key
        results.append(metrics)

    df = pd.DataFrame(results)
    if 'embedding' in df.columns:
        df = df.set_index('embedding')

    return df


def format_metrics_table(metrics_df: pd.DataFrame) -> str:
    """
    Format metrics DataFrame as a readable table string.

    Parameters
    ----------
    metrics_df : pd.DataFrame
        DataFrame with integration metrics.

    Returns
    -------
    str
        Formatted table string.
    """
    # Select key columns in preferred order
    preferred_order = [
        'batch_silhouette', 'iLISI', 'graph_connectivity', 'batch_score',
        'label_silhouette', 'cLISI', 'NMI', 'ARI', 'bio_score',
        'overall_score'
    ]

    cols = [c for c in preferred_order if c in metrics_df.columns]
    df_display = metrics_df[cols].copy()

    # Round for display
    df_display = df_display.round(4)

    lines = ["\n  Integration Metrics Summary"]
    lines.append("  " + "=" * 60)

    # Header
    header = "  {:12s}".format("Method")
    for col in cols:
        header += " {:>12s}".format(col[:12])
    lines.append(header)
    lines.append("  " + "-" * 60)

    # Data rows
    for idx, row in df_display.iterrows():
        row_str = "  {:12s}".format(str(idx)[:12])
        for col in cols:
            val = row[col]
            if pd.isna(val):
                row_str += " {:>12s}".format("N/A")
            else:
                row_str += " {:>12.4f}".format(val)
        lines.append(row_str)

    lines.append("  " + "=" * 60)

    return "\n".join(lines)
