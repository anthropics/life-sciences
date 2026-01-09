# Cellular Microenvironment Modeling with scVIVA

## Overview

scVIVA (single-cell Variational Inference with Vicinity Attention) is a deep generative model that leverages both cell-intrinsic and neighboring gene expression profiles to produce:

- **Stochastic embeddings** reflecting internal cell state AND surrounding tissue context
- **Fine-grained cell partitions** based on both intrinsic and environmental factors
- **Niche-aware differential expression** testing how environment influences gene expression

This enables discovery of cell subtypes that emerge only in specific tissue contexts (e.g., tumor vs stromal endothelial cells).

---

## Workflow Steps

### Step 1: Environment Setup and Data Loading

```python
import os
import random
import numpy as np
import scanpy as sc
import scvi
import torch
import matplotlib.pyplot as plt
import seaborn as sns

# Set reproducibility
scvi.settings.seed = 0
random.seed(42)
np.random.seed(42)
torch.manual_seed(42)
sc.set_figure_params(figsize=(4, 4))
sns.set_theme()

# Load spatial transcriptomics data
adata = sc.read_h5ad("path/to/spatial_data.h5ad")

print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")
print(f"Cell types: {adata.obs['cell_type'].nunique()}")
```

**Data Requirements**:
- Raw counts in `adata.X` or `adata.layers['counts']`
- Spatial coordinates in `adata.obsm['spatial']`
- Cell type annotations in `adata.obs` (recommended)
- Sample/batch information if multi-sample

---

### Step 2: Prepare Data

```python
# Ensure counts layer
if 'counts' not in adata.layers:
    adata.layers['counts'] = adata.X.copy()

# Optional: Run scANVI first to get expression embeddings
# This provides better initialization for scVIVA
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="sample")
scvi_model = scvi.model.SCVI(adata)
scvi_model.train(max_epochs=100)
adata.obsm["X_scANVI"] = scvi_model.get_latent_representation()
```

---

### Step 3: Preprocessing for scVIVA

```python
# Define setup parameters
setup_kwargs = {
    "sample_key": "sample",  # Batch/sample column
    "labels_key": "cell_type",  # Cell type column
    "cell_coordinates_key": "spatial",  # Spatial coordinates key
    "expression_embedding_key": "X_scANVI",  # Optional: pre-computed embedding
}

# Compute spatial neighborhood graph
scvi.external.SCVIVA.preprocessing_anndata(
    adata,
    k_nn=20,  # Number of spatial neighbors
    **setup_kwargs
)
```

**Key Parameter**: `k_nn` defines the neighborhood size - larger values capture broader tissue context.

---

### Step 4: Setup and Train Model

```python
# Register AnnData
scvi.external.SCVIVA.setup_anndata(
    adata,
    layer="counts",
    batch_key="sample",
    **setup_kwargs
)

# Initialize model
model = scvi.external.SCVIVA(adata)

# Train
model.train(
    max_epochs=600,
    early_stopping=True,
    check_val_every_n_epoch=1,
    batch_size=512,
    plan_kwargs={"lr": 5e-4}
)

# Plot training history
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
model.history["elbo_validation"].plot(ax=axes[0], title="ELBO")
model.history["niche_compo_validation"].plot(ax=axes[1], title="Niche Composition")
plt.tight_layout()
plt.show()
```

---

### Step 5: Extract Latent Representation

```python
# Get scVIVA embedding (combines intrinsic + environment)
adata.obsm["X_scVIVA"] = model.get_latent_representation()

# Compute neighbors and UMAP
sc.pp.neighbors(adata, use_rep="X_scVIVA", n_neighbors=30)
sc.tl.umap(adata)

# Visualize
sc.pl.umap(adata, color=["cell_type", "sample"], ncols=2)
```

---

### Step 6: Identify Tissue Niches via Clustering

```python
# Cluster to identify niches (cell types in specific contexts)
sc.tl.leiden(adata, resolution=0.5, key_added="scviva_clusters")

# Visualize clusters spatially and in UMAP
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
sc.pl.umap(adata, color="scviva_clusters", ax=axes[0], show=False, title="UMAP")
sc.pl.spatial(adata, color="scviva_clusters", spot_size=30, ax=axes[1], show=False, title="Spatial")
plt.tight_layout()
plt.show()
```

---

### Step 7: Subset Analysis (e.g., Endothelial Cells)

```python
# Focus on specific cell type to find context-dependent subtypes
cell_type_of_interest = "Endothelial"
adata_subset = adata[adata.obs["cell_type"] == cell_type_of_interest].copy()

# Re-cluster within cell type
sc.pp.neighbors(adata_subset, use_rep="X_scVIVA", n_neighbors=15)
sc.tl.leiden(adata_subset, resolution=0.3, key_added="subtype_clusters")

sc.pl.umap(adata_subset, color="subtype_clusters", title=f"{cell_type_of_interest} Subtypes")
```

---

### Step 8: Niche-Aware Differential Expression

```python
# Compare cell groups AND their neighborhoods
# Example: tumor endothelial (cluster 1) vs stromal endothelial (cluster 0)

DE_results = model.differential_expression(
    adata_subset,
    groupby="subtype_clusters",
    group1="1",  # Tumor-associated
    group2="0",  # Stromal-associated
    k_nn=6,  # Neighbors for niche analysis
    delta=[0.05, 0.15, 0.05, 0.05],  # LFC thresholds for 4 comparisons
    niche_mode=True,  # Enable niche comparisons
    n_samples_overall=1e5,
    fdr_target=0.2,
    pseudocounts=1e-4
)

print(f"Significant DE genes: {len(DE_results)}")
print(DE_results.head(20))
```

**The 4 Comparisons**:
1. G1 vs G2: Cell group comparison
2. G1 vs N1: Cells vs their own neighbors
3. N1 vs G2: Group 1 neighbors vs Group 2 cells
4. N1 vs N2: Neighbor comparison

---

### Step 9: Visualize Niche Marker Genes

```python
# Get normalized expression for visualization
sc.pp.normalize_total(adata_subset, target_sum=1e4)
sc.pp.log1p(adata_subset)

# Plot top DE genes spatially
top_genes = DE_results.head(4).index.tolist()

sc.pl.spatial(
    adata_subset,
    color=top_genes,
    spot_size=30,
    ncols=2,
    cmap="plasma",
    vmax="p99"
)
```

---

### Step 10: Save Results

```python
# Save model
model.save("scviva_model", overwrite=True)

# Save annotated data
adata.write_h5ad("spatial_niche_analyzed.h5ad")

# Save DE results
DE_results.to_csv("niche_de_results.csv")

# Reload later
# model = scvi.external.SCVIVA.load("scviva_model", adata=adata)
```

---

## Key Parameters Reference

| Parameter | Default | Description |
|-----------|---------|-------------|
| `k_nn` (preprocessing) | 20 | Spatial neighborhood size |
| `max_epochs` | 600 | Training epochs |
| `batch_size` | 512 | Mini-batch size |
| `lr` | 5e-4 | Learning rate |
| `delta` (DE) | [0.05, 0.15, 0.05, 0.05] | LFC thresholds per comparison |
| `niche_mode` | True | Enable neighborhood comparisons |
| `fdr_target` | 0.2 | FDR correction threshold |

---

## Troubleshooting

### Common Issues

1. **Poor niche separation**:
   - Increase `k_nn` to capture broader context
   - Ensure sufficient cells per niche (>200)
   - Train longer (increase `max_epochs`)

2. **Memory errors**:
   - Reduce `batch_size` to 256
   - Subset to region of interest

3. **No DE genes found**:
   - Adjust `delta` thresholds (lower = more sensitive)
   - Increase `n_samples_overall`
   - Check cluster sizes (need >50 cells per group)

4. **Slow training**:
   - Enable GPU: `scvi.settings.device = "cuda"`
   - Reduce `k_nn` for faster preprocessing

---

## References

- [scVIVA Tutorial](https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/spatial/scVIVA_tutorial.html)
- [scvi-tools Documentation](https://docs.scvi-tools.org/)
