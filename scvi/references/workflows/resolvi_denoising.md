# Spatial Transcriptomics Denoising with ResolVI

## Overview

ResolVI addresses noise and systematic biases in cellular-resolved spatial transcriptomics by correcting two primary issues:

1. **Erroneous co-expression patterns** from cellular segmentation errors (signal diffusion between adjacent cells)
2. **Unspecific background signal** contaminating expression measurements

The model decomposes observed expression into three components: true signal, diffusion artifacts, and background noise, enabling cleaner downstream analysis.

**Note**: Optimize cell segmentation before applying ResolVI - both processes are complementary.

---

## Workflow Steps

### Step 1: Environment Setup and Data Loading

```python
import os
import numpy as np
import scanpy as sc
import scvi
import matplotlib.pyplot as plt

# Set reproducibility
scvi.settings.seed = 0
sc.set_figure_params(figsize=(4, 4))

# Load spatial transcriptomics data
# Supports Xenium, MERFISH, CosMx, and other cellular-resolution platforms
adata = sc.read_h5ad("path/to/spatial_data.h5ad")

# Verify structure
print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")
print(f"Spatial coordinates available: {'X_spatial' in adata.obsm}")
```

**Data Requirements**:
- Raw counts in `adata.X` or `adata.layers['counts']`
- Spatial coordinates in `adata.obsm['X_spatial']` or `adata.obsm['spatial']`
- Cell type labels in `adata.obs` (optional, for semi-supervised mode)

---

### Step 2: Prepare Count Layer

```python
# Ensure counts layer exists
if 'counts' not in adata.layers:
    adata.layers['counts'] = adata.X.copy()

# Verify spatial coordinates key
spatial_key = 'X_spatial' if 'X_spatial' in adata.obsm else 'spatial'
print(f"Using spatial key: {spatial_key}")
```

---

### Step 3: Setup ResolVI AnnData

```python
# Setup for semi-supervised mode (with cell type labels)
scvi.external.RESOLVI.setup_anndata(
    adata,
    layer="counts",
    labels_key="cell_type"  # Column with cell type annotations
)

# Alternative: Unsupervised mode (no labels required)
# scvi.external.RESOLVI.setup_anndata(adata, layer="counts")
```

**Note**: The setup step automatically computes spatial neighbors using the coordinates.

---

### Step 4: Initialize and Train Model

```python
# Create model - semi-supervised mode
model = scvi.external.RESOLVI(
    adata,
    semisupervised=True  # Set False for unsupervised
)

# Train the model
model.train(
    max_epochs=100,  # Use 50 for quick testing, 100+ for production
    early_stopping=True
)

# Plot training history
plt.figure(figsize=(8, 4))
plt.plot(model.history['elbo_train'].values)
plt.xlabel('Epoch')
plt.ylabel('ELBO')
plt.title('ResolVI Training')
plt.show()
```

---

### Step 5: Extract Corrected Expression

```python
# Get denoised/corrected expression using posterior sampling
samples_corr = model.sample_posterior(
    model=model.module.model_corrected,
    return_sites=["px_rate"],
    summary_fun={"post_sample_q50": np.median},
    num_samples=3,
    summary_frequency=30
)

# Store corrected expression
adata.layers["corrected_expression"] = samples_corr.loc["post_sample_q50", "px_rate"]

print("Corrected expression stored in adata.layers['corrected_expression']")
```

---

### Step 6: Cell Type Prediction (Semi-supervised)

```python
# Get cell type predictions (soft probabilities)
adata.obsm["resolvi_celltypes"] = model.predict(
    adata,
    num_samples=3,
    soft=True  # Returns probability matrix
)

# Get hard predictions
adata.obs["resolvi_predicted"] = adata.obsm["resolvi_celltypes"].idxmax(axis=1)

# Compare with original labels
print(f"Original vs Predicted agreement: {(adata.obs['cell_type'] == adata.obs['resolvi_predicted']).mean():.2%}")
```

---

### Step 7: Extract Latent Representation

```python
# Get latent representation for downstream analysis
adata.obsm["X_resolVI"] = model.get_latent_representation(adata)

# Compute neighbors and UMAP
sc.pp.neighbors(adata, use_rep="X_resolVI")
sc.tl.umap(adata)

# Visualize
sc.pl.umap(adata, color=["cell_type", "resolvi_predicted"], ncols=2)
```

---

### Step 8: Analyze Noise Components

```python
# Decompose signal into true expression, diffusion, and background
samples = model.sample_posterior(
    model=model.module.model_residuals,
    return_sites=["mixture_proportions"],
    summary_fun={"post_sample_means": np.mean},
    num_samples=3,
    summary_frequency=100
)

# Extract proportions
proportions = samples.loc["post_sample_means", "mixture_proportions"]
adata.obs["true_proportion"] = proportions[:, 0]
adata.obs["diffusion_proportion"] = proportions[:, 1]
adata.obs["background_proportion"] = proportions[:, 2]

# Visualize noise components spatially
sc.pl.spatial(adata, color=["true_proportion", "diffusion_proportion", "background_proportion"],
              spot_size=30, ncols=3)
```

---

### Step 9: Transfer to New Data (Query Mapping)

```python
# Prepare query data (new spatial sample)
query_adata = sc.read_h5ad("path/to/new_sample.h5ad")
query_adata.obs["cell_type"] = "unknown"  # Mark as unlabeled

# Prepare query for transfer learning
model.prepare_query_anndata(query_adata, reference_model=model)

# Load query into model
query_model = model.load_query_data(query_adata, reference_model=model)

# Fine-tune on query
query_model.train(max_epochs=20)

# Predict cell types on query
query_adata.obs["predicted_celltype"] = query_model.predict(query_adata, soft=False)
```

---

### Step 10: Save Results

```python
# Save model
model.save("resolvi_model", overwrite=True)

# Save annotated data
adata.write_h5ad("spatial_denoised.h5ad")

# Reload model later
# model = scvi.external.RESOLVI.load("resolvi_model", adata=adata)
```

---

## Key Parameters Reference

| Parameter | Default | Description |
|-----------|---------|-------------|
| `semisupervised` | True | Use cell type labels for guided learning |
| `max_epochs` | 100 | Training epochs (50 for testing, 100+ for production) |
| `num_samples` | 3 | Posterior samples for predictions |
| `soft` | True | Return probability distributions vs hard labels |
| `summary_fun` | median | Aggregation for posterior (median recommended) |

---

## Advanced Visualization

### Before/After Denoising Comparison

```python
import matplotlib.pyplot as plt

# Select marker genes to visualize
marker_genes = ["CD3D", "MS4A1", "CD14"]  # Replace with your markers

fig, axes = plt.subplots(len(marker_genes), 2, figsize=(10, 4*len(marker_genes)))

for idx, gene in enumerate(marker_genes):
    gene_idx = adata.var_names.get_loc(gene)

    # Original expression
    adata.obs["_orig"] = np.asarray(adata.layers["counts"][:, gene_idx]).flatten()
    sc.pl.spatial(adata, color="_orig", ax=axes[idx, 0], show=False,
                  title=f"{gene} (Original)", cmap="viridis", spot_size=30)

    # Corrected expression
    adata.obs["_corr"] = adata.layers["corrected_expression"][:, gene_idx]
    sc.pl.spatial(adata, color="_corr", ax=axes[idx, 1], show=False,
                  title=f"{gene} (Corrected)", cmap="viridis", spot_size=30)

plt.tight_layout()
plt.savefig("before_after_denoising.png", dpi=150)
```

### Signal Component Breakdown

```python
# Visualize true signal vs noise proportions
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

components = ["true_proportion", "diffusion_proportion", "background_proportion"]
titles = ["True Signal", "Diffusion Artifact", "Background Noise"]
cmaps = ["Greens", "Oranges", "Reds"]

for ax, comp, title, cmap in zip(axes, components, titles, cmaps):
    sc.pl.spatial(adata, color=comp, ax=ax, show=False, title=title,
                  cmap=cmap, spot_size=30, vmax="p99")

plt.suptitle("Signal Component Decomposition", fontsize=14)
plt.tight_layout()
plt.savefig("signal_decomposition.png", dpi=150)
```

### Cell Type Prediction Confidence

```python
# Visualize prediction confidence spatially
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Predicted cell types
sc.pl.spatial(adata, color="resolvi_predicted", ax=axes[0], show=False,
              title="Predicted Cell Types", spot_size=30, palette="tab20")

# Prediction confidence (max probability)
adata.obs["pred_confidence"] = adata.obsm["resolvi_celltypes"].max(axis=1)
sc.pl.spatial(adata, color="pred_confidence", ax=axes[1], show=False,
              title="Prediction Confidence", cmap="plasma", spot_size=30)

plt.tight_layout()
plt.savefig("predictions_confidence.png", dpi=150)
```

### Uncertainty Map

```python
# Calculate prediction entropy as uncertainty measure
from scipy.stats import entropy

pred_probs = adata.obsm["resolvi_celltypes"].values
adata.obs["prediction_entropy"] = [entropy(p) for p in pred_probs]

fig, ax = plt.subplots(figsize=(8, 8))
sc.pl.spatial(adata, color="prediction_entropy", ax=ax, show=False,
              title="Prediction Uncertainty (Entropy)", cmap="YlOrRd", spot_size=30)
plt.savefig("prediction_uncertainty.png", dpi=150)
```

---

## Troubleshooting

### Common Issues

1. **No spatial neighbors computed**: Ensure coordinates are in `adata.obsm['X_spatial']` or `adata.obsm['spatial']`

2. **Poor denoising quality**:
   - Increase `max_epochs` to 150+
   - Ensure sufficient cells per type (>100 recommended)
   - Check that raw counts (not normalized) are used

3. **Memory errors**:
   - Process tissue regions separately
   - Reduce `num_samples` for posterior sampling

4. **Label mismatch in transfer**:
   - Ensure query data uses "unknown" for unlabeled cells
   - Check gene overlap between reference and query

---

## References

- [ResolVI Tutorial](https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/spatial/resolVI_tutorial.html)
- [scvi-tools Documentation](https://docs.scvi-tools.org/)
