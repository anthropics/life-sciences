# Spatial Deconvolution with Stereoscope

## Overview

Stereoscope is a probabilistic method for deconvoluting cell type compositions in spatial transcriptomics. It uses a two-stage approach:

1. **RNAStereoscope**: Learns cell-type-specific gene expression profiles from scRNA-seq reference
2. **SpatialStereoscope**: Applies learned profiles to estimate cell type proportions in spatial spots

**Key Features**:
- Robust to sequencing depth differences
- Works with any tissue type
- Compatible with 10x Visium and similar platforms

---

## Workflow Steps

### Step 1: Environment Setup

```python
import os
import numpy as np
import scanpy as sc
import scvi
import matplotlib.pyplot as plt
import seaborn as sns
import torch
from scvi.external import RNAStereoscope, SpatialStereoscope

scvi.settings.seed = 0
sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")
```

---

### Step 2: Load scRNA-seq Reference Data

```python
# Load annotated scRNA-seq reference
sc_adata = sc.read_h5ad("path/to/scrna_reference.h5ad")

print(f"Reference: {sc_adata.n_obs} cells, {sc_adata.n_vars} genes")
print(f"Cell types:\n{sc_adata.obs['cell_type'].value_counts()}")

# Ensure raw counts are available
if 'counts' not in sc_adata.layers:
    sc_adata.layers['counts'] = sc_adata.X.copy()
```

---

### Step 3: Preprocess Reference Data

```python
# Filter low-abundance genes
sc.pp.filter_genes(sc_adata, min_counts=10)

# Remove mitochondrial genes (recommended for deconvolution)
non_mito_genes = [g for g in sc_adata.var_names if not g.startswith('MT-')]
sc_adata = sc_adata[:, non_mito_genes].copy()

# Normalize and log-transform for HVG selection (keep raw counts)
sc.pp.normalize_total(sc_adata, target_sum=1e5)
sc.pp.log1p(sc_adata)

# Select highly variable genes
sc.pp.highly_variable_genes(
    sc_adata,
    n_top_genes=7000,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="batch" if "batch" in sc_adata.obs.columns else None
)

print(f"After preprocessing: {sc_adata.n_vars} genes")
```

---

### Step 4: Load Spatial Transcriptomics Data

```python
# Load Visium or other spatial data
st_adata = sc.read_h5ad("path/to/spatial_data.h5ad")

# Or from 10x Space Ranger output
# st_adata = sc.read_visium("path/to/spaceranger/outs")
# st_adata.var_names_make_unique()

print(f"Spatial: {st_adata.n_obs} spots, {st_adata.n_vars} genes")

# Basic QC
st_adata.var['mt'] = st_adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(st_adata, qc_vars=['mt'], inplace=True)
sc.pp.filter_cells(st_adata, min_counts=500)
sc.pp.filter_cells(st_adata, min_genes=500)

# Ensure counts layer
if 'counts' not in st_adata.layers:
    st_adata.layers['counts'] = st_adata.X.copy()
```

---

### Step 5: Align Gene Sets

```python
# Find shared genes between reference and spatial data
shared_genes = np.intersect1d(sc_adata.var_names, st_adata.var_names)
print(f"Shared genes: {len(shared_genes)}")

# Subset both datasets to shared genes
sc_adata = sc_adata[:, shared_genes].copy()
st_adata = st_adata[:, shared_genes].copy()
```

---

### Step 6: Train Reference Model (RNAStereoscope)

```python
# Setup reference data for Stereoscope
RNAStereoscope.setup_anndata(
    sc_adata,
    layer="counts",
    labels_key="cell_type"  # Cell type annotation column
)

# Initialize reference model
sc_model = RNAStereoscope(sc_adata)

# Train reference model
sc_model.train(max_epochs=100)

# Plot training
plt.figure(figsize=(8, 4))
plt.plot(sc_model.history['elbo_train'].values)
plt.xlabel('Epoch')
plt.ylabel('ELBO')
plt.title('Reference Model Training')
plt.show()

# Save reference model
sc_model.save("stereoscope_reference_model", overwrite=True)
```

---

### Step 7: Train Spatial Model (SpatialStereoscope)

```python
# Setup spatial data
SpatialStereoscope.setup_anndata(st_adata, layer="counts")

# Initialize spatial model from reference
spatial_model = SpatialStereoscope.from_rna_model(st_adata, sc_model)

# Train spatial deconvolution model
# Note: Requires more epochs than reference model
spatial_model.train(max_epochs=2000)

# Plot training
plt.figure(figsize=(8, 4))
plt.plot(spatial_model.history['elbo_train'].values)
plt.xlabel('Epoch')
plt.ylabel('ELBO')
plt.title('Spatial Model Training')
plt.show()

# Save spatial model
spatial_model.save("stereoscope_spatial_model", overwrite=True)
```

---

### Step 8: Extract Cell Type Proportions

```python
# Get deconvolution results
st_adata.obsm["deconvolution"] = spatial_model.get_proportions()

# Add proportions to obs for visualization
for ct in st_adata.obsm["deconvolution"].columns:
    st_adata.obs[ct] = st_adata.obsm["deconvolution"][ct]

# View proportion summary
print("Cell type proportion statistics:")
print(st_adata.obsm["deconvolution"].describe())
```

---

### Step 9: Visualize Results

```python
# Configure visualization
sc.settings.set_figure_params(
    dpi=100,
    color_map="inferno",
    dpi_save=200,
    vector_friendly=True
)

# Get cell types for plotting
cell_types = st_adata.obsm["deconvolution"].columns.tolist()

# Spatial visualization of cell type proportions
sc.pl.spatial(
    st_adata,
    img_key="hires",
    color=cell_types[:6],  # First 6 cell types
    ncols=3,
    size=1.2,
    cmap="inferno",
    vmax=1.0
)

# Heatmap of proportions across spots
plt.figure(figsize=(12, 8))
sns.clustermap(
    st_adata.obsm["deconvolution"],
    cmap="viridis",
    figsize=(12, 8),
    col_cluster=True,
    row_cluster=True
)
plt.title("Cell Type Proportions Across Spots")
plt.show()
```

---

### Step 10: Save Results

```python
# Save deconvolved spatial data
st_adata.write_h5ad("spatial_deconvolved_stereoscope.h5ad")

# Export proportions as CSV
st_adata.obsm["deconvolution"].to_csv("stereoscope_proportions.csv")

# Reload models later
# sc_model = RNAStereoscope.load("stereoscope_reference_model", adata=sc_adata)
# spatial_model = SpatialStereoscope.load("stereoscope_spatial_model", adata=st_adata)
```

---

## Key Parameters Reference

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_epochs` (Reference) | 100 | Training epochs for scRNA-seq model |
| `max_epochs` (Spatial) | 2000 | Training epochs for spatial model |
| `n_top_genes` | 7000 | HVGs for reference preprocessing |
| `labels_key` | - | Column with cell type annotations |

---

## Best Practices

### Reference Data Quality
- Include all expected cell types in reference
- Balance cell type representation when possible
- Remove low-quality cells and doublets
- Use batch correction for multi-sample references

### Gene Selection
- Remove mitochondrial genes before deconvolution
- Use 5000-7000 HVGs for best results
- Ensure marker genes are included in shared gene set

### Training Tips
- Reference model: 100-150 epochs usually sufficient
- Spatial model: 2000+ epochs for convergence
- Monitor ELBO loss for convergence

---

## Advanced Visualization

### Publication-Quality Spatial Maps

```python
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Multi-panel figure with all cell types
cell_types = st_adata.obsm["deconvolution"].columns.tolist()
n_types = len(cell_types)
ncols = 4
nrows = int(np.ceil(n_types / ncols))

fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows))
axes = axes.flatten()

for idx, ct in enumerate(cell_types):
    st_adata.obs[f"prop_{ct}"] = st_adata.obsm["deconvolution"][ct]
    sc.pl.spatial(st_adata, color=f"prop_{ct}", ax=axes[idx], show=False,
                  title=ct, cmap="viridis", vmax="p95", size=1.2)

for j in range(n_types, len(axes)):
    axes[j].axis("off")

plt.tight_layout()
plt.savefig("proportions_grid.png", dpi=200, bbox_inches="tight")
```

### Dominant Cell Type Map

```python
# Show which cell type dominates each spot
props = st_adata.obsm["deconvolution"]
st_adata.obs["dominant_celltype"] = props.idxmax(axis=1)
st_adata.obs["max_proportion"] = props.max(axis=1)

# Filter low-confidence spots
st_adata.obs.loc[st_adata.obs["max_proportion"] < 0.1, "dominant_celltype"] = "Mixed"

sc.pl.spatial(st_adata, color="dominant_celltype", title="Dominant Cell Type per Spot",
              size=1.3, palette="tab20")
```

### Cell Type Co-localization Heatmap

```python
import seaborn as sns

# Calculate spatial correlation between cell types
props = st_adata.obsm["deconvolution"]
corr_matrix = props.corr(method="spearman")

plt.figure(figsize=(10, 8))
sns.clustermap(corr_matrix, cmap="RdBu_r", center=0, vmin=-1, vmax=1,
               annot=True, fmt=".2f", figsize=(10, 8))
plt.title("Cell Type Co-localization (Spearman r)")
plt.savefig("colocalization_heatmap.png", dpi=150, bbox_inches="tight")
```

### Niche Composition Analysis

```python
# Cluster spots and analyze cell type composition per cluster
sc.pp.neighbors(st_adata, use_rep="X_pca")  # or spatial neighbors
sc.tl.leiden(st_adata, resolution=0.5, key_added="spatial_niche")

# Stacked bar chart of cell type composition per niche
props = st_adata.obsm["deconvolution"]
niche_composition = props.groupby(st_adata.obs["spatial_niche"]).mean()

niche_composition.plot(kind="bar", stacked=True, figsize=(12, 6), colormap="tab20")
plt.xlabel("Spatial Niche")
plt.ylabel("Mean Cell Type Proportion")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.savefig("niche_composition.png", dpi=150)
```

---

## Troubleshooting

### Common Issues

1. **Poor deconvolution**:
   - Increase spatial model epochs
   - Improve gene selection (more HVGs)
   - Check reference quality

2. **Missing cell types in output**:
   - Cell type not in reference
   - Low proportion below detection
   - Marker genes not captured

3. **All spots look similar**:
   - Training not converged (more epochs)
   - Gene overlap insufficient
   - Reference cell types too coarse

4. **Memory errors**:
   - Reduce number of genes
   - Process tissue regions separately
   - Use GPU acceleration

---

## References

- [Stereoscope Tutorial](https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/spatial/stereoscope_heart_LV_tutorial.html)
- [Stereoscope Paper](https://www.nature.com/articles/s42003-020-01247-y)
- [scvi-tools Documentation](https://docs.scvi-tools.org/)
