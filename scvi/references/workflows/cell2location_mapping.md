# Spatial Cell Type Mapping with Cell2location

## Overview

Cell2location is a Bayesian statistical method that integrates 10X Visium spatial transcriptomics with scRNA-seq reference data to:

- **Estimate absolute cell abundances** (not just proportions) per spatial location
- **Account for technical variations** like platform effects and background signal
- **Model uncertainty** in cell type assignments
- **Handle complex mixtures** of cell types in each spot

**Two-Stage Workflow**:
1. **Regression Model**: Estimates cell-type-specific gene expression signatures from scRNA-seq
2. **Cell2location Model**: Maps signatures to spatial locations to infer cell abundances

---

## Workflow Steps

### Step 1: Environment Setup

```python
import numpy as np
import scanpy as sc
import scvi
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import torch
from cell2location.models import Cell2location, RegressionModel
from cell2location.utils.filtering import filter_genes

scvi.settings.seed = 0
sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")
```

---

### Step 2: Load scRNA-seq Reference Data

```python
# Load annotated scRNA-seq reference
adata_ref = sc.read_h5ad("path/to/scrna_reference.h5ad")

print(f"Reference: {adata_ref.n_obs} cells, {adata_ref.n_vars} genes")
print(f"Cell types:\n{adata_ref.obs['cell_type'].value_counts()}")

# Remove raw if present (Cell2location uses X)
if adata_ref.raw is not None:
    del adata_ref.raw
```

---

### Step 3: Load Spatial Transcriptomics Data

```python
# Load Visium spatial data
adata_vis = sc.read_h5ad("path/to/spatial_data.h5ad")

# Or from Space Ranger output
# adata_vis = sc.read_visium("path/to/spaceranger/outs")

# Add sample identifier
adata_vis.obs["sample"] = "sample_1"  # Modify for your data

print(f"Spatial: {adata_vis.n_obs} spots, {adata_vis.n_vars} genes")
```

---

### Step 4: Preprocess and Remove Mitochondrial Genes

```python
# Remove mitochondrial genes from spatial data (critical step!)
# MT genes are uninformative for cell type mapping and add noise
adata_vis.var['mt'] = adata_vis.var_names.str.startswith('MT-')
adata_vis = adata_vis[:, ~adata_vis.var['mt']].copy()

# Also remove from reference
adata_ref.var['mt'] = adata_ref.var_names.str.startswith('MT-')
adata_ref = adata_ref[:, ~adata_ref.var['mt']].copy()

print(f"After MT removal - Spatial: {adata_vis.n_vars}, Reference: {adata_ref.n_vars}")
```

---

### Step 5: Filter Genes for Cell2location

```python
# Apply Cell2location's gene filtering
# Permissive filtering to retain rare cell type markers
selected_genes = filter_genes(
    adata_ref,
    cell_count_cutoff=5,
    cell_percentage_cutoff2=0.03,
    nonz_mean_cutoff=1.12
)

adata_ref = adata_ref[:, selected_genes].copy()
print(f"Selected genes: {len(selected_genes)}")
```

---

### Step 6: Train Regression Model (Reference Signatures)

```python
# Setup reference for regression model
RegressionModel.setup_anndata(
    adata=adata_ref,
    batch_key="sample" if "sample" in adata_ref.obs.columns else None,
    labels_key="cell_type",  # Cell type column
    categorical_covariate_keys=["batch"] if "batch" in adata_ref.obs.columns else None
)

# Initialize regression model
ref_model = RegressionModel(adata_ref)

# Train to estimate cell type signatures
ref_model.train(
    max_epochs=250,
    batch_size=2500,
    train_size=1,
    lr=0.002
)

# Plot training
plt.figure(figsize=(8, 4))
plt.plot(ref_model.history['elbo_train'].values)
plt.xlabel('Epoch')
plt.ylabel('ELBO')
plt.title('Reference Signature Training')
plt.show()
```

---

### Step 7: Export Reference Signatures

```python
# Export posterior distributions of signatures
adata_ref = ref_model.export_posterior(
    adata_ref,
    sample_kwargs={"num_samples": 1000, "batch_size": 2500}
)

# Save reference model and data
ref_model.save("cell2location_reference", overwrite=True)
adata_ref.write_h5ad("cell2location_reference/reference.h5ad")

# Extract signature matrix for spatial mapping
if "means_per_cluster_mu_fg" in adata_ref.varm.keys():
    inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][
        [f"means_per_cluster_mu_fg_{i}" for i in adata_ref.uns["mod"]["factor_names"]]
    ].copy()
else:
    inf_aver = adata_ref.var[
        [f"means_per_cluster_mu_fg_{i}" for i in adata_ref.uns["mod"]["factor_names"]]
    ].copy()

inf_aver.columns = adata_ref.uns["mod"]["factor_names"]
print(f"Signature matrix: {inf_aver.shape}")
```

---

### Step 8: Prepare Spatial Data and Train Cell2location

```python
# Align genes between signature and spatial data
shared_genes = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, shared_genes].copy()
inf_aver = inf_aver.loc[shared_genes, :].copy()

print(f"Shared genes for mapping: {len(shared_genes)}")

# Setup spatial data
Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")

# Initialize Cell2location model
# N_cells_per_location: tissue-specific - adjust based on expected cell density
spatial_model = Cell2location(
    adata_vis,
    cell_state_df=inf_aver,
    N_cells_per_location=30,  # Expected cells per spot (tissue-dependent)
    detection_alpha=200  # RNA detection parameter
)

# Train spatial mapping model
spatial_model.train(
    max_epochs=30000,
    batch_size=None,  # Full batch
    train_size=1
)

# Plot training convergence
spatial_model.plot_history(1000)  # Skip first 1000 epochs
plt.show()
```

---

### Step 9: Export Results and Visualize

```python
# Export posterior estimates
adata_vis = spatial_model.export_posterior(
    adata_vis,
    sample_kwargs={"num_samples": 1000, "batch_size": adata_vis.n_obs}
)

# Save results
spatial_model.save("cell2location_spatial", overwrite=True)
adata_vis.write_h5ad("cell2location_spatial/spatial_mapped.h5ad")

# Add cell abundances to obs for plotting
# q05 = 5th percentile (conservative estimate)
cell_types = adata_vis.uns["mod"]["factor_names"]
adata_vis.obs[cell_types] = adata_vis.obsm["q05_cell_abundance_w_sf"]

# Visualize cell type abundances spatially
sc.pl.spatial(
    adata_vis,
    color=cell_types[:6],  # First 6 cell types
    ncols=3,
    cmap="magma",
    size=1.3,
    img_key="hires",
    vmin=0,
    vmax="p99.2"  # Cap at 99.2nd percentile
)
```

---

### Step 10: Advanced Visualization

```python
# Multi-panel visualization with custom styling
from cell2location.plt import plot_spatial

# Select cell types of interest
cell_types_to_plot = ["T_cells", "B_cells", "Macrophages"]  # Modify for your data

with mpl.rc_context({"figure.figsize": (15, 15)}):
    fig = plot_spatial(
        adata=adata_vis,
        color=cell_types_to_plot,
        labels=cell_types_to_plot,
        show_img=True,
        style="fast",
        max_color_quantile=0.992,
        circle_diameter=6,
        colorbar_position="right"
    )
    plt.show()
```

---

## Key Parameters Reference

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N_cells_per_location` | 30 | Expected cells per spot (tissue-dependent) |
| `detection_alpha` | 200 | RNA detection normalization (use 20 for high variation) |
| `max_epochs` (Regression) | 250 | Reference model training epochs |
| `max_epochs` (Cell2location) | 30000 | Spatial model training epochs |
| `batch_size` (Regression) | 2500 | Reference training batch size |
| `lr` | 0.002 | Learning rate for regression model |

---

## Tissue-Specific N_cells_per_location

| Tissue Type | Recommended Value |
|-------------|-------------------|
| Brain | 5-15 |
| Lymph node | 20-40 |
| Tumor | 10-30 |
| Heart | 5-15 |
| Liver | 10-20 |
| Intestine | 15-30 |

---

## Output Interpretation

### Cell Abundance Estimates
- `q05_cell_abundance_w_sf`: 5th percentile (conservative)
- `q50_cell_abundance_w_sf`: Median estimate
- `q95_cell_abundance_w_sf`: 95th percentile (upper bound)

### When to Use Each
- **q05**: Publication-quality, conservative claims
- **q50**: General analysis and visualization
- **q95**: Detecting rare populations

---

## Advanced Visualization

### Abundance vs Proportion Comparison

```python
# Cell2location gives abundances - compare with proportions
abundances = adata_vis.obsm["q50_cell_abundance_w_sf"]
proportions = abundances / abundances.sum(axis=1, keepdims=True)

cell_type = "T_cells"  # Replace with your cell type

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Absolute abundance
adata_vis.obs["_abundance"] = abundances[:, cell_types.index(cell_type)]
sc.pl.spatial(adata_vis, color="_abundance", ax=axes[0], show=False,
              title=f"{cell_type} Abundance", cmap="magma", vmax="p99", size=1.3)

# Relative proportion
adata_vis.obs["_proportion"] = proportions[:, cell_types.index(cell_type)]
sc.pl.spatial(adata_vis, color="_proportion", ax=axes[1], show=False,
              title=f"{cell_type} Proportion", cmap="viridis", vmax=1.0, size=1.3)

plt.tight_layout()
plt.savefig("abundance_vs_proportion.png", dpi=150)
```

### Uncertainty Visualization

```python
# Plot uncertainty using q05-q95 range
q05 = adata_vis.obsm["q05_cell_abundance_w_sf"]
q95 = adata_vis.obsm["q95_cell_abundance_w_sf"]
uncertainty = (q95 - q05) / (q95 + 1e-6)  # Relative uncertainty

cell_type = "T_cells"
ct_idx = cell_types.index(cell_type)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Median estimate
adata_vis.obs["_median"] = adata_vis.obsm["q50_cell_abundance_w_sf"][:, ct_idx]
sc.pl.spatial(adata_vis, color="_median", ax=axes[0], show=False,
              title=f"{cell_type} (Median)", cmap="magma", vmax="p99", size=1.3)

# Uncertainty
adata_vis.obs["_uncertainty"] = uncertainty[:, ct_idx]
sc.pl.spatial(adata_vis, color="_uncertainty", ax=axes[1], show=False,
              title=f"{cell_type} Uncertainty", cmap="YlOrRd", size=1.3)

plt.tight_layout()
plt.savefig("abundance_uncertainty.png", dpi=150)
```

### Multi-Cell Type Overlay

```python
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches

# Select 3 cell types for RGB overlay
ct1, ct2, ct3 = "T_cells", "B_cells", "Macrophages"  # Replace with yours

# Normalize to [0, 1] range
def normalize(x):
    return (x - x.min()) / (x.max() - x.min() + 1e-6)

rgb = np.zeros((adata_vis.n_obs, 3))
rgb[:, 0] = normalize(adata_vis.obs[ct1])  # Red
rgb[:, 1] = normalize(adata_vis.obs[ct2])  # Green
rgb[:, 2] = normalize(adata_vis.obs[ct3])  # Blue

# Plot with scatter
coords = adata_vis.obsm["spatial"]
fig, ax = plt.subplots(figsize=(10, 10))
ax.scatter(coords[:, 0], coords[:, 1], c=rgb, s=10, alpha=0.8)
ax.set_aspect("equal")
ax.invert_yaxis()
ax.set_title("Cell Type RGB Overlay")

# Legend
handles = [
    mpatches.Patch(color="red", label=ct1),
    mpatches.Patch(color="green", label=ct2),
    mpatches.Patch(color="blue", label=ct3),
]
ax.legend(handles=handles, loc="upper right")
plt.savefig("rgb_overlay.png", dpi=200)
```

### Signature Quality Check

```python
# Visualize gene expression signatures learned from reference
import seaborn as sns

# Get signature matrix
sig = inf_aver.copy()

# Top genes per cell type
n_top = 10
fig, ax = plt.subplots(figsize=(12, 8))

# Select top variable genes across signatures
gene_var = sig.var(axis=1)
top_genes = gene_var.nlargest(50).index

sns.clustermap(sig.loc[top_genes, :], cmap="viridis",
               figsize=(10, 12), z_score=0)
plt.suptitle("Cell Type Gene Signatures (Top Variable Genes)", y=1.02)
plt.savefig("reference_signatures.png", dpi=150, bbox_inches="tight")
```

---

## Troubleshooting

### Common Issues

1. **Poor convergence**:
   - Increase `max_epochs` for spatial model
   - Check ELBO plot for plateau
   - Adjust `N_cells_per_location`

2. **Missing cell types**:
   - Check reference contains expected types
   - Verify gene overlap
   - Lower detection threshold

3. **Noisy results**:
   - Use q05 estimates (more conservative)
   - Increase `detection_alpha` for high technical variation
   - Check spatial data quality

4. **Memory errors**:
   - Reduce `num_samples` in export
   - Process samples separately
   - Use GPU

---

## References

- [Cell2location Tutorial](https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/spatial/cell2location_lymph_node_spatial_tutorial.html)
- [Cell2location Paper](https://www.nature.com/articles/s41587-021-01139-4)
- [Cell2location Documentation](https://cell2location.readthedocs.io/)
