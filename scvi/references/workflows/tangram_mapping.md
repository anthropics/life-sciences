# Spatial Mapping with Tangram

## Overview

Tangram is a spatial mapping tool that aligns single-cell and spatial transcriptomics data by learning an optimal transport mapping. It can:

- **Map cell type annotations** from scRNA-seq to spatial spots
- **Project gene expression** from single-cell to spatial resolution
- **Impute unmeasured genes** in spatial data using scRNA-seq
- **Operate in two modes**: "cells" (probabilistic) and "constrained" (fixed cell count)

---

## Workflow Steps

### Step 1: Environment Setup

```python
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import squidpy as sq
import mudata
import matplotlib.pyplot as plt
from scvi.external import Tangram

scvi.settings.seed = 0
sc.set_figure_params(figsize=(6, 6))
```

---

### Step 2: Load Data

```python
# Load spatial transcriptomics data (e.g., Visium)
adata_sp = sc.read_h5ad("path/to/spatial_data.h5ad")
# Or from Space Ranger: adata_sp = sc.read_visium("path/to/spaceranger/outs")

# Load scRNA-seq reference data
adata_sc = sc.read_h5ad("path/to/scrna_reference.h5ad")

print(f"Spatial: {adata_sp.n_obs} spots, {adata_sp.n_vars} genes")
print(f"scRNA-seq: {adata_sc.n_obs} cells, {adata_sc.n_vars} genes")
print(f"Cell types: {adata_sc.obs['cell_type'].value_counts()}")

# Filter genes
sc.pp.filter_genes(adata_sp, min_cells=1)
sc.pp.filter_genes(adata_sc, min_cells=1)
```

---

### Step 3: Select Marker Genes

```python
# Find marker genes for mapping (critical for accuracy)
sc.tl.rank_genes_groups(adata_sc, groupby="cell_type", use_raw=False)

# Get top markers per cell type
markers_df = pd.DataFrame(adata_sc.uns["rank_genes_groups"]["names"]).iloc[:100, :]
genes_sc = np.unique(markers_df.melt().value.values)

# Intersect with spatial genes
genes_st = adata_sp.var_names.values
genes = list(set(genes_sc).intersection(set(genes_st)))

print(f"Marker genes for mapping: {len(genes)}")

# Subset to shared genes
adata_sc_sub = adata_sc[:, genes].copy()
adata_sp_sub = adata_sp[:, genes].copy()
```

---

### Step 4: Calculate Density Priors

```python
# Option 1: If cell counts per spot are known
if "cell_count" in adata_sp.obs.columns:
    target_count = adata_sp.obs["cell_count"].sum()
    adata_sp.obs["density_prior"] = adata_sp.obs["cell_count"] / target_count

# Option 2: Use RNA counts as proxy for cell density
else:
    rna_count_per_spot = np.asarray(adata_sp.X.sum(axis=1)).squeeze()
    adata_sp.obs["density_prior"] = rna_count_per_spot / rna_count_per_spot.sum()
    target_count = adata_sc.n_obs  # Assume mapping all cells

# Option 3: Uniform prior (no density information)
# adata_sp.obs["density_prior"] = np.ones(adata_sp.n_obs) / adata_sp.n_obs
```

---

### Step 5: Create MuData Object

```python
# Tangram requires MuData with spatial and single-cell modalities
mdata = mudata.MuData({
    "sp": adata_sp_sub,
    "sc": adata_sc_sub
})

print(mdata)
```

---

### Step 6: Setup and Train Tangram

```python
# Setup Tangram
Tangram.setup_mudata(
    mdata,
    density_prior_key="density_prior",
    modalities={
        "density_prior_key": "sp",
        "sc_layer": "sc",
        "sp_layer": "sp",
    }
)

# Initialize model
# constrained=True: maps fixed number of cells (use target_count)
# constrained=False: probabilistic mapping ("cells" mode)
model = Tangram(
    mdata,
    constrained=True,
    target_count=target_count
)

# Train
model.train()

# Plot training
plt.figure(figsize=(8, 4))
plt.plot(model.history["train_loss"])
plt.xlabel("Iteration")
plt.ylabel("Loss")
plt.title("Tangram Training")
plt.show()
```

---

### Step 7: Get Mapper Matrix

```python
# Extract cell-to-spot mapping matrix
mapper = model.get_mapper_matrix()

# Shape: (n_cells, n_spots)
print(f"Mapper shape: {mapper.shape}")

# Store in scRNA-seq object
mdata.mod["sc"].obsm["tangram_mapper"] = mapper
```

---

### Step 8: Project Cell Type Annotations

```python
# Map cell type labels to spatial locations
labels = mdata.mod["sc"].obs["cell_type"]

mdata.mod["sp"].obsm["tangram_ct_pred"] = model.project_cell_annotations(
    mdata.mod["sc"],
    mdata.mod["sp"],
    mapper,
    labels
)

# Add to obs for visualization
ct_predictions = mdata.mod["sp"].obsm["tangram_ct_pred"]
for ct in ct_predictions.columns:
    adata_sp.obs[f"tangram_{ct}"] = ct_predictions[ct].values

# Visualize cell type predictions spatially
cell_types = ct_predictions.columns.tolist()
sc.pl.spatial(
    adata_sp,
    color=[f"tangram_{ct}" for ct in cell_types[:6]],
    ncols=3,
    cmap="viridis"
)
```

---

### Step 9: Project Gene Expression (Imputation)

```python
# Impute genes from scRNA-seq to spatial
# This can add genes not measured in spatial data
adata_sp_projected = model.project_genes(
    mdata.mod["sc"],
    mdata.mod["sp"],
    mapper
)

# Store as new modality or layer
mdata.mod["sp_projected"] = adata_sp_projected

# Visualize imputed genes
genes_to_plot = ["Gene1", "Gene2", "Gene3"]  # Replace with genes of interest
sc.pl.spatial(
    adata_sp_projected,
    color=genes_to_plot,
    ncols=3,
    cmap="plasma"
)
```

---

### Step 10: Save Results

```python
# Save MuData
mdata.write("tangram_results.h5mu")

# Save spatial data with annotations
adata_sp.write_h5ad("spatial_mapped.h5ad")

# Save mapper matrix
np.save("tangram_mapper.npy", mapper)

# Export cell type predictions
ct_predictions.to_csv("spatial_celltype_predictions.csv")

# Reload later
# mdata = mudata.read("tangram_results.h5mu")
```

---

## Key Parameters Reference

| Parameter | Default | Description |
|-----------|---------|-------------|
| `constrained` | True | Fixed cell count (True) vs probabilistic (False) |
| `target_count` | - | Total cells to map (required if constrained=True) |
| `density_prior_key` | - | Column with spot density priors |
| Marker genes | - | 100-200 markers recommended |

---

## Mode Selection

### Constrained Mode (`constrained=True`)
- Maps exact `target_count` cells to spatial locations
- Use when cell counts per spot are known
- More interpretable proportions

### Cells Mode (`constrained=False`)
- Probabilistic mapping without fixed count
- Use when cell density is uncertain
- More flexible, may overestimate dense regions

---

## Troubleshooting

### Common Issues

1. **Poor mapping quality**:
   - Increase marker genes (use top 150-200 per type)
   - Ensure marker genes are expressed in spatial data
   - Check that cell types match expected tissue composition

2. **Cell types not detected**:
   - Reference may lack relevant cell types
   - Marker genes may not be captured spatially
   - Try different marker selection methods

3. **Memory errors**:
   - Subsample scRNA-seq reference
   - Reduce number of genes
   - Process regions separately

4. **Uniform predictions**:
   - Improve density prior (use cell segmentation if available)
   - Increase training iterations

---

## References

- [Tangram Tutorial](https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/spatial/tangram_scvi_tools.html)
- [Tangram Paper](https://www.nature.com/articles/s41592-021-01264-7)
- [scvi-tools Documentation](https://docs.scvi-tools.org/)
- [Squidpy Documentation](https://squidpy.readthedocs.io/)
