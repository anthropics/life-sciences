# Spatial Deconvolution with DestVI

## Overview

DestVI (Deconvolution of Spatial Transcriptomics using Variational Inference) is a two-stage probabilistic method that:

1. **Stage 1 (scLVM)**: Learns cell-type-specific gene expression from scRNA-seq reference using CondSCVI
2. **Stage 2 (stLVM)**: Transfers this knowledge to spatial data for deconvolution

**Key Advantage**: Unlike discrete proportion methods, DestVI explicitly models **continuous variation within cell types** (gamma values), enabling detection of cell state gradients across tissue.

---

## Approach 1: Complete Analysis Pipeline (Recommended)

For standard DestVI deconvolution, use the convenience script `scripts/destvi_analysis.py`:

```bash
# Basic deconvolution
python3 scripts/destvi_analysis.py --reference scrna_ref.h5ad --spatial visium.h5ad --labels-key cell_type

# With more HVGs and longer training
python3 scripts/destvi_analysis.py --reference ref.h5ad --spatial spatial.h5ad --labels-key annotation --n-hvg 4000 --spatial-epochs 3000

# With weighted observations for imbalanced reference
python3 scripts/destvi_analysis.py --reference ref.h5ad --spatial spatial.h5ad --labels-key cell_type --weight-obs
```

**Parameters:**

Core parameters:
- `--reference` - scRNA-seq reference .h5ad (required)
- `--spatial` - Spatial transcriptomics .h5ad (required)
- `--labels-key` - Cell type column in reference (required)

Gene selection:
- `--n-hvg` - Number of HVGs (default: 3000)

Training:
- `--ref-epochs` - Reference model epochs (default: 300)
- `--spatial-epochs` - Spatial model epochs (default: 2500, min 1000 recommended)
- `--weight-obs` - Weight cells by inverse frequency

**Outputs:**
- `reference_training.png` - Reference model loss
- `spatial_training.png` - Spatial model loss
- `proportions_spatial.png` - Spatial proportion maps
- `cell_type_proportions.csv` - Proportion matrix
- `reference_model/` - Saved reference model
- `spatial_model/` - Saved spatial model
- `*_deconvolved.h5ad` - Processed spatial data

---

## Approach 2: Modular Building Blocks (For Custom Workflows)

For custom workflows, use functions from `scripts/destvi_core.py`:

```python
import sys
sys.path.append('scripts/')
from destvi_core import (
    setup_reference_data,
    train_reference_model,
    setup_spatial_data,
    train_spatial_model,
    get_proportions,
    get_gamma_values,
    get_celltype_expression,
    filter_shared_genes,
)

# Custom workflow
sc_adata, st_adata = filter_shared_genes(sc_adata, st_adata)
setup_reference_data(sc_adata, layer='counts', labels_key='cell_type')
ref_model = train_reference_model(sc_adata, max_epochs=300)

setup_spatial_data(st_adata, layer='counts')
spatial_model = train_spatial_model(st_adata, ref_model, max_epochs=2500)

# Get cell-type-specific expression
proportions = get_proportions(spatial_model)
gamma = get_gamma_values(spatial_model)
```

**Available functions:**
- `setup_reference_data()` / `setup_spatial_data()` - Register data
- `train_reference_model()` / `train_spatial_model()` - Train models
- `get_proportions()` - Cell type proportions per spot
- `get_gamma_values()` - Intra-cell-type variation
- `get_celltype_expression()` - Cell-type-specific gene expression
- `filter_shared_genes()` - Filter to common genes

---

## Workflow Steps (Manual)

### Step 1: Environment Setup

```python
import numpy as np
import scanpy as sc
import scvi
from scvi.model import CondSCVI, DestVI
import matplotlib.pyplot as plt
import torch

scvi.settings.seed = 0
sc.set_figure_params(figsize=(6, 6), frameon=False)
torch.set_float32_matmul_precision("high")
```

---

### Step 2: Load Reference scRNA-seq Data

```python
# Load annotated scRNA-seq reference
sc_adata = sc.read_h5ad("path/to/scrna_reference.h5ad")

print(f"Reference: {sc_adata.n_obs} cells, {sc_adata.n_vars} genes")
print(f"Cell types: {sc_adata.obs['cell_type'].value_counts()}")

# Ensure raw counts
if 'counts' not in sc_adata.layers:
    sc_adata.layers['counts'] = sc_adata.X.copy()
```

**Reference Requirements**:
- Raw counts (not normalized)
- Cell type annotations in `adata.obs`
- Ideally from same tissue/species as spatial data

---

### Step 3: Load Spatial Transcriptomics Data

```python
# Load Visium or other spatial data
st_adata = sc.read_h5ad("path/to/spatial_data.h5ad")

# Or load from 10x Space Ranger output
# st_adata = sc.read_visium("path/to/spaceranger/outs")

print(f"Spatial: {st_adata.n_obs} spots, {st_adata.n_vars} genes")

# Ensure raw counts
if 'counts' not in st_adata.layers:
    st_adata.layers['counts'] = st_adata.X.copy()
```

---

### Step 4: Preprocessing and Gene Selection

```python
# Filter genes present in both datasets
shared_genes = np.intersect1d(sc_adata.var_names, st_adata.var_names)
print(f"Shared genes: {len(shared_genes)}")

sc_adata = sc_adata[:, shared_genes].copy()
st_adata = st_adata[:, shared_genes].copy()

# Optional: Select highly variable genes for efficiency
sc.pp.highly_variable_genes(
    sc_adata,
    n_top_genes=3000,
    subset=False,
    layer="counts",
    flavor="seurat_v3"
)

# Use HVGs for both datasets
hvg = sc_adata.var_names[sc_adata.var.highly_variable]
sc_adata = sc_adata[:, hvg].copy()
st_adata = st_adata[:, hvg].copy()

print(f"After HVG selection: {len(hvg)} genes")
```

---

### Step 5: Train Reference Model (CondSCVI)

```python
# Setup reference data
CondSCVI.setup_anndata(
    sc_adata,
    layer="counts",
    labels_key="cell_type"  # Cell type annotation column
)

# Initialize and train scLVM (reference model)
sc_model = CondSCVI(
    sc_adata,
    weight_obs=False  # Equal weight to all cells
)

sc_model.train(max_epochs=300)

# Plot training
plt.figure(figsize=(8, 4))
plt.plot(sc_model.history['elbo_train'].values)
plt.xlabel('Epoch')
plt.ylabel('ELBO')
plt.title('Reference Model Training')
plt.show()

# Save reference model
sc_model.save("destvi_reference_model", overwrite=True)
```

---

### Step 6: Train Spatial Model (DestVI)

```python
# Setup spatial data
DestVI.setup_anndata(st_adata, layer="counts")

# Initialize stLVM from reference model
st_model = DestVI.from_rna_model(st_adata, sc_model)

# Train spatial deconvolution model
st_model.train(max_epochs=2500)  # Minimum 1000 recommended

# Plot training
plt.figure(figsize=(8, 4))
plt.plot(st_model.history['elbo_train'].values)
plt.xlabel('Epoch')
plt.ylabel('ELBO')
plt.title('Spatial Model Training')
plt.show()

# Save spatial model
st_model.save("destvi_spatial_model", overwrite=True)
```

---

### Step 7: Extract Cell Type Proportions

```python
# Get deconvolution results
st_adata.obsm["proportions"] = st_model.get_proportions()

# Add to obs for easy plotting
for ct in st_adata.obsm["proportions"].columns:
    st_adata.obs[f"prop_{ct}"] = st_adata.obsm["proportions"][ct]

# Visualize proportions spatially
cell_types = st_adata.obsm["proportions"].columns.tolist()
sc.pl.spatial(
    st_adata,
    color=[f"prop_{ct}" for ct in cell_types[:6]],  # First 6 types
    ncols=3,
    cmap="viridis",
    vmax=1.0
)
```

---

### Step 8: Extract Intra-Cell-Type Variation (Gamma Values)

```python
# Get continuous cell state variation per cell type
gamma_dict = st_model.get_gamma()

# Store gamma values in obsm
for ct, gamma in gamma_dict.items():
    st_adata.obsm[f"{ct}_gamma"] = gamma
    # First principal component of gamma for visualization
    st_adata.obs[f"{ct}_gamma_PC1"] = gamma[:, 0]

# Visualize cell state gradients
sc.pl.spatial(
    st_adata,
    color=[f"{ct}_gamma_PC1" for ct in list(gamma_dict.keys())[:4]],
    ncols=2,
    cmap="coolwarm"
)
```

---

### Step 9: Get Cell-Type-Specific Expression

```python
# Get deconvolved expression for specific cell type
cell_type = "Macrophage"  # Replace with your cell type

# Only analyze spots with sufficient proportion
proportion_threshold = 0.1
mask = st_adata.obs[f"prop_{cell_type}"] > proportion_threshold
indices = np.where(mask)[0]

# Get cell-type-specific expression
ct_expression = st_model.get_scale_for_ct(
    cell_type,
    indices=indices
)

# Store for analysis
st_adata.layers[f"{cell_type}_expression"] = np.zeros((st_adata.n_obs, st_adata.n_vars))
st_adata.layers[f"{cell_type}_expression"][indices, :] = ct_expression.values
```

---

### Step 10: Save Results

```python
# Save deconvolved spatial data
st_adata.write_h5ad("spatial_deconvolved.h5ad")

# Export proportions as CSV
st_adata.obsm["proportions"].to_csv("cell_type_proportions.csv")

# Reload models later
# sc_model = CondSCVI.load("destvi_reference_model", adata=sc_adata)
# st_model = DestVI.load("destvi_spatial_model", adata=st_adata)
```

---

## Key Parameters Reference

| Parameter | Default | Description |
|-----------|---------|-------------|
| `weight_obs` | False | Weight cells by inverse frequency (use for imbalanced data) |
| `max_epochs` (CondSCVI) | 300 | Reference model training epochs |
| `max_epochs` (DestVI) | 2500 | Spatial model epochs (min 1000) |
| `vamp_prior_p` | - | Controls gamma discretization |
| `l1_sparsity` | - | Increases proportion sparsity |

---

## Interpretation Guidelines

### Proportion Thresholds
- **>0.1**: Likely present
- **>0.3**: Dominant cell type
- **<0.05**: May be noise

### Gamma Values
- Capture continuous variation WITHIN a cell type
- Use PCA or clustering on gamma to identify cell states
- Compare gamma patterns across tissue regions

---

## Troubleshooting

### Common Issues

1. **Poor deconvolution (all spots similar)**:
   - Train longer (>2500 epochs)
   - Use more HVGs (4000-5000)
   - Check reference quality and cell type annotations

2. **Missing cell types**:
   - Ensure reference includes all expected types
   - Check gene overlap between datasets

3. **Noisy proportions**:
   - Increase `l1_sparsity` for sparser results
   - Apply proportion threshold before interpretation

4. **Memory errors**:
   - Reduce number of genes
   - Train on GPU

---

## References

- [DestVI Tutorial](https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/spatial/DestVI_tutorial.html)
- [DestVI Paper](https://www.nature.com/articles/s41587-022-01272-8)
- [scvi-tools Documentation](https://docs.scvi-tools.org/)
