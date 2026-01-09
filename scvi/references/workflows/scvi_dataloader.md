# scVI Data Loader

This skill provides guidance for loading and preparing single-cell RNA-seq and multi-modal data for scvi-tools models (SCVI, SCANVI, TOTALVI, MULTIVI, CellAssign).

## When to Use This Skill

This skill should be used when:
- Loading single-cell data for scVI model training
- Preparing CITE-seq data for TOTALVI (RNA + protein)
- Setting up multi-modal data for MULTIVI (3+ modalities)
- Combining multiple batches or experimental datasets
- Troubleshooting "raw counts not found" or data registration errors
- Converting data formats (10x, CSV, loom) to AnnData for scvi-tools

---

## Quick Start: Workflow Decision Guide

| Data Type | Model | Key Setup |
|-----------|-------|-----------|
| Single RNA-seq | SCVI | `scvi.model.SCVI.setup_anndata(adata, layer='counts')` |
| RNA + labels | SCANVI | `scvi.model.SCANVI.setup_anndata(adata, labels_key='cell_type')` |
| RNA + protein (CITE-seq) | TOTALVI | `scvi.model.TOTALVI.setup_anndata(adata, protein_expression_obsm_key='protein')` |
| RNA + ATAC | MULTIVI | `scvi.model.MULTIVI.setup_mudata(mdata)` |
| RNA + marker genes | CellAssign | `CellAssign.setup_anndata(adata, size_factor_key='size_factor')` |

---

## Critical Data Requirements

### 1. Raw Counts (Most Important!)

All scvi-tools models require **raw, unnormalized integer counts**.

```python
# Check if data is raw counts
import numpy as np

if hasattr(adata.X, 'toarray'):
    sample = adata.X[:100].toarray()
else:
    sample = adata.X[:100]

is_integer = np.allclose(sample, sample.astype(int))
is_nonneg = sample.min() >= 0

print(f"Is raw counts: {is_integer and is_nonneg}")
```

**Common mistakes:**
- Using log-normalized data (values like 0.5, 1.2)
- Using scaled/z-score data (negative values)
- Using TPM/FPKM (non-integer)

### 2. Preserve Raw Counts in Layers

```python
# Always do this first after loading!
adata.layers['counts'] = adata.X.copy()

# Now you can normalize .X for visualization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# scvi-tools will use the preserved counts
scvi.model.SCVI.setup_anndata(adata, layer='counts')
```

---

## Complete Workflow: Basic RNA-seq → SCVI

```python
import scanpy as sc
import scvi
import numpy as np

# 1. Load data
adata = sc.read_h5ad('data.h5ad')
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

# 2. Preserve raw counts FIRST
adata.layers['counts'] = adata.X.copy()

# 3. Quality control
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# 4. Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# 5. Select highly variable genes (batch-aware if applicable)
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor='seurat_v3',
    layer='counts',
    batch_key='batch'  # if you have batches
)

# 6. Subset to HVGs (optional but recommended)
adata = adata[:, adata.var['highly_variable']].copy()

# 7. Setup for scVI
scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',
    batch_key='batch'  # for batch correction
)

# 8. Train model
model = scvi.model.SCVI(adata, n_latent=30, n_layers=2)
model.train(max_epochs=400, early_stopping=True)

# 9. Get latent representation
adata.obsm['X_scvi'] = model.get_latent_representation()

# 10. Use for downstream analysis
sc.pp.neighbors(adata, use_rep='X_scvi')
sc.tl.umap(adata)
sc.tl.leiden(adata)
```

---

## Complete Workflow: CITE-seq → TOTALVI

```python
import scanpy as sc
import scvi
import numpy as np
import pandas as pd

# 1. Load RNA data
adata = sc.read_h5ad('rna_data.h5ad')

# 2. Load protein data
protein_df = pd.read_csv('protein_data.csv', index_col=0)
# Ensure cells are aligned
protein_data = protein_df.loc[adata.obs_names].values

# 3. Preserve counts
adata.layers['counts'] = adata.X.copy()

# 4. Add protein to obsm
adata.obsm['protein_expression'] = protein_data

# 5. QC and HVG selection
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=4000,
    flavor='seurat_v3',
    layer='counts'
)

# 6. Setup for TOTALVI
scvi.model.TOTALVI.setup_anndata(
    adata,
    layer='counts',
    batch_key='batch',
    protein_expression_obsm_key='protein_expression'
)

# 7. Train
model = scvi.model.TOTALVI(adata)
model.train(max_epochs=400)

# 8. Get representations
adata.obsm['X_totalvi'] = model.get_latent_representation()
```

---

## Complete Workflow: Multi-modal → MULTIVI

```python
import scanpy as sc
import scvi
import mudata as md

# 1. Load modalities
rna = sc.read_h5ad('rna.h5ad')
atac = sc.read_h5ad('atac.h5ad')

# 2. Preserve counts in each
rna.layers['counts'] = rna.X.copy()
atac.layers['counts'] = atac.X.copy()

# 3. Find common cells
common_cells = list(set(rna.obs_names) & set(atac.obs_names))
rna = rna[common_cells].copy()
atac = atac[common_cells].copy()

# 4. Create MuData
mdata = md.MuData({'rna': rna, 'atac': atac})

# 5. Setup for MULTIVI
scvi.model.MULTIVI.setup_mudata(
    mdata,
    rna_layer='counts',
    atac_layer='counts',
    batch_key='batch'
)

# 6. Train
model = scvi.model.MULTIVI(mdata)
model.train()
```

---

## Complete Workflow: Combining Batches

```python
import scanpy as sc
import anndata as ad

# 1. Load datasets
datasets = []
batch_names = ['control', 'treatment1', 'treatment2']

for name in batch_names:
    adata = sc.read_h5ad(f'{name}.h5ad')
    adata.obs['batch'] = name
    adata.layers['counts'] = adata.X.copy()
    datasets.append(adata)

# 2. Check gene overlap
gene_sets = [set(d.var_names) for d in datasets]
common_genes = gene_sets[0].intersection(*gene_sets[1:])
print(f"Common genes: {len(common_genes)}")

# 3. Concatenate (inner join for common genes only)
adata = ad.concat(datasets, join='inner')
adata.var_names_make_unique()

# 4. Verify counts layer is preserved
if 'counts' not in adata.layers:
    adata.layers['counts'] = adata.X.copy()

# 5. Proceed with scVI setup
scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key='batch')
```

---

## Complete Workflow: CellAssign Annotation

```python
import scanpy as sc
import numpy as np
import pandas as pd
from scvi.external import CellAssign

# 1. Load data
adata = sc.read_h5ad('data.h5ad')

# 2. CRITICAL: Calculate size factors BEFORE subsetting
lib_size = np.array(adata.X.sum(axis=1)).flatten()
adata.obs['size_factor'] = lib_size / np.mean(lib_size)

# 3. Load marker matrix (binary: genes × cell types)
marker_matrix = pd.read_csv('markers.csv', index_col=0)
marker_matrix = (marker_matrix > 0).astype(int)

# 4. Subset to marker genes
common_genes = list(set(adata.var_names) & set(marker_matrix.index))
adata_subset = adata[:, common_genes].copy()
marker_matrix = marker_matrix.loc[common_genes]

# 5. Setup CellAssign
CellAssign.setup_anndata(adata_subset, size_factor_key='size_factor')

# 6. Train
model = CellAssign(adata_subset, marker_matrix)
model.train(max_epochs=400)

# 7. Get predictions
predictions = model.predict()
adata.obs['celltype'] = predictions.idxmax(axis=1).values
adata.obs['confidence'] = predictions.max(axis=1).values
```

---

## Common Errors and Solutions

### "raw counts not found"

```python
# Check data type
print(f"Max value: {adata.X.max()}")  # Should be large integer
print(f"Data type: {adata.X.dtype}")

# Look for counts in layers
print(f"Layers: {list(adata.layers.keys())}")

# If counts in layer, use it
scvi.model.SCVI.setup_anndata(adata, layer='counts')

# If raw attribute exists
if adata.raw is not None:
    adata = adata.raw.to_adata()
```

### "batch_key not in obs"

```python
# Check available columns
print(f"obs columns: {adata.obs.columns.tolist()}")

# Add batch if missing
adata.obs['batch'] = 'batch1'
```

### "AnnData not registered"

```python
# Must call setup_anndata BEFORE creating model
scvi.model.SCVI.setup_anndata(adata, layer='counts')
model = scvi.model.SCVI(adata)  # Now this works
```

### Size factor error (CellAssign)

```python
# Calculate on FULL data before subsetting
adata.obs['size_factor'] = adata.X.sum(1) / np.mean(adata.X.sum(1))

# Then subset to markers
adata_subset = adata[:, marker_genes].copy()
```

---

## Best Practices Summary

1. **Always preserve raw counts first** - `adata.layers['counts'] = adata.X.copy()`

2. **Use batch-aware HVG selection** - Include `batch_key` parameter

3. **Validate data before training:**
   ```python
   # Quick validation
   X = adata.layers.get('counts', adata.X)
   sample = X[:100].toarray() if hasattr(X, 'toarray') else X[:100]
   print(f"Is counts: {np.allclose(sample, sample.astype(int)) and sample.min() >= 0}")
   ```

4. **For CellAssign**: Calculate size factors on full data BEFORE gene subsetting

5. **For combining batches**: Use `join='inner'` to keep only common genes

6. **Monitor training**: Check loss convergence
   ```python
   model.history['elbo_train'].plot()
   ```

7. **Save models with data path documentation**

---

## Model Selection Guide

| Scenario | Recommended Model |
|----------|-------------------|
| Basic integration/batch correction | SCVI |
| Semi-supervised with some labels | SCANVI |
| CITE-seq (RNA + protein) | TOTALVI |
| Multi-modal (RNA + ATAC + more) | MULTIVI |
| Annotation with marker genes | CellAssign |
| Linear interpretable model | Linear SCVI |
| ATAC-seq only | PEAKVI |

---

## Quick Reference: Setup Commands

```python
# SCVI
scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key='batch')

# SCANVI
scvi.model.SCANVI.setup_anndata(adata, layer='counts', labels_key='cell_type', unlabeled_category='Unknown')

# TOTALVI
scvi.model.TOTALVI.setup_anndata(adata, layer='counts', protein_expression_obsm_key='protein')

# MULTIVI
scvi.model.MULTIVI.setup_mudata(mdata, rna_layer='counts', atac_layer='counts')

# CellAssign
CellAssign.setup_anndata(adata, size_factor_key='size_factor')
```
