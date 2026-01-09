# MrVI (Multi-resolution Variational Inference) Guide

Comprehensive reference for multi-sample single-cell analysis using MrVI.

## Table of Contents

1. [Method Overview](#method-overview)
2. [When to Use MrVI](#when-to-use-mrvi)
3. [Mathematical Foundations](#mathematical-foundations)
4. [Data Preparation](#data-preparation)
5. [Understanding the Outputs](#understanding-the-outputs)
6. [Interpreting Results](#interpreting-results)
7. [Advanced Analysis](#advanced-analysis)
8. [Troubleshooting](#troubleshooting)
9. [Comparison with Other Methods](#comparison-with-other-methods)

---

## Method Overview

MrVI (Multi-resolution Variational Inference) is a deep generative model designed for analyzing single-cell RNA-seq data from multiple samples. Unlike batch correction methods that aim to remove sample effects, MrVI explicitly models and leverages sample-level variation.

### Key Features

- **Multi-resolution inference**: Learns both sample-independent (u) and sample-specific (z) representations
- **Single-cell resolution**: Maintains cell-level granularity for all analyses
- **Sample distance computation**: Quantifies sample relationships within cell populations
- **Covariate-linked DE/DA**: Directly models how sample-level covariates affect expression

### Two Latent Spaces

| Space | Symbol | Description |
|-------|--------|-------------|
| Sample-independent | **u** | Captures broad cell states shared across samples |
| Sample-augmented | **z** | Combines u with sample-specific effects |

---

## When to Use MrVI

### Ideal Use Cases

1. **Multi-sample comparisons**
   - Comparing patients (healthy vs disease)
   - Treatment response studies
   - Developmental time courses across individuals

2. **Sample relationship analysis**
   - Clustering samples by transcriptional similarity
   - Identifying outlier samples
   - Understanding sample heterogeneity within cell types

3. **Differential analysis**
   - Cell-level differential expression
   - Compositional/abundance changes

### MrVI vs Other Methods

| Task | Recommended Method |
|------|-------------------|
| Batch correction / integration | scVI, scANVI, Harmony |
| Sample-level variation analysis | **MrVI** |
| Query-to-reference mapping | scANVI |
| Simple DE between groups | Standard methods (DESeq2, edgeR) |
| Cell-level DE with covariates | **MrVI** |

### When NOT to Use MrVI

- Single-sample datasets
- When batch effects dominate (use scVI first)
- When samples are incomparable (different tissues, species)
- Simple clustering without sample considerations

---

## Mathematical Foundations

### Generative Model

MrVI extends the scVI framework with explicit sample modeling:

```
For each cell n from sample s:
  u_n ~ Encoder(x_n)                    # Sample-independent state
  z_n = Augment(u_n, s)                 # Sample-augmented state
  x_n ~ Decoder(z_n, library_size)      # Reconstructed expression
```

### The u and z Representations

**u (Sample-independent)**:
- Captures cell-intrinsic biological variation
- Comparable across all samples
- Used for standard cell type clustering

**z (Sample-augmented)**:
- Incorporates sample-specific effects
- Used for within-sample analyses
- Enables sample distance computation

### Sample Distance Computation

For each pair of samples (i, j) within a cell population:

```
d(i, j) = E[||z_i - z_j||] for cells in population
```

This gives a sample × sample distance matrix per cell type.

---

## Data Preparation

### Required Data Structure

```python
adata.X              # Raw count matrix (cells × genes)
adata.obs['sample']  # Sample identifiers (required)
adata.obs['batch']   # Technical batch (optional, for correction)
adata.obs['condition']  # Sample-level covariate for DE/DA
```

### Preprocessing Checklist

1. **Quality Control**
   - Remove low-quality cells
   - Filter doublets
   - Standard QC metrics applied

2. **Feature Selection**
   - Select highly variable genes
   - MrVI benefits from more genes than scVI
   - Recommended: 5000-10000 HVGs

3. **Raw Counts**
   - Must have integer counts
   - Do NOT normalize before MrVI

### Example Preparation

```python
import scanpy as sc
import scvi
from scvi.external import MRVI

# Load data
adata = sc.read_h5ad('data.h5ad')

# Select HVGs (more than typical for scVI)
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=10000,
    flavor='seurat_v3',
    subset=True
)

# Setup for MrVI
MRVI.setup_anndata(
    adata,
    sample_key='patient_id',
    batch_key='sequencing_batch'  # optional
)
```

---

## Understanding the Outputs

### Latent Representation (u)

```python
u = model.get_latent_representation()
adata.obsm['X_mrvi_u'] = u

# Use for clustering and visualization
sc.pp.neighbors(adata, use_rep='X_mrvi_u')
sc.tl.umap(adata)
```

**Interpretation**:
- Similar to PCA/scVI embedding
- Cells cluster by type, not sample
- Suitable for standard downstream analysis

### Sample Distances

```python
distances = model.get_local_sample_distances(
    groupby='cell_type',
    keep_cell=False
)

# Access for specific cell type
cd8_dist = distances.loc[{"cell_type_name": "CD8 T cells"}]
```

**Output Structure**:
- xarray.DataArray with dimensions: [cell_type, sample_x, sample_y]
- Lower values = more similar samples
- Separate matrix per cell population

### Differential Expression

```python
de_results = model.differential_expression(
    sample_cov_keys=['condition'],
    store_lfc=True
)

# Effect sizes per cell
effect = de_results.effect_size.sel(covariate='condition_Disease')

# Log fold changes per gene
lfc = de_results.lfc.sel(covariate='condition_Disease')
```

**Output Structure**:
- `effect_size`: [covariate × cell] - how much each cell is affected
- `lfc`: [covariate × gene × cell] - per-gene changes per cell

### Differential Abundance

```python
da_results = model.differential_abundance(
    sample_cov_keys=['condition']
)

# Log probabilities
healthy_probs = da_results.condition_log_probs.loc[{"condition": "Healthy"}]
disease_probs = da_results.condition_log_probs.loc[{"condition": "Disease"}]

# Log ratio (positive = enriched in disease)
log_ratio = disease_probs - healthy_probs
```

---

## Interpreting Results

### Sample Distance Heatmaps

**What to look for**:
- Block structure indicates sample groups
- Cell-type-specific patterns reveal population-specific relationships
- Outliers appear as consistently distant rows/columns

**Example interpretation**:
```
CD8 T cells: Treatment samples cluster separately from control
Monocytes: No clear separation by treatment
→ Treatment effect is specific to CD8 T cells
```

### DE Effect Sizes

**Visualization**:
```python
adata.obs['effect'] = de_results.effect_size.sel(covariate='treatment_Drug').values
sc.pl.umap(adata, color='effect', cmap='viridis')
```

**Interpretation**:
- High effect = cell strongly responds to covariate
- Spatial patterns on UMAP reveal which populations are affected
- Can identify responding subpopulations within cell types

### DA Log Ratios

**Visualization**:
```python
adata.obs['da_lfc'] = log_ratio.values
sc.pl.umap(adata, color='da_lfc', cmap='coolwarm', vmin=-1, vmax=1)
```

**Interpretation**:
- Positive (red): Enriched in numerator condition
- Negative (blue): Depleted in numerator condition
- Identifies state shifts not captured by DE

---

## Advanced Analysis

### Hierarchical Clustering of Samples

```python
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

# Get distance matrix for a cell type
cd8_dist = distances.loc[{"cell_type_name": "CD8 T cells"}].values

# Hierarchical clustering
linkage_matrix = linkage(squareform(cd8_dist), method='ward')

# Plot dendrogram
fig, ax = plt.subplots(figsize=(10, 5))
dendrogram(linkage_matrix, labels=sample_names, ax=ax)
plt.title('Sample Clustering in CD8 T Cells')
```

### Identifying Cell-Type-Specific Effects

```python
# Compare DE effect sizes across cell types
cell_types = adata.obs['cell_type'].unique()

effects_by_type = {}
for ct in cell_types:
    mask = adata.obs['cell_type'] == ct
    effects_by_type[ct] = de_results.effect_size[:, mask].mean().values

# Plot
pd.DataFrame(effects_by_type).T.plot(kind='bar')
```

### Gene Set Enrichment on DE Genes

```python
# Extract top DE genes
from mrvi_core import extract_de_genes

top_genes = extract_de_genes(
    de_results,
    covariate='treatment_Drug',
    effect_threshold=0.5,
    top_n=100
)

# Run enrichment (using your favorite tool)
# e.g., gseapy, goatools, etc.
```

### Combining DE and DA

```python
# Identify cells with both expression and abundance changes
de_effect = de_results.effect_size.sel(covariate='condition_Disease').values
da_lfc = log_ratio.values

# Concordant changes (both up or both down)
concordant_up = (de_effect > 0.5) & (da_lfc > 0.5)
concordant_down = (de_effect < -0.5) & (da_lfc < -0.5)

# Discordant changes (DE without DA, or vice versa)
de_only = (np.abs(de_effect) > 0.5) & (np.abs(da_lfc) < 0.2)
da_only = (np.abs(de_effect) < 0.2) & (np.abs(da_lfc) > 0.5)
```

---

## Troubleshooting

### Problem: Samples Don't Show Expected Relationships

**Possible causes**:
- Wrong sample_key
- Sample variation dominated by batch effects
- Insufficient training

**Solutions**:
1. Verify sample_key has correct values
2. Add batch_key for technical confounders
3. Increase max_epochs
4. Check data quality per sample

### Problem: DE Shows No Significant Effects

**Possible causes**:
- Covariate has no real effect
- Insufficient samples per group
- Covariate improperly formatted

**Solutions**:
1. Check sample_info for covariate values
2. Ensure categorical covariates are properly typed
3. Reorder categories with reference first
4. Verify enough samples per condition (>3)

### Problem: Training is Very Slow

**Solutions**:
1. Use GPU (major speedup)
2. Reduce n_top_genes
3. Increase batch_size
4. Enable early_stopping

### Problem: Out of Memory for Sample Distances

**Solutions**:
1. Use `keep_cell=False` (default)
2. Reduce batch_size
3. Process cell types separately
4. Subsample large populations

---

## Comparison with Other Methods

### MrVI vs Pseudobulk DE

| Aspect | MrVI | Pseudobulk |
|--------|------|-----------|
| Resolution | Single-cell | Sample-level |
| Heterogeneity | Preserved | Lost |
| Statistical power | Per-cell | Per-sample |
| Sample size needs | Moderate | Strict (n>3) |

### MrVI vs scVI + Standard DE

| Aspect | MrVI | scVI + DE |
|--------|------|----------|
| Sample modeling | Explicit | Implicit |
| Sample distances | Built-in | Not available |
| DE interpretation | Covariate-linked | Group comparison |

### MrVI vs MILO (Differential Abundance)

| Aspect | MrVI | MILO |
|--------|------|------|
| Approach | Model-based | Graph-based |
| Cell resolution | Per-cell | Neighborhoods |
| Integration | Combined with DE | DA only |

---

## References

1. Boyeau, P., et al. (2024). Deep generative modeling of transcriptional dynamics for RNA velocity analysis in single cells. Nature Methods.

2. scvi-tools documentation: https://docs.scvi-tools.org/

3. Lopez, R., et al. (2018). Deep generative modeling for single-cell transcriptomics. Nature Methods.
