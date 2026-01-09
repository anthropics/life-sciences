# scVI/scANVI Integration Guide

Comprehensive reference for single-cell data integration using scvi-tools.

## Table of Contents

1. [Method Overview](#method-overview)
2. [Mathematical Foundations](#mathematical-foundations)
3. [Data Preparation](#data-preparation)
4. [Model Selection Guide](#model-selection-guide)
5. [Parameter Tuning](#parameter-tuning)
6. [Troubleshooting](#troubleshooting)
7. [Advanced Use Cases](#advanced-use-cases)
8. [Comparison with Other Methods](#comparison-with-other-methods)

---

## Method Overview

### scVI (Single-cell Variational Inference)

scVI is a deep generative model that learns a probabilistic representation of single-cell gene expression data. Key features:

- **Unsupervised learning**: No cell type annotations required
- **Batch correction**: Explicitly models batch effects as a latent variable
- **Scalability**: Handles millions of cells efficiently
- **Uncertainty quantification**: Provides uncertainty estimates for predictions

**Use cases:**
- Initial exploration without annotations
- Large-scale atlas construction
- Datasets where cell types are unknown or unreliable

### scANVI (Single-cell ANnotation using Variational Inference)

scANVI extends scVI by incorporating cell type annotations in a semi-supervised manner:

- **Semi-supervised**: Uses available labels while handling unlabeled cells
- **Label transfer**: Can predict labels for unlabeled cells
- **Better bio-conservation**: Supervised signal helps preserve biological variation

**Use cases:**
- When reliable cell type annotations exist (even partial)
- Query-to-reference mapping
- Cell type prediction for new datasets

---

## Mathematical Foundations

### Generative Model

scVI assumes gene expression follows a negative binomial distribution:

```
x_ng | z_n, s_n, l_n ~ NegativeBinomial(μ_ng, θ_g)

where:
  x_ng = expression of gene g in cell n
  z_n  = latent representation (learned embedding)
  s_n  = batch indicator
  l_n  = library size (total counts)
  μ_ng = mean expression (decoded from z, conditioned on s)
  θ_g  = gene-specific dispersion parameter
```

### Variational Autoencoder Architecture

**Encoder (Recognition Model):**
- Input: Gene expression vector x_n
- Output: Parameters of latent distribution q(z_n | x_n)
- Architecture: Fully connected neural network

**Decoder (Generative Model):**
- Input: Latent vector z_n, batch indicator s_n
- Output: Parameters of gene expression distribution
- Architecture: Fully connected neural network

### Loss Function

The model is trained by maximizing the Evidence Lower Bound (ELBO):

```
ELBO = E_q[log p(x|z,s)] - KL(q(z|x) || p(z))

where:
  - First term: reconstruction accuracy
  - Second term: regularization toward prior
```

---

## Data Preparation

### Required Data Structure

```python
# AnnData structure for scVI
adata.X           # Raw count matrix (cells × genes)
adata.obs         # Cell metadata (must include batch column)
adata.var         # Gene metadata

# Optional
adata.layers['counts']  # Counts if X is normalized
adata.obs['cell_type']  # For scANVI
```

### Pre-processing Checklist

1. **Quality Control (required before integration)**
   - Remove low-quality cells (high MT%, low genes)
   - Remove doublets
   - Filter genes expressed in too few cells

2. **Feature Selection**
   - Select highly variable genes (2000-4000 recommended)
   - Can use batch-aware HVG selection

3. **Data Format**
   - Must have raw counts (integers)
   - Do NOT normalize or log-transform for scVI input
   - scVI models the count distribution directly

### Example Preparation Code

```python
import scanpy as sc
import scvi

# Load data
adata = sc.read_h5ad('data.h5ad')

# Ensure counts are available
if 'counts' not in adata.layers:
    adata.layers['counts'] = adata.X.copy()

# Select highly variable genes
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    batch_key='batch',
    flavor='seurat_v3',
    layer='counts'
)
adata = adata[:, adata.var.highly_variable].copy()

# Register with scvi-tools
scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',
    batch_key='batch'
)
```

---

## Model Selection Guide

### Decision Tree

```
Do you have cell type annotations?
├── No → Use scVI
│   └── Are annotations available for some cells?
│       └── Yes → Consider scANVI with unlabeled_category
└── Yes → Are annotations reliable?
    ├── No (noisy/incomplete) → Use scVI, then scANVI cautiously
    └── Yes → Use scANVI
        └── Do you need to map new data?
            ├── No → Standard scANVI workflow
            └── Yes → Train on reference, then query mapping
```

### When to Use Each Model

| Scenario | Recommended Model |
|----------|------------------|
| No annotations | scVI |
| Partial annotations | scANVI (with unlabeled_category) |
| Full annotations | scANVI |
| Atlas building | scVI → scANVI |
| Query-to-reference | scANVI (transfer learning) |
| Uncertainty needed | scVI (for expression values) |

---

## Parameter Tuning

### Model Architecture

**n_latent (default: 30)**
- Controls embedding dimensionality
- Lower (10-20): Simpler representation, may lose detail
- Higher (50-100): More complex, risk of overfitting
- Guidance: Start with 30, increase for complex datasets

**n_layers (default: 2)**
- Number of hidden layers in encoder/decoder
- 1: Simple datasets, fewer samples
- 2: Standard choice
- 3-4: Very complex datasets, many batches

**gene_likelihood (default: 'nb')**
- 'nb': Negative binomial (standard choice)
- 'zinb': Zero-inflated NB (for highly sparse data)
- 'poisson': Simpler, less flexible

### Training Parameters

**max_epochs**
- Auto-determined based on dataset size
- Small datasets (<10k cells): ~400 epochs
- Large datasets (>100k cells): ~100-200 epochs
- Very large: May need fewer

**early_stopping (default: True)**
- Recommended to leave enabled
- Monitors validation loss
- Prevents overfitting

**batch_size (default: 128)**
- Larger: Faster, more memory
- Smaller: More updates, potentially better convergence
- GPU memory limited: Reduce to 64

### scANVI-Specific

**max_epochs for scANVI (default: 20)**
- Since initialized from scVI, needs fewer epochs
- Too many may overfit to labels

**unlabeled_category**
- String label for unannotated cells
- These cells will have labels predicted
- Set to match your data's "unknown" label

---

## Troubleshooting

### Problem: Batches Don't Mix Well

**Symptoms:**
- UMAP shows clear batch separation
- Low iLISI scores
- Batch silhouette score is high

**Solutions:**
1. Increase n_latent (more capacity for batch correction)
2. Check if batches have true biological differences
3. Try more training epochs
4. Verify batch key is correct
5. Consider if batches are too different (technology, species)

### Problem: Cell Types Split Into Multiple Clusters

**Symptoms:**
- Same cell type appears in multiple UMAP regions
- Low NMI/ARI scores
- Bio-conservation metrics are poor

**Solutions:**
1. Decrease n_latent (prevent over-separation)
2. Use scANVI with reliable annotations
3. Check annotation quality (may have subtypes)
4. Reduce HVG count
5. Verify cell types are actually the same across batches

### Problem: Training Is Slow

**Symptoms:**
- Training takes hours
- Progress bar moves slowly

**Solutions:**
1. Use GPU (10-100x speedup)
2. Increase batch_size (if memory allows)
3. Reduce max_epochs
4. Reduce dataset size for initial exploration
5. Use fewer HVGs

### Problem: Training Loss Not Decreasing

**Symptoms:**
- Loss plateaus early
- Model doesn't converge

**Solutions:**
1. Check data isn't already normalized
2. Increase max_epochs
3. Adjust learning rate via plan_kwargs
4. Check for data issues (NaN, negative values)

### Problem: Out of Memory

**Symptoms:**
- CUDA out of memory error
- Process killed

**Solutions:**
1. Reduce batch_size
2. Use fewer HVGs
3. Subsample data for testing
4. Use CPU (slower but more memory)

---

## Advanced Use Cases

### Query-to-Reference Mapping

Map new (query) data to an existing reference atlas:

```python
import scvi

# Load reference and trained model
reference = sc.read_h5ad('reference_atlas.h5ad')
scanvi_model = scvi.model.SCANVI.load('scanvi_model/', adata=reference)

# Prepare query data
query = sc.read_h5ad('query_data.h5ad')
scvi.model.SCANVI.prepare_query_anndata(query, scanvi_model)

# Train query model
query_model = scanvi_model.load_query_data(query)
query_model.train(max_epochs=100, plan_kwargs={'weight_decay': 0.0})

# Get predictions and embedding
query.obs['predicted_type'] = query_model.predict()
query.obsm['X_scANVI'] = query_model.get_latent_representation()
```

### Online Learning (Updating Models)

Update an existing model with new data:

```python
# Not directly supported, but can:
# 1. Concatenate new data with reference
# 2. Retrain from scratch
# 3. Or use query mapping approach above
```

### Differential Expression

scVI provides built-in DE analysis:

```python
# Train model
model = scvi.model.SCVI(adata, ...)
model.train()

# DE between groups
de_df = model.differential_expression(
    groupby='condition',
    group1='treatment',
    group2='control'
)

# Filter significant genes
significant = de_df[
    (de_df['is_de_fdr_0.05']) &
    (abs(de_df['lfc_mean']) > 0.5)
]
```

### Denoising and Imputation

```python
# Get denoised expression
denoised = model.get_normalized_expression(
    library_size=10000,  # Normalize to this depth
    return_numpy=True
)

# Impute dropout zeros
# (Not recommended - better to use scVI's posterior for downstream)
```

---

## Comparison with Other Methods

### scVI vs Harmony

| Aspect | scVI | Harmony |
|--------|------|---------|
| Type | Deep learning | Linear correction |
| Speed | Slower (GPU helps) | Fast |
| Scalability | Excellent | Excellent |
| Bio-conservation | Good-Excellent | Good |
| Batch correction | Excellent | Good |
| DE analysis | Built-in | Requires separate |
| Uncertainty | Yes | No |

**Use Harmony when:**
- Quick initial exploration
- Linear batch effects
- No GPU available

**Use scVI when:**
- Complex batch effects
- Need uncertainty estimates
- Want integrated DE analysis

### scVI vs BBKNN

| Aspect | scVI | BBKNN |
|--------|------|-------|
| Approach | Embedding | Graph |
| Preservation | Excellent | Good |
| Complexity | Higher | Lower |
| Interpretability | Latent space | Direct neighbors |

### scVI vs MNN (Mutual Nearest Neighbors)

| Aspect | scVI | MNN |
|--------|------|-----|
| Scalability | Better | Limited |
| Memory | Higher | Lower |
| Batch effects | All batches jointly | Pairwise |

### Benchmark Summary (from scib)

Based on scib benchmarks, scVI/scANVI typically rank among top methods for:
- Overall score (batch + bio)
- Scalability
- Biological preservation (especially scANVI)

Lower rankings for:
- Pure batch correction (some linear methods better)
- Speed (without GPU)

---

## References

1. Lopez, R., et al. (2018). Deep generative modeling for single-cell transcriptomics. Nature Methods.

2. Xu, C., et al. (2021). Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models. Molecular Systems Biology.

3. Luecken, M.D., et al. (2022). Benchmarking atlas-level data integration in single-cell genomics. Nature Methods.

4. scvi-tools documentation: https://docs.scvi-tools.org/
