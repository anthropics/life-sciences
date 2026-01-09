# Multi-Sample Analysis with MrVI

Automated workflow for multi-sample single-cell RNA-seq analysis using MrVI (Multi-resolution Variational Inference) from scvi-tools.

## When to Use This Skill

Use when users:
- Have multiple samples from comparable sources (same tissue, cell line, or experimental system)
- Want to analyze sample-level variation while preserving single-cell resolution
- Need to compare samples within specific cell populations
- Want differential expression analysis linked to sample-level covariates (e.g., disease status, treatment)
- Need differential abundance analysis to detect compositional changes
- Want to compute sample distances for hierarchical clustering or similarity analysis

**MrVI vs scVI/scANVI:**
- Use **scVI/scANVI** for batch correction and data integration
- Use **MrVI** for analyzing sample-level variation and covariate effects

**Supported input formats:**
- `.h5ad` files (AnnData format) with raw counts

**Prerequisites:**
- Data should be QC-filtered
- Sample identifiers in `adata.obs` column
- For DE/DA analysis: sample-level covariates (e.g., condition, treatment)
- Highly variable genes selected (typically 5000-10000 genes)

## Approach 1: Complete Analysis Pipeline (Recommended)

For standard MrVI analysis, use the convenience script `scripts/mrvi_analysis.py`:

```bash
# Basic analysis
python3 scripts/mrvi_analysis.py input.h5ad --sample-key patient_id

# With nuisance covariate (e.g., batch correction)
python3 scripts/mrvi_analysis.py input.h5ad --sample-key patient_id --batch-key site

# With differential expression/abundance analysis
python3 scripts/mrvi_analysis.py input.h5ad --sample-key donor --sample-cov-keys condition treatment
```

**When to use this approach:**
- Standard multi-sample analysis workflow
- User wants the "just works" solution
- Quick exploration of sample relationships

**Requirements:** scvi-tools, anndata, scanpy, scipy, matplotlib, seaborn, numpy, torch, xarray

**Parameters:**

Core parameters:
- `--sample-key` - Column in `adata.obs` containing sample identifiers (required)
- `--batch-key` - Column for nuisance covariates to correct (optional)
- `--output-dir` - Output directory (default: `<input_basename>_mrvi_results`)

Model architecture:
- `--n-latent` - Latent space dimensions (default: 30)

Training:
- `--max-epochs` - Maximum training epochs (default: 400)
- `--early-stopping` - Enable early stopping (default: True)

Downstream analysis:
- `--cell-type-key` - Column for cell type grouping in distance analysis
- `--sample-cov-keys` - Sample-level covariates for DE/DA analysis (space-separated)
- `--resolution` - Leiden clustering resolution (default: 1.0)
- `--min-dist` - UMAP minimum distance parameter (default: 0.3)

Feature selection:
- `--n-top-genes` - Number of HVGs if not pre-selected (default: 10000)

Use `--help` to see all options.

**Outputs:**

All files saved to `<input_basename>_mrvi_results/` (or `--output-dir`):
- `mrvi_umap.png` - UMAP visualization colored by sample and cell type
- `sample_distances/` - Sample distance matrices per cell population
- `sample_distance_heatmaps.png` - Hierarchical clustering of samples
- `de_results/` - Differential expression results (if covariates provided)
- `da_results/` - Differential abundance results (if covariates provided)
- `de_effect_umap.png` - UMAP colored by DE effect sizes
- `da_lfc_umap.png` - UMAP colored by DA log fold changes
- `<basename>_mrvi.h5ad` - Processed dataset with embeddings
- `mrvi_model/` - Saved MrVI model

### Workflow Steps

The script performs:

1. **Data Setup** - Register AnnData with MrVI, configure sample and batch keys
2. **Model Training** - Train MrVI to learn sample-aware representations
3. **Latent Extraction** - Extract u (sample-independent) embedding
4. **Downstream Analysis** - Build neighbor graph, compute UMAP, cluster
5. **Sample Distances** - Compute local sample distances per cell population
6. **Differential Expression** - Analyze covariate-linked expression changes (optional)
7. **Differential Abundance** - Detect compositional shifts (optional)
8. **Visualization** - Generate comprehensive plots

## Approach 2: Modular Building Blocks (For Custom Workflows)

For custom analysis workflows, use the modular functions from `scripts/mrvi_core.py`:

```python
import anndata as ad
import sys
sys.path.append('scripts/')
from mrvi_core import (
    setup_anndata_mrvi,
    train_mrvi_model,
    get_latent_representation,
    compute_sample_distances,
    run_differential_expression,
    run_differential_abundance
)

adata = ad.read_h5ad('input.h5ad')
setup_anndata_mrvi(adata, sample_key='patient_id')
# ... custom workflow
```

**When to use this approach:**
- Need custom distance analysis for specific cell populations
- Want to run DE/DA for specific comparisons only
- Integrating into a larger analysis pipeline
- Custom visualization requirements

**Available utility functions:**

From `mrvi_core.py`:
- `setup_anndata_mrvi(adata, sample_key, batch_key)` - Register data with MrVI
- `train_mrvi_model(adata, n_latent, max_epochs)` - Train MrVI model
- `get_latent_representation(model, adata)` - Extract u embedding
- `compute_sample_distances(model, groupby, keep_cell)` - Get sample distances
- `run_differential_expression(model, sample_cov_keys, store_lfc)` - DE analysis
- `run_differential_abundance(model, sample_cov_keys)` - DA analysis
- `save_model(model, path)` / `load_model(path, adata)` - Model persistence

**Example custom workflows:**

**Example 1: Basic sample distance analysis**
```python
adata = ad.read_h5ad('input.h5ad')
setup_anndata_mrvi(adata, sample_key='donor')
model = train_mrvi_model(adata, max_epochs=400)

# Get sample-independent embedding
adata.obsm['X_mrvi_u'] = get_latent_representation(model, adata)

# Compute distances grouped by cell type
distances = compute_sample_distances(model, groupby='cell_type')

# Access distances for specific population
cd8_distances = distances.loc[{"cell_type_name": "CD8 T cells"}]
```

**Example 2: Differential expression for specific comparison**
```python
# After training model...
de_results = run_differential_expression(
    model,
    sample_cov_keys=['treatment'],
    store_lfc=True
)

# Get effect sizes for treatment
treatment_effects = de_results.effect_size.sel(covariate='treatment_Treated')

# Get log fold changes for specific gene
gene_lfc = de_results.lfc.sel(gene='IL6')
```

**Example 3: Combined DE and DA analysis**
```python
# Run both analyses
de_results = run_differential_expression(model, ['disease_status'])
da_results = run_differential_abundance(model, ['disease_status'])

# Identify cells with both DE and DA changes
de_effect = de_results.effect_size.sel(covariate='disease_status_Disease').values
da_lfc = compute_da_log_ratio(da_results, 'disease_status', 'Disease', 'Healthy')

# Cells with concordant changes
concordant = (de_effect > 0.5) & (da_lfc > 0.5)
```

## Key Concepts

### MrVI Latent Representations

MrVI learns two complementary representations:

| Representation | Description | Use Case |
|---------------|-------------|----------|
| **u** | Sample-independent cell states | Cell type clustering, visualization |
| **z** | Sample-specific augmentation of u | Detailed sample comparisons |

### Sample Distances

MrVI computes distances between samples within cell populations:
- Reveals sample relationships specific to each cell type
- Can identify sample outliers or batch effects
- Useful for hierarchical clustering of samples

### Differential Expression (DE)

MrVI DE analysis differs from standard approaches:
- **Cell-level resolution**: Effect sizes computed per cell, not per cluster
- **Covariate-linked**: Directly models how covariates affect expression
- **Outputs**: Effect sizes and log fold changes per cell and gene

### Differential Abundance (DA)

Detects compositional changes between conditions:
- Identifies shifts in cell state frequencies
- Log probability ratios quantify enrichment direction
- Complements DE by capturing state-level changes

## Key Parameters Explained

### Model Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sample_key` | required | Column identifying samples (e.g., patient_id, donor) |
| `batch_key` | None | Nuisance covariate to correct (e.g., sequencing batch) |
| `n_latent` | 30 | Latent space dimensionality |

### Training Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_epochs` | 400 | Maximum training iterations |
| `early_stopping` | True | Stop when validation loss plateaus |

### Analysis Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `groupby` | None | Column for grouping sample distances |
| `keep_cell` | False | Return cell-level distances (memory intensive) |
| `sample_cov_keys` | None | Covariates for DE/DA (e.g., ['condition', 'treatment']) |
| `store_lfc` | True | Store log fold changes in DE results |

## Best Practices

1. **Data Requirements**
   - Multiple samples from comparable sources
   - Raw counts (not normalized)
   - Clear sample identifiers
   - Sample-level metadata for DE/DA

2. **Feature Selection**
   - Use more HVGs than scVI (5000-10000 recommended)
   - MrVI benefits from broader gene coverage

3. **Sample Key Selection**
   - Should represent biological replicates (patients, mice, wells)
   - Not technical replicates or batches

4. **Batch Key Usage**
   - Use for technical confounders (site, batch, plate)
   - Don't use for biological variables of interest

5. **Interpreting Results**
   - Sample distances: Lower = more similar within cell population
   - DE effect sizes: Magnitude indicates covariate impact
   - DA log ratios: Positive = enriched in numerator condition

6. **Common Issues**
   - If samples don't separate: Check sample_key correctness
   - If DE shows no signal: Verify covariate has variation
   - If training slow: Reduce n_top_genes or use GPU

## Reference Materials

For detailed methodology, parameter rationale, and advanced techniques, see `references/mrvi_guide.md`. This reference provides:
- Mathematical foundations of MrVI
- Detailed explanation of multi-resolution inference
- Advanced sample distance analysis
- Interpreting DE and DA results
- Comparison with other multi-sample methods

## Next Steps After MrVI Analysis

Typical downstream analysis:
- Hierarchical clustering of samples using distances
- Gene set enrichment on DE genes
- Cell state trajectory analysis
- Integration with clinical metadata
- Visualization of sample relationships
