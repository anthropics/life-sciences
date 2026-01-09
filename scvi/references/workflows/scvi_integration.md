# Single-Cell Data Integration with scvi-tools

Automated integration workflow for multi-sample single-cell RNA-seq data using scVI (unsupervised) and scANVI (semi-supervised) models.

## When to Use This Skill

Use when users:
- Need to integrate multiple scRNA-seq samples/batches
- Want to correct batch effects while preserving biological variation
- Ask to harmonize data from different experiments, donors, or technologies
- Want to transfer cell type annotations from a reference dataset
- Request atlas-level data integration
- Need benchmarking metrics for integration quality

**Supported input formats:**
- `.h5ad` files (AnnData format) with raw counts in `adata.X` or `adata.layers['counts']`

**Prerequisites:**
- Data should be QC-filtered (low-quality cells removed)
- Highly variable genes should be selected (typically 2000-4000 genes)
- Batch/sample information in `adata.obs` column
- For scANVI: cell type annotations in `adata.obs` column

**Default recommendation**: Use Approach 1 (complete pipeline) for standard integration workflows. Use Approach 2 for custom requirements or when more control is needed.

## Approach 1: Complete Integration Pipeline (Recommended)

For standard integration following scvi-tools best practices, use the convenience script `scripts/integration_analysis.py`:

```bash
python3 scripts/integration_analysis.py input.h5ad --batch-key batch
# With cell type annotations for scANVI:
python3 scripts/integration_analysis.py input.h5ad --batch-key batch --labels-key cell_type
```

**When to use this approach:**
- Standard multi-sample integration workflow
- User wants the "just works" solution
- Quick exploratory analysis before downstream processing

**Requirements:** scvi-tools, anndata, scanpy, scipy, matplotlib, seaborn, numpy, torch

**Parameters:**

Core parameters:
- `--batch-key` - Column in `adata.obs` containing batch/sample identifiers (required)
- `--labels-key` - Column with cell type annotations (enables scANVI; optional)
- `--output-dir` - Output directory (default: `<input_basename>_integration_results`)
- `--skip-scanvi` - Skip scANVI even if labels are available

Model architecture:
- `--n-layers` - Number of hidden layers (default: 2)
- `--n-latent` - Latent space dimensions (default: 30)
- `--gene-likelihood` - Distribution for gene expression: 'nb' (negative binomial), 'zinb', 'poisson' (default: 'nb')

Training:
- `--max-epochs-scvi` - Maximum training epochs for scVI (default: auto-determined)
- `--max-epochs-scanvi` - Maximum training epochs for scANVI (default: 20)
- `--early-stopping` - Enable early stopping (default: True)

scANVI-specific:
- `--unlabeled-category` - Label for unlabeled cells (default: 'Unknown')

Downstream analysis:
- `--resolution` - Leiden clustering resolution (default: 1.0)
- `--n-neighbors` - Number of neighbors for kNN graph (default: 15)
- `--min-dist` - UMAP minimum distance parameter (default: 0.3)

Use `--help` to see all options with current defaults.

**Outputs:**

All files saved to `<input_basename>_integration_results/` (or `--output-dir`):
- `scvi_latent_umap.png` - UMAP colored by batch and cluster (scVI embedding)
- `scanvi_latent_umap.png` - UMAP colored by batch, cluster, and cell type (scANVI embedding)
- `integration_comparison.png` - Side-by-side comparison of embeddings
- `integration_metrics.csv` - Quantitative benchmarking metrics
- `<basename>_integrated.h5ad` - Integrated dataset with embeddings
- `scvi_model/` - Saved scVI model (for reuse/transfer)
- `scanvi_model/` - Saved scANVI model (if trained)

### Workflow Steps

The script performs:

1. **Data Setup** - Register AnnData with scvi-tools, configure batch and label keys
2. **scVI Training** - Train unsupervised variational autoencoder for batch correction
3. **Latent Extraction** - Extract integrated low-dimensional representation
4. **scANVI Training** (optional) - Train semi-supervised model using cell type labels
5. **Downstream Analysis** - Build neighbor graph, compute UMAP, perform clustering
6. **Benchmarking** - Calculate integration quality metrics (batch mixing, bio-conservation)
7. **Visualization** - Generate comprehensive plots comparing embeddings

## Approach 2: Modular Building Blocks (For Custom Workflows)

For custom analysis workflows, use the modular functions from `scripts/integration_core.py` and `scripts/integration_metrics.py`:

```python
import anndata as ad
import sys
sys.path.append('scripts/')
from integration_core import (
    setup_anndata_scvi,
    train_scvi_model,
    train_scanvi_model,
    get_latent_representation,
    compute_neighbors_and_umap
)
from integration_metrics import calculate_integration_metrics

adata = ad.read_h5ad('input.h5ad')
setup_anndata_scvi(adata, batch_key='batch', layer='counts')
# ... custom workflow
```

**When to use this approach:**
- Need different training parameters for different samples
- Want to load a pre-trained model for query-to-reference mapping
- Integrating into a larger analysis pipeline
- Partial execution (only scVI, skip scANVI)
- Custom evaluation or visualization needs

**Available utility functions:**

From `integration_core.py`:
- `setup_anndata_scvi(adata, batch_key, layer, labels_key)` - Register data with scvi-tools
- `train_scvi_model(adata, n_layers, n_latent, gene_likelihood, max_epochs)` - Train scVI
- `train_scanvi_model(scvi_model, adata, labels_key, unlabeled_category, max_epochs)` - Train scANVI from scVI
- `get_latent_representation(model, adata)` - Extract integrated embedding
- `compute_neighbors_and_umap(adata, use_rep, n_neighbors, min_dist)` - Downstream processing
- `save_model(model, path)` / `load_model(path, adata)` - Model persistence
- `predict_labels(scanvi_model, adata)` - Predict cell types with scANVI

From `integration_metrics.py`:
- `calculate_integration_metrics(adata, batch_key, label_key, embed_key)` - Compute benchmarks
- `compare_embeddings(adata, embed_keys, batch_key, label_key)` - Compare multiple methods

**Example custom workflows:**

**Example 1: scVI only (no cell type annotations)**
```python
adata = ad.read_h5ad('input.h5ad')
setup_anndata_scvi(adata, batch_key='sample', layer='counts')
scvi_model = train_scvi_model(adata, n_latent=30)
adata.obsm['X_scVI'] = get_latent_representation(scvi_model, adata)
compute_neighbors_and_umap(adata, use_rep='X_scVI')
```

**Example 2: Query-to-reference mapping (label transfer)**
```python
# Load reference atlas with trained scANVI model
reference = ad.read_h5ad('reference_atlas.h5ad')
scanvi_model = scvi.model.SCANVI.load('reference_scanvi_model/', adata=reference)

# Map query to reference
query = ad.read_h5ad('query_data.h5ad')
scvi.model.SCANVI.prepare_query_anndata(query, scanvi_model)
query_model = scanvi_model.load_query_data(query)
query_model.train(max_epochs=100, plan_kwargs={'weight_decay': 0.0})

# Transfer labels
query.obs['predicted_cell_type'] = query_model.predict()
query.obsm['X_scANVI'] = query_model.get_latent_representation()
```

**Example 3: Compare multiple integration methods**
```python
# Assuming you have run multiple integration methods
# and stored results in adata.obsm['X_scVI'], adata.obsm['X_scANVI'], adata.obsm['X_harmony']
metrics_df = compare_embeddings(
    adata,
    embed_keys=['X_scVI', 'X_scANVI', 'X_harmony'],
    batch_key='batch',
    label_key='cell_type'
)
print(metrics_df)  # Shows metrics for each method
```

## Key Parameters Explained

### Model Architecture

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_layers` | 2 | Depth of encoder/decoder networks. Increase (3-4) for complex datasets |
| `n_latent` | 30 | Dimensionality of learned embedding. Higher values capture more variance |
| `gene_likelihood` | 'nb' | 'nb' (negative binomial) is standard; 'zinb' for high zero-inflation |

### Training Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_epochs` | auto | Set based on dataset size; ~400 for small, fewer for large datasets |
| `early_stopping` | True | Stop when validation loss plateaus |
| `batch_size` | 128 | Larger values speed training but use more memory |

### scANVI-Specific

| Parameter | Default | Description |
|-----------|---------|-------------|
| `unlabeled_category` | 'Unknown' | Label for cells without annotations (will be predicted) |
| `n_samples_per_label` | 100 | Cells sampled per label during training |

## Best Practices

1. **Data Preparation**
   - Filter low-quality cells before integration (QC should be done separately)
   - Select highly variable genes (2000-4000 recommended)
   - Ensure raw counts are available (not normalized data)
   - Check that batch key has meaningful groupings

2. **Model Training**
   - Start with defaults; adjust if integration quality is poor
   - Monitor training loss curves for convergence
   - Use GPU if available (10-100x faster)
   - Save models for reproducibility and reuse

3. **scANVI vs scVI Selection**
   - Use scVI when no cell type annotations exist
   - Use scANVI when reliable annotations are available (even partial)
   - scANVI typically produces better bio-conservation with same batch mixing

4. **Quality Assessment**
   - Always visualize UMAP colored by batch (should be mixed)
   - Check that known cell types cluster together
   - Use quantitative metrics (iLISI, cLISI, silhouette) for objective comparison
   - Compare with unintegrated data to assess improvement

5. **Common Issues**
   - If batches don't mix: increase `n_latent`, check for biological differences
   - If cell types split: decrease `n_latent`, check annotation quality
   - If training is slow: use GPU, reduce `max_epochs`, increase batch size

## Reference Materials

For detailed methodology, parameter rationale, and advanced techniques, see `references/scvi_integration_guide.md`. This reference provides:
- Mathematical foundations of the VAE models
- Detailed explanation of negative binomial likelihood
- Advanced transfer learning scenarios
- Troubleshooting guide for common integration problems
- Comparison with other integration methods (Harmony, BBKNN, MNN)

## Next Steps After Integration

Typical downstream analysis:
- Differential expression analysis (using scVI's built-in DE)
- Cell type annotation (manual or automated with CellTypist)
- Trajectory analysis (using integrated embedding)
- Cross-condition comparisons
- Spatial deconvolution (if spatial data available)
