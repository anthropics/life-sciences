# Reference-to-Query Label Transfer with scANVI

Automated workflow for integrating datasets and transferring cell type annotations from a labeled reference to an unlabeled query using scVI and scANVI.

## When to Use This Skill

Use when users:
- Have a fully-annotated reference dataset and unannotated query data
- Need to integrate datasets from different sequencing technologies
- Want to transfer cell type labels from reference to query
- Need de novo integration (not using pre-trained models)
- Have cross-platform data (10x, SmartSeq2, Drop-seq, etc.)

**Label Transfer vs Other Methods:**
- Use **label transfer** when you have annotated reference + unannotated query
- Use **seed labeling** when you only have marker genes
- Use **scArches** when you have a pre-trained reference model

**Supported input formats:**
- `.h5ad` files (AnnData format) with raw counts

**Prerequisites:**
- Reference dataset with cell type annotations
- Query dataset (unannotated)
- Both datasets should have overlapping genes
- Raw counts (not normalized)

## Approach 1: Complete Label Transfer Pipeline (Recommended)

For standard label transfer workflow, use the convenience script `scripts/label_transfer_analysis.py`:

```bash
# Basic label transfer
python3 scripts/label_transfer_analysis.py \
    --reference reference.h5ad \
    --query query.h5ad \
    --labels-key cell_type

# With technology batch correction
python3 scripts/label_transfer_analysis.py \
    --reference reference.h5ad \
    --query query.h5ad \
    --labels-key cell_ontology_class \
    --batch-key technology

# With SmartSeq2 gene length normalization
python3 scripts/label_transfer_analysis.py \
    --reference smartseq2_data.h5ad \
    --query 10x_data.h5ad \
    --labels-key cell_type \
    --normalize-gene-length \
    --reference-tech SS2
```

**When to use this approach:**
- Standard reference-to-query annotation
- Cross-technology integration
- User wants automated pipeline

**Requirements:** scvi-tools, anndata, scanpy, numpy, pandas, torch

**Parameters:**

Core parameters:
- `--reference` - Path to annotated reference .h5ad file (required)
- `--query` - Path to unannotated query .h5ad file (required)
- `--labels-key` - Column in reference.obs containing cell type labels (required)
- `--output-dir` - Output directory (default: `label_transfer_results`)

Batch correction:
- `--batch-key` - Column for batch/technology correction
- `--reference-batch` - Batch label for reference data
- `--query-batch` - Batch label for query data

Technology-specific:
- `--normalize-gene-length` - Apply gene length normalization (for SmartSeq2)
- `--reference-tech` - Reference technology: '10x', 'SS2', 'other'
- `--query-tech` - Query technology: '10x', 'SS2', 'other'

Model architecture:
- `--n-latent` - Latent space dimensions (default: 30)
- `--n-layers` - Number of hidden layers (default: 2)
- `--n-top-genes` - Number of HVGs to select (default: 2000)

Training:
- `--max-epochs-scvi` - Max epochs for scVI (default: 400)
- `--max-epochs-scanvi` - Max epochs for scANVI (default: 20)

Use `--help` to see all options.

**Outputs:**

All files saved to `label_transfer_results/` (or `--output-dir`):
- `integration_umap.png` - UMAP showing integrated reference + query
- `label_transfer_umap.png` - UMAP colored by transferred labels
- `confusion_matrix.png` - Prediction accuracy heatmap
- `prediction_confidence.png` - Confidence distribution
- `transferred_labels.csv` - Predicted labels for query cells
- `integrated_data.h5ad` - Combined dataset with predictions
- `scvi_model/` - Saved scVI model
- `scanvi_model/` - Saved scANVI model

### Workflow Steps

The script performs:

1. **Data Loading** - Load reference and query datasets
2. **Gene Length Normalization** - Correct for SmartSeq2 if needed
3. **Concatenation** - Merge datasets with batch labels
4. **HVG Selection** - Select highly variable genes across batches
5. **scVI Training** - Learn integrated latent representation
6. **scANVI Training** - Train semi-supervised classifier
7. **Label Prediction** - Transfer labels to query cells
8. **Evaluation** - Compute confusion matrix and confidence
9. **Visualization** - Generate integration and transfer plots

## Approach 2: Modular Building Blocks (For Custom Workflows)

For custom analysis workflows, use the modular functions from `scripts/label_transfer_core.py`:

```python
import anndata as ad
import sys
sys.path.append('scripts/')
from label_transfer_core import (
    normalize_gene_length,
    concatenate_datasets,
    setup_combined_anndata,
    train_scvi_integration,
    train_scanvi_transfer,
    predict_labels,
    evaluate_transfer
)

reference = ad.read_h5ad('reference.h5ad')
query = ad.read_h5ad('query.h5ad')
# ... custom workflow
```

**When to use this approach:**
- Need custom preprocessing
- Want to use existing scVI model
- Integration with larger pipelines
- Custom evaluation metrics

**Available utility functions:**

From `label_transfer_core.py`:
- `normalize_gene_length(adata, gene_lengths)` - SmartSeq2 normalization
- `concatenate_datasets(reference, query, batch_key)` - Merge with labels
- `select_hvg_across_batches(adata, batch_key, n_top_genes)` - Batch-aware HVG
- `setup_combined_anndata(adata, batch_key, labels_key)` - Register for scVI
- `train_scvi_integration(adata, n_latent, n_layers)` - Train scVI
- `train_scanvi_transfer(scvi_model, adata, labels_key)` - Train scANVI
- `predict_labels(scanvi_model, adata)` - Get predictions
- `get_prediction_probabilities(scanvi_model, adata)` - Get confidence
- `evaluate_transfer(true_labels, predicted_labels)` - Compute metrics
- `save_model(model, path)` / `load_model(path, adata)` - Persistence

**Example custom workflows:**

**Example 1: Basic label transfer**
```python
reference = ad.read_h5ad('reference.h5ad')
query = ad.read_h5ad('query.h5ad')

# Add batch labels
reference.obs['batch'] = 'reference'
query.obs['batch'] = 'query'

# Concatenate
adata = concatenate_datasets(reference, query, labels_key='cell_type')

# Select HVGs
select_hvg_across_batches(adata, batch_key='batch', n_top_genes=2000)

# Train and transfer
setup_combined_anndata(adata, batch_key='batch', labels_key='cell_type_transfer')
scvi_model = train_scvi_integration(adata)
scanvi_model = train_scanvi_transfer(scvi_model, adata, 'cell_type_transfer')

# Get predictions for query
query_mask = adata.obs['batch'] == 'query'
predictions = predict_labels(scanvi_model, adata)
query_predictions = predictions[query_mask]
```

**Example 2: Cross-technology with gene length normalization**
```python
# Load SmartSeq2 reference
reference = ad.read_h5ad('smartseq2_reference.h5ad')
reference = normalize_gene_length(reference, gene_lengths_file='gene_lengths.txt')

# Load 10x query (no normalization needed)
query = ad.read_h5ad('10x_query.h5ad')

# Continue with standard workflow...
```

**Example 3: Transfer with confidence filtering**
```python
# After training...
predictions = predict_labels(scanvi_model, adata)
probabilities = get_prediction_probabilities(scanvi_model, adata)

# Filter low-confidence predictions
confidence = probabilities.max(axis=1)
high_conf_mask = confidence > 0.8

# Only keep high-confidence predictions
adata.obs['predicted_type'] = predictions
adata.obs.loc[~high_conf_mask, 'predicted_type'] = 'Low_confidence'
```

## Cross-Technology Integration

### Technology-Specific Considerations

| Technology | UMI-based | Gene Length Correction | Notes |
|------------|-----------|----------------------|-------|
| 10x Genomics | Yes | No | Standard workflow |
| SmartSeq2 | No | Yes | Full-length transcripts |
| Drop-seq | Yes | No | Similar to 10x |
| CEL-seq2 | Yes | No | Plate-based UMI |
| inDrop | Yes | No | Similar to 10x |

### Gene Length Normalization

For SmartSeq2 data, read counts are proportional to gene length:

```python
# Load gene lengths
gene_lengths = pd.read_csv('gene_lengths.txt', index_col=0)

# Normalize
adata.X = adata.X / gene_lengths.values * np.median(gene_lengths.values)
adata.X = np.rint(adata.X)  # Round to integers
```

### Batch Key Setup

When integrating different technologies:

```python
# Option 1: Use technology as batch
adata.obs['tech'] = ['10x'] * n_10x + ['SS2'] * n_ss2
setup_anndata(..., batch_key='tech')

# Option 2: Use original batch + technology
adata.obs['batch_tech'] = adata.obs['original_batch'] + '_' + adata.obs['tech']
setup_anndata(..., batch_key='batch_tech')
```

## Key Parameters Explained

### HVG Selection

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_top_genes` | 2000 | Number of highly variable genes |
| `batch_key` | None | Account for batch in HVG selection |
| `flavor` | 'seurat_v3' | HVG selection method |

**Note:** Performance degrades when gene count approaches cell count. 2000 HVGs is a good balance.

### Model Architecture

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_latent` | 30 | Latent space dimensions |
| `n_layers` | 2 | Network depth |

### Training

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_epochs_scvi` | 400 | Epochs for integration model |
| `max_epochs_scanvi` | 20 | Epochs for transfer model |
| `n_samples_per_label` | 100 | Balances label sampling |

## Best Practices

1. **Data Preparation**
   - Use raw counts (not normalized)
   - Ensure gene names match between datasets
   - Apply gene length correction for SmartSeq2

2. **Gene Selection**
   - Use batch-aware HVG selection
   - 2000 genes is typically sufficient
   - More genes may help for complex datasets

3. **Batch Correction**
   - Always specify batch_key for cross-technology
   - Consider using technology as primary batch

4. **Model Training**
   - Train scVI first for good integration
   - scANVI needs fewer epochs (initialized from scVI)
   - Monitor training loss convergence

5. **Evaluation**
   - Check confusion matrix for systematic errors
   - Filter low-confidence predictions
   - Validate with known markers

6. **Common Issues**
   - Poor integration: Check batch_key, increase epochs
   - Wrong predictions: Verify label quality in reference
   - Low confidence: May indicate novel cell types

## Reference Materials

For detailed methodology and troubleshooting, see `references/label_transfer_guide.md`. This provides:
- Mathematical foundations of label transfer
- Detailed cross-technology integration strategies
- Troubleshooting common transfer problems
- Comparison with other transfer methods

## Next Steps After Label Transfer

Typical downstream analysis:
- Validate predictions with marker expression
- Compare cell type proportions across datasets
- Differential expression between conditions
- Trajectory analysis using integrated embedding
- Export predictions for further analysis
