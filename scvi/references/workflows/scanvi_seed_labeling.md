# Semi-Supervised Cell Annotation with scANVI Seed Labeling

Automated workflow for cell type annotation using marker gene signatures and scANVI semi-supervised learning.

## When to Use This Skill

Use when users:
- Have marker gene signatures for cell types of interest
- Want to annotate cells without fully labeled reference data
- Have a small set of confidently labeled cells to extend
- Need semi-supervised cell type classification
- Want to leverage known biology (marker genes) for annotation

**Seed Labeling vs Other Annotation Methods:**
- Use **seed labeling** when you have marker genes but no reference atlas
- Use **scANVI query-to-reference** when you have a labeled reference dataset
- Use **CellTypist/manual** when you need fully manual curation

**Supported input formats:**
- `.h5ad` files (AnnData format) with raw counts

**Prerequisites:**
- Raw count data (not normalized)
- Marker gene signatures for each cell type (positive and negative markers)
- Genes in signatures must be present in the dataset

## Approach 1: Complete Annotation Pipeline (Recommended)

For standard seed labeling workflow, use the convenience script `scripts/seed_labeling_analysis.py`:

```bash
# Using a marker gene JSON file
python3 scripts/seed_labeling_analysis.py input.h5ad --markers-file markers.json

# With batch correction
python3 scripts/seed_labeling_analysis.py input.h5ad --markers-file markers.json --batch-key batch
```

**Marker Gene File Format (JSON):**

```json
{
  "CD8 T cell": {
    "positive": ["CD8A", "CD8B", "CCR7"],
    "negative": ["CD4"]
  },
  "CD4 T cell": {
    "positive": ["CD4", "CCR7"],
    "negative": ["CD8A", "CD8B"]
  },
  "Monocyte": {
    "positive": ["CD14", "LYZ", "S100A8"],
    "negative": ["CD3D", "CD3E"]
  }
}
```

**When to use this approach:**
- Standard cell annotation workflow
- User wants automated marker-based labeling
- Quick annotation with known cell type markers

**Requirements:** scvi-tools, anndata, scanpy, numpy, torch

**Parameters:**

Core parameters:
- `--markers-file` - JSON file with marker gene signatures (required)
- `--batch-key` - Column for batch correction (optional)
- `--output-dir` - Output directory (default: `<input_basename>_seed_labeling_results`)

Seed selection:
- `--n-seed-cells` - Number of seed cells per type (default: 50)
- `--min-score-percentile` - Minimum percentile for seed selection (default: 95)

Model architecture:
- `--n-latent` - Latent space dimensions (default: 30)
- `--n-layers` - Number of hidden layers (default: 2)

Training:
- `--max-epochs-scvi` - Max epochs for scVI (default: 100)
- `--max-epochs-scanvi` - Max epochs for scANVI (default: 25)

Downstream:
- `--resolution` - Leiden clustering resolution (default: 1.0)
- `--min-dist` - UMAP minimum distance (default: 0.3)

Use `--help` to see all options.

**Outputs:**

All files saved to `<input_basename>_seed_labeling_results/` (or `--output-dir`):
- `seed_labels_umap.png` - UMAP showing seed labels
- `predictions_umap.png` - UMAP with predicted cell types
- `comparison_umap.png` - Side-by-side seed vs predictions
- `prediction_confidence.png` - Prediction probability distribution
- `marker_scores.csv` - Marker gene scores per cell
- `<basename>_annotated.h5ad` - Annotated dataset
- `scvi_model/` - Saved scVI model
- `scanvi_model/` - Saved scANVI model

### Workflow Steps

The script performs:

1. **Marker Scoring** - Score cells based on marker gene expression
2. **Seed Selection** - Select top-scoring cells as seeds per cell type
3. **scVI Training** - Train base variational autoencoder
4. **scANVI Training** - Train semi-supervised classifier with seeds
5. **Prediction** - Predict labels for all unlabeled cells
6. **Confidence Assessment** - Compute prediction probabilities
7. **Visualization** - Generate annotation quality plots

## Approach 2: Modular Building Blocks (For Custom Workflows)

For custom analysis workflows, use the modular functions from `scripts/seed_labeling_core.py`:

```python
import anndata as ad
import sys
sys.path.append('scripts/')
from seed_labeling_core import (
    load_marker_genes,
    score_cells_by_markers,
    select_seed_cells,
    create_seed_labels,
    train_scvi_model,
    train_scanvi_from_seeds,
    predict_labels,
    get_prediction_confidence
)

adata = ad.read_h5ad('input.h5ad')
markers = load_marker_genes('markers.json')
# ... custom workflow
```

**When to use this approach:**
- Need custom seed selection criteria
- Want to combine multiple marker sources
- Iterative refinement of annotations
- Integration with existing pipelines

**Available utility functions:**

From `seed_labeling_core.py`:
- `load_marker_genes(path)` - Load markers from JSON/CSV
- `score_cells_by_markers(adata, markers)` - Compute marker scores
- `select_seed_cells(scores, n_cells, method)` - Select top-scoring cells
- `create_seed_labels(adata, seed_masks, cell_types)` - Assign seed labels
- `train_scvi_model(adata, n_latent, n_layers, max_epochs)` - Train base model
- `train_scanvi_from_seeds(scvi_model, adata, labels_key)` - Train classifier
- `predict_labels(scanvi_model, adata)` - Get predictions
- `get_prediction_confidence(scanvi_model, adata)` - Get probabilities

**Example custom workflows:**

**Example 1: Basic seed labeling with custom markers**
```python
adata = ad.read_h5ad('input.h5ad')

# Define markers inline
markers = {
    "T cell": {"positive": ["CD3D", "CD3E"], "negative": ["CD14"]},
    "B cell": {"positive": ["CD19", "MS4A1"], "negative": ["CD3D"]},
    "Monocyte": {"positive": ["CD14", "LYZ"], "negative": ["CD3D"]}
}

# Score and select seeds
scores = score_cells_by_markers(adata, markers)
seed_masks = select_seed_cells(scores, n_cells=50)
create_seed_labels(adata, seed_masks, list(markers.keys()))

# Train and predict
scvi_model = train_scvi_model(adata)
scanvi_model = train_scanvi_from_seeds(scvi_model, adata, 'seed_labels')
adata.obs['predicted_type'] = predict_labels(scanvi_model, adata)
```

**Example 2: Iterative refinement with confidence filtering**
```python
# First pass
scanvi_model = train_scanvi_from_seeds(scvi_model, adata, 'seed_labels')
predictions = predict_labels(scanvi_model, adata)
confidence = get_prediction_confidence(scanvi_model, adata)

# Filter high-confidence predictions as new seeds
high_conf_mask = confidence.max(axis=1) > 0.9
adata.obs['refined_seeds'] = 'Unknown'
adata.obs.loc[high_conf_mask, 'refined_seeds'] = predictions[high_conf_mask]

# Retrain with expanded seeds
scanvi_model_v2 = train_scanvi_from_seeds(scvi_model, adata, 'refined_seeds')
```

**Example 3: Combining marker-based and cluster-based seeds**
```python
# Marker-based seeds
marker_seeds = select_seed_cells(scores, n_cells=30)

# Add cluster-based seeds (from manual inspection)
cluster_seeds = adata.obs['leiden'] == '5'  # Known T cell cluster

# Combine
combined_mask = marker_seeds['T cell'] | cluster_seeds
create_seed_labels(adata, {'T cell': combined_mask}, ['T cell'])
```

## Marker Gene Signature Design

### Signature Structure

Each cell type signature should include:

| Component | Description | Example |
|-----------|-------------|---------|
| `positive` | Genes highly expressed in this type | CD3D, CD3E for T cells |
| `negative` | Genes NOT expressed in this type | CD14 for T cells |

### Best Practices for Marker Selection

1. **Use established markers**
   - Literature-validated markers
   - Cell type databases (CellMarker, PanglaoDB)
   - Previous studies on same tissue

2. **Include negative markers**
   - Critical for distinguishing similar types
   - Helps avoid false positives
   - Example: CD4 negative for CD8 T cells

3. **Balance specificity and sensitivity**
   - Too few markers: Low confidence
   - Too many markers: Over-restrictive
   - Recommended: 3-6 positive, 2-4 negative

4. **Verify marker presence**
   - Check markers exist in your gene list
   - Account for gene name conventions (symbols vs IDs)

### Example Marker Sets

**PBMC Cell Types:**
```json
{
  "CD4 T cell": {
    "positive": ["CD4", "CD3D", "CD3E"],
    "negative": ["CD8A", "CD8B", "CD14", "CD19"]
  },
  "CD8 T cell": {
    "positive": ["CD8A", "CD8B", "CD3D"],
    "negative": ["CD4", "CD14"]
  },
  "B cell": {
    "positive": ["CD19", "MS4A1", "CD79A"],
    "negative": ["CD3D", "CD14"]
  },
  "Monocyte": {
    "positive": ["CD14", "LYZ", "S100A8"],
    "negative": ["CD3D", "CD19"]
  },
  "NK cell": {
    "positive": ["NKG7", "GNLY", "NCAM1"],
    "negative": ["CD3D", "CD19", "CD14"]
  }
}
```

## Key Parameters Explained

### Seed Selection

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_seed_cells` | 50 | Cells per type for seeding |
| `min_score_percentile` | 95 | Minimum percentile cutoff |

**Guidance:**
- More seeds (100+): Better coverage, risk of noise
- Fewer seeds (20-30): Higher confidence, may miss subtypes
- Start with 50, adjust based on prediction confidence

### Model Architecture

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_latent` | 30 | Latent space dimensions |
| `n_layers` | 2 | Network depth |

### Training

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_epochs_scvi` | 100 | Epochs for base model |
| `max_epochs_scanvi` | 25 | Epochs for classifier |

**Note:** scANVI needs fewer epochs since it's initialized from scVI.

## Best Practices

1. **Marker Gene Quality**
   - Validate markers are expressed in your data
   - Use both positive and negative markers
   - Consider tissue-specific expression

2. **Seed Selection**
   - Start with conservative (fewer) seeds
   - Visualize seeds before training
   - Check for contamination across types

3. **Model Training**
   - Train scVI first for good initialization
   - Monitor training loss convergence
   - Use GPU for faster training

4. **Prediction Validation**
   - Check prediction confidence distribution
   - Validate with known markers
   - Compare to clustering results

5. **Iterative Refinement**
   - Use high-confidence predictions as new seeds
   - Retrain for better coverage
   - Manual curation for edge cases

## Reference Materials

For detailed methodology and troubleshooting, see `references/seed_labeling_guide.md`. This provides:
- Mathematical foundations of semi-supervised learning
- Detailed marker gene selection strategies
- Troubleshooting common annotation problems
- Comparison with other annotation methods

## Next Steps After Annotation

Typical downstream analysis:
- Validate annotations with marker expression
- Differential expression between cell types
- Cell type proportion analysis
- Trajectory inference within cell types
- Integration with other datasets
