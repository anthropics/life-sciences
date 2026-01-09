# scANVI Seed Labeling Guide

Comprehensive reference for semi-supervised cell type annotation using marker genes and scANVI.

## Table of Contents

1. [Method Overview](#method-overview)
2. [When to Use Seed Labeling](#when-to-use-seed-labeling)
3. [Marker Gene Selection](#marker-gene-selection)
4. [Seed Selection Strategies](#seed-selection-strategies)
5. [Model Training](#model-training)
6. [Interpreting Results](#interpreting-results)
7. [Troubleshooting](#troubleshooting)
8. [Comparison with Other Methods](#comparison-with-other-methods)

---

## Method Overview

Seed labeling is a semi-supervised approach that combines:
1. **Marker-based cell selection**: Identify high-confidence examples using known markers
2. **Deep learning classification**: Extend labels to all cells using scANVI

### How It Works

```
1. Define marker signatures (positive/negative genes per cell type)
2. Score all cells based on marker expression
3. Select top-scoring cells as "seeds" for each type
4. Train scVI for unsupervised representation learning
5. Train scANVI using seeds as labeled examples
6. Predict labels for all unlabeled cells
```

### Key Advantages

- **Leverages domain knowledge**: Uses established marker genes
- **Scales to large datasets**: Deep learning handles millions of cells
- **Uncertainty quantification**: Provides prediction probabilities
- **Iterative refinement**: Can use predictions to expand seeds

---

## When to Use Seed Labeling

### Ideal Scenarios

1. **Known biology, no reference**
   - You know marker genes but lack a labeled reference atlas
   - Example: Well-characterized tissue with established markers

2. **Partial annotations**
   - Some cells are confidently labeled, most are not
   - Example: FACS-sorted populations mixed with unsorted

3. **Custom cell types**
   - Standard references don't include your cell types of interest
   - Example: Disease-specific cell states

4. **Quality control**
   - Validate automated annotations with known markers
   - Example: Cross-check CellTypist predictions

### When NOT to Use

- **No marker knowledge**: Use clustering + manual annotation
- **Good reference available**: Use scANVI query-to-reference
- **Very rare cell types**: May not get enough seeds
- **Highly similar types**: Markers may not discriminate well

---

## Marker Gene Selection

### Signature Structure

Each cell type needs:

```json
{
  "Cell Type Name": {
    "positive": ["GENE1", "GENE2", "GENE3"],
    "negative": ["GENE4", "GENE5"]
  }
}
```

### Positive Markers

Genes that are:
- **Highly expressed** in the target cell type
- **Specific** (not expressed broadly)
- **Consistent** across conditions/batches

**Good examples:**
- CD3D, CD3E for T cells (T cell receptor)
- CD19, MS4A1 for B cells (B cell markers)
- CD14, LYZ for monocytes (myeloid markers)

### Negative Markers

Genes that are:
- **NOT expressed** in the target type
- **Expressed** in similar/confounding types
- Critical for **distinguishing** similar populations

**Good examples:**
- CD4 negative for CD8 T cells
- CD3D negative for B cells
- CD14 negative for T cells

### How Many Markers?

| Marker Count | Trade-off |
|--------------|-----------|
| 1-2 | Too permissive, may select wrong cells |
| 3-5 | Good balance (recommended) |
| 6-10 | More specific, may miss cells |
| >10 | Over-restrictive, few seeds |

### Marker Resources

1. **CellMarker** (http://xteam.xbio.top/CellMarker/)
2. **PanglaoDB** (https://panglaodb.se/)
3. **The Human Protein Atlas** (https://www.proteinatlas.org/)
4. **Literature review** for your specific tissue

### Validating Markers

Before running seed labeling:

```python
# Check marker expression in your data
import scanpy as sc

# Normalize for visualization
adata_norm = adata.copy()
sc.pp.normalize_total(adata_norm, target_sum=1e4)
sc.pp.log1p(adata_norm)

# Plot marker expression
markers = ['CD3D', 'CD14', 'CD19', 'NKG7']
sc.pl.violin(adata_norm, markers, groupby='leiden')
sc.pl.umap(adata_norm, color=markers)
```

---

## Seed Selection Strategies

### Strategy 1: Top N by Score (Default)

Select the N highest-scoring cells per type.

```python
seed_masks = select_seed_cells(scores, n_cells=50, method='top_n')
```

**Pros**: Consistent number of seeds
**Cons**: May include low-confidence cells if N is too high

### Strategy 2: Percentile Threshold

Select cells above a score percentile.

```python
seed_masks = select_seed_cells(scores, min_percentile=95, method='percentile')
```

**Pros**: Only high-confidence cells
**Cons**: Variable number of seeds per type

### Strategy 3: Combined Approach

Use top N with minimum percentile filter.

```python
seed_masks = select_seed_cells(
    scores,
    n_cells=50,
    min_percentile=90,
    method='top_n'
)
```

### How Many Seeds?

| Dataset Size | Recommended Seeds/Type |
|--------------|----------------------|
| <10,000 cells | 20-30 |
| 10,000-50,000 | 50 |
| 50,000-200,000 | 50-100 |
| >200,000 | 100-200 |

### Visualizing Seeds

Always check seeds before training:

```python
# Color UMAP by seed labels
sc.pl.umap(adata, color='seed_labels')

# Check for overlapping seeds (should be minimal)
seed_counts = adata.obs['seed_labels'].value_counts()
print(seed_counts)
```

---

## Model Training

### scVI Training

scVI learns a latent representation that captures biological variation:

```python
scvi.model.SCVI.setup_anndata(adata, labels_key='seed_labels')
scvi_model = scvi.model.SCVI(adata, n_latent=30, n_layers=2)
scvi_model.train(max_epochs=100)
```

**Key parameters:**
- `n_latent`: 30 works well for most datasets
- `n_layers`: 2 is standard
- `max_epochs`: 100-200 depending on size

### scANVI Training

scANVI adds a classifier on top of scVI:

```python
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    unlabeled_category='Unknown'
)
scanvi_model.train(max_epochs=25)
```

**Key points:**
- Initialize from trained scVI (important!)
- Fewer epochs needed (already has good representation)
- `unlabeled_category` must match your seed labels

### Training Monitoring

```python
# Check training loss
import matplotlib.pyplot as plt

plt.plot(scvi_model.history['elbo_train'])
plt.xlabel('Epoch')
plt.ylabel('ELBO')
plt.title('scVI Training Loss')
```

---

## Interpreting Results

### Prediction Confidence

```python
# Get probabilities
probs = scanvi_model.predict(soft=True)

# Maximum probability = confidence
confidence = probs.max(axis=1)

# Distribution
plt.hist(confidence, bins=50)
plt.xlabel('Confidence')
plt.ylabel('Cells')
```

**Interpretation:**
- >0.9: High confidence, trust prediction
- 0.7-0.9: Moderate confidence, likely correct
- 0.5-0.7: Low confidence, may need review
- <0.5: Very uncertain, manual curation needed

### Confidence by Cell Type

```python
# Per-type confidence
adata.obs['confidence'] = confidence
adata.obs['prediction'] = predictions

for ct in adata.obs['prediction'].unique():
    mask = adata.obs['prediction'] == ct
    mean_conf = adata.obs.loc[mask, 'confidence'].mean()
    print(f"{ct}: {mean_conf:.3f}")
```

### Validation with Markers

```python
# Check that predictions match marker expression
sc.pl.dotplot(
    adata,
    var_names=['CD3D', 'CD14', 'CD19'],
    groupby='prediction'
)
```

### Confusion Analysis

Compare seeds to predictions:

```python
# Confusion between seed and predicted labels
from sklearn.metrics import confusion_matrix

# Only for seeded cells
seeded = adata.obs['seed_labels'] != 'Unknown'
cm = confusion_matrix(
    adata.obs.loc[seeded, 'seed_labels'],
    adata.obs.loc[seeded, 'prediction']
)
```

---

## Troubleshooting

### Problem: Wrong Cell Types Predicted

**Symptoms:**
- Known cell types mislabeled
- Markers don't match predictions

**Solutions:**
1. Improve marker signatures (add negatives)
2. Reduce number of seeds (higher quality)
3. Check for batch effects
4. Visualize seeds before training

### Problem: Low Overall Confidence

**Symptoms:**
- Most predictions <0.7 confidence
- Uncertainty across all types

**Solutions:**
1. Add more discriminative markers
2. Increase number of seeds
3. Train scVI longer
4. Check data quality

### Problem: One Type Dominates

**Symptoms:**
- Most cells predicted as one type
- Other types have very few predictions

**Solutions:**
1. Balance seed numbers across types
2. Add negative markers for dominant type
3. Check marker specificity
4. May indicate true biology

### Problem: Seeds Overlap

**Symptoms:**
- Same cells selected for multiple types
- Conflicting seed labels

**Solutions:**
1. Add distinguishing negative markers
2. Use more specific positive markers
3. Reduce seed count
4. Accept that types may be similar

### Problem: Missing Cell Types

**Symptoms:**
- Known type not in predictions
- Type has no seeds

**Solutions:**
1. Verify markers are present in data
2. Check gene name format (symbols vs IDs)
3. Lower percentile threshold
4. Type may not exist in this dataset

---

## Comparison with Other Methods

### Seed Labeling vs CellTypist

| Aspect | Seed Labeling | CellTypist |
|--------|---------------|------------|
| Markers | User-defined | Pre-trained |
| Flexibility | High | Limited to models |
| Setup | Requires markers | Ready to use |
| Cell types | Any | Model-specific |

**Use CellTypist when:** Standard cell types, quick annotation
**Use Seed Labeling when:** Custom types, specific markers

### Seed Labeling vs Query-to-Reference

| Aspect | Seed Labeling | Query-to-Reference |
|--------|---------------|-------------------|
| Reference | Not needed | Required |
| Labels | From markers | From reference |
| Batch effects | Handled | Must match reference |

**Use Query-to-Reference when:** Good reference atlas exists
**Use Seed Labeling when:** No suitable reference

### Seed Labeling vs Clustering + Manual

| Aspect | Seed Labeling | Manual |
|--------|---------------|--------|
| Automation | High | Low |
| Consistency | High | Variable |
| Time | Fast | Slow |
| Flexibility | Moderate | High |

**Use Manual when:** Novel biology, uncertain markers
**Use Seed Labeling when:** Known markers, large datasets

---

## Best Practices Summary

1. **Marker quality over quantity** - Better to have 3 excellent markers than 10 mediocre ones

2. **Always include negatives** - Critical for distinguishing similar types

3. **Visualize seeds** - Never train without checking seed quality

4. **Start conservative** - Begin with fewer seeds, add more if needed

5. **Validate predictions** - Check marker expression matches predictions

6. **Iterate if needed** - Use high-confidence predictions as new seeds

7. **Document everything** - Record markers and parameters for reproducibility
