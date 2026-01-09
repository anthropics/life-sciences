# scVI/scANVI Label Transfer Guide

Comprehensive reference for transferring cell type labels from reference to query datasets.

## Table of Contents

1. [Method Overview](#method-overview)
2. [When to Use Label Transfer](#when-to-use-label-transfer)
3. [Cross-Technology Integration](#cross-technology-integration)
4. [Data Preparation](#data-preparation)
5. [Training Strategy](#training-strategy)
6. [Evaluating Transfer Quality](#evaluating-transfer-quality)
7. [Troubleshooting](#troubleshooting)
8. [Comparison with Other Methods](#comparison-with-other-methods)

---

## Method Overview

Label transfer uses scVI for integration and scANVI for semi-supervised classification:

### Two-Stage Approach

```
Stage 1: scVI Integration
- Learn joint latent space for reference + query
- Correct for batch/technology effects
- No labels used

Stage 2: scANVI Transfer
- Initialize from scVI
- Train classifier using reference labels
- Predict labels for query cells
```

### Why This Works

1. **scVI integration** creates a shared embedding where similar cells cluster together regardless of batch
2. **scANVI** leverages this embedding plus reference labels to learn a cell type classifier
3. Query cells are classified based on their position in the integrated space

---

## When to Use Label Transfer

### Ideal Scenarios

| Scenario | Recommendation |
|----------|---------------|
| Annotated reference + unannotated query | **Use this method** |
| Cross-technology (10x + SmartSeq2) | **Use this method** |
| Large reference, small query | **Use this method** |
| Novel cell types expected in query | Consider alternatives |
| Pre-trained reference model available | Use scArches instead |

### Label Transfer vs Alternatives

**vs scArches (reference mapping)**
- Label transfer: Retrain from scratch on combined data
- scArches: Use pre-trained model, fine-tune on query
- Use scArches when reference model is available

**vs Seed labeling**
- Label transfer: Uses labeled reference cells
- Seed labeling: Uses marker genes only
- Use seed labeling when no reference exists

**vs CellTypist**
- Label transfer: Custom reference, deep learning
- CellTypist: Pre-built models, logistic regression
- Use CellTypist for common tissue types

---

## Cross-Technology Integration

### Technology Comparison

| Technology | Method | UMI | Gene Length Bias |
|------------|--------|-----|------------------|
| 10x Genomics | Droplet | Yes | No |
| SmartSeq2 | Plate | No | Yes |
| Drop-seq | Droplet | Yes | No |
| CEL-seq2 | Plate | Yes | No |
| inDrop | Droplet | Yes | No |

### Gene Length Normalization

**Why it's needed for SmartSeq2:**
- SmartSeq2 sequences full-length transcripts
- Longer genes have more reads (proportional to length)
- UMI methods count molecules, not reads

**Normalization formula:**
```python
normalized_counts = raw_counts / gene_length * median(gene_lengths)
```

**When to apply:**
- Reference is SmartSeq2: Normalize reference
- Query is SmartSeq2: Normalize query
- Both are SmartSeq2: Normalize both
- Both are UMI-based: No normalization needed

### Batch Key Strategy

```python
# Option 1: Technology as batch (simplest)
adata.obs['batch'] = adata.obs['technology']

# Option 2: Dataset as batch
adata.obs['batch'] = ['ref'] * n_ref + ['query'] * n_query

# Option 3: Technology + sample (most granular)
adata.obs['batch'] = adata.obs['technology'] + '_' + adata.obs['sample']
```

**Recommendation:** Start with technology as batch. Add sample if integration is poor.

---

## Data Preparation

### Gene Matching

Ensure gene names match between datasets:

```python
# Check gene name overlap
ref_genes = set(reference.var_names)
query_genes = set(query.var_names)
common = ref_genes & query_genes

print(f"Reference genes: {len(ref_genes)}")
print(f"Query genes: {len(query_genes)}")
print(f"Common genes: {len(common)}")

# Subset to common genes
reference = reference[:, list(common)]
query = query[:, list(common)]
```

**Common issues:**
- Different gene name formats (symbols vs Ensembl IDs)
- Species-specific naming (GAPDH vs Gapdh)
- Version suffixes (ENSG00000111640.15)

### Count Data Requirements

scVI/scANVI require raw counts:

```python
# Check for raw counts
if hasattr(reference.X, 'toarray'):
    sample = reference.X[:100].toarray()
else:
    sample = reference.X[:100]

# Should be integers
is_counts = np.allclose(sample, sample.astype(int))
print(f"Data appears to be counts: {is_counts}")

# If normalized, look for raw counts
if 'counts' in reference.layers:
    print("Using counts from layers['counts']")
```

### HVG Selection

**Why batch-aware selection:**
- Prevents batch-specific genes from dominating
- Ensures shared biology is captured

```python
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    batch_key='batch',  # Important!
    flavor='seurat_v3',
    layer='counts'
)
```

**How many genes?**
- 2000: Standard, works well for most cases
- 3000-4000: Complex tissues, many cell types
- 1000: Simple datasets, faster training

---

## Training Strategy

### scVI Training

**Purpose:** Create integrated latent space

**Key parameters:**
```python
scvi.model.SCVI(
    adata,
    n_latent=30,      # Embedding dimension
    n_layers=2,        # Network depth
    gene_likelihood='nb'  # Negative binomial
)

model.train(
    max_epochs=400,    # Sufficient for convergence
    early_stopping=True
)
```

**Monitoring:**
- Watch ELBO loss decrease
- Early stopping prevents overfitting
- Check batch mixing in latent space

### scANVI Training

**Purpose:** Learn classifier using reference labels

**Key parameters:**
```python
scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    unlabeled_category='Unknown',  # Query cell label
    labels_key='cell_type'
)

model.train(
    max_epochs=20,     # Fewer epochs (initialized from scVI)
    n_samples_per_label=100  # Balance sampling
)
```

**Why fewer epochs?**
- scANVI inherits scVI's learned representation
- Only needs to train classifier head
- Too many epochs may overfit to reference labels

---

## Evaluating Transfer Quality

### Metrics

**Accuracy (for reference cells):**
```python
ref_true = adata.obs.loc[ref_mask, 'true_label']
ref_pred = adata.obs.loc[ref_mask, 'predicted']
accuracy = (ref_true == ref_pred).mean()
```

**F1 Score:**
```python
from sklearn.metrics import f1_score
f1 = f1_score(ref_true, ref_pred, average='weighted')
```

**Confusion Matrix:**
```python
cm = pd.crosstab(ref_true, ref_pred, normalize='index')
```

### Confidence Assessment

```python
# Get probabilities
probs = model.predict(soft=True)
confidence = probs.max(axis=1)

# Distribution
print(f"Mean confidence: {confidence.mean():.3f}")
print(f"Cells >0.9: {(confidence > 0.9).sum()}")
print(f"Cells <0.5: {(confidence < 0.5).sum()}")
```

**Interpretation:**
| Confidence | Interpretation |
|------------|---------------|
| >0.9 | High confidence, trust prediction |
| 0.7-0.9 | Moderate, likely correct |
| 0.5-0.7 | Low, may need review |
| <0.5 | Very uncertain, possible novel type |

### Visual Validation

1. **UMAP colored by predictions** - Do clusters make sense?
2. **UMAP colored by confidence** - Where are uncertain cells?
3. **Marker expression** - Do markers match predictions?

---

## Troubleshooting

### Problem: Poor Integration

**Symptoms:**
- Reference and query don't mix in UMAP
- Clear separation by dataset

**Solutions:**
1. Check batch_key is correct
2. Increase training epochs
3. Add technology as batch variable
4. Check for major biological differences
5. Verify gene names match

### Problem: All Cells Predicted as One Type

**Symptoms:**
- Single cell type dominates predictions
- Other types have near-zero cells

**Solutions:**
1. Check label distribution in reference
2. Increase n_samples_per_label
3. Verify labels_key is correct
4. Check for class imbalance

### Problem: Low Confidence Everywhere

**Symptoms:**
- Most cells <0.7 confidence
- Uncertain across all types

**Solutions:**
1. Train scVI longer
2. Increase n_latent
3. Add more HVGs
4. Check data quality
5. May indicate biological mismatch

### Problem: Wrong Predictions

**Symptoms:**
- Known cells misclassified
- Confusion between specific types

**Solutions:**
1. Check reference annotation quality
2. Similar types may be indistinguishable
3. Increase training data
4. Consider merging similar types

### Problem: Query Has Novel Types

**Symptoms:**
- Low confidence for many query cells
- Predictions don't match expected biology

**Solutions:**
1. This is expected behavior!
2. Filter by confidence
3. Label low-confidence as "Unknown"
4. Cluster and annotate manually

---

## Comparison with Other Methods

### Label Transfer vs Ingest (Scanpy)

| Aspect | scVI/scANVI | scanpy.tl.ingest |
|--------|-------------|------------------|
| Integration | Deep learning | kNN projection |
| Batch correction | Built-in | Limited |
| Scalability | Excellent | Good |
| Training time | Longer | Fast |
| Accuracy | Generally higher | Good for simple cases |

### Label Transfer vs Seurat v3/v4

| Aspect | scVI/scANVI | Seurat |
|--------|-------------|--------|
| Method | VAE | Anchors/CCA |
| Probabilistic | Yes | No |
| Confidence | Built-in | Requires additional |
| Cross-technology | Excellent | Good |

### Label Transfer vs Symphony

| Aspect | scVI/scANVI | Symphony |
|--------|-------------|----------|
| Method | VAE | Linear |
| Speed | Slower | Fast |
| Reference update | Retrain | Not needed |
| Memory | Higher | Lower |

---

## Best Practices Summary

1. **Always use raw counts** - Never normalize before scVI/scANVI

2. **Match gene names** - Ensure consistent naming between datasets

3. **Apply gene length normalization** - Required for SmartSeq2 data

4. **Use batch-aware HVG selection** - Prevents batch-specific gene bias

5. **Train scVI first** - Provides good initialization for scANVI

6. **Evaluate on reference** - Use known labels to assess accuracy

7. **Filter by confidence** - Don't trust low-confidence predictions

8. **Validate with markers** - Check predictions against known biology

9. **Consider novel types** - Query may have cell types not in reference

10. **Save models** - Enable reproducibility and reuse
