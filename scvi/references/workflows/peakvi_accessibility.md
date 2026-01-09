# Single-Cell ATAC-seq Analysis with PeakVI

## Overview

PeakVI is a variational autoencoder for single-cell ATAC-seq data that:

1. **Learns a latent representation** of chromatin accessibility patterns
2. **Handles sparsity** inherent in scATAC-seq data through probabilistic modeling
3. **Enables differential accessibility** testing between cell populations
4. **Supports batch correction** through conditional modeling

**Key Advantages**:
- Fast training with early stopping
- Interpretable latent space (default 13 dimensions)
- Probabilistic differential accessibility with Bayes factors
- Memory efficient for large datasets

---

## Approach 1: Complete Analysis Pipeline (Recommended)

For standard PeakVI analysis, use the convenience script `scripts/peakvi_analysis.py`:

```bash
# Basic scATAC-seq analysis
python3 scripts/peakvi_analysis.py input.h5ad

# With batch correction and custom resolution
python3 scripts/peakvi_analysis.py input.h5ad --batch-key batch --resolution 0.3

# With custom peak filtering threshold
python3 scripts/peakvi_analysis.py input.h5ad --min-detection 0.03
```

**Parameters:**

Core parameters:
- `--batch-key` - Column in obs for batch correction
- `--min-detection` - Minimum detection rate for peak filtering (default: 0.05)

Model architecture:
- `--n-latent` - Latent dimensions (default: 13)
- `--n-hidden` - Hidden layer size (default: 128)

Training:
- `--max-epochs` - Maximum epochs (default: 500)

Downstream:
- `--resolution` - Clustering resolution (default: 0.2)
- `--n-neighbors` - Neighbors for kNN graph (default: 15)
- `--skip-da` - Skip differential accessibility

**Outputs:**
- `training_history.png` - Training loss curves
- `peakvi_umap.png` - UMAP visualization
- `marker_peaks.csv` - DA marker peaks per cluster
- `peakvi_model/` - Saved model
- `*_peakvi.h5ad` - Processed data

---

## Approach 2: Modular Building Blocks (For Custom Workflows)

For custom workflows, use functions from `scripts/peakvi_core.py`:

```python
import sys
sys.path.append('scripts/')
from peakvi_core import (
    setup_anndata_peakvi,
    train_peakvi_model,
    get_latent_representation,
    run_differential_accessibility,
    get_cluster_markers,
    filter_peaks_by_detection,
)

# Custom workflow
adata = sc.read_h5ad('input.h5ad')
adata = filter_peaks_by_detection(adata, min_detection=0.05)
setup_anndata_peakvi(adata, batch_key='batch')
model = train_peakvi_model(adata, n_latent=15, max_epochs=300)

# Get marker peaks
markers = get_cluster_markers(model, adata, 'clusters', n_top=100)
```

**Available functions:**
- `setup_anndata_peakvi()` - Register data for PeakVI
- `train_peakvi_model()` - Train with custom parameters
- `get_latent_representation()` - Extract latent space
- `run_differential_accessibility()` - DA between groups
- `get_cluster_markers()` - Find marker peaks per cluster
- `filter_peaks_by_detection()` - Filter low-detection peaks

---

## Workflow Steps (Manual)

### Step 1: Environment Setup

```python
import numpy as np
import scanpy as sc
import scvi
import matplotlib.pyplot as plt
import torch

scvi.settings.seed = 0
sc.set_figure_params(figsize=(6, 6), frameon=False)
torch.set_float32_matmul_precision("high")
```

---

### Step 2: Load scATAC-seq Data

```python
# Load preprocessed peak matrix
adata = sc.read_h5ad("path/to/atac_peaks.h5ad")

print(f"Data: {adata.n_obs} cells, {adata.n_vars} peaks")

# Check data structure
print(f"Observation columns: {list(adata.obs.columns)}")
print(f"Variable columns: {list(adata.var.columns)}")

# Verify peak coordinates if available
if 'chr' in adata.var.columns:
    print(f"Chromosomes: {adata.var['chr'].unique()[:5]}...")
```

**Data Requirements**:
- Peak-by-cell matrix (can be binary 0/1 or counts)
- Peaks as variables (columns), cells as observations (rows)
- Optional: genomic coordinates (chr, start, end) in `adata.var`
- Optional: batch labels in `adata.obs` for multi-sample data

---

### Step 3: Quality Control and Filtering

```python
# Calculate peak detection frequency
adata.var['n_cells'] = np.array((adata.X > 0).sum(axis=0)).flatten()
adata.var['detection_rate'] = adata.var['n_cells'] / adata.n_obs

# Filter peaks: keep those detected in at least 5% of cells
min_detection = 0.05
min_cells = int(adata.n_obs * min_detection)
sc.pp.filter_genes(adata, min_cells=min_cells)

print(f"Peaks after filtering (>{min_detection*100}% detection): {adata.n_vars}")

# Optional: calculate cell-level QC
adata.obs['n_peaks'] = np.array((adata.X > 0).sum(axis=1)).flatten()
adata.obs['total_counts'] = np.array(adata.X.sum(axis=1)).flatten()

# Visualize QC metrics
sc.pl.violin(adata, ['n_peaks', 'total_counts'], jitter=0.4)
```

---

### Step 4: Setup AnnData for PeakVI

```python
# Basic setup (no batch correction)
scvi.model.PEAKVI.setup_anndata(adata)

# OR: Setup with batch correction
# scvi.model.PEAKVI.setup_anndata(adata, batch_key="batch")

# Verify setup
print("AnnData registered for PeakVI")
print(f"Using layer: X")
```

**Setup Parameters**:
- `batch_key`: Column in `adata.obs` for batch labels (optional)
- Automatically registers required fields

---

### Step 5: Initialize and Train PeakVI Model

```python
# Initialize model
model = scvi.model.PEAKVI(adata)

# View model architecture
print(model)

# Train with default settings
# - max_epochs=500 (early stopping usually triggers earlier)
# - Validation-based early stopping
model.train()

# Check training history
plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.plot(model.history['elbo_train'].values, label='Train')
if 'elbo_validation' in model.history:
    plt.plot(model.history['elbo_validation'].values, label='Validation')
plt.xlabel('Epoch')
plt.ylabel('ELBO')
plt.legend()
plt.title('Training Loss')

plt.subplot(1, 2, 2)
plt.plot(model.history['reconstruction_loss_train'].values)
plt.xlabel('Epoch')
plt.ylabel('Reconstruction Loss')
plt.title('Reconstruction Loss')
plt.tight_layout()
plt.show()
```

**Training Notes**:
- Default `max_epochs=500`, but early stopping typically stops at ~50%
- Larger datasets converge faster per epoch
- Monitor reconstruction loss for convergence

---

### Step 6: Extract Latent Representation

```python
# Get latent space (default 13 dimensions)
latent = model.get_latent_representation()
adata.obsm["X_peakvi"] = latent

print(f"Latent representation shape: {latent.shape}")
```

---

### Step 7: Clustering and Visualization

```python
# Build neighborhood graph on PeakVI latent space
sc.pp.neighbors(adata, use_rep="X_peakvi")

# UMAP embedding
sc.tl.umap(adata, min_dist=0.2)

# Leiden clustering
sc.tl.leiden(adata, key_added="clusters_peakvi", resolution=0.2)

# Visualize clusters
sc.pl.umap(adata, color="clusters_peakvi", title="PeakVI Clusters")

# If batch information exists, check batch mixing
if 'batch' in adata.obs.columns:
    sc.pl.umap(adata, color=["clusters_peakvi", "batch"], ncols=2)
```

**Clustering Parameters**:
- `resolution`: Adjust for more (higher) or fewer (lower) clusters
- `min_dist`: Controls UMAP point spread (0.2 is good default)

---

### Step 8: Differential Accessibility Analysis

```python
# Method 1: Compare specific clusters
da_results = model.differential_accessibility(
    groupby="clusters_peakvi",
    group1="0",           # Target cluster
    group2="1",           # Reference cluster (or None for vs rest)
)

# Method 2: Using cell indices
# idx1 = adata.obs["clusters_peakvi"] == "0"
# idx2 = adata.obs["clusters_peakvi"] == "1"
# da_results = model.differential_accessibility(idx1=idx1, idx2=idx2)

# View results
print(da_results.head(10))
print(f"\nSignificant peaks (prob_da > 0.9): {(da_results['prob_da'] > 0.9).sum()}")

# Key columns:
# - prob_da: Probability of differential accessibility
# - is_da_fdr: Boolean, FDR-controlled significance
# - bayes_factor: Effect size (not standardized like p-values)
# - est_prob1, est_prob2: Estimated accessibility probabilities
# - emp_prob1, emp_prob2: Empirical accessibility probabilities
```

---

### Step 9: Identify Marker Peaks

```python
# Find markers for all clusters
def get_cluster_markers(model, adata, cluster_key, n_top=50):
    """Find top DA peaks for each cluster vs rest."""
    clusters = adata.obs[cluster_key].unique()
    all_markers = {}

    for cluster in clusters:
        da = model.differential_accessibility(
            groupby=cluster_key,
            group1=str(cluster),
            group2=None,  # vs rest
        )
        # Filter and sort
        markers = da[da['prob_da'] > 0.8].nlargest(n_top, 'bayes_factor')
        all_markers[cluster] = markers
        print(f"Cluster {cluster}: {len(markers)} marker peaks")

    return all_markers

markers = get_cluster_markers(model, adata, "clusters_peakvi")

# Visualize top marker peaks for a cluster
top_peaks = markers["0"].index[:6].tolist()
sc.pl.umap(adata, color=top_peaks, ncols=3, layer=None, use_raw=False)
```

---

### Step 10: Save Results

```python
# Save trained model
model_dir = "peakvi_model"
model.save(model_dir, overwrite=True)

# Save processed AnnData
adata.write_h5ad("atac_peakvi_analyzed.h5ad")

# Export marker peaks
import pandas as pd
all_markers_df = pd.concat(
    [df.assign(cluster=k) for k, df in markers.items()],
    ignore_index=False
)
all_markers_df.to_csv("marker_peaks.csv")

# Reload model later
# model = scvi.model.PEAKVI.load(model_dir, adata=adata)
```

---

## Key Parameters Reference

### Model Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_latent` | 13 | Latent space dimensions |
| `n_hidden` | 128 | Hidden layer size |
| `n_layers_encoder` | 2 | Encoder depth |
| `n_layers_decoder` | 2 | Decoder depth |
| `dropout_rate` | 0.1 | Dropout for regularization |

### Training Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_epochs` | 500 | Maximum training epochs |
| `early_stopping` | True | Stop when validation plateaus |
| `batch_size` | 128 | Samples per batch |
| `lr` | 1e-3 | Learning rate |

### Differential Accessibility Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `groupby` | None | Column for group comparison |
| `group1` | None | Target group |
| `group2` | None | Reference group (None = vs rest) |
| `batch_correction` | False | Sample from multiple batches |
| `fdr_target` | 0.05 | FDR threshold for is_da_fdr |

---

## Interpretation Guidelines

### Differential Accessibility Results

| Column | Interpretation |
|--------|----------------|
| `prob_da` | Probability peak is differentially accessible (0-1) |
| `is_da_fdr` | Conservative: True if significant after FDR |
| `bayes_factor` | Effect size magnitude (higher = stronger) |
| `est_prob1` | Model-estimated accessibility in group 1 |
| `est_prob2` | Model-estimated accessibility in group 2 |
| `emp_prob1` | Empirical (observed) accessibility in group 1 |
| `emp_prob2` | Empirical (observed) accessibility in group 2 |

### Thresholds

- **prob_da > 0.9**: High confidence differential accessibility
- **prob_da > 0.8**: Moderate confidence
- **bayes_factor > 3**: Strong evidence
- **bayes_factor > 1**: Positive evidence

---

## Advanced Usage

### Multi-Batch Analysis with Batch Correction

```python
# Setup with batch key
scvi.model.PEAKVI.setup_anndata(adata, batch_key="batch")
model = scvi.model.PEAKVI(adata)
model.train()

# Differential accessibility with batch correction
da_results = model.differential_accessibility(
    groupby="cell_type",
    group1="T_cells",
    group2="B_cells",
    batch_correction=True  # Sample from all batches
)
```

### Custom Training Configuration

```python
model = scvi.model.PEAKVI(
    adata,
    n_latent=20,           # More latent dimensions
    n_hidden=256,          # Larger hidden layers
    dropout_rate=0.2       # More regularization
)

model.train(
    max_epochs=1000,
    early_stopping=True,
    early_stopping_patience=20,
    batch_size=256
)
```

### One-vs-Rest Differential Accessibility

```python
# Find peaks specific to each cluster
for cluster in adata.obs["clusters_peakvi"].unique():
    da = model.differential_accessibility(
        groupby="clusters_peakvi",
        group1=str(cluster),
        group2=None  # Compare to all other cells
    )
    significant = da[da['is_da_fdr']].index.tolist()
    print(f"Cluster {cluster}: {len(significant)} significant peaks")
```

---

## Parameter Tuning Guide

### When to Adjust Key Parameters

**Before running analysis, ask the user about their data:**
1. How many cells do you have?
2. How sparse is your data (% of zeros)?
3. What is your biological question?
4. Have you tried clustering before? Were results too coarse or too fine?

### Peak Filtering (`min_detection`)

| Scenario | Recommended Value | Rationale |
|----------|-------------------|-----------|
| Standard 10x data (~5-10k cells) | 0.05 (5%) | Default, balances signal/noise |
| Very sparse data (<3% mean accessibility) | 0.03 (3%) | Preserve more peaks |
| High-coverage data (>10% mean accessibility) | 0.10 (10%) | Focus on robust peaks |
| Small dataset (<1000 cells) | 0.03-0.05 | Lower threshold to retain peaks |
| Large dataset (>50k cells) | 0.05-0.10 | Can afford stricter filtering |

```python
# Adaptive filtering based on data characteristics
mean_accessibility = (adata.X > 0).mean()
if mean_accessibility < 0.03:
    min_detection = 0.03  # Sparse data
elif mean_accessibility > 0.10:
    min_detection = 0.10  # Dense data
else:
    min_detection = 0.05  # Standard
```

### Clustering Resolution (`resolution`)

| Result | Problem | Adjustment |
|--------|---------|------------|
| Too few clusters | Under-clustering | Increase resolution (0.3, 0.5, 0.8) |
| Too many clusters | Over-clustering | Decrease resolution (0.1, 0.15) |
| Known cell types not separating | Need finer resolution | Increase to 0.5-1.0 |
| Rare populations missing | Resolution too low | Try 0.5+ and merge later |

```python
# Iterative resolution finding
for res in [0.1, 0.2, 0.3, 0.5, 0.8]:
    sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}")
    n_clusters = adata.obs[f"leiden_{res}"].nunique()
    print(f"Resolution {res}: {n_clusters} clusters")
# Choose resolution that gives biologically meaningful number
```

### Model Architecture (`n_latent`, `n_hidden`)

| Dataset Size | n_latent | n_hidden | Notes |
|--------------|----------|----------|-------|
| Small (<5k cells) | 10-13 | 128 | Default works well |
| Medium (5-50k cells) | 13-20 | 128-256 | May benefit from more dimensions |
| Large (>50k cells) | 15-30 | 256 | Capture more complexity |
| Very heterogeneous | 20-30 | 256 | More cell states = more dimensions |

### Training (`max_epochs`)

| Scenario | max_epochs | Early Stopping |
|----------|------------|----------------|
| Quick exploration | 200-300 | Yes |
| Standard analysis | 500 (default) | Yes |
| Final/publication | 500-1000 | Yes, but higher patience |
| Not converging | 1000+ | Check loss curves |

### Differential Accessibility Thresholds

| Analysis Goal | prob_da Threshold | Use is_da_fdr? |
|---------------|-------------------|----------------|
| Exploratory (find candidates) | 0.7 | No |
| Standard marker finding | 0.8-0.9 | Optional |
| Publication/high confidence | 0.9+ | Yes |
| Very few results | 0.7 | No, rank by bayes_factor |

---

## Adaptation Prompts for Claude

When a user invokes this skill, consider asking:

1. **Data characteristics:**
   - "How many cells and peaks are in your dataset?"
   - "What platform generated your data (10x, sci-ATAC, etc.)?"
   - "Is this single-sample or multi-batch data?"

2. **Analysis goals:**
   - "Are you looking for broad cell type clusters or fine-grained states?"
   - "Do you have expected cell types you're looking for?"
   - "Will you need to compare specific populations?"

3. **Previous attempts:**
   - "Have you tried clustering this data before? What were the results?"
   - "Are you seeing too many or too few clusters?"

4. **Computational constraints:**
   - "Do you have GPU access?"
   - "Are you working on a laptop or cluster?"

Based on answers, adjust parameters accordingly before running.

---

## Troubleshooting

### Common Issues

1. **Poor clustering (cells not separating)**:
   - Increase training epochs
   - Try higher `n_latent` (15-20)
   - Check peak filtering threshold
   - Verify data quality

2. **Training not converging**:
   - Reduce learning rate
   - Increase batch size
   - Check for constant peaks (remove)

3. **Few significant DA peaks**:
   - Lower `prob_da` threshold
   - Ensure sufficient cells per group
   - Check if groups are biologically distinct

4. **Memory issues**:
   - Reduce number of peaks (stricter filtering)
   - Use smaller batch size
   - Train on GPU if available

5. **Batch effects dominating**:
   - Use `batch_key` in setup
   - Verify batch correction in UMAP
   - Check batch sizes are balanced

---

## Integration with Other Methods

### After PeakVI Analysis

```python
# Gene activity scoring (requires peak coordinates)
# Use external tools like Signac or ArchR

# Motif enrichment on DA peaks
significant_peaks = da_results[da_results['is_da_fdr']].index
# Export to bed format for HOMER/chromVAR

# Integration with scRNA-seq
# Use scGLUE, LIGER, or WNN approaches
```

---

## References

- [PeakVI Tutorial](https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/atac/PeakVI.html)
- [scvi-tools Documentation](https://docs.scvi-tools.org/)
- [Scanpy Documentation](https://scanpy.readthedocs.io/)
