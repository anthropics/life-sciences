# Multimodal Integration with MultiVI

## Overview

MultiVI enables joint analysis of datasets with heterogeneous modality profiles:

1. **Paired cell anchoring**: Uses fully-paired cells to learn cross-modal relationships
2. **Unpaired integration**: Maps single-modality cells using learned relationships
3. **Modality imputation**: Predicts missing RNA or ATAC for single-modality cells
4. **Unified embedding**: Single latent space regardless of measured modalities

**Typical Scenario**:
- 1/3 cells: RNA-only
- 1/3 cells: Paired RNA+ATAC (multiome)
- 1/3 cells: ATAC-only

**Key Insight**: "MultiVI requires the features to be ordered so that genes appear before genomic regions"

**Citation**: Ashuach et al. (2023). "MultiVI: deep generative model for the integration of multimodal data." Nature Methods.

---

## Workflow Steps

### Step 1: Environment Setup

```python
import numpy as np
import scanpy as sc
import scvi
import muon
import matplotlib.pyplot as plt
import torch

scvi.settings.seed = 0
sc.set_figure_params(figsize=(6, 6), frameon=False)
torch.set_float32_matmul_precision("high")
```

---

### Step 2: Load and Understand Your Data

```python
# Load MuData with RNA and ATAC modalities
mdata = muon.read_h5mu("path/to/multiome_data.h5mu")

# Check what modalities exist
print("Modalities:", list(mdata.mod.keys()))

# Access individual modalities
rna = mdata.mod["rna"]
atac = mdata.mod["atac"]

print(f"RNA: {rna.n_obs} cells, {rna.n_vars} genes")
print(f"ATAC: {atac.n_obs} cells, {atac.n_vars} peaks")

# Check for paired vs unpaired cells
# Cells may have different modality profiles
if "modality" in mdata.obs.columns:
    print("\nModality profiles:")
    print(mdata.obs["modality"].value_counts())
```

---

### Step 3: Identify Paired and Unpaired Cells

```python
# Determine which cells have which modalities
# Based on presence in each modality's AnnData

rna_cells = set(rna.obs_names)
atac_cells = set(atac.obs_names)

paired_cells = rna_cells.intersection(atac_cells)
rna_only_cells = rna_cells - atac_cells
atac_only_cells = atac_cells - rna_cells

print(f"Paired (RNA+ATAC): {len(paired_cells)}")
print(f"RNA-only: {len(rna_only_cells)}")
print(f"ATAC-only: {len(atac_only_cells)}")

# Add modality labels
mdata.obs["modality_profile"] = "unknown"
mdata.obs.loc[list(paired_cells), "modality_profile"] = "paired"
mdata.obs.loc[list(rna_only_cells), "modality_profile"] = "rna_only"
mdata.obs.loc[list(atac_only_cells), "modality_profile"] = "atac_only"
```

---

### Step 4: Prepare Features (CRITICAL: Genes Before Peaks)

```python
# CRITICAL: MultiVI requires genes ordered before peaks
# This is handled automatically with MuData setup, but verify

# Ensure raw counts
if "counts" not in rna.layers:
    rna.layers["counts"] = rna.X.copy()
if "counts" not in atac.layers:
    atac.layers["counts"] = atac.X.copy()

# HVG selection for RNA
sc.pp.highly_variable_genes(
    rna,
    n_top_genes=2000,
    flavor="seurat_v3",
    layer="counts",
    subset=True
)

# Peak filtering for ATAC (keep peaks in >5% of cells)
peak_counts = np.array((atac.X > 0).sum(axis=0)).flatten()
min_cells = int(atac.n_obs * 0.05)
atac = atac[:, peak_counts >= min_cells].copy()

print(f"Selected genes: {rna.n_vars}")
print(f"Selected peaks: {atac.n_vars}")

# Update MuData
mdata.mod["rna"] = rna
mdata.mod["atac"] = atac

# Convert to CSR format for faster training
rna.X = rna.X.tocsr() if hasattr(rna.X, 'tocsr') else rna.X
atac.X = atac.X.tocsr() if hasattr(atac.X, 'tocsr') else atac.X
```

---

### Step 5: Setup MultiVI

```python
# Setup MuData for MultiVI
scvi.model.MULTIVI.setup_mudata(
    mdata,
    modalities={
        "rna_layer": "rna",
        "atac_layer": "atac",
    }
)

# Note: The main batch annotation should correspond to modality type
# Additional batch effects can be specified with categorical_covariate_keys

print("MuData registered for MultiVI")
```

---

### Step 6: Initialize and Train MultiVI

```python
# Initialize model
model = scvi.model.MULTIVI(
    mdata,
    n_genes=rna.n_vars,
    n_regions=atac.n_vars
)

# View model summary
print(model)

# Train
model.train(max_epochs=500)

# Plot training history
plt.figure(figsize=(8, 4))
plt.plot(model.history['elbo_train'].values, label='Train')
if 'elbo_validation' in model.history:
    plt.plot(model.history['elbo_validation'].values, label='Validation')
plt.xlabel('Epoch')
plt.ylabel('ELBO')
plt.legend()
plt.title('MultiVI Training')
plt.show()
```

---

### Step 7: Extract Joint Embedding

```python
# Get unified latent representation (all cells)
latent = model.get_latent_representation()
mdata.obsm["X_multivi"] = latent

# Build neighborhood graph
sc.pp.neighbors(mdata, use_rep="X_multivi")

# UMAP
sc.tl.umap(mdata, min_dist=0.2)

# Clustering
sc.tl.leiden(mdata, key_added="leiden_multivi", resolution=0.5)

print(f"Joint embedding shape: {latent.shape}")
```

---

### Step 8: Impute Missing Modalities

```python
# Impute gene expression for ATAC-only cells
rna_imputed = model.get_normalized_expression()

# Impute accessibility for RNA-only cells
atac_imputed = model.get_accessibility_estimates()

# Store imputed values
mdata.obsm["rna_imputed"] = rna_imputed
mdata.obsm["atac_imputed"] = atac_imputed

# For ATAC-only cells, check imputed gene expression
atac_only_mask = mdata.obs["modality_profile"] == "atac_only"
print(f"\nImputed RNA for {atac_only_mask.sum()} ATAC-only cells")

# Visualize imputed marker genes
marker_genes = ["CD3G", "MS4A1", "NCAM1"]  # Example markers
for gene in marker_genes:
    if gene in rna.var_names:
        gene_idx = list(rna.var_names).index(gene)
        mdata.obs[f"imputed_{gene}"] = rna_imputed[:, gene_idx]
```

---

### Step 9: Visualize Results

```python
# Comprehensive visualization
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# Row 1: Integration quality
sc.pl.umap(mdata, color="modality_profile", ax=axes[0, 0], show=False,
           title="Modality Profile")
sc.pl.umap(mdata, color="leiden_multivi", ax=axes[0, 1], show=False,
           title="Clusters")

# Sample batch if present
if "batch" in mdata.obs.columns:
    sc.pl.umap(mdata, color="batch", ax=axes[0, 2], show=False,
               title="Sample Batch")
else:
    axes[0, 2].set_visible(False)

# Row 2: Imputed expression
for i, gene in enumerate(marker_genes[:3]):
    if f"imputed_{gene}" in mdata.obs.columns:
        sc.pl.umap(mdata, color=f"imputed_{gene}", ax=axes[1, i],
                   show=False, title=f"{gene} (imputed)", cmap="viridis")

plt.tight_layout()
plt.show()

# Check that imputation is consistent across modality profiles
print("\nImputed expression by modality profile:")
for gene in marker_genes[:3]:
    if f"imputed_{gene}" in mdata.obs.columns:
        for profile in ["paired", "rna_only", "atac_only"]:
            mask = mdata.obs["modality_profile"] == profile
            if mask.sum() > 0:
                mean_expr = mdata.obs.loc[mask, f"imputed_{gene}"].mean()
                print(f"  {gene} in {profile}: {mean_expr:.3f}")
```

---

### Step 10: Save Results

```python
# Save trained model
model_dir = "multivi_model"
model.save(model_dir, overwrite=True)

# Save processed MuData
mdata.write_h5mu("multiome_multivi_analyzed.h5mu")

# Export imputed values
import pandas as pd

# Imputed RNA
rna_imputed_df = pd.DataFrame(
    rna_imputed,
    index=mdata.obs_names,
    columns=rna.var_names
)
rna_imputed_df.to_csv("imputed_rna.csv")

# Reload model later
# model = scvi.model.MULTIVI.load(model_dir, adata=mdata)

print("Results saved")
```

---

## Key Parameters Reference

### Setup Parameters

| Parameter | Description |
|-----------|-------------|
| `modalities` | Dict mapping modality names |
| `categorical_covariate_keys` | Additional batch covariates |

### Model Parameters

| Parameter | Description |
|-----------|-------------|
| `n_genes` | Number of genes (from RNA modality) |
| `n_regions` | Number of peaks (from ATAC modality) |
| `n_latent` | Latent dimensions (default 10) |

### Training Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_epochs` | 500 | Maximum training |
| `batch_size` | 128 | Samples per batch |

---

## Parameter Tuning Guide

### When to Adjust Parameters

**Before running, ask the user:**
1. What fraction of cells are paired vs unpaired?
2. How many peaks after filtering?
3. Are there additional batch effects beyond modality?
4. What is the main goal (embedding, imputation, both)?

### Impact of Paired Cell Ratio

| Paired Ratio | Quality | Notes |
|--------------|---------|-------|
| >50% | Excellent | Strong anchoring |
| 25-50% | Good | Usually sufficient |
| 10-25% | Fair | May need more epochs |
| <10% | Challenging | Consider alternatives |

### Latent Dimensions

| Complexity | n_latent | Notes |
|------------|----------|-------|
| Simple cell types | 8-10 | Default |
| Moderate complexity | 12-15 | Standard |
| High heterogeneity | 15-20 | More capacity |

### Peak Filtering

| Scenario | min_detection | Notes |
|----------|---------------|-------|
| Standard | 0.05 (5%) | Default |
| Sparse ATAC | 0.03 | Preserve more peaks |
| Large dataset | 0.08-0.10 | Focus on robust peaks |

---

## Adaptation Prompts for Claude

When a user invokes this skill, consider asking:

1. **Data composition:**
   - "What fraction of cells have both modalities (paired)?"
   - "How were the different datasets generated?"
   - "Do paired and unpaired cells come from same tissue?"

2. **Technical details:**
   - "What ATAC peak caller was used?"
   - "Are peaks shared across all ATAC samples?"
   - "Any known batch effects beyond modality?"

3. **Analysis goals:**
   - "Is the embedding or imputation more important?"
   - "What downstream analyses are planned?"
   - "Do you need to impute RNA, ATAC, or both?"

### Data Validation Helper

```python
def validate_multivi_data(mdata):
    """Validate data for MultiVI analysis."""
    issues = []
    recommendations = []

    # Check modalities
    if "rna" not in mdata.mod:
        issues.append("Missing 'rna' modality")
    if "atac" not in mdata.mod:
        issues.append("Missing 'atac' modality")

    if issues:
        print("ERRORS:")
        for i in issues:
            print(f"  {i}")
        return False

    rna = mdata.mod["rna"]
    atac = mdata.mod["atac"]

    # Check for paired cells
    paired = set(rna.obs_names).intersection(set(atac.obs_names))
    paired_ratio = len(paired) / len(mdata.obs_names)

    if paired_ratio < 0.1:
        recommendations.append(f"Low paired ratio ({paired_ratio:.1%}). "
                               "Results may be less reliable.")
    elif paired_ratio < 0.25:
        recommendations.append(f"Moderate paired ratio ({paired_ratio:.1%}). "
                               "Consider longer training.")
    else:
        print(f"Good paired ratio: {paired_ratio:.1%}")

    # Check feature counts
    if atac.n_vars > 100000:
        recommendations.append(f"Many peaks ({atac.n_vars}). "
                               "Consider stricter filtering.")

    if rna.n_vars > 5000:
        recommendations.append(f"Many genes ({rna.n_vars}). "
                               "Consider using 2000-3000 HVGs.")

    print("\nRecommendations:")
    for rec in recommendations:
        print(f"  - {rec}")

    return True
```

---

## Troubleshooting

### Common Issues

1. **"Features must be ordered: genes before peaks"**:
   - MultiVI expects concatenated features with genes first
   - MuData setup usually handles this, but verify structure

2. **Poor integration (modalities separate in UMAP)**:
   - Check paired cell ratio
   - Verify cells share biological variation
   - Try more training epochs

3. **Imputation looks uniform**:
   - May indicate poor cross-modal learning
   - Check if paired cells span all cell types
   - Verify feature selection

4. **Memory issues**:
   - Filter more peaks
   - Use fewer HVGs
   - Convert to CSR sparse format

5. **Shared peak requirement**:
   - All ATAC datasets must use same peak set
   - Use ArchR, SnapATAC, or CellRanger merge

### Checking Modality Balance

```python
def check_modality_balance(mdata):
    """Check cell type distribution across modality profiles."""
    import pandas as pd

    if "leiden_multivi" in mdata.obs.columns and "modality_profile" in mdata.obs.columns:
        ct = pd.crosstab(mdata.obs["leiden_multivi"],
                         mdata.obs["modality_profile"])
        print("Cells per cluster by modality profile:")
        print(ct)

        # Check for modality-specific clusters
        ct_norm = ct.div(ct.sum(axis=1), axis=0)
        for profile in ["paired", "rna_only", "atac_only"]:
            dominated = (ct_norm[profile] > 0.9).sum() if profile in ct_norm.columns else 0
            print(f"\nClusters dominated by {profile}: {dominated}")
```

---

## Advanced Usage

### Additional Batch Correction

```python
# If you have batch effects beyond modality
mdata.obs["sample_batch"] = ...  # Your sample batch labels

scvi.model.MULTIVI.setup_mudata(
    mdata,
    modalities={"rna_layer": "rna", "atac_layer": "atac"},
    categorical_covariate_keys=["sample_batch"]  # Additional batches
)

model = scvi.model.MULTIVI(mdata, n_genes=rna.n_vars, n_regions=atac.n_vars)
model.train()
```

### Cross-Modal Correlation Analysis

```python
# After training, examine RNA-ATAC relationships
def analyze_cross_modal_correlation(model, mdata, gene, peaks_near_gene):
    """Analyze correlation between gene expression and nearby peaks."""

    rna_imputed = model.get_normalized_expression()
    atac_imputed = model.get_accessibility_estimates()

    gene_idx = list(mdata.mod["rna"].var_names).index(gene)
    gene_expr = rna_imputed[:, gene_idx]

    print(f"\nCorrelation of {gene} with nearby peaks:")
    for peak in peaks_near_gene:
        if peak in mdata.mod["atac"].var_names:
            peak_idx = list(mdata.mod["atac"].var_names).index(peak)
            peak_acc = atac_imputed[:, peak_idx]
            corr = np.corrcoef(gene_expr, peak_acc)[0, 1]
            print(f"  {peak}: r = {corr:.3f}")
```

---

## References

- [MultiVI Tutorial](https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/multimodal/MultiVI_tutorial.html)
- [MultiVI Paper](https://www.nature.com/articles/s41592-023-01909-9) - Ashuach et al., Nature Methods 2023
- [scvi-tools Documentation](https://docs.scvi-tools.org/)
- [MuData Documentation](https://mudata.readthedocs.io/)
