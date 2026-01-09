# CITE-Seq Analysis with TotalVI

## Overview

TotalVI (Total Variational Inference) is a probabilistic model for joint analysis of CITE-Seq data that:

1. **Joint modeling**: Learns a shared latent space from RNA and protein
2. **Denoising**: Separates true signal from background noise in protein measurements
3. **Batch correction**: Handles technical variation across samples/batches
4. **Differential expression**: Tests for DE in both RNA and protein simultaneously

**Key Advantages**:
- Accounts for protein-specific background noise
- Provides foreground probability for protein detection
- Single model for both modalities
- Unified differential expression framework

**Citation**: Gayoso et al. (2021). "Joint probabilistic modeling of single-cell multi-omic data with totalVI." Nature Methods.

---

## Approach 1: Complete Analysis Pipeline (Recommended)

For standard TotalVI analysis, use the convenience script `scripts/totalvi_analysis.py`:

```bash
# Basic CITE-Seq analysis with MuData
python3 scripts/totalvi_analysis.py input.h5mu --batch-key batch

# With AnnData format (protein in obsm)
python3 scripts/totalvi_analysis.py input.h5ad --protein-obsm-key protein_expression

# With differential expression between specific clusters
python3 scripts/totalvi_analysis.py input.h5mu --de-clusters 0 1 --de-delta 0.5
```

**Parameters:**

Core parameters:
- `--batch-key` - Column in obs for batch correction
- `--protein-obsm-key` - For AnnData: obsm key with protein data (default: protein_expression)
- `--rna-mod-key` / `--protein-mod-key` - For MuData: modality names (default: rna/prot)

Model architecture:
- `--n-latent` - Latent dimensions (default: 20)
- `--n-layers-encoder` - Encoder layers (default: 2)

Training:
- `--max-epochs` - Maximum epochs (default: 400)
- `--early-stopping` / `--no-early-stopping`

Downstream:
- `--resolution` - Clustering resolution (default: 0.5)
- `--de-clusters` - Two cluster IDs to compare for DE
- `--skip-de` - Skip differential expression

**Outputs:**
- `training_history.png` - Training loss curves
- `totalvi_umap.png` - UMAP visualization
- `differential_expression.csv` - DE results (if not skipped)
- `totalvi_model/` - Saved model
- `*_totalvi.h5mu/h5ad` - Processed data

---

## Approach 2: Modular Building Blocks (For Custom Workflows)

For custom workflows, use functions from `scripts/totalvi_core.py`:

```python
import sys
sys.path.append('scripts/')
from totalvi_core import (
    setup_mudata_totalvi,
    setup_anndata_totalvi,
    train_totalvi_model,
    get_latent_representation,
    get_denoised_expression,
    get_foreground_probability,
    run_differential_expression,
    validate_citeseq_data,
)

# Custom workflow example
mdata = muon.read_h5mu('input.h5mu')
validate_citeseq_data(mdata, 'rna', 'prot')
setup_mudata_totalvi(mdata, rna_layer='counts', batch_key='batch')
model = train_totalvi_model(mdata, n_latent=25, max_epochs=500)

# Get denoised values
rna_denoised, protein_denoised = get_denoised_expression(model)
fg_prob = get_foreground_probability(model)
```

**Available functions:**
- `setup_mudata_totalvi()` / `setup_anndata_totalvi()` - Register data
- `train_totalvi_model()` - Train model with custom parameters
- `get_latent_representation()` - Extract joint latent space
- `get_denoised_expression()` - Get background-corrected RNA and protein
- `get_foreground_probability()` - Probability protein signal is real
- `run_differential_expression()` - DE analysis between groups
- `validate_citeseq_data()` - Check data quality before training

---

## Workflow Steps (Manual)

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

### Step 2: Load CITE-Seq Data

```python
# Option 1: Load MuData (recommended format)
mdata = muon.read_h5mu("path/to/citeseq_data.h5mu")

# Check modalities
print(f"Modalities: {list(mdata.mod.keys())}")
print(f"RNA: {mdata.mod['rna'].shape}")
print(f"Protein: {mdata.mod['prot'].shape}")

# Option 2: Load AnnData with protein in obsm
# adata = sc.read_h5ad("citeseq.h5ad")
# Protein should be in adata.obsm["protein_expression"]
```

**Data Requirements**:
- RNA counts (raw, not normalized) in `layers["counts"]` or `.X`
- Protein expression in separate modality or `obsm["protein_expression"]`
- Batch labels in `.obs` if multiple samples
- Shared cell barcodes across modalities

---

### Step 3: Prepare Data for TotalVI

```python
# Access modalities
rna = mdata.mod["rna"]
protein = mdata.mod["prot"]

# CRITICAL: Convert sparse to dense (TotalVI requirement)
if hasattr(protein.X, 'toarray'):
    protein.X = protein.X.toarray()
if hasattr(rna.X, 'toarray'):
    rna.X = rna.X.toarray()

# Ensure counts layer exists for RNA
if "counts" not in rna.layers:
    rna.layers["counts"] = rna.X.copy()
else:
    rna.layers["counts"] = rna.layers["counts"].toarray() if hasattr(
        rna.layers["counts"], 'toarray') else rna.layers["counts"]

# Check protein names
print(f"Proteins measured: {list(protein.var_names)}")

# If multiple batches, verify shared proteins
if "batch" in rna.obs.columns:
    for batch in rna.obs["batch"].unique():
        batch_cells = rna.obs["batch"] == batch
        n_proteins = (protein[batch_cells].X.sum(axis=0) > 0).sum()
        print(f"Batch {batch}: {n_proteins} proteins detected")
```

---

### Step 4: Setup MuData for TotalVI

```python
# Setup for MuData format
scvi.model.TOTALVI.setup_mudata(
    mdata,
    rna_layer="counts",
    protein_layer=None,  # Uses .X from protein modality
    batch_key="batch",   # Optional: for batch correction
    modalities={
        "rna_layer": "rna",      # Name of RNA modality
        "protein_layer": "prot", # Name of protein modality
        "batch_key": "rna"       # Where batch info is stored
    }
)

# Alternative: Setup for AnnData with protein in obsm
# scvi.model.TOTALVI.setup_anndata(
#     adata,
#     layer="counts",
#     batch_key="batch",
#     protein_expression_obsm_key="protein_expression"
# )

print("Data registered for TotalVI")
```

---

### Step 5: Initialize and Train TotalVI Model

```python
# Initialize model
model = scvi.model.TOTALVI(mdata)

# View model summary
print(model)

# Train with early stopping (recommended)
model.train(early_stopping=True)

# Plot training history
plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.plot(model.history['elbo_train'].values, label='Train')
if 'elbo_validation' in model.history:
    plt.plot(model.history['elbo_validation'].values, label='Validation')
plt.xlabel('Epoch')
plt.ylabel('ELBO')
plt.legend()
plt.title('TotalVI Training')

plt.subplot(1, 2, 2)
plt.plot(model.history['reconstruction_loss_train'].values)
plt.xlabel('Epoch')
plt.ylabel('Reconstruction Loss')
plt.title('Reconstruction Loss')
plt.tight_layout()
plt.show()
```

**Training Notes**:
- Default ~400 epochs with early stopping
- ~3 minutes for 10-15k cells
- Watch for ELBO stabilization

---

### Step 6: Extract Latent Representation

```python
# Get joint latent representation
latent = model.get_latent_representation()

# Store in RNA modality for downstream analysis
rna.obsm["X_totalVI"] = latent

print(f"Latent representation shape: {latent.shape}")
```

---

### Step 7: Get Denoised Expression Values

```python
# Get denoised RNA and protein
rna_denoised, protein_denoised = model.get_normalized_expression(
    n_samples=25,      # Monte Carlo samples
    return_mean=True   # Return mean of samples
)

# Store denoised values
rna.layers["denoised_rna"] = rna_denoised
protein.layers["denoised_protein"] = protein_denoised

# Get protein foreground probability
# (probability that signal is real, not background)
protein_fg_prob = model.get_protein_foreground_probability(
    n_samples=25,
    return_mean=True
)
protein.layers["foreground_prob"] = protein_fg_prob

print("Denoised values extracted")
print(f"Protein foreground probability range: {protein_fg_prob.min():.2f} - {protein_fg_prob.max():.2f}")
```

**Key Outputs**:
- `denoised_rna`: Background-corrected RNA expression
- `denoised_protein`: Background-corrected protein expression
- `foreground_prob`: Probability protein signal is real (0-1)

---

### Step 8: Clustering and Visualization

```python
# Build neighborhood graph on joint latent space
sc.pp.neighbors(rna, use_rep="X_totalVI")

# UMAP embedding
sc.tl.umap(rna)

# Leiden clustering
sc.tl.leiden(rna, key_added="leiden_totalVI", resolution=0.5)

# Visualize
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Clusters
sc.pl.umap(rna, color="leiden_totalVI", ax=axes[0], show=False, title="Clusters")

# Batch (if present)
if "batch" in rna.obs.columns:
    sc.pl.umap(rna, color="batch", ax=axes[1], show=False, title="Batch")
else:
    axes[1].set_visible(False)

# Key protein marker (use denoised values)
# Transform for visualization
protein_for_viz = np.log1p(protein_denoised)
rna.obs["CD3_denoised"] = protein_for_viz[:, protein.var_names.get_loc("CD3")]
sc.pl.umap(rna, color="CD3_denoised", ax=axes[2], show=False,
           title="CD3 (denoised)", vmax="p99")

plt.tight_layout()
plt.show()
```

---

### Step 9: Differential Expression Analysis

```python
# Differential expression between clusters
# TotalVI tests both RNA and protein together
de_results = model.differential_expression(
    groupby="rna:leiden_totalVI",  # Format: modality:column
    group1="0",                     # Compare cluster 0
    group2="1",                     # vs cluster 1
    delta=0.5,                      # Effect size threshold
    batch_correction=True           # Account for batch effects
)

# View results
print(de_results.head(10))

# Key columns:
# - is_de_fdr: Significant after FDR correction
# - bayes_factor: Evidence strength
# - lfc_mean: Log fold change
# - raw_normalized_mean1/2: Expression in each group

# Filter for significant markers
# Protein markers
protein_markers = de_results[
    (de_results["is_de_fdr"]) &
    (de_results["bayes_factor"] > 0.7) &
    (de_results["lfc_mean"] > 0) &
    (de_results.index.isin(protein.var_names))
]

# RNA markers
rna_markers = de_results[
    (de_results["is_de_fdr"]) &
    (de_results["bayes_factor"] > 3) &
    (de_results["non_zeros_proportion1"] > 0.1) &
    (de_results["lfc_mean"] > 0) &
    (~de_results.index.isin(protein.var_names))
]

print(f"Protein markers: {len(protein_markers)}")
print(f"RNA markers: {len(rna_markers)}")
```

---

### Step 10: Save Results

```python
# Save trained model
model_dir = "totalvi_model"
model.save(model_dir, overwrite=True)

# Save processed MuData
mdata.write_h5mu("citeseq_totalvi_analyzed.h5mu")

# Export DE results
de_results.to_csv("differential_expression.csv")

# Reload model later
# model = scvi.model.TOTALVI.load(model_dir, adata=mdata)
```

---

## Key Parameters Reference

### Setup Parameters

| Parameter | Description |
|-----------|-------------|
| `rna_layer` | Layer with RNA counts |
| `protein_layer` | Layer with protein (None uses .X) |
| `batch_key` | Column for batch labels |
| `protein_expression_obsm_key` | For AnnData format |

### Model Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_latent` | 20 | Latent space dimensions |
| `n_layers_encoder` | 2 | Encoder depth |
| `n_layers_decoder` | 1 | Decoder depth |
| `latent_distribution` | "normal" | Latent prior |

### Training Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_epochs` | 400 | Maximum training epochs |
| `early_stopping` | True | Stop on validation plateau |
| `batch_size` | 128 | Samples per batch |

### Differential Expression Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `groupby` | None | Column for comparison |
| `delta` | 0.25 | Effect size threshold |
| `batch_correction` | False | Account for batch |
| `fdr_target` | 0.05 | FDR threshold |

---

## Parameter Tuning Guide

### When to Adjust Parameters

**Before running, ask the user:**
1. How many cells and proteins?
2. How many batches?
3. Quality of protein staining?
4. Looking for subtle or major differences?

### Latent Dimensions (`n_latent`)

| Dataset Complexity | n_latent | Notes |
|--------------------|----------|-------|
| Simple (<5 cell types) | 10-15 | Reduce overfitting |
| Standard | 20 | Default works well |
| Complex (>15 cell types) | 25-30 | Capture more variation |

### Clustering Resolution

| Observation | Adjustment |
|-------------|------------|
| Too few clusters | Increase to 0.8-1.0 |
| Too many clusters | Decrease to 0.3-0.4 |
| Protein markers not distinguishing | Check denoised values |

### Differential Expression Thresholds

| Analysis Goal | delta | bayes_factor (protein) | bayes_factor (RNA) |
|---------------|-------|----------------------|-------------------|
| Exploratory | 0.25 | 0.5 | 2 |
| Standard | 0.5 | 0.7 | 3 |
| Stringent | 1.0 | 0.9 | 5 |

---

## Adaptation Prompts for Claude

When a user invokes this skill, consider asking:

1. **Data characteristics:**
   - "How many cells and proteins are in your dataset?"
   - "Is this a single batch or multiple batches?"
   - "Do all batches have the same proteins measured?"

2. **Analysis goals:**
   - "Are you looking for cell type markers?"
   - "Do you need to compare conditions?"
   - "Will you be integrating with other datasets later?"

3. **Quality considerations:**
   - "How does your protein staining quality look?"
   - "Are there known problematic antibodies?"
   - "Have you done QC on the protein data?"

### Data Validation Helper

```python
def validate_citeseq_data(mdata):
    """Validate CITE-Seq data for TotalVI."""
    issues = []
    recommendations = []

    rna = mdata.mod["rna"]
    protein = mdata.mod["prot"]

    # Check for dense matrices
    if hasattr(protein.X, 'toarray'):
        issues.append("Protein matrix is sparse - convert to dense")

    # Check cell counts
    if rna.n_obs < 1000:
        recommendations.append(f"Low cell count ({rna.n_obs}). Results may be noisy.")

    # Check protein count
    if protein.n_vars < 10:
        recommendations.append(f"Few proteins ({protein.n_vars}). Consider if TotalVI is needed.")
    if protein.n_vars > 200:
        recommendations.append(f"Many proteins ({protein.n_vars}). May need longer training.")

    # Check for zeros
    zero_proteins = (protein.X.sum(axis=0) == 0).sum()
    if zero_proteins > 0:
        issues.append(f"{zero_proteins} proteins have zero counts in all cells")

    # Check batch info
    if "batch" in rna.obs.columns:
        batch_sizes = rna.obs["batch"].value_counts()
        if batch_sizes.min() < 100:
            recommendations.append("Some batches have <100 cells. May affect batch correction.")

    print("CITE-Seq Data Validation:")
    for issue in issues:
        print(f"  ERROR: {issue}")
    for rec in recommendations:
        print(f"  NOTE: {rec}")

    return len(issues) == 0
```

---

## Troubleshooting

### Common Issues

1. **"Sparse matrix" error**:
   - TotalVI requires dense matrices
   - Fix: `protein.X = protein.X.toarray()`

2. **Poor clustering**:
   - Check if proteins are informative for cell types
   - Try adjusting resolution
   - Verify batch correction is working

3. **Protein values look noisy**:
   - Use denoised values, not raw
   - Check foreground probability
   - Some proteins may have poor antibodies

4. **Batch effects persist**:
   - Verify `batch_key` is correctly specified
   - Check if batches have different protein panels
   - Subset to shared proteins

5. **DE returns few results**:
   - Lower `delta` threshold
   - Check if groups are truly different
   - Ensure sufficient cells per group

### Interpreting Foreground Probability

```python
# Foreground probability interpretation
# 0.0-0.3: Likely background noise
# 0.3-0.7: Uncertain
# 0.7-1.0: Likely true signal

# Visualize for a specific protein
protein_name = "CD3"
idx = protein.var_names.get_loc(protein_name)
fg_prob = protein.layers["foreground_prob"][:, idx]

plt.figure(figsize=(8, 4))
plt.hist(fg_prob, bins=50)
plt.xlabel(f"{protein_name} Foreground Probability")
plt.ylabel("Cells")
plt.axvline(0.5, color='r', linestyle='--', label='Threshold')
plt.legend()
plt.title(f"{protein_name}: Real vs Background Signal")
plt.show()
```

---

## Advanced Usage

### Custom Model Architecture

```python
model = scvi.model.TOTALVI(
    mdata,
    n_latent=25,           # More latent dimensions
    n_layers_encoder=3,    # Deeper encoder
    n_layers_decoder=2,    # Deeper decoder
    latent_distribution="normal",
    gene_likelihood="nb"   # Negative binomial for RNA
)

model.train(
    max_epochs=500,
    early_stopping=True,
    early_stopping_patience=30
)
```

### Compare Denoised vs Raw Protein

```python
# Compare denoised to raw for a marker
protein_name = "CD19"
idx = protein.var_names.get_loc(protein_name)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Raw
sc.pl.umap(rna, color=None, ax=axes[0], show=False)
scatter = axes[0].scatter(
    rna.obsm["X_umap"][:, 0],
    rna.obsm["X_umap"][:, 1],
    c=protein.X[:, idx],
    s=1, cmap="viridis"
)
plt.colorbar(scatter, ax=axes[0])
axes[0].set_title(f"{protein_name} (Raw)")

# Denoised
scatter = axes[1].scatter(
    rna.obsm["X_umap"][:, 0],
    rna.obsm["X_umap"][:, 1],
    c=np.log1p(protein.layers["denoised_protein"][:, idx]),
    s=1, cmap="viridis"
)
plt.colorbar(scatter, ax=axes[1])
axes[1].set_title(f"{protein_name} (Denoised, log)")

plt.tight_layout()
plt.show()
```

---

## References

- [TotalVI Tutorial](https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/multimodal/totalVI.html)
- [TotalVI Paper](https://www.nature.com/articles/s41592-020-01050-x) - Gayoso et al., Nature Methods 2021
- [scvi-tools Documentation](https://docs.scvi-tools.org/)
- [MuData Documentation](https://mudata.readthedocs.io/)
