# CITE-Seq and RNA-Only Integration with TotalVI

## Overview

TotalVI can integrate CITE-Seq with RNA-only data by:

1. **Joint modeling**: Learn from both complete (CITE-Seq) and incomplete (RNA-only) data
2. **Missing data handling**: Setting protein values to 0 signals "missing" to TotalVI
3. **Protein imputation**: Use learned RNA-protein relationships to predict missing proteins
4. **Unified embedding**: Single latent space across all cells regardless of measured modalities

**Key Insight**: "The quality of totalVI's protein imputation depends on datasets largely sharing cell subpopulations."

**Citation**: Gayoso et al. (2021). "Joint probabilistic modeling of single-cell multi-omic data with totalVI." Nature Methods.

---

## Workflow Steps

### Step 1: Environment Setup

```python
import numpy as np
import scanpy as sc
import scvi
import pandas as pd
import matplotlib.pyplot as plt
import torch

scvi.settings.seed = 0
sc.set_figure_params(figsize=(6, 6), frameon=False)
torch.set_float32_matmul_precision("high")
```

---

### Step 2: Load CITE-Seq and RNA-Only Datasets

```python
# Load CITE-Seq data (has both RNA and protein)
adata_cite = sc.read_h5ad("path/to/citeseq.h5ad")
print(f"CITE-Seq: {adata_cite.n_obs} cells")
print(f"Proteins: {adata_cite.obsm['protein_expression'].shape[1]}")

# Load RNA-only data (no protein)
adata_rna = sc.read_h5ad("path/to/rnaseq.h5ad")
print(f"RNA-only: {adata_rna.n_obs} cells")

# Add batch labels
adata_cite.obs["batch"] = "CITE"
adata_rna.obs["batch"] = "RNA_only"
adata_cite.obs["has_protein"] = True
adata_rna.obs["has_protein"] = False
```

---

### Step 3: Align and Merge Datasets

```python
# Find shared genes
shared_genes = adata_cite.var_names.intersection(adata_rna.var_names)
print(f"Shared genes: {len(shared_genes)}")

# Subset to shared genes
adata_cite = adata_cite[:, shared_genes].copy()
adata_rna = adata_rna[:, shared_genes].copy()

# Concatenate datasets
adata = sc.concat(
    [adata_cite, adata_rna],
    join="inner",
    label="batch",
    keys=["CITE", "RNA_only"]
)

# Ensure counts layer
if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()

print(f"Combined: {adata.n_obs} cells, {adata.n_vars} genes")
```

---

### Step 4: Handle Missing Proteins (CRITICAL)

```python
# Get protein info from CITE-Seq data
protein_names = adata_cite.obsm["protein_expression"].columns \
    if hasattr(adata_cite.obsm["protein_expression"], 'columns') \
    else [f"Protein_{i}" for i in range(adata_cite.obsm["protein_expression"].shape[1])]

n_proteins = len(protein_names)

# Create protein DataFrame for combined data
protein_df = pd.DataFrame(
    index=adata.obs_names,
    columns=protein_names,
    dtype=float
)

# Fill CITE-Seq cells with actual protein values
cite_mask = adata.obs["batch"] == "CITE"
protein_df.loc[cite_mask] = adata_cite.obsm["protein_expression"].values

# CRITICAL: Set RNA-only cells to 0 (signals "missing" to TotalVI)
rna_mask = adata.obs["batch"] == "RNA_only"
protein_df.loc[rna_mask] = 0.0

# Store in obsm
adata.obsm["protein_expression"] = protein_df

print(f"Protein matrix: {adata.obsm['protein_expression'].shape}")
print(f"CITE-Seq cells: {cite_mask.sum()} (with protein)")
print(f"RNA-only cells: {rna_mask.sum()} (protein set to 0 = missing)")
```

**Key Point**: Setting protein to 0 tells TotalVI "this data is missing" - it's not treated as zero expression.

---

### Step 5: Gene Selection

```python
# Select highly variable genes (batch-aware)
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=4000,
    batch_key="batch",
    flavor="seurat_v3",
    subset=True,
    layer="counts"
)

print(f"Genes after HVG selection: {adata.n_vars}")
```

---

### Step 6: Setup and Train TotalVI

```python
# Setup TotalVI
scvi.model.TOTALVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="batch",
    protein_expression_obsm_key="protein_expression"
)

# Initialize model
model = scvi.model.TOTALVI(
    adata,
    latent_distribution="normal",
    n_layers_decoder=2
)

# Train
model.train(max_epochs=400, early_stopping=True)

# Plot training
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

---

### Step 7: Get Latent Representation

```python
# Get joint latent space (all cells)
latent = model.get_latent_representation()
adata.obsm["X_totalVI"] = latent

print(f"Latent representation: {latent.shape}")
```

---

### Step 8: Impute Proteins for RNA-Only Cells

```python
# Counterfactual imputation: "What if RNA-only cells had protein measured?"
# Use transform_batch to specify source batch for protein prediction

_, protein_imputed = model.get_normalized_expression(
    transform_batch="CITE",  # Predict as if from CITE-Seq batch
    n_samples=25,
    return_mean=True
)

# Store imputed proteins
adata.obsm["protein_imputed"] = protein_imputed

# Compare imputed vs observed for CITE-Seq cells
cite_cells = adata.obs["batch"] == "CITE"
print("\nCITE-Seq cells: Imputed vs Observed correlation")
for i, prot in enumerate(protein_names[:5]):
    observed = adata.obsm["protein_expression"].values[cite_cells, i]
    imputed = protein_imputed[cite_cells, i]
    corr = np.corrcoef(observed, imputed)[0, 1]
    print(f"  {prot}: r = {corr:.3f}")
```

---

### Step 9: Cluster and Visualize

```python
# Build neighborhood graph
sc.pp.neighbors(adata, use_rep="X_totalVI")
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="leiden_totalVI", resolution=0.5)

# Visualization
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# Row 1: Dataset integration
sc.pl.umap(adata, color="batch", ax=axes[0, 0], show=False,
           title="Dataset (CITE vs RNA-only)")
sc.pl.umap(adata, color="leiden_totalVI", ax=axes[0, 1], show=False,
           title="Clusters")
sc.pl.umap(adata, color="has_protein", ax=axes[0, 2], show=False,
           title="Has Protein Measured")

# Row 2: Imputed proteins (for all cells)
for i, prot in enumerate(protein_names[:3]):
    adata.obs[f"imputed_{prot}"] = np.log1p(protein_imputed[:, i])
    sc.pl.umap(adata, color=f"imputed_{prot}", ax=axes[1, i],
               show=False, title=f"{prot} (imputed)", cmap="viridis")

plt.tight_layout()
plt.show()
```

---

### Step 10: Validate and Save Results

```python
# If you held out some proteins for validation
def validate_imputation(adata, protein_imputed, protein_names, held_out_data):
    """Validate imputation against held-out ground truth."""
    from scipy.stats import pearsonr

    print("Imputation Validation (RNA-only cells with held-out protein):")
    for i, prot in enumerate(protein_names):
        if prot in held_out_data.columns:
            true = held_out_data[prot].values
            pred = protein_imputed[~adata.obs["has_protein"], i]
            r, p = pearsonr(true, pred)
            print(f"  {prot}: Pearson r = {r:.3f}, p = {p:.2e}")

# Save results
adata.write_h5ad("integrated_cite_rna.h5ad")
model.save("totalvi_integration_model", overwrite=True)

# Export imputed proteins
imputed_df = pd.DataFrame(
    protein_imputed,
    index=adata.obs_names,
    columns=protein_names
)
imputed_df.to_csv("imputed_proteins.csv")

print("Results saved")
```

---

## Key Parameters Reference

### Data Setup

| Parameter | Description |
|-----------|-------------|
| `layer` | RNA counts layer |
| `batch_key` | Distinguishes CITE vs RNA-only |
| `protein_expression_obsm_key` | Protein matrix (0s = missing) |

### Model Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `latent_distribution` | "normal" | Latent prior |
| `n_layers_decoder` | 1 | Decoder depth (2 often helps) |
| `n_latent` | 20 | Latent dimensions |

### Imputation Parameters

| Parameter | Description |
|-----------|-------------|
| `transform_batch` | Source batch for imputation |
| `n_samples` | Monte Carlo samples |
| `include_protein_background` | Include technical noise |

---

## Parameter Tuning Guide

### When to Adjust Parameters

**Before running, ask the user:**
1. Ratio of CITE-Seq to RNA-only cells?
2. Do the datasets share cell types?
3. How many proteins in CITE-Seq panel?
4. Quality of protein staining?

### Cell Ratio Considerations

| Ratio (CITE : RNA) | Considerations |
|--------------------|----------------|
| >1:1 | Good protein learning |
| 1:1 | Balanced |
| 1:2+ | Protein imputation may be less accurate |
| 1:10+ | Consider training on CITE-Seq first |

### Decoder Depth

| Scenario | n_layers_decoder | Notes |
|----------|------------------|-------|
| Standard | 1 | Default |
| Complex relationships | 2 | Better RNA-protein mapping |
| Very large dataset | 2-3 | More capacity |

### Imputation Quality Factors

| Factor | Impact on Quality |
|--------|-------------------|
| Shared cell types | High - critical! |
| Protein panel size | More proteins = harder |
| Sequencing depth | Higher = better |
| Batch effects | Can reduce accuracy |

---

## Adaptation Prompts for Claude

When a user invokes this skill, consider asking:

1. **Dataset composition:**
   - "How many cells in each dataset (CITE vs RNA-only)?"
   - "Do the datasets come from the same tissue?"
   - "What cell types are expected in each?"

2. **Protein panel:**
   - "How many proteins in the CITE-Seq panel?"
   - "Are there specific proteins you care most about?"
   - "Can you share the protein list?"

3. **Validation:**
   - "Do you have any RNA-only cells with held-out protein for validation?"
   - "What accuracy metrics matter for your application?"

### Quality Check Helper

```python
def assess_integration_quality(adata):
    """Assess CITE-RNA integration quality."""

    # Check batch mixing
    print("Integration Assessment:")

    # Cluster composition by batch
    ct = pd.crosstab(adata.obs["leiden_totalVI"], adata.obs["batch"])
    ct_norm = ct.div(ct.sum(axis=1), axis=0)

    # Check for batch-dominated clusters
    cite_dominated = (ct_norm["CITE"] > 0.9).sum()
    rna_dominated = (ct_norm["RNA_only"] > 0.9).sum()

    print(f"  Total clusters: {len(ct_norm)}")
    print(f"  CITE-dominated clusters (>90%): {cite_dominated}")
    print(f"  RNA-dominated clusters (>90%): {rna_dominated}")

    if cite_dominated + rna_dominated > len(ct_norm) * 0.3:
        print("  WARNING: Poor batch mixing - imputation may be less reliable")
    else:
        print("  Good batch mixing - imputation should be reliable")

    return ct_norm
```

---

## Troubleshooting

### Common Issues

1. **Imputed proteins all look similar**:
   - Datasets may not share cell types
   - Check batch integration in UMAP
   - Verify protein values set to 0 for RNA-only

2. **Batches don't mix in latent space**:
   - Strong batch effects
   - Different cell type composition
   - Try more HVGs or different preprocessing

3. **Low correlation with held-out validation**:
   - Cell types may differ
   - Protein panel may be noisy
   - Consider which proteins are most reliable

4. **Memory issues**:
   - Reduce number of genes
   - Process in chunks
   - Use smaller n_samples for imputation

### Verifying Missing Data Handling

```python
# Verify zeros are interpreted as missing
def check_missing_handling(adata):
    """Verify protein zeros are set correctly."""

    protein_matrix = adata.obsm["protein_expression"]
    cite_mask = adata.obs["batch"] == "CITE"
    rna_mask = adata.obs["batch"] == "RNA_only"

    cite_zeros = (protein_matrix[cite_mask] == 0).all(axis=1).sum()
    rna_zeros = (protein_matrix[rna_mask] == 0).all(axis=1).sum()

    print(f"CITE-Seq cells with all-zero protein: {cite_zeros}")
    print(f"RNA-only cells with all-zero protein: {rna_zeros}")

    if rna_zeros != rna_mask.sum():
        print("WARNING: Not all RNA-only cells have zero protein!")
    else:
        print("All RNA-only cells correctly have zero protein (= missing)")
```

---

## Advanced Usage

### Multiple RNA-Only Batches

```python
# If you have multiple RNA-only datasets
adata_cite.obs["batch"] = "CITE"
adata_rna1.obs["batch"] = "RNA_batch1"
adata_rna2.obs["batch"] = "RNA_batch2"

# Concatenate all
adata = sc.concat([adata_cite, adata_rna1, adata_rna2])

# Set protein to 0 for all RNA-only batches
rna_batches = ["RNA_batch1", "RNA_batch2"]
for batch in rna_batches:
    mask = adata.obs["batch"] == batch
    adata.obsm["protein_expression"].loc[mask] = 0

# Continue with training...
```

### Differential Expression Across Original Batches

```python
# After integration, compare cell types across original batches
de_results = model.differential_expression(
    groupby="leiden_totalVI",
    group1="0",
    group2="1",
    batch_correction=True
)

# Filter for protein markers
protein_de = de_results[de_results.index.isin(protein_names)]
print("Protein markers between clusters:")
print(protein_de[protein_de["is_de_fdr"]].head(10))
```

---

## References

- [CITE-RNA Integration Tutorial](https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/multimodal/cite_scrna_integration_w_totalVI.html)
- [TotalVI Paper](https://www.nature.com/articles/s41592-020-01050-x) - Gayoso et al., Nature Methods 2021
- [scvi-tools Documentation](https://docs.scvi-tools.org/)
