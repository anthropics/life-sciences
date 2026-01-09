# TotalVI Reference Mapping and Protein Imputation

## Overview

TotalVI reference mapping enables:

1. **Query mapping**: Project new data onto TotalVI reference latent space
2. **Protein imputation**: Predict protein expression for RNA-only queries
3. **Label transfer**: Transfer cell type annotations via classifiers
4. **Cross-modality prediction**: Leverage reference multimodal relationships

**Key Advantage**: Even if query data lacks protein measurements, you can impute them using the learned RNA-protein relationships from the reference.

**Citation**: Gayoso et al. (2021). "Joint probabilistic modeling of single-cell multi-omic data with totalVI." Nature Methods.

---

## Workflow Steps

### Step 1: Environment Setup

```python
import numpy as np
import scanpy as sc
import scvi
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
import torch

scvi.settings.seed = 0
sc.set_figure_params(figsize=(6, 6), frameon=False)
torch.set_float32_matmul_precision("high")
```

---

### Step 2: Prepare Reference CITE-Seq Data

```python
# Load CITE-Seq reference with RNA + protein
adata_ref = sc.read_h5ad("path/to/citeseq_reference.h5ad")

print(f"Reference: {adata_ref.n_obs} cells")
print(f"Genes: {adata_ref.n_vars}")
print(f"Proteins: {adata_ref.obsm['protein_expression'].shape[1]}")

# Ensure counts layer
if "counts" not in adata_ref.layers:
    adata_ref.layers["counts"] = adata_ref.X.copy()

# Preprocessing
sc.pp.normalize_total(adata_ref, target_sum=1e4)
sc.pp.log1p(adata_ref)
adata_ref.raw = adata_ref  # Store for later

# HVG selection
sc.pp.highly_variable_genes(
    adata_ref,
    n_top_genes=4000,
    batch_key="batch",
    flavor="seurat_v3",
    subset=True
)
```

---

### Step 3: Train Reference TotalVI Model with Protein Masking

**Key Strategy**: Mask proteins in some batches during training to teach the model to impute missing proteins.

```python
# Setup with protein masking for better generalization
scvi.model.TOTALVI.setup_anndata(
    adata_ref,
    layer="counts",
    batch_key="batch",
    protein_expression_obsm_key="protein_expression"
)

# Initialize TotalVI reference
model_ref = scvi.model.TOTALVI(
    adata_ref,
    use_layer_norm="both",      # scArches compatible
    use_batch_norm="none"       # scArches compatible
)

# Train
model_ref.train(max_epochs=400)

# Get reference latent and store
adata_ref.obsm["X_totalVI"] = model_ref.get_latent_representation()

# Save for query mapping
model_ref.save("totalvi_reference", overwrite=True)

print("Reference TotalVI model saved")
```

---

### Step 4: Train Cell Type Classifier

```python
# Train classifier on reference latent space
latent_ref = model_ref.get_latent_representation()

clf = RandomForestClassifier(
    n_estimators=100,
    random_state=0,
    n_jobs=-1
)
clf.fit(latent_ref, adata_ref.obs["cell_type"])

# Store reference UMAP for consistent projection
sc.pp.neighbors(adata_ref, use_rep="X_totalVI")
sc.tl.umap(adata_ref)

print(f"Classifier trained on {len(adata_ref.obs['cell_type'].unique())} cell types")
```

---

### Step 5: Prepare Query Data (RNA-Only or CITE-Seq)

```python
# Load query data
adata_query = sc.read_h5ad("path/to/query.h5ad")

print(f"Query: {adata_query.n_obs} cells, {adata_query.n_vars} genes")

# Check if query has protein
has_protein = "protein_expression" in adata_query.obsm
print(f"Query has protein: {has_protein}")

# Ensure counts
if "counts" not in adata_query.layers:
    adata_query.layers["counts"] = adata_query.X.copy()

# Same preprocessing as reference
sc.pp.normalize_total(adata_query, target_sum=1e4)
sc.pp.log1p(adata_query)

# Prepare query for TotalVI reference
scvi.model.TOTALVI.prepare_query_anndata(
    adata_query,
    "totalvi_reference"
)

# If query lacks protein, create empty protein matrix
if not has_protein:
    n_proteins = adata_ref.obsm["protein_expression"].shape[1]
    protein_names = adata_ref.uns.get("protein_names",
                                       [f"Protein_{i}" for i in range(n_proteins)])
    # Set to zeros (TotalVI interprets as missing)
    adata_query.obsm["protein_expression"] = np.zeros(
        (adata_query.n_obs, n_proteins)
    )
    print(f"Created empty protein matrix ({n_proteins} proteins)")
```

---

### Step 6: Map Query to Reference

```python
# Load query model from reference
model_query = scvi.model.TOTALVI.load_query_data(
    adata_query,
    "totalvi_reference"
)

# Train query model (preserves reference space)
model_query.train(
    max_epochs=200,
    plan_kwargs={"weight_decay": 0.0}  # Critical!
)

# Get query latent representation
adata_query.obsm["X_totalVI"] = model_query.get_latent_representation()

print("Query mapped to reference")
```

---

### Step 7: Impute Proteins for Query

```python
# Impute proteins using counterfactual prediction
# "What would protein expression be if query came from reference batch?"

# Get normalized expression with protein imputation
_, protein_imputed = model_query.get_normalized_expression(
    adata_query,
    n_samples=25,
    return_mean=True,
    transform_batch=adata_ref.obs["batch"].unique().tolist()  # Reference batches
)

# Store imputed proteins
adata_query.obsm["protein_imputed"] = protein_imputed

# Also get denoised RNA
rna_denoised, _ = model_query.get_normalized_expression(
    adata_query,
    n_samples=25,
    return_mean=True
)
adata_query.layers["denoised_rna"] = rna_denoised

print(f"Imputed {protein_imputed.shape[1]} proteins for {protein_imputed.shape[0]} cells")
```

**Key Parameter**: `transform_batch` specifies which reference batch(es) to use as the "source" for protein imputation.

---

### Step 8: Transfer Cell Type Labels

```python
# Predict cell types using reference classifier
latent_query = model_query.get_latent_representation()
predictions = clf.predict(latent_query)
probabilities = clf.predict_proba(latent_query)

# Store predictions
adata_query.obs["predicted_cell_type"] = predictions
adata_query.obs["prediction_confidence"] = probabilities.max(axis=1)

print("Predicted cell types:")
print(adata_query.obs["predicted_cell_type"].value_counts())
```

---

### Step 9: Visualize Results

```python
# Compute query UMAP
sc.pp.neighbors(adata_query, use_rep="X_totalVI")
sc.tl.umap(adata_query)

# Visualization
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# Row 1: Clustering and labels
sc.pl.umap(adata_query, color="predicted_cell_type", ax=axes[0, 0],
           show=False, title="Predicted Cell Types")
sc.pl.umap(adata_query, color="prediction_confidence", ax=axes[0, 1],
           show=False, title="Prediction Confidence", cmap="RdYlGn")

# Row 2: Imputed proteins
protein_names = adata_ref.uns.get("protein_names", ["CD3", "CD4", "CD8"])[:3]
for i, prot in enumerate(protein_names[:3]):
    if i < protein_imputed.shape[1]:
        adata_query.obs[f"imputed_{prot}"] = np.log1p(protein_imputed[:, i])
        sc.pl.umap(adata_query, color=f"imputed_{prot}", ax=axes[1, i],
                   show=False, title=f"{prot} (imputed)", cmap="viridis")

plt.tight_layout()
plt.show()
```

---

### Step 10: Save Results

```python
# Save mapped query with predictions and imputations
adata_query.write_h5ad("query_mapped_totalvi.h5ad")

# Save query model
model_query.save("query_totalvi_model", overwrite=True)

# Export imputed proteins as CSV
import pandas as pd
protein_df = pd.DataFrame(
    protein_imputed,
    index=adata_query.obs_names,
    columns=adata_ref.uns.get("protein_names", range(protein_imputed.shape[1]))
)
protein_df.to_csv("imputed_proteins.csv")

print("Results saved")
```

---

## Key Parameters Reference

### Reference Model Setup

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `use_layer_norm` | "both" | scArches compatibility |
| `use_batch_norm` | "none" | scArches compatibility |
| `protein_expression_obsm_key` | Key for protein data |

### Query Training

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `max_epochs` | 200 | Query adaptation |
| `weight_decay` | 0.0 | Preserve reference space |

### Protein Imputation

| Parameter | Description |
|-----------|-------------|
| `n_samples` | Monte Carlo samples (25 recommended) |
| `return_mean` | Return mean of samples |
| `transform_batch` | Which batch(es) to impute from |

---

## Parameter Tuning Guide

### When to Adjust Parameters

**Before running, ask the user:**
1. Is query RNA-only or CITE-Seq?
2. How different is query from reference (tissue, platform)?
3. How many proteins in reference?
4. What batch(es) should be used for imputation?

### Transform Batch Selection

| Scenario | transform_batch | Notes |
|----------|-----------------|-------|
| Query similar to one ref batch | That specific batch | Most accurate |
| Query different from all | All reference batches | Average behavior |
| Unsure | All batches | Safe default |

### Query Training Epochs

| Scenario | max_epochs | Notes |
|----------|------------|-------|
| Same platform as reference | 100 | Fast |
| Different platform | 200 | Standard |
| Very different conditions | 300 | More adaptation |

---

## Adaptation Prompts for Claude

When a user invokes this skill, consider asking:

1. **Query data type:**
   - "Does your query have protein measurements or RNA-only?"
   - "What platform generated your query data?"
   - "How many cells in query?"

2. **Reference details:**
   - "What proteins are in the reference panel?"
   - "How many batches in reference?"
   - "What cell types are annotated?"

3. **Imputation needs:**
   - "Which proteins are most important for your analysis?"
   - "Do you need all proteins or specific ones?"

### Imputation Quality Assessment

```python
def assess_imputation_quality(protein_imputed, predicted_cell_types):
    """Assess quality of protein imputation."""

    # Check for reasonable distributions
    print("Imputed Protein Statistics:")
    for i in range(min(5, protein_imputed.shape[1])):
        prot_vals = protein_imputed[:, i]
        print(f"  Protein {i}: mean={prot_vals.mean():.2f}, "
              f"std={prot_vals.std():.2f}, "
              f"zeros={100*(prot_vals < 0.1).mean():.1f}%")

    # Check cell-type-specific patterns
    print("\nProtein patterns by cell type:")
    for ct in np.unique(predicted_cell_types)[:5]:
        mask = predicted_cell_types == ct
        mean_expr = protein_imputed[mask].mean(axis=0)
        top_protein = mean_expr.argmax()
        print(f"  {ct}: Highest protein = {top_protein} "
              f"(mean={mean_expr[top_protein]:.2f})")
```

---

## Troubleshooting

### Common Issues

1. **Imputed proteins all similar**:
   - Query may not align well with reference
   - Check latent space overlap
   - Try different transform_batch

2. **Poor cell type predictions**:
   - Cell types may not exist in reference
   - Query from very different tissue
   - Increase training epochs

3. **Gene mismatch errors**:
   - Run `prepare_query_anndata` first
   - Check gene name format (symbols vs IDs)

4. **Memory issues with large query**:
   - Process in batches
   - Reduce n_samples for imputation

### Validating Imputation

```python
# If you have held-out proteins to validate
def validate_imputation(true_protein, imputed_protein):
    """Calculate correlation between true and imputed."""
    from scipy.stats import pearsonr, spearmanr

    correlations = []
    for i in range(true_protein.shape[1]):
        r, p = pearsonr(true_protein[:, i], imputed_protein[:, i])
        correlations.append(r)

    print(f"Mean Pearson correlation: {np.mean(correlations):.3f}")
    print(f"Median Pearson correlation: {np.median(correlations):.3f}")
    return correlations
```

---

## Advanced Usage

### Batch-Specific Imputation

```python
# Impute as if from specific reference batch
# Useful when query is similar to one batch

specific_batch = "PBMC10k"
_, protein_batch_specific = model_query.get_normalized_expression(
    adata_query,
    n_samples=25,
    return_mean=True,
    transform_batch=specific_batch
)
```

### Combined Reference-Query Visualization

```python
# Project query onto reference UMAP
import anndata

adata_combined = anndata.concat(
    [adata_ref, adata_query],
    label="source",
    keys=["reference", "query"]
)

# Get combined latent through query model
adata_combined.obsm["X_totalVI"] = model_query.get_latent_representation(adata_combined)

# UMAP
sc.pp.neighbors(adata_combined, use_rep="X_totalVI")
sc.tl.umap(adata_combined)

sc.pl.umap(adata_combined, color=["source", "cell_type"],
           ncols=2, frameon=False)
```

---

## References

- [TotalVI Reference Mapping Tutorial](https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/multimodal/totalVI_reference_mapping.html)
- [TotalVI Paper](https://www.nature.com/articles/s41592-020-01050-x) - Gayoso et al., Nature Methods 2021
- [scvi-tools Documentation](https://docs.scvi-tools.org/)
