# Reference Atlas Mapping with scArches

## Overview

scArches (single-cell architecture surgery) enables:

1. **Reference building**: Train models with settings that allow future query mapping
2. **Query mapping**: Map new data to existing reference latent space
3. **Label transfer**: Transfer cell type annotations via trained classifiers
4. **Atlas-level integration**: Incrementally add new datasets to atlases

**Key Advantages**:
- Preserves reference latent space exactly
- Efficient (no need to retrain full model)
- Works with scVI, SCANVI, TotalVI, and other scvi-tools models
- Enables atlas-scale integration

**Critical Requirement**: Reference models must be trained with `encode_covariates=True` to enable scArches functionality.

**Citation**: Lotfollahi et al. (2022). "Mapping single-cell data to reference atlases by transfer learning." Nature Biotechnology.

---

## Approach 1: Complete Analysis Pipeline (Recommended)

For standard scArches mapping, use the convenience script `scripts/scarches_analysis.py`:

```bash
# Train new reference and map query
python3 scripts/scarches_analysis.py --reference atlas.h5ad --query new_data.h5ad --labels-key cell_type

# Use existing reference model
python3 scripts/scarches_analysis.py --reference-model-path ./ref_model --reference atlas.h5ad --query new_data.h5ad --labels-key cell_type

# With batch correction in reference
python3 scripts/scarches_analysis.py --reference atlas.h5ad --query new_data.h5ad --labels-key cell_type --batch-key batch
```

**Parameters:**

Input:
- `--reference` - Reference .h5ad (if training new model)
- `--reference-model-path` - Path to existing reference model
- `--query` - Query .h5ad to map (required)
- `--labels-key` - Cell type column (required)

Reference model:
- `--n-latent` - Latent dimensions (default: 30)
- `--n-layers` - Hidden layers (default: 2)
- `--ref-epochs` - Reference training epochs (default: 400)
- `--batch-key` - Batch column for reference

Query mapping:
- `--query-epochs` - Query training epochs (default: 200)
- `--n-hvg` - HVGs for reference (default: 2000)

**Outputs:**
- `combined_umap.png` - Combined reference + query visualization
- `confidence_distribution.png` - Label transfer confidence
- `predictions.csv` - Cell type predictions
- `reference_model/` - Saved reference model
- `query_model/` - Saved query model
- `*_mapped.h5ad` - Mapped query data
- `combined.h5ad` - Combined dataset

---

## Approach 2: Modular Building Blocks (For Custom Workflows)

For custom workflows, use functions from `scripts/scarches_core.py`:

```python
import sys
sys.path.append('scripts/')
from scarches_core import (
    setup_reference_data,
    train_reference_model,
    train_classifier,
    prepare_query_data,
    map_query_to_reference,
    transfer_labels,
    verify_reference_preservation,
)

# Train scArches-compatible reference
setup_reference_data(adata_ref, layer='counts', batch_key='batch')
model_ref = train_reference_model(adata_ref, n_latent=30)  # Uses scArches-compatible settings

# Train classifier
clf, accuracy = train_classifier(model_ref, adata_ref, labels_key='cell_type')

# Map query
prepare_query_data(adata_query, 'reference_model')
model_query = map_query_to_reference(adata_query, 'reference_model')

# Transfer labels
predictions, confidence = transfer_labels(clf, model_query, adata_query)

# Verify reference space preserved
verify_reference_preservation(model_ref, model_query, adata_ref)
```

**Available functions:**
- `setup_reference_data()` - Register reference data
- `train_reference_model()` - Train scArches-compatible model (uses correct settings)
- `train_classifier()` - Train label transfer classifier
- `prepare_query_data()` - Align query genes with reference
- `map_query_to_reference()` - Fine-tune on query (preserves reference space)
- `transfer_labels()` - Predict cell types for query
- `verify_reference_preservation()` - Check reference embeddings unchanged

---

## Workflow Steps (Manual)

### Step 1: Environment Setup

```python
import numpy as np
import scanpy as sc
import scvi
import anndata
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, accuracy_score
import torch

scvi.settings.seed = 0
sc.set_figure_params(figsize=(6, 6), frameon=False)
torch.set_float32_matmul_precision("high")
```

---

### Step 2: Prepare Reference Data

```python
# Load reference data with cell type annotations
adata_ref = sc.read_h5ad("path/to/reference.h5ad")

print(f"Reference: {adata_ref.n_obs} cells, {adata_ref.n_vars} genes")
print(f"Cell types: {adata_ref.obs['cell_type'].value_counts()}")

# Ensure counts are available
if "counts" not in adata_ref.layers:
    adata_ref.layers["counts"] = adata_ref.X.copy()

# Select highly variable genes (on reference only)
sc.pp.highly_variable_genes(
    adata_ref,
    n_top_genes=2000,
    batch_key="batch",  # If multiple batches in reference
    flavor="seurat_v3",
    subset=True,
    layer="counts"
)

print(f"Genes after HVG selection: {adata_ref.n_vars}")
```

---

### Step 3: Train scArches-Compatible Reference Model

**CRITICAL**: The model architecture must be configured for scArches compatibility.

```python
# Setup reference data
scvi.model.SCVI.setup_anndata(
    adata_ref,
    layer="counts",
    batch_key="batch"  # Important for batch-aware mapping
)

# Initialize model with scArches-compatible settings
model_ref = scvi.model.SCVI(
    adata_ref,
    use_layer_norm="both",      # CRITICAL for scArches
    use_batch_norm="none",      # CRITICAL for scArches
    encode_covariates=True,     # CRITICAL - enables batch mapping
    dropout_rate=0.2,
    n_layers=2,
    n_latent=30
)

# Train reference model
model_ref.train(max_epochs=400)

# Get reference latent representation
adata_ref.obsm["X_scVI"] = model_ref.get_latent_representation()

# Save model for future query mapping
model_ref.save("reference_model", overwrite=True)

print("Reference model trained and saved")
```

**Why These Settings?**
- `use_layer_norm="both"`: Stabilizes representations across batches
- `use_batch_norm="none"`: Prevents batch-specific normalization
- `encode_covariates=True`: Allows new batches to be encoded

---

### Step 4: Train Classifier for Label Transfer (Optional)

```python
# Train classifier on reference latent space
latent_ref = model_ref.get_latent_representation()

clf = RandomForestClassifier(
    n_estimators=100,
    random_state=0,
    n_jobs=-1
)
clf.fit(latent_ref, adata_ref.obs["cell_type"])

# Evaluate on reference (should be high)
ref_predictions = clf.predict(latent_ref)
ref_accuracy = accuracy_score(adata_ref.obs["cell_type"], ref_predictions)
print(f"Reference self-accuracy: {ref_accuracy:.3f}")
```

---

### Step 5: Prepare Query Data

```python
# Load query data
adata_query = sc.read_h5ad("path/to/query.h5ad")

print(f"Query: {adata_query.n_obs} cells, {adata_query.n_vars} genes")

# Ensure counts
if "counts" not in adata_query.layers:
    adata_query.layers["counts"] = adata_query.X.copy()

# CRITICAL: Subset query to reference genes
# prepare_query_anndata handles gene alignment
scvi.model.SCVI.prepare_query_anndata(
    adata_query,
    "reference_model"  # Path to saved reference model
)

print(f"Query genes after alignment: {adata_query.n_vars}")
```

**What `prepare_query_anndata` does**:
- Reorders genes to match reference
- Pads missing genes with zeros
- Validates data structure

---

### Step 6: Load and Fine-tune Query Model

```python
# Load query model from reference
model_query = scvi.model.SCVI.load_query_data(
    adata_query,
    "reference_model"
)

# Fine-tune on query data
# CRITICAL: weight_decay=0.0 preserves reference latent space
model_query.train(
    max_epochs=200,
    plan_kwargs={"weight_decay": 0.0}  # CRITICAL!
)

# Get query latent representation
adata_query.obsm["X_scVI"] = model_query.get_latent_representation()

print("Query model trained")
```

**Why `weight_decay=0.0`?**
- Ensures reference cell representations stay identical
- Only query cells get new embeddings
- Maintains consistency for combined analysis

---

### Step 7: Transfer Cell Type Labels

```python
# Predict cell types for query using reference classifier
latent_query = model_query.get_latent_representation()
predictions = clf.predict(latent_query)
probabilities = clf.predict_proba(latent_query)

# Store predictions
adata_query.obs["predicted_cell_type"] = predictions
adata_query.obs["prediction_confidence"] = probabilities.max(axis=1)

# Summarize predictions
print("Predicted cell types:")
print(adata_query.obs["predicted_cell_type"].value_counts())

# Visualize confidence distribution
plt.figure(figsize=(8, 4))
plt.hist(adata_query.obs["prediction_confidence"], bins=50)
plt.xlabel("Prediction Confidence")
plt.ylabel("Cells")
plt.title("Label Transfer Confidence")
plt.axvline(0.7, color='r', linestyle='--', label='High confidence threshold')
plt.legend()
plt.show()
```

---

### Step 8: Combined Visualization

```python
# Combine reference and query for visualization
adata_combined = anndata.concat(
    [adata_ref, adata_query],
    label="dataset",
    keys=["reference", "query"]
)

# Get latent for combined data through query model
# (This ensures reference cells have same embedding)
adata_combined.obsm["X_scVI"] = model_query.get_latent_representation(adata_combined)

# UMAP
sc.pp.neighbors(adata_combined, use_rep="X_scVI")
sc.tl.umap(adata_combined)

# Visualize
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

sc.pl.umap(adata_combined, color="dataset", ax=axes[0], show=False,
           title="Reference vs Query")

sc.pl.umap(adata_combined, color="cell_type", ax=axes[1], show=False,
           title="Cell Types (Reference)", na_color="lightgray")

# Show predictions for query
adata_combined.obs["display_cell_type"] = adata_combined.obs["cell_type"].copy()
query_mask = adata_combined.obs["dataset"] == "query"
adata_combined.obs.loc[query_mask, "display_cell_type"] = \
    adata_combined.obs.loc[query_mask, "predicted_cell_type"]

sc.pl.umap(adata_combined, color="display_cell_type", ax=axes[2], show=False,
           title="Cell Types (with predictions)")

plt.tight_layout()
plt.show()
```

---

### Step 9: Evaluate Label Transfer (If Ground Truth Available)

```python
# If query has ground truth labels
if "true_cell_type" in adata_query.obs.columns:
    # Calculate accuracy
    accuracy = accuracy_score(
        adata_query.obs["true_cell_type"],
        adata_query.obs["predicted_cell_type"]
    )
    print(f"Label transfer accuracy: {accuracy:.3f}")

    # Confusion matrix
    from sklearn.metrics import classification_report
    print("\nClassification Report:")
    print(classification_report(
        adata_query.obs["true_cell_type"],
        adata_query.obs["predicted_cell_type"]
    ))
```

---

### Step 10: Save Results

```python
# Save query data with predictions
adata_query.write_h5ad("query_mapped.h5ad")

# Save combined data
adata_combined.write_h5ad("reference_query_combined.h5ad")

# Save query model (for future incremental mapping)
model_query.save("query_model", overwrite=True)

print("Results saved")
```

---

## SCANVI for Semi-Supervised Label Transfer

If you have reference labels and want probabilistic predictions:

```python
# First train scVI reference (as above)
# Then convert to SCANVI

scanvi_ref = scvi.model.SCANVI.from_scvi_model(
    model_ref,
    unlabeled_category="Unknown",
    labels_key="cell_type"
)

# Train SCANVI reference
scanvi_ref.train(max_epochs=20, n_samples_per_label=100)
scanvi_ref.save("scanvi_reference")

# Map query with SCANVI
scvi.model.SCANVI.prepare_query_anndata(adata_query, "scanvi_reference")
scanvi_query = scvi.model.SCANVI.load_query_data(adata_query, "scanvi_reference")
scanvi_query.train(max_epochs=100, plan_kwargs={"weight_decay": 0.0})

# Get predictions (probabilistic)
predictions = scanvi_query.predict(adata_query)
adata_query.obs["scanvi_prediction"] = predictions
```

---

## Key Parameters Reference

### Reference Model Architecture (CRITICAL)

| Parameter | Required Value | Why |
|-----------|---------------|-----|
| `use_layer_norm` | "both" | Stabilizes across batches |
| `use_batch_norm` | "none" | Prevents batch-specific norm |
| `encode_covariates` | True | Enables new batch encoding |

### Training Parameters

| Parameter | Reference | Query | Notes |
|-----------|-----------|-------|-------|
| `max_epochs` | 400 | 100-200 | Query needs less |
| `weight_decay` | default | 0.0 | Preserves ref space |
| `early_stopping` | True | Optional | |

---

## Parameter Tuning Guide

### When to Adjust Parameters

**Before running, ask the user:**
1. Do you have an existing reference model to use?
2. Is your reference model scArches-compatible?
3. How different is query from reference (same tissue, different platform)?
4. Do you have ground truth labels for validation?

### Reference Model Settings

| Scenario | n_latent | n_layers | Notes |
|----------|----------|----------|-------|
| Small reference (<10k cells) | 20-30 | 2 | Standard |
| Large atlas (>100k cells) | 30-50 | 2-3 | More capacity |
| Many cell types (>20) | 30-40 | 2 | Capture diversity |

### Query Training

| Scenario | max_epochs | Notes |
|----------|------------|-------|
| Similar to reference | 100 | Fast convergence |
| Different platform | 200 | More adaptation |
| Very different | 300-400 | Careful: may diverge |

### Classifier Choice

| Scenario | Classifier | Notes |
|----------|------------|-------|
| Standard | RandomForest | Fast, robust |
| Many classes | RandomForest(n_estimators=200) | More trees |
| Probabilistic needed | SCANVI | Built-in uncertainty |
| Complex boundaries | MLP/Neural | sklearn MLPClassifier |

---

## Adaptation Prompts for Claude

When a user invokes this skill, consider asking:

1. **Reference status:**
   - "Do you have an existing reference model?"
   - "Was it trained with scArches-compatible settings?"
   - "What cell type annotations are available?"

2. **Query characteristics:**
   - "How many cells in your query?"
   - "Same tissue/organism as reference?"
   - "Same sequencing platform?"

3. **Goals:**
   - "Do you need cell type predictions?"
   - "Will you add more queries later?"
   - "Do you have ground truth for validation?"

### Check Reference Compatibility

```python
def check_scarches_compatibility(model_path):
    """Check if saved model is scArches-compatible."""
    import json
    import os

    config_path = os.path.join(model_path, "model.pt")
    # Note: Would need to load and inspect model config

    print("To verify scArches compatibility, check:")
    print("  1. encode_covariates=True")
    print("  2. use_layer_norm='both'")
    print("  3. use_batch_norm='none'")
    print("\nIf not compatible, need to retrain reference.")
```

---

## Troubleshooting

### Common Issues

1. **"encode_covariates must be True" error**:
   - Reference model wasn't trained correctly
   - Must retrain with `encode_covariates=True`

2. **Poor label transfer accuracy**:
   - Query very different from reference
   - Missing cell types in reference
   - Try longer training, different classifier

3. **Query cells cluster separately**:
   - May indicate batch effects
   - Check `weight_decay=0.0` was used
   - Verify gene alignment worked

4. **Reference latent space changed after query**:
   - Forgot `weight_decay=0.0`
   - Retrain query with correct setting

5. **Gene mismatch errors**:
   - Run `prepare_query_anndata` first
   - Check gene names match (symbols vs IDs)

### Verifying Reference Preservation

```python
# After query training, verify reference embeddings unchanged
latent_ref_original = model_ref.get_latent_representation(adata_ref)
latent_ref_through_query = model_query.get_latent_representation(adata_ref)

# Should be nearly identical
diff = np.abs(latent_ref_original - latent_ref_through_query).max()
print(f"Max difference in reference latent: {diff:.6f}")
if diff < 1e-5:
    print("Reference space preserved correctly!")
else:
    print("WARNING: Reference space may have changed")
```

---

## Advanced Usage

### Incremental Atlas Building

```python
# Add multiple queries incrementally
queries = ["query1.h5ad", "query2.h5ad", "query3.h5ad"]

current_model_path = "reference_model"

for i, query_path in enumerate(queries):
    # Load and prepare query
    adata_q = sc.read_h5ad(query_path)
    scvi.model.SCVI.prepare_query_anndata(adata_q, current_model_path)

    # Map query
    model_q = scvi.model.SCVI.load_query_data(adata_q, current_model_path)
    model_q.train(max_epochs=200, plan_kwargs={"weight_decay": 0.0})

    # Save for next iteration
    current_model_path = f"model_with_query_{i+1}"
    model_q.save(current_model_path, overwrite=True)

    print(f"Added query {i+1}")
```

### TotalVI Reference Mapping (for CITE-Seq)

```python
# Reference with protein
scvi.model.TOTALVI.setup_anndata(
    adata_ref,
    batch_key="batch",
    protein_expression_obsm_key="protein_expression"
)

totalvi_ref = scvi.model.TOTALVI(
    adata_ref,
    use_layer_norm="both",
    use_batch_norm="none"
)
totalvi_ref.train()
totalvi_ref.save("totalvi_reference")

# Query mapping (even RNA-only!)
totalvi_query = scvi.model.TOTALVI.load_query_data(
    adata_query,
    "totalvi_reference"
)
totalvi_query.train(max_epochs=200, plan_kwargs={"weight_decay": 0.0})

# Impute proteins for query
_, protein_imputed = totalvi_query.get_normalized_expression(
    adata_query,
    n_samples=25,
    return_mean=True
)
```

---

## References

- [scArches Tutorial](https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/multimodal/scarches_scvi_tools.html)
- [scArches Paper](https://www.nature.com/articles/s41587-021-01001-7) - Lotfollahi et al., Nature Biotechnology 2022
- [scvi-tools Documentation](https://docs.scvi-tools.org/)
