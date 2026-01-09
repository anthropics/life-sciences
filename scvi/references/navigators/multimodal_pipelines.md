# Multimodal CITE-Seq Analysis Pipelines

This reference provides detailed multi-stage pipelines, method comparison guidance, and scenario-based recommendations for multimodal single-cell analysis.

For method selection decision trees, see the main SKILL.md.

---

## Data Assessment

### Data Inventory Questions

**Answer these questions about your multimodal data:**

| Question | Options |
|----------|---------|
| What modalities do you have? | RNA only / RNA + Protein / RNA + ATAC / All three |
| Is protein data available for all cells? | Yes / Partial / No (want to impute) |
| Do you have a reference atlas? | Yes / No / Want to create one |
| Are your datasets paired or unpaired? | All paired / Mix / All unpaired |
| Do you need to map new data to existing model? | Yes / No |

### Data Structure Classification

```
MULTIMODAL DATA TYPES:
├── CITE-Seq (RNA + Protein)
│   ├── All cells have both modalities → TotalVI
│   ├── Some cells missing protein → CITE-RNA Integration
│   └── Map to reference → TotalVI Reference Mapping
│
├── Multiome (RNA + ATAC)
│   ├── All cells paired → scVI + scATAC separately, or MultiVI
│   └── Mix of paired/unpaired → MultiVI
│
├── Triple Modality (RNA + Protein + ATAC)
│   └── Use MultiVI
│
└── REFERENCE MAPPING SCENARIOS
    ├── Map RNA-only to CITE-Seq reference → TotalVI Reference
    ├── Map any query to scVI reference → scArches
    └── Transfer cell type labels → scArches + SCANVI
```

---

## Analysis Goals Checklist

**Core Analysis**
- [ ] Joint dimensionality reduction of RNA + protein
- [ ] Cluster cells using both modalities
- [ ] Identify cell types from protein markers
- [ ] Denoise protein measurements

**Protein Imputation**
- [ ] Impute proteins for RNA-only cells
- [ ] Transfer protein predictions across batches
- [ ] Fill in missing modalities

**Reference Mapping**
- [ ] Map new data to existing reference atlas
- [ ] Transfer cell type annotations
- [ ] Project onto reference UMAP
- [ ] Preserve reference latent space

**Integration**
- [ ] Integrate CITE-Seq with scRNA-seq
- [ ] Combine datasets with different protein panels
- [ ] Handle batch effects across modalities

---

## Common Analysis Pipelines

### Pipeline A: Standard CITE-Seq Analysis

**Goal**: Joint analysis of RNA and protein from CITE-Seq experiment

```
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 1: Data Preparation                                       │
│  ├── Load MuData with RNA and protein modalities                 │
│  ├── Convert sparse to dense matrices                            │
│  ├── Store raw counts in layers                                  │
│  └── Subset to shared proteins across batches                    │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 2: TotalVI Training                                       │
│  Workflow: totalvi_citeseq.md                                    │
│  ├── Setup MuData with batch_key                                 │
│  ├── Train model (early_stopping=True)                          │
│  └── Extract latent representation                              │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 3: Downstream Analysis                                    │
│  ├── Clustering on joint latent space                           │
│  ├── Get denoised RNA and protein                               │
│  ├── Differential expression (RNA + protein together)           │
│  └── Visualize protein markers on UMAP                          │
└─────────────────────────────────────────────────────────────────┘
```

**Code outline:**
```python
import scvi
import muon
import scanpy as sc

# Load CITE-Seq data
mdata = muon.read_h5mu("citeseq_data.h5mu")

# Setup TotalVI
scvi.model.TOTALVI.setup_mudata(
    mdata,
    rna_layer="counts",
    protein_layer=None,
    batch_key="batch",
    modalities={"rna_layer": "rna", "protein_layer": "prot", "batch_key": "rna"}
)

# Train
model = scvi.model.TOTALVI(mdata)
model.train(early_stopping=True)

# Extract outputs
rna = mdata.mod["rna"]
rna.obsm["X_totalVI"] = model.get_latent_representation()
rna_denoised, protein_denoised = model.get_normalized_expression(n_samples=25, return_mean=True)

# Cluster
sc.pp.neighbors(rna, use_rep="X_totalVI")
sc.tl.umap(rna)
sc.tl.leiden(rna, key_added="leiden_totalVI")

# Differential expression
de_results = model.differential_expression(
    groupby="leiden_totalVI",
    delta=0.5,
    batch_correction=True
)
```

---

### Pipeline B: Reference Mapping with Label Transfer

**Goal**: Map new query data to existing reference, transfer cell type labels

```
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 1: Prepare Reference Model                                │
│  Workflow: scarches_reference_mapping.md                         │
│  ├── Setup with encode_covariates=True (critical!)              │
│  ├── use_layer_norm="both", use_batch_norm="none"               │
│  ├── Train reference model                                       │
│  └── Train classifier on reference labels                       │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 2: Map Query Data                                         │
│  ├── prepare_query_anndata() to align genes                     │
│  ├── load_query_data() with reference model                     │
│  ├── Train query with weight_decay=0.0 (preserves ref space)    │
│  └── Extract query latent representation                        │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 3: Transfer Annotations                                   │
│  ├── Apply reference classifier to query latent                 │
│  ├── Visualize query on reference UMAP                          │
│  └── Assess label transfer accuracy                             │
└─────────────────────────────────────────────────────────────────┘
```

**Code outline:**
```python
import scvi
from sklearn.ensemble import RandomForestClassifier

# Reference model (with scArches-compatible settings)
scvi.model.SCVI.setup_anndata(adata_ref, batch_key="batch", layer="counts")
model_ref = scvi.model.SCVI(
    adata_ref,
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,  # Critical for scArches!
    n_layers=2
)
model_ref.train()
model_ref.save("reference_model")

# Train classifier
latent_ref = model_ref.get_latent_representation()
clf = RandomForestClassifier()
clf.fit(latent_ref, adata_ref.obs["cell_type"])

# Map query
scvi.model.SCVI.prepare_query_anndata(adata_query, "reference_model")
model_query = scvi.model.SCVI.load_query_data(adata_query, "reference_model")
model_query.train(max_epochs=200, plan_kwargs={"weight_decay": 0.0})

# Transfer labels
latent_query = model_query.get_latent_representation()
adata_query.obs["predicted_cell_type"] = clf.predict(latent_query)
```

---

### Pipeline C: Integrate CITE-Seq with RNA-Only Data

**Goal**: Combine CITE-Seq dataset with RNA-only dataset, impute proteins

```
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 1: Prepare Combined Dataset                               │
│  ├── Merge CITE-Seq and RNA-only AnnData                        │
│  ├── Set protein values to 0 for RNA-only cells                 │
│  │   (This signals "missing" to TotalVI)                        │
│  └── Select shared highly variable genes                        │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 2: TotalVI Integration                                    │
│  Workflow: citeseq_rna_integration.md                            │
│  ├── Setup with batch_key and protein_expression_obsm_key       │
│  ├── TotalVI automatically detects missing proteins             │
│  └── Train joint model                                          │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 3: Protein Imputation                                     │
│  ├── get_normalized_expression with transform_batch             │
│  ├── Impute proteins for RNA-only cells                         │
│  └── Validate imputation quality                                │
└─────────────────────────────────────────────────────────────────┘
```

**Code outline:**
```python
import scvi
import numpy as np

# Set protein to 0 for RNA-only batch (signals "missing")
batch = adata.obs["batch"].values
adata.obsm["protein_expression"].loc[batch == "rna_only_batch"] = 0

# Setup and train
scvi.model.TOTALVI.setup_anndata(
    adata,
    batch_key="batch",
    protein_expression_obsm_key="protein_expression"
)
model = scvi.model.TOTALVI(adata)
model.train()

# Impute proteins using counterfactual prediction
_, protein_imputed = model.get_normalized_expression(
    transform_batch="cite_seq_batch",  # "What if from CITE batch?"
    n_samples=25,
    return_mean=True
)

# Store imputed values
adata.obsm["protein_imputed"] = protein_imputed
```

---

### Pipeline D: Multimodal Integration with MultiVI

**Goal**: Integrate datasets with different modality combinations

```
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 1: Prepare Multimodal Data                                │
│  ├── Ensure shared peaks/genes across datasets                   │
│  ├── Order features: genes before peaks                          │
│  ├── Format as MuData                                            │
│  └── Label modality batches appropriately                       │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 2: MultiVI Training                                       │
│  Workflow: multivi_multimodal.md                                 │
│  ├── Setup with n_genes and n_regions                           │
│  ├── Specify modality annotations                                │
│  └── Train model                                                │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 3: Joint Analysis                                         │
│  ├── Extract unified latent representation                      │
│  ├── Cluster on joint embedding                                 │
│  ├── Impute missing modalities                                  │
│  └── Cross-modal analysis                                       │
└─────────────────────────────────────────────────────────────────┘
```

**Code outline:**
```python
import scvi

# Setup MultiVI
scvi.model.MULTIVI.setup_mudata(
    mdata,
    modalities={
        "rna_layer": "rna",
        "atac_layer": "atac",
    }
)

model = scvi.model.MULTIVI(
    mdata,
    n_genes=len(mdata.mod["rna"].var),
    n_regions=len(mdata.mod["atac"].var)
)
model.train()

# Joint embedding
mdata.obsm["X_multivi"] = model.get_latent_representation()
sc.pp.neighbors(mdata, use_rep="X_multivi")
sc.tl.umap(mdata, min_dist=0.2)

# Impute missing modalities
rna_imputed = model.get_normalized_expression()
atac_imputed = model.get_accessibility_estimates()
```

---

## Method Comparison Guide

### When to Use Each Method

| Feature | TotalVI | scArches | CITE-RNA Int. | MultiVI |
|---------|---------|----------|---------------|---------|
| **Primary data** | CITE-Seq | Any + query | CITE + RNA-only | Mixed modality |
| **Protein imputation** | Yes | Via TotalVI | Yes | No |
| **Reference mapping** | Via separate skill | Yes | No | No |
| **Label transfer** | Via classifier | Yes (SCANVI) | Via classifier | Via clustering |
| **Unpaired data** | No | No | Yes | Yes |
| **ATAC support** | No | scATAC models | No | Yes |

### Method Chaining Guide

| First Method | Then Use | Why |
|--------------|----------|-----|
| TotalVI | scArches | Create reference, then map queries |
| TotalVI | Classifier | Train cell type predictor on latent |
| CITE-RNA Integration | DE analysis | Compare imputed vs measured |
| MultiVI | Motif/peak analysis | Use imputed accessibility |

### Key Parameter Comparison

| Parameter | TotalVI | scArches | MultiVI |
|-----------|---------|----------|---------|
| Batch correction | `batch_key` | `batch_key` + architecture | Modality batch |
| Missing data | Zeros in protein | N/A | Mixed modality |
| Query training | N/A | `weight_decay=0.0` | N/A |
| Model surgery | N/A | `encode_covariates=True` | N/A |

---

## Quick Start by Scenario

### Scenario 1: "I have standard CITE-Seq data from one experiment"

**Recommended workflow:**
1. Load as MuData (RNA + protein modalities)
2. Run TotalVI for joint embedding
3. Get denoised values and cluster
4. Run differential expression

**Start with:** `totalvi_citeseq.md`

---

### Scenario 2: "I want to map my data to a public CITE-Seq atlas"

**Recommended workflow:**
1. Download/load reference TotalVI model
2. Prepare query with same genes
3. Use TotalVI reference mapping
4. Impute proteins if query is RNA-only

**Start with:** `totalvi_reference_mapping.md`

---

### Scenario 3: "I have CITE-Seq and want to add RNA-only samples"

**Recommended workflow:**
1. Combine datasets, set protein=0 for RNA-only cells
2. Train integrated TotalVI model
3. Impute proteins for RNA-only cells
4. Analyze unified embedding

**Start with:** `citeseq_rna_integration.md`

---

### Scenario 4: "I want to build a reference atlas and map future queries"

**Recommended workflow:**
1. Train reference with scArches-compatible settings
2. Save model with `encode_covariates=True`
3. Use scArches to map new batches
4. Optional: Train SCANVI for label transfer

**Start with:** `scarches_reference_mapping.md`

---

### Scenario 5: "I have multiome (RNA+ATAC) with some unpaired samples"

**Recommended workflow:**
1. Ensure shared peaks across datasets
2. Use MultiVI to leverage paired cells as anchors
3. Impute missing modalities
4. Unified analysis across modality types

**Start with:** `multivi_multimodal.md`

---

## Interactive Assessment

When helping users with multimodal analysis, guide them through these questions:

### Essential Questions to Ask

```
QUESTION 1: MODALITIES
"What modalities does your data include?"
- [ ] RNA + Protein (CITE-Seq)
- [ ] RNA + ATAC (Multiome)
- [ ] RNA + Protein + ATAC (Triple)
- [ ] RNA only (want to add to multimodal reference)

QUESTION 2: DATA COMPLETENESS
"Do all cells have all modalities measured?"
- [ ] Yes, all cells are fully paired
- [ ] No, some cells are missing modalities
- [ ] Mix of paired and unpaired cells

QUESTION 3: REFERENCE MAPPING
"Are you working with a reference atlas?"
- [ ] I have an existing reference to map to
- [ ] I want to create a reference for future mapping
- [ ] No reference mapping needed

QUESTION 4: PRIMARY GOAL
"What's your main analysis goal?"
- [ ] Standard joint analysis (clustering, DE)
- [ ] Impute missing proteins/modalities
- [ ] Transfer cell type labels from reference
- [ ] Integrate multiple datasets
```

### Parameter Adaptation by User Situation

| User Says | Adaptation |
|-----------|------------|
| "Small protein panel (<50 proteins)" | Standard TotalVI settings work well |
| "Large protein panel (>200 proteins)" | May need more training epochs |
| "Many batches with different proteins" | Subset to shared proteins |
| "RNA-only cells outnumber CITE-Seq" | Imputation quality may vary |
| "Want to preserve reference exactly" | Use `weight_decay=0.0` for query |
| "Query has different genes" | `prepare_query_anndata()` handles this |

---

## References

### Tutorials
- [TotalVI](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/multimodal/totalVI.html)
- [scArches](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/multimodal/scarches_scvi_tools.html)
- [TotalVI Reference Mapping](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/multimodal/totalVI_reference_mapping.html)
- [CITE-RNA Integration](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/multimodal/cite_scrna_integration_w_totalVI.html)
- [MultiVI](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/multimodal/MultiVI_tutorial.html)

### Papers
- TotalVI: Gayoso et al., Nature Methods (2021)
- scArches: Lotfollahi et al., Nature Biotechnology (2022)
- MultiVI: Ashuach et al., Nature Methods (2023)
