# Spatial Transcriptomics Analysis Pipelines

This reference provides detailed multi-stage pipelines, method chaining guidance, and scenario-based recommendations for spatial transcriptomics analysis.

For method selection decision trees, see the main SKILL.md.

---

## Data Assessment

### Platform Identification

**Answer these questions about your spatial data:**

| Question | Options |
|----------|---------|
| What platform was used? | Visium, Xenium, MERFISH, CosMx, Slide-seq, STARmap, other |
| What is the resolution? | Spot-based (~50-100μm) or Cellular (~single cell) |
| How many genes measured? | Panel (<500), Medium (500-5000), Whole transcriptome (>5000) |
| Do you have matched scRNA-seq? | Yes / No |

### Resolution Classification

```
SPOT-BASED (Multi-cell per spot):
├── 10x Visium (55μm spots, ~1-10 cells/spot)
├── Slide-seq/Slide-seqV2 (10μm beads)
├── HDST (2μm pixels, binned)
└── Stereo-seq (binned to spots)

CELLULAR-RESOLUTION (Single-cell):
├── 10x Xenium (subcellular)
├── MERFISH (subcellular)
├── CosMx (subcellular)
├── seqFISH+ (subcellular)
└── STARmap (subcellular)
```

---

## Analysis Goals Checklist

Check all that apply to your project:

**Data Quality & Preprocessing**
- [ ] Correct segmentation errors in cellular data
- [ ] Remove background noise/contamination
- [ ] Denoise expression measurements

**Cell Type Analysis**
- [ ] Identify what cell types are present
- [ ] Estimate cell type proportions per spot
- [ ] Get absolute cell counts per location
- [ ] Transfer annotations from reference data

**Spatial Biology**
- [ ] Understand how location affects gene expression
- [ ] Identify tissue niches/microenvironments
- [ ] Find spatially variable genes
- [ ] Analyze cell-cell interactions

**Integration**
- [ ] Map scRNA-seq data to spatial locations
- [ ] Impute genes not measured spatially
- [ ] Compare multiple spatial samples

---

## Common Analysis Pipelines

### Pipeline A: Cellular-Resolution Analysis (Xenium/MERFISH/CosMx)

**Goal**: Clean data, identify niches, understand microenvironment effects

```
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 1: Data Cleaning                                         │
│  Workflow: resolvi_denoising.md                                 │
│  ├── Correct segmentation errors                                │
│  ├── Remove background signal                                   │
│  └── Output: Denoised expression matrix                        │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 2: Cell Type Assignment                                  │
│  Options:                                                       │
│  ├── ResolVI semi-supervised (if labels available)             │
│  ├── Standard clustering + marker genes                        │
│  └── scANVI label transfer (if reference available)            │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 3: Niche Analysis                                        │
│  Workflow: scviva_environment.md                                │
│  ├── Model neighborhood effects                                 │
│  ├── Identify tissue niches                                     │
│  └── Niche-aware differential expression                       │
└─────────────────────────────────────────────────────────────────┘
```

**Code outline:**
```python
# Stage 1: ResolVI denoising
scvi.external.RESOLVI.setup_anndata(adata, layer="counts", labels_key="cell_type")
resolvi_model = scvi.external.RESOLVI(adata, semisupervised=True)
resolvi_model.train(max_epochs=100)
adata.layers["denoised"] = resolvi_model.get_corrected_expression()

# Stage 2: Use ResolVI predictions or refine with clustering
adata.obs["cell_type_refined"] = resolvi_model.predict(adata, soft=False)

# Stage 3: scVIVA niche analysis
scvi.external.SCVIVA.preprocessing_anndata(adata, k_nn=20, ...)
scvi.external.SCVIVA.setup_anndata(adata, layer="counts", ...)
scviva_model = scvi.external.SCVIVA(adata)
scviva_model.train(max_epochs=600)
adata.obsm["X_scVIVA"] = scviva_model.get_latent_representation()
```

---

### Pipeline B: Comprehensive Visium Deconvolution

**Goal**: Get cell type composition with multiple methods for validation

```
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 1: Reference Preparation                                 │
│  ├── Load scRNA-seq reference                                   │
│  ├── QC and filter                                              │
│  ├── Ensure cell type annotations                               │
│  └── Select shared genes with spatial data                     │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 2: Run Multiple Deconvolution Methods                    │
│  (Can run in parallel)                                          │
│  ├── Stereoscope: stereoscope_deconvolution.md                 │
│  ├── DestVI: destvi_deconvolution.md                           │
│  └── Cell2location: cell2location_mapping.md                   │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 3: Compare and Validate                                  │
│  ├── Correlate proportions across methods                       │
│  ├── Check against known tissue structure                       │
│  ├── Identify consensus cell types                              │
│  └── Report confidence based on agreement                      │
└─────────────────────────────────────────────────────────────────┘
```

**Code outline:**
```python
# Run all three methods (can parallelize)
# Method 1: Stereoscope
RNAStereoscope.setup_anndata(sc_adata, layer="counts", labels_key="cell_type")
stereo_sc = RNAStereoscope(sc_adata)
stereo_sc.train(max_epochs=100)
SpatialStereoscope.setup_anndata(st_adata, layer="counts")
stereo_st = SpatialStereoscope.from_rna_model(st_adata, stereo_sc)
stereo_st.train(max_epochs=2000)
st_adata.obsm["stereoscope"] = stereo_st.get_proportions()

# Method 2: DestVI
CondSCVI.setup_anndata(sc_adata, layer="counts", labels_key="cell_type")
destvi_sc = CondSCVI(sc_adata)
destvi_sc.train(max_epochs=300)
DestVI.setup_anndata(st_adata, layer="counts")
destvi_st = DestVI.from_rna_model(st_adata, destvi_sc)
destvi_st.train(max_epochs=2500)
st_adata.obsm["destvi"] = destvi_st.get_proportions()

# Method 3: Cell2location
RegressionModel.setup_anndata(sc_adata, labels_key="cell_type")
c2l_ref = RegressionModel(sc_adata)
c2l_ref.train(max_epochs=250)
# ... continue with Cell2location spatial model

# Compare results
import pandas as pd
comparison = pd.DataFrame({
    'stereoscope': st_adata.obsm["stereoscope"]["T_cells"],
    'destvi': st_adata.obsm["destvi"]["T_cells"],
    'cell2location': st_adata.obsm["cell2location"]["T_cells"]
})
print(comparison.corr())  # Check correlation between methods
```

---

### Pipeline C: Visium + Gene Imputation

**Goal**: Map cell types AND impute genes not measured in spatial data

```
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 1: Tangram Mapping                                       │
│  Workflow: tangram_mapping.md                                   │
│  ├── Map scRNA-seq cells to spatial spots                      │
│  ├── Project cell type annotations                              │
│  └── Impute unmeasured genes                                   │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 2: Deconvolution on Original Data                        │
│  Workflow: destvi_deconvolution.md                              │
│  ├── Get refined proportions                                    │
│  └── Extract cell state variation (gamma)                      │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 3: Analysis with Imputed Genes                           │
│  ├── Spatial differential expression                            │
│  ├── Pathway analysis using imputed genes                       │
│  └── Cell-cell communication (e.g., CellChat, LIANA)           │
└─────────────────────────────────────────────────────────────────┘
```

---

### Pipeline D: Multi-Sample Spatial Analysis

**Goal**: Compare spatial patterns across conditions/samples

```
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 1: Process Each Sample                                   │
│  ├── QC and preprocessing per sample                            │
│  ├── Deconvolution per sample (same reference)                  │
│  └── Store results with sample labels                          │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 2: Integration                                           │
│  ├── Concatenate samples with batch annotation                  │
│  ├── Batch-correct if needed (Harmony, scVI)                    │
│  └── Joint embedding for comparison                            │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  STAGE 3: Differential Analysis                                 │
│  ├── Compare cell type proportions across conditions            │
│  ├── Identify condition-specific spatial patterns               │
│  └── Statistical testing (MiloR, etc.)                         │
└─────────────────────────────────────────────────────────────────┘
```

---

## Method Chaining Guide

### When to Chain Methods

| First Method | Then Use | Why |
|--------------|----------|-----|
| ResolVI | scVIVA | Clean data before modeling environment |
| ResolVI | Any clustering | Denoised data gives better clusters |
| Tangram | DestVI | Impute genes, then get refined proportions |
| Tangram | Cell-cell communication | Need ligand/receptor genes imputed |
| DestVI | Downstream on gamma | Cell states for trajectory/pseudotime |
| Any deconvolution | Squidpy | Spatial statistics on cell type maps |

### Incompatible Combinations

| Don't Do This | Why |
|---------------|-----|
| ResolVI on Visium | ResolVI is for cellular-resolution only |
| scVIVA on Visium | scVIVA needs cellular coordinates |
| DestVI without reference | Requires scRNA-seq reference |
| Cell2location on panels | Needs sufficient gene overlap |

---

## Quick Start by Scenario

### Scenario 1: "I have Visium data and a matched scRNA-seq reference"

**Recommended workflow:**
1. Preprocess both datasets (QC, filter, normalize for visualization)
2. Choose deconvolution based on needs:
   - **Just proportions?** → Stereoscope (fastest)
   - **Need cell states?** → DestVI (get gamma values)
   - **Need cell counts?** → Cell2location (absolute abundances)
3. Optionally run Tangram for gene imputation

**Start with:** `destvi_deconvolution.md` (most versatile)

---

### Scenario 2: "I have Xenium data with some noise issues"

**Recommended workflow:**
1. Run ResolVI to denoise and correct segmentation
2. Extract corrected expression for downstream use
3. Run scVIVA for niche analysis
4. Perform niche-aware differential expression

**Start with:** `resolvi_denoising.md`

---

### Scenario 3: "I want to compare multiple deconvolution methods"

**Recommended workflow:**
1. Prepare reference once
2. Run Stereoscope, DestVI, and Cell2location in parallel
3. Compare results (correlation, spatial patterns)
4. Report consensus findings

**Use:** `scripts/compare_deconvolution_methods.py`

---

### Scenario 4: "I need to impute genes for cell-cell communication analysis"

**Recommended workflow:**
1. Run Tangram to map scRNA-seq to spatial
2. Project all genes (including ligands/receptors)
3. Use imputed data for CellChat/LIANA/CellPhoneDB

**Start with:** `tangram_mapping.md`

---

## Visualization Utilities

Use `scripts/spatial_viz_utils.py` for publication-quality figures:

```python
from spatial_viz_utils import (
    plot_proportions_grid,
    plot_dominant_celltype,
    plot_celltype_comparison,
    plot_niche_composition,
    plot_spatial_correlation,
    create_summary_figure,
)

# Grid of all cell type proportions
plot_proportions_grid(adata, proportions_key="destvi", ncols=4, save_path="proportions.png")

# Dominant cell type per spot
plot_dominant_celltype(adata, threshold=0.15, save_path="dominant.png")

# Compare methods for specific cell type
plot_celltype_comparison(
    adata,
    cell_type="T_cells",
    method_results={"DestVI": destvi_props, "Stereoscope": stereo_props},
    save_path="tcell_comparison.png"
)

# Niche composition bar chart
plot_niche_composition(adata, cluster_key="leiden", save_path="niche_composition.png")

# Cell type co-localization heatmap
plot_spatial_correlation(adata, proportions_key="destvi", save_path="colocalization.png")

# Comprehensive summary figure
create_summary_figure(adata, cluster_key="leiden", save_path="summary.png")
```

---

## Troubleshooting

### "My deconvolution results look wrong"

```
Check these in order:
1. Is your reference from the same tissue/species?
   └── NO → Find better reference or use public atlas
2. Are enough genes shared between datasets?
   └── <1000 genes → Filter less aggressively
3. Did the model converge?
   └── NO → Train longer, check loss curves
4. Are cell types in reference complete?
   └── NO → Results will miss those types
```

### "I don't have a scRNA-seq reference"

```
Options:
1. Use public atlases (Human Cell Atlas, Tabula Sapiens, etc.)
2. Use Tangram with existing annotated datasets
3. For Xenium: Use panel-specific references from 10x
4. Cluster and annotate manually (no deconvolution)
```

### "My data is too large for memory"

```
Solutions:
1. Process tissue regions separately
2. Reduce gene count (use HVGs)
3. Use GPU acceleration
4. Subsample for initial exploration
```
