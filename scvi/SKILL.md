---
name: scvi
description: Deep learning for single-cell analysis using scvi-tools. This skill should be used when users need (1) batch integration/correction with scVI/scANVI, (2) label transfer with scANVI or seed labeling, (3) reference mapping with scArches, (4) multi-sample analysis with MrVI, (5) CITE-seq/multimodal analysis with TotalVI/MultiVI, (6) ATAC-seq analysis with PeakVI, (7) spatial transcriptomics deconvolution with DestVI/Stereoscope/Cell2location/Tangram, (8) spatial denoising with ResolVI, (9) microenvironment modeling with scVIVA, or (10) any deep learning-based single-cell method. Triggers include mentions of scVI, scANVI, TotalVI, PeakVI, MultiVI, DestVI, MrVI, Stereoscope, Cell2location, Tangram, ResolVI, scVIVA, scArches, variational autoencoder, VAE, batch correction, data integration, multi-modal, CITE-seq, multiome, spatial deconvolution, reference mapping, label transfer, latent space.
---

# scvi-tools Deep Learning Skill

This skill provides guidance for deep learning-based single-cell analysis using scvi-tools, the leading framework for probabilistic models in single-cell genomics.

## How to Use This Skill

1. Identify the appropriate workflow from the model selection table or decision trees below
2. Read the corresponding reference file in `references/workflows/` for detailed steps
3. Use scripts in `scripts/` to avoid rewriting common code
4. For installation or GPU issues, consult `references/common/environment_setup.md`
5. For debugging, consult `references/common/troubleshooting.md`

## Model Selection Guide

| Data Type | Model | Reference File |
|-----------|-------|----------------|
| scRNA-seq | **scVI/scANVI** | `references/workflows/scvi_integration.md` |
| scRNA-seq + labels | **scANVI** | `references/workflows/scanvi_label_transfer.md` |
| scRNA-seq + marker signatures | **scANVI seed** | `references/workflows/scanvi_seed_labeling.md` |
| Query to reference | **scArches** | `references/workflows/scarches_reference_mapping.md` |
| Multi-sample comparison | **MrVI** | `references/workflows/mrvi_multisample.md` |
| CITE-seq (RNA+protein) | **TotalVI** | `references/workflows/totalvi_citeseq.md` |
| CITE-seq reference mapping | **TotalVI** | `references/workflows/totalvi_reference_mapping.md` |
| CITE-seq + RNA-only | **TotalVI** | `references/workflows/citeseq_rna_integration.md` |
| Multiome (RNA+ATAC) | **MultiVI** | `references/workflows/multivi_multimodal.md` |
| scATAC-seq | **PeakVI** | `references/workflows/peakvi_accessibility.md` |
| Spatial + scRNA ref | **DestVI** | `references/workflows/destvi_deconvolution.md` |
| Spatial proportions | **Stereoscope** | `references/workflows/stereoscope_deconvolution.md` |
| Spatial cell counts | **Cell2location** | `references/workflows/cell2location_mapping.md` |
| Spatial gene imputation | **Tangram** | `references/workflows/tangram_mapping.md` |
| Cellular spatial denoising | **ResolVI** | `references/workflows/resolvi_denoising.md` |
| Spatial microenvironment | **scVIVA** | `references/workflows/scviva_environment.md` |

## Quick Decision Tree

```
Need to integrate scRNA-seq data?
├── Have cell type labels? → scANVI (scanvi_label_transfer.md)
├── Have marker gene signatures? → scANVI seed labeling (scanvi_seed_labeling.md)
└── No labels? → scVI (scvi_integration.md)

Need reference mapping?
├── Map query to scVI/scANVI reference? → scArches (scarches_reference_mapping.md)
└── Map query to CITE-seq reference? → TotalVI reference (totalvi_reference_mapping.md)

Have multi-sample data?
└── Compare samples, differential abundance? → MrVI (mrvi_multisample.md)

Have multi-modal data? → See Multimodal section below

Have spatial data? → See Spatial section below
```

---

## Multimodal CITE-Seq Analysis

**Decision Tree:**

```
What type of multimodal data do you have?

CITE-Seq (RNA + Protein):
├── All cells have protein measurements?
│   ├── YES → Do you have a reference atlas?
│   │   ├── YES → TotalVI Reference Mapping (totalvi_reference_mapping.md)
│   │   └── NO → TotalVI (totalvi_citeseq.md)
│   └── NO (some cells RNA-only) → CITE-RNA Integration (citeseq_rna_integration.md)

Multiome (RNA + ATAC):
├── All cells paired? → MultiVI or separate scVI + PeakVI
└── Mix of paired/unpaired → MultiVI (multivi_multimodal.md)

Triple Modality (RNA + Protein + ATAC):
└── MultiVI (multivi_multimodal.md)
```

**Quick Reference:**

| Situation | Method | Key Output |
|-----------|--------|------------|
| Standard CITE-Seq | TotalVI | Joint embedding, denoised protein |
| CITE-Seq + RNA-only cells | CITE-RNA Integration | Imputed proteins |
| Map to CITE-Seq atlas | TotalVI Reference | Transferred labels, imputed proteins |
| Multiome (paired/unpaired) | MultiVI | Unified embedding, imputed modalities |

---

## Spatial Transcriptomics Analysis

**Decision Tree:**

```
What is your spatial data resolution?

CELLULAR-RESOLUTION (Xenium, MERFISH, CosMx):
├── Is there segmentation noise or background issues?
│   ├── YES → ResolVI (resolvi_denoising.md) → then continue
│   └── NO → Continue
└── Do you want to understand microenvironment effects?
    ├── YES → scVIVA (scviva_environment.md)
    └── NO → Standard scanpy/scVI workflow

SPOT-BASED (Visium, Slide-seq):
├── Do you have a matched scRNA-seq reference?
│   ├── YES → What's your primary goal?
│   │   ├── Cell type proportions only → Stereoscope (stereoscope_deconvolution.md)
│   │   ├── Proportions + cell state variation → DestVI (destvi_deconvolution.md)
│   │   ├── Absolute cell abundances → Cell2location (cell2location_mapping.md)
│   │   └── Gene imputation + mapping → Tangram (tangram_mapping.md)
│   └── NO → Limited options (clustering only, or find public atlas)
```

**Quick Reference:**

| Situation | Method | Key Output |
|-----------|--------|------------|
| Xenium/MERFISH noise | ResolVI | Denoised expression |
| Niche/microenvironment | scVIVA | Niche-aware embedding |
| Visium proportions | Stereoscope | Cell type proportions |
| Visium proportions + states | DestVI | Proportions + gamma (cell states) |
| Visium cell counts | Cell2location | Absolute abundances |
| Visium gene imputation | Tangram | Imputed genes |

**Comparing Deconvolution Methods:**

Use `scripts/compare_deconvolution_methods.py` to run multiple methods and compare:
```bash
python scripts/compare_deconvolution_methods.py \
    --spatial visium.h5ad \
    --reference scrna_ref.h5ad \
    --labels-key cell_type
```

---

## Critical Requirements

1. **Raw counts required**: scvi-tools models require integer count data
   ```python
   adata.layers["counts"] = adata.X.copy()  # Before normalization
   scvi.model.SCVI.setup_anndata(adata, layer="counts")
   ```

2. **HVG selection**: Use 2000-4000 highly variable genes
   ```python
   sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="batch", layer="counts", flavor="seurat_v3")
   adata = adata[:, adata.var['highly_variable']].copy()
   ```

3. **Batch information**: Specify batch_key for integration
   ```python
   scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
   ```

---

## CLI Scripts

Model-specific scripts for common workflows:

| Script | Purpose |
|--------|---------|
| `integration_analysis.py` | scVI/scANVI batch integration |
| `label_transfer_analysis.py` | scANVI label transfer |
| `seed_labeling_analysis.py` | scANVI semi-supervised with markers |
| `scarches_analysis.py` | scArches reference mapping |
| `mrvi_analysis.py` | MrVI multi-sample analysis |
| `totalvi_analysis.py` | TotalVI CITE-seq analysis |
| `peakvi_analysis.py` | PeakVI ATAC-seq analysis |
| `destvi_analysis.py` | DestVI spatial deconvolution |
| `compare_deconvolution_methods.py` | Compare spatial methods |
| `validate_adata.py` | Check data compatibility |

**Quick Examples:**
```bash
# Validate data before running any model
python scripts/validate_adata.py data.h5ad

# Basic scVI integration
python scripts/integration_analysis.py data.h5ad --batch-key batch

# scANVI with cell type labels
python scripts/integration_analysis.py data.h5ad --batch-key batch --labels-key cell_type

# Label transfer to query data
python scripts/label_transfer_analysis.py query.h5ad --reference-model scanvi_model/

# Spatial deconvolution
python scripts/destvi_analysis.py spatial.h5ad --reference scrna_ref.h5ad --labels-key cell_type
```

**Shared utilities:** `model_utils.py` provides importable functions for custom workflows.

---

## Environment & Troubleshooting

- **Installation, GPU setup**: `references/common/environment_setup.md`
- **Data preparation**: `references/common/data_preparation.md`
- **Common issues**: `references/common/troubleshooting.md`

---

## Key Resources

- [scvi-tools Documentation](https://docs.scvi-tools.org/)
- [scvi-tools Tutorials](https://docs.scvi-tools.org/en/stable/tutorials/index.html)
- [Model Hub](https://huggingface.co/scvi-tools)
- [GitHub Issues](https://github.com/scverse/scvi-tools/issues)
