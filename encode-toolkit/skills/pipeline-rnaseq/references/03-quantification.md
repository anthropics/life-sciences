# Stage 3: Gene and Transcript Quantification

## Tools
- **RSEM 1.3.3**: Gene/transcript quantification (Li & Dewey 2011, ~6,000 citations)
- **Kallisto 0.50.1**: Fast pseudoalignment quantification (Bray et al. 2016, ~4,000 citations)

## RSEM Quantification (Primary)

RSEM uses an expectation-maximization (EM) algorithm to probabilistically assign
multi-mapped reads to genes and transcripts, providing accurate quantification even
for overlapping gene families and repetitive elements.

### RSEM Index Preparation (one-time prep, outside the workflow)

```bash
# Prepare RSEM reference (run once)
rsem-prepare-reference --gtf gencode.v38.primary_assembly.annotation.gtf \
  --star GRCh38.primary_assembly.genome.fa rsem_index/GRCh38
```

The last argument is a **prefix**, not a directory, and it is what you pass as
`--rsem_index` (here `rsem_index/GRCh38`). The workflow stages every file whose name
starts with that prefix. The annotation is fixed at this point; there is no `--gtf`
parameter downstream.

### RSEM Quantification

```bash
rsem-calculate-expression \
  --paired-end \
  --bam \
  --no-bam-output \
  --estimate-rspd \
  --strandedness reverse \
  --num-threads 8 \
  Aligned.toTranscriptome.out.bam \
  rsem_index/GRCh38 \
  sample_name
```

### RSEM Output Files

| File | Contents | Key Columns |
|------|----------|-------------|
| `rsem/<sample>.genes.results` | Gene-level quantification | gene_id, transcript_id(s), length, effective_length, expected_count, TPM, FPKM |
| `rsem/<sample>.isoforms.results` | Transcript-level quantification | transcript_id, gene_id, length, effective_length, expected_count, TPM, FPKM, IsoPct |
| `rsem/<sample>.stat/` | Model and read statistics | Directory; parsed by MultiQC |

`--no-bam-output` is passed, so RSEM writes no BAM of its own.

### RSEM Strandedness Flags

The workflow passes its `--strandedness` value straight through to RSEM, so these are the
same three values:

| Library Type | Pipeline and RSEM value | Description |
|-------------|-----------|-------------|
| dUTP / rf-stranded | `--strandedness reverse` | ENCODE standard, the pipeline default |
| fr-stranded | `--strandedness forward` | Directional ligation |
| Unstranded | `--strandedness none` | SMARTer, SMART-Seq2, older protocols |

## Kallisto Quantification (Optional Fast Alternative)

Kallisto uses pseudoalignment (k-mer matching without full alignment) for ultra-fast
transcript quantification. It runs 10-100x faster than STAR+RSEM but does not produce
BAM files or support fusion detection.

In this workflow Kallisto runs on the trimmed FASTQ files, in parallel with (not instead
of) STAR and RSEM. Skip it with `--skip_kallisto`, which also stops `--kallisto_index`
from being read.

### Kallisto Index (one-time prep, outside the workflow)

```bash
# Build Kallisto index from transcriptome FASTA (run once)
kallisto index -i kallisto_index.idx gencode.v38.transcripts.fa
```

Pass the resulting file as `--kallisto_index`. Build it with kallisto 0.50.1, the version
in the image: 0.50.1 writes index version 13, and an index built with kallisto 0.48 or
earlier is rejected at load time. Rebuild rather than reuse an older `.idx`.

### Kallisto Quantification

```bash
kallisto quant \
  -i kallisto_index.idx \
  -o sample/ \
  --rf-stranded \
  -t 8 \
  R1.fq.gz R2.fq.gz
```

The strand flag follows `--strandedness`: `--rf-stranded` for `reverse`, `--fr-stranded`
for `forward`, and no flag for `none`. With `--single_end` the workflow substitutes the
fixed `--single -l 200 -s 20`; if your fragment length distribution differs, rerun
kallisto by hand with the right values.

### Kallisto Output

| File | Contents | Published |
|------|----------|-----------|
| `kallisto/<sample>/abundance.tsv` | transcript_id, length, effective_length, est_counts, tpm | Yes |
| `kallisto/<sample>/run_info.json` | Run metadata and statistics | Yes |
| `kallisto/<sample>/abundance.h5` | Binary HDF5 format (for sleuth) | Only when the kallisto build has HDF5 support; the process declares it as an optional output |

## TPM vs FPKM vs Raw Counts

| Metric | Definition | Cross-Sample Comparable | Use Case |
|--------|-----------|------------------------|----------|
| **Raw counts** | Number of reads/fragments mapped to gene | No | Input for DESeq2/edgeR differential expression |
| **TPM** (Transcripts Per Million) | Counts normalized by gene length then library size | Yes | Cross-sample expression comparison |
| **FPKM** (Fragments Per Kilobase per Million) | Counts normalized by library size then gene length | No | Legacy; avoid for cross-sample comparison |

### When to Use Each

- **Differential expression**: Use raw `expected_count` from RSEM with DESeq2 or edgeR.
  These tools apply their own normalization (median-of-ratios or TMM).
- **Cross-sample comparison**: Use TPM. It sums to 1M per sample, enabling direct comparison.
- **Single-gene reporting**: TPM is appropriate for reporting expression of individual genes.
- **Avoid FPKM**: FPKM does not sum to a constant across samples, making it unreliable
  for cross-sample comparison. TPM is strictly preferred.

## QC Checkpoints

None of these are computed by the workflow; each is a check to run on the published
quantifications (see `05-qc-metrics.md` for the detected-gene command).

| Check | Threshold | Action if Failed |
|-------|-----------|------------------|
| Detected genes (TPM>1) | >12,000 (human) | Check sequencing depth, RNA quality |
| RSEM mapping rate (from `<sample>.stat/`) | >70% of transcriptome BAM reads | Check strandedness setting |
| TPM correlation between replicates | r > 0.95 (Pearson) | Check batch effects, sample swap |
| Gene count distribution | Log-normal shape expected | Skewed distribution suggests degradation |
