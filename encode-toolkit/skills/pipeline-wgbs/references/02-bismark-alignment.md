# Bismark Alignment for WGBS

Bisulfite-converted reads require a specialized aligner that accounts for
C-to-T conversion. Bismark handles this by aligning to both C-to-T and G-to-A
converted genomes simultaneously.

## Genome Preparation (one-time prep, outside the workflow)

Build the bisulfite-converted genome index. The directory you prepare is what you pass as
`--genome_dir`, and it must still contain the genome `.fa`, because MethylDackel reads the
FASTA from the same directory later.

```bash
bismark_genome_preparation \
    --bowtie2 \
    --parallel 4 \
    --verbose \
    /ref/genome/
```

This creates two converted genomes:
- C-to-T converted (for original top strand and complementary to bottom)
- G-to-A converted (for original bottom strand and complementary to top)

Requires approximately 12 GB disk space for human genome.

## Bismark Alignment

This is what the workflow runs, in the task directory:

```bash
mkdir -p tmp
bismark \
    --genome /ref/genome/ \
    --bowtie2 \
    --parallel 4 \
    --score_min L,0,-0.2 \
    --no_mixed \
    --no_discordant \
    --maxins 1000 \
    --temp_dir $PWD/tmp \
    -1 sample_R1_val_1.fq.gz \
    -2 sample_R2_val_2.fq.gz
```

Of the two outputs, only the `*_PE_report.txt` is published (to
`bismark/alignments/`). The raw unsorted BAM stays in the work directory; what reaches
`bismark/alignments/<sample>.sorted.bam` is the deduplicated, coordinate-sorted BAM from
the later step.

### Key Parameters

| Parameter | Value | Reason |
|-----------|-------|--------|
| `--parallel 4` | 4 instances | Each uses 2 Bowtie2 threads = 8 total threads |
| `--score_min L,0,-0.2` | Linear penalty | ENCODE default; tolerates bisulfite mismatches |
| `--no_mixed` | Discard | Both mates must align |
| `--no_discordant` | Discard | Mates must be properly paired |
| `--maxins 1000` | 1000 bp | Maximum insert size for paired reads |

### Memory Requirements

Bismark parallel mode runs several Bowtie2 instances, each holding its own copy of a
converted genome index. `nextflow.config` allocates 48 GB to `BISMARK_ALIGN` with
`--parallel 4` on a human genome, and multiplies it by the attempt number on a retry,
capped by `--max_memory` (64 GB by default, so a retry gets 64 GB). Budget for the 48 GB
figure rather than a per-instance estimate.

## Alternative: bwa-meth

bwa-meth is not part of the bundled workflow or container image. The commands below are for
running it by hand when Bismark is too slow or memory-hungry for your genome.

bwa-meth is faster than Bismark for large genomes and uses less memory:

```bash
# Index genome (one-time)
bwameth.py index /ref/genome/genome.fa

# Align
bwameth.py \
    --threads 8 \
    --reference /ref/genome/genome.fa \
    sample_R1_trimmed.fq.gz \
    sample_R2_trimmed.fq.gz \
    | samtools sort -@ 4 -o sample_bwameth.bam

samtools index sample_bwameth.bam
```

### Bismark vs bwa-meth Comparison

| Feature | Bismark | bwa-meth |
|---------|---------|----------|
| Speed | Slower (~2x) | Faster |
| RAM | 48 GB allocated (parallel) | ~16 GB |
| Accuracy | Gold standard | Comparable |
| This workflow | The only aligner it runs | Not available; manual only |
| Methylation calling | Built-in | Requires MethylDackel |

Recommendation: Use Bismark for ENCODE compatibility. Use bwa-meth outside this workflow
when processing many samples and speed is critical.

## Lambda/pUC19 Spike-in Alignment (not run by this workflow)

If spike-in DNA was used, align to the spike-in genome by hand to measure the conversion
rate. The workflow does not do this and reports no conversion rate:

```bash
bismark \
    --genome /ref/lambda/ \
    --bowtie2 \
    -1 sample_R1_trimmed.fq.gz \
    -2 sample_R2_trimmed.fq.gz \
    --output_dir lambda_out/ \
    --unmapped
```

The lambda genome is fully unmethylated, so any detected methylation represents
incomplete bisulfite conversion. Expect conversion rate ≥98%.

## Alignment QC Checks

After alignment, verify in the Bismark report (`bismark/alignments/*_PE_report.txt`,
also parsed into the MultiQC report):
- **Mapping efficiency**: >70% for WGBS (lower than standard WGS due to conversion)
- **Unique alignments**: Should dominate over multimappers
- **C methylated in CpG context**: Typically 70-85% for mammalian somatic tissue
- **C methylated in CHG context**: Should be <1% in somatic tissue (>5% in plants/ESCs)
- **C methylated in CHH context**: Should be <1% in somatic tissue
