# Stage 2: Alignment

## Tools
- **BWA-MEM 0.7.18** (image version): Primary aligner for ChIP-seq (Li & Durbin, 2009)
- **Samtools 1.19** (image version): BAM conversion, sorting, indexing, and statistics

## Reference Genome

The workflow does not download or build an index. It expects a directory (`--bwa_index`,
default `./<genome>_index`) that already contains the FASTA and its BWA index files,
because BWA is invoked as `bwa mem ... <dir>/<genome>.fa`:

```
GRCh38_index/
  GRCh38.fa
  GRCh38.fa.amb
  GRCh38.fa.ann
  GRCh38.fa.bwt
  GRCh38.fa.pac
  GRCh38.fa.sa
```

Build it once from the ENCODE "no alt" analysis-set FASTA:

| Organism | Assembly | FASTA source |
|----------|----------|--------------|
| Human | GRCh38 (hg38) | `https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/` |
| Mouse | mm10 (GRCm38) | `https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/` |

```bash
mkdir -p GRCh38_index
cp GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta GRCh38_index/GRCh38.fa
bwa index GRCh38_index/GRCh38.fa

# Chromosome sizes for --chrom_sizes (Stage 5 needs them)
samtools faidx GRCh38_index/GRCh38.fa
cut -f1,2 GRCh38_index/GRCh38.fa.fai > GRCh38.chrom.sizes
```

**Important**: Use the ENCODE "no alt" analysis set which excludes alternate haplotype
contigs and decoy sequences. This prevents ambiguous multi-mapping to alternate loci.

## Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| BWA algorithm | mem | Recommended for reads >70bp |
| BWA threads (`-t`) | 8 | `BWA_MEM` is configured with 8 CPUs in `nextflow.config` |
| BWA `-M` flag | yes | Mark shorter split hits as secondary (Picard compatible) |
| MAPQ filter | 30 | Remove multi-mappers (MAPQ<30) |
| Sort order | coordinate | Required for Picard and peak calling |
| Sort memory | 2G per thread | `samtools sort -m 2G`; 8 threads x 2G fits the 32 GB process limit |

## Commands

The workflow runs the equivalent of:

```bash
# Paired-end alignment
bwa mem -t 8 -M GRCh38_index/GRCh38.fa trimmed_R1.fq.gz trimmed_R2.fq.gz | \
  samtools view -@ 8 -bS -q 30 - | \
  samtools sort -@ 8 -m 2G -o aligned.bam -
samtools index aligned.bam

# Single-end alignment (--single_end)
bwa mem -t 8 -M GRCh38_index/GRCh38.fa trimmed.fq.gz | \
  samtools view -@ 8 -bS -q 30 - | \
  samtools sort -@ 8 -m 2G -o aligned.bam -
samtools index aligned.bam

# Alignment statistics (the workflow runs flagstat only)
samtools flagstat aligned.bam > flagstat.txt
```

`samtools idxstats` and `samtools stats` are useful manual follow-ups but are not run by
the workflow.

## Expected Output
- `aligned/<sample>.bam` + `.bam.bai` -- coordinate-sorted, indexed BAM
- `aligned/<sample>.flagstat.txt` -- alignment summary (total, mapped, paired, properly paired)

When `--control` is given, the control libraries appear here too, prefixed `CONTROL_`.

## QC Checkpoints

| Check | Threshold | Action if Failed |
|-------|-----------|------------------|
| Mapping rate | >80% | Check genome build match, contamination |
| Mapped reads | >=20M (TF), >=45M (histone) | Sequence more or pool |
| MAPQ>=30 fraction | >70% of mapped reads | Check for repetitive regions enrichment |
| Properly paired (PE) | >90% of mapped | Check library preparation |
| Mitochondrial reads | <5% (manual, via `samtools idxstats`) | Not measured or filtered by this workflow |

## Notes

- The MAPQ 30 filter removes reads mapping to multiple locations. For repetitive element
  analysis, consider relaxing this threshold (it is fixed in `main.nf`).
- BWA-MEM is preferred over BWA-ALN for reads longer than 70bp. For older datasets with
  shorter reads, BWA-ALN may be more appropriate.
- The `-M` flag ensures compatibility with Picard MarkDuplicates in Stage 3.
- This workflow performs no mitochondrial filtering at any stage.
