# Stage 2: Alignment

## Tools
- **Bowtie2 2.5.4** (image version): Primary aligner for ATAC-seq (Langmead & Salzberg 2012, ~30,000 citations)
- **Samtools 1.19** (image version): BAM conversion, sorting, indexing, and statistics

## Why Bowtie2 Instead of BWA
Bowtie2 is preferred for ATAC-seq because:
1. Better handling of short fragments (NFR <150 bp)
2. `--very-sensitive` mode provides optimal alignment for ATAC-seq read characteristics
3. Concordant paired-end alignment with fragment size constraints
4. Standard in the ENCODE ATAC-seq pipeline

## Reference Genome

The workflow does not download or build an index. It expects a directory
(`--bowtie2_index`, default `./<genome>_bowtie2_index`) whose index files are named after
the genome, because Bowtie2 is invoked as `bowtie2 ... -x <dir>/<genome>`:

```
GRCh38_bowtie2_index/
  GRCh38.1.bt2  GRCh38.2.bt2  GRCh38.3.bt2  GRCh38.4.bt2
  GRCh38.rev.1.bt2  GRCh38.rev.2.bt2
```

Build it once from the ENCODE "no alt" analysis-set FASTA:

| Organism | Assembly | FASTA source |
|----------|----------|--------------|
| Human | GRCh38 (hg38) | `https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/` |
| Mouse | mm10 (GRCm38) | `https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/` |

```bash
mkdir -p GRCh38_bowtie2_index
bowtie2-build --threads 8 GRCh38.fa GRCh38_bowtie2_index/GRCh38
```

Use the ENCODE "no alt" analysis set (excludes alternate haplotypes and decoys).

## Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| Bowtie2 mode | --very-sensitive | Optimal sensitivity for ATAC-seq |
| Threads | 8 | `BOWTIE2_ALIGN` is configured with 8 CPUs in `nextflow.config` |
| Max fragment size (-X) | 2000 | Accommodate di/tri-nucleosomal fragments |
| No mixed (--no-mixed) | yes | Require both mates to align |
| No discordant (--no-discordant) | yes | Require concordant alignment |
| MAPQ filter | 30 | Remove multi-mappers |
| Proper pairs only | `-f 2` | Applied in the same samtools step |
| Sort memory | 2G per thread | `samtools sort -m 2G`; 8 threads x 2G fits the 32 GB process limit |

## Commands

The workflow runs the equivalent of:

```bash
# Paired-end alignment with Bowtie2
bowtie2 --very-sensitive -X 2000 --no-mixed --no-discordant \
  --threads 8 -x GRCh38_bowtie2_index/GRCh38 \
  -1 trimmed_R1.fq.gz -2 trimmed_R2.fq.gz 2> sample.bowtie2.log | \
  samtools view -@ 8 -bS -q 30 -f 2 - | \
  samtools sort -@ 8 -m 2G -o aligned.bam -
samtools index aligned.bam

# Alignment statistics (the workflow runs flagstat only)
samtools flagstat aligned.bam > flagstat.txt
```

`samtools idxstats` is run on this BAM in Stage 3 and published as
`qc/<sample>.idxstats.txt`, which is where the mitochondrial fraction comes from.

## Expected Output
- `aligned/<sample>.bam` + `.bam.bai` -- coordinate-sorted, indexed BAM
- `aligned/<sample>.flagstat.txt` -- alignment summary
- `aligned/<sample>.bowtie2.log` -- Bowtie2 alignment rate summary (parsed by MultiQC)

## QC Checkpoints

| Check | Threshold | Action if Failed |
|-------|-----------|------------------|
| Overall alignment rate | >80% | Check genome build, contamination |
| Concordant pair rate | >90% | Check library prep |
| MAPQ>=30 fraction | >60% of mapped | Expected lower than ChIP due to open chromatin |
| Mitochondrial reads | <20% | Measured and removed in Stage 3 |
| Total mapped reads | >=25M after mito removal | May need deeper sequencing |

## Notes
- The `-f 2` flag retains only properly paired reads; combined with `-q 30` this is the
  only BAM flag filtering the workflow performs.
- `-X 2000` allows Bowtie2 to map large fragments from di/tri-nucleosomal DNA.
- Mitochondrial read fraction varies by cell type and preparation; it is measured and
  filtered in Stage 3.
