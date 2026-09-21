# BWA `-SP5M` Alignment for Hi-C

Hi-C reads are chimeric: each mate can originate from a different genomic
location due to proximity ligation. Standard paired-end alignment fails
because mates are not from a contiguous fragment. The solution is `bwa mem
-SP5M`, which skips pairing and mate rescue while still writing both mates to
one BAM.

## Genome Index Preparation

```bash
# Build BWA index (one-time step, ~1 hour for human genome)
bwa index -a bwtsw genome.fa
```

Requires ~8 GB disk and ~8 GB RAM for human genome.

## Alignment Strategy Used by This Workflow

Align both mates in a single `bwa mem -SP5M` call and write one BAM:

```bash
bwa mem -t 8 -SP5M genome.fa \
    sample_R1.fastq.gz sample_R2.fastq.gz \
    | samtools view -@ 4 -bhS - \
    > sample.paired.bam
```

`main.nf` publishes this file as `alignment/{sample}.paired.bam`. It keeps the
pair information in one file, which pairtools parses directly, and it is the
ENCODE-recommended approach.

### BWA-MEM Flags for Hi-C

| Flag | Meaning | Reason |
|------|---------|--------|
| `-S` | Skip mate rescue | Mates are on different chromosomes |
| `-P` | Skip pairing | Do not try to pair mates |
| `-5` | Split alignment: primary = 5' end | Ensures primary alignment is the 5' portion |
| `-M` | Mark shorter split as secondary | Compatible with downstream tools |

These flags are CRITICAL for Hi-C. Without `-SP`, BWA will try to pair
mates expecting a standard insert size, which fails for Hi-C contacts.

## Per-Mate Alignment (Alternative, not used here)

Some pipelines align R1 and R2 independently as single-end reads and merge the
two BAMs before parsing:

```bash
bwa mem -t 8 -SP5M genome.fa sample_R1.fastq.gz \
    | samtools view -@ 4 -bS - > sample_R1.bam

bwa mem -t 8 -SP5M genome.fa sample_R2.fastq.gz \
    | samtools view -@ 4 -bS - > sample_R2.bam
```

This workflow does not do this and never writes per-mate BAMs.

## Alignment QC (manual)

`main.nf` runs no flagstat on the Hi-C BAM; run it yourself on the published
file if you need it:

```bash
samtools flagstat alignment/sample.paired.bam
```

Expected metrics:
- **Mapped**: >80% (both mates combined)
- **Properly paired**: Low percentage is NORMAL for Hi-C (mates are on different chromosomes)
- **Supplementary alignments**: 5-20% (chimeric reads split across junction)

### Chimeric Read Handling

BWA-MEM produces supplementary alignments for reads spanning ligation junctions.
The `-5` flag ensures the 5' portion is reported as primary:

```
Read spans junction:
  [ChromA:100-200]---GATCGATC---[ChromB:500-600]

BWA reports:
  Primary:      ChromA:100-200  (5' end)
  Supplementary: ChromB:500-600  (3' end)
```

pairtools uses the primary alignment for contact calling.

## MAPQ Filtering

Filter low-quality alignments after pairing (in pairtools step), not here.
Keeping all alignments allows pairtools to classify pair types correctly.

Typical MAPQ threshold: 30 (applied during pairtools parse).

## Memory and Time Estimates

| Genome | Threads | RAM | Time (500M reads) |
|--------|---------|-----|---------------------|
| Human (hg38) | 8 | 12 GB | 3-5 hours |
| Mouse (mm10) | 8 | 10 GB | 2-4 hours |

BWA-MEM memory usage scales with index size, not read count. 8 GB is
sufficient for human genome alignment.
