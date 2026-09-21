# QC and Trimming for Hi-C Data

Hi-C reads contain chimeric sequences from ligation junctions. Ligation
junctions are handled downstream by pairtools.

**This workflow does not trim.** `main.nf` runs FastQC on the raw reads and
passes the same raw reads to `bwa mem -SP5M`, which soft-clips adapter
sequence. Trim Galore is present in the container image but no process calls
it. The trimming commands below are a manual, optional pre-processing step.

## Pre-Alignment QC with FastQC

```bash
fastqc --threads 4 --outdir fastqc_raw/ sample_R1.fastq.gz sample_R2.fastq.gz
```

Key checks:
- Per-base quality (expect Phred >28 across most positions)
- Adapter content (Illumina adapters)
- Sequence length distribution (Hi-C reads are typically 50-150 bp)
- Duplication level (Hi-C libraries often show high duplication at low depths)

**Note**: Hi-C data will show unusual insert size distributions because reads
originate from ligation junctions, not contiguous fragments. This is expected.

## Optional: Adapter Trimming with Trim Galore (not run by this workflow)

If you choose to trim before running the workflow, keep it light: do NOT
aggressively trim quality, since chimeric reads may have lower quality at the
junction. Feed the resulting `*_val_{1,2}.fq.gz` files to `--reads`.

```bash
trim_galore \
    --paired \
    --quality 15 \
    --phred33 \
    --length 30 \
    --cores 4 \
    --fastqc \
    sample_R1.fastq.gz \
    sample_R2.fastq.gz
```

### Parameter Rationale

| Parameter | Value | Reason |
|-----------|-------|--------|
| `--quality 15` | Phred 15 | Lenient -- chimeric reads have junction artifacts |
| `--length 30` | 30 bp | Short reads still carry valid contact information |
| `--cores 4` | 4 | Trim Galore uses ~3x threads internally |

### Default here: No Trimming

Like Juicer, this workflow skips trimming entirely and relies on the aligner
to handle adapter contamination via soft-clipping, which BWA-MEM does:

```bash
# No trim_galore step; raw FASTQs go straight to bwa mem -SP5M,
# which soft-clips adapter sequence.
```

## Restriction Enzyme Verification

Before processing, verify which restriction enzyme was used:

| Enzyme | Recognition Site | Ligation Junction | Average Fragment |
|--------|-----------------|-------------------|------------------|
| MboI/DpnII | GATC | GATCGATC | ~256 bp |
| HindIII | AAGCTT | AAGCTAGCTT | ~4 kb |
| Arima | Two sites | Multiple | ~160 bp |
| NcoI | CCATGG | CCATGCATGG | ~2 kb |

```bash
# Check for ligation junction sequence in reads (MboI example)
zcat sample_R1.fastq.gz | head -10000 | grep -c 'GATCGATC'
```

If the junction sequence appears frequently (>1% of reads), the enzyme
assignment is confirmed.

## Post-Trimming Summary (only if you trimmed manually)

After trimming, verify:
- >95% reads pass filters
- Median read length >40 bp
- Adapter detection rate is reasonable (typically 5-20% for Hi-C)
