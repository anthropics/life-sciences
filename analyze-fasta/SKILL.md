---
name: analyze-fasta
description: Analyzes FASTA files (nucleotide or protein) with automatic sequence type detection, comprehensive bioinformatics metrics, and biological interpretation. Use when users want to analyze DNA/RNA sequences, protein sequences, assess assembly quality, find ORFs, or compute physicochemical properties.
---

# FASTA Sequence Analyzer

Automated analysis workflow for FASTA files supporting both nucleotide and protein sequences.

## When to Use This Skill

Use when users:
- Want to analyze a FASTA file (nucleotide or protein)
- Need sequence composition metrics (GC%, base composition, amino acid composition)
- Want to find open reading frames (ORFs) in a nucleotide sequence
- Need physicochemical properties of a protein (molecular weight, pI, GRAVY, stability)
- Want an HTML report with visualizations of sequence properties
- Need assembly quality metrics (N50, length distribution)
- Ask for biological interpretation of sequence data

**Supported input formats:**
- `.fasta`, `.fa`, `.fna`, `.faa` files
- Single or multi-sequence FASTA files
- Nucleotide (DNA/RNA) and protein sequences (auto-detected)

## Quick Start

```bash
# Text report to terminal
python3 scripts/analyze_fasta.py input.fasta

# JSON output (for pipelines or programmatic use)
python3 scripts/analyze_fasta.py input.fasta --json

# HTML report (auto-generates in output/ directory)
python3 scripts/analyze_fasta.py input.fasta --html
```

**Requirements:** biopython, numpy

Install with: `pip install biopython numpy`

## What It Analyzes

### For Nucleotide Sequences

| Metric | Description |
|--------|-------------|
| Length | Sequence length in base pairs |
| GC Content | Percentage of G+C bases |
| AT Content | Percentage of A+T bases |
| Base Composition | Count and percentage of each nucleotide (A, T, G, C, N) |
| Top Dinucleotides | 5 most frequent two-base combinations |
| ORFs | Open reading frames in 3 forward frames (minimum 100 aa / 300 bp) |
| Molecular Weight | Estimated molecular weight of the DNA strand |
| N50 | Assembly quality metric (for multi-sequence files) |

### For Protein Sequences

| Metric | Description |
|--------|-------------|
| Length | Sequence length in amino acids |
| Molecular Weight | In Daltons, via Biopython ProteinAnalysis |
| Isoelectric Point (pI) | pH at which the protein has no net charge |
| Instability Index | Classifies as stable (<40) or unstable (>=40) |
| GRAVY | Grand Average of Hydropathy. Negative = hydrophilic/soluble, positive = hydrophobic/membrane |
| Aromaticity | Fraction of aromatic amino acids (Phe, Trp, Tyr) |
| Secondary Structure | Predicted helix, sheet, and turn fractions |
| Charged Residues | Percentage of Arg, Lys, Asp, Glu |
| Aromatic Residues | Percentage of Phe, Trp, Tyr |
| Amino Acid Composition | Count and percentage of each amino acid |

### Summary Statistics (Multi-Sequence Files)

| Metric | Description |
|--------|-------------|
| Total Sequences | Number of sequences in the file |
| Total Residues | Sum of all sequence lengths |
| Min / Max / Avg Length | Length distribution |
| N50 | Minimum length such that 50% of total residues are in sequences of that length or longer |
| Average GC | Mean GC content across all sequences (nucleotide only) |
| Total ORFs | Sum of all ORFs found (nucleotide only) |

## Output Modes

### Text Report (default)
Human-readable summary printed to stdout. Shows summary statistics followed by per-sequence details (up to 20 sequences).

### JSON (`--json`)
Structured JSON output suitable for programmatic consumption or piping to other tools. Contains the complete analysis including all metrics, compositions, and ORF coordinates.

### HTML Report (`--html`)
Interactive HTML report with:
- Summary statistics table
- Per-sequence detail cards with visual elements:
  - **Nucleotides**: colored base composition bar (A/T/G/C), dinucleotide frequency charts, ORF table with frame and coordinates
  - **Proteins**: metric cards (MW, pI, GRAVY, stability), secondary structure prediction bars
- Collapsible JSON data viewer
- "How it works" section explaining the analysis pipeline and libraries used
- Responsive design, works on desktop and mobile

When `--html` is used without a filename, the report is auto-generated at `output/output_<fasta_name>.html` relative to the input file.

## Biological Interpretation

After running the script, Claude should interpret the results in biological context:

### For Nucleotide Sequences
- **Sequence type inference**: Based on length and GC content, suggest whether it is genomic DNA, a plasmid, rRNA gene, coding region, etc.
- **GC content context**: ~41% is typical for human, ~50% for many bacteria, >60% suggests thermophiles, <35% suggests endosymbionts or AT-rich organisms.
- **ORF analysis**: Evaluate gene density. Long ORFs in a single frame suggest coding regions. Multiple short ORFs across frames suggest non-coding or intergenic regions.
- **Organism inference**: If the header contains organism info, provide biological context. If not, use GC% and composition to suggest possible taxonomic groups.

### For Protein Sequences
- **Protein family inference**: Based on MW, pI, amino acid composition, and secondary structure predictions, suggest likely protein families or functions.
- **Localization**: pI and GRAVY together suggest localization (cytoplasmic, membrane, secreted).
- **Stability assessment**: Instability index in biological context (in vivo half-life implications).
- **Functional clues**: High charged residue content suggests DNA-binding; high hydrophobicity suggests membrane association; specific amino acid biases suggest metabolic roles.

### Always
- Be quantitative: compare values against known references.
- If the FASTA header contains gene/organism identifiers, use them for richer context.
- Suggest follow-up analyses: BLAST for homology search, multiple sequence alignment, phylogenetic inference, structure prediction with AlphaFold, domain search with InterPro.

## ORF Detection Details

The script searches for open reading frames in the 3 forward reading frames (+1, +2, +3):
- Start codon: `ATG`
- Stop codons: `TAA`, `TAG`, `TGA`
- Minimum length: 300 bp (100 amino acids)

**Limitations**: Does not search the reverse complement strand. Does not model introns or splicing. For complete gene prediction, use dedicated tools like Augustus, GenScan, or Prodigal.

## How It Works

```
FASTA File
    |
    v
detect_sequence_type()     Classifies as nucleotide or protein
    |                      (threshold: >85% ATCGUN chars = nucleotide)
    v
analyze_nucleotide()  OR  analyze_protein()
    |                          |
    |  GC%, composition,       |  MW, pI, GRAVY,
    |  dinucleotides, ORFs,    |  stability, secondary
    |  molecular weight        |  structure, composition
    v                          v
analyze_fasta()            Orchestrates analysis,
    |                      computes summary stats (N50, etc.)
    v
format_text_report()  OR  json.dumps()  OR  generate_html_report()
```

**Libraries used:**
- `Bio.SeqIO` (Biopython) for FASTA parsing
- `Bio.SeqUtils.gc_fraction` for GC content
- `Bio.SeqUtils.molecular_weight` for DNA molecular weight
- `Bio.SeqUtils.ProtParam.ProteinAnalysis` for all protein physicochemical properties
- `numpy` (Biopython dependency)
- Python standard library: `collections.Counter`, `json`, `argparse`, `pathlib`

## Author

Created by Santiago Rodriguez ([@santiago-rodriguezs](https://github.com/santiago-rodriguezs)) as part of a bioinformatics seminar project exploring Claude's capabilities for omics analysis.
