#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ENCODE RNA-seq Pipeline — Nextflow DSL2
// FASTQ -> QC -> STAR 2-pass -> RSEM -> Kallisto (optional) -> Signal -> RSeQC

params.reads          = null
params.genome         = 'GRCh38'
params.outdir         = 'results'
params.single_end     = false
params.strandedness   = 'reverse'   // reverse (dUTP), forward, or none
params.skip_kallisto  = false
params.star_index     = null
params.rsem_index     = null   // RSEM reference prefix (rsem-prepare-reference output), not a directory
params.kallisto_index = null   // kallisto index file; not read with --skip_kallisto
params.rseqc_bed      = null   // BED12 gene model for RSeQC
params.chrom_sizes    = null

// Genome-specific defaults. Kept in a function because scripts that declare processes
// cannot also hold top-level variables. The annotation is part of the STAR and RSEM
// references, so it is chosen when those are built rather than passed to this workflow.
def genomeDefaults() {
    return [
    'GRCh38': [
        star_index:  'GRCh38_star_index',
        rsem_index:  'GRCh38_rsem_index/GRCh38',
        kallisto_idx: 'gencode.v38.kallisto.idx',
        rseqc_bed:   'hg38_RefSeq.bed'
    ],
    'mm10': [
        star_index:  'mm10_star_index',
        rsem_index:  'mm10_rsem_index/mm10',
        kallisto_idx: 'gencode.vM27.kallisto.idx',
        rseqc_bed:   'mm10_RefSeq.bed'
    ]
    ]
}

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*.{html,zip}"), emit: reports

    script:
    """
    fastqc -t ${task.cpus} --outdir . ${reads}
    """
}

process TRIM_GALORE {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: '*{.fq.gz,trimming_report.txt}'
    publishDir "${params.outdir}/fastqc",  mode: 'copy', pattern: '*_fastqc.{html,zip}'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*{val_1.fq.gz,val_2.fq.gz,trimmed.fq.gz}"), emit: trimmed
    path("*trimming_report.txt"),                                             emit: log
    path("*_fastqc.{html,zip}"),                                              emit: fastqc

    script:
    if (params.single_end)
        """
        trim_galore --quality 20 --length 36 --fastqc --cores ${task.cpus} ${reads}
        """
    else
        """
        trim_galore --paired --quality 20 --length 36 --fastqc --cores ${task.cpus} ${reads[0]} ${reads[1]}
        """
}

process STAR_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/star", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path(star_index)

    output:
    tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam"),
                          path("${sample_id}.Aligned.sortedByCoord.out.bam.bai"),    emit: genome_bam
    tuple val(sample_id), path("${sample_id}.Aligned.toTranscriptome.out.bam"),      emit: transcriptome_bam
    path("${sample_id}.Log.final.out"),                                               emit: log
    path("${sample_id}.SJ.out.tab"),                                                  emit: junctions
    path("${sample_id}.ReadsPerGene.out.tab"),                                        emit: gene_counts
    // STAR always names the first pair of signal files ".str1", also with --outWigStrand
    // Unstranded, which writes str1 only (STAR 2.7.11b, source/signalFromBAM.cpp lines 45-49).
    tuple val(sample_id), path("${sample_id}.Signal.UniqueMultiple.str*.out.bg"),    emit: bedgraph
    path("${sample_id}.Signal.Unique.str*.out.bg"),                                   emit: bedgraph_unique

    script:
    def input_reads = params.single_end ? "${reads}" : "${reads[0]} ${reads[1]}"
    def wig_strand  = params.strandedness == 'none' ? 'Unstranded' : 'Stranded'
    """
    STAR --genomeDir ${star_index} \\
      --readFilesIn ${input_reads} \\
      --readFilesCommand zcat \\
      --runThreadN ${task.cpus} \\
      --outSAMtype BAM SortedByCoordinate \\
      --outSAMunmapped Within \\
      --outFilterMultimapNmax 20 \\
      --alignSJoverhangMin 8 \\
      --alignSJDBoverhangMin 1 \\
      --outFilterMismatchNmax 999 \\
      --outFilterMismatchNoverReadLmax 0.04 \\
      --alignIntronMin 20 \\
      --alignIntronMax 1000000 \\
      --alignMatesGapMax 1000000 \\
      --quantMode TranscriptomeSAM GeneCounts \\
      --twopassMode Basic \\
      --outWigType bedGraph \\
      --outWigStrand ${wig_strand} \\
      --outFileNamePrefix ${sample_id}.

    samtools index ${sample_id}.Aligned.sortedByCoord.out.bam
    """
}

process RSEM_QUANT {
    tag "$sample_id"
    publishDir "${params.outdir}/rsem", mode: 'copy'

    input:
    tuple val(sample_id), path(transcriptome_bam)
    path(rsem_index)

    output:
    path("${sample_id}.genes.results"),    emit: genes
    path("${sample_id}.isoforms.results"), emit: isoforms
    path("${sample_id}.stat"),             emit: stats

    script:
    def pe_flag = params.single_end ? "" : "--paired-end"
    // The reference files are staged into the task directory, so use the prefix basename.
    def rsem_prefix = file(params.rsem_index ?: genomeDefaults()[params.genome].rsem_index).name
    """
    rsem-calculate-expression \\
      ${pe_flag} \\
      --bam \\
      --no-bam-output \\
      --estimate-rspd \\
      --strandedness ${params.strandedness} \\
      --num-threads ${task.cpus} \\
      ${transcriptome_bam} \\
      ${rsem_prefix} \\
      ${sample_id}
    """
}

process KALLISTO_QUANT {
    tag "$sample_id"
    publishDir "${params.outdir}/kallisto", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path(kallisto_index)

    output:
    path("${sample_id}/abundance.tsv"), emit: abundance
    path("${sample_id}/abundance.h5"),  emit: abundance_h5, optional: true   // only from kallisto builds with HDF5
    path("${sample_id}/run_info.json"), emit: run_info

    when:
    !params.skip_kallisto

    script:
    def strand_flag = params.strandedness == 'reverse' ? '--rf-stranded' :
                      params.strandedness == 'forward' ? '--fr-stranded' : ''
    def input_reads = params.single_end ?
        "--single -l 200 -s 20 ${reads}" : "${reads[0]} ${reads[1]}"
    """
    kallisto quant \\
      -i ${kallisto_index} \\
      -o ${sample_id}/ \\
      ${strand_flag} \\
      -t ${task.cpus} \\
      ${input_reads}
    """
}

process SIGNAL_TRACKS {
    tag "$sample_id"
    publishDir "${params.outdir}/signal", mode: 'copy'

    input:
    tuple val(sample_id), path(bedgraphs)
    path(chrom_sizes)

    output:
    path("${sample_id}_*.bw"), emit: bigwig

    script:
    // STAR writes str1 for fragments whose read 1 maps to the + strand. Read 1 of a
    // reverse-stranded (dUTP) library is antisense, so str1 is minus-strand signal and str2 is
    // plus-strand signal, as in ENCODE's STAR_RSEM.sh (str[1]=-; str[2]=+). Forward-stranded
    // libraries are the opposite, and unstranded libraries give a single str1 track.
    def tracks = params.strandedness == 'reverse' ? [str1: 'minus', str2: 'plus'] :
                 params.strandedness == 'forward' ? [str1: 'plus', str2: 'minus'] :
                                                    [str1: 'unstranded']
    def commands = tracks.collect { str, strand ->
        "sort -k1,1 -k2,2n ${sample_id}.Signal.UniqueMultiple.${str}.out.bg > ${strand}_sorted.bg\n" +
        "    bedGraphToBigWig ${strand}_sorted.bg ${chrom_sizes} ${sample_id}_${strand}.bw"
    }.join('\n    ')
    """
    ${commands}
    """
}

process RSEQC {
    tag "$sample_id"
    publishDir "${params.outdir}/qc/rseqc", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path(rseqc_bed)

    output:
    path("${sample_id}.infer_experiment.txt"),     emit: strandedness
    path("${sample_id}.read_distribution.txt"),    emit: read_dist
    path("${sample_id}.geneBody_coverage*"),        emit: gene_body
    path("${sample_id}.inner_distance*"),           emit: inner_dist, optional: true

    script:
    def inner_dist_cmd = params.single_end ? "" :
        "inner_distance.py -r ${rseqc_bed} -i ${bam} -o ${sample_id}.inner_distance"
    """
    infer_experiment.py -r ${rseqc_bed} -i ${bam} > ${sample_id}.infer_experiment.txt
    read_distribution.py -r ${rseqc_bed} -i ${bam} > ${sample_id}.read_distribution.txt
    geneBody_coverage.py -r ${rseqc_bed} -i ${bam} -o ${sample_id}.geneBody_coverage
    ${inner_dist_cmd}
    """
}

process MULTIQC {
    publishDir "${params.outdir}/qc/multiqc", mode: 'copy'

    input:
    path('*')

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data"),        emit: data

    script:
    """
    multiqc . -o . -f
    """
}

workflow {
    // ---- Parameter validation ----
    if (!params.reads) { error "Missing required parameter: --reads" }
    if (!genomeDefaults().containsKey(params.genome)) {
        error "Unsupported --genome '${params.genome}': expected one of ${genomeDefaults().keySet().join(', ')}"
    }
    if (!(params.strandedness in ['reverse', 'forward', 'none'])) {
        error "Invalid --strandedness '${params.strandedness}': expected 'reverse', 'forward', or 'none'"
    }

    // Google Batch and AWS Batch stage every task through object storage, so the matching
    // profile cannot run without a bucket work directory and a project or job queue.
    def active_profiles = workflow.profile.tokenize(',')
    def work_uri        = workflow.workDir.toUriString()
    if (active_profiles.contains('gcp') && !(params.gcp_project && work_uri.startsWith('gs://'))) {
        error "-profile gcp requires --gcp_project <project-id> and --gcp_workdir gs://<bucket>/work"
    }
    if (active_profiles.contains('aws') && !(params.aws_queue && work_uri.startsWith('s3://'))) {
        error "-profile aws requires --aws_queue <job-queue> and --aws_workdir s3://<bucket>/work"
    }

    // ---- Input channels ----
    def defaults      = genomeDefaults()[params.genome]
    def star_idx_path = params.star_index ?: defaults.star_index
    def rsem_idx_path = params.rsem_index ?: defaults.rsem_index

    ch_reads      = channel.fromFilePairs(params.reads, size: params.single_end ? 1 : 2, checkIfExists: true)
    ch_star_idx   = channel.fromPath(star_idx_path, type: 'dir', checkIfExists: true)
    // rsem_idx_path is a prefix such as .../GRCh38, so stage every file that starts with it
    ch_rsem_idx   = channel.fromPath("${rsem_idx_path}*", checkIfExists: true)
    ch_rseqc_bed  = channel.fromPath(params.rseqc_bed ?: defaults.rseqc_bed, checkIfExists: true)

    // Stage 1: QC and Trimming
    FASTQC(ch_reads)
    TRIM_GALORE(ch_reads)

    // Stage 2: STAR Alignment (2-pass mode)
    STAR_ALIGN(TRIM_GALORE.out.trimmed, ch_star_idx.collect())

    // Stage 3: Quantification
    RSEM_QUANT(STAR_ALIGN.out.transcriptome_bam, ch_rsem_idx.collect())

    // Kallisto (optional)
    if (!params.skip_kallisto) {
        ch_kallisto_idx = channel.fromPath(params.kallisto_index ?: defaults.kallisto_idx, checkIfExists: true)
        KALLISTO_QUANT(TRIM_GALORE.out.trimmed, ch_kallisto_idx.collect())
    }

    // Stage 4: Signal Tracks
    ch_chrom_sizes = params.chrom_sizes ?
        channel.fromPath(params.chrom_sizes, checkIfExists: true) :
        channel.fromPath("${star_idx_path}/chrNameLength.txt", checkIfExists: true)
    SIGNAL_TRACKS(STAR_ALIGN.out.bedgraph, ch_chrom_sizes.collect())

    // Stage 5: QC Metrics
    RSEQC(STAR_ALIGN.out.genome_bam, ch_rseqc_bed.collect())

    // MultiQC aggregation
    ch_multiqc = FASTQC.out.reports
        .mix(TRIM_GALORE.out.log)
        .mix(TRIM_GALORE.out.fastqc)
        .mix(STAR_ALIGN.out.log)
        .mix(RSEM_QUANT.out.stats)
        .mix(RSEQC.out.strandedness)
        .mix(RSEQC.out.read_dist)
        .collect()
    MULTIQC(ch_multiqc)
}
