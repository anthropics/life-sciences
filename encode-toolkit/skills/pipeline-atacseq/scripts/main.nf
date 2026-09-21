#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ENCODE ATAC-seq Pipeline — Nextflow DSL2
// FASTQ -> QC -> Bowtie2 -> Tn5 Shift -> Filter -> Peaks -> IDR -> Signal

params.reads         = null
params.genome        = 'GRCh38'
params.outdir        = 'results'
params.single_end    = false
params.skip_idr      = false
params.blacklist     = null   // defaults to the ENCODE blacklist v2 for --genome
params.bowtie2_index = null   // directory holding the Bowtie2 index <genome>.*.bt2; default ./<genome>_bowtie2_index
params.mito_name     = 'chrM'
params.nfr_max       = 150

// Genome-specific defaults. Kept in a function because scripts that declare processes
// cannot also hold top-level variables.
def genomeDefaults() {
    return [
    'GRCh38': [
        index:     'GRCh38_bowtie2_index',
        blacklist: 'https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz',
        gsize:     'hs'
    ],
    'mm10': [
        index:     'mm10_bowtie2_index',
        blacklist: 'https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz',
        gsize:     'mm'
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
        trim_galore --nextera --quality 20 --length 20 --fastqc --cores ${task.cpus} ${reads}
        """
    else
        """
        trim_galore --paired --nextera --quality 20 --length 20 --fastqc --cores ${task.cpus} ${reads[0]} ${reads[1]}
        """
}

process BOWTIE2_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path(genome_index)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam
    path("${sample_id}.flagstat.txt"),                                              emit: flagstat
    path("${sample_id}.bowtie2.log"),                                               emit: log

    script:
    def input_reads = params.single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    bowtie2 --very-sensitive -X 2000 --no-mixed --no-discordant \\
      --threads ${task.cpus} -x ${genome_index}/${params.genome} \\
      ${input_reads} 2> ${sample_id}.bowtie2.log | \\
      samtools view -@ ${task.cpus} -bS -q 30 -f 2 - | \\
      samtools sort -@ ${task.cpus} -m 2G -o ${sample_id}.bam -
    samtools index ${sample_id}.bam
    samtools flagstat ${sample_id}.bam > ${sample_id}.flagstat.txt
    """
}

process MITO_FILTER {
    tag "$sample_id"
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.idxstats.txt'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.no_mito.bam"), emit: bam
    path("${sample_id}.idxstats.txt"),                       emit: stats

    script:
    """
    # Reads per chromosome before filtering: MultiQC's samtools module turns this into the
    # mitochondrial fraction. It also works when the genome has no ${params.mito_name} contig.
    samtools idxstats ${bam} > ${sample_id}.idxstats.txt

    # Remove mitochondrial reads
    samtools view -@ ${task.cpus} -b ${bam} \$(samtools idxstats ${bam} | \\
      awk '\$1 != "${params.mito_name}" && \$1 != "*" {print \$1}' | tr '\\n' ' ') > ${sample_id}.no_mito.bam
    """
}

process MARK_DUPLICATES {
    tag "$sample_id"
    publishDir "${params.outdir}/filtered", mode: 'copy', pattern: '*_metrics.txt'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bam.bai"), emit: bam
    path("${sample_id}.dup_metrics.txt"),                                                       emit: metrics

    script:
    """
    samtools sort -@ ${task.cpus} -o sorted.bam ${bam}
    picard MarkDuplicates \\
      INPUT=sorted.bam OUTPUT=${sample_id}.dedup.bam \\
      METRICS_FILE=${sample_id}.dup_metrics.txt \\
      REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
    samtools index ${sample_id}.dedup.bam
    """
}

process TN5_SHIFT {
    tag "$sample_id"
    publishDir "${params.outdir}/filtered/shifted", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)   // alignmentSieve needs the index

    output:
    tuple val(sample_id), path("${sample_id}.shifted.bam"), path("${sample_id}.shifted.bam.bai"), emit: bam

    script:
    """
    alignmentSieve --bam ${bam} --outFile shifted_unsorted.bam \\
      --ATACshift --numberOfProcessors ${task.cpus}
    samtools sort -@ ${task.cpus} -o ${sample_id}.shifted.bam shifted_unsorted.bam
    samtools index ${sample_id}.shifted.bam
    """
}

process BLACKLIST_FILTER {
    tag "$sample_id"
    publishDir "${params.outdir}/filtered", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path(blacklist_bed)

    output:
    tuple val(sample_id), path("${sample_id}.final.bam"), path("${sample_id}.final.bam.bai"), emit: bam
    path("${sample_id}.final.flagstat.txt"),                                                    emit: flagstat

    script:
    """
    bedtools intersect -v -abam ${bam} -b ${blacklist_bed} > ${sample_id}.final.bam
    samtools index ${sample_id}.final.bam
    samtools flagstat ${sample_id}.final.bam > ${sample_id}.final.flagstat.txt
    """
}

process NFR_SELECTION {
    tag "$sample_id"
    publishDir "${params.outdir}/filtered/nfr", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.nfr.bam"), path("${sample_id}.nfr.bam.bai"), emit: nfr
    tuple val(sample_id), path("${sample_id}.mononuc.bam"),                                emit: mononuc

    script:
    """
    alignmentSieve --bam ${bam} --outFile nfr_unsorted.bam \\
      --maxFragmentLength ${params.nfr_max} --numberOfProcessors ${task.cpus}
    samtools sort -@ ${task.cpus} -o ${sample_id}.nfr.bam nfr_unsorted.bam
    samtools index ${sample_id}.nfr.bam

    alignmentSieve --bam ${bam} --outFile mononuc_unsorted.bam \\
      --minFragmentLength ${params.nfr_max} --maxFragmentLength 300 --numberOfProcessors ${task.cpus}
    samtools sort -@ ${task.cpus} -o ${sample_id}.mononuc.bam mononuc_unsorted.bam
    """
}

process MACS2_CALLPEAK {
    tag "$sample_id"
    publishDir "${params.outdir}/peaks/narrow", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}*narrowPeak"), emit: peaks
    path("${sample_id}*.bdg"),                              emit: bdg
    path("${sample_id}*.xls"),                              emit: xls
    path("${sample_id}*_summits.bed"),                      emit: summits

    script:
    def gsize = genomeDefaults()[params.genome].gsize
    """
    macs2 callpeak -t ${bam} \\
      -f BAMPE -g ${gsize} -n ${sample_id} \\
      --nomodel --keep-dup all --call-summits \\
      --qvalue 0.05 -B
    """
}

process IDR_ANALYSIS {
    tag "${rep1_id}_vs_${rep2_id}"
    publishDir "${params.outdir}/peaks/idr", mode: 'copy'

    input:
    tuple val(rep1_id), path(rep1_peaks), val(rep2_id), path(rep2_peaks)

    output:
    path("${rep1_id}_vs_${rep2_id}.idr_peaks.txt"),     emit: peaks
    path("${rep1_id}_vs_${rep2_id}.idr_peaks.txt.png"), emit: plot, optional: true

    script:
    """
    idr --samples ${rep1_peaks} ${rep2_peaks} \\
      --input-file-type narrowPeak --rank p.value \\
      --output-file ${rep1_id}_vs_${rep2_id}.idr_peaks.txt --plot --idr-threshold 0.05
    """
}

process SIGNAL_TRACKS {
    tag "$sample_id"
    publishDir "${params.outdir}/signal", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path("${sample_id}.signal.bw"), emit: bw

    script:
    """
    bamCoverage -b ${bam} -o ${sample_id}.signal.bw \\
      --normalizeUsing RPKM --binSize 10 \\
      --numberOfProcessors ${task.cpus} --extendReads
    """
}

process FRIP {
    tag "$sample_id"
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(peaks)

    output:
    path("${sample_id}.frip_mqc.tsv"), emit: frip

    script:
    // Fraction of reads in peaks: alignments of the final BAM that overlap a called peak,
    // over all alignments of that BAM. One row per peak file, in a table MultiQC picks up.
    """
    total=\$(samtools view -c ${bam})
    {
        echo "# id: 'frip'"
        echo "# section_name: 'Fraction of reads in peaks'"
        echo "# description: 'Alignments of the final BAM that overlap a called peak, over all alignments of that BAM.'"
        echo "# plot_type: 'table'"
        echo "# pconfig:"
        echo "#     id: 'frip_table'"
        echo "#     namespace: 'FRiP'"
        printf 'Peak set\\tFRiP\\treads_in_peaks\\ttotal_reads\\n'
        for peak_file in ${peaks}; do
            in_peaks=\$(bedtools intersect -u -a ${bam} -b "\$peak_file" | samtools view -c -)
            frip=\$(awk -v a="\$in_peaks" -v b="\$total" 'BEGIN { printf "%.4f", (b > 0) ? a / b : 0 }')
            printf '%s\\t%s\\t%s\\t%s\\n' "\$peak_file" "\$frip" "\$in_peaks" "\$total"
        done
    } > ${sample_id}.frip_mqc.tsv
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
    if (params.single_end) {
        error "This ATAC-seq workflow needs paired-end reads: Tn5 shifting, nucleosome-free selection, and BAMPE peak calling all depend on fragment length"
    }
    if (!genomeDefaults().containsKey(params.genome)) {
        error "Unsupported --genome '${params.genome}': expected one of ${genomeDefaults().keySet().join(', ')}"
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
    def defaults  = genomeDefaults()[params.genome]
    def blacklist = params.blacklist ?: defaults.blacklist

    ch_reads  = channel.fromFilePairs(params.reads, size: params.single_end ? 1 : 2, checkIfExists: true)
    ch_genome = channel.fromPath(params.bowtie2_index ?: defaults.index, type: 'dir', checkIfExists: true)
    ch_black  = channel.fromPath(blacklist)

    // Stage 1: QC and Trimming
    FASTQC(ch_reads)
    TRIM_GALORE(ch_reads)

    // Stage 2: Alignment (Bowtie2)
    BOWTIE2_ALIGN(TRIM_GALORE.out.trimmed, ch_genome.collect())

    // Stage 3: Tn5 Shift and Filtering
    MITO_FILTER(BOWTIE2_ALIGN.out.bam)
    MARK_DUPLICATES(MITO_FILTER.out.bam)
    TN5_SHIFT(MARK_DUPLICATES.out.bam)
    BLACKLIST_FILTER(TN5_SHIFT.out.bam, ch_black.collect())
    NFR_SELECTION(BLACKLIST_FILTER.out.bam)

    // Stage 4: Peak Calling on NFR fragments
    MACS2_CALLPEAK(NFR_SELECTION.out.nfr)

    // IDR (optional, with 2+ replicates)
    if (!params.skip_idr) {
        // IDR compares two replicates at a time, so every pair of samples is compared:
        // two samples give one comparison, three give three. One sample gives none.
        ch_idr_pairs = MACS2_CALLPEAK.out.peaks
            .toSortedList { a, b -> a[0] <=> b[0] }
            .flatMap { samples ->
                [samples, samples].combinations()
                    .findAll { pair -> pair[0][0] < pair[1][0] }
                    .collect { pair -> [pair[0][0], pair[0][1], pair[1][0], pair[1][1]] }
            }
        IDR_ANALYSIS(ch_idr_pairs)
    }

    // Stage 5: Signal Tracks and QC. FRiP counts all filtered fragments (the final BAM)
    // against the peaks called on the nucleosome-free fragments.
    SIGNAL_TRACKS(BLACKLIST_FILTER.out.bam)
    FRIP(BLACKLIST_FILTER.out.bam.join(MACS2_CALLPEAK.out.peaks))

    // MultiQC
    ch_multiqc = FASTQC.out.reports
        .mix(TRIM_GALORE.out.log)
        .mix(TRIM_GALORE.out.fastqc)
        .mix(BOWTIE2_ALIGN.out.log)
        .mix(MITO_FILTER.out.stats)
        .mix(MARK_DUPLICATES.out.metrics)
        .mix(BLACKLIST_FILTER.out.flagstat)
        .mix(FRIP.out.frip)
        .collect()
    MULTIQC(ch_multiqc)
}
