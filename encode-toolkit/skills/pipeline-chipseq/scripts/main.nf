#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ENCODE ChIP-seq Pipeline — Nextflow DSL2
// FASTQ -> QC -> Align -> Filter -> Peaks -> IDR -> Signal

params.reads       = null
params.control     = null
params.genome      = 'GRCh38'
params.peak_type   = 'narrow'
params.outdir      = 'results'
params.single_end  = false
params.skip_idr    = false
params.blacklist   = null   // defaults to the ENCODE blacklist v2 for --genome
params.chrom_sizes = null   // required: two-column chromosome sizes for the bigWig tracks
params.bwa_index   = null   // directory holding <genome>.fa and its BWA index; default ./<genome>_index

// Genome-specific defaults. Kept in a function because scripts that declare processes
// cannot also hold top-level variables.
def genomeDefaults() {
    return [
    'GRCh38': [
        blacklist:  'https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz',
        gsize:      'hs'
    ],
    'mm10': [
        blacklist:  'https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz',
        gsize:      'mm'
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

process BWA_MEM {
    tag "$sample_id"
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path(genome_index)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam
    path("${sample_id}.flagstat.txt"),                                              emit: flagstat

    script:
    def input_reads = params.single_end ? "${reads}" : "${reads[0]} ${reads[1]}"
    """
    bwa mem -t ${task.cpus} -M ${genome_index}/${params.genome}.fa ${input_reads} | \\
      samtools view -@ ${task.cpus} -bS -q 30 - | \\
      samtools sort -@ ${task.cpus} -m 2G -o ${sample_id}.bam -
    samtools index ${sample_id}.bam
    samtools flagstat ${sample_id}.bam > ${sample_id}.flagstat.txt
    """
}

process FILTER_SORT {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.filtered.bam"), emit: bam

    script:
    def flags = params.single_end ? '-F 1028' : '-F 1804'
    """
    samtools view -@ ${task.cpus} -b ${flags} -q 30 ${bam} | \\
      samtools sort -@ ${task.cpus} -o ${sample_id}.filtered.bam -
    """
}

process MARK_DUPLICATES {
    tag "$sample_id"
    publishDir "${params.outdir}/filtered", mode: 'copy', pattern: '*_metrics.txt'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), emit: bam
    path("${sample_id}.dup_metrics.txt"),                  emit: metrics

    script:
    """
    picard MarkDuplicates \\
      INPUT=${bam} OUTPUT=${sample_id}.dedup.bam \\
      METRICS_FILE=${sample_id}.dup_metrics.txt \\
      REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
    """
}

process BLACKLIST_FILTER {
    tag "$sample_id"
    publishDir "${params.outdir}/filtered", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
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

process MACS2_CALLPEAK {
    tag "$sample_id"
    publishDir "${params.outdir}/peaks/${params.peak_type}", mode: 'copy'

    input:
    tuple val(sample_id), path(treatment_bam), path(treatment_bai)
    path(control_bams)   // all control BAMs (MACS2 pools them); empty when no --control is given

    output:
    tuple val(sample_id), path("${sample_id}*Peak"),    emit: peaks
    tuple val(sample_id), path("${sample_id}_peaks.{narrowPeak,broadPeak}"), emit: main_peaks
    tuple val(sample_id), path("${sample_id}*.bdg"),    emit: bdg
    path("${sample_id}*.xls"),                           emit: xls
    path("${sample_id}*_summits.bed"),                   emit: summits, optional: true   // narrow peaks only

    script:
    def format_flag  = params.single_end ? 'BAM' : 'BAMPE'
    def broad_flags  = params.peak_type == 'broad' ? '--broad --broad-cutoff 0.1' : '--call-summits'
    def control_flag = control_bams ? "-c ${control_bams}" : ''
    def gsize        = genomeDefaults()[params.genome].gsize
    """
    macs2 callpeak \\
      -t ${treatment_bam} ${control_flag} \\
      -f ${format_flag} -g ${gsize} -n ${sample_id} \\
      --qvalue 0.05 --nomodel --keep-dup all \\
      ${broad_flags} -B
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
    tuple val(sample_id), path(bdg_files)
    path(chrom_sizes)

    output:
    path("${sample_id}.fc.bw"),   emit: fc_bw
    path("${sample_id}.pval.bw"), emit: pval_bw

    script:
    """
    macs2 bdgcmp -t ${sample_id}_treat_pileup.bdg -c ${sample_id}_control_lambda.bdg -o fc.bdg -m FE
    sort -k1,1 -k2,2n fc.bdg > fc.sorted.bdg
    bedGraphToBigWig fc.sorted.bdg ${chrom_sizes} ${sample_id}.fc.bw

    macs2 bdgcmp -t ${sample_id}_treat_pileup.bdg -c ${sample_id}_control_lambda.bdg -o pval.bdg -m ppois
    sort -k1,1 -k2,2n pval.bdg > pval.sorted.bdg
    bedGraphToBigWig pval.sorted.bdg ${chrom_sizes} ${sample_id}.pval.bw
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
    if (!params.reads)       { error "Missing required parameter: --reads" }
    if (!params.chrom_sizes) { error "Missing required parameter: --chrom_sizes" }
    if (!genomeDefaults().containsKey(params.genome)) {
        error "Unsupported --genome '${params.genome}': expected one of ${genomeDefaults().keySet().join(', ')}"
    }
    if (!(params.peak_type in ['narrow', 'broad'])) {
        error "Invalid --peak_type '${params.peak_type}': expected 'narrow' or 'broad'"
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
    def n_reads   = params.single_end ? 1 : 2
    def blacklist = params.blacklist ?: genomeDefaults()[params.genome].blacklist

    ch_treatment = channel.fromFilePairs(params.reads, size: n_reads, checkIfExists: true)
    ch_genome    = channel.fromPath(params.bwa_index ?: "${params.genome}_index", type: 'dir', checkIfExists: true).collect()
    ch_black     = channel.fromPath(blacklist).collect()
    ch_chromsz   = channel.fromPath(params.chrom_sizes, checkIfExists: true).collect()

    // A process can be called only once per workflow, so control libraries travel through the
    // same steps as the ChIP samples. They are tagged with a CONTROL_ prefix and split off
    // again before peak calling.
    ch_control = params.control
        ? channel.fromFilePairs(params.control, size: n_reads, checkIfExists: true)
              .map { sample_id, reads -> ["CONTROL_${sample_id}", reads] }
        : channel.empty()
    ch_reads = ch_treatment.mix(ch_control)

    // Stage 1: QC and Trimming
    FASTQC(ch_reads)
    TRIM_GALORE(ch_reads)

    // Stage 2: Alignment
    BWA_MEM(TRIM_GALORE.out.trimmed, ch_genome)

    // Stage 3: Filtering
    FILTER_SORT(BWA_MEM.out.bam)
    MARK_DUPLICATES(FILTER_SORT.out.bam)
    BLACKLIST_FILTER(MARK_DUPLICATES.out.bam, ch_black)

    ch_final = BLACKLIST_FILTER.out.bam.branch { sample_id, _bam, _bai ->
        control: sample_id.startsWith('CONTROL_')
        treatment: true
    }
    ch_control_bams = params.control
        ? ch_final.control.map { _sample_id, bam, _bai -> bam }.collect()
        : []

    // Stage 4: Peak Calling (controls are pooled, as MACS2 does for multiple -c files)
    MACS2_CALLPEAK(ch_final.treatment, ch_control_bams)

    // IDR (optional, narrow peaks with 2+ replicates)
    if (!params.skip_idr && params.peak_type == 'narrow') {
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

    // Stage 5: Signal Tracks and FRiP
    SIGNAL_TRACKS(MACS2_CALLPEAK.out.bdg, ch_chromsz)
    FRIP(ch_final.treatment.join(MACS2_CALLPEAK.out.main_peaks))

    // MultiQC
    ch_multiqc = FASTQC.out.reports
        .mix(TRIM_GALORE.out.log)
        .mix(TRIM_GALORE.out.fastqc)
        .mix(BWA_MEM.out.flagstat)
        .mix(MARK_DUPLICATES.out.metrics)
        .mix(BLACKLIST_FILTER.out.flagstat)
        .mix(FRIP.out.frip)
        .collect()
    MULTIQC(ch_multiqc)
}
