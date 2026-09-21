#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================================
// ENCODE CUT&RUN Pipeline — FASTQ to Peaks and Signal Tracks
// Tools: Bowtie2, SEACR, MACS2, deepTools
// ============================================================================

params.reads          = null
params.bowtie2_index  = null
params.spikein_index  = null
params.chrom_sizes    = null
params.blacklist      = null
params.outdir         = './results'
params.seacr_mode     = 'stringent'  // 'stringent', 'relaxed', or 'both'
params.seacr_norm     = 'norm'       // 'norm' or 'non'; only used when --control is given
params.seacr_threshold = 0.01        // top fraction of signal kept when no --control is given
params.control        = null         // IgG control BAM (filtered, deduplicated)
params.peak_caller    = 'seacr'      // 'seacr', 'macs2', or 'both'
params.macs2_gsize    = 'hs'         // MACS2 effective genome size ('hs', 'mm', or a number)
params.skip_spikein   = false

// ---- Processes ----

process FASTQC_RAW {
    tag "${sample_id}"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    cpus 2
    memory { 4.GB * task.attempt }

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*.{html,zip}"), emit: reports

    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${reads}
    """
}

process TRIM_GALORE {
    tag "${sample_id}"
    publishDir "${params.outdir}/trim_galore", mode: 'copy'
    cpus 4
    memory { 4.GB * task.attempt }

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_val_{1,2}.fq.gz"), emit: trimmed
    path("*_trimming_report.txt"), emit: reports
    path("*_fastqc.{html,zip}"), emit: fastqc

    script:
    """
    trim_galore \\
        --paired \\
        --quality 20 \\
        --phred33 \\
        --length 20 \\
        --cores ${task.cpus} \\
        --nextera \\
        --fastqc \\
        ${reads[0]} ${reads[1]}
    """
}

process BOWTIE2_ALIGN {
    tag "${sample_id}"
    cpus 8
    memory { 8.GB * task.attempt }

    input:
    tuple val(sample_id), path(reads)
    path bt2_idx

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam
    path("${sample_id}.bowtie2.log"), emit: log

    script:
    // The index files are staged into the task directory, so use the prefix basename.
    def idx_prefix = file(params.bowtie2_index).name
    """
    bowtie2 \\
        --very-sensitive \\
        --no-mixed \\
        --no-discordant \\
        --dovetail \\
        --phred33 \\
        -I 10 -X 700 \\
        --threads ${task.cpus} \\
        -x ${idx_prefix} \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        2> ${sample_id}.bowtie2.log \\
        | samtools view -@ 4 -bS - \\
        | samtools sort -@ 4 -o ${sample_id}.sorted.bam

    samtools index ${sample_id}.sorted.bam
    """
}

process SPIKEIN_ALIGN {
    tag "${sample_id}"
    publishDir "${params.outdir}/spikein", mode: 'copy'
    cpus 4
    memory { 4.GB * task.attempt }

    input:
    tuple val(sample_id), path(bam), path(bai)
    path spikein_idx

    output:
    path("${sample_id}.spikein_counts.txt"), emit: counts

    script:
    def idx_prefix = file(params.spikein_index).name
    """
    # Extract read pairs that did not map to the primary genome
    samtools view -b -f 12 -F 256 ${bam} | samtools sort -n -@ 2 -o unmapped.bam
    bedtools bamtofastq -i unmapped.bam -fq unmap_R1.fq -fq2 unmap_R2.fq

    # Align to the spike-in genome
    bowtie2 \\
        --very-sensitive \\
        --no-mixed --no-discordant --dovetail \\
        -I 10 -X 700 \\
        --threads ${task.cpus} \\
        -x ${idx_prefix} \\
        -1 unmap_R1.fq -2 unmap_R2.fq \\
        2> spikein.log \\
        | samtools view -bS -q 10 -F 1804 -f 2 - \\
        | samtools sort -o spikein.bam

    spikein_count=\$(samtools view -c spikein.bam)
    printf '%s\\t%s\\n' "${sample_id}" "\${spikein_count}" > ${sample_id}.spikein_counts.txt
    """
}

process FILTER_DEDUP {
    tag "${sample_id}"
    publishDir "${params.outdir}/alignment", mode: 'copy'
    cpus 4
    memory { 8.GB * task.attempt }

    input:
    tuple val(sample_id), path(bam), path(bai)
    path blacklist

    output:
    tuple val(sample_id), path("${sample_id}.filtered.bam"), path("${sample_id}.filtered.bam.bai"), emit: bam
    path("${sample_id}.dup_metrics.txt"), emit: dup_metrics
    path("${sample_id}.flagstat.txt"), emit: flagstat

    script:
    """
    # Quality filter
    samtools view -b -h -q 10 -F 1804 -f 2 ${bam} \\
        | samtools sort -@ ${task.cpus} -o qfilt.bam

    # Mark and remove duplicates
    picard MarkDuplicates \\
        INPUT=qfilt.bam \\
        OUTPUT=dedup.bam \\
        METRICS_FILE=${sample_id}.dup_metrics.txt \\
        REMOVE_DUPLICATES=true \\
        VALIDATION_STRINGENCY=LENIENT \\
        ASSUME_SORTED=true

    # Remove blacklist regions
    bedtools intersect -a dedup.bam -b ${blacklist} -v \\
        > ${sample_id}.filtered.bam

    samtools index ${sample_id}.filtered.bam
    samtools flagstat ${sample_id}.filtered.bam > ${sample_id}.flagstat.txt

    rm qfilt.bam dedup.bam
    """
}

process COMPUTE_SCALE_FACTOR {
    publishDir "${params.outdir}/spikein", mode: 'copy'
    cpus 1
    memory { 1.GB * task.attempt }

    input:
    path(counts)

    output:
    path("scale_factors.txt"), emit: factors

    script:
    """
    cat ${counts} > all_counts.txt

    # Scale every sample to the smallest non-zero spike-in count (factor = min / count).
    # A sample with no spike-in reads cannot be calibrated and is left unscaled (factor 1).
    min_count=\$(awk -F'\\t' '\$2 > 0 {print \$2}' all_counts.txt | sort -n | head -1)
    awk -F'\\t' -v min="\${min_count:-0}" 'BEGIN {OFS="\\t"} {
        factor = (\$2 > 0 && min > 0) ? min / \$2 : 1
        print \$1, \$2, factor
    }' all_counts.txt > scale_factors.txt
    """
}

process FRAGMENT_BEDGRAPH {
    tag "${sample_id}"
    // Publish each sample's fragment BED; the IgG control's copy is an intermediate
    publishDir "${params.outdir}/signal", mode: 'copy', pattern: '*.fragments.bed',
        saveAs: { name -> sample_id == '__control__' ? null : name }
    cpus 2
    memory { 4.GB * task.attempt }

    input:
    tuple val(sample_id), path(bam)
    path chrom_sizes

    output:
    tuple val(sample_id), path("${sample_id}_fragments.bedGraph"), emit: bedgraph
    path("${sample_id}.fragments.bed"), emit: fragments

    script:
    """
    # bamtobed -bedpe needs mates on adjacent lines, so name-sort first
    samtools sort -n -@ ${task.cpus} -o namesorted.bam ${bam}

    # Keep properly paired fragments on one chromosome and shorter than 1 kb
    bedtools bamtobed -bedpe -i namesorted.bam \\
        | awk 'BEGIN {OFS="\\t"} \$1 == \$4 && \$6 - \$2 < 1000 {print \$1, \$2, \$6}' \\
        | sort -k1,1 -k2,2n -k3,3n \\
        > ${sample_id}.fragments.bed

    bedtools genomecov -bg -i ${sample_id}.fragments.bed -g ${chrom_sizes} \\
        > ${sample_id}_fragments.bedGraph

    rm namesorted.bam
    """
}

process SEACR_PEAKS {
    tag "${sample_id}"
    publishDir "${params.outdir}/peaks", mode: 'copy'
    cpus 2
    memory { 4.GB * task.attempt }

    input:
    tuple val(sample_id), path(bedgraph)
    path control_bedgraph   // empty when no --control is given

    output:
    tuple val(sample_id), path("${sample_id}.seacr.*.bed"), emit: peaks

    script:
    // SEACR takes either a control bedGraph or a numeric threshold as its second argument.
    // Normalization to the control only applies when a control bedGraph is supplied.
    def ctrl  = control_bedgraph ? "${control_bedgraph}" : "${params.seacr_threshold}"
    def norm  = control_bedgraph ? params.seacr_norm : 'non'
    def modes = params.seacr_mode == 'both' ? ['stringent', 'relaxed'] : [params.seacr_mode]
    def calls = modes.collect { mode ->
        "SEACR_1.3.sh ${bedgraph} ${ctrl} ${norm} ${mode} ${sample_id}.seacr"
    }.join('\n    ')
    """
    ${calls}
    """
}

process MACS2_PEAKS {
    tag "${sample_id}"
    publishDir "${params.outdir}/peaks", mode: 'copy'
    cpus 2
    memory { 4.GB * task.attempt }

    input:
    tuple val(sample_id), path(bam), path(bai)
    path control_bam   // empty when no --control is given

    output:
    tuple val(sample_id), path("${sample_id}.macs2_peaks.narrowPeak"), emit: peaks

    script:
    def ctrl_flag = control_bam ? "-c ${control_bam}" : ""
    """
    macs2 callpeak \\
        -t ${bam} \\
        ${ctrl_flag} \\
        -f BAMPE \\
        -g ${params.macs2_gsize} \\
        -n ${sample_id}.macs2 \\
        --nomodel \\
        --keep-dup all \\
        -q 0.05 \\
        --outdir .
    """
}

process SIGNAL_TRACK {
    tag "${sample_id}"
    publishDir "${params.outdir}/signal", mode: 'copy'
    cpus 4
    memory { 8.GB * task.attempt }

    input:
    tuple val(sample_id), path(bam), path(bai), val(scale_factor)

    output:
    path("${sample_id}.normalized.bw"), emit: bigwig

    script:
    // With a spike-in scale factor the track is spike-in calibrated; otherwise fall back to RPKM.
    def normalization = scale_factor
        ? "--scaleFactor ${scale_factor} --normalizeUsing None"
        : "--normalizeUsing RPKM"
    """
    bamCoverage \\
        --bam ${bam} \\
        --outFileName ${sample_id}.normalized.bw \\
        --binSize 10 \\
        ${normalization} \\
        --extendReads \\
        --numberOfProcessors ${task.cpus}
    """
}

process FRAGMENT_SIZES {
    tag "${sample_id}"
    publishDir "${params.outdir}/qc", mode: 'copy'
    cpus 1
    memory { 4.GB * task.attempt }

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path("${sample_id}.fragment_sizes.txt"), emit: sizes

    script:
    """
    samtools view -f 2 -F 1804 ${bam} | \\
        awk '{if(\$9 > 0 && \$9 < 1000) print \$9}' | \\
        sort -n | uniq -c | \\
        awk '{print \$2, \$1}' > ${sample_id}.fragment_sizes.txt
    """
}

process FRIP {
    tag "$sample_id"
    publishDir "${params.outdir}/qc", mode: 'copy'
    cpus 2
    memory { 4.GB * task.attempt }

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
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    cpus 1
    memory { 4.GB * task.attempt }

    input:
    path('*')

    output:
    path("multiqc_report.html"), emit: report

    script:
    """
    multiqc --title "ENCODE CUT&RUN Pipeline" --filename multiqc_report --force .
    """
}

// ---- Workflow ----

workflow {
    // ---- Parameter validation ----
    if (!params.reads)         { error "Missing required parameter: --reads" }
    if (!params.bowtie2_index) { error "Missing required parameter: --bowtie2_index" }
    if (!params.chrom_sizes)   { error "Missing required parameter: --chrom_sizes" }
    if (!params.blacklist)     { error "Missing required parameter: --blacklist" }
    if (!(params.seacr_mode in ['stringent', 'relaxed', 'both'])) {
        error "Invalid --seacr_mode '${params.seacr_mode}': expected 'stringent', 'relaxed', or 'both'"
    }
    if (!(params.seacr_norm in ['norm', 'non'])) {
        error "Invalid --seacr_norm '${params.seacr_norm}': expected 'norm' or 'non'"
    }
    if (!(params.peak_caller in ['seacr', 'macs2', 'both'])) {
        error "Invalid --peak_caller '${params.peak_caller}': expected 'seacr', 'macs2', or 'both'"
    }

    def use_spikein = !params.skip_spikein && params.spikein_index
    def use_seacr   = params.peak_caller in ['seacr', 'both']
    def use_macs2   = params.peak_caller in ['macs2', 'both']

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

    // ---- Channels ----
    ch_reads       = channel.fromFilePairs(params.reads, checkIfExists: true)
    ch_bt2_index   = channel.fromPath("${params.bowtie2_index}*", checkIfExists: true).collect()
    ch_chrom_sizes = channel.fromPath(params.chrom_sizes, checkIfExists: true).collect()
    ch_blacklist   = channel.fromPath(params.blacklist, checkIfExists: true).collect()
    ch_control_bam = params.control
        ? channel.fromPath(params.control, checkIfExists: true).collect()
        : []

    // ---- Alignment and filtering ----
    FASTQC_RAW(ch_reads)
    TRIM_GALORE(ch_reads)
    BOWTIE2_ALIGN(TRIM_GALORE.out.trimmed, ch_bt2_index)
    FILTER_DEDUP(BOWTIE2_ALIGN.out.bam, ch_blacklist)

    // ---- Spike-in calibration ----
    if (use_spikein) {
        ch_spikein_index = channel.fromPath("${params.spikein_index}*", checkIfExists: true).collect()
        SPIKEIN_ALIGN(BOWTIE2_ALIGN.out.bam, ch_spikein_index)
        COMPUTE_SCALE_FACTOR(SPIKEIN_ALIGN.out.counts.collect())

        // scale_factors.txt columns: sample, spike-in count, scale factor
        ch_factors = COMPUTE_SCALE_FACTOR.out.factors
            .splitCsv(sep: '\t')
            .map { row -> [row[0], row[2]] }
        ch_signal_in = FILTER_DEDUP.out.bam.join(ch_factors)
    } else {
        ch_signal_in = FILTER_DEDUP.out.bam.map { sample_id, bam, bai -> [sample_id, bam, bai, ''] }
    }

    // ---- Fragment bedGraphs (samples, plus the IgG control when given) ----
    ch_fragment_in = FILTER_DEDUP.out.bam.map { sample_id, bam, _bai -> [sample_id, bam] }
    if (params.control) {
        ch_fragment_in = ch_fragment_in.mix(ch_control_bam.map { bams -> ['__control__', bams[0]] })
    }
    FRAGMENT_BEDGRAPH(ch_fragment_in, ch_chrom_sizes)

    ch_bedgraph = FRAGMENT_BEDGRAPH.out.bedgraph.branch { sample_id, _bedgraph ->
        control: sample_id == '__control__'
        sample: true
    }
    ch_control_bedgraph = params.control
        ? ch_bedgraph.control.map { _sample_id, bedgraph -> bedgraph }.collect()
        : []

    // ---- Peak calling ----
    if (use_seacr) {
        SEACR_PEAKS(ch_bedgraph.sample, ch_control_bedgraph)
    }
    if (use_macs2) {
        MACS2_PEAKS(FILTER_DEDUP.out.bam, ch_control_bam)
    }

    // ---- Signal and QC ----
    SIGNAL_TRACK(ch_signal_in)
    FRAGMENT_SIZES(FILTER_DEDUP.out.bam)

    // FRiP for every peak set that was called for a sample (SEACR modes and/or MACS2)
    ch_sample_peaks = (use_seacr ? SEACR_PEAKS.out.peaks : channel.empty())
        .mix(use_macs2 ? MACS2_PEAKS.out.peaks : channel.empty())
        .groupTuple()
        .map { sample_id, peak_sets -> [sample_id, peak_sets.flatten()] }
    FRIP(FILTER_DEDUP.out.bam.join(ch_sample_peaks))

    ch_multiqc = FASTQC_RAW.out.reports
        .mix(TRIM_GALORE.out.reports)
        .mix(TRIM_GALORE.out.fastqc)
        .mix(BOWTIE2_ALIGN.out.log)
        .mix(FILTER_DEDUP.out.flagstat)
        .mix(FILTER_DEDUP.out.dup_metrics)
        .mix(FRIP.out.frip)
        .collect()

    MULTIQC(ch_multiqc)
}
