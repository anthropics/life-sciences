#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================================
// ENCODE DNase-seq Pipeline — FASTQ to Hotspots and Footprints
// Tools: BWA-MEM, Hotspot2, HINT (RGT)
// ============================================================================

params.reads                = null
params.bwa_index            = null
params.chrom_sizes          = null
params.hotspot_center_sites = null   // center_sites.starch made once per genome by extractCenterSites.sh
params.hotspot_mappable     = null   // optional: the mappable-regions BED used to make the center sites
params.blacklist            = null
params.outdir               = './results'
params.fdr                  = 0.05
params.skip_footprint       = false
params.organism             = 'hg38' // genome name registered in the RGT data directory (HINT)
params.rgt_data             = null   // populated RGT data directory; required unless --skip_footprint

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
        --fastqc \\
        ${reads[0]} ${reads[1]}
    """
}

process BWA_ALIGN {
    tag "${sample_id}"
    cpus 8
    memory { 16.GB * task.attempt }

    input:
    tuple val(sample_id), path(reads)
    path bwa_idx

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam

    script:
    // The index files are staged into the task directory, so use the prefix basename.
    def idx_base = file(params.bwa_index).name
    """
    bwa mem -t ${task.cpus} -M \\
        ${idx_base} \\
        ${reads[0]} ${reads[1]} \\
        | samtools view -@ 4 -bS - \\
        | samtools sort -@ 4 -o ${sample_id}.sorted.bam

    samtools index ${sample_id}.sorted.bam
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
    # Remove chrM and filter quality
    samtools idxstats ${bam} | \\
        awk '\$1 != "chrM" && \$1 != "*" {print \$1}' > chroms.txt

    samtools view -b -h -q 30 -F 1804 -f 2 \\
        ${bam} \$(cat chroms.txt | tr '\\n' ' ') \\
        | samtools sort -@ ${task.cpus} -o nuclear.bam

    # Mark and remove duplicates
    picard MarkDuplicates \\
        INPUT=nuclear.bam \\
        OUTPUT=dedup.bam \\
        METRICS_FILE=${sample_id}.dup_metrics.txt \\
        REMOVE_DUPLICATES=true \\
        VALIDATION_STRINGENCY=LENIENT \\
        ASSUME_SORTED=true

    # Remove blacklist regions
    bedtools intersect \\
        -a dedup.bam \\
        -b ${blacklist} \\
        -v \\
        > ${sample_id}.filtered.bam

    samtools index ${sample_id}.filtered.bam
    samtools flagstat ${sample_id}.filtered.bam > ${sample_id}.flagstat.txt

    rm nuclear.bam dedup.bam
    """
}

process HOTSPOT2 {
    tag "${sample_id}"
    publishDir "${params.outdir}/hotspots", mode: 'copy'
    cpus 4
    memory { 8.GB * task.attempt }

    input:
    tuple val(sample_id), path(bam), path(bai)
    path chrom_sizes
    path center_sites
    path mappable   // empty when --hotspot_mappable is not given

    output:
    tuple val(sample_id), path("${sample_id}.hotspots.fdr${params.fdr}.bed"), emit: hotspots
    tuple val(sample_id), path("${sample_id}.peaks.narrowPeak"), emit: peaks
    path("${sample_id}.SPOT.txt"), emit: spot
    path("${sample_id}.allcalls.bed"), emit: allcalls

    script:
    def mappable_opt = mappable ? "-M ${mappable}" : ''
    // The site-calling threshold (-F) may not be stricter than the hotspot threshold (-f).
    // hotspot2.sh names the hotspots, peaks and SPOT files after -f (HOTSPOT_FDR_THRESHOLD),
    // never after -F, so the names below always use params.fdr.
    def sitecall_fdr = Math.max(params.fdr as double, 0.05d)
    // hotspot2.sh names every output after the BAM basename
    def base = "hotspot2_out/${bam.baseName}"
    """
    # hotspot2.sh wants chromosome sizes as a sorted BED file with column 2 set to 0
    awk 'BEGIN {OFS="\\t"} {print \$1, 0, \$2}' ${chrom_sizes} | sort-bed - > chrom_sizes.bed

    # Usage: hotspot2.sh [options] in.bam outdir   (-c and -C are mandatory)
    hotspot2.sh \\
        -c chrom_sizes.bed \\
        -C ${center_sites} \\
        ${mappable_opt} \\
        -f ${params.fdr} \\
        -F ${sitecall_fdr} \\
        ${bam} hotspot2_out

    # Convert starch to BED
    unstarch ${base}.hotspots.fdr${params.fdr}.starch > ${sample_id}.hotspots.fdr${params.fdr}.bed
    unstarch ${base}.peaks.fdr${params.fdr}.narrowpeaks.starch > ${sample_id}.peaks.narrowPeak
    cp ${base}.SPOT.fdr${params.fdr}.txt ${sample_id}.SPOT.txt
    unstarch ${base}.allcalls.starch > ${sample_id}.allcalls.bed
    """
}

process SIGNAL_TRACK {
    tag "${sample_id}"
    publishDir "${params.outdir}/hotspots", mode: 'copy'
    cpus 2
    memory { 4.GB * task.attempt }

    input:
    tuple val(sample_id), path(bam), path(bai)
    path chrom_sizes

    output:
    path("${sample_id}.density.bw"), emit: bigwig

    script:
    """
    # Generate RPM-normalized signal
    total=\$(samtools view -c -F 1804 -f 2 ${bam})
    scale=\$(echo "scale=10; 1000000 / \$total" | bc)

    bedtools genomecov \\
        -ibam ${bam} \\
        -bg \\
        -pc \\
        -g ${chrom_sizes} \\
        | awk -v s=\$scale 'BEGIN{OFS="\\t"} {\$4=\$4*s; print}' \\
        | sort -k1,1 -k2,2n \\
        > ${sample_id}_rpm.bedGraph

    bedGraphToBigWig ${sample_id}_rpm.bedGraph ${chrom_sizes} ${sample_id}.density.bw
    """
}

process FOOTPRINTING {
    tag "${sample_id}"
    publishDir "${params.outdir}/footprints", mode: 'copy'
    cpus 4
    memory { 8.GB * task.attempt }

    input:
    tuple val(sample_id), path(bam), path(bai), path(peaks)
    path rgt_data

    output:
    path("${sample_id}.footprints.bed"), emit: footprints

    script:
    """
    # RGT reads its genome data from \$RGTDATA (default ~/rgtdata), which does not exist for the
    # unprivileged user the container runs as, so point it at the staged directory.
    export RGTDATA="\$PWD/${rgt_data}"

    rgt-hint footprinting \\
        --dnase-seq \\
        --paired-end \\
        --organism ${params.organism} \\
        --output-location fp_out/ \\
        --output-prefix ${sample_id} \\
        ${bam} \\
        ${peaks}

    cp fp_out/${sample_id}.bed ${sample_id}.footprints.bed
    """
}

process INSERT_SIZES {
    tag "${sample_id}"
    publishDir "${params.outdir}/qc", mode: 'copy'
    cpus 1
    memory { 4.GB * task.attempt }

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path("${sample_id}.insert_sizes.txt"), emit: metrics

    script:
    """
    # Picard CollectInsertSizeMetrics needs R for its mandatory histogram, which the image does
    # not ship. samtools stats reports the same insert-size distribution and MultiQC reads it.
    samtools stats ${bam} > ${sample_id}.insert_sizes.txt
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
    multiqc --title "ENCODE DNase-seq Pipeline" --filename multiqc_report --force .
    """
}

// ---- Workflow ----

workflow {
    // ---- Parameter validation ----
    if (!params.reads)                { error "Missing required parameter: --reads" }
    if (!params.bwa_index)            { error "Missing required parameter: --bwa_index" }
    if (!params.chrom_sizes)          { error "Missing required parameter: --chrom_sizes" }
    if (!params.hotspot_center_sites) { error "Missing required parameter: --hotspot_center_sites" }
    if (!params.blacklist)            { error "Missing required parameter: --blacklist" }
    if (!params.skip_footprint && !params.rgt_data) {
        error "Footprinting needs --rgt_data (an RGT data directory set up for --organism '${params.organism}'); pass it or use --skip_footprint"
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

    // ---- Channels ----
    ch_reads        = channel.fromFilePairs(params.reads, checkIfExists: true)
    ch_bwa_index    = channel.fromPath("${params.bwa_index}*", checkIfExists: true).collect()
    ch_chrom_sizes  = channel.fromPath(params.chrom_sizes, checkIfExists: true).collect()
    ch_center_sites = channel.fromPath(params.hotspot_center_sites, checkIfExists: true).collect()
    ch_blacklist    = channel.fromPath(params.blacklist, checkIfExists: true).collect()
    ch_mappable     = params.hotspot_mappable
        ? channel.fromPath(params.hotspot_mappable, checkIfExists: true).collect()
        : []

    FASTQC_RAW(ch_reads)
    TRIM_GALORE(ch_reads)
    BWA_ALIGN(TRIM_GALORE.out.trimmed, ch_bwa_index)
    FILTER_DEDUP(BWA_ALIGN.out.bam, ch_blacklist)
    HOTSPOT2(FILTER_DEDUP.out.bam, ch_chrom_sizes, ch_center_sites, ch_mappable)
    SIGNAL_TRACK(FILTER_DEDUP.out.bam, ch_chrom_sizes)
    INSERT_SIZES(FILTER_DEDUP.out.bam)

    if (!params.skip_footprint) {
        ch_rgt_data = channel.fromPath(params.rgt_data, type: 'dir', checkIfExists: true).collect()
        // join on sample_id so each BAM is footprinted against its own peaks
        FOOTPRINTING(FILTER_DEDUP.out.bam.join(HOTSPOT2.out.peaks), ch_rgt_data)
    }

    ch_multiqc = FASTQC_RAW.out.reports
        .mix(TRIM_GALORE.out.reports)
        .mix(TRIM_GALORE.out.fastqc)
        .mix(FILTER_DEDUP.out.flagstat)
        .mix(FILTER_DEDUP.out.dup_metrics)
        .mix(INSERT_SIZES.out.metrics)
        .mix(HOTSPOT2.out.spot)
        .collect()

    MULTIQC(ch_multiqc)
}
