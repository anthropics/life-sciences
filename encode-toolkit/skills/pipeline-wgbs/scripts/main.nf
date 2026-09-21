#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================================
// ENCODE WGBS Pipeline — FASTQ to bedMethyl
// Tools: Trim Galore, Bismark (Bowtie2), MethylDackel
// ============================================================================

params.reads         = null
params.genome_dir    = null    // Bismark genome folder (bismark_genome_preparation output + the .fa)
params.outdir        = './results'
params.min_coverage  = 5
params.merge_context = true    // merge the two strands of each CpG/CHG into one record
params.skip_dedup    = false

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
        --length 36 \\
        --cores ${task.cpus} \\
        --clip_R2 10 \\
        --three_prime_clip_R1 1 \\
        --fastqc \\
        ${reads[0]} ${reads[1]}
    """
}

process BISMARK_ALIGN {
    tag "${sample_id}"
    publishDir "${params.outdir}/bismark/alignments", mode: 'copy', pattern: '*_report.txt'
    cpus 8
    memory { 48.GB * task.attempt }

    input:
    tuple val(sample_id), path(reads)
    path genome_dir

    output:
    tuple val(sample_id), path("*.bam"), emit: bam
    path("*_report.txt"), emit: report

    script:
    """
    mkdir -p tmp
    bismark \\
        --genome ${genome_dir} \\
        --bowtie2 \\
        --parallel 4 \\
        --score_min L,0,-0.2 \\
        --no_mixed \\
        --no_discordant \\
        --maxins 1000 \\
        --temp_dir \$PWD/tmp \\
        -1 ${reads[0]} \\
        -2 ${reads[1]}
    """
}

process DEDUPLICATE {
    tag "${sample_id}"
    publishDir "${params.outdir}/bismark/dedup_reports", mode: 'copy', pattern: '*.txt'
    cpus 2
    memory { 16.GB * task.attempt }

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("*.deduplicated.bam"), emit: bam
    path("*.deduplication_report.txt"), emit: report

    when:
    !params.skip_dedup

    script:
    """
    deduplicate_bismark --bam --paired ${bam}
    """
}

process SAMTOOLS_SORT_INDEX {
    tag "${sample_id}"
    publishDir "${params.outdir}/bismark/alignments", mode: 'copy'
    cpus 4
    memory { 8.GB * task.attempt }

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam ${bam}
    samtools index ${sample_id}.sorted.bam
    """
}

process METHYLDACKEL_MBIAS {
    tag "${sample_id}"
    publishDir "${params.outdir}/bismark/mbias", mode: 'copy'
    cpus 2
    memory { 8.GB * task.attempt }

    input:
    tuple val(sample_id), path(bam), path(bai)
    path genome_dir

    output:
    path("*.svg"), emit: plots
    path("*.txt"), emit: report

    script:
    """
    GENOME_FA=\$(ls ${genome_dir}/*.fa | head -1)
    MethylDackel mbias \\
        --CHG --CHH \\
        \$GENOME_FA \\
        ${bam} \\
        ${sample_id}_mbias > ${sample_id}_mbias_report.txt 2>&1
    """
}

process METHYLDACKEL_EXTRACT {
    tag "${sample_id}"
    publishDir "${params.outdir}/bismark/methylation", mode: 'copy'
    cpus 4
    memory { 8.GB * task.attempt }

    input:
    tuple val(sample_id), path(bam), path(bai)
    path genome_dir

    output:
    tuple val(sample_id), path("*.bedGraph"), emit: bedgraph
    tuple val(sample_id), path("*.bedMethyl.gz"), emit: bedmethyl
    path("*.bedMethyl.gz.tbi"), emit: index

    script:
    // MethylDackel never counts both mates of an overlapping pair, so there is no overlap
    // switch. --mergeContext is a separate choice: per-CpG/CHG records instead of per-cytosine.
    // No read positions are excluded here: TRIM_GALORE already clips the end-repair bias at
    // the 5' end of read 2. Check bismark/mbias/ before trusting the calls.
    def merge = params.merge_context ? '--mergeContext' : ''
    """
    GENOME_FA=\$(ls ${genome_dir}/*.fa | head -1)
    MethylDackel extract \\
        ${merge} \\
        --CHG --CHH \\
        --opref ${sample_id} \\
        \$GENOME_FA \\
        ${bam}

    # Convert each context to ENCODE bedMethyl (https://www.encodeproject.org/data-standards/wgbs/):
    # column 5 = score, the read count capped at 1000; column 10 = coverage; column 11 = percent
    # methylated. Sites below --min_coverage are left out. MethylDackel bedGraphs start with a
    # "track" header line and carry no strand, so the header is skipped and strand is ".".
    for context in CpG CHG CHH; do
        awk -v min_cov=${params.min_coverage} 'BEGIN {OFS="\\t"} !/^track/ {
            cov = \$5 + \$6
            if (cov == 0 || cov < min_cov) next
            score = (cov > 1000) ? 1000 : cov
            print \$1, \$2, \$3, ".", score, ".", \$2, \$3, "0,0,0", cov, int((\$5 / cov) * 100 + 0.5)
        }' ${sample_id}_\${context}.bedGraph \\
            | sort -k1,1 -k2,2n \\
            | bgzip > ${sample_id}.\${context}.bedMethyl.gz

        tabix -p bed ${sample_id}.\${context}.bedMethyl.gz
    done
    """
}

process COVERAGE_STATS {
    tag "${sample_id}"
    publishDir "${params.outdir}/coverage", mode: 'copy'
    cpus 2
    memory { 4.GB * task.attempt }

    input:
    tuple val(sample_id), path(bedgraph)

    output:
    path("${sample_id}.coverage_stats.txt"), emit: stats

    script:
    """
    awk -v min_cov=${params.min_coverage} '!/^track/ {
        cov = \$5 + \$6; sum += cov; n++;
        if (cov >= 5)       c5++;
        if (cov >= 10)      c10++;
        if (cov >= min_cov) cmin++
    } END {
        if (n == 0) { print "Covered CpGs: 0"; exit }
        printf "Covered CpGs (>=1x): %d\\n", n;
        printf "Mean coverage of covered CpGs: %.1f\\n", sum/n;
        printf "Covered CpGs >=5x: %d (%.1f%%)\\n", c5, c5/n*100;
        printf "Covered CpGs >=10x: %d (%.1f%%)\\n", c10, c10/n*100;
        printf "Covered CpGs >=%dx (--min_coverage, kept in bedMethyl): %d (%.1f%%)\\n", min_cov, cmin, cmin/n*100
    }' ${sample_id}_CpG.bedGraph > ${sample_id}.coverage_stats.txt
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
    multiqc --title "ENCODE WGBS Pipeline" --filename multiqc_report --force .
    """
}

// ---- Workflow ----

workflow {
    // ---- Parameter validation ----
    if (!params.reads)      { error "Missing required parameter: --reads" }
    if (!params.genome_dir) { error "Missing required parameter: --genome_dir" }

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
    ch_reads  = channel.fromFilePairs(params.reads, checkIfExists: true)
    ch_genome = channel.fromPath(params.genome_dir, type: 'dir', checkIfExists: true).collect()

    FASTQC_RAW(ch_reads)
    TRIM_GALORE(ch_reads)

    BISMARK_ALIGN(TRIM_GALORE.out.trimmed, ch_genome)

    if (!params.skip_dedup) {
        DEDUPLICATE(BISMARK_ALIGN.out.bam)
        SAMTOOLS_SORT_INDEX(DEDUPLICATE.out.bam)
        ch_dedup_report = DEDUPLICATE.out.report
    } else {
        SAMTOOLS_SORT_INDEX(BISMARK_ALIGN.out.bam)
        ch_dedup_report = channel.empty()
    }

    METHYLDACKEL_MBIAS(SAMTOOLS_SORT_INDEX.out.bam, ch_genome)
    METHYLDACKEL_EXTRACT(SAMTOOLS_SORT_INDEX.out.bam, ch_genome)
    COVERAGE_STATS(METHYLDACKEL_EXTRACT.out.bedgraph)

    ch_multiqc = FASTQC_RAW.out.reports
        .mix(TRIM_GALORE.out.reports)
        .mix(TRIM_GALORE.out.fastqc)
        .mix(BISMARK_ALIGN.out.report)
        .mix(ch_dedup_report)
        .collect()

    MULTIQC(ch_multiqc)
}
