#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================================
// ENCODE Hi-C Pipeline — FASTQ to Contact Matrices and Loops
// Tools: BWA-MEM, pairtools, Juicer, cooler, HiCCUPS
// ============================================================================

// Contacts are processed at read-pair level with pairtools, which does not depend on the
// restriction enzyme, so the same workflow applies to MboI/DpnII, HindIII, Arima, and Micro-C.

params.reads            = null
params.bwa_index        = null
params.chrom_sizes      = null
params.outdir           = './results'
params.resolutions      = '1000,5000,10000,25000,50000,100000,250000,500000,1000000'   // smallest = base bin
params.hiccups_resolutions = '5000,10000,25000'   // loop-calling resolutions; each must also be in --resolutions
params.min_mapq         = 30
params.assembly         = 'hg38'
params.hiccups_gpu      = false   // HiCCUPS runs its CPU mode unless an NVIDIA GPU and CUDA are available

// "--resolutions 5000" arrives as a number and "1000, 5000" carries spaces, so every use of the
// two resolution lists goes through this parser and renders the normalized list.
def parseResolutions(value) {
    return value.toString().tokenize(',').collect { r -> r.trim() as long }
}

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

process BWA_ALIGN {
    tag "${sample_id}"
    publishDir "${params.outdir}/alignment", mode: 'copy'
    cpus 8
    memory { 16.GB * task.attempt }

    input:
    tuple val(sample_id), path(reads)
    path bwa_idx

    output:
    tuple val(sample_id), path("${sample_id}.paired.bam"), emit: bam

    script:
    // The index files are staged into the task directory, so use the prefix basename.
    def idx_base = file(params.bwa_index).name
    """
    bwa mem -t ${task.cpus} -SP5M \\
        ${idx_base} \\
        ${reads[0]} ${reads[1]} \\
        | samtools view -@ 4 -bhS - \\
        > ${sample_id}.paired.bam
    """
}

process PAIRTOOLS_PARSE_SORT {
    tag "${sample_id}"
    // The parse statistics hold the pair-type breakdown (UU, NU, MM, WW, ...) used for QC
    publishDir "${params.outdir}/pairs", mode: 'copy', pattern: '*.parse_stats.txt'
    cpus 4
    memory { 16.GB * task.attempt }

    input:
    tuple val(sample_id), path(bam)
    path chrom_sizes

    output:
    tuple val(sample_id), path("${sample_id}.sorted.pairs.gz"), emit: pairs
    path("${sample_id}.parse_stats.txt"), emit: stats

    script:
    """
    mkdir -p tmp   # pairtools sort hands --tmpdir to GNU sort, which does not create it

    pairtools parse \\
        --chroms-path ${chrom_sizes} \\
        --min-mapq ${params.min_mapq} \\
        --walks-policy mask \\
        --max-inter-align-gap 30 \\
        --nproc-in ${task.cpus} \\
        --nproc-out ${task.cpus} \\
        --output-stats ${sample_id}.parse_stats.txt \\
        ${bam} \\
        | pairtools sort \\
            --nproc ${task.cpus} \\
            --tmpdir \$PWD/tmp \\
            -o ${sample_id}.sorted.pairs.gz
    """
}

process PAIRTOOLS_DEDUP {
    tag "${sample_id}"
    publishDir "${params.outdir}/pairs", mode: 'copy'
    cpus 4
    memory { 16.GB * task.attempt }

    input:
    tuple val(sample_id), path(pairs)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.pairs.gz"), emit: pairs
    path("${sample_id}.dedup_stats.txt"), emit: stats

    script:
    """
    pairtools dedup \\
        --nproc-in ${task.cpus} \\
        --nproc-out ${task.cpus} \\
        --mark-dups \\
        --output-stats ${sample_id}.dedup_stats.txt \\
        -o ${sample_id}.dedup.pairs.gz \\
        ${pairs}
    """
}

process PAIRTOOLS_SELECT {
    tag "${sample_id}"
    cpus 2
    memory { 8.GB * task.attempt }

    input:
    tuple val(sample_id), path(pairs)

    output:
    tuple val(sample_id), path("${sample_id}.valid.pairs.gz"), emit: pairs
    path("${sample_id}.pair_stats.txt"), emit: stats

    script:
    """
    pairtools select \\
        '(pair_type == "UU")' \\
        ${pairs} \\
        -o ${sample_id}.valid.pairs.gz

    pairtools stats \\
        ${sample_id}.valid.pairs.gz \\
        -o ${sample_id}.pair_stats.txt
    """
}

process JUICER_HIC {
    tag "${sample_id}"
    publishDir "${params.outdir}/matrices", mode: 'copy'
    cpus 4
    memory { 64.GB * task.attempt }

    input:
    tuple val(sample_id), path(pairs)
    path chrom_sizes

    output:
    tuple val(sample_id), path("${sample_id}.hic"), emit: hic

    script:
    // The JVM needs memory beyond its heap, so the heap gets 85% of the task's allocation
    def heap_gb     = Math.max(1, (task.memory.toGiga() * 0.85) as int)
    def resolutions = parseResolutions(params.resolutions).join(',')
    """
    # Convert to Juicer short format: str1 chr1 pos1 frag1 str2 chr2 pos2 frag2.
    # No restriction-site file is used, so the fragment fields carry the dummy values
    # 0 and 1 that Juicer's `pre` command documents for this case.
    zcat ${pairs} | awk 'BEGIN{OFS="\\t"} !/^#/ {
        s1 = (\$6 == "+") ? 0 : 16;
        s2 = (\$7 == "+") ? 0 : 16;
        print s1, \$2, \$3, 0, s2, \$4, \$5, 1
    }' > juicer_medium.txt

    # --threads speeds up the normalization step. Building the matrix itself stays
    # single-threaded (pre says so on stderr) because no --mndindex is supplied.
    java -Xmx${heap_gb}g -jar /opt/juicer_tools.jar pre \\
        --threads ${task.cpus} \\
        -r ${resolutions} \\
        -k KR,VC,VC_SQRT \\
        juicer_medium.txt \\
        ${sample_id}.hic \\
        ${chrom_sizes}

    rm juicer_medium.txt
    """
}

process COOLER_MCOOL {
    tag "${sample_id}"
    publishDir "${params.outdir}/matrices", mode: 'copy'
    cpus 4
    memory { 16.GB * task.attempt }

    input:
    tuple val(sample_id), path(pairs)
    path chrom_sizes

    output:
    tuple val(sample_id), path("${sample_id}.mcool"), emit: mcool

    script:
    // Bin at the smallest requested resolution; zoomify coarsens that into the others
    def resolutions = parseResolutions(params.resolutions)
    def base_bin    = resolutions.min()
    """
    cooler cload pairs \\
        --chrom1 2 --pos1 3 --chrom2 4 --pos2 5 \\
        --assembly ${params.assembly} \\
        ${chrom_sizes}:${base_bin} \\
        ${pairs} \\
        ${sample_id}_base.cool

    cooler zoomify \\
        --balance \\
        --resolutions ${resolutions.join(',')} \\
        --nproc ${task.cpus} \\
        ${sample_id}_base.cool \\
        -o ${sample_id}.mcool
    """
}

process HICCUPS {
    tag "${sample_id}"
    publishDir "${params.outdir}/loops", mode: 'copy'
    cpus 4
    memory { 16.GB * task.attempt }

    input:
    tuple val(sample_id), path(hic)

    output:
    path("${sample_id}.hiccups_loops.bedpe"), emit: loops

    script:
    // The image has no CUDA runtime, so use HiCCUPS' CPU mode unless a GPU is requested.
    // CPU mode only searches near the diagonal (8 Mb by default).
    def cpu_flag = params.hiccups_gpu ? '' : '--cpu'
    def heap_gb  = Math.max(1, (task.memory.toGiga() * 0.85) as int)
    // Juicer's published settings per resolution: peak width (-p), window width (-i) and the
    // radius for merging nearby pixels (-d). HiCCUPS wants one value per resolution for each
    // of them (HiCCUPSConfiguration reads every list with resolutions.length), whatever its
    // usage text says about -d.
    def settings = [
        5000L:  [peak: 4, window: 7, radius: 20000],
        10000L: [peak: 2, window: 5, radius: 20000],
        25000L: [peak: 1, window: 3, radius: 50000],
    ]
    def loop_res = parseResolutions(params.hiccups_resolutions)
    def fdr      = loop_res.collect { _r -> '0.1' }.join(',')
    def peaks    = loop_res.collect { r -> settings[r].peak }.join(',')
    def windows  = loop_res.collect { r -> settings[r].window }.join(',')
    def radii    = loop_res.collect { r -> settings[r].radius }.join(',')
    """
    java -Xmx${heap_gb}g -jar /opt/juicer_tools.jar hiccups \\
        ${cpu_flag} \\
        --threads ${task.cpus} \\
        -k KR \\
        -r ${loop_res.join(',')} \\
        -f ${fdr} \\
        -p ${peaks} \\
        -i ${windows} \\
        -d ${radii} \\
        ${hic} \\
        hiccups_out/

    cp hiccups_out/merged_loops.bedpe ${sample_id}.hiccups_loops.bedpe
    """
}

process CONTACT_STATS {
    tag "${sample_id}"
    publishDir "${params.outdir}/qc", mode: 'copy'
    cpus 1
    memory { 4.GB * task.attempt }

    input:
    tuple val(sample_id), path(pairs)

    output:
    path("${sample_id}.contact_stats.txt"), emit: stats

    script:
    """
    pairtools stats ${pairs} -o ${sample_id}.contact_stats.txt
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
    multiqc --title "ENCODE Hi-C Pipeline" --filename multiqc_report --force .
    """
}

// ---- Workflow ----

workflow {
    // ---- Parameter validation ----
    if (!params.reads)       { error "Missing required parameter: --reads" }
    if (!params.bwa_index)   { error "Missing required parameter: --bwa_index" }
    if (!params.chrom_sizes) { error "Missing required parameter: --chrom_sizes" }

    // Reject empty lists, zero, negative and non-numeric values before anything does arithmetic
    // on them (the smallest resolution is a divisor below and cooler's bin size).
    [resolutions: params.resolutions, hiccups_resolutions: params.hiccups_resolutions].each { name, value ->
        // split(',', -1) keeps empty fields, so "5000," and "5000,,10000" are rejected too
        def tokens = value.toString().split(',', -1).collect { r -> r.trim() }
        if (tokens.any { r -> !(r ==~ /[1-9][0-9]*/) }) {
            error "--${name} must be a comma-separated list of positive integers (got '${value}')"
        }
    }

    // cooler bins once at the smallest resolution and coarsens from there, and HiCCUPS reads
    // its resolutions from the .hic file, so both lists have to be consistent.
    def resolutions = parseResolutions(params.resolutions)
    def off_grid    = resolutions.findAll { r -> r % resolutions.min() != 0 }
    if (off_grid) {
        error "--resolutions must all be multiples of the smallest one (${resolutions.min()}): ${off_grid.join(', ')} are not"
    }
    def loop_res    = parseResolutions(params.hiccups_resolutions)
    def unsupported = loop_res.findAll { r -> !(r in [5000L, 10000L, 25000L]) }
    if (unsupported) {
        error "--hiccups_resolutions accepts 5000, 10000 and 25000 only (got ${unsupported.join(', ')})"
    }
    def not_built   = loop_res.findAll { r -> !(r in resolutions) }
    if (not_built) {
        error "--hiccups_resolutions ${not_built.join(', ')} must also be listed in --resolutions"
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
    ch_reads       = channel.fromFilePairs(params.reads, checkIfExists: true)
    ch_bwa_index   = channel.fromPath("${params.bwa_index}*", checkIfExists: true).collect()
    ch_chrom_sizes = channel.fromPath(params.chrom_sizes, checkIfExists: true).collect()

    FASTQC_RAW(ch_reads)
    BWA_ALIGN(ch_reads, ch_bwa_index)
    PAIRTOOLS_PARSE_SORT(BWA_ALIGN.out.bam, ch_chrom_sizes)
    PAIRTOOLS_DEDUP(PAIRTOOLS_PARSE_SORT.out.pairs)
    PAIRTOOLS_SELECT(PAIRTOOLS_DEDUP.out.pairs)

    JUICER_HIC(PAIRTOOLS_SELECT.out.pairs, ch_chrom_sizes)
    COOLER_MCOOL(PAIRTOOLS_SELECT.out.pairs, ch_chrom_sizes)
    HICCUPS(JUICER_HIC.out.hic)
    CONTACT_STATS(PAIRTOOLS_SELECT.out.pairs)

    ch_multiqc = FASTQC_RAW.out.reports
        .mix(PAIRTOOLS_PARSE_SORT.out.stats)
        .mix(PAIRTOOLS_DEDUP.out.stats)
        .mix(PAIRTOOLS_SELECT.out.stats)
        .collect()

    MULTIQC(ch_multiqc)
}
