nextflow.enable.dsl=2

genome_dir = file(params.hisat2_index_dir)
annotation = file(params.annotation)
samplesheet = file(params.samplesheet)
edger_script = file('src/gene_expression.R')

raw_reads = Channel
    .fromPath(samplesheet)
    .ifEmpty { error "Samplesheet not found at ${samplesheet}" }
    .splitCsv(header: true)
    .map { row ->
        def sampleId = row.sample_id?.toString()
        def fastq = row.fastq ? file(row.fastq) : null
        if (!sampleId) {
            error "Samplesheet row missing sample_id: ${row}"
        }
        if (!fastq?.exists()) {
            error "FASTQ file ${row.fastq} for sample ${sampleId} does not exist."
        }
        tuple(sampleId, fastq)
    }

genome_dir_ch = Channel.value(genome_dir)
annotation_ch = Channel.value(annotation)
samplesheet_ch = Channel.value(samplesheet)
edger_script_ch = Channel.value(edger_script)

process FastQC {
    tag "$sample_id"
    publishDir "${params.results_dir}/qc", mode: 'copy'
    container 'staphb/fastqc:0.11.9'
    cpus 2

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.html"), emit: html
    tuple val(sample_id), path("*_fastqc.zip"), emit: zip

    script:
    """
    fastqc ${reads} --outdir .
    """
}

process FastQCTrimmed {
    tag "$sample_id"
    publishDir "${params.results_dir}/qc/trimmed", mode: 'copy'
    container 'staphb/fastqc:0.11.9'
    cpus 2

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.html"), emit: html
    tuple val(sample_id), path("*_fastqc.zip"), emit: zip

    script:
    """
    fastqc ${reads} --outdir .
    """
}

process Fastp {
    tag "$sample_id"
    publishDir "${params.results_dir}/trimmed", mode: 'copy'
    container 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'
    cpus 4

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), emit: trimmed
    tuple val(sample_id), path("${sample_id}_fastp.html"), emit: html
    tuple val(sample_id), path("${sample_id}_fastp.json"), emit: json

    script:
    """
    fastp \\
        -i ${reads} \\
        -o ${sample_id}_trimmed.fastq.gz \\
        --html ${sample_id}_fastp.html \\
        --json ${sample_id}_fastp.json \\
        --thread ${task.cpus}
    """
}

process HISAT2_Align {
    tag "$sample_id"
    publishDir "${params.results_dir}/alignment", mode: 'copy'
    container 'quay.io/biocontainers/hisat2:2.2.1--h87f3376_3'
    cpus 8

    input:
    tuple val(sample_id), path(reads)
    path genome_dir

    output:
    tuple val(sample_id), path("${sample_id}.sam"), emit: sam
    tuple val(sample_id), path("${sample_id}_alignment_summary.txt"), emit: summary

    script:
    """
    hisat2 \\
        -p ${task.cpus} \\
        -x ${genome_dir}/${params.hisat2_index_prefix} \\
        -U ${reads} \\
        -S ${sample_id}.sam \\
        --summary-file ${sample_id}_alignment_summary.txt
    """
}

process SortBAM {
    tag "$sample_id"
    publishDir "${params.results_dir}/alignment", mode: 'copy'
    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
    cpus 4

    input:
    tuple val(sample_id), path(sam)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: bam
    path "${sample_id}_sorted.bam.bai", emit: bai

    script:
    """
    samtools view -bS ${sam} | samtools sort -o ${sample_id}_sorted.bam
    samtools index ${sample_id}_sorted.bam
    """
}

process FeatureCounts {
    tag "featureCounts"
    publishDir "${params.results_dir}/counts", mode: 'copy'
    container 'quay.io/biocontainers/subread:2.0.6--hed695b0_1'
    cpus 4

    input:
    path annotation
    path bam_files

    output:
    path "gene_counts.txt", emit: counts
    path "gene_counts.txt.summary", emit: summary

    script:
    """
    featureCounts \\
        -T ${task.cpus} \\
        -a ${annotation} \\
        -o gene_counts.txt \\
        ${bam_files}
    """
}

process EdgeRAnalysis {
    tag "edgeR"
    publishDir "${params.results_dir}/edger", mode: 'copy'
    container params.analysis_container
    cpus 2

    input:
    path counts
    path metadata
    path script_file

    output:
    path "edger_results.csv"
    path "deg_summary.csv"
    path "volcano_plot.png"
    path "pca_plot.png"
    path "top50_heatmap.png"
    path "significant_genes.txt"
    path "session_info.txt"

    script:
    """
    Rscript ${script_file} \\
        --counts ${counts} \\
        --metadata ${metadata} \\
        --outdir . \\
        --lfc 1 \\
        --pvalue 0.05
    """
}

workflow {
    raw_qc = FastQC(raw_reads)
    fastp_results = Fastp(raw_reads)
    trimmed_reads = fastp_results.trimmed
    trimmed_qc = FastQCTrimmed(trimmed_reads)
    aligned = HISAT2_Align(trimmed_reads, genome_dir_ch)
    sorted = SortBAM(aligned.sam)

    bam_files = sorted.bam.map { it[1] }.collect()
    fc = FeatureCounts(annotation_ch, bam_files)

    EdgeRAnalysis(fc.counts, samplesheet_ch, edger_script_ch)
}
