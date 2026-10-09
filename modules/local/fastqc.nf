process FASTQC {
    tag "${meta.id}:${stage}"
    label 'process_low'
    publishDir path: { "${params.outdir}/fastqc/${stage}" }, mode: params.publish_dir_mode

    conda 'bioconda::fastqc=0.13.0'
    container "${ workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/816bda6dd014e810753f642611232443ddeabb2b47f9681ad43f9044816f00ff/data'
        : 'community.wave.seqera.io/library/fastqc:0.13.0--b342db7b694c9af4' }"

    input:
    tuple val(meta), path(reads, stageAs: 'input/*')
    val stage    // 'raw' or 'filtered' (after primer removal and DADA2 filterAndTrim)

    output:
    tuple val(meta), path('*.html'), emit: html
    tuple val(meta), path('*.zip') , emit: zip

    script:
    // Uniform names: "<sample>_R1" (raw) / "<sample>_filtered_R1"; keep .gz only if the input is compressed
    def ext = reads[0].name.endsWith('.gz') ? 'fastq.gz' : 'fastq'
    def prefix = stage == 'raw' ? meta.id : "${meta.id}_${stage}"
    """
    ln -s ${reads[0]} ${prefix}_R1.${ext}
    ln -s ${reads[1]} ${prefix}_R2.${ext}
    fastqc --quiet --threads ${task.cpus} ${prefix}_R1.${ext} ${prefix}_R2.${ext}
    """

    stub:
    def prefix = stage == 'raw' ? meta.id : "${meta.id}_${stage}"
    """
    touch ${prefix}_R1_fastqc.html ${prefix}_R2_fastqc.html ${prefix}_R1_fastqc.zip ${prefix}_R2_fastqc.zip
    """
}
