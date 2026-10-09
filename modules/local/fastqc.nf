process FASTQC {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/fastqc", mode: params.publish_dir_mode

    conda 'bioconda::fastqc=0.13.0'
    container "${ workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/816bda6dd014e810753f642611232443ddeabb2b47f9681ad43f9044816f00ff/data'
        : 'community.wave.seqera.io/library/fastqc:0.13.0--b342db7b694c9af4' }"

    input:
    tuple val(meta), path(reads, stageAs: 'input/*')

    output:
    tuple val(meta), path('*.html'), emit: html
    tuple val(meta), path('*.zip') , emit: zip

    script:
    // Uniform names -> MultiQC shows "<sample>_R1" / "<sample>_R2"
    """
    ln -s ${reads[0]} ${meta.id}_R1.fastq.gz
    ln -s ${reads[1]} ${meta.id}_R2.fastq.gz
    fastqc --quiet --threads ${task.cpus} ${meta.id}_R1.fastq.gz ${meta.id}_R2.fastq.gz
    """

    stub:
    """
    touch ${meta.id}_R1_fastqc.html ${meta.id}_R2_fastqc.html ${meta.id}_R1_fastqc.zip ${meta.id}_R2_fastqc.zip
    """
}
