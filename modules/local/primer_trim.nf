process PRIMER_TRIM {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/primers/stats", mode: params.publish_dir_mode, pattern: '*.primer_stats.tsv'

    conda "bioconda::bioconductor-dada2=1.38.0 conda-forge::r-base=4.5.2 conda-forge::r-digest=0.6.39 conda-forge::tbb=2022.3.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/81153df5d53322e6d91b2c4c9bc4da50774fb1d101ead002fe75bb75fc6f036c/data'
        : 'community.wave.seqera.io/library/bioconductor-dada2_r-base_r-digest_tbb:38acac09bac46f36' }"

    input:
    tuple val(meta), path(reads, stageAs: 'input/*')
    path primers

    output:
    tuple val(meta), path("${meta.id}_R{1,2}.trim.fastq.gz"), emit: reads
    path "${meta.id}.primer_stats.tsv"                        , emit: stats

    script:
    """
    primer_trim.R --sample ${meta.id} --r1 ${reads[0]} --r2 ${reads[1]} --primers ${primers}
    """

    stub:
    """
    echo "" | gzip > ${meta.id}_R1.trim.fastq.gz; echo "" | gzip > ${meta.id}_R2.trim.fastq.gz
    printf 'sample\\tinput\\twith_primers\\n${meta.id}\\t1\\t1\\n' > ${meta.id}.primer_stats.tsv
    """
}
