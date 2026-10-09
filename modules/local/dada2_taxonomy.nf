process DADA2_TAXONOMY {
    label 'process_high'
    publishDir "${params.outdir}/taxonomy", mode: params.publish_dir_mode

    conda "bioconda::bioconductor-dada2=1.38.0 conda-forge::r-base=4.5.2 conda-forge::r-digest=0.6.39 conda-forge::tbb=2022.3.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/81153df5d53322e6d91b2c4c9bc4da50774fb1d101ead002fe75bb75fc6f036c/data'
        : 'community.wave.seqera.io/library/bioconductor-dada2_r-base_r-digest_tbb:38acac09bac46f36' }"

    input:
    path seqtab
    path silva_db

    output:
    path 'taxonomy.rds', emit: taxa
    path 'taxonomy.tsv', emit: tsv

    script:
    """
    dada2_taxonomy.R --seqtab ${seqtab} --db ${silva_db} --min_boot ${params.min_boot} \\
        --cpus ${task.cpus} --seed ${params.seed}
    """

    stub:
    """
    touch taxonomy.rds taxonomy.tsv
    """
}
