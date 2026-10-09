process DADA2 {
    label 'process_high'
    publishDir "${params.outdir}/dada2", mode: params.publish_dir_mode, saveAs: { fn -> (fn.contains('_mqc.') || fn.startsWith('filt/')) ? null : fn }

    conda "bioconda::bioconductor-dada2=1.38.0 conda-forge::r-base=4.5.2 conda-forge::r-digest=0.6.39 conda-forge::tbb=2022.3.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/81153df5d53322e6d91b2c4c9bc4da50774fb1d101ead002fe75bb75fc6f036c/data'
        : 'community.wave.seqera.io/library/bioconductor-dada2_r-base_r-digest_tbb:38acac09bac46f36' }"

    input:
    tuple val(ids), path(r1, stageAs: 'reads/r1_?/*'), path(r2, stageAs: 'reads/r2_?/*')
    path primer_stats

    output:
    path 'seqtab_nochim.rds'          , emit: seqtab
    path 'dada2_read_tracking.tsv'    , emit: tracking
    path 'dada2_truncation.tsv'       , emit: truncation
    path '*_mqc.{tsv,png}'            , emit: mqc
    path '*.{png,pdf}'                , emit: plots
    path 'dada2.log'                  , emit: log
    path '*_sessionInfo.txt'          , emit: session
    path 'filt/*_filt.fastq.gz'       , emit: filtered, optional: true   // for FastQC after filtering

    script:
    def r1s = r1 instanceof List ? r1 : [r1]
    def r2s = r2 instanceof List ? r2 : [r2]
    def tf = params.trunc_len_f ?: 'NA'
    def tr = params.trunc_len_r ?: 'NA'
    // filterAndTrim and removeBimeraDenovo fork one R process per thread (~1-2 GB each):
    // use at most one thread per 2 GB of the task's memory
    def threads = Math.max(1, Math.min(task.cpus as int, (task.memory.toMega() / 2048) as int))
    """
    ids=(${ids.join(' ')})
    r1=(${r1s.join(' ')})
    r2=(${r2s.join(' ')})
    {
        printf 'sample\tr1\tr2\n'
        for i in "\${!ids[@]}"; do printf '%s\t%s\t%s\n' "\${ids[\$i]}" "\${r1[\$i]}" "\${r2[\$i]}"; done
    } > manifest.tsv

    dada2.R \\
        --manifest manifest.tsv \\
        --trunc_len_f ${tf} \\
        --trunc_len_r ${tr} \\
        --trunc_qmin ${params.trunc_qmin} \\
        --max_ee ${params.max_ee} \\
        --min_overlap ${params.min_overlap} \\
        --cpus ${threads} \\
        --seed ${params.seed} \\
        2>&1 | tee dada2.log
    """

    stub:
    """
    touch seqtab_nochim.rds dada2_read_tracking.tsv dada2_truncation.tsv dada2.log dada2_sessionInfo.txt
    touch dada2_read_tracking_mqc.tsv dada2_truncation_mqc.tsv dada2_quality_truncation_mqc.png
    touch dada2_quality_truncation.png dada2_quality_truncation.pdf
        mkdir filt; echo "" | gzip > filt/S1_F_filt.fastq.gz; echo "" | gzip > filt/S1_R_filt.fastq.gz
    """
}
