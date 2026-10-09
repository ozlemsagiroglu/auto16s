process BETA_DIVERSITY {
    label 'process_low'
    publishDir "${params.outdir}/beta_diversity", mode: params.publish_dir_mode, saveAs: { fn -> fn.endsWith('_mqc.png') || fn.endsWith('_mqc.tsv') ? null : fn }

    conda "${projectDir}/envs/r_analysis.yml"
    container "${ params.r_container }"

    input:
    path ps_rare

    output:
    path '*.{png,pdf}'      , emit: plots, optional: true
    path '*.{tsv,csv,txt,rds,fasta}', emit: tables, optional: true
    path '*_mqc.{tsv,png}'  , emit: mqc, optional: true


    script:
    """
    beta_diversity.R --ps ${ps_rare} --group_col '${params.group_col}' --permutations ${params.permutations} --seed ${params.seed}
    """

    stub:
    """
    touch pcoa_bray.png beta_diversity_tests.tsv pcoa_bray_mqc.png
    """
}
