process ALPHA_DIVERSITY {
    tag "${cmp.id}"
    label 'process_low'
    publishDir path: { cmp.publish ? "${params.outdir}/${cmp.publish}/alpha_diversity" : "${params.outdir}/alpha_diversity" }, mode: params.publish_dir_mode, saveAs: { fn -> fn.endsWith('_mqc.png') || fn.endsWith('_mqc.tsv') ? null : fn }

    conda "${moduleDir}/../../envs/r_analysis.yml"
    container "${ params.r_container }"

    input:
    tuple val(cmp), path(ps_rare)

    output:
    path '*.{png,pdf}'      , emit: plots, optional: true
    path '*.{tsv,csv,txt,rds,fasta}', emit: tables, optional: true
    tuple val(cmp), path('*_mqc.{tsv,png}'), emit: mqc, optional: true


    script:
    def args = "--group_col '${cmp.group}' --subset_col '${cmp.subset_col}' --subset_value '${cmp.subset_value}' " +
               "--covariates '${cmp.covariates}' --interaction ${cmp.interaction} --prefix '${cmp.prefix}' --label '${cmp.label}'"
    """
    alpha_diversity.R --ps ${ps_rare} ${args}
    """

    stub:
    """
    touch ${cmp.prefix}alpha_diversity.png ${cmp.prefix}alpha_diversity.tsv ${cmp.prefix}alpha_diversity_mqc.png
    """
}
