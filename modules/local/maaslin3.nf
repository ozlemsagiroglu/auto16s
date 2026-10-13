process MAASLIN3 {
    tag "${cmp.id}"
    label 'process_medium'
    publishDir path: { cmp.publish ? "${params.outdir}/${cmp.publish}/maaslin3" : "${params.outdir}/maaslin3" }, mode: params.publish_dir_mode, saveAs: { fn -> fn.endsWith('_mqc.png') || fn.endsWith('_mqc.tsv') ? null : fn }

    conda "${moduleDir}/../../envs/r_analysis.yml"
    container "${ params.r_container }"

    input:
    tuple val(cmp), path(ps)

    output:
    path '*.{png,pdf}'      , emit: plots, optional: true
    path '*.{tsv,csv,txt,rds,fasta}', emit: tables, optional: true
    tuple val(cmp), path('*_mqc.{tsv,png}'), emit: mqc, optional: true
    path '*maaslin3*_output', emit: raw

    script:
    def ref = cmp.reference ?: 'NA'
    def args = "--group_col '${cmp.group}' --subset_col '${cmp.subset_col}' --subset_value '${cmp.subset_value}' " +
               "--covariates '${cmp.covariates}' --interaction ${cmp.interaction} --prefix '${cmp.prefix}' --label '${cmp.label}'"
    """
    maaslin3.R --ps ${ps} ${args} --reference '${ref}' \\
        --min_prevalence ${params.maaslin_min_prevalence} --qval ${params.maaslin_qval}
    """

    stub:
    """
    touch ${cmp.prefix}maaslin3_volcano.png ${cmp.prefix}maaslin3_significant.csv; mkdir ${cmp.prefix}maaslin3_output
    """
}
