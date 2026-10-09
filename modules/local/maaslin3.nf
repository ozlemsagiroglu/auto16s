process MAASLIN3 {
    label 'process_medium'
    publishDir "${params.outdir}/maaslin3", mode: params.publish_dir_mode, saveAs: { fn -> fn.endsWith('_mqc.png') || fn.endsWith('_mqc.tsv') ? null : fn }

    conda "${moduleDir}/../../envs/r_analysis.yml"
    container "${ params.r_container }"

    input:
    path ps

    output:
    path '*.{png,pdf}'      , emit: plots, optional: true
    path '*.{tsv,csv,txt,rds,fasta}', emit: tables, optional: true
    path '*_mqc.{tsv,png}'  , emit: mqc, optional: true
    path 'maaslin3_output'  , emit: raw

    script:
    """
    maaslin3.R --ps ${ps} --group_col '${params.group_col}' --reference ${params.maaslin_reference ?: 'NA'} \\
        --min_prevalence ${params.maaslin_min_prevalence} --qval ${params.maaslin_qval}
    """

    stub:
    """
    touch maaslin3_volcano.png maaslin3_significant.csv; mkdir maaslin3_output
    """
}
