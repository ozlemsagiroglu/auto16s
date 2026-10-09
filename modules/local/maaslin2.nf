process MAASLIN2 {
    label 'process_medium'
    publishDir "${params.outdir}/maaslin2", mode: params.publish_dir_mode, saveAs: { fn -> fn.endsWith('_mqc.png') || fn.endsWith('_mqc.tsv') ? null : fn }

    conda "${projectDir}/envs/r_analysis.yml"
    container "${ params.r_container }"

    input:
    path ps

    output:
    path '*.{png,pdf}'      , emit: plots, optional: true
    path '*.{tsv,csv,txt,rds,fasta}', emit: tables, optional: true
    path '*_mqc.{tsv,png}'  , emit: mqc, optional: true
    path 'maaslin2_output'                , emit: raw

    script:
    """
    maaslin2.R --ps ${ps} --group_col '${params.group_col}' --reference ${params.maaslin_reference ?: 'NA'} \\
        --min_prevalence ${params.maaslin_min_prevalence} --qval ${params.maaslin_qval} --cpus ${task.cpus}
    """

    stub:
    """
    touch maaslin2_volcano.png maaslin2_significant.csv; mkdir maaslin2_output
    """
}
