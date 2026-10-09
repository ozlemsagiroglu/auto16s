process COMPOSITION {
    label 'process_low'
    publishDir "${params.outdir}/composition", mode: params.publish_dir_mode, saveAs: { fn -> fn.endsWith('_mqc.png') || fn.endsWith('_mqc.tsv') ? null : fn }

    conda "${moduleDir}/../../envs/r_analysis.yml"
    container "${ params.r_container }"

    input:
    path ps

    output:
    path '*.{png,pdf}'      , emit: plots, optional: true
    path '*.{tsv,csv,txt,rds,fasta}', emit: tables, optional: true
    path '*_mqc.{tsv,png}'  , emit: mqc, optional: true


    script:
    """
    composition.R --ps ${ps} --group_col '${params.group_col}' --top_n ${params.top_n}
    """

    stub:
    """
    touch genus_barplot_samples.png genus_heatmap_mqc.png
    """
}
