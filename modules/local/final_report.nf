process FINAL_REPORT {
    label 'process_single'
    publishDir "${params.outdir}/report", mode: params.publish_dir_mode

    conda "${moduleDir}/../../envs/r_analysis.yml"
    container "${ params.r_container }"

    input:
    path files, stageAs: 'inputs/*'
    path samplesheet

    output:
    path 'auto16s_report.html', emit: report

    script:
    """
    final_report.R --inputs inputs --samplesheet ${samplesheet} --out auto16s_report.html \\
        --group_col '${params.group_col}' --qval ${params.maaslin_qval} --min_prevalence ${params.maaslin_min_prevalence} \\
        --max_ee ${params.max_ee} --min_overlap ${params.min_overlap} --permutations ${params.permutations} --seed ${params.seed} \\
        --silva_db '${file(params.silva_db).name}' --pipeline_version '${workflow.manifest.version}' \\
        --nextflow_version '${workflow.nextflow.version}' --run_name '${workflow.runName}'
    """

    stub:
    """
    touch auto16s_report.html
    """
}
