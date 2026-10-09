process MULTIQC {
    label 'process_single'
    publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode

    conda 'bioconda::multiqc=1.35'
    container "${ workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8e346f4f6080eadf1253505e6ff09ef004454fc18e8d672006fd7b222cc412e/data'
        : 'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc' }"

    input:
    path files, stageAs: 'inputs/*'
    path config

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_report_data', emit: data

    script:
    """
    multiqc --force --config ${config} --filename multiqc_report.html inputs
    """

    stub:
    """
    touch multiqc_report.html; mkdir multiqc_report_data
    """
}
