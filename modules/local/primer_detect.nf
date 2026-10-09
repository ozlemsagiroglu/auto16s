process PRIMER_DETECT {
    label 'process_low'
    publishDir "${params.outdir}/primers", mode: params.publish_dir_mode, saveAs: { fn -> fn.contains('_mqc.') ? null : fn }

    conda "bioconda::bioconductor-dada2=1.38.0 conda-forge::r-base=4.5.2 conda-forge::r-digest=0.6.39 conda-forge::tbb=2022.3.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/81153df5d53322e6d91b2c4c9bc4da50774fb1d101ead002fe75bb75fc6f036c/data'
        : 'community.wave.seqera.io/library/bioconductor-dada2_r-base_r-digest_tbb:38acac09bac46f36' }"

    input:
    tuple val(ids), path(r1, stageAs: 'reads/r1_?/*'), path(r2, stageAs: 'reads/r2_?/*')
    path primer_db

    output:
    path 'primers.tsv'           , emit: primers
    path '*_mqc.{tsv,png}'       , emit: mqc
    path 'primer_detection.{png,pdf}', emit: plots

    script:
    def r1s = r1 instanceof List ? r1 : [r1]
    def r2s = r2 instanceof List ? r2 : [r2]
    """
    ids=(${ids.join(' ')})
    r1=(${r1s.join(' ')})
    r2=(${r2s.join(' ')})
    {
        printf 'sample\\tr1\\tr2\\n'
        for i in "\${!ids[@]}"; do printf '%s\\t%s\\t%s\\n' "\${ids[\$i]}" "\${r1[\$i]}" "\${r2[\$i]}"; done
    } > manifest.tsv

    primer_detect.R --manifest manifest.tsv --primer_db ${primer_db} \\
        --fw_primer ${params.fw_primer ?: 'NA'} --rv_primer ${params.rv_primer ?: 'NA'}
    """

    stub:
    """
    printf 'key\\tvalue\\nstatus\\tdetected\\nfw_primer\\tGTGYCAGCMGCCGCGGTAA\\nrv_primer\\tGGACTACNVGGGTWTCTAAT\\n' > primers.tsv
    touch primer_detection_mqc.tsv primer_detection_mqc.png primer_detection.png primer_detection.pdf
    """
}
