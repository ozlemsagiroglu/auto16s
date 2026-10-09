process PHYLOSEQ_BUILD {
    label 'process_medium'
    publishDir "${params.outdir}/phyloseq", mode: params.publish_dir_mode, saveAs: { fn -> fn.endsWith('_mqc.png') || fn.endsWith('_mqc.tsv') ? null : fn }

    conda "${projectDir}/envs/r_analysis.yml"
    container "${ params.r_container }"

    input:
    path seqtab
    path taxa
    path samplesheet

    output:
    path '*.{png,pdf}'      , emit: plots, optional: true
    path '*.{tsv,csv,txt,rds,fasta}', emit: tables, optional: true
    path '*_mqc.{tsv,png}'  , emit: mqc, optional: true
    path 'phyloseq_object.rds'          , emit: ps
    path 'phyloseq_object_rarefied.rds' , emit: ps_rare

    script:
    """
    phyloseq_build.R --seqtab ${seqtab} --taxa ${taxa} --samplesheet ${samplesheet} \\
        --group_col '${params.group_col}' --rarefy_depth ${params.rarefy_depth ?: 'NA'} --seed ${params.seed}
    """

    stub:
    """
    touch phyloseq_object.rds phyloseq_object_rarefied.rds asv_table.tsv sample_depths_mqc.tsv rarefaction_curves_mqc.png
    """
}
