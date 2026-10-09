#!/usr/bin/env nextflow
/*
 * auto16s — one fixed 16S rRNA workflow, raw paired-end FASTQ to statistics and figures
 *
 *   FastQC
 *   Primers: detected automatically (library of common 16S primers) and removed per read
 *   DADA2: truncLen (auto) + maxEE -> learnErrors -> dada -> mergePairs -> chimeras
 *   SILVA taxonomy (assignTaxonomy)
 *   phyloseq: remove unassigned / Eukaryota / chloroplast / mitochondria; rarefy for alpha/beta
 *   composition (unrarefied) | alpha + beta diversity (rarefied) | MaAsLin2 genus level (unrarefied)
 *   MultiQC: QC + key figures and test results in one report
 */
nextflow.enable.dsl = 2

include { FASTQC          } from './modules/local/fastqc'
include { PRIMER_DETECT   } from './modules/local/primer_detect'
include { PRIMER_TRIM     } from './modules/local/primer_trim'
include { DADA2           } from './modules/local/dada2'
include { DADA2_TAXONOMY  } from './modules/local/dada2_taxonomy'
include { PHYLOSEQ_BUILD  } from './modules/local/phyloseq_build'
include { COMPOSITION     } from './modules/local/composition'
include { ALPHA_DIVERSITY } from './modules/local/alpha_diversity'
include { BETA_DIVERSITY  } from './modules/local/beta_diversity'
include { MAASLIN2        } from './modules/local/maaslin2'
include { MULTIQC         } from './modules/local/multiqc'

def helpMessage() {
    log.info """
    auto16s ${workflow.manifest.version}

    nextflow run main.nf -profile docker \\
        --input samplesheet.csv

    Required
      --input            CSV: sample,fastq_1,fastq_2,${params.group_col}[,other metadata columns]

    Optional (sensible defaults; change only if needed)
      --group_col        Samplesheet column with the groups to compare   [${params.group_col}]
      --silva_db         DADA2 SILVA training set (local path or URL)     [SILVA 138.1, Zenodo]
      --outdir           Output folder                                    [${params.outdir}]
      --maaslin_reference  Reference group for MaAsLin2                   [alphabetically first]
      --fw_primer / --rv_primer   Primer sequences, only if yours are not detected automatically
      --trunc_len_f / --trunc_len_r   Override automatic truncLen
      --rarefy_depth     Override automatic rarefaction depth
    """.stripIndent()
}

def readSamplesheet(path) {
    def sheet = file(path, checkIfExists: true)
    def rows  = sheet.splitCsv(header: true, strip: true)
    if (!rows) error "Samplesheet ${path} is empty."
    def absent = ['sample', 'fastq_1', 'fastq_2', params.group_col] - rows[0].keySet()
    if (absent) error "Samplesheet is missing column(s): ${absent.join(', ')}"

    def resolve = { String p ->
        (p.startsWith('/') || p.contains('://')) ? file(p, checkIfExists: true)
                                                 : file(sheet.parent.resolve(p).toString(), checkIfExists: true)
    }
    def seen = [] as Set
    def out = rows.collect { row ->
        if (!(row.sample ==~ /[A-Za-z][A-Za-z0-9._-]*/))
            error "Sample name '${row.sample}': use letters, digits, '.', '_' or '-', starting with a letter"
        if (!seen.add(row.sample)) error "Duplicate sample name: ${row.sample}"
        if (!row.fastq_1 || !row.fastq_2) error "Sample ${row.sample}: fastq_1 and fastq_2 are both required"
        if (!row[params.group_col]) error "Sample ${row.sample}: empty '${params.group_col}'"
        [[id: row.sample, group: row[params.group_col]], [resolve(row.fastq_1), resolve(row.fastq_2)]]
    }
    def groups = out.collect { it[0].group }.unique()
    if (groups.size() < 2) error "'${params.group_col}' must contain at least 2 groups (found: ${groups.join(', ')})"
    groups.each { g ->
        def n = out.count { it[0].group == g }
        if (n < 3) log.warn "Group '${g}' has only ${n} sample(s); statistical tests will have very little power."
    }
    return out
}

workflow {
    if (params.help) { helpMessage(); return }
    def missing = ['input'].findAll { params[it] == null }
    if (missing) error "Missing required parameter(s): ${missing.collect { '--' + it }.join(', ')}. Use --help."

    ch_reads = Channel.fromList(readSamplesheet(params.input))
    ch_sheet = Channel.value(file(params.input))
    ch_silva = Channel.value(file(params.silva_db, checkIfExists: true))

    // QC
    FASTQC(ch_reads)

    // Primers: detect once from all samples, then remove per sample
    def toManifest = { ch -> ch.toSortedList { a, b -> a[0].id <=> b[0].id }
        .map { list -> [list.collect { it[0].id }, list.collect { it[1][0] }, list.collect { it[1][1] }] } }
    PRIMER_DETECT(toManifest(ch_reads), file("${projectDir}/assets/primers_16s.tsv", checkIfExists: true))
    PRIMER_TRIM(ch_reads, PRIMER_DETECT.out.primers)

    // ASVs: all samples together (one error model for the run)
    DADA2(toManifest(PRIMER_TRIM.out.reads), PRIMER_TRIM.out.stats.collect())
    DADA2_TAXONOMY(DADA2.out.seqtab, ch_silva)

    // phyloseq + rarefaction
    PHYLOSEQ_BUILD(DADA2.out.seqtab, DADA2_TAXONOMY.out.taxa, ch_sheet)

    // statistics and figures
    COMPOSITION(PHYLOSEQ_BUILD.out.ps)
    ALPHA_DIVERSITY(PHYLOSEQ_BUILD.out.ps_rare)
    BETA_DIVERSITY(PHYLOSEQ_BUILD.out.ps_rare)
    MAASLIN2(PHYLOSEQ_BUILD.out.ps)

    // one report
    ch_mqc = FASTQC.out.zip.map { meta, z -> z }.flatten()
        .mix(PRIMER_DETECT.out.mqc.flatten(), DADA2.out.mqc.flatten(), PHYLOSEQ_BUILD.out.mqc.flatten(), COMPOSITION.out.mqc.flatten(),
             ALPHA_DIVERSITY.out.mqc.flatten(), BETA_DIVERSITY.out.mqc.flatten(), MAASLIN2.out.mqc.flatten())
        .collect()
    MULTIQC(ch_mqc, file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true))
}

workflow.onComplete {
    log.info(workflow.success
        ? "\nDone. Open ${params.outdir}/multiqc/multiqc_report.html\n"
        : "\nFailed. See .nextflow.log and the work directory of the failed task.\n")
}
