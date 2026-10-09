#!/usr/bin/env nextflow
/*
 * auto16s — one fixed 16S rRNA workflow, raw paired-end FASTQ to statistics and figures
 *
 *   FastQC
 *   Primers: detected automatically (library of common 16S primers) and removed per read
 *   DADA2: truncLen (auto) + maxEE -> learnErrors -> dada -> mergePairs -> chimeras
 *   SILVA taxonomy (assignTaxonomy)
 *   phyloseq: remove unassigned / Eukaryota / chloroplast / mitochondria; rarefy for alpha/beta
 *   composition (unrarefied) | alpha + beta diversity (rarefied) | MaAsLin3 genus level (unrarefied)
 *   MultiQC: QC + key figures and test results in one report
 */
nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW      } from './modules/local/fastqc'
include { FASTQC as FASTQC_FILTERED } from './modules/local/fastqc'
include { PRIMER_DETECT   } from './modules/local/primer_detect'
include { PRIMER_TRIM     } from './modules/local/primer_trim'
include { DADA2           } from './modules/local/dada2'
include { DADA2_TAXONOMY  } from './modules/local/dada2_taxonomy'
include { PHYLOSEQ_BUILD  } from './modules/local/phyloseq_build'
include { COMPOSITION     } from './modules/local/composition'
include { ALPHA_DIVERSITY } from './modules/local/alpha_diversity'
include { BETA_DIVERSITY  } from './modules/local/beta_diversity'
include { MAASLIN3        } from './modules/local/maaslin3'
include { MULTIQC         } from './modules/local/multiqc'

def helpMessage() {
    log.info """
    auto16s ${workflow.manifest.version}

    nextflow run main.nf -profile docker \\
        --input samplesheet.csv

    Required
      --input            Samplesheet (CSV: sample,fastq_1,fastq_2,${params.group_col}[,other columns]),
                         or a folder of FASTQ files: writes a samplesheet to fill in and stops

    Optional (sensible defaults; change only if needed)
      --group_col        Samplesheet column with the groups to compare   [${params.group_col}]
      --silva_db         DADA2 SILVA training set (local path or URL)     [SILVA 138.1, Zenodo]
      --outdir           Output folder                                    [${params.outdir}]
      --maaslin_reference  Reference group for MaAsLin3                   [alphabetically first]
      --fw_primer / --rv_primer   Primer sequences, only if yours are not detected automatically
      --trunc_len_f / --trunc_len_r   Override automatic truncLen
      --rarefy_depth     Override automatic rarefaction depth
    """.stripIndent()
}

// ---------------------------------------------------------------------------------------------
// Samplesheet helpers (written for Nextflow's strict syntax: plain functions, no top-level statements)
// ---------------------------------------------------------------------------------------------
def isFastq(name) {
    return name ==~ /(?i).*\.(fastq|fq)(\.gz)?$/
}

def csvQuote(v) {
    return '"' + (v ?: '').toString().replace('"', '""') + '"'
}

// Comma, semicolon (Turkish/European Excel) or tab, whichever the header line uses most
def detectSeparator(f) {
    def first = (f.readLines()[0] ?: '').replace('﻿', '')
    def counts = [',': first.count(','), ';': first.count(';'), '\t': first.count('\t')]
    def best = counts.max { e -> e.value }
    return best.value > 0 ? best.key : ','
}

def resolveFastq(sheet, p) {
    def f = (p.startsWith('/') || p.contains('://')) ? file(p) : file(sheet.parent.resolve(p).toString())
    if (!f.exists()) error "File not found: ${p} (relative paths are read from the samplesheet's folder: ${sheet.parent})"
    if (!isFastq(f.name)) error "Not a FASTQ file (.fastq, .fq, .fastq.gz or .fq.gz): ${p}"
    return f
}

def readSamplesheet(path) {
    def sheet = file(path, checkIfExists: true)
    def sep   = detectSeparator(sheet)
    def rows  = sheet.splitCsv(header: true, strip: true, sep: sep, quote: '"')
        .collect { r -> r.collectEntries { k, v -> [(k.replace('﻿', '').trim()): v?.trim()] } }   // Excel BOM
        .findAll { r -> r.values().any { v -> v } }                                                     // empty lines
    if (!rows) error "Samplesheet ${path} is empty."
    def absent = ['sample', 'fastq_1', 'fastq_2', params.group_col] - rows[0].keySet()
    if (absent) error "Samplesheet is missing column(s): ${absent.join(', ')}. The first line must be: sample,fastq_1,fastq_2,${params.group_col}"

    def seen = [] as Set
    def samples = rows.collect { row ->
        if (!(row.sample ==~ /[A-Za-z][A-Za-z0-9._-]*/))
            error "Sample name '${row.sample}': use English letters, digits, '.', '_' or '-', starting with a letter (no spaces)"
        if (!seen.add(row.sample)) error "Duplicate sample name: ${row.sample}"
        if (!row.fastq_1 || !row.fastq_2) error "Sample ${row.sample}: fastq_1 and fastq_2 are both required"
        if (!row[params.group_col]) error "Sample ${row.sample}: empty '${params.group_col}'"
        [[id: row.sample, group: row[params.group_col]], [resolveFastq(sheet, row.fastq_1), resolveFastq(sheet, row.fastq_2)]]
    }
    def groups = samples.collect { s -> s[0].group }.unique()
    if (groups.size() < 2) error "'${params.group_col}' must contain at least 2 groups (found: ${groups.join(', ')})"
    def sizes = groups.collect { g -> [g, samples.count { s -> s[0].group == g }] }
    sizes.findAll { gs -> gs[1] < 3 }.each { gs -> log.warn "Group '${gs[0]}' has only ${gs[1]} sample(s); statistical tests will have very little power." }
    log.info "Samplesheet: ${samples.size()} samples; groups: ${sizes.collect { gs -> "${gs[0]} (${gs[1]})" }.join(', ')}"

    // comma-separated copy with absolute paths for the R steps
    def cols = rows[0].keySet() as List
    def lines = [cols.collect { c -> csvQuote(c) }.join(',')]
    rows.eachWithIndex { r, i ->
        lines << cols.collect { c ->
            c == 'fastq_1' ? csvQuote(samples[i][1][0].toString()) : (c == 'fastq_2' ? csvQuote(samples[i][1][1].toString()) : csvQuote(r[c]))
        }.join(',')
    }
    def norm = file("${workflow.workDir}/samplesheet_normalized.csv")
    norm.parent.mkdirs()
    norm.text = lines.join('\n') + '\n'
    return [samples, norm]
}

// Sample name from a file-name prefix: Illumina suffix removed, accents and Turkish letters to ASCII
def cleanSampleName(p) {
    def n = p.replaceAll(/_S\d+(_L\d{3})?$/, '')
    n = java.text.Normalizer.normalize(n, java.text.Normalizer.Form.NFD).replaceAll(/\p{M}/, '')
    n = n.replace('ı', 'i').replace('İ', 'I').replaceAll(/[^A-Za-z0-9._-]/, '_')
    return (n ==~ /[A-Za-z].*/) ? n : 'S' + n
}

def uniqueName(n, used) {
    def u = used.contains(n) ? (2..100000).collect { k -> "${n}_${k}".toString() }.find { c -> !used.contains(c) } : n
    used.add(u)
    return u
}

// --input <folder>: pair R1/R2 files by name and write a samplesheet to fill in
def draftSamplesheet(dir) {
    def fq = dir.listFiles().findAll { f -> f.isFile() && isFastq(f.name) }.sort { f -> f.name }
    if (!fq) error "No FASTQ files (.fastq, .fq, .fastq.gz, .fq.gz) in ${dir}"
    if (fq.any { f -> f.name.contains('�') })
        error "Some file names contain characters that cannot be read with the current system language setting. Run with a UTF-8 locale (export LANG=C.UTF-8) or rename the files using English letters only."
    def pairs = [:]
    def unmatched = []
    fq.each { f ->
        def m = f.name =~ /^(.+?)[._-]R?([12])(_001)?\.(?i:fastq|fq)(\.gz)?$/
        if (m.matches()) {
            pairs.get(m.group(1), [:])[m.group(2)] = f
        } else {
            unmatched << f.name
        }
    }
    def complete = pairs.findAll { _k, v -> v['1'] && v['2'] }
    unmatched += pairs.findAll { _k, v -> !(v['1'] && v['2']) }.collect { _k, v -> v.values()*.name }.flatten()
    if (!complete) error "Could not pair any R1/R2 files in ${dir}. Expected names like S1_R1.fastq.gz / S1_R2.fastq.gz or S1_1.fq.gz / S1_2.fq.gz."

    def used = [] as Set
    def lines = ["sample,fastq_1,fastq_2,${params.group_col}".toString()]
    complete.each { k, v -> lines << "${uniqueName(cleanSampleName(k), used)},${v['1']},${v['2']},".toString() }
    def target = file("${workflow.launchDir}/samplesheet.csv")
    if (target.exists()) target = file("${workflow.launchDir}/samplesheet_draft.csv")
    target.text = lines.join('\n') + '\n'

    log.info "\n${complete.size()} samples found in ${dir}\nSamplesheet written to: ${target}\n\n" +
             "Next step: open it, fill in the '${params.group_col}' column for every sample (e.g. control / patient),\n" +
             "then run the pipeline again with --input ${target.name}\n"
    if (unmatched) log.warn "Files without a matching R1/R2 partner (not included): ${unmatched.join(', ')}"
}

// DADA2's filtered reads (<sample>_F_filt.fastq.gz / _R_) back to [meta, [R1, R2]]
def filteredPairs(ch) {
    return ch
        .flatten()
        .map { f -> [f.name.replaceAll(/_[FR]_filt\.fastq\.gz$/, ''), f] }
        .groupTuple(size: 2)
        .map { id, files -> [[id: id], files.sort { f -> f.name }] }      // _F_ sorts before _R_
}

def toManifest(ch) {
    return ch
        .toSortedList { a, b -> a[0].id <=> b[0].id }
        .map { list -> [list.collect { x -> x[0].id }, list.collect { x -> x[1][0] }, list.collect { x -> x[1][1] }] }
}

workflow {
    if (params.help) {
        helpMessage()
    } else if (params.input == null) {
        error "Missing required parameter: --input. Use --help."
    } else if (file(params.input).isDirectory()) {
        draftSamplesheet(file(params.input))
    } else {
        def sheet = readSamplesheet(params.input)
        log.info "Results will be written to ${params.outdir}/; the report will be ${params.outdir}/multiqc/multiqc_report.html"
        ch_reads = channel.fromList(sheet[0])
        ch_sheet = channel.value(sheet[1])
        ch_silva = channel.value(file(params.silva_db, checkIfExists: true))

        // QC of the raw reads
        FASTQC_RAW(ch_reads, 'raw')

        // Primers: detect once from all samples, then remove per sample
        PRIMER_DETECT(toManifest(ch_reads), file("${projectDir}/assets/primers_16s.tsv", checkIfExists: true))
        PRIMER_TRIM(ch_reads, PRIMER_DETECT.out.primers)

        // ASVs: all samples together (one error model for the run)
        DADA2(toManifest(PRIMER_TRIM.out.reads), PRIMER_TRIM.out.stats.collect())
        DADA2_TAXONOMY(DADA2.out.seqtab, ch_silva)

        // QC again after primer removal and DADA2 filtering: did the trimming work?
        FASTQC_FILTERED(filteredPairs(DADA2.out.filtered), 'filtered')

        // phyloseq + rarefaction
        PHYLOSEQ_BUILD(DADA2.out.seqtab, DADA2_TAXONOMY.out.taxa, ch_sheet)

        // statistics and figures
        COMPOSITION(PHYLOSEQ_BUILD.out.ps)
        ALPHA_DIVERSITY(PHYLOSEQ_BUILD.out.ps_rare)
        BETA_DIVERSITY(PHYLOSEQ_BUILD.out.ps_rare)
        MAASLIN3(PHYLOSEQ_BUILD.out.ps)

        // one report
        ch_mqc = FASTQC_RAW.out.zip.map { _meta, z -> z }.flatten()
            .mix(FASTQC_FILTERED.out.zip.map { _meta, z -> z }.flatten())
            .mix(PRIMER_DETECT.out.mqc.flatten(), DADA2.out.mqc.flatten(), PHYLOSEQ_BUILD.out.mqc.flatten(),
                 COMPOSITION.out.mqc.flatten(), ALPHA_DIVERSITY.out.mqc.flatten(), BETA_DIVERSITY.out.mqc.flatten(),
                 MAASLIN3.out.mqc.flatten())
            .collect()
        MULTIQC(ch_mqc, file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true))
    }

}
