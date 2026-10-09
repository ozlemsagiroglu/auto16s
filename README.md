# auto16s

[![CI](https://github.com/ozlemsagiroglu/auto16s/actions/workflows/ci.yml/badge.svg)](https://github.com/ozlemsagiroglu/auto16s/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://apptainer.org/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**auto16s** analyses paired-end 16S rRNA amplicon data, from raw FASTQ files to statistics and publication-ready figures, in one fixed workflow.

It is built for the common case: one Illumina run and a comparison between groups. There are no alternative routes to choose from. The steps and their settings follow the standard DADA2-based workflow used in the microbiome literature. The pipeline determines from the data what usually has to be looked up by hand: which primers are in the reads, where to truncate, and how deep to rarefy. Every decision is shown in the report.

The only required input is a samplesheet.

![auto16s overview](docs/images/auto16s_metro_map.png)

## Contents

- [Pipeline summary](#pipeline-summary)
- [Quick start](#quick-start)
- [Input](#input)
- [Parameters and defaults](#parameters-and-defaults)
- [What the pipeline decides from the data](#what-the-pipeline-decides-from-the-data)
- [Quality filtering](#quality-filtering)
- [Outputs](#outputs)
- [Statistics](#statistics)
- [Figures](#figures)
- [Reading the report](#reading-the-report)
- [Scope and limitations](#scope-and-limitations)
- [Testing](#testing)
- [Citations](#citations)
- [License](#license)

## Pipeline summary

1. **Read QC** ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. **Primer detection**: the read starts are matched against a library of common 16S primers ([`assets/primers_16s.tsv`](assets/primers_16s.tsv)).
3. **Primer removal**, per read and by sequence. This handles variable-length spacers before the primer and pairs in reverse orientation. Pairs without both primers are dropped.
4. **ASV inference** ([DADA2](https://benjjneb.github.io/dada2/)):
   - quality filtering with automatic `truncLen` and `maxEE`;
   - error models for R1 and R2;
   - denoising;
   - `mergePairs`;
   - consensus chimera removal.
5. **Taxonomy**: DADA2 `assignTaxonomy` against [SILVA 138.1](https://doi.org/10.5281/zenodo.4587955).
6. **phyloseq object** ([phyloseq](https://joey711.github.io/phyloseq/)). Unassigned, eukaryotic, chloroplast and mitochondrial ASVs are removed, then rarefaction at an automatically chosen depth.
7. **Composition**: phylum, family and genus bar plots per sample and per group, plus a genus heatmap (all reads).
8. **Alpha diversity**: Observed, Shannon and Simpson, with Wilcoxon or Kruskal-Wallis tests (rarefied).
9. **Beta diversity**: Bray-Curtis PCoA, PERMANOVA, betadisper and pairwise PERMANOVA ([vegan](https://github.com/vegandevs/vegan); rarefied).
10. **Differential abundance** at genus level ([MaAsLin2](https://huttenhower.sph.harvard.edu/maaslin/); all reads).
11. **Report** ([MultiQC](https://multiqc.info/)): QC, warnings, all key figures and test results in one HTML file.

## Quick start

1. Install [Nextflow](https://www.nextflow.io/docs/latest/install.html) (≥ 24.04, Java 17+) and **one** of the software stacks below. On Windows, use [WSL](https://learn.microsoft.com/en-us/windows/wsl/install).

2. Test the installation (8 small samples, a few minutes). Nextflow downloads the pipeline from GitHub:

   ```bash
   nextflow run ozlemsagiroglu/auto16s -profile test,docker      # or test,singularity / test,conda
   ```

3. Run your data:

   ```bash
   nextflow run ozlemsagiroglu/auto16s -profile docker --input samplesheet.csv --outdir results
   ```

4. Open `results/multiqc/multiqc_report.html`. Read the **"DADA2 summary and warnings"** table first.

Add `-resume` to rerun after a change; finished steps are reused. Add `-r v0.3.0` to run a fixed release.

### Software stacks: Docker is not required

All tools come either from containers or from Conda. Choose the profile that matches what is installed:

| Profile | Needs | Notes |
|---|---|---|
| `docker` | Docker | Simplest on a laptop or workstation |
| `singularity` / `apptainer` | Singularity or Apptainer | **No Docker needed.** Usual choice on HPC clusters; the same images are pulled and converted automatically |
| `conda` | Conda / Mamba / Miniforge | **No containers at all.** Environments are created on the first run (this takes a while once) |

The tools of the R analysis steps (phyloseq, vegan, MaAsLin2, ggplot2) are in the image `ghcr.io/ozlemsagiroglu/auto16s-r`, built from [`containers/r-analysis`](containers/r-analysis) and published automatically. Nothing has to be built by hand. For offline clusters, the image can be built without Docker:

```bash
cd containers/r-analysis && apptainer build auto16s-r.sif auto16s-r.def
nextflow run ozlemsagiroglu/auto16s -profile singularity --r_container $PWD/auto16s-r.sif --input samplesheet.csv
```

**Without any container or Conda**, the pipeline also runs on a machine where FastQC, MultiQC and R are installed. It needs R ≥ 4.3 with dada2, phyloseq, vegan, MaAsLin2 and ggplot2; the R packages can be installed with `Rscript containers/r-analysis/install_r_packages.R` (plus `BiocManager::install("dada2")`). Then run without `-profile`.

## Input

One CSV file holds the reads **and** the sample metadata:

```csv
sample,fastq_1,fastq_2,group,age
Patient01,raw/Patient01_R1.fastq.gz,raw/Patient01_R2.fastq.gz,control,34
Patient02,raw/Patient02_R1.fastq.gz,raw/Patient02_R2.fastq.gz,case,41
```

| Column | Description |
|---|---|
| `sample` | Unique name. Starts with a letter; letters, digits, `.`, `_` and `-` only |
| `fastq_1`, `fastq_2` | R1 and R2 files (`.fastq.gz`). Absolute paths, or relative to the CSV's folder |
| `group` | The groups to compare (at least 2). Another column can be used with `--group_col` |
| any other column | Kept as sample metadata in the phyloseq object |

The samplesheet is checked before anything runs. The pipeline stops on duplicate names, missing files, empty groups or fewer than 2 groups, and warns about groups with fewer than 3 samples.

## Parameters and defaults

**Required:** `--input`.

**Optional**, only if needed:

| Parameter | Default | Description |
|---|---|---|
| `--outdir` | `results` | Output folder |
| `--group_col` | `group` | Samplesheet column with the groups |
| `--maaslin_reference` | alphabetically first group | Reference group for MaAsLin2 coefficients |
| `--silva_db` | SILVA 138.1 training set, downloaded from Zenodo | Local copy, to avoid downloading in every new run |
| `--fw_primer`, `--rv_primer` | detected | Primer sequences (5'→3', IUPAC), only if your primers are not in the library |
| `--trunc_len_f`, `--trunc_len_r` | automatic | Fixed truncation lengths |
| `--rarefy_depth` | automatic | Fixed rarefaction depth |
| `--max_cpus`, `--max_memory`, `--max_time` | `8`, `32.GB`, `24.h` | Upper limits per task |

**Settings applied in every run.** These are literature-standard values, stored in [`nextflow.config`](nextflow.config):

| Step | Setting | Value |
|---|---|---|
| Primer removal | primer start | within the first 12 bases of the read |
| | mismatches allowed | 1 (primer < 20 nt) or 2 (≥ 20 nt), IUPAC codes honoured |
| DADA2 `filterAndTrim` | `truncLen` | automatic (see below) |
| | `maxEE` | 2 (R1 and R2) |
| | `truncQ`, `maxN`, `rm.phix` | 2, 0, `TRUE` |
| DADA2 `mergePairs` | `minOverlap`, `trimOverhang` | 12, `TRUE` |
| Chimera removal | `removeBimeraDenovo` | `consensus` |
| Taxonomy | `assignTaxonomy` `minBoot` | 50 |
| Rarefaction | depth, seed | automatic (see below), 711 |
| PERMANOVA / betadisper | permutations | 999 |
| MaAsLin2 | model, normalisation, transform | LM, TSS, LOG |
| | `min_prevalence`, significance | 0.1, q < 0.05 (Benjamini-Hochberg) |
| Bar plots | taxa shown | top 12 by mean relative abundance + "Other" |

## What the pipeline decides from the data

### Primers

The first reads of every sample are matched against [`assets/primers_16s.tsv`](assets/primers_16s.tsv). The list has 20 widely used 16S primers, among them 27F, 341F, 515F, 515F-Y, 799F, 805R, 806R, 806RB, 926R and 1193R. The forward/reverse pair found in most read pairs is selected and removed from every read by sequence, which also handles:

- **spacers**: variable-length bases before the primer ("heterogeneity spacers", phased primers);
- **mixed orientation**: pairs whose R1 starts with the reverse primer, as in some ligation-based libraries, are swapped.

Pairs that do not carry both primers are removed; they are usually off-target products. The detected primers, the share of read pairs that carry them and a per-sample summary are in `primers/` and in the report.

If no known primer pair is found, the pipeline assumes the primers were already removed by the sequencing provider and uses the reads as they are. The report then shows the most common read start, so an unlisted primer can be recognised. Give such a primer with `--fw_primer`/`--rv_primer`, or add one line to the TSV.

### Truncation length (`truncLen`)

The truncation length is computed from the reads themselves. FastQC is only for viewing. On the first reads of every sample, after primer removal, the pipeline:

1. computes the median quality per position and truncates where it drops below Q25;
2. estimates the amplicon length by aligning R1 to the reverse complement of R2;
3. checks that R1 and R2 still overlap by at least 20 bp after truncation; if not, it truncates less.

Step 3 matters because a too-short truncation is the most common reason why read pairs silently fail to merge in DADA2. The quality profile, the chosen positions and the expected overlap are shown in `dada2/dada2_quality_truncation.png`.

### Rarefaction depth

Rarefaction is used for alpha and beta diversity only. The depth is the smallest depth among samples with at least 1,000 reads and at least 10 % of the median depth. One failed sample therefore cannot force all others down to a few hundred reads. Samples below the depth are left out of alpha and beta diversity and listed in the report. Composition and MaAsLin2 use all samples and all reads.

## Quality filtering

Quality control is done by DADA2 rather than a separate trimmer, so it may help to know how it compares with Trimmomatic.

| | Trimmomatic | DADA2 (this pipeline) |
|---|---|---|
| Low-quality 3' end | cut per read (sliding window) → reads of different lengths | cut at the same position in every read (`truncLen`) |
| Read with many errors | kept if no window falls below the threshold | **removed** if its expected number of errors exceeds 2 (`maxEE`) |
| Single low-quality bases | cut the read there | kept; handled by the error model learned from the run |

The expected number of errors of a read is the sum of the error probabilities of its bases (Q30 = 0.001, Q20 = 0.01, Q10 = 0.1). It measures the error load of the whole read, which a sliding window does not. Reads of equal length are a requirement of DADA2's denoising, not a cosmetic choice. The few single errors that remain are corrected by the denoising itself; correcting them is what an ASV method does.

On the test data (sample S01, 2,500 read pairs):

| | Kept (72 %) | Removed (28 %) |
|---|---|---|
| Median expected errors R1 / R2 | 0.24 / 0.22 | 3.84 / 1.52 |
| Mean quality R1 / R2 | Q36.4 / Q35.9 | Q28.7 / Q31.5 |
| Bases below Q20 (median, R1 / R2) | 1.3 % / 1.4 % | 19.2 % / 11.9 % |

None of the kept pairs exceeds 2 expected errors. Every removed pair would also have been trimmed by a `SLIDINGWINDOW:4:20` check. About half of the kept pairs contain a short dip below Q20 that Trimmomatic would have cut. DADA2 keeps them because their overall error load is low (median 0.8 expected errors).

## Outputs

All results are in `--outdir`. Besides the figures, every number behind a figure is delivered as a table, and the full data as phyloseq objects, so you can make your own figures or analyses.

| Folder | Content |
|---|---|
| `multiqc/` | `multiqc_report.html`: everything in one report |
| `phyloseq/` | `phyloseq_object.rds` (all reads), `phyloseq_object_rarefied.rds`, ASV tables (raw and rarefied), taxonomy, ASV sequences (FASTA), sample metadata, depths, rarefaction curves |
| `composition/` | Counts and relative abundances per phylum, family and genus; bar plots; genus heatmap |
| `alpha_diversity/` | Values per sample, tests, figure |
| `beta_diversity/` | Bray-Curtis distance matrix, PCoA coordinates, PERMANOVA, betadisper, interpretation, figures |
| `maaslin2/` | All results and significant genera (CSV), volcano, coefficient plot, boxplots, raw MaAsLin2 output |
| `primers/` | Detected primers, detection plot, per-sample primer statistics |
| `dada2/` | Summary and warnings, read tracking, truncation, error models, ASV lengths, sequence table |
| `taxonomy/` | SILVA assignment per ASV sequence |
| `fastqc/` | FastQC report per file |
| `pipeline_info/` | Nextflow execution report, timeline, trace |

Each file is described in [docs/output.md](docs/output.md). Figures are saved as PNG (300 dpi) and PDF.

**Continuing in R**

```r
library(phyloseq)
ps      <- readRDS("results/phyloseq/phyloseq_object.rds")            # all reads
ps_rare <- readRDS("results/phyloseq/phyloseq_object_rarefied.rds")   # rarefied

sample_data(ps)                       # your samplesheet columns
plot_bar(tax_glom(ps, "Phylum"), fill = "Phylum")

genus <- read.delim("results/composition/genus_relative_abundance.tsv", row.names = 1, check.names = FALSE)
```

## Statistics

| Analysis | Data | Method |
|---|---|---|
| Alpha diversity | rarefied | 2 groups: Wilcoxon rank-sum. More than 2: Kruskal-Wallis and pairwise Wilcoxon, Benjamini-Hochberg |
| Beta diversity | rarefied | Bray-Curtis; PERMANOVA (`adonis2`) and `betadisper`. More than 2 groups: pairwise PERMANOVA, Benjamini-Hochberg |
| Differential abundance | all reads, genus level | MaAsLin2 (LM, TSS, LOG), Benjamini-Hochberg q-values |

PERMANOVA cannot tell a shift in composition from a difference in within-group spread. The pipeline therefore reads it together with betadisper and writes the interpretation into the report. Rarefying for alpha and beta diversity, but not for differential abundance, follows current recommendations (Schloss 2024).

## Figures

All figures are made with ggplot2. Colours are fixed, so a group or a taxon looks the same in every figure:

| Element | Colours |
|---|---|
| Groups | Okabe-Ito palette (colour-blind safe), in alphabetical order of group names |
| Taxa in bar plots | 20 distinct colours by abundance rank; "Other" in grey |
| Heatmap | viridis, log10 relative abundance |
| MaAsLin2 | orange = higher, blue = lower than the reference group |

## Reading the report

The report starts with **DADA2 summary and warnings**. Possible warnings:

| Warning | Most likely cause |
|---|---|
| No known 16S primer pair found | Primers were already removed (fine), or a primer that is not in the library. Check the reported read start |
| Only x % of read pairs carry the primers | Off-target products, or a primer variant not in the library |
| < 80 % of reads survived chimera removal | Primers still in the reads |
| Median merge rate < 70 % | Reads too short to overlap after truncation |

**Sanity check: ASV length.** `dada2/dada2_asv_length.png` should show one main peak at the length of your region without primers (e.g. V4 ≈ 253 bp, V3-V4 ≈ 400–430 bp). Small peaks far from it are usually off-target products.

## Scope and limitations

The workflow is deliberately fixed and assumes:

- **One Illumina paired-end run.** DADA2's error model is run-specific; samples from different runs should not be combined in one analysis.
- **Primers from the built-in list**, or given with `--fw_primer`/`--rv_primer`.
- **No negative controls.** If extraction or PCR blanks were sequenced, contaminant removal (e.g. [decontam](https://github.com/benjjneb/decontam)) is standard and is not included.
- **No phylogenetic tree,** so UniFrac and Faith's PD are not computed.
- **A comparison between groups.** Covariates and paired or repeated-measures designs are not modelled; use the phyloseq objects in R for those.
- **One differential abundance method.** Differential abundance methods can give different results on the same data (Nearing et al. 2022). Report MaAsLin2 results as such.
- **Small groups.** Permutation tests cannot produce small p-values with very few samples. With 3 vs 2 samples there are only 10 distinct permutations, so p ≥ 0.1.

For other designs, data types or analyses, [nf-core/ampliseq](https://nf-co.re/ampliseq) offers many more options.

## Testing

`-profile test` runs real reads: the V4 example data shipped with DADA2 (2×250 bp), with 515F-Y/806RB primers added and their degenerate bases resolved per read. It has 8 samples in 2 groups with a built-in difference in composition. Two other scenarios can be generated with the same script:

```bash
Rscript tests/make_test_data.R tests/data_spacer spacer     # 0-7 nt spacers, half of the pairs reverse-oriented
Rscript tests/make_test_data.R tests/data_noprimer none     # primers already removed
nextflow run main.nf -profile docker --input tests/data_spacer/samplesheet.csv \
    --silva_db tests/data/example_train_set.fa.gz
```

In all three, the primers are found (or correctly reported as absent), the amplicon length is estimated at 253 bp (V4), and no reads are lost to chimeras.

## Citations

If you use auto16s, please cite the tools it relies on:

- **Nextflow**: Di Tommaso P. et al. (2017). Nextflow enables reproducible computational workflows. *Nat Biotechnol* 35:316–319.
- **FastQC**: Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data.
- **DADA2**: Callahan B.J. et al. (2016). DADA2: high-resolution sample inference from Illumina amplicon data. *Nat Methods* 13:581–583.
- **SILVA**: Quast C. et al. (2013). The SILVA ribosomal RNA gene database project. *Nucleic Acids Res* 41:D590–D596.
- **SILVA for DADA2**: McLaren M.R., Callahan B.J. (2021). SILVA 138.1 prokaryotic SSU taxonomic training data formatted for DADA2. Zenodo. doi:10.5281/zenodo.4587955
- **phyloseq**: McMurdie P.J., Holmes S. (2013). phyloseq. *PLoS ONE* 8:e61217.
- **vegan**: Oksanen J. et al. vegan: Community Ecology Package.
- **betadisper**: Anderson M.J. (2006). Distance-based tests for homogeneity of multivariate dispersions. *Biometrics* 62:245–253.
- **MaAsLin2**: Mallick H. et al. (2021). Multivariable association discovery in population-scale meta-omics studies. *PLoS Comput Biol* 17:e1009442.
- **MultiQC**: Ewels P. et al. (2016). MultiQC. *Bioinformatics* 32:3047–3048.
- **ggplot2**: Wickham H. (2016). *ggplot2: Elegant Graphics for Data Analysis*. Springer.

Background for design choices:

- Klindworth A. et al. (2013). Evaluation of general 16S ribosomal RNA gene PCR primers. *Nucleic Acids Res* 41:e1.
- Schloss P.D. (2024). Rarefaction is currently the best approach to control for uneven sequencing effort in amplicon sequence analyses. *mSphere*.
- Nearing J.T. et al. (2022). Microbiome differential abundance methods produce different results across 38 datasets. *Nat Commun* 13:342.

## License

auto16s is released under the [MIT License](LICENSE). The tools it runs keep their own licenses.
