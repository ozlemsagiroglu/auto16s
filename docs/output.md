# auto16s: output

All paths are relative to `--outdir` (default `results/`). Figures are saved as PNG (300 dpi) and PDF with the same name; only the PNG is listed. Tables are tab-separated (`.tsv`) unless they are `.csv`.

## Start here: `multiqc/`

| File | Content |
|---|---|
| `multiqc_report.html` | One self-contained report: DADA2 summary and warnings, primer detection, FastQC, read tracking, filtering, depths, and all key figures and tests |
| `multiqc_report_data/` | The report's data in machine-readable form |

## `primers/`

| File | Content |
|---|---|
| `primers.tsv` | Detected (or given) primers: names, sequences, 16S region, share of read pairs with both primers, share in reverse orientation, warnings. `status` is `detected`, `user` or `none` |
| `primer_detection.png` | Share of R1 and R2 reads starting with each candidate primer; the selected pair is highlighted |
| `stats/<sample>.primer_stats.tsv` | Per sample: read pairs in, pairs with both primers, % reverse-oriented (swapped), longest spacer seen, and `primer_left_R1_pct`/`primer_left_R2_pct`: % of reads that still contain the primer after removal (should be ≈ 0) |

## `dada2/`

| File | Content |
|---|---|
| `dada2_truncation.tsv` | Chosen `truncLen` and why, estimated amplicon length, expected R1/R2 overlap, share of read pairs expected to pass filtering, merge rate, chimera rate, **QC warnings** |
| `dada2_read_tracking.tsv` | Read pairs per sample after each step: input, primers found, filtered, denoised, merged, non-chimeric; merge rate and % retained |
| `dada2_quality_truncation.png` | Median quality (and IQR) per position for R1 and R2 with the chosen truncation |
| `dada2_read_tracking_plot.png` | % of input read pairs retained per step, one line per sample |
| `dada2_errors_R1.png`, `dada2_errors_R2.png` | Learned error models (points: observed; black line: fitted). The line should follow the points |
| `dada2_asv_length.png` (also in the report) | Reads per ASV length; one main peak at the region length is expected |
| `seqtab_nochim.rds` | DADA2 sequence table (samples × ASV sequences), before taxonomic filtering |
| `dada2.log` | Full log of the DADA2 step |

## `taxonomy/`

| File | Content |
|---|---|
| `taxonomy.tsv` / `.rds` | SILVA assignment (Kingdom to Genus) for every ASV sequence; `NA` = not assigned at that rank (bootstrap < 50) |

## `phyloseq/`

| File | Content |
|---|---|
| `phyloseq_object.rds` | phyloseq object: all reads after taxonomic filtering, sample metadata from the samplesheet. Used for composition and MaAsLin 3 |
| `phyloseq_object_rarefied.rds` | The same, rarefied. Used for alpha and beta diversity |
| `asv_table.tsv` | ASV × sample read counts (all reads) |
| `asv_table_rarefied.tsv` | ASV × sample read counts (rarefied) |
| `asv_taxonomy.tsv` | ASV id, taxonomy and sequence |
| `asv_sequences.fasta` | ASV sequences; ids `ASV_1`, `ASV_2`, … ordered by total abundance |
| `sample_metadata.tsv` | Sample metadata as used in the analysis |
| `taxonomic_filtering.tsv` | ASVs and reads removed as unassigned, Eukaryota, chloroplast or mitochondria |
| `sample_depths.tsv` | Reads per sample after filtering; whether the sample is used for alpha/beta diversity |
| `sample_depths_plot.png` | Reads per sample with the rarefaction depth |
| `rarefaction_curves.tsv` / `.png` | Observed ASVs vs. reads sampled, per sample, with the rarefaction depth |

## `composition/` (all reads)

| File | Content |
|---|---|
| `<rank>_counts.tsv` | Read counts per phylum / family / genus and sample. Unassigned ranks are named after the nearest assigned rank, e.g. `Lachnospiraceae (unclassified)` |
| `<rank>_relative_abundance.tsv` | Relative abundance (0-1) per phylum / family / genus and sample |
| `<rank>_barplot_samples.png` | Stacked bars per sample, faceted by group: top 12 taxa + Other |
| `<rank>_barplot_groups.png` | Mean composition per group |
| `genus_heatmap.png` | Top 25 genera, log10 relative abundance; genera clustered, samples clustered within group |

## `alpha_diversity/` (rarefied)

| File | Content |
|---|---|
| `alpha_diversity.tsv` | Observed, Shannon and Simpson (1 − D) per sample |
| `alpha_diversity_tests.tsv` | Wilcoxon (2 groups) or Kruskal-Wallis + pairwise Wilcoxon with BH-adjusted p (> 2 groups) |
| `alpha_diversity.png` | Box plots with all samples; test p-values in the panel titles |

## `beta_diversity/` (rarefied)

| File | Content |
|---|---|
| `bray_curtis_distance.tsv` | Sample × sample Bray-Curtis distance matrix |
| `pcoa_bray_coords.tsv` | PCoA coordinates (axes 1 and 2) per sample |
| `beta_diversity_tests.tsv` | PERMANOVA (R², F, p) and betadisper (F, p) |
| `permanova_pairwise.tsv` | Pairwise PERMANOVA with BH-adjusted p (only with > 2 groups) |
| `betadisper_distances.tsv` | Distance of every sample to its group centroid |
| `beta_diversity_interpretation.txt` | PERMANOVA and betadisper read together, in one sentence |
| `pcoa_bray.png` | PCoA with 95 % ellipses (≥ 4 samples per group) and group centroids |
| `betadisper_bray.png` | Distances to centroid per group |

## `maaslin3/` (all reads, genus level)

| File | Content |
|---|---|
| `maaslin3_all_results.csv` | Every tested genus, once per model: `model` (`abundance` or `prevalence`), comparison, coefficient, standard error, p, q (`qval_individual`), `qval_joint` (either model), number of samples / non-zero samples, `error` (model did not fit) |
| `maaslin3_significant.csv` | Rows with q < 0.05 and no fitting error |
| `maaslin3_volcano.png` | Coefficient vs. −log10(q), one panel per model; significant genera labelled |
| `maaslin3_coefficients.png` | Coefficients ± SE of the significant genera, per model (positive = higher than the reference group) |
| `maaslin3_abundance_boxplots.png` | Abundance hits: relative abundance in the samples where the genus was found (up to 12 genera) |
| `maaslin3_prevalence_bars.png` | Prevalence hits: share of samples per group in which the genus was found (up to 12 genera) |
| `maaslin3_output/` | Unmodified MaAsLin 3 output |

Genera without a SILVA genus assignment are named after the nearest assigned rank, e.g. `Lachnospiraceae (unclassified)`; the column `resolved` is `FALSE` for them. Do not report them as named genera.

## `fastqc/` and `pipeline_info/`

| File | Content |
|---|---|
| `fastqc/raw/<sample>_R1_fastqc.html`, `_R2` | FastQC report per raw read file |
| `fastqc/filtered/<sample>_filtered_R1_fastqc.html`, `_R2` | FastQC report of the reads after DADA2 quality filtering and truncation (the reads that go into denoising) |
| `pipeline_info/report.html`, `timeline.html`, `trace.tsv`, `dag.html` | Nextflow execution report: resources and run time per task |

Every R step also writes `<step>_sessionInfo.txt` with the exact R and package versions used.
