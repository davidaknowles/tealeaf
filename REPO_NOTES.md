# tealeaf Repository Notes

## Purpose

`tealeaf` is a Python package and command-line toolkit for Transcript
Expression Augmented LEAFcutter. It combines transcript-level abundance
estimates with annotation-derived transcript structure to produce
LeafCutter-like intron clusters and PSI-style intron usage ratios.

The main use case is alternative splicing analysis when direct junction
coverage is limited, including low-coverage bulk, pseudobulk, and
single-cell RNA-seq settings.

## Package Entry Points

Console commands are exposed through `tealeaf/__main__.py`:

- `tealeaf-map`: builds isoform-to-intron and isoform-to-exon maps from a GTF.
- `tealeaf-cluster`: runs the bulk or pseudobulk Salmon `quant.sf` pipeline.
- `tealeaf-sc`: runs the single-cell alevin-fry pipeline.
- `tealeaf-ggsashimi`: plots tealeaf intron/exon counts in a sashimi-style view.

Packaging is defined in both `pyproject.toml` and `setup.py`. Note that
`pyproject.toml` appears to misspell the ggsashimi entry point as
`tealeaf_ggshashimi`; `setup.py` uses the correct `tealeaf_ggsashimi`.

## Main Data Flow

### 1. Reference Map Generation

Implemented in `tealeaf/map_gen/tealeaf_map_gen.py`.

`compute_transcript_intron_map()` reads a GTF with `pyranges`, filters unwanted
annotations, and walks through genes and transcripts to infer the introns and
exons supported by each transcript.

Main outputs:

- `isoform_intron_map.tsv`: transcript to supported introns and exons.
- `intron_exon_connectivity.tsv`: intron to neighboring exons and strand.
- `intron_source_map.tsv`: intron to transcript source type for Gencode input.
- `isoform_intron_matrix.npz`: sparse isoform by intron matrix for single-cell.
- `isoform_exon_matrix.npz`: sparse isoform by exon matrix for single-cell.
- `isoform_rows.txt`, `intron_cols.txt`, `exon_cols.txt`: sparse matrix labels.

Optional virtual first and last introns can be added to capture alternative
first exon and alternative last exon usage, but this is marked as experimental
in the code and README.

### 2. Bulk and Pseudobulk Counting

Implemented in `tealeaf/clustering/tealeaf_clustering.py`.

The bulk path reads a list of Salmon `quant.sf` files. If `--use_TPM` is not
set, it first creates a `normalized_count` column using one of three modes:

- `junction`: simulate junction support from reads, effective transcript length,
  read length, overhang, paired-end status, and a sizing factor.
- `global`: globally scale TPM by the relationship between total reads and total
  TPM.
- `local`: scale TPM within each gene using transcript-to-gene mapping from the
  GTF.

`count_introns()` then projects transcript abundance onto all introns and exons
supported by each transcript in the reference map.

Main outputs:

- `{prefix}count_intron`
- `{prefix}count_exon`
- `{prefix}dropped_trancripts`

### 3. Single-Cell Pipeline

Implemented in `tealeaf/sc/tealeaf_sc.py` and `tealeaf/sc/sc_utils.py`.

The single-cell path starts from alevin-fry equivalence-class output. It can
generate pseudobulk samples from barcode-to-cell-type assignments in two ways:

- `metacells`: shuffle barcodes within a type and group fixed-size chunks.
- `bootstrapping`: repeatedly sample fixed-size barcode sets within a type.

`pseudo_eq_conversion()` loads barcode by equivalence-class counts, aggregates
barcodes into pseudobulk samples, filters low-count equivalence classes and
transcripts, then runs EM to estimate transcript abundance. The EM
implementation is `sc_utils.EM()`, which estimates transcript proportions from
equivalence-class counts using transcript effective-length weights.

After transcript quantification, `sc_intron_count()` multiplies the pseudobulk
transcript matrix by the reference isoform-intron and isoform-exon sparse
matrices. This produces the same `{prefix}count_intron` and
`{prefix}count_exon` files consumed by the shared clustering code.

### 4. Shared Clustering and PSI Calculation

Implemented in `tealeaf/shared_functions.py`.

This file contains the common LeafCutter-like clustering logic used by both the
bulk and single-cell pipelines.

Important functions:

- `build_init_cluster()`: sorts introns by chromosome, strand, start, and end,
  then groups overlapping introns into initial clusters.
- `process_clusters()`: filters and refines clusters.
- `filter_introns()`: removes low-support introns by absolute count and fraction
  of cluster support.
- `cluster_intervals()`: reclusters overlapping intervals after filtering.
- `refine_links()`: refines clusters by shared splice sites and optionally
  shared neighboring exon sites.
- `compute_ratio()`: converts refined intron counts into per-cluster intron
  usage ratios.

The internal intron representation used in this file is:

```python
[start, end, total_count, name, exon_set]
```

Cluster refinement modes are intended to be:

- `1`: overlap only.
- `2`: overlap plus shared intron splice site.
- `3`: overlap plus shared intron splice site plus shared exon splice site.

The final clustering outputs are:

- `{prefix}refined_cluster`
- `{prefix}ratio_count`

`ratio_count` divides each intron count by the total count of its cluster for
each sample, yielding PSI-style values compatible with downstream LeafCutter
tools.

## Plotting

Implemented in `tealeaf/ggsashimi/tealeaf_ggsashimi.py`.

This is an adapted ggsashimi-style plotter. Instead of reading BAM files, it
reads tealeaf intron and exon count files. It reconstructs exon coverage and
junction arcs over a requested genomic interval, reads GTF annotation, generates
an R script using `ggplot2`, `gridExtra`, and `data.table`, and runs R to write
the plot.

Important inputs:

- `--intron`: tealeaf intron count file.
- `--exon`: tealeaf exon count file.
- `--strand_info`: `intron_exon_connectivity.tsv` from `tealeaf-map`.
- `--coordinates`: genomic interval, `chr:start-end`.
- `--gtf`: annotation GTF.

## Extra Scripts

`extra_scripts/ccp_gen.py` generates cluster edge representations and compares
cluster structure across samples using graph operations. It uses NetworkX, which
is not listed in the package dependencies.

`extra_scripts/unify_definition.py` appears to be unfinished or untested. It
calls `result_df.save_csv`, which should likely be `result_df.to_csv`.

## Important File Formats

Introns are usually named:

```text
chr:start-end
```

Rows in refined clusters are named:

```text
chr:start:end:clu_N_strand
```

Bulk and single-cell pipelines converge on the same count files:

```text
Name Chr Start End [Gene] sample_1 sample_2 ...
```

`compute_ratio()` expects refined cluster row names to contain the cluster ID in
the fourth colon-delimited field.

## Notes and Potential Issues

- `pyproject.toml` has a likely typo in the `tealeaf-ggsashimi` script target.
- `cluster_def` is parsed in the bulk and single-cell CLIs, but the current calls
  to `process_clusters()` do not pass `mode=options.cluster_def`, so mode `3`
  is always used.
- `sample_data/` is currently untracked in this working tree.
- Some CLI help text and function names contain typos, for example
  `dropped_trancripts` and several help-string spelling mistakes.
- `extra_scripts/ccp_gen.py` imports `networkx`, but `networkx` is not listed in
  `pyproject.toml` dependencies.
