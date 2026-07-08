# Lab Notebook

## 2026-07-08

Created `REPO_NOTES.md` with a codebase overview for tealeaf.

Key components documented:

- `tealeaf-map` builds annotation-derived isoform-to-intron and isoform-to-exon
  references.
- `tealeaf-cluster` maps Salmon transcript abundance to intron and exon counts,
  then runs shared cluster refinement.
- `tealeaf-sc` aggregates barcodes into pseudobulk samples, estimates transcript
  abundance from alevin-fry equivalence classes using EM, then reuses the shared
  intron counting and clustering path.
- `tealeaf/shared_functions.py` contains the common cluster construction,
  filtering, refinement, and PSI-ratio calculation logic.
- `tealeaf-ggsashimi` plots tealeaf intron/exon count files with an adapted
  ggsashimi workflow.

Implementation notes:

- Intron representation in the shared clustering code is
  `[start, end, total_count, name, exon_set]`.
- Bulk and single-cell pipelines both converge on `{prefix}count_intron`,
  `{prefix}count_exon`, `{prefix}refined_cluster`, and `{prefix}ratio_count`.
- Potential issues observed: `pyproject.toml` misspells the ggsashimi entry
  point, `cluster_def` is parsed but not passed to `process_clusters()`, and
  `extra_scripts/ccp_gen.py` uses `networkx` without declaring it as a package
  dependency.

## 2026-07-08 Microglia Split-Pool Smoke Test

Input data:

- `/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/salmon_spliceu`
- Transcript-level alevin-fry output:
  `out_permit_known/quant_t2t_dedup/alevin`
- Salmon spliceu reference FASTA:
  `/gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu/spliceu.fa`

Reference build:

- The spliceu `clean_gtf.gtf` is gene-only and cannot be used by
  `tealeaf-map`.
- Built a tealeaf reference from
  `/gpfs/commons/home/daknowles/knowles_lab/index/kallisto/mus_musculus/with_precursor/gencode.vM32.basic.annotation.gtf`.
- Output prefix:
  `/scratch/daknowles/tealeaf_microglia_test/ref/vM32_`.
- `tealeaf-map` finished in 17m 11s and produced isoform-intron/exon maps,
  sparse matrices, connectivity, and source-map files.

Smoke-test run:

- Streamed the first 1,000 rows from the transcript-level
  `quants_mat.mtx` and aggregated them into 10 pseudobulk samples of 100 cells.
- Wrote preprocessed tealeaf-sc inputs under
  `/scratch/daknowles/tealeaf_microglia_test/run/smoke_*`.
- Ran `tealeaf/sc/tealeaf_sc.py --preprocessed` with `--ref_prefix vM32_`,
  `--introncutoff 5`, and `--minclucounts 10`.
- The run completed in 4m 23s.

Outputs:

- `smoke_count_intron`: 180,934 lines.
- `smoke_count_exon`: 213,661 lines.
- `smoke_refined_cluster`: 29,574 lines.
- `smoke_ratio_count`: 29,574 lines.

Implementation note:

- `tealeaf_sc()` had a preprocessed-mode bug: `out_prefix` was only assigned
  inside the non-preprocessed branch. Moved `out_prefix = options.outprefix` to
  the top of the function so `--preprocessed` can use existing pseudo matrices.
