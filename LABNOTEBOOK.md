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

## 2026-07-08 Microglia Cluster-Name Pseudobulk Run

The `calcutta` notebook `analyze_microglialess.ipynb` identified the relevant
combined AnnData object:

- `/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/salmon_per_sublib/t2g_permit_known_combined.h5ad`

The AnnData `.obs` table contains `CB_polydT`, `CB_ranhex`, `class`, and
`cluster_name`. Both barcode columns overlap the transcript-level alevin row
barcodes in:

- `/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/salmon_spliceu/out_permit_known/quant_t2t_dedup/alevin/quants_mat_rows.txt`

Generated barcode-to-cluster mapping:

- `/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/barcodes_to_cluster_name.csv`
- 88,920 barcode/group rows.
- 20 cluster-name groups.
- Barcodes with conflicting labels were dropped.

Moved the tealeaf reference from local scratch to a GPFS-visible run directory:

- `/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/ref/vM32_*`

Added helper scripts:

- `extra_scripts/aggregate_pseudobulk_sparse.py`: aggregates a cell-by-transcript
  sparse matrix into tealeaf-sc preprocessed pseudobulk transcript matrices.
- `extra_scripts/run_microglia_cluster_tealeaf.sbatch`: Slurm recipe for this
  microglia-less cluster-name run.

Run details:

- Slurm job: `18728233`.
- Partition: `bigmem`.
- Requested memory: 500G.
- MaxRSS reported: 11,317,496K.
- Elapsed time: 4m 23s.
- Aggregated transcript matrix:
  - Input shape: 744,997 cells by 116,918 transcripts.
  - Input nonzeros: 797,002,207.
  - Output shape: 20 cluster-name pseudobulks by 116,918 transcripts.
  - Output nonzeros: 1,445,078.

Final outputs:

- `/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/cluster_count_intron`: 167,800 lines.
- `/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/cluster_count_exon`: 196,748 lines.
- `/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/cluster_refined_cluster`: 27,414 lines.
- `/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/cluster_ratio_count`: 27,414 lines.

## 2026-07-08 Microglia Cluster x DX x Mouse Pseudobulk Run

Reran tealeaf with pseudobulks defined by:

- `cluster_name`
- `DX`
- `sample`

Here `sample` is treated as the mouse/individual identifier. The metadata join
used `sublibrary`, `rnd1_well`, `rnd2_well`, and `rnd3_well` as keys between
the combined AnnData object and `GSM5693472_cell_metadata.txt.gz`.

Input mapping:

- `/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/barcodes_to_cluster_dx_mouse.csv`
- 51,294 annotated cells joined to metadata.
- 936 `cluster_name x DX x sample` groups.
- 88,920 barcode/group rows after dropping 6,016 barcodes with conflicting
  group labels.
- Group sizes range from 2 to 1,592 barcode rows, median 38.

Added Slurm recipe:

- `extra_scripts/run_microglia_cluster_dx_mouse_tealeaf.sbatch`

Run details:

- Slurm job: `18729972`.
- Partition: `bigmem`.
- Requested memory: 500G.
- MaxRSS reported: 12,598,008K.
- Elapsed time: 41m 30s.
- Aggregated transcript matrix:
  - Input shape: 744,997 cells by 116,918 transcripts.
  - Input nonzeros: 797,002,207.
  - Output shape: 936 pseudobulks by 116,918 transcripts.
  - Output nonzeros: 26,144,309.

Final outputs:

- `/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/cluster_dx_mouse_count_intron`: 167,800 lines.
- `/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/cluster_dx_mouse_count_exon`: 196,748 lines.
- `/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/cluster_dx_mouse_refined_cluster`: 27,414 lines.
- `/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/cluster_dx_mouse_ratio_count`: 27,414 lines.

## 2026-07-09 Filtered Pseudobulk EM

Created a filtered `cluster_name x DX x sample` pseudobulk mapping requiring
at least 20 cells and at least 100,000 UMI per pseudobulk.

Filtering summary:

- Started with 936 pseudobulks and 88,920 barcode/group rows.
- Kept 551 pseudobulks and 83,858 barcode/group rows.
- Dropped 385 low-support pseudobulks.

Added Slurm recipe:

- `extra_scripts/run_microglia_cluster_dx_mouse_filtered_em_tealeaf.sbatch`

Run details:

- Slurm job: `18732826`.
- Partition: `bigmem`.
- Requested memory: 500G.
- MaxRSS reported: 27,185,384K.
- Elapsed time: 23m 15s.
- Exit code: 0.

## 2026-07-09 NNLS Quantification Backends

Implemented transcript quantification alternatives for the single-cell
equivalence-class path described in `docs/glm.tex`.

Key implementation choices:

- Added `--quant_method em|nnls|nnls_nucnorm` to `tealeaf-sc`; the default
  remains `em`.
- `nnls` solves independent non-negative least squares problems for each
  pseudobulk using SciPy's bounded sparse least-squares solver.
- `nnls_nucnorm` is a reference many-pseudobulk proximal-gradient
  implementation with singular-value thresholding and nonnegative projection.
- Transcript effective-length factors are computed from the Salmon reference
  with the existing `get_feature_weights()` convention.
- Alevin-fry `--dump-eqclasses` does not include EC/region effective lengths,
  so the GLM factor `u_s` currently defaults to `1.0` for all equivalence
  classes.

Validation:

- `python -m py_compile tealeaf/sc/sc_utils.py tealeaf/sc/tealeaf_sc.py`
  succeeded under `~/venv/jax` after loading
  `Python/3.12.3-GCCcore-13.3.0`.
- A tiny synthetic EC/transcript system produced nonnegative, normalized
  outputs for EM, NNLS, and nuclear-norm NNLS.

## 2026-07-09 Salmon/Alevin Pipeline Recipes

Pulled the useful microglia-less Salmon/Alevin pipeline pieces from the older
Calcutta workflow into this repo.

Added recipes:

- `extra_scripts/run_microglia_salmon_alevin_align.sbatch`: original
  `salmon alevin --splitseqV2 --sketch` RAD-generation step.
- `extra_scripts/run_microglia_alevin_fry_t2t_quant.sbatch`: original
  permit-list, collate, and transcript-level `alevin-fry quant
  --dump-eqclasses` step.
- `extra_scripts/run_microglia_salmon_dump_eq_weights.sbatch`: new Salmon
  `quant --dumpEq --dumpEqWeights` pass for the weighted equivalence-class
  design matrix described in `docs/glm.tex`.
- `extra_scripts/alevin_make_t2t.py` and
  `extra_scripts/alevin_dedup_t2g.py`: mapping helpers used to build
  transcript-level alevin-fry mappings.

Important caveat: `--dumpEqWeights` is exposed by `salmon quant`, not by the
documented `salmon alevin` options. The new weighted-EC recipe is therefore a
bulk Salmon pass over the same FASTQs used for Alevin RAD generation. Before
using those weighted EC rows with the existing alevin-fry cell-by-EC matrix, we
need to verify that the EC definitions can be matched or add a conversion step.
