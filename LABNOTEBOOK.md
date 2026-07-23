# Lab Notebook

## 2026-07-22 Full Factorized FISTA Rerun

The pre-FISTA all-cell rank-CV outputs are not reused. Added a fresh binary and
fixed-weighted chain over all 169,533 cells with at least 500 raw UMIs. Each
three-fold rank path starts at ranks 1, 2, 4, 8, 16, 32, 64, and 128, expands
through 512 if the selected rank reaches the upper boundary, uses 32 exact
inner FISTA steps, and allows up to 4,096 convergence-controlled epochs per
candidate. Per-fold/rank completion records are flushed to the job logs.

Selected full fits use a separate output tag, deterministic FISTA, the same
UMI threshold, and an 8,192-epoch ceiling with objective-patience stopping.
Dependent scoring uses library-normalized log1p gene expression followed by
PCA and the existing mouse-group-held-out cell-type benchmark. Binary CV and
fit jobs `19305462 -> 19305463`, weighted jobs `19305464 -> 19305465`, and
joint scoring job `19305466` all completed successfully. Removed stale scoring
job `19299644`, whose dependency could never be satisfied. The independent
pre-existing Frank--Wolfe CV jobs remain active because they do not use the
corrected factorized update.

All 24 fold-by-rank fits converged. Rank 1 failed the profile-variance
nondegeneracy criterion for both designs, while every rank from 2 through 128
passed. Mean held-out loss decreased through rank 128, but the
one-standard-error rule selected rank 32 for both binary and weighted designs.
The selected ranks were not on the search boundary, so no grid expansion was
needed. Binary rank CV completed in 9 minutes and weighted rank CV in 59
minutes, compared with the pre-FISTA runs that took more than a day and either
failed convergence selection or reached the maximum tested rank.

The selected all-cell fits contain 169,533 cells and 99,679 transcripts. The
binary rank-32 fit converged by objective patience after 204 epochs at objective
1,278.199 and completed in 2 minutes. The weighted rank-32 fit converged after
465 epochs at objective 690.490 and completed in 21 minutes. Both have all
cells active; normalized profile relative variance is 0.0305 for binary and
0.00294 for weighted.

Of the filtered cells, 78,552 have reference cell-type labels. Group-held-out
multinomial prediction from 30-dimensional log1p gene PCA gave binary accuracy
0.502, balanced accuracy 0.471, and macro-F1 0.410. Weighted gave accuracy
0.363, balanced accuracy 0.386, and macro-F1 0.309. Reference-label
silhouettes were -0.212 and -0.305, and k-means cluster silhouettes were 0.297
and 0.286, respectively. Every labeled fitted profile was active and finite.

Expanded `docs/glm.tex` with the paired-primer observation model and the exact
projected FISTA recurrences, curvature choices, factor balancing, warm starts,
and stopping rule. Added full-data validation tables for the count-selected
rank-32 fits and the paired-primer fits. Rebuilt `docs/glm.pdf` with Tectonic;
the final compile has no box, reference, or TeX warnings. The full Python test
suite passes 46 tests.

## 2026-07-22 Accelerated Factorized Convergence

Diagnosed the paired rank-64 fits that reached 2,048 iterations without
converging. Minibatch and deterministic epochs had nearly identical costs
(about 0.27 seconds binary and 0.56 seconds weighted), but the runs spent 2,016
epochs in minibatch mode and only 32 in deterministic polishing. More
importantly, the transcript-factor update divided its gradient by the number
of cells and then capped the normalized step at 0.05. This made the effective
raw-gradient step approximately 0.05 divided by the cell count rather than the
scale-correct reciprocal-curvature step.

Replaced the deterministic update with alternating nonnegative FISTA solves.
The cell subproblem uses the exact rank-by-rank Hessian, and the transcript
subproblem uses matrix-free products with \(A'A\), the cell-factor Gram matrix,
and a safety-adjusted spectral estimate. Added configurable inner steps,
deterministic fitting as the default, and factor-file continuation that checks
both cell and transcript identifiers before loading a warm start.

On the full paired binary data, 32 inner steps reduced the objective by 37.18
in 10.1 seconds, compared with 0.378 in 9.27 seconds for one inner step. This is
about 90-fold more objective reduction per second. Continuing the saved binary
fit converged after 1,800 exact epochs and 9.5 minutes of GPU optimization,
with final objective 60.019 versus 143.140 before continuation. The complete
46-test suite passes.

The weighted fit decreased from objective 151.950 to 58.911 in the first 2,048
accelerated epochs. It was still making progress near the strict \(10^{-5}\)
relative-objective threshold, so a second warm-started block was run without
relaxing the tolerance. It satisfied the ten-epoch patience rule after another
264 epochs at objective 58.753. The two weighted blocks took about 12 minutes
end to end. An obsolete weighted factor-rank CV process using the superseded
solver was stopped; the independent Frank--Wolfe CV processes remained active.

Both converged representations were evaluated on all 27,383 paired cells with
cell-type labels. Features were transcript estimates aggregated to genes,
library-normalized to 10,000, transformed with \(\log(1+x)\), and reduced to 30
principal components. Group-held-out multinomial prediction gave the binary
fit accuracy 0.688, balanced accuracy 0.636, and macro-F1 0.567. The weighted
fit gave accuracy 0.487, balanced accuracy 0.514, and macro-F1 0.439. Reference
label silhouettes were -0.239 and -0.110, respectively; k-means cluster
silhouettes were 0.352 and 0.252. Both representations had finite active
profiles for every labeled cell. These results replace the scores from the
unconverged 2,048-iteration fits.

## 2026-07-21 Paired Primer GLM

Implemented a paired-observation GLM for the poly(dT)- and random-hexamer-
primed halves of one biological cell. Each half is normalized independently,
then the response and design are stacked with equal one-half weights. The two
halves therefore share one genome-wide coefficient vector while using separate
primer design matrices. The initial implementation retains complete pairs only
and requires at least 500 UMIs in each half.

Added a reusable pair-table interface to the package and kept the AnnData
translation in a dataset-specific script. The AnnData table has 109,275 source
records but only 77,218 primer pairs occur in one sublibrary record; 14,857
pair keys are reused across sublibraries and are excluded because the pooled
alevin rows cannot separate those biological cells. Separate weighted designs
are built in one streaming pass
over the alevin-fry probability sidecar by grouping rows according to the
half-cell barcode. Binary and weighted paired factorization launchers use the
same response and differ only in their EC designs.

After excluding reused barcode pairs and requiring at least 500 UMIs in each
half, the full-data preparation retains 48,568 biological cells. The paired
response has 475,656 columns (two blocks of 237,828 retained ECs), 209,569,988
nonzeros, and row sums within 3e-7 of one. The paired theta design has 95,842
transcripts. Binary preparation completed in 66 seconds on a CPU node.

The weighted primer designs remain an empirical approximation. The patched RAD
and alevin-fry path preserves alignment score and transcript-position evidence,
but it does not fit separate primer-specific fragment-start distributions.
Consequently, it can capture differences in within-EC alignment evidence but
is not a complete correction for the strong poly(dT) 3-prime coverage bias. A
full correction requires primer-separated Salmon positional-bias models or a
RAD-position export consumed by a custom primer-specific design builder.

The complete 45-test suite passed on a Slurm CPU node. Tests cover grouped
probability designs, column normalization, complete-pair filtering, equal
primer weighting, and the existing scalable solvers and CV paths.

Execution jobs: primer-design streaming build `19300273`; paired binary fit
`19300274`; paired weighted fit `19300275`; dependent label/PCA/silhouette
scoring `19300276`.

All four jobs completed. The primer-specific design build took 7:22:26 and
produced 948,916 by 116,918 matrices with 784,615 nonzeros for poly(dT) and
883,909 for random hexamer. Binary and weighted rank-64 fits took 10:21 and
35:22, respectively. Both used 2,048 iterations, including 32 deterministic
polishing iterations, and both hit the iteration limit without satisfying the
1e-5 objective tolerance. Their final relative objective changes were about
9e-5, so the representations are provisional.

Scoring used 27,383 labeled retained pairs. Binary achieved accuracy 0.0828,
balanced accuracy 0.0908, macro-F1 0.0451, and reference-label silhouette
-0.395. Weighted achieved accuracy 0.0994, balanced accuracy 0.1205,
macro-F1 0.0635, and reference-label silhouette -0.302. Weighted therefore
improves all supervised metrics and label silhouette relative to paired binary,
but both are weak and neither fit has converged. Rank 64 was not selected by
paired-data CV.

## 2026-07-21 Primer Half-Cell Audit

The current alevin cell-by-EC matrix contains poly(dT)-primed and random-
hexamer-primed half-cells as separate rows. Its 547,077 row names do not match
the biological-cell index in the standard AnnData object; 91,940 match its
annotated `CB_polydT` values and 91,475 match `CB_ranhex`. The GLM does not
currently pair these rows or include primer type as a covariate. It normalizes
each retained half-cell response to unit mass, fits both primer types in one
shared low-rank model, and partitions molecule counts within each half-cell for
CV.

At the 500-UMI threshold, the 169,533 retained GLM rows include 79,423
annotated poly(dT) barcodes, 71,110 annotated random-hexamer barcodes, and
19,000 rows not matched to either annotation barcode set. Across 109,275
annotated biological-cell metadata records, 78,141 have both halves passing,
18,166 have only poly(dT), 9,220 have only random hexamer, and 3,748 have
neither. Median observed depth is 2,057 UMIs for poly(dT) and 1,935 for random
hexamer; among retained halves the medians are 2,566 and 2,971, respectively.

The label files assign paired barcodes the same cell-type and
cluster-by-diagnosis-by-mouse group. Grouped classification therefore keeps
paired halves in the same mouse fold, but scoring still treats each retained
half as a separate observation. Earlier pseudobulk aggregation sums both
barcodes into the same cluster-by-diagnosis-by-mouse group; the current
single-cell GLM does not perform that merge.

## 2026-07-21 Cell-Minibatch Solver Optimization

Implemented a reusable normalized sparse-data context for the scalable GLM
solvers. Cell-by-EC response rows are normalized once, split into cached sparse
blocks, and either retained on CUDA or staged from pinned host memory. Warm-start
CV paths reuse the same prepared training context instead of rebuilding and
transferring a fold for every candidate.

The unregularized direct factorization now performs seeded random-reshuffling
cell-minibatch updates of both cell and transcript factors. It reserves a final
deterministic full-gradient polishing phase, and only that phase can satisfy the
objective-patience convergence rule. Factorized ADMM keeps its global
transcript split exact: cell primal/copy/dual rows are updated in blocks,
sufficient statistics are accumulated across the complete epoch, and the
transcript primal/copy/dual and adaptive rho are updated once per epoch.

Both solvers now accumulate `C.T @ V` in EC space and apply `A.T` once per
exact epoch instead of once per cell block. Exact losses are evaluated from
`V.T @ V` and `C.T @ V`, repeated `A @ U` products are reused, and ADMM
residual reductions remain on-device until the epoch boundary. Adaptive rank
CV retains fold warm starts on CPU and evaluates only newly added larger ranks.

CPU validation passes the full 43-test suite, including normalized
cache reuse, sufficient-statistic loss equivalence, prepared-context fitting,
deterministic polishing, and incremental rank expansion. CUDA validation and
representative batch-size profiling are queued until an existing GPU analysis
releases an allocation; the four running analyses were left undisturbed.

The single-cell fitting entry point now applies an optional raw-UMI threshold
before constructing solver state and writes factors/barcodes only for retained
cells. The selected-fit launcher uses the same 500-UMI population as CV,
supports either a selected direct-factor rank or a selected factorized-ADMM
penalty, and propagates the benchmark-selected cell batch size.

Queued CUDA validation as job `19297650`. Binary and weighted representative
old-versus-new benchmarks are jobs `19297651` and `19297652`; both depend on
CUDA validation. Dependent job `19297657` will select the fastest batch size
within 1% of the best short-run objective, submit four optimized all-cell CV
jobs, then submit selected full fits and log-gene-PCA label and silhouette
scoring. Existing long-running CV jobs remain active and were not cancelled.

The first CUDA validation reached both resident solver fits, then failed when
the pinned backend called `record_stream` on a sparse CUDA wrapper. PyTorch
does not implement that operation for the sparse backend. Stream lifetime
tracking now records the COO index and value tensors plus the dense nonempty
mask instead. Focused CPU validation passed 37 tests, and replacement CUDA job
`19299162` completed successfully with identical resident and pinned objectives
for direct factorization and factorized ADMM. Binary and weighted jobs
`19299163` and `19299164` completed the optimized benchmark phase and wrote
valid results, but their legacy phase could not import the old worktree because
it was local to the login node. The stale selectors and failed-dependency jobs
were cancelled. The old commit is now on shared storage; legacy-only jobs
`19299625` and `19299626` are running, with replacement selector `19299627`
waiting on their successful completion.

## 2026-07-20 UMI-Filtered Single-Cell CV

Added raw per-cell UMI filtering to the reusable GLM CV preparation and
sampling path. The threshold is evaluated from the unfiltered deduplicated
cell-by-equivalence-class matrix, before globally rare equivalence classes are
removed. This keeps the cell-depth definition independent of the GLM design
filter.

The tuning command now accepts `--min-cell-umis`; `--cells 0` uses every cell
meeting that threshold. Reports distinguish total, threshold-eligible, and
sampled CV cells, and convert the selected dimensionless multiplier using the
scale of the full eligible population. The dataset launcher exposes these
settings through environment variables and gives large CV runs a four-day
ceiling.

For the microglia-less data, 169,533 of 547,077 cells have at least 500 raw
deduplicated UMIs. The planned rerun uses all 169,533 cells, three molecule
count folds, no cell-type labels during selection, and the existing adaptive
grids for binary and fixed-weighted penalized Frank--Wolfe and factorized ADMM.

Validation: the full test suite passed 33 tests.

Submitted four independent L40S jobs:

- `19282598`: binary penalized Frank--Wolfe.
- `19282599`: fixed-weighted penalized Frank--Wolfe.
- `19282601`: binary factorized ADMM.
- `19282600`: fixed-weighted factorized ADMM.

Reports use the suffix `_minumi500_all.json`; all four jobs started running.
These first submissions were cancelled after about two minutes when inspection
showed that CV candidates were initialized independently. Fold-local
continuation is now implemented: factorized ADMM traverses lambda multipliers
from largest to smallest and carries nonnegative primal factors while resetting
split copies and duals; penalized Frank--Wolfe traverses tau multipliers from
smallest to largest and carries its nonzero atoms into the enlarged nuclear
ball. Adaptive expansion replays the full candidate path to preserve this
ordering when a new endpoint is added. Reports record the path and warm-start
rank for each fold/candidate fit.

The continuation implementation and full suite passed 36 tests. Submitted the
corrected all-eligible-cell jobs with separate output suffix
`_minumi500_all_warmstart.json`:

- `19282628`: binary penalized Frank--Wolfe.
- `19282629`: fixed-weighted penalized Frank--Wolfe.
- `19282631`: binary factorized ADMM.
- `19282630`: fixed-weighted factorized ADMM.

All four corrected jobs started running.

Removed factor L2 regularization from the separate direct `factorized` solver.
It now minimizes only EC reconstruction loss subject to nonnegative factors and
a rank cap. Prediction-preserving column-norm balancing controls factor scale.
Added molecule-count CV for rank with an ascending warm-start path: existing
columns are retained and small positive columns are appended. Selection uses
the smallest converged, nondegenerate rank within one standard error of the
minimum held-out loss. The initial rank grid is 1, 2, 4, 8, 16, 32, 64, and
128; a selected upper boundary expands to 256 and replays the path. Cell-type
labels are not used.

Validation: the full suite passed 39 tests; the rank-CV launcher and Python
entry point also passed shell and bytecode checks.

Submitted all-eligible-cell direct-factorization rank CV:

- `19282971`: binary EC design.
- `19282972`: fixed-weighted EC design.

Both jobs started running and write reports with suffix
`_minumi500_all_rankcv.json`.

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

## 2026-07-11 Genome-Wide Torch GLMs

Added reusable Torch-based genome-wide EC GLM solvers in
`tealeaf/sc/glm_solvers.py`:

- `admm` is a dense, bounded convex nuclear-norm reference solver.
- `admm_factorized` uses factor-sized ADMM state and is non-convex.
- `frank_wolfe` stores a convex combination of nonnegative rank-one atoms.
- `factorized` fits a nonnegative fixed-rank variational nuclear-norm model.

The scalable methods stream sparse cell-by-EC blocks and retain only
transcript-by-rank and cell-by-rank factors. They avoid forming a global
transcript-by-cell abundance matrix during fitting. `single_cell` mode writes
thresholded sparse output chunks, factors, and diagnostics; it intentionally
does not invoke the pseudobulk intron clustering path.

Torch is an optional package dependency (`tealeaf[glm]`). The GPU execution
path uses the project Torch environment and the provided GPU submission
wrapper. Dense ADMM is capped and is not included in the large single-cell
submission default.

Static validation and the Torch synthetic suite passed in the project Torch
environment. The site PyTorch module was unsuitable because it referenced a
missing system RDMA library. The selected environment's CUDA runtime was
checked with an allocated-device tensor operation before full-data submission.

Updated docs/glm.tex to distinguish the bounded dense ADMM reference from the
streamed genome-wide factorized solvers, and to document raw-cell output
chunks, factors, diagnostics, and the current EC compatibility construction.

Added a phi-or-theta regularization target. The response remains the EC-count
proportion in each cell. Phi uses the inverse-effective-length design; theta
applies the linear reparameterization A_theta = A_phi diag(effective_length)
and regularizes the resulting molecular-abundance coefficient matrix without
a simplex constraint.

Source inspection clarified that gene_eqclass_probs.tsv.gz contains
cell/UMI-specific alignment-likelihood vectors computed from alignment scores
and transcript-end positions. Alevin-fry uses these vectors inside weighted
EM. Tealeaf does not yet consume them because its GLM currently assumes one
shared EC-by-transcript design matrix across cells.

Validation covered the identity A_theta = A_phi diag(effective_length), invalid
target rejection, and a theta-parameterized nuclear-norm fit recovering an
identity-design normalized response when the penalty is zero.

Added fixed EC design choices for regularized methods. `binary` uses
column-normalized EC membership. `weighted` averages the dumped per-UMI
likelihood vectors within each global EC, fills ECs without vectors uniformly,
then column-normalizes over ECs for each transcript. The weighted matrix is
cached because its sidecar is large. The comparison launcher schedules cache
construction once before weighted GPU fits and runs binary and weighted
designs as separate jobs.

The full fixed weighted design has shape 948,916 by 116,918 and 600,236
nonzeros, compared with 14,974,929 nonzeros in the binary design. Its support
is a subset retaining 4.0% of binary-compatible entries, reflecting exact-zero
alignment likelihoods. After column normalization, weighted versus binary has
cosine similarity 0.919 and relative Frobenius difference 0.408. This is a
large enough design change to compare fitted genome-wide models rather than
treat the two constructions as interchangeable.

The first full single-cell output attempt exposed two scaling issues. Writing
reconstructed transcript chunks produced multi-gigabyte files for individual
cell blocks, so genome-wide runs now write compact factors first and skip
transcript chunks unless explicitly requested. Factor manifests report finite
status, factor norms, and the fraction of cells with an active representation.
The initial factors are scaled by both rank and transcript count so their
predicted mass matches the unit-sum response; the previous rank-only scaling
caused factorized ADMM cell factors to collapse to zero.
The first rescaled full-data retry showed that a fixed factor learning rate
could still cross the nonnegative boundary. Factorized and factorized-ADMM
updates now cap each step by a Lipschitz bound. Cell updates use the exact
rank-sized curvature matrix, while transcript updates combine cell-factor
curvature with a sparse power estimate of the design spectral norm.
The comparison runs use all 64 configured iterations. An initial Frank-Wolfe
fit stopped after two iterations because the float32 objective rounded to an
unchanged value; its held-out macro-F1 was 0.029 and adjusted Rand index was
near zero. Objective-change stopping is therefore disabled for the full
method comparison so rounded losses cannot terminate an unfitted model.

All six 64-iteration fits and the dependent label benchmark completed. Every
fit had finite factors and active representations for all matched cells. The
fixed weighted Frank-Wolfe fit recovered the reference labels best: held-out
accuracy 0.154, balanced accuracy 0.156, macro-F1 0.110, adjusted Rand index
0.024, and normalized mutual information 0.075. Binary Frank-Wolfe was second
among the Frank-Wolfe designs by balanced accuracy (0.072), while all four
factorized or factorized-ADMM fits were near the 0.05 balanced-accuracy chance
level (0.053--0.056).

Reference-label silhouettes were negative for every fit. Binary and weighted
factorized fits were closest to zero (-0.024 to -0.020); weighted Frank-Wolfe
was -0.212 and binary Frank-Wolfe was -0.252. The corresponding k-means
silhouettes were positive (0.108--0.528), but ARI remained at most 0.024. Thus
the factor spaces contain cluster structure, especially for Frank-Wolfe, but
that structure aligns weakly with the standard cell-type labels. The weighted
EC design improves supervised recovery for Frank-Wolfe but does not yet
recapitulate the standard cell types well in absolute terms.

The 64-iteration comparison is retained as an intermediate result but is being
superseded by convergence-controlled fits. Factorized methods now require at
least 50 iterations and ten consecutive relative-objective changes below
`1e-6`; factorized ADMM also checks relative primal and dual residuals.
Frank-Wolfe uses its relative duality gap with the same minimum and patience,
ten power iterations per approximate oracle, and a separate 256-atom capacity.
The initial convergence-controlled run used a 512-iteration ceiling. Binary
factorized and both factorized-ADMM fits reached that ceiling without meeting
their stopping criteria: the binary factorized objective was still changing
by about `1.7e-5`, while ADMM relative residuals remained between about
`2e-4` and `8e-4`. These three fits are therefore rerun with a configurable
2,048-iteration ceiling. Diagnostics record objective changes,
ADMM residuals, Frank-Wolfe gaps and line-search steps, and an explicit
convergence reason.

Added reference-label scoring for completed fits. The benchmark aligns cell
factors to the standard-analysis `cluster_name` labels, uses five-fold
stratified splits grouped by mouse, and reports accuracy, balanced accuracy,
and macro-F1. It also reports unsupervised adjusted Rand index and normalized
mutual information from mini-batch k-means with the reference number of cell
types. Active finite cell coverage is reported so a collapsed fit cannot
receive an apparently valid score.
Each factor representation is also standardized and projected into up to 30
principal components. The benchmark reports silhouette coefficients for both
the standard-analysis labels and fitted k-means clusters in this PCA space.
Silhouette uses a fixed 10,000-cell subsample because exact pairwise distances
over all labeled cells are quadratic, while PCA and the other metrics continue
to use all matched cells.

The preceding label benchmark used PCA of standardized cell factors, not
normalized log expression. Those scores are retained as factor-space
diagnostics but are not the final standard-analysis comparison. Added a
streamed fitted-expression benchmark that sums transcript loadings to the gene
IDs in the standard AnnData, normalizes each fitted cell to 10,000, applies
`log1p`, selects 2,000 variable genes per fit, and computes 30-component
incremental PCA. All supervised and unsupervised label metrics now use this
log-gene PCA. Reconstruction uses GPU matrix multiplication in cell blocks;
only gene moments and PCA blocks move to CPU memory.

A real-factor probe found that 62,344 fitted transcripts map into the selected
standard gene universe and contribute nonzero loadings to 21,107 of 21,134
genes. A 100-cell blocked reconstruction produced finite PCA coordinates and
matched dense synthetic normalization, moments, and PCA geometry in tests.

The convergence-controlled binary Frank-Wolfe fit is rank one in practice:
only one stored atom has nonzero mass. Its cell factor changes only the fitted
library size, so normalization to 10,000 cancels all cell-to-cell differences.
Across sampled normalized profiles, the largest gene range was about `1e-7`,
which is float32 noise. Incremental PCA divided its similarly tiny component
variance by an unstable near-zero total variance and reported an impossible
explained-variance ratio. Scoring now combines block moments with a centered,
numerically stable variance calculation, rejects normalized representations
whose total between-cell variance is negligible relative to their mean-square
abundance, and computes the PCA ratio against that stable variance. The binary
Frank-Wolfe label scores are therefore invalid rather than evidence of a weak
but nonzero cell-type signal. The solver's zero approximate dual gap after its
first atom remains a limitation of its nonnegative linear oracle, not reliable
evidence that the full optimization problem is solved.

Implemented a separate `frank_wolfe_penalized` method rather than changing the
legacy fit. It optimizes squared reconstruction loss plus a smooth squared
penalty on negative coefficient entries over an ordinary nuclear-norm ball.
The signed leading-singular-vector oracle is applied to the true penalized
negative gradient. Matrix-vector products use sparse count/design operations
and reconstruct only cell blocks needed for the negative part; they do not
create the dense transcript-by-cell gradient blocks used by the legacy oracle.
The penalty multiplier is scaled by the estimated squared design spectral
norm. Output diagnostics include negative Frobenius and absolute-mass
fractions, oracle residuals, candidate gaps, and an explicit flag that the
power-iteration candidate gap is not an exact certificate. Signed fitted gene
abundances are rectified before standard-like library normalization, with row
totals computed after rectification.

The first full penalized-Frank--Wolfe comparison used 128 atoms, ten power
iterations, and a negative-mass multiplier of one. Both binary and weighted
fits completed in about 2.5 hours and retained variable normalized gene
profiles. Neither converged: both reached atom capacity, with final relative
candidate gaps of 0.326 (binary) and 0.643 (weighted). The binary fit retained
0.74% negative absolute mass and the weighted fit retained 0.014%, so the
penalty is adequate for weighted coefficients but only approximate for binary.

Despite incomplete optimization, both fits improved recovery of the standard
cell labels. Binary penalized Frank--Wolfe achieved mean accuracy 0.290,
balanced accuracy 0.309, macro-F1 0.231, ARI 0.0167, and NMI 0.0557. Weighted
penalized Frank--Wolfe achieved accuracy 0.224, balanced accuracy 0.226,
macro-F1 0.164, ARI 0.0152, and NMI 0.0646. Reference-label silhouettes remained
negative at -0.156 and -0.200. Binary therefore performed better on supervised
label recovery, while weighted had slightly higher NMI and much tighter
nonnegativity. These results support the signed nuclear-ball oracle over the
legacy clipped oracle, but more atoms or a continuation fit are required before
comparing converged objectives.

Added residual-balanced adaptive `rho` to dense and factorized ADMM. Every ten
iterations, a primal residual more than ten times the dual residual doubles
`rho`, while the reverse imbalance halves it. Scaled dual variables are
rescaled by the old/new `rho` ratio, and manifests retain the full update
history. Statistical regularization remains controlled only by `lambda`.

Added reusable count-fold cross-validation for `lambda` and `tau`. Integer UMI
counts are partitioned into three folds, models are fit to two folds, and
normalized EC reconstruction is evaluated on the held-out molecules using the
same design and cells. Lambda candidates are fractions of
`||A.T @ C||_2`. Tau candidates are multiples of the best rank-one line-search
scale from zero, `sigma_max(A.T @ C) / ||A @ u||^2`. These matrix scales are
estimated with streamed power products. Tuning uses a reproducible cell subset;
the selected dimensionless multiplier is combined with a scale recomputed on
all cells before the final fit. Dataset-specific Slurm scripts run binary and
weighted tuning for adaptive factorized ADMM and penalized Frank--Wolfe, launch
the four dependent all-cell fits, and then run the common log-gene PCA scoring
benchmark.

Corrected the factorized solvers so cell averaging applies to the complete
left-factor gradient, including regularization and ADMM terms. Previously only
the data term was averaged, making the effective penalty grow with the number
of cells and preventing a subset-selected lambda fraction from transferring to
the full dataset. Production fitting refuses a CV optimum on a grid boundary
that does not have a theoretical endpoint so that the grid must first be
widened. Zero is an explicit baseline, and lambda/lambda-max equal to one is a
valid terminal candidate because it is the zero-solution threshold.
Validation loss ignores cells with no molecules in that held-out count fold,
rather than treating their unavailable response as an observed all-zero row.
The all-cell binary smoke calculation estimated lambda-max at approximately
2.18e5. The production ADMM grid therefore includes zero and fractions from
1e-9 through 1e-5; larger initial fractions overwhelm the current low-rank
initialization. Automatic expansion starts if 1e-5 wins.

CV now expands an open grid boundary automatically. It evaluates only one new
candidate per round, multiplying an upper endpoint or dividing a lower endpoint
by a configurable factor, and merges those fold losses with the existing
results. Lambda fractions are capped at the closed lambda-max endpoint. Reports
record each expansion and whether the configured expansion limit was exhausted;
the all-cell fit is blocked only if the optimum remains on an open boundary.

For constrained penalized Frank--Wolfe, selection now uses the
one-standard-error rule: among candidates within one standard error of the
minimum fold-mean loss, choose the smallest tau, which is the most regularized
model. Candidates must satisfy the stopping rule in every fold. Because the
penalized signed oracle's candidate gap is not a valid certificate, it is no
longer part of stopping; convergence uses patient relative objective change.
CV uses a 1e-4 tolerance and up to 1,024 atoms, while the selected all-cell fit
allows 2,048 atoms with the same tolerance.

The convergence-aware CV and dependent full fits completed. All seven FW
candidates converged in every fold under objective-patience stopping. The raw
minimum reconstruction loss remained at multiplier 256 for both designs, but
the one-standard-error rule selected 0.25 for binary and 0 for weighted. Both
selected all-cell fits converged in 50 iterations. The binary representation
collapsed to relative between-cell gene variance 4.91e-14; weighted had no
active cells because tau was zero. Neither could be scored for cell-type label
recovery.

The label benchmark therefore rejects reconstruction-loss one-standard-error
selection for this application. Earlier non-CV penalized FW remains strongest:
binary achieved balanced accuracy 0.309 and macro-F1 0.231, while weighted
achieved 0.226 and 0.164. CV-selected ADMM was valid but weak (balanced
accuracy 0.0547 binary and 0.0560 weighted) and both all-cell fits stopped at
4,096 iterations without satisfying residual convergence. Selecting tau must
include a non-collapse constraint or a cell-label-aware secondary criterion;
fold-to-fold baseline variation makes the conventional one-standard-error
threshold too permissive here.

Added a scale-invariant non-collapse constraint to CV. For each candidate, the
solver reconstructs the 512 transcripts with the largest loading norms in cell
blocks, rectifies abundance, normalizes every active cell to a common library
size, applies log1p, and computes relative between-cell profile variance. A
candidate is eligible only if every fold has at least 90% active cells and
relative variance above 1e-6. This rejects zero fits and rank-one solutions that
encode only library size before applying the one-standard-error rule. The same
diagnostics are written to full-fit manifests.

The non-collapse FW rerun completed. Binary excluded multipliers 0 and 0.25
and selected 1.0; weighted excluded 0 and selected 0.25. Both all-cell fits
converged in 50 iterations with all cells active. Their normalized profile
relative variances were 0.00324 and 1.51e-5, and downstream log-gene relative
variances were 0.00414 and 2.19e-5, so both passed scoring. Binary achieved
balanced accuracy 0.182 and macro-F1 0.111; weighted achieved 0.168 and 0.112.
This fixes degenerate selection but remains below the earlier larger-radius
penalized FW fits (0.309/0.231 binary and 0.226/0.164 weighted). Reconstruction
one-standard-error selection therefore still over-regularizes relative to the
cell-type recovery objective; a label-aware secondary selection criterion is
needed if label recovery is the target.

Rather than tune directly on those labels, which would compromise their role as
an external evaluation, selection now uses profile-variance retention within
the reconstruction one-standard-error set. It finds the maximum mean normalized
profile variance among statistically admissible, converged, nondegenerate
candidates; retains candidates with at least 90% of that variance; and chooses
the smallest tau among them. This is a two-objective rule: near-optimal held-out
count reconstruction plus preservation of cell-to-cell structure. Applied to
the completed folds, it selects multiplier 16 for binary and 64 for weighted.

The profile-variance-retention rerun selected those multipliers. Both all-cell
fits converged in 50 iterations with all cells active. Binary achieved balanced
accuracy 0.280, macro-F1 0.201, and reference-label silhouette -0.168; weighted
achieved 0.217, 0.155, and -0.180. These recover most of the earlier
larger-radius penalized FW performance (0.309/0.231 binary and 0.226/0.164
weighted) while choosing parameters without using cell labels. Weighted label
silhouette improved relative to the earlier fit (-0.180 versus -0.200). This
supports variance retention over strict reconstruction minimization and
conventional one-standard-error selection for representation-oriented fitting.

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

## 2026-07-09 Alevin-Fry Weighted EC Source Inspection

Checked upstream `alevin-fry` and Salmon source to see whether the GLM weights
can be recovered directly from the existing single-cell output.

Findings:

- `alevin-fry quant --dump-eqclasses` writes `geqc_counts.mtx` and
  `gene_eqclass.txt.gz` from a global gene/transcript EC map, but the writer
  only stores label sets plus EC ids.
- Standard short-read single-cell records dispatch through
  `BasicEqClassPayload`, so the EC payload has counts only.
- `OptionalAlignmentExtras` returns no alignment-score/position extras for
  `AlevinFryReadRecordT`, `AlevinFryReadRecordWithPositionT`, and
  `MultiBarcodeReadRecordT`. Probability vectors are currently only available
  for `ScLongReadRecordT`.
- Salmon v1.10 `--dumpEqWeights` is implemented in the bulk quantification
  writer; the separate Alevin single-cell writer dumps compatibility sets and
  barcode/UMI count details, not the bulk combined weights.

Consequence: for the existing short-read Split-seq alevin-fry output, there is
no hidden weighted EC sidecar to request from alevin-fry. A conversion/matching
step against a separately generated Salmon weighted EC dump is required, unless
Salmon/Alevin is modified upstream to write the relevant conditional weights
into the RAD or single-cell EC output.

Added `extra_scripts/check_eqweight_compatibility.py` to compare Salmon
weighted EC rows with alevin-fry EC rows by transcript-set key and report
missing or ambiguous matches.

## 2026-07-09 Weighted RAD Patch Probe

Implemented and built a local Salmon v1.10.3 patch that writes the fields
needed by alevin-fry's weighted single-cell record path during RAD generation.

Patch:

- `external_patches/salmon-v1.10.3-weighted-rad.patch`
- Applied locally in `/gpfs/commons/home/daknowles/projects/salmon-v1.10.3`.
- Built binary:
  `/gpfs/commons/home/daknowles/projects/salmon-v1.10.3/build/src/salmon`

Behavior:

- Use `salmon alevin --rad --splitseqV2`, not `--sketch`.
- The selective-alignment path computes alignment scores and positions; the
  sketch path does not.
- The patched RAD records declare and write five alignment-level tags:
  `compressed_ori_refid`, `as`, `start`, `end`, and `tlen`.
- Those tags match current alevin-fry's `ScLongReadRecord` parser, which routes
  quantification through `LongReadEqClassPayload`.

Validation:

- `salmon --version` reports `salmon 1.10.3`.
- A 10,000 read-pair Split-seq probe generated RAD successfully with
  `--rad` and no `--sketch`.
- The probe RAD header parsed as:
  file tags `cblen`, `ulen`; read tags `b`, `u`; alignment tags
  `compressed_ori_refid`, `as`, `start`, `end`, `tlen`.
- Current alevin-fry source built with Rust 1.97.0 and recognized the patched
  RAD as `long read single-cell RNA-seq`.
- `alevin-fry generate-permit-list`, `collate`, and `quant --dump-eqclasses`
  completed on the probe.

Remaining gap:

- Current alevin-fry uses the probability payload internally for the long-read
  route, but its standard `--dump-eqclasses` output still writes only EC labels
  and count matrices. If the GLM needs the per-molecule/per-UMI probability
  vectors on disk, add an alevin-fry dump from the in-memory
  `LongReadEqClassPayload` before collapsing the per-cell EC map.

## 2026-07-09 Alevin-Fry Probability Dump Patch

Implemented a local alevin-fry v0.16.0 patch for the weighted RAD route.

Patch:

- `external_patches/alevin-fry-v0.16.0-dump-eq-probs.patch`

Behavior:

- When `alevin-fry quant --dump-eqclasses` runs on a long-read/probability
  payload record type, it writes `gene_eqclass_probs.tsv.gz` next to
  `gene_eqclass.txt.gz` and `geqc_counts.mtx`.
- The sidecar columns are `cell_idx`, `eqid`, `umi_rank`, and `probs`.
- `gene_eqclass.txt.gz` remains the source of `eqid -> transcript/gene label`
  definitions; the probability sidecar stores only the probability vector.
- The small-cell shortcut remains enabled for ordinary count-only records, but
  is bypassed for probability payloads so the EC payload is populated.

Validation:

- Patched alevin-fry built successfully with Rust 1.97.0.
- A 10,000 read-pair weighted-RAD probe produced nonempty probability rows with
  `-r parsimony-em`.
- The same probe with `-r cr-like-em` produced EC counts but no probability
  rows, because that resolution path does not use the probability-aware
  parsimony graph machinery. Full weighted runs therefore use
  `parsimony-em`.

Full-data plan:

- Generate weighted RAD with the patched Salmon v1.10.3 binary.
- Run patched alevin-fry with known barcodes, `--min-reads 100`,
  `-r parsimony-em`, `--use-mtx`, and `--dump-eqclasses`.
- Run tealeaf on filtered cluster x diagnosis x mouse pseudobulks for all
  current GLM variants: `em`, `nnls`, and `nnls_nucnorm`.

Full-data status:

- Weighted RAD job `18755367` completed after 16:57:32 and wrote `map.rad`.
- The first alevin-fry job failed immediately with `Illegal instruction`
  because the binary had been built with `target-cpu=native` on a newer CPU
  node and then scheduled on an older CPU node.
- Switched the dataset script default to a portable alevin-fry binary path
  under `target-portable/release/alevin-fry`.
- Added `extra_scripts/build_alevin_fry_portable.sbatch` and submitted retry
  chain `18759219 -> 18759220 -> 18759221` for portable build, weighted
  alevin-fry quantification, and all-variant tealeaf GLM.
- Retry quant job `18759220` generated the permit list and collated RAD, then
  hung after worker panics in long-read EM. The panic occurred when some
  `LongReadEqClassPayload` entries had fewer probability rows than molecule
  counts, which can happen when large connected components fall back to
  count-only `cr-like` resolution.
- Patched long-read EM to use available probability rows and fall back to
  uniform compatibility weights for missing rows, avoiding an indexing panic.
- Moved the partial failed quant output to scratch and submitted retry chain
  `18772912 -> 18772913 -> 18772914` for portable rebuild, weighted
  alevin-fry quantification, and all-variant tealeaf GLM.

## 2026-07-22 GSE233208 Parse Replication

Selected GSE233208 as a public human-brain Parse Evercode WT v1 replication
dataset. ENA exposes 40 paired snRNA-seq runs in five batches of eight
sublibraries, totaling 1.000 TB of compressed FASTQs. A streaming probe found
74 nt read-1 cDNA and 86 nt read-2 barcode/UMI reads, matching Salmon's
Split-seq v2 geometry.

Added reusable acquisition and processing components:

- `tealeaf.data.ena` queries ENA, validates paired runs, writes normalized
  manifests, and performs resumable size- and MD5-verified FASTQ downloads.
- Parameterized Salmon weighted-RAD and alevin-fry scripts replace
  microglia-specific paths.
- `tealeaf.data.alevin.merge_alevin_quantifications` merges independently
  processed sublibraries. It aligns transcript features and ECs by
  transcript identity, prefixes cell barcodes with batch, and remaps the
  patched per-UMI probability sidecars into the merged EC coordinate system.
- `tealeaf.data.parse` derives batch-aware poly(dT)/random-hexamer half-cell
  pairs from the ordered Parse RT barcode list using indexed suffix lookup,
  making pairing linear in the number of called half cells.

Each of the 40 sublibraries is quantified independently. Inspection of the
published metadata showed that combinatorial barcodes recur between
sublibraries, including sublibraries in the same biological batch. Pooling the
eight sublibraries in a batch would therefore collapse distinct nuclei before
cell calling. Merging only after alevin-fry preserves cell identity by
prefixing it with the ENA run accession without introducing separate
sublibrary-specific observation blocks into the genome-wide GLM. The merger
streams one sublibrary structure or count matrix at a time to bound peak
memory and accepts either an alevin-fry output root or its nested `alevin`
matrix directory. The resulting merged count matrix and fixed weighted design
support binary, weighted, and paired-primer genome-wide fits using the same
reusable CV and scoring code as the initial dataset.

Probability remapping precomputes transcript permutations once per local EC.
Vectors whose local and merged transcript order already agree are copied
without numeric parsing; only reordered ECs are parsed and permuted. The
merged temporary sidecar uses fast gzip compression because it is consumed
into sparse fixed-design caches.

The public processed Seurat object and case table will supply reference cell
types, subjects, and diagnoses for external representation scoring. Labels are
reserved for evaluation and are not used to select rank or regularization.

The authors' public analysis repository confirms that the processed object
contains `cell_barcode`, `Batch`, `annotation`, `subtype`, diagnosis, and case
metadata. Published batch numbering follows the original sequencing lanes,
not GEO accession order: the five internal manifest groups map to published
Batch1, Batch4, Batch2, Batch3, and Batch5, respectively. Evaluation therefore
joins each `(published batch, sublibrary, cell_barcode)` key to its ENA run
accession after applying this map. Each annotated published barcode is
expanded to its poly(dT) and random-hexamer RT barcode, producing 792,554
primer-specific reference IDs from 396,277 annotated nuclei. Fine
`annotation` labels are the primary classification and silhouette target;
case identity is the cross-validation grouping variable.

Genome-wide jobs are chained after preprocessing for binary and weighted
designs using rank-CV factorization, adaptive-rho ADMM, and penalized
Frank-Wolfe. Paired-primer binary and weighted factorized fits use independent
rank CV. All CV runs include every eligible observation with at least 500 raw
UMIs (and at least 500 UMIs in each primer half for paired fits), use held-out
molecule reconstruction rather than labels, and refit the selected value on
all eligible observations. After all selected fits finish, the common scoring
launcher reconstructs gene-level expression, applies library-size
normalization and `log1p`, computes a high-variable-gene PCA representation,
and measures donor-held-out label prediction and silhouette scores. When no
external gene universe is supplied, scoring uses all genes represented by the
transcript-to-gene map.

The corrected execution graph was submitted with 40 independent Salmon tasks,
40 dependent alevin-fry tasks, a streaming merge, all eight CV/full-fit
analysis pairs, and a final joint representation-scoring task. Job and fit
manifests are written for each submission so status checks and evaluation do
not depend on parsing console output.

A post-merge validation gate checks sparse matrix dimensions, unique cell and
feature identifiers, all expected run prefixes, the weighted probability
sidecar, complete primer pairs, cells above the UMI threshold, and overlap with
published reference labels. Downstream design construction does not start
unless this gate passes.

The fixed-weight cache stage computes the overall, poly(dT), and
random-hexamer EC probability averages in one streaming pass over the merged
per-UMI sidecar. Cells outside the primer groups are still included in the
overall design. This avoids a second decompression and parse of the largest
post-quantification file.

Representation scoring is submitted independently for each selected fit and
depends only on that fit, rather than serializing all reconstructions after the
slowest method. A final CPU task combines the per-fit summaries after all eight
scoring jobs finish.
