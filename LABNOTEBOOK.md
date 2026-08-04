# Lab Notebook

## 2026-07-27 Hierarchical Gene--Isoform Proposal

Added `docs/hierarchical.tex`, a proposed count model that decomposes
transcript abundance into gene abundance and within-gene isoform allocation.
A single cell state controls both gene log abundance and within-gene
transcript logits through separate decoder loadings. There is no independent
splicing state and no gene-specific dictionary of isoform programs. This
avoids independent cell-level EM: pooled EM initializes baseline within-gene
isoform logits, while the one low-dimensional cell state and shared decoder
loadings are fit jointly under the paired-primer multinomial EC likelihood.

The implementation sequence is gene-factor initialization from
gene-unambiguous ECs, pooled isoform-logit initialization, fitting the isoform
decoder against the complete EC likelihood, and joint refinement. The exact
likelihood evaluates transcript logits in cell minibatches but never stores a
genome-wide cell-by-transcript matrix. The note also makes explicit that an
unrestricted hierarchy is equivalent to a flat multinomial model; statistical
sharing comes from the one low-dimensional cell state and shared gene and
transcript decoder loadings.

## 2026-07-27 Factorization Representation Audit

Audited the paired positional rank-32 factorization on the same 157,006
annotated cells and donor-held-out folds used for the standard-analysis
comparison. Scoring the saved cell factor directly gave accuracy 0.578,
balanced accuracy 0.644, macro F1 0.574, and label silhouette -0.119. The
previous transcript-to-gene reconstruction followed by library normalization,
log1p, HVG selection, and PCA gave accuracy 0.493, balanced accuracy 0.564,
macro F1 0.508, and silhouette -0.303. Row-normalized and log-row-normalized
factors scored between these two results. Thus the fit contains biological
signal that the nonlinear reconstruction path discards, but postprocessing is
not the primary performance gap.

Continued the selected rank-32 fit to rank 64. Its objective decreased from
1278.731 to 1278.223, but direct-factor accuracy was 0.578, balanced accuracy
0.647, macro F1 0.575, and label silhouette -0.112. Donor-centering the
rank-32 factors also reduced accuracy to 0.560. Increasing rank and subtracting
a simple donor mean are therefore not useful fixes.

Exported a matched gene-expression PCA from the published RNA assay without
using cell-type labels. A sparse randomized PCA implementation avoids
materializing the scaled cell-by-gene matrix. On the same cells and folds,
gene PCA achieved accuracy 0.896, balanced accuracy 0.918, macro F1 0.877,
and label silhouette 0.105. The published scVI representation remains higher
at 0.953 accuracy. The large gene-PCA advantage indicates that the EC
squared-error objective spends factor capacity on transcript ambiguity instead
of preserving marker-gene totals; the negligible rank-64 gain shows that rank
alone does not explain this.

Implemented an optional multi-resolution paired-primer objective. It appends
normalized gene totals derived from ECs whose compatible transcripts all map
to one gene, weighted relative to the original EC blocks, while retaining the
same transcript-level theta and nonnegative factors. Primer-specific sampling
exposures are applied to the gene design. Count-fold cross-validation splits
raw molecules before deriving either EC or gene summaries, preventing fold
leakage. Added direct raw, row-normalized, and log-row-normalized factor
evaluation to the reusable scorer and exposed both features through the
single-cell and CV command interfaces.

The fixed \(\alpha=1\) gene-auxiliary rank-32 fit converged by objective
patience after 2,139 epochs. Gene-unambiguous assignment covered 26.7% of
retained EC identities. Direct-factor donor-held-out accuracy increased to
0.721, balanced accuracy to 0.779, and macro F1 to 0.712. ARI increased from
0.035 to 0.118 and NMI from 0.058 to 0.262, although label silhouette remained
negative at -0.138. Explicit gene reconstruction therefore recovers a large
part of the missing cell-type signal, but does not close the gap to the
standard gene PCA. The auxiliary weight is a fixed label-blind diagnostic,
not a label-selected hyperparameter.

Post-fit transforms produced a supervised-versus-clustering tradeoff. L1
normalization gave accuracy 0.712 and label silhouette -0.054. Applying
log1p after L1 normalization gave accuracy 0.688, but improved label silhouette
to 0.005, ARI to 0.188, and NMI to 0.434. Reconstructed log-gene PCA from the
same fit scored 0.698 accuracy and -0.222 silhouette, with relative gene
variance only \(7.15\times10^{-4}\). Thus post-fit gene collapsing still loses
information: the gain comes from gene-aware training, not from applying the
standard pipeline to an already compressed transcript fit.

The practical current strategy is to use gene expression for annotation and
retain transcript theta for isoform analysis. A unified next model should
factor gene abundance and within-gene isoform allocation separately, ideally
under a count likelihood with donor and chemistry covariates, rather than ask
one low-rank Gaussian EC-proportion model to serve both tasks.

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
The initial all-cell binary smoke calculation estimated lambda-max at
approximately 2.18e5 and motivated fractions from 1e-9 through 1e-5. That
range became obsolete after correcting cell scaling and changing to the
primer-specific positional design. Production CV now spans zero, 1e-4, 1e-3,
1e-2, 1e-1, and the closed lambda-max endpoint. Strong-to-weak continuation
starts from the zero-solution threshold, and zero remains the unregularized
endpoint.

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

Salmon and alevin-fry stages write completion markers only after checking their
required outputs. A repeated task reuses marked output and rejects an
unmarked directory, permitting restart-safe arrays without treating partial
tool output as complete.

While the remaining public FASTQs were downloading, early arrays were
submitted for the 30 sublibraries whose read pairs had already passed size and
MD5 verification. Preprocessing was split into disjoint 30-run and 10-run
Salmon/fry branches after scheduler estimates showed that serializing the
remaining runs behind the early branch would delay completion. The second
branch depends on completion of the download, and merge depends on both fry
branches. Early fry can overlap Salmon processing of the remaining
sublibraries without duplicate writers.

Within each branch, fry uses Slurm's corresponding-array dependency rather
than waiting for every Salmon task in the branch. Each sublibrary can
therefore enter permit-list generation, collation, and quantification as soon
as its own RAD output is complete; the merge gate still requires all fry
tasks.

Historical weighted full-dataset runs peaked near 28 GB RAM for both Salmon
and alevin-fry. Independent sublibrary tasks therefore request 48 GB for
Salmon and 40 GB for fry rather than the earlier 164/96 GB reservations.
Common CPU nodes have 28 cores: eight-core tasks allow three sublibraries per
node, while the old 20-core request allowed one. Salmon used its assigned cores
but gains aggregate throughput from run-level parallelism; fry historically
used about eight cores despite a 20-core allocation.

Per-sublibrary Salmon limits are six hours and fry limits are two hours. The
earlier four-day limit was based on whole-dataset processing and prevented
useful backfill. A prior weighted run processed roughly 5.6 billion fragments
in 17 hours on 20 CPUs; the corrected fry stage for the same data completed in
28 minutes. The first public Salmon tasks sustained millions of fragments per
minute, and even the slowest initial task projected below three hours. Six
hours leaves Salmon margin at eight CPUs, while two hours is conservative for
the smaller per-sublibrary fry tasks and improves backfill placement.

The first production Salmon completion was SRR24710697 after 28 minutes 45
seconds. It processed 387,223,905 fragments, mapped 143,163,109, and wrote a
9.0 GB RAD plus valid metadata and a completion marker. Its corresponding fry
array element was released independently as intended.

That first fry task completed quantification in 4 minutes 41 seconds but then
failed because the validation launcher relied on a `python` command inherited
through the batch environment. Resolving the module interpreter by path was
not sufficient because that worker environment lacked NumPy. Dataset launchers
now consistently use the configured analysis virtual environment, while the
reusable fry script continues to accept any explicit interpreter path. The fry
script also recognizes complete matrix and probability output without a
marker, validates it, and writes the marker. This makes failures between tool
completion and validation restartable without repeating quantification.

A reusable pre-merge gate reads the ENA manifest and requires nonempty Salmon
and fry completion markers, RAD and metadata files, and per-run validation
reports. It independently rechecks positive cells, ECs, weighted
compatibilities, and molecule counts no larger than Salmon mapped fragments.
The merge depends on this content-based gate after all array attempts finish,
so a recovered post-quantification launcher failure does not invalidate
validated data and any unrecovered task still stops the graph.

Slurm snapshots a submitted batch script, so updating the launcher did not
change already-pending fry array elements. Those elements can still complete
quantification but fail when their stored launcher reaches validation. Two
one-time recovery arrays run the corrected launcher after the original early
and remaining fry arrays finish. Complete quantifications are only validated;
missing ones are processed normally. The live content gate depends on these
recovery arrays rather than the historical array exit states.

The 30-run early branch completed with all per-run count and structure checks
passing. It contains 1,738,306 called primer half-cells and 1,799,846,896
deduplicated molecules from 4,898,240,110 Salmon-mapped fragments. The slowest
Salmon task took 2 hours 41 minutes, confirming the six-hour limit has useful
margin. The final fry task was run in the existing compute allocation rather
than waiting several days for a 40 GB backfill slot; it completed in about 12
minutes and retained 60,985,855 molecules from 206,595,909 mapped fragments.
After all 30 outputs had validated markers, the stale early fry and recovery
jobs were cancelled. The content gate now waits only for the remaining
10-run branch.

All 80 FASTQs are now present with their manifest sizes and checksums. ENA
limited individual HTTPS streams to roughly 1 MB/s. The final 13.6 GB file was
therefore fetched as eight HTTP byte ranges, assembled in byte order, and
accepted only after its total size and published MD5 matched exactly. The
single-stream partial was moved aside rather than deleted. This reduced the
last transfer from a projected tens of minutes to about 20 minutes and
released the remaining 10-run Salmon branch.

A bounded real-data integration job uses one million paired reads from a
checksum-verified GSE233208 sublibrary. It runs the same patched Salmon,
alevin-fry, probability dump, and validation components as production. This is
an interface check alongside, not a replacement for, the full 40-sublibrary
analysis.

The fixed-weight cache stage computes the overall, poly(dT), and
random-hexamer EC probability averages in one streaming pass over the merged
per-UMI sidecar. Cells outside the primer groups are still included in the
overall design. This avoids a second decompression and parse of the largest
post-quantification file.

Representation scoring is submitted independently for each selected fit and
depends only on that fit, rather than serializing all reconstructions after the
slowest method. A final CPU task combines the per-fit summaries after all eight
scoring jobs finish.

The first million-read weighted integration run exposed a count-conservation
failure: fry reported 362,310 processed records, but the dumped EC matrix
contained 3,570,371 molecules. The patched probability-bearing record path
disabled the small-cell shortcut, while the inherited cleanup still cleared
the reusable per-worker EC map only for cells above the shortcut size
threshold. Small cells therefore retained preceding cells' EC state and
re-emitted cumulative counts.

The alevin-fry patch now clears the EC map whenever `used_fast_path` is false,
independent of cell size. Rebuilding and repeating the same integration run
produced 336,260 deduplicated molecules from 362,310 processed records, with
309,373 rather than 2,956,669 nonzero cell-by-EC entries. The reusable
quantification validator accepts a maximum total molecule count, and the
bounded real-data integration job sets that maximum to its input read-pair
count. This provides a regression check for stale state or other count
inflation before a patched binary is used for full data.

Every production fry task now applies the same count-conservation gate using
Salmon's authoritative `num_mapped` value. Matrix structure, weighted
probability output, nonnegative counts, nonempty cells, and total molecule
mass are validated before the per-run completion marker is written.

The corrected real-data output was also carried through primer pairing,
weighted-design caching, paired rank CV, a bounded genome-wide factorized fit,
and representation scoring. This integration pass found two downstream
assumptions that unit-level synthetic data had not exposed:

- The design-cache command required merge-generated sparse `.npz` caches even
  though raw alevin-fry output only guarantees Matrix Market and EC text
  files. It now uses the shared format-aware alevin structure loader and works
  on either raw or merged quantifications. Paired preparation uses the same
  public structure and count loaders, including nested alevin output roots.
- Paired-primer responses have fractional entries because each primer half is
  normalized to mass 0.5. Generic molecule-fold CV correctly rejected those
  entries as non-integer but left both submitted paired CV jobs unable to run.
  Paired CV now partitions the original integer molecules within each primer
  half, then independently normalizes each training and validation partition
  to equal primer mass. Raw fold counts live in a dedicated, non-serialized
  preparation field so fit manifests remain JSON-compatible.

With these corrections, a two-fold paired weighted rank path completed on 50
real cells and warm-started rank two from rank one. A bounded rank-two fit
wrote compact factors and the scorer aligned 2,993 published annotations,
aggregated transcript loadings to genes, formed a `log1p` gene-PCA embedding,
and computed donor-held-out classification and silhouette metrics. The
bounded fit used two iterations and is only an interface test; its scores are
not biological results.

The production Python environment is built against Python 3.12. Dataset
launchers previously mixed that interpreter with a Python 3.10 module, and the
preprocessing gate loaded no Python module. This did not affect Salmon or fry,
but it caused the historical fry jobs to fail when they reached the Python
count validator despite having complete quantifications. The Python module is
now declared once in the dataset configuration and loaded consistently by
preprocessing, merge, metadata, validation, and GLM submission scripts. The
configured interpreter was checked by importing NumPy, SciPy, and PyTorch
after loading that module.

All 40 public GSE233208 sublibraries passed the content gate. Together they
contain 2,520,302 called primer half-cells and 2,185,504,840 deduplicated
molecules from 6,050,007,444 Salmon-mapped fragments. Because the queued merge
would have waited several days for a large-memory slot, the same merge launcher
was run in an existing compute allocation with sufficient memory. The GLM
submitter now supports either a pending merge dependency or an already
materialized merged count matrix, so a completed external merge can root the
validation and fitting graph without a placeholder job.

The merged quantification has 2,520,302 cell rows, 9,048,099 genome-wide ECs,
243,927 transcript features, and 1,262,210,898 nonzero cell-by-EC entries. At
the 500-UMI threshold, 568,887 cells remain; 754,660 rows overlap published
cell-type labels, and 650,297 poly(dT)/random-hexamer pairs have both halves.
The merged weighted sidecar is larger than the compressed inputs because
reordered probability vectors use full precision and global indices require
more digits.

Validation and fixed-weight cache reports can now root GLM submission directly
when both were produced outside Slurm. The submitter still uses dependencies
when merge, validation, or cache jobs are pending, and it refuses to reuse
pre-existing artifacts for a new pending merge. This avoids large-memory queue
delays without allowing stale caches to bypass a new preprocessing run.

## 2026-07-24 Primer-Specific Positional Bias

Inspection of the patched alevin-fry implementation showed that its weighted
sidecar applies one hard-coded transcript-end model to every read. It does not
learn oligo(dT) coverage and applies the same 3-prime preference to
random-hexamer cells. The grouped sidecar designs remain useful alignment-
weight comparators, but they are not a valid primer-bias correction.

The corrected path classifies Parse reads from the RT barcode in read 2 and
streams transcript read 1 to separate poly(dT) and random-hexamer Salmon
processes. Exact RT barcodes and unique one-substitution corrections are
retained. Each sublibrary is processed independently in a Slurm array, so
Salmon learns positional models from every run and the demultiplexing bottleneck
is parallel. No primer-split FASTQ copies are materialized.

Each Salmon run uses positional-bias learning and rich EC-weight output.
Salmon models observed and expected positions within transcript-length
classes, and its effective lengths include the resulting primer-specific bias
correction. Rich range-factorized rows sharing a transcript set are collapsed
across sublibraries by fragment-count-weighted averaging, joined to alevin ECs
by transcript-set identity, and normalized by transcript column. Effective
lengths are mapped-fragment-weighted across runs. Missing ECs use inverse
primer-specific effective lengths rather than uniform weights. The builder
reports exact-join coverage by both EC identity and primer-specific UMI mass.

The paired GLM now accepts a `positional` design. Phi fits use the two
column-normalized primer designs directly. Theta fits multiply each design by
its own Salmon effective-length vector before stacking, while retaining one
shared coefficient row per biological cell. This avoids reintroducing a
shared reference-length approximation after fitting primer-specific bias.

For the primary paired theta analysis, oligo(dT) UMIs are now assumed
proportional to transcript molar abundance (TPM): anchored priming supplies no
additional transcript-length factor. Random-hexamer UMIs retain relative
Salmon effective-length exposure, normalized by its median over retained
transcripts so the two separately normalized primer losses have comparable
units. The reusable API also exposes the previous effective-length model for
both primers and a TPM model for both. CV reports and fit diagnostics record
the selected model and scale constants. The paired-data API and both command
line interfaces default to this primary model. This also ensures previously
submitted Slurm launchers that did not include the newer explicit option
resolve the intended model from the live Python implementation.

Because each primer response is normalized separately, the exact ranhex mean
contains a cell-specific denominator equal to the abundance-weighted effective
length. The median-length scaling is a linear approximation to this
denominator. The oligo(dT) block anchors the coefficient scale and keeps
nuclear-norm penalties comparable; an exact treatment would require a
compositional objective or alternating primer-by-cell exposure factors.

The full-data preparation preflight showed that host memory is dominated by
overlapping sparse copies rather than the bounded solver tail. Paired
preparation now extracts the two filtered integer response blocks and releases
the 2.5-million-row merged matrix before loading primer designs. After stacking
the final compatibility matrix, it similarly releases the membership and
primer-specific design intermediates before normalizing responses. This keeps
the same matrices and values while reducing peak overlap in every CV and fit
job.

Paired count-fold CV no longer retains all raw and normalized
training/validation matrices. A deterministic reusable fold plan partitions
integer molecules with the same sequential binomial draws and regenerates one
fold at a time, which also supports adaptive rank expansion without storing
three full copies of each response. The tuning launcher avoids row-slicing
copies when the requested CV population is already the complete prepared
population.

The merged paired raw matrix stores integer counts directly. Fold-plan
iteration now copies those values straight into its one mutable `int64`
remainder array rather than rounding and comparing another full-size array on
every adaptive pass. Non-integer inputs remain validated at plan creation.

Long sparse phases now emit JSON progress records. The tuning launcher reports
completed data preparation, each paired fold reports start, elapsed time, and
training/validation nonzeros, and regularized paths report each completed
fold-candidate fit. These events distinguish slow host preparation from GPU
optimization in Slurm logs.

Penalized Frank--Wolfe CV now releases the temporary CUDA allocator cache used
to estimate its scale before constructing validation and training operators.
The requested response backend is propagated into the FW solver instead of
being resolved again internally. FW diagnostics now include backend,
sparse-cache time, epoch throughput, and peak CUDA allocation.

The full-scale binary-design CUDA preflight exercised all 201,023 biological
cells passing 500 UMIs in each primer half, 915,119,948 response nonzeros, and
two streamed count folds. It completed in 69 minutes with 64.9 GB peak host
RSS and 25.7 GB peak CUDA allocation under the 192 GB host and B6K GPU
production requests. After sparse cache construction, the two deliberately
truncated rank-one fits processed about 75,000 and 237,000 cells per second.
They ran only two iterations and were not expected to converge; this was a
scale and interface gate, not rank selection.

All 40 primer-specific positional runs are complete. The final run conserved
all 334,102,139 input pairs, assigning 267,766,447 to one primer and
classifying 66,335,692 as unknown or ambiguous. Both Salmon outputs reached
their normal inference stopping criterion and wrote completion markers,
releasing the validated 40-run design reducer.

A 100,000-pair real-data check assigned 79,022 transcript reads: 52,110 to
poly(dT), 26,912 to random hexamer, and 20,978 to unknown or ambiguous RT
barcodes. Separate Salmon runs completed with 59.3 and 60.5 percent mapping.
When joined to the existing million-read alevin smoke output, exact Salmon
rows covered 16.6 and 11.1 percent of distinct ECs but 71.0 and 61.9 percent
of molecule mass; all remaining compatible entries received effective-length
fallbacks. The full test suite passes with 68 tests, and `docs/glm.tex`
compiles after documenting the correction and parameterization.

The full public dataset run was submitted as a 40-task positional-quantification
array with at most ten concurrent tasks. Its reducer feeds separate
factorization, adaptive-rho ADMM, and penalized Frank-Wolfe count-fold CV and
full-fit branches. Each completed fit has an independent cell-type label and
PCA-silhouette scoring job, followed by one aggregate report.

The initial array requested 32 CPUs and 112 GB per task and remained pending
without starting. Historical full-sublibrary Salmon tasks used 4.2--5.7 GB
RSS and completed within 2 hours 15 minutes on eight CPUs. Before any new task
started, the untouched graph was replaced with 16-CPU, 32-GB tasks: two
six-thread Salmon processes and two two-thread decompressors run concurrently.

The right-sized array started immediately but most tasks found no pre-existing
user directory on node-local scratch and stopped before Salmon. One task
created only 8 KB of partial output before cancellation; it was moved to
scratch. The launcher now creates its node-local scratch parent. It also calls
the pigz executable by absolute path rather than loading its GCC 12 module,
which had downgraded dependencies after loading the GCC 13 Python module.

Before reduction, every run now passes an independent positional-output gate.
It requires conservation of the original raw-pair count, conservation across
the two primer assignments and rejected RT barcodes, agreement between
assigned reads and each Salmon `num_processed`, positive mapped counts, the
full reference target count, nonempty rich ECs, all four readable observed and
expected positional models, and completion markers. Reduction starts only
after all 40 reports pass. The reducer also writes one aggregate validation
summary with run IDs, conserved read totals, assignment fraction, primer-
specific processed and mapped totals, run-level mapping-rate bounds, rich-EC
totals, and reference target counts.

Bulk Salmon auto-detection consistently chose unstranded for poly(dT) but
forward-stranded for nine of ten initial random-hexamer runs. The original
alevin RAD metadata reports an unstranded library (`IU`), so accepting the
ranhex auto-detection would remove alignments that remain in the EC response.
The first wave was stopped before completing any run and its 168 KB of partial
metadata was moved to scratch. Both bulk primer runs now explicitly use
single-end unstranded type `U`, and the validation gate requires that type.

The first completed unstranded run conserved all 350,606,188 input pairs and
assigned 257,601,910 (73.5 percent) to a primer. Poly(dT) and random hexamer
mapping rates were 43.3 and 43.5 percent. Both outputs contained all 243,927
targets and valid rich EC and positional-model products. The compressed rich
EC files contained 3,451,554 poly(dT) rows and 1,269,938 random-hexamer rows.
To keep the 40-run reduction bounded, the parser now retains only transcript
sets present in the 948,916-row alevin EC universe. On the first poly(dT)
file, this reduced 3,451,554 range-factorized rows to 555,814 relevant unique
sets. Skipping weight-vector allocation for irrelevant rows reduced a
single-file profile from 220 to 142 seconds and from 5.4 to 4.9 GB peak
memory. A Slurm-requeued task left 32 KB of partial output, which was moved to
scratch before its clean retry.

All 40 validated runs conserved 13,614,523,624 input read pairs and assigned
10,712,743,461 (78.69 percent) to one primer. Salmon processed 6.89 billion
poly(dT) and 3.82 billion random-hexamer reads, mapping 48.46 and 50.41
percent, respectively. Every run retained all 243,927 reference targets and
the required unstranded library type.

The genome-wide reducer produced two finite, nonnegative
9,048,099-by-243,927 primer designs with 73,276,245 nonzeros each and 238,282
active transcript columns. Exact transcript-set joins covered 23.31 percent of
distinct poly(dT) ECs and 16.15 percent of random-hexamer ECs, but represented
97.92 and 98.86 percent of the corresponding molecule mass. The unmatched
tail therefore uses the documented effective-length fallback without
controlling the fitted signal.

For oligo(dT), the primary theta model treats expected UMI abundance as
proportional to latent transcript molecule fraction, or TPM scale, with no
additional transcript-length exposure. Anchored priming provides
approximately one UMI-generating opportunity per captured polyadenylated
molecule. Primer-specific effective lengths remain part of Salmon's
alignment, positional-bias, and fallback-weight calculations; they are not
multiplied into the oligo(dT) response model again. This does not assert equal
capture efficiency across transcripts: poly(A)-tail state, degradation,
internal priming, end annotation, and mappability can still introduce
transcript-specific effects. Random-hexamer counts retain relative
effective-length exposure because longer molecules provide more priming
opportunities.

The first complete positional-design CV branch used all 201,023 paired cells
passing 500 UMIs in each primer half. Direct nonnegative factorization
converged for ranks 1 through 128 in every count fold. Mean held-out loss
decreased from 0.006700 at rank 1 to 0.006663 at rank 128; the
one-standard-error rule selected rank 32 without using cell-type labels. The
selected all-cell rank-32 fit converged by objective patience after 273
iterations, with finite factors, effectively complete active-cell coverage,
and 26.4 GB peak CUDA allocation.

The corresponding log1p gene-PCA benchmark scored all 157,006 labeled cells.
Five-fold held-out label prediction had mean accuracy 0.4935, balanced
accuracy 0.5641, and macro F1 0.5077. The reference-label silhouette was
-0.3029 and the k-means silhouette was 0.2876. These are intermediate
factorization results; ADMM and penalized Frank--Wolfe CV, their selected
fits, and the aggregate standard-analysis comparison remain in progress.

The published Seurat object contains a 30-dimensional scVI representation
rather than a PCA reduction. A reusable reference-embedding scorer maps the
published nucleus IDs to internal poly(dT) IDs and applies the same donor-held-
out classifier and silhouette implementation after restricting to the exact
paired-GLM cell set. On the same 157,006 labeled nuclei, published scVI
achieved 0.9525 accuracy, 0.9731 balanced accuracy, 0.9528 macro F1, and
0.1920 label silhouette. This establishes a matched standard-analysis
baseline without using labels in GLM fitting or hyperparameter selection.

The first full penalized Frank--Wolfe CV attempt exhausted CUDA memory after
three fold-zero candidates. Its 16,384-cell block created a
236,543-by-16,384 fitted-value workspace of approximately 14.4 GiB; the
nonnegativity calculations held multiple such workspaces while the sparse
response cache was resident. The solver now clamps temporary negative blocks
in place and multiplies the direction workspace in place. Frank--Wolfe
launchers use 4,096-cell blocks by default, while factorization and ADMM retain
their 16,384-cell throughput setting. The CV launcher also exposes the sparse
data backend for explicit fallback if the cache-bounded rerun requires it.

Scoring failures now remain inspectable as `status=invalid` artifacts but
produce a nonzero process exit. Aggregation independently rejects every
non-`ok` row, preventing a failed fit from satisfying the evaluation graph.

The first positional ADMM signal-bearing candidate showed why one tolerance
should not control both optimization and split feasibility. At 4,096
iterations its relative objective change was 6.7e-6 and primal residual
5.6e-5, while the dual residual was 7.0e-4. The objective was stable, but the
candidate was excluded under a shared 1e-4 cutoff. Factorized ADMM now retains
a 1e-4 patient objective tolerance and uses a separate 1e-3 split-residual
tolerance for production CV and full fitting. Candidate progress records both
tolerances, final residuals, final rho, and the number of adaptive-rho updates.

Loosening only the feasibility tolerance exposed a second issue: a candidate
warm-started from the nearly zero lambda-max factors reported no float32
objective change and stopped after 50 iterations, but its held-out loss and
profile variance showed that it had not activated signal. This is the
bilinear zero saddle: both factor gradients vanish when both factors are zero.
ADMM continuation now retains the leading singular vectors computed for
lambda-max. At each lambda it forms the exact best rank-one step from zero,
with magnitude `(sigma - lambda) / ||A u||^2`, adds small random remaining
columns, and compares that seed's training objective with the preceding
candidate. The lower-objective initializer is used, so genuine continuation
is retained without allowing the strong endpoint to trap the path.

The initial selected all-cell ADMM launcher still cold-started only the
CV-selected endpoint, even though CV itself used continuation. This was
detected while the unregularized rank-64 refit remained active far longer than
the corresponding fold endpoint. The run was stopped rather than accepted.
Selected ADMM fitting now recomputes lambda-max on all eligible cells and
replays a geometric strong-to-selected path. At every stage it compares the
previous factors with the spectral descent seed, resets ADMM split state, and
records the initializer, objective, residuals, rho updates, and convergence
status in the final manifest. This makes the full refit follow the same
warm-start policy used for model selection.

The corrected continuation run completed in 10 minutes 37 seconds and all six
stages met the objective and split-residual rules, but downstream evaluation
showed that this was not sufficient evidence of optimization. Its final
objective was 1326.64 versus 1278.73 for direct factorization, mean
reconstructed library mass was 0.014 rather than approximately one, and
donor-held-out label accuracy was 0.169. The shared transcript-factor gradient
had been divided by 201,023 cells while its mathematically corresponding step
was still capped at 0.05. Split feasibility could therefore converge while
the data-fitting update remained microscopic.

Factorized ADMM now solves the cell-factor primal block exactly using one
rank-by-rank Cholesky factorization per epoch. Its shared-factor step uses the
reciprocal of the Lipschitz curvature after applying the same cell-count
scaling as the gradient. The low-accuracy fit is diagnostic only and is not
accepted as the final ADMM result; CV, selected continuation, and scoring must
be rerun with the corrected updates.

The corrected ADMM rerun selected the closed zero-penalty endpoint in all-cell
CV. Its mean held-out loss was 0.006746, compared with 0.006783 at
lambda/lambda-max 0.1 and 0.007023 at lambda-max. The selected continuation
reached objective 1294.79 with mean reconstructed library mass 0.873, fixing
the earlier mass collapse. Biological recovery remained weak: donor-held-out
accuracy was 0.164, balanced accuracy 0.251, macro F1 0.196, and label
silhouette -0.297. This is now a valid negative solver result rather than an
artifact of the cell-count step cap.

Penalized Frank--Wolfe exposed a separate adaptive-grid cost. Each upper
radius expansion replayed every completed fold and candidate, so reaching
multipliers 1024 and 4096 repeated hours of unchanged fits. Adaptive GLM CV
now retains each fold's endpoint factors and result rows when an expansion
continues the natural path. It evaluates only the new larger Frank--Wolfe
radius or smaller ADMM penalty. A path is replayed only when an expansion
changes its starting endpoint.

The same audit found that the selected all-cell Frank--Wolfe launcher
cold-started only the final radius even though each CV endpoint accumulated
atoms along the increasing-radius path. Selected Frank--Wolfe fitting now
reads the exact path from the CV report, truncates it at the selected
multiplier, and replays it on all eligible cells. The final manifest records
the warm-start rank, final rank, objective, gap, and stopping state for every
radius, making the selected refit comparable to the CV candidate.

Profiling the first incremental rerun showed the GPU dropping to zero
utilization between candidates while the same fold-specific sparse training
cache was reconstructed on one CPU core. Frank--Wolfe CV now constructs one
`SparseGLM` training context per fold, uses it for radius-scale estimation and
all candidates, and retains one separate validation context. A test verifies
that three radii over two folds create four contexts rather than rebuilding a
training context for every candidate.

The cached production Frank--Wolfe CV completed in 27 hours 29 minutes with
78.5 GB peak host RSS and approximately 51 GB peak CUDA allocation. All 33
fold-candidate fits met their configured objective-patience rule. The adaptive
path extended the original radius multipliers through 1024, 4096, 16384, and
65536 without replaying completed candidates. Mean held-out loss decreased
monotonically from 0.007023 at radius zero to 0.006762 at 65536. Mean
normalized profile variance also increased through the final radius, from
0.0833 at 4096 to 0.0955 at 16384 and 0.1074 at 65536.

The one-standard-error threshold at the final radius was 0.0067633, so no
smaller radius was eligible. The selected multiplier therefore remained the
open upper boundary after all four expansions
(`best_on_boundary=true`, `grid_exhausted=true`). This is not a valid
regularization selection. The progressively larger radius is still changing
Frank--Wolfe step size and optimization progress: every continuation stage
stopped after another 50 atoms even though its accumulated validation
improvement exceeded fold uncertainty. The pending full fit, label score, and
aggregate jobs were cancelled rather than treating the boundary endpoint as a
selected biological model. The complete CV report is retained as a negative
solver diagnostic.

The completed GSE233208 comparison therefore contains three valid,
label-blind fits or references on the same 157,006 labeled cells. The direct
rank-32 nonnegative factorization achieved 0.493 accuracy, 0.564 balanced
accuracy, 0.508 macro F1, and -0.303 reference-label silhouette. Corrected
rank-64 ADMM achieved 0.164 accuracy, 0.251 balanced accuracy, 0.196 macro F1,
and -0.297 silhouette. Published scVI achieved 0.953 accuracy, 0.973 balanced
accuracy, 0.953 macro F1, and 0.192 silhouette. Thus the primer-aware
factorization recovers substantial cell-type signal but remains far below the
standard analysis, while ADMM is a valid weak result. Penalized
Frank--Wolfe has no accepted label score because its unsupervised
hyperparameter selection did not close.

## 2026-07-27: one-state hierarchical gene--isoform model

The hierarchical branch implements a genome-wide paired-primer multinomial
model with one latent cell state. A gene decoder defines positive gene
abundance, and a transcript decoder defines a within-gene softmax over
isoforms. Their product is passed through the fixed primer-specific
positional EC designs. The implementation evaluates only observed sparse
cell--EC events and obtains each multinomial normalizer from design column
sums; it does not instantiate a cell-by-transcript matrix.

Fitting is staged. Sparse log-normalized gene PCA initializes the shared
state from gene-unambiguous EC counts. A gene multinomial stage updates the
state and gene decoder, a fixed-state EC stage learns isoform intercepts and
loadings, and joint refinement minimizes paired EC loss plus an equally
weighted gene multinomial loss. Labels are not used in fitting.

The first full rank-32 run included 201,023 cells passing 500 UMIs in each
primer half, 54,227 genes, and 236,543 transcripts. Gene loss decreased from
8.092 to 8.075. Paired EC loss decreased from 13.631 to 13.336 with the cell
state fixed, then to 12.929 during joint refinement. Ten epochs per stage
completed in 23 minutes with approximately 48 GB peak host memory.

On the matched 157,006 labeled cells with donor-held-out folds, the initial
gene PCA state achieved 0.893 accuracy, 0.911 balanced accuracy, 0.873 macro
F1, and 0.109 label silhouette. Gene multinomial pretraining improved these
to 0.915, 0.936, 0.906, and 0.144. Joint paired-primer refinement gave 0.915,
0.936, 0.906, and 0.146. The matched standard gene pipeline gave 0.896,
0.918, 0.877, and 0.105, respectively. The hierarchical state therefore
exceeds the standard pipeline on this external assessment, though the
published scVI representation remains higher at 0.953 accuracy.

During assessment, raw signed factor states were incorrectly filtered by
requiring a positive row sum. This excluded roughly half of centered PCA
states. Raw-factor scoring now defines active rows by finite nonzero norm;
the positive-total criterion remains in place for nonnegative L1 transforms.

## 2026-07-28: count-aware flat GLM ablations

The poor direct EC-factorization result did not distinguish transcript
ambiguity from differences in response scale and feature preprocessing. A
reusable flat transcript GLM now uses one signed cell state, an exponential
transcript decoder, the exact paired positional EC likelihood, and one of
three gene objectives: fixed-concentration negative binomial, Poisson, or
squared error on standardized library-normalized log counts. Gene features
are selected by transformed variance after a prevalence filter. The NB
implementation optimizes the exact likelihood with minibatch Adam; its
log-mean Fisher weight saturates at the fixed concentration in the same way
as an NB IRLS update.

The full rank-32 NB run used 201,023 cells, 2,000 selected genes,
concentration 10, and ten epochs each of gene fitting, fixed-state EC fitting,
and joint refinement. NB gene loss decreased from 0.583 to 0.575. EC loss
decreased from 13.488 to 12.832 with the state fixed and to 12.814 during
joint refinement. The run completed in 22 minutes with approximately 48 GB
peak host memory.

On 157,006 labeled cells with donor-held-out folds, the old direct EC
factorization had 0.493 accuracy, 0.564 balanced accuracy, 0.508 macro F1,
and -0.303 label silhouette. Standardized log-gene PCA initialization alone
gave 0.893, 0.911, 0.873, and 0.109. The selected-gene NB stage gave 0.910,
0.929, 0.897, and 0.133; paired-EC joint refinement gave 0.912, 0.931, 0.899,
and 0.136. The hierarchical joint model remained slightly higher at 0.915,
0.936, 0.906, and 0.146.

Two label-blind ablations showed that the gain was not specific to the chosen
NB concentration or HVG cutoff. Poisson fitting on the same 2,000 genes gave
0.911 accuracy and 0.136 silhouette. NB fitting on all 24,002 genes observed
in at least 200 cells gave 0.910 accuracy and 0.131 silhouette. A direct
standardized-log joint objective reached 0.900 accuracy and 0.123 silhouette.
Random initialization followed by 30 NB epochs reached 0.908 accuracy but
only 0.115 silhouette, compared with 0.910 and 0.133 from the transformed PCA
initialization.

The evidence therefore rejects EC ambiguity as the dominant explanation for
the old GLM gap. Normalized log-count initialization, variance scaling, and a
gene-count training objective account for most of the recovery. The explicit
hierarchical gene/isoform decomposition contributes a smaller remaining
improvement.

## 2026-07-28: microglialess count-aware comparison

The count-aware and hierarchical models were also run on the microglialess
paired-primer data. Both used the primer-specific fixed weighted EC designs,
the oligo(dT)-TPM sampling model, rank 32, and the same requirement of at
least 500 UMIs in each primer half. This retained 48,568 complete pairs. The
count model used 2,000 variance-selected genes observed in at least 50 cells,
fixed NB concentration 10, and ten epochs per stage. The matched hierarchical
fit used the same rank and epoch schedule.

The count-model gene loss decreased from 0.634 to 0.615, and its fixed-state
EC loss decreased from 3.279 to 3.153. Joint refinement ended at gene loss
0.625 and EC loss 3.174. The hierarchical gene loss decreased from 6.966 to
6.914; its EC loss decreased from 3.288 to 3.221 with the state fixed and to
3.189 during joint fitting. Each full fit completed in under four minutes.

On all 27,383 labeled fitted pairs, the shared log-gene PCA initialization
gave 0.817 accuracy, 0.767 balanced accuracy, 0.685 macro F1, and -0.010
label silhouette. NB gene fitting improved these to 0.847, 0.814, 0.733, and
0.053. Count-model joint refinement reached 0.854, 0.824, 0.744, and 0.056.
The hierarchical gene stage reached 0.894, 0.893, 0.824, and 0.068, while its
joint stage reached 0.899, 0.902, 0.835, and 0.070.

For context, the earlier direct weighted EC factorization achieved only
0.487 accuracy and -0.110 silhouette, while the earlier binary fit achieved
0.688 and -0.239. The count-aware model therefore recovers much more
cell-type information than direct EC fitting, but the hierarchical
multinomial model is clearly better on this dataset. In particular, its
larger gains in balanced accuracy and macro F1 indicate better preservation
of the less frequent cell types rather than an improvement confined to
majority classes.

## 2026-07-28: hierarchical stage ablation

The microglialess hierarchical fit was ablated to determine whether separate
gene and fixed-state isoform stages are needed before joint optimization.
Every schedule started from the same rank-32 log-gene PCA state and used a
matched budget of 30 epochs. The original schedule used 10 gene, 10 isoform,
and 10 joint epochs. The alternatives used 10 gene plus 20 joint epochs, 10
isoform plus 20 joint epochs, or 30 direct joint epochs. Direct joint fitting
was tested at learning rates \(10^{-3}\), \(3\times10^{-3}\), and \(10^{-2}\).
No labels were used during fitting.

The original three-stage schedule ended at total loss 10.1029. Gene plus
joint fitting reached 10.1018, while isoform plus joint reached 10.1214.
Direct joint endpoints were 10.1152, 10.1044, and 10.1076 at increasing
learning rates. Thus fixed-state isoform pretraining is unnecessary for
optimizing the joint objective. Gene pretraining is useful relative to a
direct joint fit at the original \(10^{-3}\) refinement rate, but direct
joint optimization is nearly as effective when its learning rate is raised.

On the 27,383 labeled pairs, the original schedule gave 0.899 accuracy,
0.902 balanced accuracy, 0.835 macro F1, and 0.070 label silhouette. Gene
plus joint improved these to 0.902, 0.907, 0.842, and 0.072. Isoform plus
joint fell to 0.845, 0.812, 0.733, and 0.011. Direct joint fitting at
\(10^{-3}\), \(3\times10^{-3}\), and \(10^{-2}\) achieved accuracies 0.858,
0.896, and 0.911, respectively; the corresponding label silhouettes were
0.021, 0.062, and 0.082. At \(10^{-2}\), direct joint fitting also achieved
0.919 balanced accuracy and 0.861 macro F1.

The three-stage procedure is therefore not necessary. Within a label-blind
workflow, gene plus joint fitting is the conservative replacement because it
has the lowest final objective and slightly improves every reported label
metric without using labels for selection. The stronger post-hoc label score
of direct joint fitting at \(10^{-2}\) shows that gene pretraining is not
intrinsically required either; the earlier apparent requirement mostly came
from applying a refinement-scale learning rate from the initial PCA state.
That learning rate should not be selected from these label scores. A direct
joint default would need an unsupervised validation or convergence criterion.

## 2026-07-29: subisoform covariance and differential splicing

A reusable differential-splicing module now converts GTF exon chains into
local transcript paths between constitutive exonic regions. Transcripts with
the same path through a block are collapsed to one subisoform. Terminal
blocks represent alternative first and last regions, and genes without a
shared constitutive region receive one whole-gene block. Path usage is
represented in orthonormal log-ratio coordinates.

For each mouse-by-condition-by-cell-type pseudobulk, cell-level hierarchical
transcript estimates are decoded and UMI-weighted. The exact primer-specific
weighted EC designs are retained. Expected Fisher information for transcript
log abundance is calculated from EC responsibilities and propagated to path
log ratios. Null-space projection distinguishes estimable path contrasts from
the misleading zero variance that a bare pseudoinverse would return.

A second conditional estimate perturbs all transcripts assigned to the same
path while fixing the within-path transcript mixture and total spliced mass.
The paired multinomial likelihood is optimized directly. Two-path fits use a
bounded log-ratio grid for global initialization because different primer
normalizers can make the likelihood non-concave. Conditional and
transcript-profiled covariance are both written.

Downstream testing uses independent mouse pseudobulks. A scalar biological
variance is estimated by REML under the no-condition design, followed by a
multivariate GLS Wald test. Condition permutations reuse the null variance
and diagnose calibration. Simulations with known inferential covariance gave
2.7--5.0% type-I error at a nominal 5% over the biological-variance settings
tested.

The first real-data bootstrap exposed two limits of the local Gaussian
approximation. Information-null contrasts were already detected correctly,
but weak information and fitted paths near zero produced mode switching and
heavy-tailed log-ratio errors. Estimates are therefore retained for
inspection but enter GLS only when maximum log-ratio covariance eigenvalue is
at most \(1/8\) and every fitted path has at least 15% usage. In the two-path
smoke analysis, 1,257 fits all converged; 95.5% had identifiable conditional
covariance. Among bootstrap fits passing both reliability rules, mean squared
standardized error was 1.09 and component-wise 95% coverage was 94.0%.
Permuted-label asymptotic p-values were conservative but close to uniform:
3.3% were below 0.05 for conditional covariance.

Rare-path tests need a profile-likelihood or parametric-bootstrap procedure
rather than relaxed Gaussian thresholds. The current reliability filter
avoids claiming covariance-calibrated inference for those cases while
preserving their point estimates in the output.

The genome-wide microglialess run completed in 24 minutes with approximately
6.4 GB peak memory. Blocks with up to eight paths were fit. The gene-count
filter retained 2,820 blocks and 46,393 pseudobulk estimates, all of which
converged. Conditional covariance was identifiable for 73.5% of estimates
and 12.0% passed both reliability checks; profile covariance was identifiable
for 29.8% and 7.9% was reliable. This yielded 170 conditional and 113 profile
cell-type/block tests. Five conditional and three profile tests had nominal
\(p<0.05\), but no result passed 5% FDR.

The production parametric bootstrap comprised 3,000 refits. Mean squared
standardized error was 1.52 and component-wise 95% coverage was 94.6%.
Condition permutations were conservative: 2.9% of conditional and 2.5% of
profile asymptotic p-values fell below 0.05. The retained Gaussian tests are
therefore calibrated conservatively. Power is limited by mouse replication
and by the fraction of pseudobulks with enough path-specific information,
not by local optimizer convergence.

### Gene-coverage cutoff sweep

The genome-wide analysis was repeated with minimum pseudobulk gene counts of
25, 50, 200, and 500 UMIs and compared with the original 100-UMI run. The
number of fitted blocks was 10,100, 5,896, 2,820, 1,034, and 195 at increasing
cutoffs. The corresponding numbers of reliable differential tests were
1,335, 755, 283, 69, and 10. No cutoff produced a discovery at 5% FDR.

Calibration did not deteriorate at low coverage because the downstream
information and path-proportion filters remove estimates for which the
Gaussian approximation is unreliable. Conditional permutation rejection
rates at a nominal 5% were 2.7%, 3.1%, 2.9%, 2.8%, and 3.0%; profile rates
were 2.7%, 3.0%, 2.5%, 2.6%, and 5.0%. The final estimate is based on only
100 permutations. Parametric-bootstrap 95% coverage was 94.5%, 93.6%, 94.6%,
94.8%, and 94.6%, respectively.

For reliable two-condition, two-path conditional tests, mean asymptotic power
for an ILR effect of 0.5 was 0.49, 0.51, 0.44, 0.56, and 0.59. The last two
values are selected from only eight and two eligible tests. Among 103 test
identities shared by the 25- and 50-UMI runs, power was nearly unchanged
(0.453 versus 0.449). Lowering the cutoff therefore adds testable blocks
without reducing calibrated precision among shared tests. The 25-UMI cutoff
is preferred for this dataset when paired with the existing reliability
rules.

The initial 25-UMI run also exposed numerical underflow when an unconstrained
optimizer drove one path abundance to zero. Path-logratio perturbations are
now bounded to \([-20,20]\). The repeated run completed 359,711 of 359,713
local fits; the two remaining nonconvergences were excluded by the existing
fit-status check.

## 2026-07-30: paired cell-type differential splicing

The original differential table tested condition effects separately within
each cell type. It did not test the expected larger differences among cell
types. A second hypothesis family now tests each block across cell types
using the same reliable pseudobulk path estimates.

For a block, the cell-type model includes a fixed effect for each biological
mouse and treatment-coded cell-type effects. The mouse identifier includes
condition, so mouse fixed effects absorb condition and all other shared
mouse-level shifts. The omnibus Wald statistic tests all cell-type
coefficients. This is a paired comparison because multiple cell types from
the same mouse share its fixed effect. Cell-type labels are permuted within
mouse for calibration, preserving the observed missingness pattern.

At the 25-gene-UMI cutoff, 369 blocks had a full-rank conditional-covariance
test across at least two cell types; 218 also had a profile-covariance test.
The represented set spanned 12 cell types, and each included cell type
required at least three mice. Conditional covariance gave 57 nominal and 22
FDR-significant blocks. Profile covariance gave 28 nominal and five
FDR-significant blocks, all five among the conditional discoveries. The
union therefore contains 22 FDR-significant blocks with cell-type-dependent
path usage.

Within-mouse label permutations were conservative: 2.0% of conditional and
1.4% of profile asymptotic p-values were below 0.05. The absence of
condition-within-cell-type discoveries and the presence of 22 paired
cell-type discoveries are therefore consistent: they are different
hypotheses, and the larger cell-type effects are detected as expected.

## 2026-07-30: cell-resolved splice-path likelihood

To test whether pseudobulk EC aggregation was limiting power, I implemented a
cell-resolved local likelihood. Cells retain distinct hierarchical
transcript baselines and paired-primer EC count vectors, but all cells in one
mouse-by-condition-by-cell-type group share a single path-logratio
perturbation. The objective is the sum of cell-level multinomial likelihoods.
Analytic Fisher information is summed across cells and primers. The reported
group path response averages fitted cell path masses using decoded expected
gene molecules, and its covariance is obtained by propagating the shared
perturbation covariance. Downstream REML--GLS still operates on independent
mouse groups, so cells are not used as biological replicates.

The reusable implementation is in `tealeaf/sc/differential.py`. The dataset
runner supports block sharding, and a separate comparison script reconstructs
both methods on exactly matched reliable samples. A homogeneous-baseline test
verifies equality to aggregated pseudobulk fitting; a heterogeneous-baseline
test verifies recovery of an imposed shared shift.

The full assessment targeted the 527 blocks represented in the primary
25-gene-UMI condition analysis. Eight shards completed 41,242 local fits in
about 7.5 minutes wall time; every optimizer converged. The cell-resolved
method retained 21,624 reliable estimates and 836 condition tests, versus
21,010 estimates and 829 tests for pseudobulk conditional covariance on the
same block set. It produced no 5% FDR discoveries. Across 16,720
condition-label permutations, 2.49% of null p-values were below 0.05.

There were 807 downstream tests with exactly matched reliable samples. The
median cell-resolved/pseudobulk local covariance-trace ratio was 0.984 and
the median path-response RMS difference was 0.033 ILR units. P-value
Spearman correlation was 0.958, with 24 versus 25 nominal tests. For 259
matched two-condition, two-path tests, mean power for a 0.5 ILR effect was
0.5038 for pseudobulk and 0.5038 for cell resolution; the median
condition-coefficient SE ratio was 0.995.

The result is negative: retaining cells provides almost no usable precision
gain here. EC aggregation and the UMI-weighted decoder baseline already
preserve nearly all local information within a cell type, and estimated
between-mouse variation dominates the small inferential covariance
reduction. The pseudobulk likelihood remains the default. The cell-resolved
implementation is retained for datasets with stronger within-group
heterogeneity.

## 2026-07-30: compositional and joint condition tests

I evaluated two routes to greater differential-splicing power:
Dirichlet--multinomial regression for path composition and joint tests across
cell types or blocks. All methods preserve mouse-level condition assignment.

For each reliable path estimate, an effective multinomial count is obtained
by matching the theoretical multinomial ILR covariance to the EC-derived
conditional covariance in Frobenius norm, capped by observed gene UMIs.
Median effective count was 32 molecules, or 62% of gene UMIs; 0.6% reached
the cap. Fractional counts are fitted with a reference-logit
Dirichlet--multinomial regression. Concentration is estimated under the null
and fixed under the alternative. Fits at the infinite-concentration boundary
switch to an exact multinomial likelihood.

The annotation contains repeated genomic blocks that induce the same
transcript-to-path partition. These have identical statistics and cannot be
treated as independent hypotheses. Canonicalizing the unordered transcript
sets per path reduced 829 conditional block-by-cell-type rows to 505 unique
tests over 317 path partitions. Outputs retain all equivalent interval IDs.

Across the 505 matched tests, conditional GLS gave 12 nominal results and the
Dirichlet--multinomial gave 18; neither gave a 5% FDR discovery. The
Dirichlet--multinomial was conservative, with 3.76% of 10,100 permutation
p-values below 0.05. Pure multinomial regression was invalid: it gave 197
nominal and 155 FDR results, while 39.1% of permutation p-values were below
0.05. Between-mouse overdispersion is much larger than UMI sampling
variation.

ACAT combination across cell types produced 8 nominal GLS and 12 nominal
Dirichlet--multinomial partitions among 317 tests, with no FDR results.
Gene-level ACAT across distinct partitions produced 7 and 11 nominal genes
among 276, also with no FDR results.

I also fit a common condition effect across cell types with cell-type
intercepts. Naive Dirichlet--multinomial pooling produced 43 nominal and one
FDR partition, but synchronized mouse permutations showed 13.2% type-I
rejection. Its apparent discovery had asymptotic p-value 0.000064 and
permutation p-value 0.143. Naive shared GLS was also inflated at 8.4%.

The corrected joint Gaussian model includes a shared mouse random intercept
and a separate cell-type residual variance, both estimated by REML under the
null. It was conservative under synchronized mouse permutations: 2.19% of
7,169 null p-values were below 0.05. It produced 11 nominal partitions among
359 and no FDR discoveries; the smallest p-value was 0.0126 with FDR 0.86.
Gene-level ACAT produced 9 nominal genes among 313 and no discoveries.

The only calibrated sensitivity increase was therefore 12 to 18 nominal
within-cell-type results from the Dirichlet--multinomial. No tested strategy
produced an FDR-significant condition effect. The reusable effective-count,
multinomial, Dirichlet--multinomial, and clustered-GLS implementations are in
`tealeaf/sc/differential.py`; the dataset-independent comparison runner is
`extra_scripts/run_compositional_splicing.py`.

## 2026-07-30: exact compositional cell-type tests

The cell-type analysis was extended from Gaussian ILR Wald tests to joint
path-composition likelihood-ratio tests. Equivalent annotation intervals are
canonicalized to one transcript-to-path partition. For each partition,
EC-derived ILR covariance is converted to fractional effective path counts.
The null Dirichlet--multinomial model has mouse fixed effects; the alternative
adds cell-type effects. This tests all path log-ratios and represented cell
types jointly while absorbing condition and shared mouse effects.

The asymptotic reference was not adequate. Across 32,291 within-mouse
permutation fits in an initial sweep, 15.6% of asymptotic null p-values were
below 0.05. A pooled empirical null restored an overall 5.17% rejection rate
but was overly conservative because the null distribution varied strongly
with effective path count. The production analysis therefore used 500
within-mouse label permutations separately for every event and calculated
the finite-permutation rank p-value. Null parameters were reused, while every
permuted alternative was optimized. All 646 observed tests converged; valid
permutation counts ranged from 476 to 500.

A label-effect-independent minimum of 30 pseudobulk observations retained
406 tests for the primary FDR family. There were 149 exact p-values below
0.05 and 101 canonical splice-path partitions in 95 genes passed 5% BH FDR.
The count was 101 at minimum-observation thresholds 29, 30, 31, and 32, then
99 at 33. Median maximum cell-type logit effect among discoveries was 0.84;
median effective path count was 21.7. Thirteen of 14 canonical discoveries
from the earlier Gaussian analysis were recovered.

Pairwise Dirichlet--multinomial testing produced 97 significant
block-by-pair hypotheses but only 32 distinct blocks after correcting across
3,708 hypotheses. Permutation-calibrated multinomial testing produced 43
significant pair hypotheses over 15 blocks. A label-independent projection
onto the two molecule-weighted dominant paths was implemented by exact ILR
covariance transformation. It yielded 58 discoveries with a 5% path floor,
88 with a 1% floor, and 91 without a floor under four-bin cross-fitted
empirical-null calibration. None exceeded the full-path exact omnibus
result.

Reusable likelihood fitting, fitted-null reuse, convergence retries, and
composition projection are in `tealeaf/sc/differential.py`. The reusable
runner and shard merger are
`extra_scripts/run_celltype_compositional_splicing.py` and
`extra_scripts/merge_celltype_compositional_splicing.py`; the dataset wrapper
only supplies paths and run parameters.

## 2026-07-31: free alternative concentration

I tested re-estimating one global Dirichlet--multinomial concentration under
the cell-type alternative instead of fixing it at the null estimate. The
null and alternative still contain the same number of concentration nuisance
parameters, so the likelihood-ratio degrees of freedom continue to count
only cell-type composition coefficients. The implementation now evaluates
both finite Dirichlet--multinomial fits and the exact multinomial endpoint.
This was necessary because a finite cap at kappa (10^6) caused 25
alternatives at the boundary to be reported as nonconverged.

The corrected eight-shard run completed all 646 observed fits and all 323,000
permutation alternatives. With 500 within-mouse permutations, the 5% path
floor, and the 30-pseudobulk inference filter, the free-concentration model
found 101 partitions in 95 genes at 5% FDR. This is exactly the fixed-
concentration discovery set. In the matched, tolerance-aware reruns, exact
p-values had Spearman correlation 0.9995; 37 were lower, 556 equal, and 53
higher under the free model. Both models had 148 nominal events. Thus the
larger observed likelihood-ratio statistic does not imply uniformly greater
permutation-calibrated significance.

The null optimum was exact multinomial for 499 of 646 events and the
alternative optimum for 580. Of 147 events with finite null concentration,
81 moved to the multinomial endpoint under the alternative. For the 66 with
finite concentration under both models, median alternative/null kappa was
1.54. Cell-type mean effects therefore explain part of the residual
overdispersion, but allowing this change does not increase FDR discoveries.

This run exposed numerical ties in permutation statistics. Strict comparison
made observed statistics larger than mathematically tied permutations by
roughly (10^{-8}) to (10^{-7}), creating spuriously minimal ranks. Exact
ranks now count statistics within relative tolerance (10^{-6}) as ties.
Recomputing from saved event-specific null draws changed 207 p-values and
removed four apparent free-kappa discoveries; the final 101-event set was
unchanged from the fixed-kappa analysis. Future null tables store statistics
directly, and the shard merger recomputes tolerant ranks when null draws are
available.

## 2026-07-31: differential-splicing methods specification

I rewrote the differential-splicing methods document so the inferential
workflow is self-contained. The document now introduces the biological
questions and analysis sequence before the model, defines ECs, splice paths,
subisoforms, PSI, ILR coordinates, and the Helmert basis, and gives the
domain or dimension of every model variable when it is introduced.
Overloaded notation was removed, and matrix orientations were checked
explicitly.

Seven pseudocode algorithms now cover block construction, pseudobulk and
cell-resolved local likelihoods, Gaussian condition and paired cell-type
tests, fixed- and free-concentration compositional tests, pairwise and
dominant-path sensitivities, and ACAT or clustered shared-effect condition
tests. This was a documentation-only change: model implementations and
reported numerical results were not altered. The standalone document
compiles without warnings.

The notation was subsequently tightened to distinguish gene-level model
isoforms from block-specific splice paths. The local EC design contains all
retained isoforms assigned to the block's gene; this set has size \(T_g\)
and is shared across that gene's blocks. Only retained annotated spliced
isoforms receive nonzero entries in a block's isoform-to-path map. Other
columns, including modeled unspliced precursors, remain nuisance isoforms so
their EC probability is represented in the likelihood normalization.

The differential-splicing document now distinguishes the conditional path
covariance from the full-isoform nuisance covariance previously called only
the profile calculation. The latter is not an unconstrained refit. It
evaluates full isoform expected Fisher information at the constrained local
estimate, propagates its pseudoinverse through the path-ILR Jacobian, and
rejects responses that overlap an information-null isoform direction. The
document gives the Jacobian entrywise, the eigenspace thresholds, and the
identifiability test used by the implementation.

All method pseudocode was moved from prose enumerations into numbered
algorithm and algpseudocode environments. The full-isoform covariance is now
a separate algorithm, alongside block construction, local pseudobulk and
cell-resolved fits, Gaussian and compositional tests, sensitivity tests, and
joint condition tests.

## 2026-07-31: formal GLMM and Liang--Zeger GEE comparison

I implemented two formal mouse-clustered models for the joint condition test.
Both use cell-type-specific baselines and one condition effect shared across
cell types. The Liang--Zeger GEE fits a reference-category multinomial mean,
an exchangeable within-mouse working correlation, Pearson scale and
correlation updates, and the canonical cluster sandwich covariance. The
GLMM uses a full-path multinomial logit with an isotropic Gaussian
mouse-level random intercept in the non-reference path coordinates. It fits
the marginal likelihood by a Laplace approximation, with posterior modes
found by damped Newton updates and the outer objective differentiated by JAX.

The covariance-matched effective path counts are fractional. GEE can use
these directly as a quasi-score. A formal multinomial likelihood cannot, so
the GLMM integerizes each row by largest remainder while preserving its
rounded effective total. Shape-specialized JAX functions are cached, fitted
nulls are reused across condition permutations, and GEE permutations warm
start from the observed coefficients.

The complete comparison contained 359 canonical path partitions. Every
method used the same 50 synchronized mouse-level condition permutations.
For 88 two-condition partitions, 50 semi-synthetic replicates at ILR effects
0.25 and 0.5 yielded 4,369 full-rank replicate designs per effect. There were
no fit exceptions. Clustered GLS, Dirichlet--multinomial, and GLMM converged
for all 359 observed tests and all 17,917 null fits. GEE converged for 355
observed tests and 17,702 null fits.

The asymptotic references were poorly calibrated. Permuted-null rejection at
0.05 was 2.55% for clustered GLS, 13.79% for Dirichlet--multinomial, 41.49%
for GEE, and 10.76% for GLMM. GEE's 116 apparent asymptotic FDR discoveries
are therefore artifacts of severe small-cluster sandwich inflation. The
Dirichlet--multinomial had one asymptotic FDR call; GLS and GLMM had none.
After event-matched finite-permutation calibration, no method had a 5% FDR
discovery. Nominal permutation counts were 17, 18, 8, and 19, respectively.

Event-matched power at ILR effects 0.25 and 0.5 was 17.1% and 32.2% for
clustered GLS, 18.8% and 36.4% for Dirichlet--multinomial, 13.9% and 26.8%
for GEE, and 19.9% and 39.0% for GLMM. The GLMM was best, improving over
clustered GLS by 2.8 and 6.9 percentage points, while canonical GEE was worse
after calibration. The gain is modest and does not remove the need for
event-specific permutation inference.

Reusable GEE, GLMM, integerization, and stable permutation-rank code is in
`tealeaf/sc/differential.py`. The dataset-independent comparison and merger
are `extra_scripts/run_clustered_compositional_models.py` and
`extra_scripts/merge_clustered_compositional_models.py`; the analysis wrapper
only supplies inputs and execution parameters.

## 2026-08-01: GLMMs on raw equivalence-class counts

I added a paired-primer GLMM whose observations are the integer EC counts,
not effective splice-path counts. For each gene, fixed condition or cell-type
effects and mouse random intercepts act on reference-isoform logits. The EC
probability is the compatibility-weighted isoform mass divided by the sum of
that mass over gene-assigned ECs. Oligo(dT) and random-hexamer observations
share biological coefficients and mouse effects while retaining separate
weighted compatibility matrices and normalizers. Both multinomial and direct
EC Dirichlet--multinomial likelihoods are implemented.

Three inference families are compared. Laplace inference differentiates
through damped Newton posterior-mode updates. Its Hessian calculation was
changed from one dense mouse-by-isoform matrix to mouse-specific blocks built
by accumulating per-observation Hessians; conditional independence makes
these formulations identical and removes cubic scaling in the total number
of mice. The multinomial variational method uses diagonal Gaussian mouse
effects, analytic latent EC-to-isoform responsibilities, and the tilted
softmax expectation bound from Knowles and Minka (2011). The cited algorithm
is NCVMP, not conjugate CAVI; the implementation alternates optimized local
coordinates and bounded global updates. Exact-likelihood Monte Carlo VI fits
either the ELBO or the Li and Turner variational Renyi bound at alpha 0.5.
The Renyi fit uses 128 fixed antithetic draws and warm starts from the ELBO
fit. An initial global importance ratio over all mice had ESS 1/128 and drove
the random-effect variance to its lower bound. Because the posterior
factorizes by mouse conditional on the global parameters, the corrected
objective sums mouse-specific Renyi log-averages. A simulation then gave
median ESS above 127/128 and converged for both likelihood families. The
direct EC Dirichlet--multinomial does not admit the tilted coordinate
bound because its gamma functions contain nonlinear EC probabilities, so it
uses only Laplace and Monte Carlo VI.

The reusable implementation is `tealeaf/sc/ec_glmm.py`. Ambiguous-EC
simulation tests verify the tilted expectation bound, coefficient recovery,
and Dirichlet--multinomial support. The genome-wide runner keeps genes with
two to six EC-supported isoforms, at most 128 gene-assigned ECs, and mean
gene depth of at least 100 UMIs per retained pseudobulk. It precomputes EC
totals for eligibility scanning and uses a locked 160 MB preparation cache;
the initial implementation's repeated full sparse-matrix slice per gene was
the dominant startup cost. The real-data and simulation comparisons are
running. Final convergence, calibration, power, and timing results will be
added when those jobs complete.

The simulation comparison used 96 null and 96 effect-0.5 replicates per
likelihood. Empirical 95th-percentile null thresholds gave 5.21% null
rejection for every method. Multinomial Laplace, tilted ELBO, and Renyi VI
all had 54.17% calibrated power and condition-coefficient RMSE 0.223.
Dirichlet--multinomial Laplace, Monte Carlo ELBO, and Renyi VI all had 41.67%
power and RMSE 0.236. Median nested-fit times at the nonzero effect were
9.87, 1.87, and 1.58 seconds for the multinomial methods and 21.41, 3.06,
and 3.35 seconds for the Dirichlet--multinomial methods. Laplace converged in
every replicate. VI convergence was 97--100%, and median mouse-level ESS was
127.5--127.7 of 128. Thus the approximation choice did not affect calibrated
power or coefficient error in this setting, while VI was five- to seven-fold
faster. The genome-wide real-data arrays cover 65 of 66 depth-eligible genes;
the omitted gene has 3,905 ECs, while all retained genes have at most six
supported isoforms. Both condition and cell-type arrays completed for all 65
genes without fit exceptions.

The first real-data pass initialized each approximation independently. Its
pooled multinomial Laplace-versus-Renyi evidence-gain correlation was only
0.41. A second, authoritative pass initialized diagonal Gaussian VI from the
Laplace posterior means and variance diagonal and initialized Renyi VI from
the ELBO fit; relative objective tolerance was relaxed to 1e-8 and the cap
raised to 500 iterations. This raised convergence to 92--98% for condition
and 94--100% for cell type, with no fit exceptions. It did not remove the
dimension-dependent discrepancy: Laplace-versus-Renyi gain correlation was
0.974 (condition) and 0.9996 (cell type) for 38 two-isoform genes, but 0.49
and 0.59 for 20 three-isoform genes. This points to the diagonal variational
covariance, which cannot represent within-mouse posterior correlation among
isoform logits. A full-covariance Gaussian extension is needed before using
mean-field VI for multi-isoform evidence comparisons.

For two-isoform genes, the exact-likelihood ELBO and alpha-0.5 Renyi gains
evaluated at the same tilted-bound posterior were almost identical: Spearman
correlation 0.998 for condition and 1.0 for cell type. Median condition gains
were 5.54 for the tilted training bound, 4.95 for exact Monte Carlo ELBO,
4.93 for Renyi evaluation, 4.93 for Laplace, and 4.95 for separately
optimized Renyi VI. Cell-type medians were 89.06, 74.38, 74.36, 72.12, and
74.36. The training bound is differentially looser under the alternative and
inflates nested gains, so it must not be treated as a likelihood-ratio
statistic. One-dimensional simulation thresholds apply only to two-isoform
condition tests, not multi-isoform or cell-type omnibus tests.

## 2026-08-01: full-covariance EC GLMM and logistic-normal overdispersion

I implemented a mouse-factorized full-covariance Gaussian posterior for the
EC-count multinomial GLMM. Without extra observation noise, each mouse block
contains all nonreference isoform-logit random intercepts. With observation
noise, it jointly contains the shared mouse intercept and one isotropic
logistic-normal residual per observed pseudobulk. The generative covariance
remains isotropic for both effect classes; only the posterior Cholesky is
unrestricted. A correlated version of the Knowles--Minka tilted log-sum-exp
bound uses the full row covariance, and the Gaussian KL is analytic. Exact
ELBO and alpha-0.5 Renyi objectives use mouse-specific antithetic importance
weights. Full VI warm-starts from diagonal VI by copying means and marginal
variances and initializing off-diagonal covariance to zero.

The ordinary Laplace implementation now exposes its full posterior covariance
and supports the same observation-level residual model through a general
latent-block selection matrix. Cholesky warm starts symmetrize covariance and
add the minimum eigenvalue-based jitter needed for positive definiteness. A
fixed-draw Monte Carlo objective can overfit the much larger residual-model
Cholesky. The authoritative residual fit therefore uses the analytic tilted
objective and evaluates exact ELBO and Renyi bounds with fresh draws. The
reusable full-posterior implementation is `tealeaf/sc/ec_glmm_full.py`; the
dataset runner and simulation benchmark expose all new variants.

In 96 null and 96 effect simulations with three isoforms, full covariance did
not materially change multinomial VI. Diagonal and full tilted gains had
Spearman correlation 0.9999 and median absolute difference 0.073; diagonal
and full Renyi gains had correlation 0.99996 and median difference 0.044.
Every VI variant had calibrated power 43.75% and coefficient RMSE 0.237--0.238.
The full tilted fit took 2.43 seconds per nested effect fit versus 2.00 seconds
for diagonal tilted VI.

The direct EC Dirichlet--multinomial gave the same result. In 96 three-isoform
null and effect simulations, diagonal versus full ELBO gains had correlation
0.9997 and median absolute difference 0.150; diagonal versus full Renyi gains
had correlation 0.9999 and median difference 0.092. Both posterior families
had calibrated power 21.88%. Across 65 real-data genes, diagonal versus full
ELBO gains had correlations 0.998 for condition and 0.9998 for cell type,
with median absolute differences 0.0035 and 0.0050. Full DM ELBO fits
converged for every gene. Fixed-draw full Renyi fits were less stable, as in
the multinomial analysis.

The genome-wide result was consistent. Fresh-sample condition ELBO gains from
diagonal and full tilted posteriors had correlation 0.991 and median absolute
difference 0.002 across 65 genes. For the 20 three-isoform genes these were
0.938 and 0.024; cell-type gains had correlation 1.000 and median difference
0.023 in that subset. The earlier low Laplace-versus-VI correlation therefore
was not caused by omitted variational posterior covariance. It more likely
reflects the Laplace outer objective or optimization; at least one
three-isoform real-data Laplace fit remains at its all-zero outer start while
reporting convergence.

I also simulated logistic-normal overdispersion with true residual standard
deviation 0.6. Tilted full-covariance VI recovered median 0.559 and had 35.42%
calibrated power, close to direct Dirichlet--multinomial Laplace at 36.46% and
above the misspecified ordinary multinomial at 31.25%. The high-dimensional
logistic-normal Laplace approximation had only 14.58% power. In real data,
the tilted residual model converged under both nested models for 80.0% of
condition genes and 86.2% of cell-type genes. Median residual standard
deviations were 0.040 and 0.029. It changed condition ELBO gains modestly
(median absolute difference 0.335) but reduced the much larger cell-type
gains by a median absolute 10.68. Given weaker convergence and many residual
scales near zero, it is an informative overdispersion sensitivity analysis,
not a replacement for the primary model.

## 2026-08-01: block-specific EC GLMM tests

The original EC GLMM runner tested unrestricted gene-wide isoform effects. It
therefore could not report splice-block events, even though it retained the EC
counts. I added a block-specific nested model that continues to evaluate every
EC and represented isoform in the gene likelihood. Under the null, condition
has unrestricted nonreference-isoform coefficients. Under the alternative,
cell-type coefficients are added only in the Helmert contrast space of the
tested block's paths. The null also includes cell-type coefficients spanning
every isoform-logit direction orthogonal to those block effects; without
these terms, unrelated isoform changes elsewhere in the gene would inflate
the block statistic. The alternative combines the nuisance and block bases
and is unrestricted over gene isoform logits. Isoforms outside the block
remain in the likelihood normalizer as nuisance isoforms. Mouse random
intercepts remain unrestricted. For `J` cell types and `S_b` represented
paths, the tested dimension is `(J - 1)(S_b - 1)`.

Cell-type label permutation is invalid for this conditional block null: it
would destroy the real orthogonal cell-type effects that the null explicitly
permits. Calibration instead uses a parametric bootstrap from each fitted
block null. It draws mouse and optional observation effects, preserves every
observed primer-row total, generates EC counts under the selected likelihood,
and refits both nested models. Per-event bootstrap counts are kept modest,
and each observed test is calibrated against a leave-event-out pool matched
by likelihood, tested dimension, and gene-depth quantile. Rare tested
dimensions are pooled across depth within likelihood. The merger also
calibrates every bootstrap
statistic against the corresponding leave-event-out null pool. It withholds
FDR values for a likelihood unless the aggregate empirical type-I error at
the 0.05 threshold lies between 0.025 and 0.075.

There are 406 canonical coverage-eligible blocks in the primary cell-type
table. Of these, 399 have at least two EC-represented paths and at most ten
supported gene isoforms and can be tested. The other seven have only one
represented path in the EC design, so their block effect is not identifiable.
A real-data smoke fit of the corrected multinomial model converged for both
nested models. Its observed statistic was 0.00035 with 12 degrees of freedom,
and one parametric-bootstrap statistic was 0.000034. This contrasts with the
much larger statistic from the model that omitted orthogonal cell-type nuisance
effects and confirms that the omission could create false block attribution.
The full-covariance Dirichlet--multinomial objective is memory-sensitive to
its fixed Monte Carlo sample count; the analysis uses 16 synchronized
antithetic samples and validates the approximation by empirical null
calibration.

The final pooled calibration used five parametric-bootstrap replicates per
block, giving 1,830--1,960 converged null statistics per likelihood after
failed refits were excluded. The ordinary multinomial converged for 395 of
399 observed blocks, had empirical null rejection 0.0480, and found 111
nominal and 73 FDR-significant blocks in 66 genes. The logistic-normal model
converged for 362 blocks, had null rejection 0.0486, and found 104 nominal and
70 FDR-significant blocks in 63 genes. The direct Dirichlet--multinomial
converged for 389 blocks, had null rejection 0.0494, and found 70 nominal and
23 FDR-significant blocks in 21 genes. All methods passed the prespecified
calibration gate. Multinomial and logistic-normal results shared 66 blocks;
16 blocks were significant under all three EC likelihoods. Their overlaps
with the 101-event primary path-count analysis were 29, 29, and seven blocks,
respectively.

## 2026-08-02: direct EC-count block eligibility

The 399-block EC analysis inherited its test universe from the compositional
path-count analysis. That was too restrictive for an EC likelihood: its
upstream filters required identifiable per-pseudobulk path estimates, every
retained path proportion to be at least 0.05, and at least 30 such estimates.
Only the requirement that at least two block paths occur in the EC design is
intrinsic to the EC GLMM.

The EC runner now screens coverage directly from the primer-specific EC rows
that enter each gene likelihood. For each gene it retains pseudobulks meeting
a gene-UMI cutoff, retains cell types represented in at least three mice, and
requires at least two such cell types. The fixed-effect design and degrees of
freedom are then built from these gene-specific rows and cell types. An event
table can still be supplied as an optional restriction, but no path-estimate
table is required.

An initial diagnostic that summed every gene-associated EC column suggested
2,205 blocks at 25 gene UMIs plus 30 covered pseudobulks and 6,022 blocks at
10 gene UMIs without the sample-count requirement. An end-to-end smoke test
showed that this total included primer-specific EC rows with no compatibility
to the retained gene isoforms; those rows are absent from the likelihood.
Recomputing coverage from exactly the modeled EC rows gives 1,926 and 5,650
blocks for the same two cutoff rules.

All six likelihood-by-coverage runs completed. Five null replicates per block
were pooled by likelihood and tested dimension using leave-block-out ranks.
Splitting each dimension by depth quartile made the minimum attainable
p-values too coarse for BH over thousands of hypotheses. Pooling by dimension
alone restored tail resolution, while an independent calibration audit in
four depth quartiles showed no count-dependent inflation. Aggregate null
rejection rates were 0.0491--0.0498 for the strict analyses and 0.0497 for all
three permissive analyses. Across depth quartiles, rates ranged from 0.0439 to
0.0582 for strict tests and 0.0428 to 0.0541 for permissive tests. Every method
passed the prespecified 0.025--0.075 calibration interval both overall and in
each depth quartile.

For the strict universe, multinomial, logistic-normal multinomial, and direct
Dirichlet--multinomial convergence rates were 1.000, 0.991, and 0.988. They
gave 415, 377, and 145 nominal discoveries and 174, 157, and zero BH
discoveries. For the permissive universe, convergence was 0.999, 0.984, and
0.989; nominal discoveries were 997, 974, and 421; and BH discoveries were
294, 193, and zero. MN and LN shared 152 strict discoveries and 192
permissive discoveries.

The permissive cutoff therefore improves absolute empirical discovery yield
for MN and LN, but less than proportionally to the 2.9-fold larger hypothesis
family. On the 1,926 shared blocks, refitting with the permissive rows and
applying BH only within that common family gives 316 MN, 265 LN, and 40 DM
discoveries, compared with 174, 157, and zero from the strict fits. Lower-count
rows do add signal; the genome-wide BH burden spends much of that gain. For MN,
90 blocks are significant under both genome-wide analyses, 84 only under the
strict analysis, 94 only under the permissive fit among shared blocks, and 110
only in the expanded universe. The analogous LN counts are 62, 95, 53, and
78. Statistic rank correlations between strict and permissive fits on shared
blocks are 0.781, 0.764, and 0.760 for MN, LN, and DM.

## 2026-08-02: one block EC GLMM for cell-type and condition DS

The block EC runner now supports both biological contrasts through the same
nested logistic-normal multinomial GLMM. The cell-type test uses condition as
a nuisance fixed effect and tests block-path by cell-type interactions. The
condition test is run separately within each cell type, uses an intercept-only
nuisance design, and tests block-path by condition interactions. Both retain
mouse-level random effects, paired primer-specific EC likelihoods, and pooled
leave-test-out parametric-bootstrap calibration. A test identifier combines
the block and cell type for condition tests so calibration and multiplicity
operate on block-by-cell-type hypotheses rather than collapsing them by block.

For condition DS, each retained condition must have EC-covered observations
from at least three mice. At a cutoff of 10 modeled gene UMIs this gives
17,685 block-by-cell-type tests over 6,550 distinct blocks. Screening is now
written to a compact candidate manifest and reused by every shard; this avoids
repeating genome-wide annotation and EC-support screening in each array task.
A production-criteria smoke fit converged under both nested models in 132 and
179 iterations. The full logistic-normal condition analysis will use five
parametric-bootstrap replicates per test and the same aggregate, depth, and
sample-count calibration audits as the cell-type analysis.

The completed condition run initially contained 17,685 annotated
block-by-cell-type tests. Inspection of tied top statistics showed that
distinct annotated blocks can induce the same isoform partition after
unsupported isoforms are removed. Collapsing tests by gene, retained isoforms,
canonical path partition, covered rows, and contrast levels leaves 14,833
nonredundant tests over 5,608 representative blocks and 4,080 genes. Of 14,827
converged and calibratable tests, 1,404 have nominal p-values below 0.05, but
none pass BH at 0.05; the minimum FDR is 0.163. Empirical null rejection is
0.0500 overall, 0.0462--0.0611 across depth quartiles, and 0.0475--0.0555
across sample-count quartiles.

Applying the identical supported-partition collapse to the completed
cell-type logistic-normal fits reduces 5,650 annotated blocks to 4,825 tests.
There are 800 nominal results among 4,744 converged and calibratable tests,
but no BH discoveries under raw empirical ranks; the minimum FDR is 0.0531.
This apparent loss was a tail-resolution artifact. After equivalent tests
were removed, each degrees-of-freedom pool had a minimum p-value near
$4.5\times10^{-4}$, while BH required values near $10^{-5}$ at the first
ranks.

The corrected calibration cross-fits a generalized Pareto
peaks-over-threshold model within each degrees-of-freedom stratum, holding out
whole block tests. It then maps both observed and held-out-null tail scores
through a global leave-block-out empirical null CDF. The second stage removes
small GPD tail bias while retaining pooled-null resolution. At a 90% tail
threshold, held-out null rejection for the cell-type logistic-normal test is
0.04994, 0.00999, and 0.000977 at p-value thresholds 0.05, 0.01, and 0.001.
The multinomial values are 0.04999, 0.01000, and 0.000982, and the
Dirichlet--multinomial values are 0.05001, 0.00997, and 0.001002.

After nonredundant BH correction, multinomial, logistic-normal multinomial,
and Dirichlet--multinomial identify 272, 243, and 59 cell-type DS partitions.
The logistic-normal result is stable across 80%, 85%, 90%, and 95% GPD
thresholds, which give 243, 242, 243, and 239 calls. Their intersection has
237 partitions; 229 are also multinomial calls, and 44 are called by all
three EC likelihoods.

For condition within cell type, the 90% threshold remains nonsignificant
with minimum FDR 0.0637. Thresholds of 80%, 85%, 90%, and 95% give 27, 27,
zero, and zero calls, with an empty four-way intersection. Thus no condition
partition is robust to the tail threshold. The primary DS procedure is one
nonredundant logistic-normal block EC GLMM parameterized for either contrast:
it gives strong cell-type DS but no robust condition DS in this dataset.

## 2026-08-03: differential-splicing manuscript reorganization

The differential-splicing manuscript now presents the implemented
multinomial and logistic-normal EC GLMMs before every supporting method. The
main text gives the paired-primer EC probability model, mouse and pseudobulk
random effects, orthogonal block-path and nuisance tensor construction for
both biological contrasts, full-covariance tilted variational fit,
null-to-alternative warm start, parametric bootstrap, cross-fitted GPD
scoring, global leave-block-out recalibration, supported-partition
deduplication, and separate BH families. It includes pseudocode for fitting
and calibration and then reports the primary cell-type and condition results.

The local path-estimate covariance model, Gaussian and compositional tests,
direct EC Dirichlet--multinomial, GEE, Laplace, CAVI, Monte Carlo and
R\'enyi variants, power studies, implementation details, and supporting
results now follow the appendix marker. The document title and notation table
were also revised to match the primary observed-EC analysis.

## 2026-08-03: matched junction-method benchmark on two datasets

Implemented a reusable STARsolo junction benchmark for the microglialess and
GSE233208 datasets. Split-seq barcode extraction uses the three eight-base
segments from read 2 and the ten-base UMI. STAR writes complex cell barcodes
with segment separators; the reader and BAM tagger now remove those separators
before matching reference cell IDs. A 200,000-read GSE233208 smoke run retained
3,576 junction UMIs across 2,717 cells and recovered all 18,040 labeled cell
barcodes from that run.

The shared representation is a sparse subject-by-cell-type junction UMI
matrix plus junction and sample metadata. Spliced BAM evidence is restricted
to unique alignments and deduplicated by cell, UMI, reference, start, and
CIGAR before reads are tagged by pseudobulk. This supplies biological
replicates rather than treating cells as replicates. Reusable adapters export
LeafCutter junction files, annotate scQuint intron groups once from the GTF,
and create pairwise condition-within-cell-type and subject-matched cell-type
contrasts. LeafCutter uses subject as a paired-test covariate. scQuint lacks a
covariate interface, so it receives the same pseudobulks but cannot use the
pairing.

MAJIQ academic v3.0.23 was built against the available HTSlib and Zstandard
libraries. An annotation smoke test found that a VM-built extension required
a newer glibc than the compute nodes provide, so a reproducible Slurm installer
now builds the environment on a compute node. The rebuilt executable applied
the academic license and parsed 58,870 human and 56,953 mouse genes in
annotation smoke tests. Its workflow
converts the annotation to a splice graph, extracts one SJ object per
pseudobulk BAM, builds a common splice graph, creates a multi-sample
PSI-coverage object, and runs Heterogen for every matched contrast. The common
comparison uses the 95th percentile of the posterior-bootstrap Mann--Whitney
p-value, which is more conservative than the raw plug-in p-value. MAJIQ also
lacks a paired-test design, a limitation recorded in the manuscript.

Both full alignment and comparator dependency chains have been submitted.
The GSE233208 primary multinomial and logistic-normal block EC GLMMs are also
scheduled with the same cell-type and within-cell-type condition definitions.
Results will be merged into a common per-feature schema with within-contrast
BH values and Simes cell-type omnibus scores after all chains complete.

The first full mouse alignment array completed STAR successfully for every
run but failed in the post-alignment tagger because the batch scripts had been
submitted before the package-path fix. The aligned BAMs and STARsolo junction
matrices were intact. Added a reusable tag-only array and replaced the stale
mouse aggregation and comparator dependencies, avoiding any realignment.

## 2026-08-03: junction benchmark recovery

GSE233208 reused numeric case identifiers across diagnoses. Grouping only by
the published case number therefore combined cells from different biological
conditions in 84 pseudobulk samples. Reference-label preparation now supports
an optional group-prefix column, and this dataset uses diagnosis-qualified
case identifiers. Regeneration produced 792,554 labeled cells from 40 runs,
51 subjects, and no subject assigned to more than one condition. The completed
STAR alignments are being retagged and reaggregated; no realignment is needed.

Two sparse-contrast failure modes were found in the mouse comparator run.
LeafCutter can write a valid cluster-significance table with no finite
p-values and then fail during its optional effect-size calculation. The
adapter now accepts that table and records zero testable clusters. scQuint can
abort a full contrast when one intron-group optimizer does not improve on its
null fit. The adapter now retains the other groups and assigns the failed
group a conservative p-value of one with an explicit failure flag. Focused
adapter and reference-label tests pass, and the affected LeafCutter and
scQuint arrays have been restarted while the independent MAJIQ run continues.

## 2026-08-04: terminal comparator fixes

The mouse junction and MAJIQ build stages completed for all 936 pseudobulk
samples. MAJIQ Heterogen then completed each statistical calculation but
failed while writing its table because v3.0.23 unconditionally dropped an
optional `mannwhitneyu-stats_extra` column that was not present. A small
external patch makes that drop tolerant of an absent optional column, and the
cluster installer now applies the patch reproducibly. An end-to-end contrast
then completed and produced a normalized result, so the 490 mouse contrasts
were resubmitted without rebuilding the splice graph or PSI coverage.

Of the 490 repaired mouse scQuint jobs, 480 completed. Nine sparse contrasts
failed after filtering removed every intron, and one failed when process-based
joblib workers could not reconstruct their semaphores. The adapter now allows
an empty intermediate intron matrix to reach scQuint's existing empty-result
return and uses thread-based group parallelism with one Torch thread per fit.
The ten affected contrasts were resubmitted; the sparse cases now complete.

GSE233208 retagging, junction aggregation, and pseudobulk BAM splitting all
completed. The first comparator submission processes exceeded Slurm's default
memory while loading the aggregate junction bundle. They were rerun with an
explicit allocation, producing 434 replicate-supported contrasts. The MAJIQ
submission helper also now initializes Python before reading an already
existing contrast plan. LeafCutter and scQuint are running, and the MAJIQ
annotation, 1,449-sample junction extraction, build, tests, and common summary
are linked in a new dependency chain.

The completed mouse summary initially failed because deliberately empty
no-test outputs contain a newline but no columns. The common summarizer now
skips both zero-byte and no-column tables before concatenation. This preserves
zero testable features for those contrasts without treating them as failed
analyses.
