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

## 2026-08-04: GSE233208 EC block identifier mismatch

The first GSE233208 EC block arrays completed all 192 shards but performed no
fits. The quantification and hierarchical feature files use unversioned
Ensembl gene and transcript identifiers, whereas the GENCODE block cache uses
versioned identifiers. Exact matching therefore rejected every block. Gene
and transcript matching now removes only the terminal Ensembl version suffix;
block IDs remain unchanged. Corrected screening retained 11,597 nonredundant
cell-type block tests and 73,144 condition-within-cell-type tests across 9,191
blocks. The runner now fails early when the genome-wide
screen returns zero candidates, and the merger skips legitimate empty shards
but rejects a collection with no fitted results.

The initial corrected fit layout used 96 shards per contrast family. After
about eight hours, array-wide progress logs projected median shard runtimes of
62 hours for cell type and 38 hours for condition, exceeding the submitted
16-hour limit. Running-job limits could not be extended. The arrays were
stopped and resubmitted as 768 shards per family, with at most 96 concurrent
tasks and limits of three days for cell type and two days for condition. This
reuses both candidate manifests and the prepared EC cache while reducing each
cell-type shard from about 121 to 15 tests and each condition shard from about
762 to 95 tests.
## 2026-08-05: nonbootstrap EC-block inference

The parametric bootstrap requires one observed null/alternative pair plus a
null/alternative pair for every simulated replicate. I implemented a
nonbootstrap sensitivity analysis using the actual Laplace marginal
likelihood rather than treating an ELBO difference as a likelihood ratio.
The Laplace fitter now accepts the block-specific fixed-effect tensor used by
the primary GLMM. The alternative tensor begins with every null column and
appends only the block-path interaction columns, so the tested degrees of
freedom are the number of appended coefficients. The alternative is
warm-started by copying the null parameters and setting those coefficients to
zero.

For maximized negative Laplace objectives \(L_0\) and \(L_1\), the runner
reports \(D=2\max(0,L_0-L_1)\), the nominal chi-square p-value with the tested
dimension, and the BIC approximation
\(\log \widehat{\mathrm{BF}}_{10}=(D-q\log M)/2\), where \(M\) is the number
of independent subjects. The latter is an approximation, not a
prior-integrated Bayes factor. Both results require two fits instead of two
fits per bootstrap replicate. A merger deduplicates EC-equivalent block
partitions and applies BH within method.

I also tested a nuisance-adjusted score calculation based on automatic
differentiation of the Laplace objective. A full observed Hessian used more
than a gigabyte for a two-isoform unit test, and differentiating clusterwise
objective contributions had similar scaling because both differentiate
through the iterative random-effect mode. That route was removed rather than
used genome-wide. Hessian-vector or analytic implicit-differentiation methods
remain possible future work, but they are not needed by the LRT/BIC path.

Twenty-two focused EC-GLMM tests pass. One real mouse block required about 49
seconds for null and alternative fits under both multinomial and
logistic-normal Laplace models. Its six-degree-of-freedom LRT p-value was
0.708 under both models, compared with completed bootstrap/GPD p-values 0.594
and 0.536. Its observation-noise estimate reached the numerical lower bound,
so the full mouse comparison must assess both calibration agreement and the
frequency of variance-boundary fits before this can replace bootstrap
calibration. That validation is running over the same direct-coverage
candidate family used by the completed mouse benchmark.

The first unscaled validation batch exposed line-search failures in about 5\%
of fits. These stopped after only a few iterations with gradients in the
thousands. The Laplace objective is a sum over all observed molecules, and
its magnitude varied enough across genes to make one set of L-BFGS stopping
and line-search scales unreliable. The optimizer now divides both objective
and gradient by the absolute initial objective while retaining the unscaled
objective for the LRT. This positive constant leaves every optimum and
likelihood ratio unchanged. A previously failed block then converged in 37
null and 17 alternative iterations for the multinomial model, and in 36 and
34 iterations for the logistic-normal model. The scaled full validation was
restarted; the partial unscaled results are not used for final inference.

As a separate calibration diagnostic, applying the nominal chi-square
reference to the 26,373 converged multinomial statistics already simulated
under the variational-fit null gave rejection rates 0.0389, 0.00781, and
0.00091 at levels 0.05, 0.01, and 0.001. The corresponding 21,955
logistic-normal null statistics gave 0.0374, 0.00824, and 0.00087. Thus the
chi-square reference is mildly conservative for those simulated null
statistics rather than anti-conservative. This is informative but not a
substitute for comparing the completed Laplace results with the empirical
bootstrap analysis, because the fitting objectives differ.

The final strict mouse run completed all 4,825 nonredundant cell-type tests.
The multinomial Laplace null and alternative both met the scaled-gradient and
mode-score criteria for 94.8% of tests; logistic-normal convergence was 93.6%.
Among tests converged under both Laplace and the completed bootstrap fit, the
Laplace and variational statistics had Spearman correlations 0.982 for
multinomial and 0.964 for logistic-normal. P-value rank correlations were
0.950 and 0.922. The multinomial LRT gave 703 nominal calls and 309 BH calls
among jointly converged tests, compared with 753 and 261 for bootstrap/GPD;
243 BH calls overlapped. Logistic-normal gave 633 nominal and 266 BH calls,
compared with 707 and 230; 204 BH calls overlapped. The differing nominal and
BH ordering reflects the more extreme analytic chi-square tail despite fewer
p-values below 0.05.

The BIC threshold \(\widehat{\mathrm{BF}}_{10}>10\) selected 155 multinomial
and 128 logistic-normal tests in the jointly converged comparison. Of these,
149 and 122 were also LRT BH calls, while 140 and 116 were bootstrap BH calls.
This is a conservative high-confidence subset rather than an MHT-corrected
decision rule. In the logistic-normal model, 758 of 4,438 jointly converged
alternative fits put the observation-noise SD at its lower bound. The
ordinary multinomial LRT is therefore the cleaner nonbootstrap sensitivity;
the logistic-normal Wilks approximation is less regular and remains a
secondary result. Both likelihoods were retained for the independent human
cell-type analysis; its larger within-cell-type condition analysis uses the
recommended multinomial model.

Before the human run, the logistic-normal nested fits were changed to
warm-start from the corresponding multinomial null and alternative fits by
appending only the observation-noise scale. On a mouse smoke block this
reduced logistic-normal outer iterations from 38 and 34 to 14 and 18, while
changing the LRT statistic from 3.7691 to 3.7701. This warm start is used only
when both Laplace methods are requested in multinomial-first order.

The human analysis completed 11,597 cell-type tests under both likelihoods
and 73,144 within-cell-type condition tests under the recommended multinomial
model. Cell-type convergence was 0.949 for multinomial and 0.943 for
logistic-normal. The methods gave 2,580 and 2,161 nominal calls, 1,487 and
1,125 BH calls, and 237 and 148 converged-fit BIC BF values above 10. Their BH
sets overlap at 1,113 tests, so only 12 logistic-normal calls are absent from
the multinomial result.

The condition analysis converged for 0.941 of tests and gave 5,503 nominal
calls, 191 BH calls, and 521 converged-fit BIC BF values above 10. The 191 BH
calls represent 137 blocks and 135 genes; 175 also have BIC BF above 10. They
are distributed across cell types: EX1, EX2, EX5, EX3, and EX4 contribute 63,
31, 21, 19, and 14, with smaller counts in other neuronal and glial groups.
The nominal rate among converged condition tests is 0.0799, compared with
0.2345 for multinomial cell-type tests. The merger initially summarized BF
counts without its convergence filter; the per-test table and BH values were
unaffected. The summary now applies the same nested-convergence requirement
to BF counts, and all mouse and human summaries were regenerated.

## 2026-08-05: matched external-method result tables

The differential-splicing manuscript now reports matched Tealeaf,
LeafCutter, MAJIQ Heterogen, and scQuint results for both datasets. Each table
gives the number of tests, BH discoveries, and discovery fraction separately
for cell-type and condition effects. The mouse Tealeaf row uses the primary
empirically calibrated logistic-normal EC GLMM (243 cell-type and zero
condition discoveries); the human row uses the recommended multinomial
Laplace LRT (1,487 and 191 discoveries). After correcting the cell-type
omnibus summary, the corresponding LeafCutter, MAJIQ, and scQuint cell-type
discovery counts are 229, 27, and 828 in mouse and 238, 342, and 2,239 in
human. Condition counts are 40, 34, and 16 in
mouse and 475, 278, and 36 in human.

The tables explicitly avoid interpreting raw discovery count or fraction as
power. Tealeaf, LeafCutter, MAJIQ, and scQuint test EC blocks, intron clusters,
LSV edges, and intron groups, respectively, and method-specific filtering creates
substantially different test universes. The datasets have no ground-truth DS
labels. External-method cell-type entries are Simes omnibus summaries of
pairwise tests, whereas the Tealeaf cell-type hypothesis is fit directly as
an omnibus model.

## 2026-08-05: comparator audit and split-subject power benchmark

The original external-method summary counted every feature-by-pairwise-cell-
type row even though the manuscript described a feature-level Simes omnibus.
The omnibus table was computed but not used to construct `summary.tsv`. This
inflated mouse scQuint cell-type tests from 16,081 unique intron groups to
586,476 pairwise rows and human tests from 27,097 to 2,034,311. The same issue
affected LeafCutter and MAJIQ. The summarizer now uses the omnibus table for
cell type and retains pairwise rows for within-cell-type condition tests.

MAJIQ normalization had selected `raw_pvalue` through a fallback while the
manuscript claimed the 95th-percentile posterior statistic. MAJIQ's own
documentation says that posterior quantile is conservative but not calibrated
as a p-value, so BH must use the raw per-LSV-edge test p-value. MAJIQ edge IDs
are now built from genomic coordinates and gene metadata rather than unstable
within-gene row numbers. The corrected mouse cell-type omnibus counts are
1,020/229 for LeafCutter, 46/27 for MAJIQ, and 16,081/828 for scQuint
(tests/BH discoveries). Human counts are 2,155/238, 582/342, and 27,097/2,239.

scQuint's remaining larger universe is real but permissive. The mouse and
human junction bundles contain 80,986 and 107,912 annotated three-prime groups
with at least two junctions; 46,236 and 57,280 contain exactly two. Running
scQuint on pseudobulks with only three nonzero samples required per arm allows
16,081 and 27,097 groups into at least one cell-type contrast. Tealeaf instead
requires EC-identifiable block paths and block-level gene coverage.

A two-fold mouse reproducibility benchmark was implemented and submitted.
Each fold contains 24 mice, exactly four from each of six conditions. Tealeaf,
LeafCutter, scQuint, and MAJIQ are rerun independently in each half. Features
are mapped to genes, within-gene p-values are combined by Simes, and scoring
uses a pairwise shared universe: each comparator and Tealeaf are evaluated on
the intersection of genes eligible for both methods in both folds. This avoids
letting MAJIQ's small mouse universe determine every comparison. The
primary endpoint applies BH to the conjunction p-value
`max(p_fold0, p_fold1)`; secondary endpoints are held-fold nominal replication,
log-p rank correlation, and top-K overlap. A junction-coordinate map assigns
6,568 LeafCutter clusters uniquely to genes. Subject-label permutations remain
necessary to establish calibration before replicated discovery count is
described as power.

The first LeafCutter split arrays inherited a malformed Lmod module table from
the submission shell, causing most tasks to fail before R started. A direct
smoke contrast completed after resetting the module state. Both LeafCutter
arrays and their dependent summaries were resubmitted with that fix; the
failed logs were retained.

The completed split fits were initially held up at comparator summarization:
the summary launcher invoked the Python virtual environment without loading
its matching Python module, so both jobs exited before reading results. The
launcher now resets inherited Lmod state and loads the configured Python
module, matching the reproducibility scorer. No model fits need to be rerun.
The scorer's launcher also accepts `PYTHON_BIN`; plotting requires the
`benchmark` optional dependencies, which now explicitly include plotnine.
Metric tables are written before plotting, and plotting is skipped with a
warning when plotnine is unavailable. This keeps optional visualization
dependencies from blocking the statistical benchmark.

The repaired split-subject benchmark completed on 2026-08-06. On each
method-specific shared gene universe, replicated discoveries based on BH of
the conjunction p-value were 26 Tealeaf versus 30 LeafCutter (166 genes), 34
Tealeaf versus 66 scQuint (608 genes), and 2 Tealeaf versus 4 MAJIQ (10
genes). The MAJIQ comparison is too small to interpret. Held-out nominal
replication was 0.696 versus 0.700 for LeafCutter, 0.666 versus 0.767 for
scQuint, and 0.833 versus 0.786 for MAJIQ. These results do not support a
claim of improved Tealeaf power over the alternatives. Calibration by
subject-label permutation and investigation of Tealeaf's restrictive shared
gene universe remain necessary.

## 2026-08-06: paired cell-type EC-GLMM

The initial split comparison was not hypothesis-matched. LeafCutter, scQuint,
and MAJIQ were run for each eligible pair of cell types using only subjects
represented in both members of the pair; feature-level p-values were then
combined by Simes. Tealeaf instead fit one omnibus cell-type effect using all
covered pseudobulks. Its degrees of freedom ranged from 1 to 48, so sparse
high-dimensional alternatives could lose power relative to pairwise tests.

`run_ec_block_glmm.py` now supports `cell_type_pairwise`. For each gene and
cell-type pair, it retains subjects with at least the requested gene EC count
in both cell types, includes condition as a nuisance effect, and clusters the
random effect by the shared subject. Each block LRT then has only
`n_paths - 1` tested degrees of freedom. Pairwise p-values are combined by
Simes within block, followed by Simes across blocks within gene, matching the
comparator aggregation hierarchy. The split candidate manifests contain
7,671 tests, 1,293 blocks, and 935 genes in fold 0, and 5,140 tests, 1,036
blocks, and 750 genes in fold 1. The two folds have 34 and 31 eligible
cell-type pairs, respectively.

The first gene-level scorer matched only genes, which can still compare
different biological hypotheses when gene-specific coverage changes the set
of eligible cell-type pairs. The paired scorer now first intersects
`(gene, unordered cell-type pair)` keys across both methods and both folds.
It combines feature tests within each matched gene-pair and then combines the
matched pair p-values within gene. Thus neither method can gain discoveries
from a cell-type contrast unavailable to the other method.

An audit of the base candidate universe found that the complexity caps are not
the main source of attrition: 10,681 genes pass the limits of 10 supported
isoforms and 128 ECs, while only 34 additional genes would pass if those caps
were raised to 50 isoforms and 512 ECs. The per-pseudobulk gene-count cutoff
is consequential. Reducing it from 25 to 10 increases the paired split
manifests from 7,671 to 31,004 tests (935 to 2,503 genes) in fold 0 and from
5,140 to 23,055 tests (750 to 2,238 genes) in fold 1. The lower threshold must
be evaluated by held-out reproducibility and null calibration rather than by
raw discovery count alone.

For a fixed gene and tested cell-type pair, every block-specific alternative
is the same unrestricted cell-type effect on the isoform logits. Blocks differ
only in which path contrast is removed under the null. The implementation now
expresses this shared alternative in a canonical unrestricted tensor basis,
maps the first block-null warm start into that basis by least squares, and
caches the fitted alternative by gene, retained rows, and method. This avoids
refitting identical alternatives for genes with multiple tested blocks while
leaving each block-specific null and LRT degrees of freedom unchanged.
Candidate sharding now treats all tests sharing `(gene, retained rows)` as an
indivisible group and greedily balances those groups by test count. The prior
stride partition could place reusable alternatives in different workers and
defeat the cache.

For inexpensive calibration checks, `--pairwise-null-seed` now performs a
deterministic label swap within each paired subject. The swap is shared across
genes for a given cell-type pair and preserves one observation from each level
per subject. `--max-candidates` selects a reproducible genome-wide subset
before grouped sharding, allowing the null p-value distribution to be checked
without a full permutation run.

The first 1,000-test fold-0 label-swap diagnostic rejects the plain
multinomial asymptotic calibration. Of 1,000 fits, 984 converged; 10.98% had
`p <= 0.05`, 2.03% had `p <= 0.001`, and 21 survived BH at 0.05 under the
permuted labels. The uniform KS p-value was `1.43e-17`. The excess tail is
consistent with subject-specific cell-type heterogeneity that is not captured
by a shared subject random intercept. Observed replicated discoveries from
this model cannot be interpreted as calibrated power. A matched null run with
the logistic-normal observation-noise model was submitted to test whether its
per-observation isotropic logit random effect controls this heterogeneity.
In the first 908 converged observation-noise null fits, 7.49% had
`p <= 0.05`, 0.44% had `p <= 0.001`, and none survived BH at 0.05. This is a
substantial tail correction but still shows mild nominal inflation. Full
observed split fits are required to determine whether the better-calibrated
model retains enough replicated power.

The label-independent coverage sweep shows why filtering must be tied to
calibration. At median gene UMI at least 150, the uncalibrated plain model has
9 versus 2 replicated genes for scQuint and 6 versus 4 for LeafCutter, but its
permuted nominal rate remains 12.7%, so this is not valid evidence of higher
power. For the observation-noise model, a median gene UMI threshold of 100 has
a 5.8% null nominal rate, no null p-values below 0.001, and no BH null
discoveries in the calibration subset. The scorer now supports this
label-independent threshold; 100 UMIs is the calibration-selected threshold
for the final noise-model comparison.

The calibrated noise model at 100 UMIs ties scQuint at 11 replicated genes
and remains behind LeafCutter. At 150 UMIs it exceeds scQuint (8 versus 2) and
LeafCutter (4 versus 3), with higher held-out replication, but choosing that
cutoff from the observed benchmark would overstate the evidence. The
calibration-compatible 75-UMI threshold still trails both comparators.

A joint gene-level EC test was added to address block-level multiplicity. For
each gene and eligible cell-type pair, its null allows condition effects but
no cell-type effect on any isoform logit; its alternative allows an
unrestricted cell-type effect with `n_isoforms - 1` degrees of freedom.
Block-level tests remain the localization stage after a joint discovery. This
produces 5,710 gene-pair tests across 935 genes in fold 0 and 3,930 tests
across 750 genes in fold 1, replacing 7,671 and 5,140 block-pair tests.

A null-BIC hybrid selector was also implemented. It selects the
logistic-normal model when twice its improvement in the null objective exceeds
`log(n_subjects)` and otherwise uses the plain multinomial result. On the
label-swap subset it selects noise for 11.5% of tests, has an 8.1% nominal
rate, and has no BH null discoveries. It does not improve the 100-UMI split
benchmark over the noise model, so it is retained as an ablation rather than
the recommended test.

The observation-noise model was rescored on the independently selected
high-confidence universe of at least eight paired subjects and at least 75
median gene UMIs. It still trails LeafCutter (5 versus 11 replicated genes)
and scQuint (12 versus 20), although its held-out nominal replication against
LeafCutter is higher (0.889 versus 0.780). Filtering alone therefore does not
establish a power advantage.

## 2026-08-07: retire unused paired EC extension

The unused two-cell-type EC extension and its dedicated CLI methods were
removed. No reported TeX result used this model. Pairwise cell-type contrasts,
the random-intercept EC GLMM, and the observation-noise model remain available.

## 2026-08-06: calibrated paired junction benchmark

The cross-method cell-type benchmark now includes a Tealeaf paired junction
log-ratio test on exactly the filtered junction groups used by scQuint. For
each subject and two-cell-type contrast, the test adds a 0.5 pseudocount,
takes the within-sample CLR of the junction counts, differences the two CLR
vectors within subject, and projects to Helmert coordinates. A two-junction
group uses a paired one-sample t test; a larger group uses Hotelling T-squared
with its finite-sample F reference. Tests require at least eight paired
subjects, more subjects than log-ratio dimensions, and a full-rank sample
covariance.

The eight-subject threshold is calibration selected. Across 12,760 eligible
tests from 32 paired-label-swap families, the empirical rates at p-value
thresholds 0.05, 0.01, 0.001, and 0.0001 are 0.05086, 0.00964, 0.000627, and
0.0000784. Recomputing BH within each filtered null family produces two calls
across 32 families. The paired-null generator was corrected to exclude both
the all-unswapped and all-swapped patterns: for a two-sided paired test either
pattern reproduces the observed contrast and is not a null randomization.

Subjects were divided into two condition-balanced halves. Every method was
rerun in each half, feature tests were combined to genes by Simes, and each
Tealeaf-comparator comparison was restricted to the genes and cell-type pairs
available to both methods in both halves. BH was applied to the conjunction
p-value `max(p_fold0, p_fold1)`. With native comparator p-values, Tealeaf has
193 versus 36 replicated genes for LeafCutter, 646 versus 154 for scQuint,
and 11 versus 6 for MAJIQ. Tealeaf's held-half nominal replication and
fold-to-fold Spearman correlation are respectively 0.906 and 0.685 against
LeafCutter, 0.788 and 0.642 against scQuint, and 0.958 and 0.713 against
MAJIQ.

The corrected null nominal 0.05 rates are 0.0620 for LeafCutter, 0.0837 for
scQuint, and 0.0260 for MAJIQ. Mapping each comparator's p-values through its
pooled method-specific null CDF reduces its replicated discoveries to 21,
134, and 0, respectively; Tealeaf remains at 193, 646, and 11. Therefore the
power advantage is present under native p-values and grows after a
label-independent calibration sensitivity analysis.

Several lower-power or invalid approaches remain useful negative controls.
The paired dominant-two-path Dirichlet-multinomial model produced at most six
replicated Tealeaf genes in preliminary split fits. A junction multinomial
GEE had a 0.13 null rate at nominal 0.05 and extreme zero p-values, so it is
not valid for discovery. A Laplace Dirichlet-multinomial EC smoke benchmark
was instead too conservative, with no p-values below 0.05 in either observed
or null fits. These approaches are not recommended for the benchmark.

## 2026-08-07: expanded MAJIQ sensitivity analysis

The stringent split-subject MAJIQ comparison was limited to 12 common genes.
This was not caused by Tealeaf's eight-paired-subject filter: zero-, four-,
and eight-subject Tealeaf thresholds all produced the same MAJIQ universe.
The two original MAJIQ folds contained only 700 and 1,096 feature-contrast
rows, with 46 unique LSV edges from 13 genes. The common build had supplied
all 936 pseudobulks as one experiment group while retaining MAJIQ's default
50% build support. Together with PSI-coverage defaults of 10 reads and three
nonzero bins per pseudobulk and 50% Heterogen support per arm, this selected
only ubiquitous, highly expressed genes.

The MAJIQ wrappers now expose build support, PSI-coverage read/bin thresholds,
and Heterogen support as parameters. The original defaults remain unchanged.
The relaxed UMI sensitivity analysis uses absolute support in three
pseudobulks at build, at least three UMIs in one nonzero bin per quantified
pseudobulk, and at least three quantified pseudobulks in each Heterogen arm.
Strict/nonredundant LSV selection is retained. Existing SJ files are reusable
through an independently configurable SJ directory, so changing build
filters does not require BAM reprocessing.

The relaxed IMM-versus-PER smoke contrast expanded from zero to 985 edge tests
across 250 genes. Genome-wide relaxed runs completed all 190 contrasts in
each subject fold plus 32 non-identity paired-label null permutations. The
folds contain 20,954 and 18,914 unique LSV features. After matching genes and
cell-type pairs in both folds and requiring eight Tealeaf paired subjects,
the common universe expands from 12 to 1,630 genes and from 133 to 24,722
gene-pair tests.

On this expanded universe Tealeaf has 490 replicated genes, held-half
replication 0.791, and fold Spearman correlation 0.638. MAJIQ has 50
replicated genes, held-half replication 0.723, and correlation 0.369 using
its native p-values. Across 178,793 relaxed MAJIQ null edge tests, nominal
rates at 0.05, 0.01, and 0.001 are 0.03197, 0.00456, and 0.000151, and no
permutation family has a BH discovery. Mapping observed MAJIQ p-values through
this conservative empirical null raises its replicated count to 89, still
well below Tealeaf's 490. Thus Tealeaf's advantage is not an artifact of the
original small MAJIQ universe.

## 2026-08-07: pseudobulk terminology and block contrasts

The differential-splicing documentation now calls each biological observation
a pseudobulk rather than a row. Matrix rows, design rows, and result-table rows
retain their algebraic or tabular meaning. Primer-specific observations are
described as EC count vectors within a pseudobulk, because the two primer
protocols share one latent pseudobulk composition.

The block-contrast description now separates the two components of the test.
The Helmert-derived basis encodes independent local path log-ratios, while the
tested design identifies the cell-type or condition differences applied to
those directions. The null retains factor-associated isoform changes
orthogonal to the local path balance; the alternative adds only the
factor-by-path directions.

The primary-method text now gives the implemented construction explicitly.
It lifts Helmert path rows to isoforms, centers the resulting full-logit
matrix, and applies the reference-logit map. Centering cancels under this map
because reference logits are invariant to a common shift. The nuisance-space
orthogonality is imposed in centered full-logit coordinates before both the
tested and nuisance bases are converted to reference logits.

## 2026-08-07: TeX source formatting

Hard wrapping was removed from `docs/differential.tex`. Prose, captions, list
items, algorithm statements, and table rows now occupy one logical source line;
display equations retain their semantic multiline layout. A whitespace-free
comparison confirmed that the document content did not change.

## 2026-08-07: Bouchard CAVI for the EC GLMM

I added a deterministic quadratic variational alternative for the multinomial
EC GLMM. EC-to-isoform responsibilities lower-bound each positive EC
numerator. Bouchard's softmax majorization followed by the Jaakkola--Jordan
softplus bound makes every subtracted primer normalizer quadratic in the
isoform logits. This gives closed-form updates for the full within-mouse
Gaussian factors and weighted least-squares updates for fixed effects. Mouse
and pseudobulk variance components use empirical-Bayes second-moment updates,
so the implementation is coordinate-ascent variational EM rather than fully
Bayesian CAVI for every parameter.

The reusable implementation is in `tealeaf/sc/ec_glmm_full.py`. It supports
the random-intercept and logistic-normal variants, returns the same packed
posterior representation as tilted full-covariance VI, and exposes both pure
CAVI and CAVI followed by tilted L-BFGS refinement. Both gene-level and
block-level EC runners expose these methods. The derivation and variable
dimensions are in `docs/cavi.tex`; `docs/differential.tex` gives the high-level
description and comparison.

In six null and six effect simulations per setting, CAVI initialization did
not improve tilted convergence. In the three-isoform logistic-normal setting,
direct and hybrid fitting each converged in four of six null replicates;
direct fitting converged in four of six effect replicates and the hybrid in
three. In the random-intercept setting, the hybrid reduced median tilted
iterations by two to three but added CAVI time and left final objectives and
errors unchanged. Standalone CAVI optimized a substantially looser bound and
had larger effect error.

The exact 65-gene logistic-normal screens gave a stronger negative result.
The established Laplace-initialized tilted procedure converged under both
nested models for 52 condition genes and 56 cell-type genes. Replacing its
initializer with 25 CAVI sweeps reduced those counts to 11 and 10 and rescued
none of the historical failures. Pure CAVI converged under both models for
five genes in each contrast by 200 sweeps. On genes where both tilted fits
converged, their final objectives were close, indicating that the CAVI start
impedes reaching the same optimum rather than identifying a better one. CAVI
therefore remains an opt-in reference method, not the default initializer.

## 2026-08-07: tilted L-BFGS convergence and numerical precision

I profiled the full-covariance logistic-normal optimizer on the two exact
65-gene real-data screens. Eight fixed-point updates for the correlated
tilted normalizer matched 25 updates, and differentiating through the local
updates changed the gradient by less than \(10^{-8}\). The production path now
uses eight updates and the envelope gradient.

The main conditioning problem was a ridge between latent posterior moments
and their estimated prior standard deviations. L-BFGS now optimizes each
latent mean and Cholesky row after scaling by the corresponding prior standard
deviation to the 0.75 power. This is an invertible optimizer
reparameterization; returned parameters and the variational family are
unchanged. A fit that exhausts the strict iteration limit receives at most 50
continuation iterations with a relaxed relative-objective tolerance.

In the contemporaneous full screens, the original 25-update differentiated
schedule converged for both nested models in 40 of 65 genes in each contrast.
The new schedule converged in 64 of 65 before a separate numerical edge-case
fix. Median time per nested pair decreased from 6.47 to 5.86 seconds for cell
type and from 6.76 to 5.97 seconds for condition. The new objectives were not
traded for convergence, median changes were positive in all four nested-model
comparisons, and the worst decrease was 0.33 objective units.

The remaining gene had no ECs for one primer. The old JAX likelihood evaluated
its zero primer total times a negative-infinite normalizer and returned NaN.
Empty primer matrices now contribute zero in the shared likelihood and both
tilted implementations. The previously failing gene then converged under both
nested models in both contrasts, giving corrected joint convergence of 65 of
65. A regression test exercises the empty-primer case.

I also compared numerical precision on 24 difficult fits. fp32
reduced median time from 11.38 to 8.58 seconds and reported convergence for
all genes, but its null and alternative objectives were lower than fp64 by
median 0.22 and 0.68, respectively. The nested gain changed by median absolute
0.45 and by as much as 9.27. On one hard gene, bfloat16 produced coarse,
wrong objectives, an abnormal null termination, and longer runtime. A
25-iteration fp32 warm-up followed by fp64 refinement restored objective
accuracy and found one alternative optimum 7.52 units above the direct fp64
fit, but increased median time by 34%. This mixed path is a useful sensitivity
restart, not a default initializer. Neither lower-precision format is suitable
for final fitting.

## 2026-08-07, comparator rerun and resource accounting

The affected mouse Tealeaf analysis was rerun after the tilted L-BFGS convergence update. The cell-type family retained 4,825 tests, 4,822 nested fits converged, 771 tests were nominally significant, and 252 passed BH at the primary 90% GPD threshold. Thresholds of 80%, 85%, 90%, and 95% gave 242, 240, 252, and 245 calls, with 238 shared by all four thresholds. Held-out null rejection rates were 0.04999, 0.00999, and 0.000995 at nominal levels 0.05, 0.01, and 0.001. The condition family retained 14,833 tests, 14,827 nested fits converged, 1,397 were nominally significant, and none passed BH. Its corresponding null rates were 0.05001, 0.00999, and 0.000998.

Tilted variational fitting now caches the compiled value-and-gradient function within each shard when the design and tensor shapes are unchanged, while passing the observed or simulated EC counts as a dynamic argument. The same executable is reused across the observed fit, continuation attempts, and five calibration replicates without changing the mathematical objective. A numerical comparison with altered counts matched an uncached fit in objective and fitted parameters. This reduced aggregate worker time for the two mouse families from 185.5 to 36.1 hours, a 5.1-fold improvement, with 1.77 GiB peak worker memory.

Resource accounting was added for the matched Tealeaf, LeafCutter, MAJIQ, and scQuint analyses. From prepared method inputs, mouse aggregate worker hours and peak worker GiB were 36.1 and 1.77 for Tealeaf, 11.7 and 1.28 for LeafCutter, 6.6 and 1.38 for MAJIQ, and 15.3 and 3.99 for scQuint. The final human values were 463.5 and 10.63 for Tealeaf, 38.6 and 3.26 for LeafCutter, 9.5 and 1.92 for MAJIQ, and 94.1 and 4.76 for scQuint. Tealeaf time includes every failed resource attempt and its successful retry. Relative to the previous human Tealeaf run, fitting only the recommended multinomial model with adaptive mode steps reduced aggregate worker time 3.1-fold and peak memory by 25%.

A controlled benchmark on 20 human cell-type tests compared 30 and 12 Newton updates for each random-effect posterior mode inside the Laplace objective. Both settings converged for every nested fit and had Spearman statistic correlation 1.0. Median absolute differences in the null objective, alternative objective, statistic, and p-value were $5.7\times10^{-14}$, $1.8\times10^{-10}$, $3.7\times10^{-10}$, and $9.5\times10^{-12}$. The largest p-value change was 0.012 in a nonsignificant test. Twelve updates reduced elapsed time from 2,745 to 1,019 seconds, a 2.7-fold improvement, and reduced peak memory from 7.00 to 6.70 GiB. A complete 12-step cell-type run nevertheless reduced joint convergence from the historical 0.949 to 0.921, showing that the 20-test benchmark underrepresented difficult fits. The analysis runner therefore uses 12 mode updates on the first attempt and 30 on continuation attempts, while the lower-level fitting API retains its conservative 30-update default. The final cell-type analysis used a 30-step continuation for at least one nested model in 1,354 of 11,597 tests and achieved 0.946 joint convergence, 2,557 nominal tests, 1,458 BH discoveries, and 236 BIC Bayes factors above 10. The condition analysis used a continuation in 6,746 of 73,002 tests and achieved 0.939 joint convergence, 5,475 nominal tests, 195 BH discoveries, and 521 Bayes factors above 10.

## 2026-08-08, harmonized Laplace optimization

The Laplace fitter now offers implicit differentiation of each converged random-effect posterior mode. For outer parameter vector $\theta\in\mathbb R^p$, cluster mode $u_m\in\mathbb R^d$, mode score $s_m=\partial\ell_m/\partial u_m\in\mathbb R^d$, and positive mode Hessian $H_m=-\partial s_m/\partial u_m\in\mathbb S_{++}^d$, the stationarity equation $s_m(\widehat u_m,\theta)=0$ gives $\partial\widehat u_m/\partial\theta=H_m^{-1}\partial s_m/\partial\theta\in\mathbb R^{d\times p}$. The reverse pass therefore solves a $d\times d$ system and applies the score pullback instead of retaining the computation graph for all Newton updates. This changes the derivative implementation, not the fitted Laplace objective.

Counts and tensor-valued fixed effects are dynamic arguments to the compiled objective. The runner reuses a compiled function among contiguous contrasts for the same gene, pseudobulk set, and shape, then clears the cache before the next group. An initial unbounded cache exceeded available memory, so its incomplete artifacts were archived outside the repository and the cache lifetime was bounded. A regression test with changed counts and fixed-effect tensor matched uncached parameters and objective to numerical precision. On a 72-test human condition shard, bounded gene-level caching and per-test compilation produced identical results for all tests and 69 jointly converged fits; caching reduced elapsed time from 2,273 to 2,228 seconds with similar peak memory, about 6.65 and 6.66 GiB. Early termination of the Newton mode loop produced identical 20-test results but did not improve elapsed time, so it was not retained.

On 20 human cell-type tests, implicit and unrolled mode derivatives both converged for every nested fit. Objectives, statistics, and p-values were close, with a largest p-value difference of 0.015 and no change among significant tests. Implicit differentiation reduced elapsed time from 1,119 to 752 seconds, aggregate CPU time from 39:01 to 21:35, and peak memory from 6.49 to 6.30 GiB. On 20 human condition tests it reduced elapsed time from 1,498 to 651 seconds and aggregate CPU time from 50:15 to 17:38. A complete 72-test human condition shard took 2,228 seconds and 6.66 GiB with implicit differentiation versus 4,781 seconds and 6.88 GiB with unrolled differentiation, improvements of 2.15-fold in elapsed time and 2.56-fold in aggregate CPU time. Both produced 69 jointly converged tests. Median absolute differences in the statistic and p-value were below $10^{-12}$, the largest statistic difference was 0.221, and the largest p-value difference among significant tests was $9.0\times10^{-6}$. The optimized procedure uses an ordinary multinomial EC GLMM, a 12-update posterior-mode solve, implicit mode differentiation, a bounded gene-and-shape compilation cache, a Laplace likelihood-ratio test, and BIC evidence; a 30-update continuation is used when a nested fit fails convergence. The same model, optimizer, inference rule, and continuation are used for mouse and human, while dataset-specific minimum-subject screening remains unchanged.

The four harmonized production arrays completed without failed shards. Mouse cell-type results contained 4,825 tests, 94.22% joint convergence, 708 nominal tests, 323 BH discoveries, and 160 BIC Bayes factors above 10. Mouse condition results contained 14,833 tests, 92.89% joint convergence, 677 nominal tests, three BH discoveries, and 51 BIC Bayes factors above 10. Human cell-type results contained 11,597 tests, 94.69% joint convergence, 2,557 nominal tests, 1,459 BH discoveries, and 237 BIC Bayes factors above 10. Human condition results contained all 73,144 screened tests, 94.02% joint convergence, 5,497 nominal tests, 193 BH discoveries, and 522 BIC Bayes factors above 10. The 193 human condition discoveries represent 137 blocks and 135 genes, and 176 also exceed the BIC evidence threshold.

The common Laplace procedure used 118.8 aggregate mouse worker hours with 4.74 GiB peak worker memory and 729.4 aggregate human worker hours with 6.83 GiB peak worker memory. The corresponding scQuint totals were 15.3 and 94.1 hours, so Tealeaf remained about 7.8-fold slower in both datasets after the tested optimization. The new human peak memory is lower than the previous 10.63 GiB, but the aggregate time is higher than the previous 463.5-hour adaptive campaign because condition shards have a heavy runtime tail. The representative-shard benchmark established that implicit differentiation is still faster than unrolled differentiation for the same implementation. Relative to mouse, human has 4.3-fold more tests; its median test has five rather than three isoforms, 26 rather than nine ECs, and five rather than two block paths. These dimensions enlarge every likelihood, gradient, and mode-Hessian calculation and explain why human is slower. Early mode stopping did not reduce real-data time, fp32 and bfloat16 did not preserve final-fit accuracy, and CAVI initialization reduced convergence, so no tested alternative justified another production rerun.

## 2026-08-08, pre-collapsed local-path comparator

A lower-dimensional comparator now estimates a label-independent pooled isoform composition for each gene and covered pseudobulk set, normalizes those weights within each represented local path, and replaces the path's compatibility column by the weighted mean of its isoform columns. Isoforms outside the tested block map remain separate normalization nuisances. The comparator continues to fit the raw paired-primer EC vectors with a multinomial random-intercept GLMM, so it retains EC multinomial covariance, ambiguity, and subject correlation. It does not retain isoform-specific changes within one path, and it is exact only when the conditional isoform mixture within each path is constant and correctly estimated.

The exactness audit supports treating this as an approximation rather than an algebraic reparameterization. In the two mouse subject-fold manifests, only 0.22% and 0.24% of multi-isoform path-by-primer compatibility groups had identical columns. The mean latent-dimension reduction was 11.0% and 11.3%. In the full discovery manifests, only 29.7% of mouse tests and 12.8% of human tests reduce at all; the mean fractional reductions are 10.6% and 4.4%, respectively. Median dimension changes from three to three in mouse and five to four in human. A matched four-test smoke fit took 66.9 seconds in isoform space and 64.5 seconds in pooled path space, with every nested fit converged. This small pilot does not indicate a material speed improvement.

The complete mouse discovery comparison contained the same 4,825 tests in each latent space. Isoform space had 94.22% joint convergence, 323 BH discoveries among 4,546 converged tests, and a 7.11% discovery rate. Path space had 94.90% convergence, 354 BH discoveries among 4,579 converged tests, and a 7.73% discovery rate. Restricting both methods to 4,496 jointly converged tests gave 319 and 348 discoveries, of which 299 were shared; the discovery-set Jaccard index was 0.813. Thus the path model adds calls, but discovery yield alone cannot distinguish added power from effects of its fixed-mixture assumption.

Both subject folds were rerun for isoform and path space with the current implicit-gradient Laplace procedure. Reproducibility scoring required each block test to converge in both latent spaces and both folds, leaving 1,009 block tests and 744 genes. Isoform space had 37 genes passing BH on the conjunction p-value, held-half nominal replication 0.676, and fold Spearman correlation 0.433. Path space had 39 replicated genes, held-half replication 0.682, and correlation 0.395. The two-count difference and 0.006 held-half difference are small, while the rank correlation moves against path space. The extra full-data discoveries therefore do not produce a clear replication improvement.

Path collapse did not improve resource use. The complete mouse discovery run used 33.23 aggregate worker hours and 4.04 GiB peak memory in path space versus 31.18 hours and 3.97 GiB in isoform space, making path space 6.6% slower. In fold 0 the totals were 8.40 versus 8.17 hours and 3.64 versus 3.55 GiB; in fold 1 they were 8.26 versus 7.60 hours and 3.61 versus 3.25 GiB. The block-specific collapsed maps prevent reuse of a gene-level unrestricted alternative across blocks, and the pooled within-path fit adds preprocessing. These costs outweigh the small latent-dimension reduction.

A human sensitivity analysis used the 359 completed matched shards, containing 4,282 tests, before stopping the remaining queued tasks. Isoform and path convergence were 94.89% and 94.86%. Path space had 506 BH calls among converged tests versus 477 for isoform space, with discovery rates 12.46% and 11.74%. On 3,968 common-converged tests the counts were 496 and 472, with 457 shared calls and discovery-set Jaccard index 0.894. The exact matched shards used 36.55 path-space worker hours versus 33.69 isoform-space hours, an 8.5% slowdown, and peak memory was 6.65 versus 6.39 GiB. This subset is not a genome-wide human discovery analysis, but it supports the mouse conclusion that the small human dimension reduction does not improve speed and that path space tends to add calls.

## 2026-08-08, estimate-once path Wald ablation

The comparator-like ablation estimates each pseudobulk's local-path ILR response once from its paired-primer EC vectors, retains the full conditional Fisher covariance among ILR coordinates, and then runs a precision-weighted regression with a mouse-clustered sandwich Wald test. The pooled isoform baseline is label-independent and cached across blocks from the same gene and covered pseudobulk set. Nuisance and tested fixed-effect columns and mouse clusters match the primary EC GLMM. Unlike the EC GLMM, the ablation freezes observation covariance and does not revisit EC ambiguity during fixed-effect inference. Pseudobulks with unidentifiable path covariance are omitted, and the block fails when the remaining rows do not identify the alternative design or when the number of mice does not exceed the number of design columns.

A 20-block mouse smoke run completed in 33.2 seconds. Eighteen blocks were testable; their median pseudobulk use rate was 1.0 and mean use rate was 0.973. The two failures lacked enough identifiable pseudobulks to span the regression design. Full mouse discovery, both mouse subject folds, and full human discovery were submitted as jobs 19650603--19650606.

The initial REML second stage retained 4,451 of 4,825 mouse tests, made 132 BH calls among eligible tests, and had a 2.97% discovery rate. On 4,170 blocks eligible for isoform-space EC GLMM, collapsed-path EC GLMM, and estimate-once REML Wald, the methods made 303, 336, and 126 calls; 79 estimate-once calls overlapped isoform-space calls. The complete three-method replication universe contained 936 blocks and 713 genes. Isoform-space, collapsed-path, and estimate-once REML Wald had 35, 40, and eight replicated genes. Their held-half replication rates were 0.668, 0.671, and 0.811, while fold rank correlations were 0.416, 0.375, and 0.677. This version selected fewer effects but its selected effects and full p-value ranking were more stable across halves.

The initial REML version used 3.30 aggregate mouse worker hours and 0.55 GiB peak memory, versus 31.18 hours and 3.97 GiB for isoform-space Tealeaf. Subject folds used 0.52 and 0.56 hours, versus 8.17 and 7.60 hours for isoform space. Thus resolving EC ambiguity once removed most of Tealeaf's runtime and memory cost, confirming that repeated joint EC likelihood and latent-Gaussian inference, not storing path covariance itself, creates the primary gap.

An attempted partial-identification sensitivity replaced infinite variance in unidentified ILR directions with diffuse finite variance while retaining inverse-Fisher covariance in identified directions. On high-dimensional human blocks this produced condition numbers around $10^{16}$--$10^{19}$, repeated linear-algebra warnings, and much slower REML fits. The sensitivity was stopped and not retained; the reported implementation omits pseudobulks whose complete path coordinate vector is unidentifiable.

Seven of 96 human REML shards remained in their first 20 tests after 15 minutes. These tail shards had 500--900 pseudobulks and up to 27 cell types, making repeated variance-component profiling the remaining bottleneck. The second stage was replaced by one precision-weighted fit with a mouse-clustered sandwich covariance and degrees-of-freedom correction. This keeps the estimated ILR responses and their off-diagonal conditional Fisher covariance, but represents mouse dependence marginally instead of estimating a random intercept. A 20-test human smoke run completed in 94.9 seconds with no variance-component runtime tail. The matched mouse discovery, mouse folds, and human discovery reruns are jobs 19650886--19650889.

In 1,000 paired null simulations with 30 clusters, the cluster-corrected sandwich statistic with a chi-square reference rejected 5.1%, 2.0%, and 0.2% at nominal 0.05, 0.01, and 0.001. Using the finite-cluster reference $W/q\sim F_{q,M-k}$ gave 4.4%, 1.0%, and 0.0%. The production ablation therefore uses the F reference. The initial sandwich arrays were stopped before completion and replaced after this calibration check.

The final audit found eight mouse discovery blocks, three fold blocks, and five human blocks for which the mouse count did not exceed the number of design columns. The initial implementation had used a chi-square fallback when the finite-cluster denominator degrees of freedom were zero. The code now makes these tests ineligible, a regression test covers this condition, and only the affected shards were rerun. All final reported tests have positive denominator degrees of freedom.

The final mouse cell-type universe contained 4,825 tests. Estimate-once Wald was eligible for 4,499 and made 3,069 BH calls, a 68.2% discovery rate, versus 323 calls among 4,546 eligible isoform-space EC GLMM tests and 354 among 4,579 collapsed-path EC GLMM tests. On 4,211 common-eligible blocks, the call counts were 303, 335, and 2,876. Estimate-once Wald overlapped 276 of the isoform-space calls but had a discovery-set Jaccard index of 0.095 because its selected set was much larger.

The three-method split-half universe contained 948 blocks and 724 genes. Isoform-space EC GLMM, collapsed-path EC GLMM, and estimate-once Wald made 35, 40, and 276 gene-level BH calls in both halves. Held-half nominal replication was 0.669, 0.685, and 0.751, while fold Spearman correlation was 0.429, 0.393, and 0.426. The extra estimate-once calls therefore produced substantially more replicated genes and did not reduce held-half replication, but the sixfold increase in discovery rate is large enough that realistic structured-null calibration is required before adopting the working-model test as the primary procedure.

The human cell-type universe contained 11,597 tests. Estimate-once Wald was eligible for 8,062 and made 5,906 BH calls, a 73.3% discovery rate. The isoform-space EC GLMM was eligible for 11,003 and made 1,487 calls, a 13.5% rate. On 7,685 common-eligible blocks the methods selected 5,633 and 1,170, with 873 overlapping calls and a Jaccard index of 0.147. The 69.5% human eligibility rate is lower than the 93.2% mouse rate because the ablation requires the complete local-path ILR vector and covariance to be identifiable in every retained pseudobulk.

Including the targeted denominator-degrees-of-freedom correction reruns, mouse discovery used 2.78 aggregate worker hours and 0.54 GiB peak memory; the two mouse folds used 0.40 and 0.30 hours with 0.33 and 0.36 GiB. Human discovery used 14.18 hours and 3.58 GiB. The matched cell-type isoform-space EC GLMM used 31.18 mouse hours and 3.97 GiB and 98.39 human hours and 6.55 GiB. Estimate-once Wald is therefore 11.2-fold faster on mouse and 6.9-fold faster on human. Its mouse time is below the complete LeafCutter, scQuint, and MAJIQ campaign totals of 11.71, 15.28, and 6.61 hours. Its human time is below LeafCutter's 38.57 and scQuint's 94.09 hours but above MAJIQ's 9.55 hours. Comparator feature units and campaign scopes differ, so these latter comparisons establish the runtime range rather than a per-test speed ratio. The result confirms that storing and using the off-diagonal path covariance is not the main expense; repeated raw-EC likelihood evaluation and latent-Gaussian inference under two nested models create the primary Tealeaf runtime gap.

## 2026-08-08, structured-null audit of estimate-once Wald

The 68--73% discovery rates prompted a design-preserving audit. For each block, local path estimates and their full conditional Fisher covariances were computed once. A label null then independently shuffled the observed multi-level cell-type labels within each mouse, preserving the mouse's observed rows, EC-derived response, covariance, condition, and cell-type counts. A second restricted wild-cluster null fit the nuisance-only precision-weighted model and multiplied each mouse's complete multivariate residual vector by one Rademacher sign. Thirty-two null families were used for full mouse discovery and 16 per fold. The earlier balanced two-level Gaussian simulation did not reproduce the real tested dimensions, cluster leverage, incomplete cell-type patterns, or covariance heterogeneity and was therefore not an adequate calibration check.

The CR1 finite-cluster F reference failed severely. Across 143,460 successful label-permuted tests, rejection rates were 0.658, 0.583, and 0.517 at nominal 0.05, 0.01, and 0.001. Across 143,904 restricted wild-cluster tests, rates were 0.603, 0.532, and 0.468. A label-permuted family produced a mean 2,844 BH calls, range 2,794--2,896, while a wild family produced a mean 2,581, range 2,539--2,616. The observed CR1 analysis made 3,069 calls, so most of that discovery count is explained by null miscalibration. Per-test empirical label-permutation ranks placed 5.16% of observed blocks at or below 0.05, but 32 draws give a minimum p-value of 0.0303 and cannot resolve genome-wide BH tails. Pooling structured-null p-values within exact numerator-degree and binned cluster-degree strata gave 215 nominal label-calibrated tests and 286 wild-calibrated tests, with zero BH calls under either calibration.

Cluster leverage is a major contributor but not a complete explanation. CR2 premultiplies each whitened mouse residual by the inverse square root of its cluster hat-matrix complement and changed little in the 20-block smoke audit. CR3 uses the full inverse complement. In that smoke set CR3 reduced label-null rejection from 0.465 to 0.125 and wild-null rejection from 0.347 to 0.0625. Genome-wide, however, CR3 still rejected 0.367 of label permutations and 0.338 of wild nulls at 0.05. Its null families averaged 1,430 and 1,277 BH calls. CR3 reduced the observed mouse count from 3,069 to 1,906, but matched null calibration again left zero BH discoveries. Miscalibration increases with the joint tested dimension, yet one-degree-of-freedom CR3 tests still rejected roughly 0.12--0.24 depending on cluster degrees of freedom, so restricting to small blocks is insufficient.

Path direction was compared in a parameterization invariant to each fold's reference cell type. Each fold-specific ILR coefficient matrix was transformed to centered-log path effects, then all cell-type-pair contrasts shared by both folds were compared. The common universe contained 985 blocks and 14,059 pair-by-path coordinates. Overall sign agreement was 0.568 and fold effect Spearman correlation was 0.151. Coordinates in blocks passing uncalibrated BH in both folds had agreement 0.577 and correlation 0.138. The corresponding null sign agreement was 0.510 for label permutations and 0.502 for wild nulls. Requiring absolute z-scores of at least 1.96 in both CR1 folds retained 1,035 coordinates from 279 blocks, with agreement 0.616 and correlation 0.173. Under CR3 this filter retained 217 coordinates from 69 blocks, with agreement 0.806 and correlation 0.415. Thus a small strong-effect subset has directional support, but neither analytic procedure supplies calibrated p-values for selecting it.

An evenly spaced human sensitivity used 11 of 96 shards, 912 eligible observed blocks, and 16 null families. Label-permutation rejection was 0.753 at 0.05, 0.723 at 0.01, and 0.687 at 0.001; wild-null rejection was 0.723, 0.688, and 0.652. Null families averaged 666 and 652 BH calls. Matched stratified null calibration left zero BH calls. The twelfth selected shard was stopped after becoming a runtime tail; the completed subset is already large enough to establish that the mouse failure generalizes.

The estimate-once analysis remains useful as a computational ablation: resolving EC ambiguity once moves runtime into the comparator range and shows that dense path covariance storage is not the main cost. Its discovery and replication counts are invalid for inference, however, and it must not replace the primary EC GLMM. A future fast alternative would need direct randomization inference with enough draws to resolve the BH tail, or a model whose cluster-level reference distribution is validated under these structured nulls.

## 2026-08-09, Gaussian mixed-model LRT and score Bayes-factor path ablation

The estimate-once path responses were coupled to a Gaussian marginal mixed model with their conditional EC covariance treated as known. Each block has an isotropic path-coordinate random intercept for mouse and an isotropic extra residual component. Null and alternative models separately maximize the ordinary likelihood, profile their fixed effects, and explicitly consider the variance-component boundaries. The LRT compares the profiled fits. A BIC approximation uses the independent mouse count as its sample size. A second Bayes factor integrates the null efficient score under normal effect priors with path-ILR standard deviations 0.1, 0.25, 0.5, and 1.0 and averages the four alternatives. The contrast prior is induced from independent latent cell-type effects, giving covariance $s^2\{(I+\bm1\bm1^\prime)\otimes I_d\}$ and making the Bayes factor invariant to the reference cell type.

The covariance solver uses a block-diagonal observation factor and a path-dimensional Woodbury correction within each mouse. Dense and Woodbury calculations on identical estimated responses agree to numerical precision. Woodbury reduced a 20-block high-dimensional human benchmark from 193.6 to 168.0 seconds, 13.2%, but increased aggregate mouse worker time from 4.46 to 4.83 hours while reducing the slowest shard from 502 to 418 seconds. It remains the default because it avoids a cubic factorization in pseudobulks per mouse and is safer for the human runtime tail.

The full mouse observed analysis retained 4,472 of 4,825 blocks. The analytic LRT made 200 BH calls, 4.47% of eligible tests. The score-mixture Bayes factor exceeded 10 for 277 blocks and 100 for 163; the BIC Bayes factor exceeded those thresholds for 83 and 51. The corresponding isoform-space EC GLMM made 323 calls among 4,546 eligible blocks. On 4,232 common-eligible blocks, the raw EC model selected 305 and the path LRT selected 195, with 105 shared calls and Jaccard index 0.266.

Mouse subject halves retained 963 common block tests and 725 common genes when compared with the isoform-space EC GLMM. The raw EC model produced 35 genes significant in both halves, held-half nominal replication 0.668, and fold log-p correlation 0.432. Gaussian path LRT produced 17 genes significant in both halves, held-half nominal replication 0.724, and fold correlation 0.716. Among all 14,014 shared pair-by-path coordinates, sign agreement was 0.573 and effect correlation was 0.148. Requiring absolute z-scores of at least 1.96 in both folds retained 171 coordinates from 54 blocks, with sign agreement 0.947 and correlation 0.932. Strong Gaussian path effects are therefore more directionally stable than the sandwich-Wald selections, although calibration remains a separate requirement.

An evenly spaced human observed sensitivity used 11 of 96 shards and retained 935 of 1,329 candidates. It made 153 analytic BH calls. The score-mixture Bayes factor exceeded 10 for 196 blocks and the BIC Bayes factor exceeded 10 for 19. Across eight within-subject label-permutation families, pooled LRT rejection was 0.0608, 0.0337, and 0.0209 at nominal 0.05, 0.01, and 0.001. Thus the bulk null is close to nominal at 0.05 but the extreme tail is inflated. Empirical tail calibration retained 37 LRT, 32 score-mixture, and 19 BIC-evidence calls at estimated FDR 0.05, but eight null families are too few to regard these counts as definitive.

All 96 mouse shards completed 32 within-mouse label permutations. Pooled LRT rejection was 0.0416, 0.0151, and 0.00478 at nominal 0.05, 0.01, and 0.001, and one null family produced a mean 11.0 analytic BH calls. The bulk is conservative but the deep tail is inflated. The restricted wild-cluster sensitivity covered 92 shards and gave rejection 0.0224, 0.00638, and 0.00150, with a mean 1.56 BH calls. Empirical tail-FDR calibration retained 12 observed LRT calls under label permutations and 194 under the wild null. The score-mixture Bayes factor retained 12 and 269, while the BIC approximation retained 27 and 99. The conservative label/wild intersections therefore contain 12 LRT, 12 score-mixture, and 27 BIC-evidence blocks.

Fixed Bayes-factor thresholds are not calibrated either. At score-mixture BF at least 10, the observed count was 277 and the label-permutation expected null count at the observed family size was 28.3; at BF at least 100 the counts were 163 and 12.3. BIC BF at least 10 gave 83 observed and 7.94 expected null calls, while BF at least 100 gave 51 and 3.37. The score Bayes factor does not solve the near-boundary path misspecification, while BIC's stronger tested-dimension penalty yields a somewhat larger empirically supported set. The estimate-once Gaussian model is a useful ranking and computational ablation, not an analytically calibrated replacement for the raw-EC GLMM.

## 2026-08-09, EC-GLMM outer optimization and batching audit

The production Laplace fit was profiled on matched mouse and human block sets. The existing runner already deduplicates supported partitions, groups candidates by gene and pseudobulk rows, reuses one unrestricted isoform-space alternative across blocks in that group, and reuses compiled objectives within a bounded gene-and-shape group. The 4,825 mouse candidates form 3,508 gene-row groups, so unrestricted-alternative reuse avoids 1,317 fits. The 11,597 human candidates form 10,175 groups and avoid 1,422 fits. In the matched baseline, objective evaluations occupied 92.8 of 211 seconds for mouse and 94.4 of 288 seconds for human, while the SciPy optimizer occupied 30.2 and 41.6 seconds. Compilation, argument construction, posterior summaries, and other per-fit work therefore remain material targets.

The Laplace fitter now exposes projected Adam and Adam followed by L-BFGS as experimental alternatives and records objective-evaluation and optimizer timings. Pure Adam used 300 bounded steps with a cosine learning-rate schedule. Learning rates 0.01 and 0.03 produced zero of 12 jointly converged mouse tests; on 11 human tests they produced zero and two jointly converged tests. Maximum absolute LRT-statistic differences from direct L-BFGS were 424 and 295 for mouse and 1,169 and 113 for human. Adam is easy to vectorize, but these results rule out replacing the line-searched optimizer with it. Adam followed by L-BFGS recovered all matched convergence checks, but maximum statistic differences were 1.74 for mouse and 3.79 for human because the warm-up sometimes selected a different local optimum. It remains a sensitivity restart rather than a default.

Muon was considered but not implemented. The outer parameters comprise block-contrast coefficients, positive scalar variance parameters, and cluster-specific posterior modes eliminated by the inner solve. Reshaping contrast coefficients into a matrix would make Muon's orthogonalized update depend on the arbitrary contrast basis, while treating mice as a matrix axis would couple independent random-effect blocks. There is no invariant matrix geometry analogous to a neural-network weight matrix. Momentum SGD lacks adaptive curvature, while RMSProp, Adagrad, AdamW, and Lion retain fixed first-order schedules and no line search; given Adam's large residual gradients and objective errors, none has a credible accuracy advantage here. A batched Newton, natural-gradient, or L-BFGS update is better aligned with the likelihood geometry, but independent fits require active masks and separate acceptance or line-search state.

Warm starts were also tested. Projecting the already fitted unrestricted alternative into later null-model bases reduced matched wall time from 220 to 183 seconds for mouse and from 303 to 216 seconds for human, but moved some local optima and changed the maximum statistic by 1.74 and 0.20. A previous-null warm start was not faster on the human subset and was not fully matched after the deterministic benchmark cache was rebuilt. Neither warm start is a safe production default unless multiple starts are evaluated and the best objective retained, which removes much of the speed gain.

Two independent one-core L-BFGS workers reduced matched wall time from 220 to 137 seconds for mouse and from 303 to 230 seconds for human. Aggregate CPU time changed from 397 to 350 seconds and from 654 to 644 seconds, while peak memory rose from 2.13 to 2.68 GiB and from 4.63 to 6.97 GiB. Mouse statistics were identical. One difficult human statistic differed by 0.279 under independent execution, with the same ten of eleven tests jointly converged, which is consistent with the local numerical sensitivity seen elsewhere. Existing production arrays already provide this coarse-grained parallelism; increasing worker concurrency trades memory for wall time but does not materially reduce total compute.

The genome-wide shape census found that only 813 mouse and 434 human candidates belong to exact tensor-shape buckets containing at least two tests. If pseudobulk and mouse counts may vary while isoform count, two primer EC dimensions, and fixed-effect dimension remain fixed, 3,726 mouse and 7,971 human candidates become batchable, with maximum bucket sizes 249 and 150. If a fully ragged optimizer groups only by isoform and fixed-effect dimension, 4,777 and 11,523 candidates become batchable. Cost-sorted padding at batch size eight inflates an EC-likelihood work proxy by 1.97-fold for mouse and 2.11-fold for human; at batch size 64 the factors are 3.43 and 3.42.

An exact likelihood-and-gradient prototype compared padded tensors with JAX ragged dot products. In the largest low-padding bucket, padding ratios were 1.13 for mouse and 1.34 for human; padded median evaluation times were 0.00065 and 0.00137 seconds, versus 0.00162 and 0.00296 seconds for ragged evaluation. In the deliberately worst-padding buckets, ratios were 3.15 and 7.09; ragged evaluation was 1.09-fold and 1.86-fold faster and compilation was 2.1-fold and 3.0-fold faster. Values and gradients agreed numerically. The appropriate design is a hybrid scheduler with small, cost-sorted padded microbatches for ordinary buckets and ragged kernels for high-padding tails, not universal ragged storage. Completing this requires batched posterior-mode solves and independent L-BFGS histories and line searches; a shared scalar objective with Adam would be simpler, but the optimizer audit shows it is not statistically reliable.

A conservative screening audit ranked path-eligible mouse tests by the estimate-once Gaussian path LRT while forcing every ineligible block through the exact EC GLMM. Fitting 53.7% of exact tests recovered 86.0% of the 315 full-analysis discoveries, and fitting 76.8% recovered 94.9%. The path screen itself used 4.46 aggregate worker hours. Score-mixture and BIC Bayes-factor ranks recovered fewer discoveries. Screening therefore offers only modest net savings if high discovery recall is required and was not adopted.

### Independent batched L-BFGS

A separable bounded L-BFGS implementation now retains an independent curvature history, projected gradient, convergence state, and Armijo backtracking line search for every fit. A failed or converged fit is frozen by its active mask while the remaining fits continue. Objective evaluation is vectorized, but JAX still evaluates inactive lanes in the compiled padded batch, so masking protects optimizer state without eliminating all wasted arithmetic. Unit tests cover different per-fit Hessians, bound optima, different iteration counts, and a line search in which one fit must backtrack while another is already converged.

The full benchmark pads pseudobulks, mice, and primer-specific EC dimensions within a fixed isoform, null-coefficient, and alternative-coefficient bucket. It evaluates the complete multinomial Laplace objective, including 12 posterior-mode Newton updates, implicit mode differentiation, and separate null and alternative optimization. Empty primers and padded EC rows require conditional likelihood evaluation; without it, second derivatives of the all-zero padded mapping are undefined. The benchmark also deduplicates the unrestricted alternative by gene and pseudobulk set, matching the production runner.

For a representative three-isoform mouse bucket with 17 null and 22 alternative coefficients, eight nested tests took 24.2 seconds in batch and 132.6 seconds sequentially, a 5.48-fold speedup. Batched optimization converged all eight tests and scalar L-BFGS-B converged seven. On the seven commonly converged tests, median and maximum absolute LRT-statistic differences were 0.00191 and 0.00898. For a representative three-isoform human bucket with four null and six alternative coefficients, eight tests took 19.1 versus 123.1 seconds, a 6.43-fold speedup. Both implementations converged every test; median and maximum statistic differences were 0.000066 and 0.000584.

The gain shrinks with model dimension. A four-test mouse bucket with seven isoforms, 76 null coefficients, and 84 alternative coefficients took 27.0 seconds in batch and 41.5 seconds sequentially after reducing four identical unrestricted alternatives to one fit, a 1.53-fold speedup. Every nested fit converged and the LRT statistics were identical. A four-test human bucket with ten isoforms, 18 null coefficients, 252 alternative coefficients, roughly 900 pseudobulks, and up to 97 ECs in one primer took 268.7 seconds in batch and 165.4 seconds sequentially, so batching was 1.62-fold slower. Batched and scalar optimization jointly converged on two tests, whose absolute statistic differences were 6.88 and 14.80. The batched fits found lower null objectives on all four tests, but one alternative objective was 6.67 units worse, confirming sensitivity to distinct local optima rather than a simple arithmetic mismatch.

The tested implementation is therefore promising only as a hybrid path. Cost-sorted batches of about eight are useful for low-padding, moderate-parameter fits, while the high-dimensional human tail should remain scalar or use smaller or ragged batches. Switching the optimizer changes trajectories and occasionally optima even when both convergence criteria pass, so any production use requires structured-null calibration and a comparison of discovery calls. No production result was rerun with batched L-BFGS.

## 2026-08-09, full-data selective batching audit

The selective scheduler was tested on all 4,825 mouse and 11,597 human cell-type block tests. It admitted cost-sorted batches of eight only when at least four fits shared isoform and coefficient dimensions, the padded likelihood-work ratio was at most 2, and both nested models had at most 128 outer coefficients. The high-dimensional tail remained scalar. Unrestricted alternatives were deduplicated by gene and covered pseudobulk set. Failed batched fits received scalar continuation. This routed 2,349 mouse nulls and 1,956 mouse alternatives through an initial batch, including continuations, and routed 2,147 human nulls and 5,708 human alternatives through an initial batch.

The mouse hybrid converged under both nested models for 4,540 tests versus 4,546 for the scalar production run. It made 338 BH discoveries versus 323, with 313 shared, 25 hybrid-only, 10 scalar-only, and discovery-set Jaccard index 0.899. Among 4,483 tests converged by both procedures, median absolute LRT-statistic drift was 0.00750, Spearman correlation was 0.953, 484 statistics moved by more than one, 98 moved by more than five, and the maximum absolute change was 54.0. Elapsed worker time was 30.48 hours versus 31.18, only 2.2% lower, while peak worker memory fell from 3.97 to 1.40 GiB.

The human hybrid converged under both nested models for 10,889 tests versus 10,981 for scalar production. It made 1,490 BH discoveries versus 1,459, with 1,441 shared, 49 hybrid-only, 18 scalar-only, and Jaccard index 0.956. Among 10,846 common-converged tests, median absolute statistic drift was 0.0381, Spearman correlation was 0.997, 1,485 statistics moved by more than one, 363 moved by more than five, and the maximum change was 674.9. Elapsed worker time was 125.21 hours, including one transient failed attempt and its successful retry, versus 98.39 for scalar production, a 27.3% increase. Peak memory fell from 6.55 to 4.19 GiB.

The full-data failure has two implementation causes that the homogeneous microbenchmarks did not exercise. Production scalar fitting keeps candidates from the same gene and covered pseudobulk set together, reuses compiled objectives within a bounded gene-and-shape cache, and initializes each first unrestricted alternative from its fitted null. The global batching schedule groups by tensor shape and cost instead, clears the JAX cache after every work unit, and fits deduplicated alternatives independently from the default initialization. The missing compilation reuse explains much of the human slowdown. The initialization change moves local optima even for scalar-to-scalar routes, one common-converged human scalar-to-scalar test changed by 674.9 statistic units because the hybrid alternative objective was 337.5 units lower.

Selective global batching is therefore rejected as a production replacement. Its mouse speed gain is negligible, its human run is slower, and its discovery and statistic drift require structured-null recalibration. Calibration was not run because the method failed the preceding runtime and numerical-agreement gates. The memory reduction is useful evidence that short cache lifetimes bound resident memory, but production results continue to use the scalar gene-grouped fitter with bounded compilation reuse and null-to-alternative initialization.

## 2026-08-09, follow-up EC GLMM speed options

The scalar Laplace implementation now passes counts, primer-specific EC maps, fixed-effect designs or tensors, cluster membership, selection arrays, and cluster indices as dynamic JAX arguments. A static shape signature permits bounded compilation reuse across genes, and regression tests change mappings, counts, tensors, and cluster assignments while matching uncached objectives and parameters. On the same 20 randomly selected real mouse cell-type tests, gene-local caching took 617.0 seconds and 3.97 GiB peak memory, while a 16-entry global shape cache took 553.4 seconds and 4.41 GiB. The global cache reduced elapsed time by 10.3% and summed objective-evaluation time from 198.5 to 179.4 seconds, with 19 joint convergences and identical LRT statistics in both runs. Exact cross-gene shape reuse is sparse and memory rose by 11.0%, so global caching remains an explicit option rather than the production default.

The posterior-mode solver can reuse the mode from the last accepted L-BFGS iterate and terminate when the mode score is below a selected tolerance. At tolerance $10^{-5}$ this reduced summed optimizer time from 14.3 to 6.43 seconds on the 20-test audit, while elapsed time fell from 617.0 to 589.5 seconds because compilation and initial evaluations dominate. All 19 common fits remained jointly converged, but the median and maximum absolute LRT-statistic changes were $6.1\times10^{-5}$ and 0.319. Tightening the tolerance to $10^{-8}$ reduced these differences to $1.5\times10^{-6}$ and 0.187 but did not improve elapsed time. Accepted-mode warm starts are therefore retained as an experimental option, not the default.

A plug-in raw-EC score test fits the null, maps it into the minimal nested alternative, computes the Laplace outer gradient and observed fixed-effect Hessian, and uses the fixed-effect Schur-complement quadratic while holding variance components at their null estimates. On the same 20 tests it took 903.4 seconds, compared with 617.0 seconds for the full nested LRT, because second-order differentiation and compilation cost more than alternative refinement. Score and LRT statistics had Spearman correlation 0.910, score ranking recovered three of the LRT top five and eight of the top ten, and the score made five nominal calls versus three for the LRT. Seven fixed-effect Hessians were indefinite and the largest condition number was $1.7\times10^{13}$. One block had score statistic 2,911 but LRT statistic 9.88 and LRT p-value 0.873; its minimum Hessian eigenvalue was -0.0132. The raw-EC observed-Hessian score is neither a speed improvement nor sufficiently stable for screening or calibration.

The batched benchmark now projects every fitted null into its corresponding unrestricted alternative before optimization, matching production initialization. For eight moderate mouse fits, the seven tests converged by both batched and scalar procedures had median and maximum absolute statistic differences 0.000303 and 0.0192; the remaining scalar null failed and generated the large unmatched outlier. For 32 cost-sorted fits from the same three-isoform, 17-null-coefficient, 22-alternative-coefficient family, padding inflated the likelihood-work proxy 2.31-fold. CPU batching took 76.2 seconds and GPU batching took 24.3 seconds, a 3.13-fold device speedup, while 29 batched fits converged on both devices. Their median and maximum absolute CPU--GPU statistic differences were below $10^{-11}$ and 0.00254. Larger GPU microbatches are promising, but a production scheduler must preserve gene-local alternative reuse, null-projected initialization, scalar treatment of high-dimensional tails, and structured-null calibration before a full rerun.

## 2026-08-09, production-initialized GPU scheduler

The full-data scheduler now saves every fitted null parameter vector and requires completed null shards before fitting alternatives. Each alternative, whether routed through batched or scalar optimization, begins by projecting its representative fitted null logits and variance parameter into the unrestricted alternative basis. GPU-batched and CPU-scalar work are partitioned into separate routes so slow scalar GPU fits do not erase the vectorization gain. Scalar routes use the bounded global shape cache, while batched routes clear compiled state after each padded work unit. A four-test end-to-end smoke run loaded saved nulls and completed null-projected alternatives. The planned full audit uses batches of 32, minimum batch size four, padding ratio at most three, and at most 128 outer coefficients; this routes 2,534 mouse nulls, 2,333 mouse alternatives, 1,344 human nulls, and 4,871 human alternatives through GPU batches, with all remaining fits on CPU.

## 2026-08-09, full production-initialized GPU audit

An early mouse null run exposed 25 batched-route failures for which scalar production had converged. Failed batch lanes now receive scalar continuation and then a clean restart from the production initialization if necessary. On the rerun, 2,455 of 2,534 batched-route mouse nulls converged, only one scalar-converged null was lost, and 52 nulls that failed scalar production converged. The corresponding human route converged 1,303 of 1,344 nulls. Alternative rescue restarts use the saved fitted-null projection rather than a default alternative vector. Tests cover route construction, saved parameter loading, multilevel label permutation, route cost summaries, and structured-null summaries.

The complete mouse scheduler converged both nested models for 4,597 of 4,825 tests, compared with 4,546 for harmonized scalar production. It made 324 BH calls versus 323, with 319 shared and Jaccard index 0.973. Among 4,544 tests converged by both procedures, median absolute LRT drift was 0.000131, Spearman correlation was 0.996, 281 statistics differed by more than 0.1, 67 by more than 1, four by more than 5, and the maximum was 423.7. Worker time fell from 31.18 to 24.55 hours, 21.3%, while peak memory rose from 3.97 to 5.31 GiB. Nulls used 16.01 hours and alternatives 8.54; GPU-batched portions used only 1.03 and 0.79 hours, while scalar tails used 14.98 and 7.76.

The complete human scheduler converged both nested models for 11,004 of 11,597 tests, compared with 10,981 for harmonized scalar production. It made 1,460 BH calls versus 1,459, with 1,453 shared and Jaccard index 0.991. Among 10,978 common-converged tests, median absolute LRT drift was 0.000191, Spearman correlation was 0.9997, 1,560 statistics differed by more than 0.1, 261 by more than 1, 41 by more than 5, and the maximum was 411.6. Worker time increased from 98.39 to 134.74 hours, 36.9%, and peak memory increased from 6.55 to 8.14 GiB. The 10,253 scalar nulls consumed 85.23 hours, compared with 1.03 hours for 1,344 batch-routed nulls. Alternative GPU batches covered 4,871 of 10,175 unique alternatives in 3.04 hours, but the remaining scalar alternatives required 45.44 hours. Heterogeneous high-dimensional nulls and alternatives therefore dominate human runtime even though the same mathematical fit, tolerances, retries, and routing gates are used in both datasets.

Four full mouse structured-null families independently permuted multilevel cell-type labels within each subject for every gene-row design. They held EC counts, fixed nuisance columns, covered pseudobulks, and each subject's observed label counts fixed. Across 18,379 converged hybrid null tests, rejection rates were 0.0983, 0.0385, and 0.0141 at nominal 0.05, 0.01, and 0.001. Every family made BH calls, 64, 58, 51, and 56. Mean hybrid family cost was 23.17 worker-hours and peak memory was 5.61 GiB.

An exact original gene-grouped scalar rerun used the identical permutation for family 101. It rejected 0.1047, 0.0432, and 0.0148 at the same levels and made 62 BH calls, compared with 64 hybrid calls. The procedures shared 60 calls, with Jaccard index 0.909. Among 4,533 common-converged tests, median statistic drift was $6.24\times10^{-5}$, Spearman correlation was 0.994, 242 statistics differed by more than 0.1, 73 by more than 1, and one local-optimum outlier differed by 422.7. Scalar cost was 31.21 worker-hours and 4.83 GiB, compared with 24.81 hours for hybrid family 101. The multilevel omnibus chi-square LRT is anti-conservative under the design-preserving null regardless of optimizer; batching changes some local optima but does not create the aggregate inflation.

The production-initialized scheduler is rejected as a common production replacement. Its repeatable mouse speedup does not offset the human slowdown, higher memory, and residual local-optimum drift. Gene-grouped scalar fitting remains the shared reported procedure. The ordinary multilevel omnibus BH results are descriptive and not FDR-controlled; this does not alter the separately bootstrap-calibrated primary mouse result or the paired-label pairwise reproducibility benchmark.

## 2026-08-10, analytic EC derivatives and second-order outer optimization

Profiling showed that the latent Laplace Hessian is already block diagonal by subject and each block has only the latent isoform-logit dimension, so sparse factorization has little remaining structure to exploit. The avoidable nested work was instead the automatic differentiation of an automatically differentiated row Hessian inside every outer likelihood gradient. For each primer, the multinomial EC score and Hessian are now evaluated in closed form from the compatibility-weighted abundance matrix, EC masses, total mass, and counts. Finite differences agreed below $10^{-9}$, direct analytic-versus-autodiff Laplace objectives and gradients passed, empty-primer handling is tested, and the original autodiff path remains selectable as a reference.

Matched benchmarks used identical data, candidates, initialization, tolerances, retries, implicit mode gradients, and L-BFGS settings. Across 12 mouse tests, both implementations jointly converged on 11; analytic derivatives reduced end-to-end time from 226.59 to 152.55 seconds, a factor of 1.485, and objective-evaluation time by a factor of 1.239. Across 10 human tests, both jointly converged on nine; elapsed time fell from 283.96 to 192.27 seconds, a factor of 1.477, and objective-evaluation time by a factor of 1.479. Median absolute LRT drift was $2.41\times10^{-11}$ in mouse and $2.27\times10^{-13}$ in human, with maxima 0.0404 and 0.1243. Peak benchmark memory fell from 2.35 to 1.33 GiB in mouse and from 4.81 to 4.32 GiB in human.

Complete analytic-derivative production reruns used 96 shards to reduce repeated compilation, compared with 256 mouse and 1,024 human shards in harmonized autodiff production. Aggregate worker accounting, including result merging, fell from 31.18 to 12.57 hours in mouse, a factor of 2.48, with 4,546 joint convergences under both implementations. Analytic derivatives made 322 BH calls versus 323, all 322 shared, and median absolute LRT drift was $3.64\times10^{-12}$; the maximum local-optimum difference was 6.83. Peak memory was 4.00 versus 3.97 GiB. Human worker time fell from 98.39 to 55.17 hours, a factor of 1.78. Analytic derivatives jointly converged for 10,982 tests versus 10,981, made 1,463 BH calls versus 1,459, and shared 1,458 calls, with Jaccard index 0.996. Median absolute LRT drift was $1.40\times10^{-9}$ and the maximum was 26.68. Human peak memory increased from 6.55 to 8.25 GiB. The matched benchmark isolates the derivative gain, while the larger full-data gain combines that kernel with lower compilation overhead. Closed-form derivatives with gene-grouped L-BFGS replace nested-autodiff derivatives as scalar production without changing the biological conclusions.

An exact Hessian--vector product implementation enabled bounded trust-region Newton--CG outer optimization. A 20-iteration mouse fit took 90.41 seconds versus 40.52 for converged L-BFGS; neither trust nested model converged, and 261 Hessian--vector products consumed 67.28 seconds. The human fit took 151.36 versus 54.34 seconds; neither trust model converged, and 251 Hessian--vector products consumed 125.38 seconds. Larger matched trust jobs remained unfinished after more than three times the analytic L-BFGS benchmark duration and were stopped. Second-order outer optimization is rejected because differentiating through the implicit posterior mode and Laplace log determinant dominates its cost.

## 2026-08-10, retired experimental EC GLMM paths

The codebase was reduced to the supported EC GLMM path after tracing package functions, command-line choices, analysis launchers, and tests. The retained Laplace implementation uses fixed posterior-mode Newton steps, implicit or reference unrolled differentiation, bounded L-BFGS, analytic multinomial derivatives with the autodiff reference available for validation, null-to-alternative initialization, and bounded gene-level or global shape caching. The global cache remains optional because it gave a measured 10.3% elapsed-time reduction with identical statistics, despite its 11.0% memory cost.

Removed speed experiments include projected Adam, Adam followed by L-BFGS, accepted-mode warm starts with adaptive termination, projected alternative and previous-null starts, trust-region outer optimization, independent batched L-BFGS, padded and ragged batching prototypes, and the production-initialized GPU scheduler. These paths either failed convergence, changed local optima or test statistics, slowed the human analysis, or increased memory without a common-dataset runtime benefit. Their benchmark result tables and notebook findings remain as the audit record, while their package modules, runners, launchers, and dedicated tests were retired.

Removed inferential experiments include the raw-EC observed-Hessian score screen, estimate-once clustered Wald test, estimate-once Gaussian mixed-model LRT and score Bayes factors, and the screening policy built from those rankings. The score screen was slower than the nested LRT and had indefinite Hessians. The path Wald reference was severely anti-conservative under design-preserving nulls, while the Gaussian LRT and fixed Bayes-factor thresholds retained deep-tail inflation and required expensive empirical calibration. Local path estimation and isoform-to-path collapse remain because they support the calibrated path analyses and the retained latent-space ablation; only the rejected second-stage tests were removed.

Bouchard CAVI and CAVI-initialized tilted fitting were also retired from the package and runners. The derivation remains in `docs/cavi.tex`, and the high-level manuscript now labels it as a negative result. CAVI increased total time in simulation, reduced real-data joint convergence from 80.0% and 86.2% under the established condition and cell-type initializations to 16.9% and 15.4%, and rescued no historical failures.

## 2026-08-10, analytic-derivative condition production reruns

The ordinary multinomial condition-within-cell-type families were rerun with the same closed-form row derivatives, gene-grouped scalar L-BFGS, 12-step implicit posterior-mode solve, 30-step continuation, initialization, bounds, and convergence tolerances used for the completed analytic cell-type reruns. Mouse used 96 shards. A 96-shard human attempt placed about 762 heterogeneous tests in each process and aborted as compiled shapes accumulated beyond the worker memory budget. The successful human campaign used 768 shards, about 95 tests each, without changing a fit or inference setting. The aborted coarse campaign is excluded from the production throughput comparison because it produced no reportable complete analysis; its failure establishes the viable sharding bound.

Mouse condition fitting completed all 14,833 tests. Analytic and nested-autodiff production each jointly converged on 13,779 tests and made the same three BH calls. The median absolute LRT difference among common converged tests was $2.27\times10^{-13}$ and the maximum was 0.229. Worker time including merging fell from 87.60 to 38.41 hours, a factor of 2.28, while peak worker memory changed from 4.74 to 4.81 GiB.

Human condition fitting completed all 73,144 tests. Analytic derivatives jointly converged on 68,774 tests versus 68,771 under nested autodiff, with 68,767 common converged tests. Both runs made the same 193 BH calls. The median absolute LRT difference was $9.09\times10^{-13}$ and the maximum was 3.03. Worker time including merging fell from 631.05 to 258.10 hours, a factor of 2.44, and peak memory fell from 6.83 to 6.16 GiB.

Combining the updated cell-type and condition families, broad Tealeaf worker time is 50.97 hours for mouse and 313.27 hours for human, versus 118.79 and 729.44 hours before closed-form derivatives. Peak broad worker memory is 4.81 GiB for mouse and 8.25 GiB for human; the latter peak comes from cell-type fitting. Relative to the matched scQuint campaigns, Tealeaf is now 3.34-fold slower in mouse and 3.33-fold slower in human, rather than 7.8-fold slower. Condition convergence, discovery counts, and the manuscript's biological interpretation are unchanged.

## 2026-08-11, cell-type event structures and functional literature

The 238 mouse cell-type partitions that pass BH under all four GPD thresholds were annotated from the supported path signatures used by each fitted test, rather than every GTF path. Their mutually exclusive structural classes are 154 alternative-first regions, 31 alternative-last regions, 28 cassette exons, 22 whole-gene complex blocks, two compound internal blocks, and one multi-exon-skipping block. Terminal categories are called regions because a pseudo-anchor block can combine promoter, terminal-exon, polyadenylation, and internal differences. The reusable classifier removes exact intervals shared by two tested paths before recognizing canonical differences, which correctly identifies the App and Grin1 discoveries as cassette events even though their paths share other exons. The complete catalog records coordinates, tested signatures, nested canonical components, empirical p-values, and FDR values. The audit is restricted to the calibrated mouse set because ordinary multilevel LRT yields in mouse and human are descriptive under the structured-null results.

Coordinate-level literature review found five interpretable examples. Ntrk2 alternative last-region paths include the kinase-deficient TrkB.T1 family, whose astrocyte-specific loss impairs morphogenesis and support of excitatory synapses. The Grin1 event is the 111-nt C1 cassette, whose retention motifs regulate GluN1 surface trafficking. The App event is the 54-nt mouse exon 15 cassette; its exclusion produces L-APP, changes polarized secretion and amyloidogenic processing, and has a cell-specific neuronal versus non-neuronal distribution. The Enah event is the exact 63-nt exon 11a cassette, with characterized effects on actin polymerization and cell motility outside the brain. The Dlg4 N-terminal block includes the PSD-95 alpha/beta system, but its six-path omnibus does not isolate that pair. These examples establish functional precedents without assigning the omnibus statistic to one cell-type pair or effect direction.

## 2026-08-24, differential-splicing presentation

Added `docs/differential_slides.tex`, a self-contained 16:9 Beamer presentation of the EC differential-splicing analysis. The deck follows the forward model from local splice-path effects through isoform composition and primer-specific EC probabilities to observed pseudobulk counts. It summarizes the GSE189033 mouse and GSE233208 human sources, distinguishes the calibrated primary mouse analysis from descriptive ordinary-LRT sensitivity results, compares matched Tealeaf, LeafCutter, MAJIQ Heterogen, and scQuint yields and resources, and shows the structural classes and selected functional precedents among the 238 robust mouse cell-type calls. All diagrams and charts are vector-native TikZ objects in the TeX source.

Added Figure `fig:ec-generative-model` to `docs/differential.tex`. The left-to-right diagram shows block-level fixed and latent effects producing a shared isoform composition, which branches through separate oligo(dT) and random-hexamer EC maps into primer-specific EC proportions and observed multinomial count vectors. The caption states that the branch represents observation chemistry rather than independently fitted isoform compositions.

## 2026-08-28, differential manuscript cleanup

Reorganized `docs/differential.tex` so all production and comparator methods precede results. The main Results now begin with the independent-subject split reproducibility comparison on shared gene and cell-type-pair universes, followed by the full calibrated Tealeaf mouse results, structural event classes, and functional literature precedents. Raw genome-wide discovery-yield tables for Tealeaf, LeafCutter, MAJIQ, and scQuint were removed because their feature universes and filters differ and the counts do not measure comparative power. Detailed permutation calibration for the split benchmark, nonbootstrap sensitivities, computational benchmarks, alternative models, and supporting algorithms now appear only after the appendix boundary.

Defined the generalized Pareto distribution (GPD) at first use as the peaks-over-threshold tail model, replaced visible “EC block” terminology with “block” or “block-level EC-count model,” and removed the deprecated Bouchard/CAVI production discussion from this manuscript. The separate `docs/cavi.tex` derivation remains the provenance record. Recompilation completed without missing or duplicate references; no new statistical analysis was required because the primary comparison uses the existing completed subject splits and paired-label permutation families.

## 2026-08-28, separate EC fitting from the paired junction baseline

Revised `docs/differential.tex` to distinguish the two supported fitting options for the same EC GLMM. The main methods now describe both full-covariance tilted VI, which supplies the primary bootstrap-calibrated result, and analytic Laplace fitting, which supplies conventional approximate marginal likelihood, LRT, and BIC summaries. The text records that tilted VI was faster than the original Laplace implementation in controlled simulations, while closed-form EC derivatives and compilation reuse made the later Laplace implementation practical for large nonbootstrap screens. It also states that the ordinary multilevel Laplace chi-square reference is anti-conservative under structured mouse null permutations and therefore does not supersede the calibrated VI result.

Moved the paired junction centered-log-ratio method, its subject-split reproducibility table, and its permutation calibration out of the main methods and results into the exploratory appendix. Relabeled it throughout as the in-house paired CLR baseline rather than Tealeaf. The appendix now states that this test was introduced to match the pairwise hypotheses and filtered junction representation used by the external comparators, and that its reproducibility advantage cannot be attributed to EC ambiguity resolution, primer-specific compatibility, splice blocks, or latent EC fitting. The main Results now begin with the full calibrated EC-GLMM analysis. Recompiled `docs/differential.pdf`, 49 pages, without missing references, duplicate labels, or overfull boxes.

## 2026-08-28, restore the comparator result using the actual Tealeaf model

The preceding edit correctly relabeled the paired junction CLR method but incorrectly removed the requested main-text comparator result. Rescored the completed pairwise logistic-normal EC-GLMM subject-half fits against LeafCutter, relaxed-coverage MAJIQ Heterogen, and scQuint. This is the actual Tealeaf likelihood, retaining paired-primer EC counts, primer maps, isoform ambiguity, condition nuisance effects, mouse random intercepts, and pseudobulk logistic-normal effects; only the inference option is Laplace rather than tilted VI. Each comparison is restricted to genes and cell-type pairs eligible for both methods in both condition-balanced subject halves.

The primary split benchmark uses the calibration-selected minimum median gene depth of 100 modeled UMIs and at least four paired mice. In the corresponding label-swap subset, 77 converged tests reject at rate 0.0519 at nominal 0.05, none has p-value below 0.01, and none passes BH. Against LeafCutter, the shared universe has 36 genes and Tealeaf has 4 replicated genes versus 9, held-half replication 0.757 versus 0.766, and fold rank correlation 0.397 versus 0.661. Against relaxed MAJIQ, the shared universe has 63 genes and the corresponding values are 11 versus 17, 0.782 versus 0.823, and 0.494 versus 0.700. Against scQuint, 61 genes are shared and the values are 11 versus 11, 0.782 versus 0.801, and 0.457 versus 0.616. Thus Tealeaf ties scQuint in replicated discoveries but does not outperform the comparators overall.

An eight-paired-mouse sensitivity retains 30, 54, and 52 shared genes. Tealeaf has 3 versus 8 replicated genes against LeafCutter, 9 versus 14 against relaxed MAJIQ, and 9 versus 11 against scQuint. Added both scoring outputs under `analyses/ec_glmm_split_reproducibility_min100_min4` and `analyses/ec_glmm_split_reproducibility_min100_min8`. The main Methods now define the actual EC-GLMM split benchmark, and the first main Results subsection reports its comparator table. The paired CLR method and its separate result remain appendix-only and are not labeled Tealeaf.

## 2026-08-28, scalar production EC-GLMM CPU--GPU benchmark

Benchmarked the unchanged gene-grouped scalar Laplace EC-GLMM on CPU and GPU within one interactive GPU allocation. Both routes used the same JAX 0.9.2 environment with 64-bit arithmetic, data, candidate order, analytic multinomial derivatives, 12 posterior-mode Newton updates, implicit gradients, L-BFGS settings, and retries. The matched host supplied two Intel Xeon Gold 5115 CPU cores and one 16 GB NVIDIA V100. Timing excludes input loading and candidate construction and includes objective compilation and fitting.

Across 12 random mouse tests, CPU required 128.15 seconds and GPU 187.16 seconds, so GPU was 1.46-fold slower. Across nine gene-grouped dimension-stratified mouse tests, CPU and GPU required 90.59 and 90.22 seconds. Across eight grouped human tests they required 95.65 and 94.26 seconds, a 1.5% GPU gain. A production-like human sensitivity with one block from each of four genes required 70.31 seconds on CPU and 71.33 seconds on GPU, a 1.4% GPU slowdown. Every nested fit converged on both devices. Median absolute LRT differences were at most $3.0 \times 10^{-11}$, and the maximum was 0.0084 in the random mouse sample.

Direct GPU execution therefore does not materially accelerate the scalar production fitter. Some grouped human objectives run faster on device, but per-shape compilation, host-side SciPy control, synchronization, and small matrices offset the kernel gain. Shape-changing workloads can be substantially slower on GPU. Retain CPU as the production device. The existing batched GPU prototype addresses kernel occupancy but remains unsuitable as the common production scheduler because the full human audit was slower and showed residual local-optimum drift. Saved the device comparison under `analyses/ec_glmm_gpu_device_benchmark`.

## 2026-08-28, restore the full EC-GLMM comparator test universe

Audited the unexpectedly small comparator table and found that the 100-UMI sensitivity filter, not model fitting or identifier matching, caused the collapse. The completed pairwise logistic-normal EC GLMM contains 7,671 and 5,140 block-by-cell-type-pair fits in the two subject halves, with 7,436 and 4,949 jointly converged fits. The 100-UMI filter reduced these to 1,101 tests from 202 genes and 564 tests from 94 genes before comparator intersection. This sensitivity cannot serve as the main genome-wide benchmark.

Rescored the native-screened converged fits, requiring at least four paired mice but no post hoc depth filter. The shared LeafCutter universe contains 739 gene-pair tests from 138 genes; Tealeaf has 10 replicated genes versus 25, held-half nominal replication 0.839 versus 0.744, and fold rank correlation 0.561 versus 0.537. The relaxed MAJIQ universe contains 2,471 tests from 475 genes; the corresponding values are 17 versus 38, 0.849 versus 0.692, and 0.394 versus 0.432. The scQuint universe contains 2,519 tests from 494 genes; the values are 17 versus 57, 0.830 versus 0.773, and 0.366 versus 0.464.

In the 1,000-test paired-label-swap audit, 967 fits jointly converge and reject at rates 0.0745, 0.0186, and 0.00414 at nominal levels 0.05, 0.01, and 0.001, with no BH discoveries. The native asymptotic p-values are mildly anti-conservative, so the main table is explicitly an empirical split-reproducibility comparison and its BH counts are descriptive rather than FDR-controlled power estimates. Moved the 100-UMI result to an appendix sensitivity and saved the full-universe scoring under `analyses/ec_glmm_split_reproducibility_full_min4`.

## 2026-08-29, gene-level EC-GLMM aggregation and low-depth audit

Audited the restored split benchmark and found that its cached manifests use a 25-modeled-gene-UMI inclusion threshold, although the manuscript and analysis README called it 10. Corrected the documentation to 25. The two production block manifests contain 7,671 and 5,140 block-by-cell-type-pair tests from 935 and 750 genes, with 7,436 and 4,949 jointly converged fits. This correction changes no reported statistic.

Tested an unrestricted joint-gene LRT at the 25-UMI screen. The two manifests contain 5,710 and 3,930 gene-by-cell-type-pair tests from 935 and 750 genes. Against LeafCutter, relaxed MAJIQ, and scQuint, Tealeaf has 14, 30, and 27 replicated genes compared with 24, 38, and 57, with held-half nominal replication 0.817, 0.870, and 0.856. In a 1,000-test paired-label-swap audit, 967 fits jointly converge and reject at rates 0.0931, 0.0310, and 0.00931 at nominal levels 0.05, 0.01, and 0.001. The joint alternative recovers discoveries but worsens calibration.

Implemented and fully evaluated a gene model in which annotated blocks contribute independent fixed-effect subspaces. The model forms the union of all block path-contrast bases for a gene and design, removes linear dependencies by singular-value decomposition, and retains tested-label effects orthogonal to that union in the null. It reduces the tested rank relative to the unrestricted joint model for 1,864 of 5,710 fold-0 tests and 1,445 of 3,930 fold-1 tests, usually by one degree of freedom. Its replicated-gene counts are 12 versus 24 for LeafCutter, 19 versus 38 for relaxed MAJIQ, and 18 versus 57 for scQuint, with held-half replication 0.844, 0.843, and 0.842. Among 966 converged null fits, rejection rates are 0.0704, 0.0217, and 0.00414, and two tests survive BH. The independent-block restriction loses annotation-orthogonal signal without improving calibration over production, so its command-line option, tensors, and dedicated tests were removed after the audit.

Lowered the screen to 10 modeled gene UMIs and reran the unrestricted joint-gene model on full data splits. The manifests expand to 22,453 tests from 2,503 genes and 16,753 tests from 2,238 genes. The lower threshold increases the number of included pseudobulks from a median 12 to 24 among tests shared with the 25-UMI run. Against LeafCutter, Tealeaf has 23 replicated genes versus 29, held-half replication 0.780 versus 0.684, and fold rank correlation 0.442 versus 0.433 on 1,302 shared gene-pair tests from 208 genes. Against relaxed MAJIQ, the corresponding values are 66 versus 35, 0.839 versus 0.659, and 0.402 versus 0.354 on 8,613 tests from 1,243 genes. Against scQuint they are 53 versus 104, 0.842 versus 0.776, and 0.369 versus 0.401 on 9,031 tests from 1,319 genes. The low-depth model beats MAJIQ in replicated discoveries and all comparators in held-half replication, but only 946 of 1,000 paired-label-swap fits jointly converge; 86, 21, and 6 fall below 0.05, 0.01, and 0.001, and three survive BH. The 10-UMI joint model is therefore retained only as a documented sensitivity. The production blockwise 25-UMI model remains the least anti-conservative broad EC-GLMM split configuration tested, although its asymptotic p-values are still descriptive and it does not outperform the junction comparators in replicated discovery count.

## 2026-08-29, local splice-path power diagnosis and corrected EC-GLMM null

The previous blockwise isoform test did not implement the intended null of constant aggregate local-path usage. For a block with $T\in\mathbb N$ supported gene isoforms and $S\in\mathbb N$ represented local paths, the old construction lifted each path contrast into the $T-1$ dimensional reference-logit space and treated its Euclidean orthogonal complement as tested-label nuisance variation. This is incorrect when isoforms assigned to one path have unequal baseline abundance. If $I_s\subseteq\{1,\ldots,T\}$ is the isoform set for path $s\in\{1,\ldots,S\}$, $\boldsymbol\theta\in\Delta^{T-1}$ is the gene isoform composition, and $\psi_s=\sum_{t\in I_s}\theta_t$, then the derivative of a local path log-ratio is weighted by $\theta_t/\psi_s$. Orthogonality to an unweighted path indicator does not imply zero derivative of $\log(\psi_s/\psi_{s'})$. The old null could therefore absorb the path effect it was meant to test.

This leakage is common in the real cassette-exon benchmark. Fold 0 has 1,081 of 1,353 cassette tests with nonzero leakage, with median nuisance-to-tested derivative norm 0.341, 90th percentile 0.522, and maximum 0.817. Fold 1 has 653 of 859 affected tests, with median 0.330, 90th percentile 0.523, and maximum 0.816. In fold 0, leakage has Spearman correlation 0.236 with the corrected model's gain in $-\log_{10}p$, and the highest-leakage quartile has median gain 0.213.

The corrected global-isoform tangent test first estimates a label-independent pooled baseline $\bar{\boldsymbol\theta}\in\Delta^{T-1}$. Let $\boldsymbol J\in\mathbb R^{(S-1)\times(T-1)}$ be the Jacobian of the block path ILR with respect to the gene reference logits at that baseline. The tested basis is the right inverse $\boldsymbol H=\boldsymbol J^{\mathsf T}(\boldsymbol J\boldsymbol J^{\mathsf T})^{-1}\in\mathbb R^{(T-1)\times(S-1)}$, so $\boldsymbol J\boldsymbol H=\boldsymbol I_{S-1}$. The tested-label nuisance basis $\boldsymbol N\in\mathbb R^{(T-1)\times(T-S)}$ spans $\ker(\boldsymbol J)$, so $\boldsymbol J\boldsymbol N=\boldsymbol0$. The null includes the tested-label design multiplied by $\boldsymbol N$, while the alternative adds the same design multiplied by $\boldsymbol H$. A pairwise two-path cassette therefore tests one coefficient. The likelihood still uses raw primer-specific EC count vectors, both compatibility maps, every supported isoform, multinomial covariance, mouse random effects, and pseudobulk logistic-normal effects.

A second implementation audit found that the unrestricted alternative cached for tests sharing a gene and covered pseudobulk set could be a poorer local optimum than a later block null. This violated numerical nesting in about 18% of original cassette fits. Reused alternatives are now checked against each fitted null, refined from the null projection when unconverged or worse, and replaced only when the refined objective improves. This repairs nesting, although it does not by itself explain the cassette power gain.

An exact local-path sensitivity collapses compatibility columns within each represented path using pooled, label-independent within-path isoform weights and retains any background isoform separately. It fits the raw EC counts directly rather than estimating global isoforms and then running a second-stage test. Among 732 cassette tests converged in both subject halves, 27 are nominal in both halves and six pass BH on the conjunction p-value, compared with 11 and two under the former global model on its common universe. All six conjunction calls agree in effect sign, 25 of 27 both-nominal calls agree, and their effect-size Spearman correlation is 0.728. Among 971 converged cassette paired-label swaps, rejection rates are 0.0587, 0.0113, and 0.00206 at nominal levels 0.05, 0.01, and 0.001, with no BH calls.

The tangent test preserves within-path isoform nuisance variation and nearly duplicates the exact local-path result. It has 1,319 and 824 jointly converged cassette tests in the two halves. Among 730 tests converged in both halves, 28 are nominal in both and six pass conjunction BH; 26 of 28 nominal effects and all six BH effects agree in sign, with effect-size Spearman correlation 0.722. Its 975 converged paired-label swaps reject at rates 0.0615, 0.0133, and 0.00205, with no BH calls. The tangent model is therefore the production correction, while the exact fixed-mixture local-path model remains a sensitivity.

Exonic reads are useful only when their compatibility patterns differ across paths. In a 100-cassette audit, 24 events lack an EC unique to at least one path. The 76 feasible events have median 96.5 path-unique UMIs, only 10.7% of crude gene UMI exposure, and a pure-EC test makes fewer nominal calls than the full local EC likelihood. This rejects simple counting of path-unique ECs. The power advantage comes from carrying all EC observations through the compatibility likelihood, not from treating every exonic read as path-identifying.

## 2026-08-29, full tangent-null split rerun

Completed the production tangent-null EC GLMM over the full native-screened mouse split manifests. Fold 0 jointly converged for 7,471 of 7,671 block-by-cell-type-pair tests, made 1,014 nominal and 316 fold-specific BH calls, and used 18.06 aggregate worker hours. Fold 1 jointly converged for 4,967 of 5,140 tests, made 558 nominal and 136 fold-specific BH calls, and used 13.06 aggregate worker hours. The corresponding former-basis split fits used 51.79 and 40.43 worker hours, so analytic likelihood derivatives and safe unrestricted-alternative reuse reduce elapsed worker cost 2.87-fold and 3.10-fold.

Scored the corrected fits against the same independent-subject comparator outputs. Against LeafCutter, 742 shared gene-pair tests from 139 genes give 12 Tealeaf versus 25 comparator replicated genes, held-half nominal replication 0.867 versus 0.744, and fold rank correlation 0.562 versus 0.543. Against MAJIQ Heterogen, 2,481 tests from 473 genes give 19 versus 38 replicated genes, held-half replication 0.828 versus 0.697, and correlation 0.362 versus 0.434. Against scQuint, 2,531 tests from 492 genes give 18 versus 57 replicated genes, held-half replication 0.828 versus 0.777, and correlation 0.358 versus 0.461.

The corrected null raises Tealeaf's replicated-gene counts from 10, 17, and 17 to 12, 19, and 18 and improves the cassette-specific result more clearly, but Tealeaf still does not exceed comparator replicated discovery counts genome-wide. Its higher held-half nominal replication indicates a smaller, more stable discovery set. The expected one-degree-of-freedom advantage is exact for a two-path cassette, not for every compound block or every external feature after gene aggregation. Shared exonic reads also do not provide direct path-composition information when their compatibility is equal across paths. The full-data 238-call bootstrap result predates the tangent correction and is now labeled provisional in the manuscript rather than being represented as a corrected-model result.

## 2026-08-30, direct local-path paired test and moderated inference

Completed the corrected tangent-null full-data bootstrap reruns. The cell-type family contains 4,825 tests and 522 nominal p-values below 0.05, while the condition family contains 14,833 tests and 922 nominal p-values below 0.05; neither family has a BH discovery after cross-fitted generalized-Pareto tail scoring and leave-block-out empirical recalibration. The calibrated null distributions are exact at the audited nominal levels, but the smallest resolved p-values, 0.00216 for cell type and 0.000148 for condition, do not reach their genome-wide BH thresholds. The correction therefore removes the former-basis discovery set rather than confirming it.

Implemented a direct local-path estimator for the pairwise subject-split benchmark. For each subject and cell-type level, it aggregates the two primer-specific EC vectors separately, perturbs the represented path logits around a label-independent pooled isoform baseline, and evaluates the raw EC multinomial likelihood through both primer maps. Within-path isoform ratios, the total mass assigned to represented paths, and unrepresented nuisance isoforms remain fixed during each local fit. A configurable symmetric path pseudocount prevents boundary estimates; the selected setting is two. Subject-level ILR differences are tested by a paired t-test for two paths and Hotelling's test for more paths. The procedure therefore retains ambiguous junction and exonic reads through the compatibility likelihood but avoids repeatedly fitting a global latent GLMM for a pairwise contrast.

The unregularized paired-path estimator had well-behaved pooled empirical calibration but unstable analytic tails and only 1, 0, and 0 replicated genes against LeafCutter, MAJIQ, and scQuint. A 0.5 path pseudocount restored analytic null calibration and raised the replicated counts to 8, 20, and 21, still below 26, 38, and 56 for the comparators. Dropping nuisance isoforms and nuisance-only ECs, combining global and locally conditioned p-values by Cauchy aggregation, restricting to two-path blocks, and alternative gene-level Cauchy aggregation did not improve calibrated replication and are rejected.

With the conventional add-one path pseudocount, the native paired tests made 502 and 332 fold-specific BH calls. Their pooled sign-flip null rejection rates were 0.0488, 0.0104, and 0.00124 in fold 0 and 0.0472, 0.00983, and 0.00132 in fold 1 at nominal 0.05, 0.01, and 0.001; no null family made more than one BH call. On matched universes, native replicated-gene counts were 24 versus 26 for LeafCutter, 46 versus 38 for MAJIQ, and 42 versus 56 for scQuint. Held-half replication was higher for Tealeaf in all three comparisons.

Most two-path tests have few paired subjects, so their separate sample variances are noisy. Added robust empirical-Bayes variance moderation for the scalar paired tests. The fitted scaled inverse-chi-squared prior uses log-variance moment matching after removal of each test's chi-squared sampling moments and 5% symmetric Winsorization. Each sign-flip null family refits its own variance prior before its moderated p-values are computed. A leave-one-event-out empirical null CDF within path dimension and subject-count strata then corrects residual tail error. This yields 803 and 510 fold-specific BH calls, exact pooled null rates of 0.0500, 0.0100, and approximately 0.0010, and zero BH calls in every one of the 64 fold-specific null families.

The strictly calibrated moderated add-one result has 31 replicated genes versus 26 for LeafCutter, 53 versus 38 for MAJIQ, and 56 versus 56 for scQuint. Tealeaf held-half replication rates are 0.865, 0.808, and 0.799, compared with 0.750, 0.698, and 0.775, and its fold rank correlations are 0.687, 0.550, and 0.574, compared with 0.542, 0.437, and 0.458. Propagating all sign-flip p-values through Simes aggregation, fold conjunction, and BH produces zero replicated-gene calls in every one of 32 null families on each matched comparator universe.

The add-two sensitivity is also exactly calibrated after empirical correction. It makes 1,163 and 770 fold-specific BH calls, with pooled null rejection rates 0.05004, 0.01000, and 0.00100 in fold 0 and 0.05001, 0.01003, and 0.00100 in fold 1; no fold-specific null family makes a BH call. On the matched intersections, Tealeaf has 50 replicated genes versus 26 for LeafCutter, 103 versus 38 for MAJIQ, and 101 versus 56 for scQuint. Its held-half replication rates are 0.799, 0.810, and 0.819, compared with 0.750, 0.698, and 0.775, and fold rank correlations are 0.711, 0.622, and 0.656, compared with 0.542, 0.437, and 0.458. The full Simes-plus-conjunction null audit again yields zero replicated-gene BH calls in all 32 null families on every comparator-matched universe. Add-two smoothing is selected for the pairwise production benchmark because it exceeds all three comparators in both independently replicated discoveries and held-half replication while retaining calibration at the reported endpoint. Because the smoothing strength was selected after examining these splits, the comparison is internally tuned and requires confirmation in an independent dataset.

The selected two-fold path fitting run used 2.516 aggregate worker hours across 192 one-core shards, with maximum shard elapsed time 78 seconds and peak worker resident memory 0.288 GiB. Merge-time variance moderation and empirical calibration add less than two minutes of serial work. Direct local-path fitting is therefore substantially cheaper than the repeated nested EC GLMM while retaining the paired-primer compatibility likelihood.

## 2026-08-30, empirical-Bayes smoothing and conditional path uncertainty

Implemented baseline-centered smoothing for the direct local-path likelihood. For block baseline path composition $\bar{\boldsymbol\psi}_b\in\Delta^{S_b-1}$ and scalar path pseudocount $\alpha_b>0$, the Dirichlet prior counts are $\boldsymbol a_b=\alpha_bS_b\bar{\boldsymbol\psi}_b\in\mathbb R_+^{S_b}$. The pooled isoform baseline is fitted without differential labels. A deterministic candidate sample is subject-cross-fitted, its held-out EC counts are scored by a Laplace marginal likelihood, and hyperparameters are cross-fitted over five gene folds. A one-component EB model was driven to the maximum precision by the null majority. A spike-and-slab evidence model with a spike at $\alpha=128$ selected finite slab values from 0.25 through 8; split and full-data estimates were usually 4 or 8, with null-spike weights from 0.91 to nearly 1.

Added conditional Laplace covariance propagation. If $I_{bmc}\in\mathbb S_+^{S_b-1}$ is the expected EC information and $C(\boldsymbol\psi)=\operatorname{diag}(\boldsymbol\psi)-\boldsymbol\psi\boldsymbol\psi^{\mathsf T}$, then $V_{bmc}=\{I_{bmc}+\alpha_bS_bB_b^{\mathsf T}C(\widehat{\boldsymbol\psi}_{bmc})B_b\}^{-1}\in\mathbb S_+^{S_b-1}$. Paired differences use $W_{bm}=V_{bm1}+V_{bm2}$ in a Gaussian measurement-error model with isotropic biological variance fitted by REML. Independent condition tests use the corresponding scaled heteroskedastic multivariate GLS. The shared covariance scale is selected outside each held-out gene fold by REML profile evidence.

Neither extension improved the split discovery endpoint. Unit-scale conditional covariance produced 684 and 524 fold-specific BH calls but only 28, 10, and 11 independently replicated genes against LeafCutter, MAJIQ Heterogen, and scQuint. Cross-fitted scale profiles over 0 through 3 selected zero in every split, full cell-type, and full condition stratum. The predictive spike-and-slab EB alpha reduced fold-specific BH calls to 106 and 42 and independently replicated genes to 3, 0, and 0, while held-half nominal replication remained 0.885, 0.889, and 0.783. Predictive evidence is dominated by the large null mass and its finite slab still over-shrinks alternatives needed for discovery. Both mechanisms remain implemented and tested as sensitivities. Production retains uniform add-two smoothing, robust scalar variance moderation for two-path paired tests, and design-matched empirical null calibration.

Replaced the deprecated full-data EC GLMM with the production local-path analysis. The cell-type family contains 11,462 pairwise block hypotheses, of which 11,315 are eligible, 3,537 have calibrated p-values below 0.05, and 2,222 pass BH at 0.05. These calls represent 585 distinct blocks and 440 genes. The pooled sign-flip null rejection rates are 0.050008, 0.010025, and 0.000991 at nominal levels 0.05, 0.01, and 0.001; none of 32 null families has a BH call. The condition-within-cell-type family contains 14,833 eligible hypotheses, 1,111 nominal events, and no BH event. Its pooled label-permutation null rates are 0.049962, 0.009952, and 0.000946, with zero BH calls in every null family.

Classified the 585 significant cell-type blocks from their supported fitted paths. The mutually exclusive classes are 347 alternative-first regions, 110 cassette exons, 66 alternative-last regions, 24 whole-gene complex blocks, 16 compound internal blocks, and 22 other internal structures. Four literature-audited structures recur in the production set: Ntrk2 in 12 cell-type pairs, Grin1 in five, Dlg4 in seven, and App in 17. Enah does not pass production BH. The catalog and summaries are under `analyses/local_path_full_event_annotation`.

The full cell-type fit used 1.0348 aggregate worker hours, a maximum shard fitting time of 151.6 seconds, and 0.384 GiB peak worker memory. The condition null design is intercept-only and therefore invariant to label permutations. Reusing its fitted REML biological variance instead of reoptimizing the identical null for every permutation is exact and reduced the optimized full condition run to 2.2459 aggregate worker hours, a maximum shard fitting time of 283.2 seconds, and 0.410 GiB peak worker memory. A regression test verifies equality with refitting after label permutation. The former global EC GLMM method and results now appear only in the appendix of `docs/differential.tex`.

## 2026-08-31, mouse-blocked omnibus cell-type regression

Implemented a block-level multivariate cell-type test as an alternative to separate cell-type-pair tests. EC counts are first summed for each mouse--cell-type combination and quantified with the same uniform add-two local-path likelihood. For a block with $N$ observed mouse--cell-type combinations, $M$ mice, $C$ cell types, and $d=S-1$ path ILRs, the response is $Y\in\mathbb R^{N\times d}$. Its design contains an intercept, $M-1$ mouse indicators, and $C-1$ cell-type indicators. The joint Wald test of the $(C-1)\times d$ cell-type coefficient block retains within-mouse pairing. Thirty-two null families permute cell-type labels within mouse and use the same leave-one-event-out empirical calibration as production.

The full omnibus manifest contains 1,785 candidate blocks. Four become rank deficient after local-fit filtering, leaving 1,781 fitted tests. There are 728 calibrated p-values below 0.05, 584 BH blocks, and 434 BH genes. The pooled null rejection rates are 0.050025, 0.009966, and 0.000965 at nominal levels 0.05, 0.01, and 0.001, and none of 32 null families has a BH call. The omnibus and event-level pairwise BH sets contain 584 and 585 blocks, with 492 blocks shared.

Compared the methods after aggregating pairwise p-values by Simes to one p-value per block. On the 1,624 full-data blocks shared by both methods, omnibus and pairwise-Simes testing call 599 and 437 blocks, with 412 shared. On 988 blocks shared by both methods and both subject halves, omnibus testing independently replicates 212 blocks and 159 genes, compared with 112 blocks and 98 genes for pairwise-Simes testing. Held-half nominal replication is 0.722 versus 0.687, and split rank correlation is 0.613 versus 0.620. Both block-level and gene-level conjunction pipelines produce zero BH calls in all 32 null families. The omnibus test is preferred for the global question of whether any cell type differs; pairwise tests remain necessary to localize the differing pairs and to match the external comparator hypotheses.

When conditional path covariance is zero, scalar-variance multivariate GLS reduces to ordinary least squares. Added a closed-form path that estimates the null variance with the same lower bound as the numerical solver and reuses it for label permutations. The final full omnibus run uses 1.181 aggregate one-core worker hours, compared with 1.445 hours before this optimization and 1.035 hours for production pairwise testing. The closed form leaves the omnibus BH set unchanged. The omnibus model remains slightly slower because it quantifies more cell types per block and fits higher-dimensional designs despite testing fewer hypotheses.

## 2026-08-31, differential-method documentation cleanup

Reorganized the non-appendix manuscript around the current Tealeaf method. Methods now contains a Tealeaf subsection divided into local-path estimation and statistical testing, followed by a separate comparator subsection. The subject-blocked multivariate omnibus regression is the primary cell-type test. The paired cell-type test is secondary, providing pair-specific localization and the hypothesis geometry required for comparisons with LeafCutter, MAJIQ, and scQuint. All non-appendix method descriptions use subject rather than dataset-specific mouse terminology.

Removed appendix-only symbols from the main notation table, including the conditional covariance, EC-probability vector, path Jacobian, right inverse, and tangent nuisance basis. The primer-specific EC probability vector is now defined when the deprecated global EC GLMM first uses it in the appendix. The main full-data table now reports the primary omnibus cell-type family, 1,785 candidates, 1,781 eligible blocks, 728 nominal blocks, 584 BH blocks from 434 genes, minimum calibrated p-value $9.53\times10^{-5}$, pooled null rejection 0.05002, and zero null-family BH calls. The 2,222 pairwise calls across 585 blocks and their structural and literature annotations are explicitly labeled as secondary localization results.

Rebuilt `docs/differential.pdf` with two clean LaTeX passes. There are no undefined references or fatal warnings; the remaining underfull-box notices are confined to narrow cells in the literature table.

## 2026-08-31, Monte Carlo null-replicate sensitivity

Defined $R_{\rm null}\in\mathbb N$ as the Monte Carlo null-replicate budget, with default $R_{\rm null}=32$. For block $b$ in calibration stratum $\mathcal S_b$, calibration pools the $R_{\rm null}$ null reference p-values from every other block in the stratum, giving approximately $R_{\rm null}|\mathcal S_b|$ calibration values rather than only $R_{\rm null}$. In contrast, the complete-family BH audit has only $R_{\rm null}$ realizations. The manuscript now gives the leave-one-block-out empirical-CDF formula and states that 32 is a computational default rather than a theoretically distinguished value.

Reran the primary omnibus, condition, full pairwise, and both pairwise subject-split analyses with 128 deterministic null replicate indices. Added `--max-null-replicates` to the paired-path merger and recalibrated the same outputs using nested prefixes of 32, 64, and 128 replicates. Across all five analyses, calibrated $-\log_{10}p$ ranks at 64 and 128 correlate at least 0.9994 with the 32-replicate values, median absolute p-value changes are at most 0.00221, and BH-set Jaccard similarity is 0.972--1.000. Omnibus BH counts are 583, 586, and 585; condition counts are zero at every budget; full pairwise counts are 2,225, 2,226, and 2,228. The independent production omnibus run calls 584 blocks, within this Monte Carlo variation.

Repeated the matched comparator scoring. Tealeaf replicated-gene counts at 32, 64, and 128 are 50, 50, and 50 against LeafCutter; 103, 104, and 103 against MAJIQ Heterogen; and 101, 98, and 99 against scQuint. The headline ordering is unchanged. Propagating all 128 null families through within-pair and across-pair Simes aggregation, two-fold conjunction, and BH gives zero replicated genes on every comparator-matched universe.

The larger whole-family audit reveals rare calls missed by the first 32 replicates. Four of 128 omnibus families make one BH call, one of 128 full pairwise families makes two, and 15 and two families make at least one fold-level call in pairwise folds 0 and 1, with maxima of two and four. No condition family makes a BH call. The omnibus complete-null frequency of any BH call is 3.1%, below 5%, but the fold-0 intermediate pairwise endpoint is anti-conservative at 11.7%. This does not propagate to the reported two-fold replication endpoint, which has zero calls in all 128 null families. Results and reusable audit scripts are under `analyses/null_replicate_sensitivity`.

## 2026-08-31, omnibus event classification and literature follow-up

Regenerated the biological event summary from the 584 BH-significant blocks in the current production omnibus family rather than the secondary pairwise localization family. The classes are 344 alternative-first regions, 97 cassette exons, 74 alternative-last regions, 30 whole-gene complex blocks, 21 compound internal blocks, and 18 other internal structures. Four of five pre-specified literature-audited events pass omnibus BH, `Ntrk2`, `Grin1`, `Dlg4`, and `App`, each at FDR 0.00152; `Enah` does not. Pairwise tests are no longer used to select Table 4 or the literature follow-up. The generated catalog and audit intersection are under `analyses/local_path_omnibus_event_annotation`.

## 2026-08-31, expanded omnibus literature screen

Screened the 584 BH-significant omnibus blocks for exact correspondence to experimentally studied splice events, using block coordinates, event structure, and exon length before reviewing primary functional evidence. Added four high-confidence matches to the literature audit. `Gria2` block B1 is the 115-nt flip/flop mutually exclusive pair, `Scn2a` block B2 is the 92-nt 5N/5A mutually exclusive pair, `Rbfox3` block B3 is the 93-nt RRM cassette, and `Mef2d` block B3 is the 21-nt beta microexon. Their omnibus FDR values are 0.00152, 0.0211, 0.0267, and 0.0323. Together with the four earlier matches, eight of nine audited events pass omnibus BH.

Excluded weaker candidates from the main functional-evidence table. The `Syt7` compound block plausibly overlaps known linker splice variants, but the strongest functional report is a recent preprint and the five-path fitted block does not uniquely isolate one published switch. The `Cacna1d` 8a/8b structure and `Grm5` 96-nt cassette are exact known events, but the located primary studies did not establish a clean isoform-specific functional consequence appropriate to the table. The short `Nrxn1` event does not map to canonical splice site 4. These exclusions avoid treating gene-level familiarity as event-level functional validation.

## 2026-08-31, pooled isoform baseline documentation

Clarified the construction of the production local-path baseline in `docs/differential.tex`. For each gene and included observation set, EC counts are summed across observations separately within each primer. One isoform composition is then estimated by maximizing the joint primer-specific conditional multinomial likelihood in reference-logit coordinates. Cell-type and condition labels do not enter this fit, so the baseline is unchanged under matched null-label transformations. It supplies fixed within-path isoform ratios, nuisance-isoform abundances, and represented mass for the subsequent block-specific path perturbations.

## 2026-08-31, uniform path-pseudocount sensitivity

Re-merged the saved \(\alpha=0\) and \(\alpha=0.5\) subject-split path fits using the final robust variance moderation and leave-one-block-out empirical calibration, making them comparable to the existing production \(\alpha=1\) and \(\alpha=2\) results. Fold-specific BH calls increase from 22 and 24 at \(\alpha=0\), to 465 and 258 at 0.5, 803 and 510 at one, and 1,163 and 770 at two. Replicated genes on the LeafCutter, MAJIQ Heterogen, and scQuint matched universes are respectively 1/0/0, 11/14/17, 31/53/56, and 50/103/101. Pooled fold-specific null rejection remains approximately 0.05 for every setting, and no setting makes a BH call in any of its 64 fold-level sign-flip null families. The unsmoothed fit retains only 7,295 and 4,911 eligible tests because of boundary instability, compared with 7,523 and approximately 5,057 for smoothed fits.

Added the common-pipeline sensitivity table to the main reproducibility results in `docs/differential.tex` and recorded its machine-readable values under `analyses/path_alpha_sensitivity`. Add-two is retained because it maximizes independently replicated genes on all three matched universes, but its selection used these subject splits and remains internally tuned.

Extended the uniform path-pseudocount sensitivity through \(\alpha=128\). Replicated genes on the LeafCutter, MAJIQ Heterogen, and scQuint matched universes increase from 50/103/101 at \(\alpha=2\), to 56/132/136 at 4, 70/141/147 at 8, 75/157/159 at 16, and an exploratory maximum of 78/161/161 at 32. They then decrease to 78/150/154 at 64 and 77/147/148 at 128. Fold-specific BH calls peak at 1,471 and 1,059 around \(\alpha=16\) and decrease to 1,250 and 1,004 at 128. Pooled sign-flip rejection remains approximately 0.05 and all 64 fold-level null families have zero BH calls at every setting, so the turnover reflects reproducible signal rather than loss of calibration. Fits at \(\alpha=256\) became substantially slower before completing a shard and were stopped because the decline was already established. The production benchmark remains at \(\alpha=2\), which was selected on the original grid; \(\alpha=32\) is explicitly reported as a post-hoc sensitivity maximum rather than substituted into the headline comparison.

## 2026-09-01, fixed-total path-smoothing concentration

Added an explicit concentration scaling to local-path fitting. The existing `per_path` rule retains total uniform prior concentration \(S_b\alpha_{\rm path}\) for a block with \(S_b\) paths, while the new `total` rule fixes concentration \(A\) and assigns \(A/S_b\) pseudo-observations to each path. The command-line runner records the selected scaling, the merger reports it, and the existing rule remains the default so saved production results are reproducible. Relevant differential and EC-block tests pass, 65 in total.

Reran the complete calibrated pairwise subject-split benchmark at fixed total concentrations \(A\in\{2,4,8,16,32,64,128,256\}\). The comparison matches each value to the previous \(\alpha_{\rm path}=A/2\), making the penalties identical for the 5,295 of 7,671 fold-0 candidates with two paths. Fixed-total peak fold discovery rates are 19.56% and 21.18%, compared with 19.55% and 20.95% under current scaling. Both curves reach 78, 161, and 161 replicated genes on the LeafCutter, MAJIQ Heterogen, and scQuint matched universes. At \(A=64\), fixed total gives 76/161/161 versus 78/161/161 under matched current scaling. At \(A=128\), fixed total gives 78/156/157 versus 78/150/154. Held-half replication is also comparable, and all fixed-total settings retain pooled null rejection near 0.05 with zero BH calls in every fold-level null family. Fixed total concentration is preferred on dimensional grounds without a material empirical penalty, but the headline and full-data results retain their original per-path fit because the scaling comparison is post-hoc on the evaluation splits. Machine-readable results are under `analyses/path_alpha_scaling`.

## 2026-09-01, corrected global empirical-Bayes smoothing

Replaced the former high-finite-concentration spike with an exact point null at zero path perturbation and removed the arbitrary finite-slab cap at eight. The candidate grid is now 0.5 through 128. Selection pools block path dimensions with equal weighting by gene. The fitted point-null weight only robustifies magnitude selection; downstream differential testing uses the selected finite scalar and does not use spike membership or a posterior probability of differential splicing. Corrected the Laplace slab evidence across dimensions by including (+\tfrac12\log S_b), the ILR-to-simplex Jacobian constant for an (S_b)-path composition.

The first corrected benchmark still made a conceptual error. It applied both the selected scalar and the baseline-centered evidence prior to downstream testing, whereas the alpha sensitivity used a uniform-centered penalty. Baseline-centered alpha 4--32 is not comparable with uniform-centered alpha 4--32. The resulting 4/0/0 replicated-gene counts were invalid evidence against EB smoothing. Corrected alpha-only transfer applies the EB-selected scalar to the existing uniform-centered production test.

With five-fold gene cross-fitting, baseline evidence selects alpha 4--8 per path or fixed-total concentration 8--32. Alpha-only transfer gives 70/137/142 replicated genes under per-path scaling and 65/135/139 under fixed-total scaling, compared with 50/103/101 at production alpha two. It retains all eligible fits, pooled null rejection near 0.05, and zero fold-level null-family BH calls.

Also tested uniform-centered evidence and skipping gene-fold cross-validation. Without cross-validation, baseline evidence selects per-path alpha four and eight in subject halves 0 and 1 and fixed-total concentration 16 in both. The fixed-total global result gives 1,442/1,017 fold BH calls, 68/139/145 replicated genes, held-half replication 0.845/0.854/0.861, and fold rank correlations 0.771/0.767/0.785. Uniform evidence selects much weaker values, alpha 0.5--2 per path or total concentration two to four. Its global results are 19/31/34 and 40/65/62 replicated genes under the two scalings; cross-fitting gives 20/41/45 and 40/65/67. This initially motivated baseline-evidence alpha transfer, but the center inconsistency is resolved by the testing-level analysis below. Results are under `analyses/path_smoothing_eb`.

## 2026-09-01, center-consistent testing-level empirical Bayes

Investigated the mismatch between baseline-centered count evidence and uniform-centered testing. Applying the same baseline-centered prior at both stages does not retain power. A pooled-baseline isotropic Gaussian prior on ILR perturbations selected precision 32 and produced 23/11 fold BH calls among 6,900/4,595 eligible tests. A heavier-tailed symmetric Dirichlet prior on relative path fold changes selected total concentration 64 and produced 31/14 calls among 7,423/4,961 eligible tests. Their pooled null rejection rates remained approximately 0.05 and no null family made a BH call. The shared biological center attenuates cell-type contrasts even when the log-ratio curvature is isotropic or the prior has linear tails. These failed prior-family options were removed from the code after the diagnostic runs.

Replaced the mismatched count-evidence transfer with testing-level empirical Bayes over the existing uniform-centered fixed-total fits. For each concentration, fitted a gene-equal beta--uniform mixture to the calibrated p-values, with one uniform null component and one decreasing beta alternative. The separate assessment halves select concentrations 32 and 64; pooling their evidence selects one global concentration 32. Repeating the complete selector on every matched null family selects concentration two in all 32 replicates, gives mean fold rejection rates 0.04995, and makes no BH call. At selected concentration 32, Tealeaf has 1,471/1,047 fold BH calls and 73/156/160 replicated genes against LeafCutter, MAJIQ Heterogen, and scQuint, compared with 68/139/145 for the earlier baseline-evidence transfer. The selector uses no fitting-data split and uses the same uniform-centered estimator for selection and testing. Its assessment remains internal because the same benchmark dataset estimated the global hyperparameter.

## 2026-09-02, promoted testing-level EB in the manuscript

Corrected `docs/differential.tex` so the main Tealeaf method uses the uniform-centered fixed-total estimator and selects one global concentration with the gene-equal beta--uniform testing-level EB objective. The main comparator result now reports the selected $A=32$ configuration, 1,471/1,047 fold BH calls, 73/156/160 replicated genes, held-half replication 0.878/0.888/0.875, and fold rank correlations 0.790/0.791/0.813. The full selector is rerun inside each matched null family, all 32 null families selected $A=2$, their mean nominal-0.05 rejection was 0.04995, and none made a BH call.

Removed claims that add-two per path is the production comparator method and moved duplicated testing-level EB derivation out of the exploratory appendix. The existing full-data omnibus, condition, event-type, and literature tables have not been rerun at fixed-total $A=32$; the manuscript and their analysis README now identify them explicitly as archived add-two fits rather than silently relabeling them as EB results.

## 2026-09-02, full-data testing-level EB rerun

Reran the full-data omnibus cell-type, pairwise cell-type, and condition-within-cell-type analyses across fixed-total concentrations (A\in\{2,4,8,16,32,64,128,256\}). An initial omnibus launch used an obsolete broad 10-UMI, three-subject candidate cache. I caught the mismatch before reporting it and reconstructed the production screen with at least 25 modeled-gene UMIs, four subjects per represented cell type, and equivalence deduplication. This gives the expected 1,785 omnibus candidates, of which 1,781 are eligible after local fitting.

The gene-equal beta--uniform testing-level selector independently chooses (A=64) for all three observed analysis families. The evidence curves turn over after 64, omnibus mean log evidence is 1.5453 at 64 and 1.5253 at 128, pairwise evidence is 0.48517 and 0.47951, and condition evidence is 0.006252 and 0.006076. The selected omnibus analysis has 894 nominal and 791 BH blocks from 572 genes. The secondary pairwise family has 4,233 nominal and 2,950 BH contrasts across 716 blocks and 531 genes. The condition family has 1,243 nominal tests and no BH event among 14,833 eligible hypotheses.

Reran the complete selector inside all 32 design-matched null families rather than auditing only the selected observed-data fit. Mean nominal-0.05 rejection is 0.05043 for omnibus, 0.05000 for pairwise, and 0.04995 for condition, with no BH call in any selected null family. The selected omnibus, pairwise, and condition fits require 0.488, 1.240, and 1.670 aggregate one-core worker hours; their complete eight-value grids require 5.020, 11.211, and 12.801 worker hours. Peak worker memory is at most 0.500 GiB. Regenerated structural annotations from the 791 omnibus discoveries, including 462 alternative-first regions, 139 cassette exons, 94 alternative-last regions, 41 whole-gene complex blocks, 27 compound internal blocks, 14 multi-exon-skipping blocks, nine mutually exclusive-exon blocks, three complex internal blocks, and two alternative-acceptor blocks. The eight previously audited literature matches remain significant.

## 2026-09-02, primer-separated alignment support

Generated sashimi-style alignment plots for three strong, structurally matched omnibus discoveries, \emph{App} exon 15, \emph{Gria2} flip/flop, and \emph{Grin1} C1. The workflow streams all 16 production spliced-pseudobulk BAMs, which contain unique genomic mappings deduplicated by cell barcode, UMI, locus, start, and CIGAR. Half-cell barcodes assign each alignment to poly(dT) or random-hexamer priming. Coverage and junction support are pooled within the localized cell-type comparison and primer, then reported per 1,000 primer-specific half-cells. Long introns are compressed and the annotation shows only the local paths tested by Tealeaf.

App contains 2,210/5,122 cortical-excitatory and 531/1,699 midbrain-inhibitory spliced UMIs under poly(dT)/random-hexamer priming. Gria2 contains 1,039/904 and 167/189. Grin1 contains 555/384 medium-spiny and 148/330 midbrain-inhibitory spliced UMIs. Both primer strategies support the tested structures, with the expected differences in positional coverage and absolute support. The plots remain descriptive because they omit unspliced, exonic-only, and ambiguous EC observations used by the production likelihood.

## 2026-09-02, sashimi horizontal alignment correction

Corrected the bundled ggsashimi grob assembly so every primer/cell-type coverage track, the shared genomic x-axis, and the local-path annotation receive the same complete vector of gtable column widths. The previous legacy adjustment modified only selected y-label columns, which shifted plotting-panel origins when track-label widths differed under the current ggplot2 release. Regenerated the App, Gria2, and Grin1 vector PDFs and rebuilt the manuscript.

## 2026-09-02, effect-size-selected sashimi contrasts

Changed the primer-separated sashimi workflow to select contrasts algorithmically from the selected fixed-total A=64 production pairwise fit. For each plotted block, the workflow first restricts to pairwise tests with BH FDR below 0.05 and then maximizes mean_difference_norm, the Euclidean norm of the across-subject mean paired difference in local-path ILR coordinates. This selects EX_cortical versus EX_hippocampus_pyramidal for App, effect norm 0.92698 and FDR 0.02109, EX_cerebellum_granule versus EX_cortical for Gria2, effect norm 0.45543 and FDR 0.00832, and INH_medium_spiny versus INH_midbrain for Grin1, effect norm 0.08811 and FDR 0.00577. Regenerated primer-specific alignment summaries and figures, added selected_contrasts.tsv as an audit table, and updated the manuscript.

## 2026-09-02, fitted path-usage panels

- Increased the three primer-separated sashimi track heights from 0.8 to 1.1 and annotation height from 0.7 to 0.9 while retaining the shared horizontal layout.
- Exported subject-level fitted path proportions for the three strongest significant pairwise localizations by rerunning the production local-path likelihood with the selected EB prior setting. All three fits completed without failure.
- Added plotnine grouped barplots of mean path usage by cell type with standard errors across paired subjects. App used four paired subjects, while Gria2 and Grin1 each used nine.
- Added the bars below the read-support panels in `docs/differential.tex` and clarified that primer-separated tracks are descriptive while path usage is jointly fitted from both primer likelihoods.

## 2026-09-02, cell-type evidence heatmaps

- Added plotnine heatmaps of fitted path usage, observed junction coverage, and tested-exon coverage for each cell type in the three displayed contrasts.
- Kept junction and exon evidence separate by priming strategy. Coverage is normalized per 1,000 primer-specific half-cells and displayed on a log1p scale.
- Summarized exon evidence as mean normalized base coverage over each unique tested exon and retained all observed junctions. Event-level feature-coverage tables map J and E labels to genomic coordinates and untransformed values.
- Integrated the heatmaps into `docs/differential.tex`; reusable feature summaries are in `tealeaf/sc/sashimi.py`.

## 2026-09-02, path-ordered block evidence heatmaps

Replaced the separate pairwise path-usage, all-junction, and unique-exon heatmaps with two combined block matrices per event, one for poly(dT) and one for random-hexamer priming. Each matrix includes every cell type represented in the production omnibus fit, nine for App, eight for Gria2, and five for Grin1. A path marker contains the across-subject mean fitted path proportion from the selected fixed-total A=64 omnibus likelihood. Its following columns alternate exon and exact annotated-junction evidence along that tested path. Shared features are deliberately duplicated under every path containing them; off-path and out-of-block features are excluded.

Component order is 5'-to-3' in RNA space. All three displayed genes are on the minus strand, so their heatmap coordinates decrease from left to right within each path. Unit tests cover both strands, exact exon and junction normalization, shared-feature duplication, and marker/component interleaving. Matrix validation confirmed two primers for every event, 22/12/14 total columns for App/Gria2/Grin1, and 8/2/4 genomic coordinates repeated across paths.

Junction evidence is annotated-junction UMI count per 1,000 primer-specific half-cells, while exon evidence is mean aligned-base depth over the exact exon per 1,000 half-cells. Since these and fitted proportions have different units, the display uses a within-feature z-score across cell types, computed separately by primer and truncated to [-2,2]. The event-level block-heatmap tables retain raw values, coordinates, feature types, path membership, and ordering. The superseded panels were archived outside the repository rather than deleted.

## 2026-09-02, gene-wide evidence labels and coverage scale

Revised the combined block matrices so E and J labels are assigned once from the strand-aware union of displayed features for each gene. Repeated coordinates now retain the same label across paths, and paths skip labels for features they do not contain. Within-path column order remains 5'-to-3' in RNA space.

Replaced within-feature z-score colors with raw primer-specific coverage. Exon tiles show mean aligned-base depth and junction tiles show exact annotated-junction UMIs, both per 1,000 primer-specific half-cells and on one raw scale within each primer panel. Path proportions are not coverage, so colored path markers use an independent 0-to-1 scale and colorbar while retaining fitted percentage labels. Regenerated the six vector panels and retained the raw values and gene-wide labels in the event matrices.

## 2026-09-02, isoform-switch search and PSI evidence

Exported subject-level path usage for all 1,781 eligible selected fixed-total \(A=64\) full-data omnibus fits and ranked the 791 BH-significant blocks by the largest across-cell-type range in mean fitted usage. Ninety-two blocks change the identity of the dominant path, but no path spans from at most 0.1 to at least 0.9. The manuscript therefore reports that no near-complete fitted switch was found and plots the three largest non-whole-gene dominant-path changes that contain a splice junction, Agtpbp1 B1, Dlg4 B1, and Ntrk2 B1. Their largest path ranges are 0.326--0.639, 0.179--0.430, and 0.255--0.480.

Replaced raw coverage in these switch panels with empirical exon and junction PSI. Variable-exon PSI is pooled exon depth divided by the mean depth of the nearest path-shared exon on each available side, clipped to \([0,1]\); shared exons equal one by definition. Junction PSI is a junction's pooled UMI count divided by the total count of unique observed junctions sharing either donor or acceptor. BAM collection now retains competitors whose other splice site lies outside the displayed interval, and a unit test covers that case. Path usage and empirical PSI use separate fixed 0--1 colorbars. The complete ranking, selected estimates, empirical matrices, manifest, and six vector plots are under analyses/isoform_switch_examples.

## 2026-09-03, switch-panel scale and smoothing diagnostic

Changed the switch evidence panels to one shared 0--1 viridis scale for fitted path usage, exon PSI, and junction PSI. Exons and junctions present in every displayed path are omitted because they carry no information about the local path choice; shared exons remain internal normalization references. Explicit path-order records preserve path markers when removing all evidence features from one path.

Audited whether the absence of a fitted 0.1--0.9 switch under selected fixed-total \(A=64\) was suspicious. Over the 791 production-significant blocks, \(A=64\) fits all but four blocks, finds 92 dominant-path changes and no near-complete cell-type mean switch. Fixed-total \(A=2\) similarly finds 88 dominant changes and no near-complete mean switch. The unsmoothed fit yields three apparent near-complete mean switches, but loses 99 of 1,785 candidate fits and puts the 1st/5th individual-usage percentiles at zero and the 95th/99th percentiles at one. Two of its three apparent switches use a low-usage cell type represented by only one subject. The production \(A=64\) values are therefore suitable as regularized estimates associated with the calibrated discovery procedure but cannot establish the biological absence of complete switches. The diagnostic is recorded in `analyses/isoform_switch_examples/smoothing_sensitivity.tsv`.
