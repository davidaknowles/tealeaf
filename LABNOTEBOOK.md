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
