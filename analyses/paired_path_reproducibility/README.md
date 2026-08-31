# Paired local-path reproducibility benchmark

This directory records the selected add-two paired local-path benchmark. Each subject-by-cell-type path composition is fit directly to both primer-specific EC count vectors around a pooled isoform baseline. Scalar paired tests use robust empirical-Bayes variance moderation. All tests are calibrated using \(R_{\rm null}\) paired sign-flip Monte Carlo replicates by a leave-one-event-out empirical CDF within path-dimension and paired-subject-count strata. The default is \(R_{\rm null}=32\); calibration pools null values across other blocks in the same stratum, whereas the whole-family BH audit has only \(R_{\rm null}\) realizations.

`gene_metrics.tsv` contains the matched split-subject comparison reported in `docs/differential.tex`. `calibration.tsv` summarizes fold-level block-test calibration. `replicated_null_audit.tsv` propagates each sign-flip family through the reported Simes aggregation, two-fold conjunction, and BH procedure. `resources.tsv` records fitting-only Slurm accounting for the two subject folds.

The nested 64- and 128-replicate checks are recorded in `analyses/null_replicate_sensitivity`. Replicated-gene counts change by at most three, and none of 128 complete null families produces a replicated gene on any comparator-matched universe.

The smoothing strength was selected after examining these subject splits. The table therefore demonstrates internally replicated performance on this dataset, while independent-data confirmation remains necessary.
