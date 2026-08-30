# Paired local-path reproducibility benchmark

This directory records the selected add-two paired local-path benchmark. Each subject-by-cell-type path composition is fit directly to both primer-specific EC count vectors around a pooled isoform baseline. Scalar paired tests use robust empirical-Bayes variance moderation. All tests are calibrated against 32 paired sign-flip null families by a leave-one-event-out empirical CDF within path-dimension and paired-subject-count strata.

`gene_metrics.tsv` contains the matched split-subject comparison reported in `docs/differential.tex`. `calibration.tsv` summarizes fold-level block-test calibration. `replicated_null_audit.tsv` propagates each sign-flip family through the reported Simes aggregation, two-fold conjunction, and BH procedure. `resources.tsv` records fitting-only Slurm accounting for the two subject folds.

The smoothing strength was selected after examining these subject splits. The table therefore demonstrates internally replicated performance on this dataset, while independent-data confirmation remains necessary.
