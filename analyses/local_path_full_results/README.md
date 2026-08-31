# Full-data local-path results

These tables summarize the production direct local-path analysis reported in `docs/differential.tex`. `summary.tsv` contains observed and design-matched null results after empirical calibration and BH correction. `runtime.tsv` reports aggregate one-core fitting time and peak worker memory from prepared analysis inputs. Structural annotations for the significant cell-type blocks are in `analyses/local_path_full_event_annotation`.

Production uses two pseudo-observations per path, a label-independent pooled isoform baseline, paired mouse-level testing for cell types, independent mouse-level multivariate regression for condition within cell type, and 32 design-matched null families. Predictive empirical-Bayes smoothing and conditional Laplace covariance propagation are retained as appendix sensitivities because they did not improve subject-split replicated discovery.
