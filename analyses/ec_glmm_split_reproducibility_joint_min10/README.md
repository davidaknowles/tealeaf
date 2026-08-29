# Low-depth joint-gene EC-GLMM split sensitivity

This directory records the exploratory subject-split comparison that lowers the modeled-gene UMI threshold from 25 to 10 and replaces separate annotated-block tests with one unrestricted isoform-composition LRT per gene and cell-type pair. The fold manifests contain 22,453 tests from 2,503 genes and 16,753 tests from 2,238 genes.

The larger hypothesis universe raises reproducible discovery yield and Tealeaf exceeds relaxed MAJIQ, but it remains below LeafCutter and scQuint. Tealeaf has higher held-half nominal replication than all three comparators. The paired-label-swap audit is anti-conservative, only 946 of 1,000 fits jointly converge, rejection rates are 0.0909, 0.0222, and 0.00634 at nominal levels 0.05, 0.01, and 0.001, and three tests pass BH at 0.05. This configuration is therefore a rejected sensitivity, not a production analysis.

`gene_metrics.tsv` contains the matched discovery and reproducibility summaries. `null_calibration.tsv` records the complete label-swap audit.
