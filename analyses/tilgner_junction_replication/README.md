# Junction-method replication in independent long reads

This directory compares full-data cell-type discoveries from Tealeaf, LeafCutter, relaxed-coverage MAJIQ Heterogen, scQuint, and a paired junction-CLR ablation with the P56 mouse-brain ScISOr-Seq2 transcript-UMI matrix from Joglekar et al. (2024). Every method is evaluated by directional agreement with the same independent long-read source under common mapping and depth rules. These are conditional validation rates, not a matched comparison of power.

The MAJIQ result uses a combined build requiring support in at least three pseudobulks and PSI coverage of at least three reads in at least three bins. Heterogen requires coverage in at least half the pseudobulks within each compared group. Across the 28 contrasts among the eight mapped cell types, MAJIQ tests 10,297 edges and calls 2,657 edges from 466 genes. Significant edges are restored to multijunction LSVs before estimating and validating direction.

`method_comparison_summary.tsv` contains the headline comparison, and `method_comparison.pdf` visualizes it. `junction_replication.tsv.gz` contains the junction-method audit rows, `junction_replication_summary.tsv` summarizes them, `paired_clr_discoveries.tsv.gz` caches the CLR discoveries, `majiq_relaxed_contrasts.json` records the rerun contrast manifest, and `manifest.json` records analysis choices.
