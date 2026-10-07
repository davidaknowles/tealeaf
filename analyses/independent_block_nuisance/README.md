# Independent splice-block nuisance challenge

This controlled simulation tests whether a B-block splice change can bias inference for a separate A block in the same gene. Four transcript states are the Cartesian product of two binary blocks. True A usage is subject dependent, but equal between cell types under the A null. B usage can change independently. Counts are sampled from known primer-specific transcript-to-read-origin maps. There are 12 subjects, two types, and 100/200 reads for the two primers. A label-blind generating-mean transcript baseline is used as a diagnostic oracle, not an available real-data estimator.

Every read selection uses the same experimental random-subject Beta EC likelihood-ratio test. The all-gene selection fixes within-A transcript shares. Target-local selections retain only A-origin junction/exon classes or A-origin junction classes. These synthetic classes preserve read origin even where compatible transcript sets coincide. Conventional transcript EC counts do not retain this origin distinction. All 64 requested trials remain in each denominator; optimization or integration failures would receive p=1 with no direction.

| Scenario | All gene, fixed shares | Local junction + exon | Local junction |
| --- | --- | --- | --- |
| Both blocks null | 1/64 calls | 1/64 calls | 1/64 calls |
| Only B changes, A remains null | 64/64 false calls | 2/64 calls | 2/64 calls |
| A and B change | 64/64 calls | 64/64 calls | 64/64 calls |

All fits pass the numerical gate. With only B changing, the all-gene model's mean A-effect bias is -0.110, compared with 0.0053 and 0.0038 for target-local selections. With both changing, mean biases are -0.0506, 0.0033, and 0.0053; all directions are correct. This strong-effect control does not establish a power advantage for exonic reads.

These results identify a structural failure of fixed within-path transcript mixtures. They do not establish the prevalence of this failure in real genes or a strategy unbeaten on both split discovery and independent long-read agreement. No production defaults or main comparator results change. The all-gene model includes off-block read classes with opportunities depending on both blocks; a correctly specified joint model could use these reads without the frozen-mixture bias.

The reusable generator is `tealeaf.sc.path_simulation.simulate_independent_binary_blocks`. `extra_scripts.audit_independent_block_nuisance` runs the three read selections, and `extra_scripts.summarize_independent_block_nuisance` validates complete draw/strategy families before summarizing. The separate `extra_scripts.audit_independent_block_estimators` sensitivity tests existing fixed-share and free-share quantifiers with one common unmoderated paired ILR test. Its test is not the production pooled calibration.

The estimator sensitivity is in `estimators_linear_cov`. With only B changing, existing fixed-share fits at total concentration 1 and 32 make 64/64 false calls. Free-share fits at concentration 1 make 4/64, with mean effect bias 0.0097. All 64 trials converge for all three estimators; all have 1/64 calls when both blocks are null. These are strict target-null diagnostics, not certification under residual biological variation.

The initial `estimators` sensitivity preserves a numerical failure for inspection, all free-share abundance fits optimized successfully but their covariance rank checks failed near a nuisance transcript boundary. It is not evidence of a zero false-positive rate. The free-share fitter now computes the same working covariance in linear simplex-tangent coordinates, without changing the objective or point estimates. Independent interior-coordinate and boundary regression tests cover the repair. The coordinate change does not validate an exact posterior uncertainty approximation.
