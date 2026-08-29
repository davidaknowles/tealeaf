# Local splice-path power audit

This analysis investigates why the pairwise Tealeaf EC-count GLMM did not produce more split-replicated discoveries than junction-based comparators, with cassette exons used as the diagnostic subset.

## Failure in the global-isoform null

The previous block test represented both the tested path effect and its nuisance complement in the Euclidean geometry of global isoform logits. That complement is not the null space of the aggregate local path log-ratio when two isoforms on one path have unequal baseline abundance. A nuisance coefficient can therefore change the quantity that the test is meant to hold fixed.

For path $s$, let $I_s$ be its isoform set and let $\boldsymbol\theta\in\Delta^{T-1}$ be the $T$-dimensional gene isoform composition. The local abundance is $\psi_s=\sum_{t\in I_s}\theta_t$. The derivative of a two-path log-ratio with respect to isoform logit $\eta_t$ is weighted by the conditional abundance $\theta_t/\psi_s$, not by an unweighted path indicator. Euclidean orthogonality to the lifted path indicator therefore does not imply a zero derivative of $\log(\psi_1/\psi_2)$.

In fold 0, 1,081 of 1,353 cassette tests have nonzero leakage of the path log-ratio into the old global null. The median nuisance-derivative norm relative to the tested derivative norm is 0.341, the 90th percentile is 0.522, and the maximum is 0.817. Fold 1 gives 653 of 859 tests with nonzero leakage, median 0.330, 90th percentile 0.523, and maximum 0.816. In fold 0, leakage and the corrected model's gain in −log10 p-value have Spearman correlation 0.236; the highest-leakage quartile has median gain 0.213. `summarize_null_geometry.py` reproduces this audit from a prepared data cache and candidate manifest.

A separate fitting audit found that reuse of an unrestricted gene alternative could occasionally return an objective worse than a later block null. This affected about 18% of cassette fits in the original subject halves. Refining a reused alternative from the fitted block null removes almost every nesting violation, although the correction alone does not improve cassette replication enough to explain the power gap.

## Local-path likelihood

The corrected model estimates a label-independent pooled isoform mixture, normalizes those weights within each represented path, and uses the weighted compatibility column for that path. It retains the original pseudobulk-by-EC count vectors, separate oligo(dT) and random-hexamer compatibility maps, EC ambiguity, multinomial covariance, the mouse random intercept, and the pseudobulk logistic-normal effect. It does not estimate global isoform proportions and then run a second-stage test.

The latent components are the represented local paths plus any isoform not represented in the block map. Mouse cassette and full pairwise manifests contain at most one such background precursor per test. A two-path cassette therefore has one path log-ratio and a one-degree-of-freedom pairwise LRT. The nuisance direction changes represented-versus-background abundance while preserving the represented path ratio exactly.

The one-degree-of-freedom statement is exact for the Tealeaf cassette test. A universal two- or three-degree-of-freedom statement is not exact for the external outputs because LeafCutter, MAJIQ Heterogen, and scQuint report different cluster, local-variation, and intron-group features, and the benchmark combines their feature p-values within a gene and cell-type pair. Tealeaf avoids testing several junction coordinates for one annotated path choice, but the size of that advantage is comparator- and event-dependent.

Conditioning on ECs compatible with only one represented path is not a viable replacement. In a 100-cassette audit, 24 events lack a path-specific EC for at least one path. Among the 76 feasible events, directly path-identifying ECs contain a median 96.5 UMIs and only 10.7% of crude gene UMI exposure. Shared-exon and precursor-compatible reads are not path-identifying merely because they are exonic. The local-path likelihood retains these observations and lets the primer-specific compatibility maps determine how much information they carry.

## Cassette validation

The corrected model was fit to every eligible cassette-exon pairwise test in two condition-balanced subject halves and to 1,000 paired-label swaps. Fold 0 has 1,311 jointly converged tests among 1,353 candidates, with 143 nominal tests and 29 fold-specific BH calls. Fold 1 has 829 jointly converged tests among 859 candidates, with 69 nominal tests and 11 fold-specific BH calls.

There are 732 exact event-by-cell-type-pair tests converged in both halves. Twenty-seven have (p < 0.05) in both halves and six pass BH on the conjunction p-value max(ρ_0, ρ_1). All six BH events have the same fitted effect sign across halves. Among the 27 both-nominal events, 25, or 92.6%, have the same sign and the effect-size Spearman correlation is 0.728. On the 711-test universe also converged under the old global model, the corrected model has 26 both-nominal and four conjunction-BH events, compared with 11 and two under the old model.

Among 971 converged cassette paired-label nulls, rejection rates are 0.0587, 0.0113, and 0.00206 at nominal levels 0.05, 0.01, and 0.001, with no BH calls. The former global model's completed native-pairwise audit, which was not restricted to cassettes, gave 0.0745, 0.0186, and 0.00414. The local-path result is close to nominal and provides no evidence that its replication gain comes from a looser reference distribution.

A global-isoform tangent formulation uses the pooled-baseline Jacobian of the local path ILR. Tested-label nuisance effects lie in its kernel and the tested basis is a right inverse. This retains within-path isoform nuisance changes while correcting the former null to first order. It has 1,319 and 824 converged cassette tests in the two halves. Among 730 exact common tests, 28 are nominal in both halves and six pass conjunction BH; 26 of 28 nominal effects agree in sign and their Spearman correlation is 0.722. Among 975 converged cassette paired-label nulls, rejection rates are 0.0615, 0.0133, and 0.00205 at nominal levels 0.05, 0.01, and 0.001, with no BH calls. Because this matches the exact local-path result without fixing within-path mixtures and remains close to nominal under label swaps, the tangent model is the production correction and the exact local model is retained as a sensitivity.

The external-method cassette-only score matches genes and cell-type pairs rather than genomic event coordinates, so a comparator can receive credit for a different event in the same gene. It is retained as a sensitivity, not as an event-specific power comparison. The native full test-universe comparison is the relevant primary benchmark.

## Full split benchmark

The production tangent model was rerun over all 7,671 and 5,140 native-screened block-by-cell-type-pair hypotheses in the two subject halves. It jointly converged for 7,471 and 4,967 tests. The two fits used 18.06 and 13.06 aggregate worker hours, compared with 51.79 and 40.43 hours for the former-basis split fits, a 2.87-fold and 3.10-fold reduction from analytic likelihood derivatives and safe reuse of the unrestricted gene alternative.

On the pairwise shared universes, Tealeaf has 12 versus 25 replicated genes against LeafCutter, 19 versus 38 against MAJIQ Heterogen, and 18 versus 57 against scQuint. The corresponding Tealeaf versus comparator held-half nominal replication rates are 0.867 versus 0.744, 0.828 versus 0.697, and 0.828 versus 0.777. Fold rank correlations are 0.562 versus 0.543, 0.362 versus 0.434, and 0.358 versus 0.461. The tangent correction improves Tealeaf's replicated counts from 10, 17, and 17, but it does not make Tealeaf more powerful than the junction methods on this genome-wide gene-aggregated benchmark. Scored outputs are under `analyses/ec_glmm_split_reproducibility_tangent_full_min4`.

The expected advantage is conditional rather than automatic. A cassette path comparison has one tested Tealeaf coefficient, but reads shared by both paths contribute depth and nuisance information rather than direct path-composition information after conditioning on the primer total. External tools also define and filter features differently before p-values are combined to genes. The exact event-level cassette analysis best isolates the intended one-degree-of-freedom comparison and shows a gain, whereas the genome-wide gene score mixes simple and compound blocks and can credit a comparator for any feature in the same gene and cell-type pair.

`summarize_split_validation.py` computes exact-event conjunction, direction, and null-calibration summaries. Intermediate fit directories and candidate caches are not part of the repository because they are reproducible analysis artifacts rather than source inputs.
