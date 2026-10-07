# Real-design within-path nuisance sensitivity

This analysis samples actual EC counts using observed gene designs and primer/depth totals for 32 covered production-eligible paired hypotheses, selected without discovery or long-read outcome screening. Each of two count draws retains the same true target path usage across cell types within each subject. Type-specific within-path transcript ratios are tilted with log-normal standard deviation one and renormalized to preserve every path mass and outside-block abundance. No-tilt controls use the same hypothesis family. Each test uses the existing paired ILR statistic and its 32 sign-flip pooled-calibration families; this is a count-null diagnostic, not a full split/LR benchmark or joint-gene FDR certification.

| Model | Concentration | No-tilt null calls | Within-path-change null calls |
| --- | --- | --- | --- |
| Fixed within-path shares | 1 | 9/64 | 11/64 |
| Free within-path shares | 1 | 8/64 | 11/64 |
| Fixed within-path shares | 32 | 19/64 | 19/64 |
| Free within-path shares | 32 | 19/64 | 19/64 |

All trials converge. Failed trials would remain p=1, and collation validates every requested test/draw/strategy combination. Only 11 hypotheses have multiple transcripts in at least one path; 21 have singleton paths, where this tilt cannot change transcript composition. Among the 22 multiplet-path trials, concentration-1 calls increase from 2 to 4 for fixed shares and from 2 to 5 for free shares. These small samples do not establish a difference between estimators. No-tilt controls already have excess rejection, showing that within-path misspecification is not the sole source of the quantification/testing problem.

The controlled four-transcript challenge in `analyses/independent_block_nuisance` demonstrates an identifiable structural failure, but its free-share repair does not suffice in this real-design screen. No production test, main comparator result, or full-data interpretation is replaced. Hyperparameters are not selected using these outcomes or long-read agreement.

The reusable simulator is `tealeaf.sc.path_simulation.simulate_counts`, with the optional `within_path_type_scale`. The existing count-null runner and summarizer are reused. `extra_scripts.compare_within_path_type_stress` validates matched families and retains per-test results and multiplet/singleton structural strata.
