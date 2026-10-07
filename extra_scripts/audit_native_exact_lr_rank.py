"""Cross cached SUPPA2 rankings and effects without assuming equal quantification."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.audit_event_information import cross_ranking_effects


def eligible(path, method):
    table = pd.read_csv(path, sep="\t")
    table = table.loc[table.method.eq(method)]
    valid = table.mapping_complete.astype(str).str.lower().eq("true") & table.minimum_pooled_depth.ge(20) & table.pooled_replicated.notna()
    return table.loc[valid].copy()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--published-mapping", type=Path, required=True)
    parser.add_argument("--exact-mapping", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    published = eligible(args.published_mapping, "SUPPA2 (full data)")
    exact = eligible(args.exact_mapping, "SUPPA2 (full data); hybrid-exact paired Wilcoxon")
    summary, shared, _, _ = cross_ranking_effects(published, exact)
    labels = {"native_p": "published normal-approximation ranking", "hybrid_p": "cached exact-test ranking", "native_effect": "published expression effect", "hybrid_effect": "cached exact expression effect"}
    summary["ranking"] = summary.ranking.map(labels)
    summary["direction"] = summary.direction.map(labels)
    summary["method"] = summary.method.str.replace("native_p", "published_p", regex=False).str.replace("native_effect", "published_effect", regex=False).str.replace("hybrid_p", "exact_p", regex=False).str.replace("hybrid_effect", "exact_effect", regex=False)
    summary["scope"] = summary.scope.replace({"native own universe": "published own universe", "hybrid own universe": "cached exact own universe"})
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output_dir / "rank_decomposition.tsv", sep="\t", index=False)
    delta = np.abs(shared.native_effect - shared.hybrid_effect)
    equal = bool(np.allclose(shared.native_effect, shared.hybrid_effect, rtol=1e-10, atol=1e-10))
    settings = {"n_shared": len(shared), "maximum_effect_difference": float(delta.max()), "effects_equal_at_1e10": equal, "external_truth": "identical event membership and source-count truth checked before crossing", "eligibility": "all finite nonzero short directions at pooled LR depth20; exact-zero LR effects remain nonagreement", "scope": "diagnostic; cached exact source includes a broader event family and potentially different quantification, not silently replacing the published baseline", "interpretation": "not a pure tail-calculation ablation when cached effects differ", "production_changes": False}
    shared.columns = [column.replace("native_", "published_").replace("hybrid_", "exact_") for column in shared.columns]
    shared.to_csv(args.output_dir / "shared_events.tsv.gz", sep="\t", index=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(settings, indent=2) + "\n")
    print(summary.loc[summary.cutoff.eq(100)].to_string(index=False))
    print(json.dumps(settings, indent=2))


if __name__ == "__main__":
    main()
