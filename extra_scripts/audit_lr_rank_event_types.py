"""Describe event types on existing method-own LR rank curves, without reranking."""

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from tealeaf.sc.replication_audit import ranked_category_summary, ranked_direction_summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--published-ranks", type=Path, required=True)
    parser.add_argument("--control-ranks", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    published = pd.read_csv(args.published_ranks, sep="\t", low_memory=False)
    # These sources carry explicit known event labels. Junction clusters and
    # path blocks require a separate classifier and are not silently called SE.
    methods = ["SUPPA2 (full data)", "SUPPA2 (primer aware)", "Tealeaf/SUPPA2 hybrid", "rMATS paired JCEC"]
    published = published.loc[published.method.isin(methods)].copy()
    if set(published.method) != set(methods):
        raise ValueError("published event-defined method curves are incomplete")
    control = pd.read_csv(args.control_ranks, sep="\t")
    if control.method.nunique() != 1 or control.event_type.isna().any():
        raise ValueError("control must contain one labeled method-own ranked family")
    table = pd.concat([published, control], ignore_index=True, sort=False)
    cutoffs = (100, 200)
    categories = pd.DataFrame(ranked_category_summary(table, cutoffs=cutoffs))
    overall = pd.DataFrame(ranked_direction_summary(table, cutoffs=cutoffs))
    if not categories.complete_prefix.all():
        raise ValueError("audit requires complete top-100 and top-200 prefixes")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    categories.to_csv(args.output_dir / "event_type_composition.tsv", sep="\t", index=False)
    overall.to_csv(args.output_dir / "overall_rank_summary.tsv", sep="\t", index=False)
    paths = (args.published_ranks, args.control_ranks)
    (args.output_dir / "manifest.json").write_text(json.dumps(dict(inputs=[dict(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for path in paths], selection="existing method-own global top-100/top-200 LR-evaluable prefixes, no event-type reranking or significance cutoff", scope="completed historical curves, not the pending sequence-map candidate", caveats="event categories are descriptive, different quantification and tested families remain confounded; the original binary control failed count-null checks, its discoveries are not valid power", production_changes=False), indent=2) + "\n")
    print(overall.to_string(index=False), flush=True)
    print(categories.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
